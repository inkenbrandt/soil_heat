"""soil_ground_heat_flux.py
====================================================
Python utilities implementing the equations from:

    Wang, Z.-H., & Bou-Zeid, E. (2012). *A novel approach for the estimation of
    soil ground heat flux*. *Agricultural and Forest Meteorology*, 154-155,
    214-221.

All equations appearing in that paper are reproduced as vectorised Python
functions, together with a few convenience helpers.  The library is fully
NumPy-aware and can be used on scalars or on time-series arrays.

Units
-----
All functions assume SI units **throughout**:

* depth ``z``            ― m  (positive downward)
* time ``t``             ― s
* temperature ``T``      ― K  (or °C provided the 0-offset is handled
  consistently)
* heat flux ``G, H, LE`` ― W m-2
* radiative flux ``Rn``  ― W m-2
* soil properties
    * thermal conductivity ``k`` ― W m-1 K-1
    * heat capacity ``rho_c``    ― J m-3 K-1
    * thermal diffusivity ``kappa = k / rho_c`` ― m² s-1

Dependencies
------------
>>> pip install numpy scipy

Example
-------
>>> import numpy as np, soil_ground_heat_flux as sghf
>>> # 30-minute series (dt = 1800 s) of flux-plate measurements at z = 0.08 m
>>> Gz = np.loadtxt('Gz_8cm.txt')
>>> G0 = sghf.estimate_G0_from_Gz(Gz, z_r=0.08, kappa=0.7e-6, dt=1800)

"""

from __future__ import annotations

import math
from typing import Callable, Sequence

import numpy as np
from scipy.special import erfc, gammaincc

__all__ = [
    # Energy balance & residuals
    "energy_balance_residual",
    "surface_energy_residual",
    # Conventional ground-heat-flux estimator (Eq. 2)
    "ground_heat_flux_conventional",
    # Heat-conduction fundamentals
    "green_function_temperature",
    "temperature_convolution_solution",
    "soil_heat_flux_from_G0",
    "estimate_G0_from_Gz",
    # Sinusoidal analytical solutions (Eqs. 13–15)
    "sinusoidal_boundary_flux",
    "soil_temperature_sinusoidal",
    "soil_heat_flux_sinusoidal",
    # Soil-property parameterisations (Eqs. 16–18)
    "heat_capacity_moist_soil",
    "pf_from_theta",
    "thermal_conductivity_moist_soil",
    "thermal_diffusivity",
]

# -----------------------------------------------------------------------------
# 1. Energy-balance bookkeeping (Eqs. 1 & 19)
# -----------------------------------------------------------------------------


def energy_balance_residual(
    Rn: float | np.ndarray,
    H: float | np.ndarray,
    LE: float | np.ndarray,
    G0: float | np.ndarray,
) -> float | np.ndarray:
    """Return the residual of the surface energy balance (Eq. 19).

    ``Res = Rn − G0 − H − LE``
    """
    return Rn - G0 - H - LE


# Alias used later in the module.
surface_energy_residual = energy_balance_residual

# -----------------------------------------------------------------------------
# 2. Conventional ground-heat-flux estimator (gradient + calorimetry) – Eq. 2
# -----------------------------------------------------------------------------


def ground_heat_flux_conventional(
    k: float,
    dT_dz_at_zr: float,
    rho_c: float,
    dT_dt_profile: Sequence[float],
    z_profile: Sequence[float],
) -> float:
    """Conventional estimate of *surface* ground heat flux, Eq. (2).

    Parameters
    ----------
    k
        Effective soil thermal conductivity **(W m-1 K-1)**.
    dT_dz_at_zr
        Vertical temperature gradient evaluated at the flux-plate depth ``z_r``
        **(K m-1)**.  A *negative* gradient means temperature decreases with
        depth.
    rho_c
        Volumetric heat capacity of the soil **(J m-3 K-1)**.
    dT_dt_profile
        Time derivatives ∂T/∂t for *each* node between the surface and ``z_r``
        **(K s-1)**.  Any iterable (list, ndarray, …).  Must align with
        ``z_profile``.
    z_profile
        Depth of each node in the temperature profile **(m)**.  Increasing,
        positive downward, **excluding** the surface (z = 0) but *including*
        ``z_r`` (last element).

    Returns
    -------
    G0 : float
        Ground heat flux at the surface **(W m-2)**.  Positive *into* the soil.
    """
    # Fourier conduction term (gradient method)
    G_conduction = -k * dT_dz_at_zr

    # Heat-storage (calorimetry) – numerical integration by trapezoid
    z = np.asarray(z_profile)
    dT_dt = np.asarray(dT_dt_profile)
    if z.shape != dT_dt.shape:
        raise ValueError("z_profile and dT_dt_profile must have same length")

    # Integration bounds: surface (0) to z_r (last node) – prepend surface
    z_nodes = np.concatenate(([0.0], z))
    dT_dt_nodes = np.concatenate(([dT_dt[0]], dT_dt))
    storage = np.trapezoid(dT_dt_nodes, x=z_nodes)
    G_storage = rho_c * storage

    return G_conduction + G_storage  # type: ignore


# -----------------------------------------------------------------------------
# 3. Heat-conduction fundamentals (Eqs. 3–12)
# -----------------------------------------------------------------------------


def green_function_temperature(z: float, t: float, kappa: float) -> float:
    """Green’s function **g_z(t)** for the semi-infinite 1-D heat equation (Eq. 7).

    Notes
    -----
    Returns **0** when *t ≤ 0* (causality).
    """
    if t <= 0:
        return 0.0
    return 2.0 / math.sqrt(math.pi) * math.sqrt(kappa * t) * math.exp(
        -(z**2) / (4.0 * kappa * t)
    ) - z * erfc(z / (2.0 * math.sqrt(kappa * t)))


def temperature_convolution_solution(
    z: float, t_series: np.ndarray, f_series: np.ndarray, kappa: float, Ti: float = 0.0
) -> np.ndarray:
    """Temperature time-series at depth *z* via Duhamel convolution (Eq. 6).

    ``T(z,t) = Ti + ∫ f(t-τ) d g_z(τ)``

    The integral becomes a discrete convolution where *f* is the boundary
    heat-flux series (W m-2  → ∂T/∂z via Fourier).
    """
    if t_series.ndim != 1 or f_series.ndim != 1:
        raise ValueError("t_series and f_series must be 1-D arrays of equal length")
    if t_series.size != f_series.size:
        raise ValueError("t_series and f_series must be the same length")

    dt = np.diff(t_series)
    if not np.allclose(dt, dt[0]):
        raise ValueError("Time vector must be uniformly spaced")
    dt = dt[0]

    g = np.array([green_function_temperature(z, t, kappa) for t in t_series])
    dg = np.diff(g, prepend=0.0)  # discrete derivative → Stieltjes measure

    # Convolution implementation: cumulative sum of f * dg (causal)
    T = Ti + np.cumsum(f_series[::-1] * dg)[::-1]  # reversed for (t-τ)
    return T


def soil_heat_flux_from_G0(
    z: float, t_series: np.ndarray, G0_series: np.ndarray, kappa: float
) -> np.ndarray:
    """Compute *G(z,t)* from a known surface flux series *G0* (Eq. 9)."""
    if t_series.ndim != 1 or G0_series.ndim != 1:
        raise ValueError("t_series and G0_series must be 1-D")
    if t_series.size != G0_series.size:
        raise ValueError("t_series and G0_series must align")

    dt = np.diff(t_series)
    if not np.allclose(dt, dt[0]):
        raise ValueError("Time vector must be uniformly spaced")
    dt = dt[0]

    # Build F_z(t) = erfc(z / 2√(κ t))  (with F_z(0) = 0 by limit)
    with np.errstate(divide="ignore", invalid="ignore"):
        Fz = erfc(z / (2.0 * np.sqrt(kappa * t_series)))
    Fz[0] = 0.0

    dF = np.diff(Fz, prepend=0.0)
    # Convolution similar to temperature_convolution_solution
    Gz = np.cumsum(G0_series[::-1] * dF)[::-1]
    return Gz


def estimate_G0_from_Gz(
    Gz_series: np.ndarray, z_r: float, kappa: float, dt: float
) -> np.ndarray:
    """Estimate *surface* ground heat flux *G0* from plate measurements *Gz*.

    Implements discretised Eq. (11) – the recursion proposed by Wang & Bou-Zeid
    (2012).  Time-series must be *regularly* sampled.

    Parameters
    ----------
    Gz_series : np.ndarray
        Soil heat-flux measurements at depth *z_r* **(W m-2)**.
    z_r : float
        Plate depth **(m)**.
    kappa : float
        Thermal diffusivity **(m² s-1)**.
    dt : float
        Sampling interval **(s)**.

    Returns
    -------
    G0 : np.ndarray
        Estimated surface heat-flux series **(W m-2)**.
    """
    Gz_series = np.asarray(Gz_series, dtype=float)
    n_steps = Gz_series.size

    # Pre-compute ΔF_z(j) for j = 1 … n-1 (Eq. 10)
    j = np.arange(n_steps)  # 0 … n-1
    t_j = j * dt
    with np.errstate(divide="ignore", invalid="ignore"):
        Fz = erfc(z_r / (2.0 * np.sqrt(kappa * t_j)))
    Fz[0] = 0.0
    dF = np.diff(Fz, prepend=0.0)

    G0 = np.zeros_like(Gz_series)
    for n in range(1, n_steps):
        # J_{n-1} term (Eq. 12)
        J = 0.0
        for j in range(1, n):
            J += 0.5 * (G0[n - j] + G0[n - j - 1]) * dF[j]
        G0[n] = (2.0 * Gz_series[n] - J) / dF[1]
        # By construction dF[1] > 0 (t = dt)

    G0[0] = Gz_series[0]  # first guess – no history available
    return G0


# -----------------------------------------------------------------------------
# 4. Sinusoidal analytical solutions (Eqs. 13–15)
# -----------------------------------------------------------------------------


def sinusoidal_boundary_flux(
    t: float | np.ndarray, A: float, omega: float, epsilon: float
) -> float | np.ndarray:
    """Sinusoidal surface heat flux forcing (Eq. 13)."""
    return A * np.sin(omega * t + epsilon)


def soil_temperature_sinusoidal(
    z: float,
    t: float | np.ndarray,
    A: float,
    omega: float,
    epsilon: float,
    Ti: float,
    kappa: float,
) -> float | np.ndarray:
    """Analytical temperature under sinusoidal forcing (Eq. 14)."""
    r = np.sqrt(omega / (2.0 * kappa))
    exp_term = np.exp(-z * r)
    phase = omega * t + epsilon - z * r - math.pi / 4.0
    steady = A / (kappa * np.sqrt(omega)) * exp_term * np.sin(phase)

    # Transient integral (third term) – numeric quadrature (vectorised)
    def _integrand(x):
        return (
            (kappa * x**2 * math.sin(epsilon) - omega * math.cos(epsilon))
            * np.cos(x * z)
            / (omega**2 + (kappa**2) * x**4)
        )

    if np.isscalar(t):
        # Scalar: quad via np.trapz on a finite domain
        xi = np.linspace(0.0, 50.0 / z if z else 50.0, 2000)
        transient = (
            -2
            * A
            * kappa
            / math.pi
            * np.trapezoid(_integrand(xi) * np.exp(-kappa * xi**2 * t), xi)  # type: ignore
        )
    else:
        transient = np.zeros_like(t, dtype=float)
        xi = np.linspace(0.0, 50.0 / z if z else 50.0, 2000)
        integ = _integrand(xi)[:, None] * np.exp(-kappa * xi[:, None] ** 2 * t[None, :])
        transient = -2 * A * kappa / math.pi * np.trapezoid(integ, xi, axis=0)

    return Ti + steady + transient


def soil_heat_flux_sinusoidal(
    z: float,
    t: float | np.ndarray,
    A: float,
    omega: float,
    epsilon: float,
    kappa: float,
) -> float | np.ndarray:
    """Analytical heat-flux solution under sinusoidal forcing (Eq. 15)."""
    r = np.sqrt(omega / (2.0 * kappa))
    exp_term = np.exp(-z * r)
    phase = omega * t + epsilon - z * r
    steady = A * exp_term * np.sin(phase)

    # Transient integral similar to temperature – numeric quadrature
    def _integrand(x):
        return (
            (kappa * x**2 * math.sin(epsilon) - omega * math.cos(epsilon))
            * x
            * np.sin(x * z)
            / (omega**2 + (kappa**2) * x**4)
        )

    if np.isscalar(t):
        xi = np.linspace(0.0, 50.0 / z if z else 50.0, 2000)
        transient = (
            -2
            * A
            * kappa
            / math.pi
            * np.trapezoid(_integrand(xi) * np.exp(-kappa * xi**2 * t), xi)  # type: ignore
        )
    else:
        xi = np.linspace(0.0, 50.0 / z if z else 50.0, 2000)
        integ = _integrand(xi)[:, None] * np.exp(-kappa * xi[:, None] ** 2 * t[None, :])
        transient = -2 * A * kappa / math.pi * np.trapezoid(integ, xi, axis=0)

    return steady + transient


# -----------------------------------------------------------------------------
# 5. Soil-property parameterisations (Eqs. 16–18)
# -----------------------------------------------------------------------------


def heat_capacity_moist_soil(
    theta_v: float | np.ndarray,
    theta_s: float,
    rho_c_s: float = 1.26e6,
    rho_c_w: float = 4.20e6,
) -> float | np.ndarray:
    """Volumetric heat capacity of moist soil, Eq. (16).

    Parameters
    ----------
    theta_v
        Volumetric water content **(m³ m-3)**.
    theta_s
        Porosity (saturated volumetric water content) **(m³ m-3)**.
    rho_c_s, rho_c_w
        Heat capacity of dry soil / water **(J m-3 K-1)**.
    """
    return theta_v * rho_c_w + (1.0 - theta_s) * rho_c_s


def pf_from_theta(
    theta_v: float | np.ndarray, theta_s: float, psi_s: float, b: float
) -> float | np.ndarray:
    """Return Pf (Eq. 18) from volumetric water content."""
    return np.log10(100.0 * psi_s * (theta_s / theta_v) ** b)


def thermal_conductivity_moist_soil(
    theta_v: float | np.ndarray, theta_s: float, psi_s: float, b: float
) -> float | np.ndarray:
    """Thermal conductivity parameterisation, Eq. (17)."""
    Pf = pf_from_theta(theta_v, theta_s, psi_s, b)
    k = np.where(Pf <= 5.1, 0.420 * np.exp(-Pf - 2.7), 0.1744)
    return k


def thermal_diffusivity(
    k: float | np.ndarray, rho_c: float | np.ndarray
) -> float | np.ndarray:
    """Return κ = k / (ρ c)."""
    return k / rho_c
