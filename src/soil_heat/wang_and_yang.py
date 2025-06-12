"""
Yang & Wang (2008) – Equation Set Implementation
=================================================
Python translation of every numbered equation in:

> Yang, K., & Wang, J. (2008). *A temperature prediction‑correction method for
> estimating surface soil heat flux from soil temperature and moisture data.*
> *Science in China Series D: Earth Sciences, 51*(5), 721‑729.
> https://doi.org/10.1007/s11430‑008‑0036‑1

The paper introduces a *Temperature‑Diffusion plus Error‑Correction* (TDEC)
approach for estimating surface soil heat flux.  This module provides a direct
one‑to‑one mapping from each equation in the paper (Eqs. 1–12) to a Python
function.  Helper utilities required by those equations—matrix creators, grid
stretching, tridiagonal solvers, etc.—are also included.

All functions accept **NumPy arrays** or scalars where applicable and are fully
type‑annotated.  Docstrings use the **NumPy docstring standard** so they can be
rendered by *Sphinx‑napoleon*.
"""

from __future__ import annotations

import math
from typing import Sequence, Tuple

import numpy as np

# -----------------------------------------------------------------------------
# Constants & type aliases
# -----------------------------------------------------------------------------
SIGMA_SB: float = 5.670e-8  # Stefan‑Boltzmann constant (W m‑2 K‑4)
RHO_WATER: float = 1_000.0  # Density of liquid water (kg m‑3)
CP_WATER: float = 4_200_000.0  # Volumetric heat capacity of water (J m‑3 K‑1)

ArrayLike = np.ndarray | Sequence[float]

# -----------------------------------------------------------------------------
# Equation (1) & (2) –  the 1‑D thermal diffusion equation & Fourier’s law
# -----------------------------------------------------------------------------


def soil_heat_flux(
    Tz: ArrayLike, dz: ArrayLike, lambda_s: ArrayLike | float
) -> np.ndarray:  # Eq. (2)
    """Compute heat flux *G* at cell interfaces using Fourier’s law.

    Parameters
    ----------
    Tz : ArrayLike
        Temperature at layer **centres** (K).
    dz : ArrayLike
        Thickness of each layer (m).
    lambda_s : ArrayLike or float
        Thermal conductivity for each layer (W m‑1 K‑1).

    Returns
    -------
    np.ndarray
        Heat flux at *interfaces* (positive downward) with shape ``len(Tz)+1``.
    """
    Tz = np.asarray(Tz, dtype=float)
    dz = np.asarray(dz, dtype=float)
    lam = np.asarray(lambda_s, dtype=float) if np.ndim(lambda_s) else float(lambda_s)  # type: ignore

    # Conduction flux between layer i (upper) and i+1 (lower)
    dT = np.diff(Tz)
    if np.ndim(lam):
        lam_int = 0.5 * (lam[:-1] + lam[1:])  # type: ignore
    else:
        lam_int = lam
    G_int = lam_int * dT / dz[:-1]  # W m‑2, positive if T decreases with depth

    # Extrapolate surface & bottom fluxes (Neumann BC → same gradient)
    G_surface = lam_int[0] * (Tz[0] - Tz[1]) / dz[0]  # type: ignore
    G_bottom = lam_int[-1] * (Tz[-2] - Tz[-1]) / dz[-1]  # type: ignore
    return np.concatenate(([G_surface], G_int, [G_bottom]))


# -----------------------------------------------------------------------------
# Equation (3) & (5) –  integral form for *G(z)*
# -----------------------------------------------------------------------------


def integrated_soil_heat_flux(
    rho_c: ArrayLike,
    T_before: ArrayLike,
    T_after: ArrayLike,
    dz: ArrayLike,
    dt: float,
    G_ref: float = 0.0,
) -> np.ndarray:  # Eq. (5)
    """Discrete integration of Eq. (3)/(5) to obtain heat‑flux profile.

    Parameters
    ----------
    rho_c : ArrayLike
        Volumetric heat capacity ``ρ_s c_s`` for each layer (J m‑3 K‑1).
    T_before, T_after : ArrayLike
        Temperatures at two successive timesteps (K).
    dz : ArrayLike
        Layer thicknesses (m).
    dt : float
        Timestep (s).
    G_ref : float, default 0
        Heat flux at the lower reference depth *z_ref* (W m‑2).  Often ≈ 0.

    Returns
    -------
    np.ndarray
        Heat flux at the *upper* interface of every layer, size ``len(dz)``.
    """
    rho_c = np.asarray(rho_c)
    dT = np.asarray(T_after) - np.asarray(T_before)
    storage_term = np.cumsum(rho_c * dT * dz) / dt  # W m‑2
    return G_ref + storage_term


# -----------------------------------------------------------------------------
# Equation (4a–c) – volumetric heat capacity
# -----------------------------------------------------------------------------


def volumetric_heat_capacity(
    theta: ArrayLike, theta_sat: float | ArrayLike
) -> np.ndarray:  # Eq. 4
    """Volumetric heat capacity of moist soil.

    Parameters
    ----------
    theta : ArrayLike
        Volumetric water content (m³ m‑3).
    theta_sat : float or ArrayLike
        Soil porosity.

    Returns
    -------
    np.ndarray
        ``ρ_s c_s`` (J m‑3 K‑1).
    """
    theta = np.asarray(theta, dtype=float)
    theta_sat = np.asarray(theta_sat, dtype=float)

    # Eq. 4b & 4c
    rho_c_dry = (1.0 - theta_sat) * 2.1e6  # J m‑3 K‑1
    rho_c_water = CP_WATER  # 4.2 MJ m‑3 K‑1

    return rho_c_dry + rho_c_water * theta


# -----------------------------------------------------------------------------
# Equation (6a–b) – stretched vertical grid
# -----------------------------------------------------------------------------


def stretched_grid(n: int, D: float, xi: float) -> np.ndarray:  # Eq. 6
    """Generate *n* layer thicknesses following the exponential stretching rule.

    Parameters
    ----------
    n : int
        Number of layers.
    D : float
        Total domain depth (m).
    xi : float
        Stretching parameter; 0 → uniform grid.

    Returns
    -------
    np.ndarray
        Thickness ``Δz_i`` for each layer *i* (m).
    """
    if xi == 0:
        return np.full(n, D / n)
    delta_z0 = D * (math.exp(xi) - 1) / (math.exp(n * xi) - 1)
    dz = delta_z0 * np.exp(xi * np.arange(n))
    return dz


# -----------------------------------------------------------------------------
# Equation (7) –  implicit TDE discretisation (tridiagonal system)
# -----------------------------------------------------------------------------


def tridiagonal_coeffs(
    dz: ArrayLike,
    rho_c: ArrayLike,
    lambda_s: ArrayLike | float,
    dt: float,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Build *A*, *B*, *C* diagonals for the implicit TDE (Eq. 7b).

    Returns three 1‑D arrays representing the sub‑, main‑, and super‑diagonals.
    """
    dz = np.asarray(dz)
    lam = np.asarray(lambda_s) if np.ndim(lambda_s) else float(lambda_s)  # type: ignore
    n = len(dz)

    if np.ndim(lam):
        lam_up = np.concatenate(([lam[0]], 0.5 * (lam[:-1] + lam[1:])))  # type: ignore
        lam_dn = np.concatenate((0.5 * (lam[:-1] + lam[1:]), [lam[-1]]))  # type: ignore
    else:
        lam_up = lam_dn = np.full(n, lam)

    A = lam_up / dz / dz  # upper diagonal (i‑1)
    C = lam_dn / dz / dz  # lower diagonal (i+1)
    B = rho_c / dt + A + C  # type: ignore # main diagonal
    return A[1:], B, C[:-1]


def solve_tde(
    T_prev: ArrayLike,
    dz: ArrayLike,
    rho_c: ArrayLike,
    lambda_s: ArrayLike | float,
    Tsfc: float,
    Tbot: float,
    dt: float,
) -> np.ndarray:
    """Implicit Crank‑Nicholson (θ = 1) solve of Eq. (7).

    Boundary conditions (Eq. 7a, 7c) are Dirichlet.
    """
    from scipy.linalg import solve_banded  # lean import

    T_prev = np.asarray(T_prev)
    A, B, C = tridiagonal_coeffs(dz, rho_c, lambda_s, dt)
    n = len(B)

    # Assemble RHS (Eq. 7b – D term)
    D_vec = rho_c * T_prev / dt
    # Apply BCs
    D_vec[0] += A[0] * Tsfc
    D_vec[-1] += C[-1] * Tbot

    # Construct banded matrix for SciPy (3 × n)
    ab = np.zeros((3, n))
    ab[0, 1:] = C[:-1]  # super
    ab[1] = B  # main
    ab[2, :-1] = A[1:]  # sub

    T_new = solve_banded((1, 1), ab, D_vec)
    return np.concatenate(([Tsfc], T_new, [Tbot]))


# -----------------------------------------------------------------------------
# Temperature‑profile correction (Section 2.2)
# -----------------------------------------------------------------------------


def correct_profile(
    T_model: ArrayLike, depths_model: ArrayLike, T_obs: ArrayLike, depths_obs: ArrayLike
) -> np.ndarray:
    """Add linear‑interpolated bias (Eq. ΔT_k) to model profile."""
    bias = np.interp(
        depths_model, depths_obs, T_obs - np.interp(depths_obs, depths_model, T_model)
    )
    return T_model + bias


# -----------------------------------------------------------------------------
# Equation (8) –  surface temperature from long‑wave radiation
# -----------------------------------------------------------------------------


def surface_temperature_longwave(
    R_lw_up: float, R_lw_dn: float, emissivity: float = 0.98
) -> float:
    """Convert upward/downward long‑wave radiation to surface temperature (Eq. 8)."""
    return ((R_lw_up - (1.0 - emissivity) * R_lw_dn) / (emissivity * SIGMA_SB)) ** 0.25


# -----------------------------------------------------------------------------
# Equation (9a–c) –  thermal conductivity parameterisation
# -----------------------------------------------------------------------------


def thermal_conductivity_yang2008(
    theta: ArrayLike, theta_sat: float, rho_dry: float | ArrayLike
) -> np.ndarray:
    """Estimate soil thermal conductivity following Yang et al. (2005) (Eq. 9)."""
    theta = np.asarray(theta, dtype=float)
    rho_dry = np.asarray(rho_dry, dtype=float)

    lam_dry = (170.0 + 64.7 * rho_dry) / (2700.0 - 947.0 * rho_dry)  # Eq. 9b
    lam_sat = 2.0  # Eq. 9c (W m‑1 K‑1)
    lam = lam_dry + (lam_sat - lam_dry) * np.exp(
        0.36 * (theta / theta_sat - 1.0)
    )  # Eq. 9a
    return lam


# -----------------------------------------------------------------------------
# Equation (10–11) –  flux error for linear interpolation (diagnostic)
# -----------------------------------------------------------------------------


def flux_error_linear(
    rho_c: ArrayLike, S2_minus_S1: ArrayLike, dt: float
) -> np.ndarray:  # Eq. 11
    """Error introduced when using a LINEAR temperature profile (diagnostic)."""
    return rho_c * S2_minus_S1 / dt  # type: ignore


# -----------------------------------------------------------------------------
# Equation (12) –  surface energy budget closure
# -----------------------------------------------------------------------------


def surface_energy_residual(R_net: float, H: float, LE: float, G0: float) -> float:
    """Return the residual *ΔE* in Eq. (12)."""
    return R_net - (H + LE + G0)


# -----------------------------------------------------------------------------
# High‑level helper: one TDEC timestep
# -----------------------------------------------------------------------------


def tdec_step(
    T_prev: ArrayLike,
    dz: ArrayLike,
    theta: ArrayLike,
    theta_sat: float,
    rho_dry: float,
    lambda_const: float,
    Tsfc: float,
    Tbot: float,
    dt: float,
    depths_model: ArrayLike,
    T_obs: ArrayLike,
    depths_obs: ArrayLike,
) -> Tuple[np.ndarray, np.ndarray]:
    """One integration step of the TDEC scheme.

    Returns
    -------
    T_corr : np.ndarray
        Corrected temperature profile at *t + dt*.
    G_prof : np.ndarray
        Heat‑flux profile (W m‑2) at layer interfaces.
    """
    rho_c = volumetric_heat_capacity(theta, theta_sat)

    # 1. Predict with constant λ
    T_pred = solve_tde(
        T_prev=T_prev,
        dz=dz,
        rho_c=rho_c,
        lambda_s=lambda_const,
        Tsfc=Tsfc,
        Tbot=Tbot,
        dt=dt,
    )

    # 2. Correct with observed bias
    T_corr = correct_profile(T_pred, depths_model, T_obs, depths_obs)

    # 3. Compute heat‑flux profile (surface downward positive)
    G_prof = integrated_soil_heat_flux(
        rho_c=rho_c,
        T_before=T_prev,
        T_after=T_corr[1:-1],  # centre nodes
        dz=dz,
        dt=dt,
        G_ref=0.0,
    )
    return T_corr, G_prof


__all__ = [
    "soil_heat_flux",
    "integrated_soil_heat_flux",
    "volumetric_heat_capacity",
    "stretched_grid",
    "solve_tde",
    "correct_profile",
    "surface_temperature_longwave",
    "thermal_conductivity_yang2008",
    "flux_error_linear",
    "surface_energy_residual",
    "tdec_step",
]
