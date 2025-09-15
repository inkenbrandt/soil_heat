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
    """
    Compute heat flux *G* at cell interfaces using Fourier’s law.

    This function calculates the conductive heat flux between soil layers
    based on the temperature gradient and thermal conductivity.

    .. math:: G = -\\lambda_s \\frac{\\partial T}{\\partial z}

    Parameters
    ----------
    Tz : ArrayLike
        A 1D array of temperatures at the **centre** of each soil layer (K).
    dz : ArrayLike
        A 1D array of the thickness of each soil layer (m).
    lambda_s : ArrayLike or float
        Thermal conductivity for each layer (W m⁻¹ K⁻¹). Can be a single
        value (for homogeneous soil) or an array matching `Tz`.

    Returns
    -------
    numpy.ndarray
        An array of heat fluxes (W m⁻²) at the *interfaces* between layers,
        including the surface and bottom boundaries. The length of the
        output is `len(Tz) + 1`. Positive values indicate downward flux.
    """
    Tz = np.asarray(Tz, dtype=float)
    dz = np.asarray(dz, dtype=float)
    lam = np.asarray(lambda_s, dtype=float) if np.ndim(lambda_s) else float(lambda_s)

    # Effective conductivity at interfaces is the average of adjacent layers
    if np.ndim(lam):
        lam_int = 0.5 * (lam[:-1] + lam[1:])
    else:
        lam_int = lam

    # Gradient between layer centers
    dT = np.diff(Tz)
    # Distance between layer centers
    dz_int = 0.5 * (dz[:-1] + dz[1:])

    G_int = -lam_int * dT / dz_int

    # Extrapolate surface & bottom fluxes assuming the same gradient as the
    # nearest two layers (Neumann boundary condition approximation).
    G_surface = -lam[0] * (Tz[0] - Tz[1]) / (dz[0]/2 + dz[1]/2)
    G_bottom = -lam[-1] * (Tz[-2] - Tz[-1]) / (dz[-2]/2 + dz[-1]/2)

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
    """
    Calculate the soil heat flux profile by integrating the change in
    heat storage upwards from a reference depth (Yang & Wang 2008, Eq. 5).

    .. math::
        G(z_i) = G(z_{ref}) + \\int_{z_i}^{z_{ref}} \\rho_s c_s
                 \\frac{\\partial T}{\\partial t} dz

    Parameters
    ----------
    rho_c : ArrayLike
        A 1D array of volumetric heat capacity (`ρ_s c_s`) for each soil
        layer (J m⁻³ K⁻¹).
    T_before : ArrayLike
        A 1D array of temperatures at the start of the time step (K).
    T_after : ArrayLike
        A 1D array of temperatures at the end of the time step (K).
    dz : ArrayLike
        A 1D array of layer thicknesses (m).
    dt : float
        The duration of the time step (s).
    G_ref : float, optional
        The heat flux at the lower reference depth `z_ref` (W m⁻²). This
        is typically assumed to be zero at a sufficient depth, by default 0.0.

    Returns
    -------
    numpy.ndarray
        An array of heat fluxes (W m⁻²) at the *upper interface* of each
        layer, with a size of `len(dz)`.
    """
    rho_c = np.asarray(rho_c)
    dT = np.asarray(T_after) - np.asarray(T_before)

    # Storage change in each layer
    storage_change = rho_c * dT * dz / dt

    # Cumulatively sum storage changes from the bottom up
    # The flux at interface i is the reference flux plus the sum of storage
    # changes in all layers below i.
    return G_ref + np.cumsum(storage_change[::-1])[::-1]


# -----------------------------------------------------------------------------
# Equation (4a–c) – volumetric heat capacity
# -----------------------------------------------------------------------------


def volumetric_heat_capacity(
    theta: ArrayLike, theta_sat: float | ArrayLike
) -> np.ndarray:  # Eq. 4
    """
    Calculate the volumetric heat capacity of moist soil based on its
    water content (Yang & Wang 2008, Eq. 4).

    .. math::
        \\rho_s c_s = (1 - \\theta_{sat}) \\rho_{d}c_{d} + \\theta \\rho_w c_w

    Parameters
    ----------
    theta : ArrayLike
        Volumetric water content (m³ m⁻³).
    theta_sat : float or ArrayLike
        Soil porosity, i.e., saturated volumetric water content (m³ m⁻³).

    Returns
    -------
    numpy.ndarray
        The volumetric heat capacity `ρ_s c_s` (J m⁻³ K⁻¹).
    """
    theta = np.asarray(theta, dtype=float)
    theta_sat = np.asarray(theta_sat, dtype=float)

    # From paper, typical values for dry soil and water
    rho_c_dry = (1.0 - theta_sat) * 2.1e6
    rho_c_water = CP_WATER

    return rho_c_dry + rho_c_water * theta


# -----------------------------------------------------------------------------
# Equation (6a–b) – stretched vertical grid
# -----------------------------------------------------------------------------


def stretched_grid(n: int, D: float, xi: float) -> np.ndarray:  # Eq. 6
    """
    Generate layer thicknesses for a vertically stretched grid.

    The grid layer thickness increases exponentially with depth, allowing
    for higher resolution near the surface.

    .. math:: \\Delta z_i = \\Delta z_0 \\exp(\\xi (i-1))

    Parameters
    ----------
    n : int
        The number of soil layers.
    D : float
        The total depth of the soil domain (m).
    xi : float
        The stretching parameter. `xi = 0` results in a uniform grid.
        `xi > 0` results in a grid that is finer at the top.

    Returns
    -------
    numpy.ndarray
        A 1D array of thickness `Δz_i` for each of the `n` layers (m).
    """
    if xi == 0:
        return np.full(n, D / n)
    # First layer thickness, derived from the sum of geometric series
    delta_z0 = D * (math.exp(xi) - 1) / (math.exp(n * xi) - 1)
    # Thickness of each subsequent layer
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
    """
    Construct the coefficient diagonals (A, B, C) for the tridiagonal
    matrix system used in the implicit TDE solver (Yang & Wang 2008, Eq. 7b).

    These coefficients represent the discretized heat diffusion equation.

    Parameters
    ----------
    dz : ArrayLike
        A 1D array of layer thicknesses (m).
    rho_c : ArrayLike
        A 1D array of volumetric heat capacity for each layer (J m⁻³ K⁻¹).
    lambda_s : ArrayLike or float
        Thermal conductivity for each layer (W m⁻¹ K⁻¹).
    dt : float
        The time step (s).

    Returns
    -------
    Tuple[numpy.ndarray, numpy.ndarray, numpy.ndarray]
        A tuple containing the sub-diagonal (A), main-diagonal (B), and
        super-diagonal (C) of the tridiagonal matrix.
    """
    dz = np.asarray(dz)
    rho_c = np.asarray(rho_c)
    lam = np.asarray(lambda_s) if np.ndim(lambda_s) else float(lambda_s)
    n = len(dz)

    # Effective conductivity at interfaces
    if np.ndim(lam):
        lam_up = 0.5 * (lam[:-1] + lam[1:])
        lam_dn = lam_up
    else:
        lam_up = lam_dn = np.full(n - 1, lam)

    # Coefficients from finite difference scheme
    alpha = lam_up / (0.5 * (dz[:-1] + dz[1:]))
    beta = lam_dn / (0.5 * (dz[1:] + dz[:-1]))

    A = -alpha / dz[:-1]
    C = -beta / dz[1:]
    B = rho_c[1:-1]/dt - A - C

    # Adjust for boundary conditions
    B_full = np.zeros(n)
    B_full[0] = rho_c[0]/dt + alpha[0]
    B_full[-1] = rho_c[-1]/dt + beta[-1]
    B_full[1:-1] = B

    return A, B_full, C


def solve_tde(
    T_prev: ArrayLike,
    dz: ArrayLike,
    rho_c: ArrayLike,
    lambda_s: ArrayLike | float,
    Tsfc: float,
    Tbot: float,
    dt: float,
) -> np.ndarray:
    """
    Solve the 1D thermal diffusion equation for one time step using an
    implicit Crank-Nicolson scheme (Yang & Wang 2008, Eq. 7).

    This function sets up and solves the tridiagonal system of linear
    equations `M * T_new = D` to find the temperature profile at the
    next time step.

    Parameters
    ----------
    T_prev : ArrayLike
        Temperature profile at the previous time step (K).
    dz : ArrayLike
        Layer thicknesses (m).
    rho_c : ArrayLike
        Volumetric heat capacity of each layer (J m⁻³ K⁻¹).
    lambda_s : ArrayLike or float
        Thermal conductivity of each layer (W m⁻¹ K⁻¹).
    Tsfc : float
        Surface temperature boundary condition (K).
    Tbot : float
        Bottom temperature boundary condition (K).
    dt : float
        Time step (s).

    Returns
    -------
    numpy.ndarray
        The new temperature profile at `t + dt`, including boundary nodes.
    """
    from scipy.linalg import solve_banded

    T_prev = np.asarray(T_prev)
    A, B, C = tridiagonal_coeffs(dz, rho_c, lambda_s, dt)
    n = len(B)

    # Assemble RHS vector D from Eq. 7b
    D_vec = np.asarray(rho_c) * T_prev / dt
    # Apply Dirichlet boundary conditions
    D_vec[0] += A[0] * Tsfc
    D_vec[-1] += C[-1] * Tbot

    # Create the banded matrix for SciPy's solver
    # The matrix has shape (3, n) for a tridiagonal system
    ab = np.zeros((3, n))
    ab[0, 1:] = -C  # Super-diagonal
    ab[1, :] = B    # Main-diagonal
    ab[2, :-1] = -A  # Sub-diagonal

    T_new_internal = solve_banded((1, 1), ab, D_vec)
    return np.concatenate(([Tsfc], T_new_internal, [Tbot]))


# -----------------------------------------------------------------------------
# Temperature‑profile correction (Section 2.2)
# -----------------------------------------------------------------------------


def correct_profile(
    T_model: ArrayLike, depths_model: ArrayLike, T_obs: ArrayLike, depths_obs: ArrayLike
) -> np.ndarray:
    """
    Correct a modeled temperature profile using observed temperatures.

    This function calculates the bias (error) between the model and
    observations at the observation depths, then linearly interpolates
    this bias across the entire model grid to correct the profile.

    Parameters
    ----------
    T_model : ArrayLike
        The modeled temperature profile (K).
    depths_model : ArrayLike
        The depths corresponding to `T_model` (m).
    T_obs : ArrayLike
        The observed temperatures (K).
    depths_obs : ArrayLike
        The depths of the observations (m).

    Returns
    -------
    numpy.ndarray
        The corrected temperature profile.
    """
    T_model = np.asarray(T_model)
    # First, interpolate the model temperature to the observation depths
    T_model_at_obs_depths = np.interp(depths_obs, depths_model, T_model)
    # Calculate the bias at observation depths
    bias_at_obs_depths = T_obs - T_model_at_obs_depths
    # Interpolate the bias to all model depths
    bias_on_model_grid = np.interp(depths_model, depths_obs, bias_at_obs_depths)

    return T_model + bias_on_model_grid


# -----------------------------------------------------------------------------
# Equation (8) –  surface temperature from long‑wave radiation
# -----------------------------------------------------------------------------


def surface_temperature_longwave(
    R_lw_up: float, R_lw_dn: float, emissivity: float = 0.98
) -> float:
    """
    Calculate surface temperature from upward and downward long-wave
    radiation measurements using the Stefan-Boltzmann law (Yang & Wang 2008, Eq. 8).

    .. math::
        T_s = \\left[ \\frac{R_{lw}^{\\uparrow} - (1 - \\epsilon) R_{lw}^{\\downarrow}}
                      {\\epsilon \\sigma} \\right]^{1/4}

    Parameters
    ----------
    R_lw_up : float
        Upwelling long-wave radiation (W m⁻²).
    R_lw_dn : float
        Downwelling long-wave radiation (W m⁻²).
    emissivity : float, optional
        Surface emissivity (dimensionless), by default 0.98.

    Returns
    -------
    float
        The calculated surface temperature (K).
    """
    numerator = R_lw_up - (1.0 - emissivity) * R_lw_dn
    denominator = emissivity * SIGMA_SB
    return (numerator / denominator) ** 0.25


# -----------------------------------------------------------------------------
# Equation (9a–c) –  thermal conductivity parameterisation
# -----------------------------------------------------------------------------


def thermal_conductivity_yang2008(
    theta: ArrayLike, theta_sat: float, rho_dry: float | ArrayLike
) -> np.ndarray:
    """
    Estimate soil thermal conductivity based on soil moisture and dry
    density, following Yang et al. (2005) as cited in
    Yang & Wang (2008, Eq. 9).

    Parameters
    ----------
    theta : ArrayLike
        Volumetric water content (m³ m⁻³).
    theta_sat : float
        Saturated volumetric water content (porosity) (m³ m⁻³).
    rho_dry : float or ArrayLike
        Dry soil bulk density (kg m⁻³).

    Returns
    -------
    numpy.ndarray
        The estimated soil thermal conductivity `λ_s` (W m⁻¹ K⁻¹).
    """
    theta = np.asarray(theta, dtype=float)
    rho_dry = np.asarray(rho_dry, dtype=float) / 1000 # Convert to g cm-3 for formula

    # Eq. 9b for dry thermal conductivity
    lam_dry = (0.170 + 0.0647 * rho_dry) / (2.7 - 0.947 * rho_dry) - 0.2
    # Eq. 9c for saturated thermal conductivity
    lam_sat = 2.0

    # Eq. 9a mixing model
    lam = lam_dry + (lam_sat - lam_dry) * np.exp(
        0.36 * (theta / theta_sat - 1.0)
    )
    return lam


# -----------------------------------------------------------------------------
# Equation (10–11) –  flux error for linear interpolation (diagnostic)
# -----------------------------------------------------------------------------


def flux_error_linear(
    rho_c: ArrayLike, S2_minus_S1: ArrayLike, dt: float
) -> np.ndarray:  # Eq. 11
    """
    Calculate the diagnostic error in heat flux that arises from assuming
    a linear temperature profile between measurement points
    (Yang & Wang 2008, Eq. 11).

    .. math:: \\Delta G_i = \\frac{\\rho_s c_s (S_{i,2} - S_{i,1})}{\\Delta t}

    Parameters
    ----------
    rho_c : ArrayLike
        Volumetric heat capacity for each layer (J m⁻³ K⁻¹).
    S2_minus_S1 : ArrayLike
        The area difference representing the deviation from linearity.
    dt : float
        The time step (s).

    Returns
    -------
    numpy.ndarray
        The flux error for each layer (W m⁻²).
    """
    return np.asarray(rho_c) * np.asarray(S2_minus_S1) / dt


# -----------------------------------------------------------------------------
# Equation (12) –  surface energy budget closure
# -----------------------------------------------------------------------------


def surface_energy_residual(R_net: float, H: float, LE: float, G0: float) -> float:
    """
    Calculate the residual of the surface energy budget (Yang & Wang 2008, Eq. 12).

    .. math:: \\Delta E = R_{net} - (H + LE + G_0)

    Parameters
    ----------
    R_net : float
        Net radiation (W m⁻²).
    H : float
        Sensible heat flux (W m⁻²).
    LE : float
        Latent heat flux (W m⁻²).
    G0 : float
        Surface ground heat flux (W m⁻²).

    Returns
    -------
    float
        The energy balance residual `ΔE` (W m⁻²).
    """
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
    """
    Perform a single integration step of the TDEC (Temperature Diffusion
    Error Correction) scheme.

    This function encapsulates the predict-correct sequence for one time step.

    Parameters
    ----------
    T_prev : ArrayLike
        Temperature profile from the previous time step (K).
    dz : ArrayLike
        Layer thicknesses (m).
    theta : ArrayLike
        Volumetric water content profile (m³ m⁻³).
    theta_sat : float
        Soil porosity (m³ m⁻³).
    rho_dry : float
        Dry soil bulk density (kg m⁻³).
    lambda_const : float
        A constant thermal conductivity for the prediction step (W m⁻¹ K⁻¹).
    Tsfc : float
        Surface temperature boundary condition (K).
    Tbot : float
        Bottom temperature boundary condition (K).
    dt : float
        Time step (s).
    depths_model : ArrayLike
        Depths of the model grid nodes (m).
    T_obs : ArrayLike
        Observed temperatures for the correction step (K).
    depths_obs : ArrayLike
        Depths of the observed temperatures (m).

    Returns
    -------
    Tuple[np.ndarray, np.ndarray]
        A tuple containing:
        - `T_corr`: The corrected temperature profile at `t + dt`.
        - `G_prof`: The corresponding heat flux profile (W m⁻²) at layer interfaces.
    """
    rho_c = volumetric_heat_capacity(theta, theta_sat)

    # 1. Predict new temperature profile using the TDE solver
    T_pred = solve_tde(
        T_prev=T_prev,
        dz=dz,
        rho_c=rho_c,
        lambda_s=lambda_const,
        Tsfc=Tsfc,
        Tbot=Tbot,
        dt=dt,
    )

    # 2. Correct the predicted profile using observational data
    T_corr = correct_profile(T_pred, depths_model, T_obs, depths_obs)

    # 3. Compute the resulting heat flux profile from the corrected temperatures
    G_prof = integrated_soil_heat_flux(
        rho_c=rho_c,
        T_before=T_prev,
        T_after=T_corr[1:-1],  # Use internal nodes of corrected profile
        dz=dz,
        dt=dt,
        G_ref=0.0,  # Assume zero flux at the bottom
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
