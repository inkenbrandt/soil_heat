"""liebethal_and_folken.py
==================================
A collection of Python functions that implement every numbered
equation from Liebethal & Foken (2006) *Evaluation of six
parameterization approaches for the ground heat flux*.

Each public function is named after the paper section and
equation number for easy cross‑referencing.  Helper utilities
for finite‑difference gradients and unit handling are provided
at the end of the module.

References
----------
Liebethal, C., & Foken, T. (2006). Evaluation of six parameterization
approaches for the ground heat flux. *Theoretical and Applied Climatology*.
DOI:10.1007/s00704‑005‑0234‑0
"""

from __future__ import annotations

import numpy as np
from typing import Sequence, Tuple

# ---------------------------------------------------------------------
#  Utility finite‑difference helpers
# ---------------------------------------------------------------------


def _central_gradient(y: np.ndarray, x: np.ndarray) -> np.ndarray:
    """Central finite‑difference gradient with edge‑order 1.

    Parameters
    ----------
    y : ndarray
        Dependent variable samples.
    x : ndarray
        Independent variable samples (monotonic).

    Returns
    -------
    ndarray
        dy/dx evaluated at *x* using second‑order central differences
        and first‑order forward/backward differences at the edges.
    """
    y = np.asarray(y)
    x = np.asarray(x)
    if y.shape != x.shape:
        raise ValueError("y and x must have identical shape")
    return np.gradient(y, x, edge_order=1)


def _pad_nan_like(arr: np.ndarray) -> np.ndarray:
    """Return a NaN array with the same shape and dtype=float64."""
    return np.full_like(arr, np.nan, dtype=float)


# ---------------------------------------------------------------------
#  Eq. (1) – Reference ground‑heat‑flux (gradient + calorimetry)
# ---------------------------------------------------------------------


def reference_ground_heat_flux(
    temp_profile: np.ndarray,
    depths: Sequence[float],
    times: Sequence[float],
    cv: float,
    thermal_conductivity: float,
    gradient_depth: float = 0.20,
) -> np.ndarray:
    """
    Compute the reference ground‑heat flux *G₀,M* using the
    gradient method combined with calorimetry for heat storage
    (Liebethal & Foken 2006, Eq. 1).

    The method combines the conductive heat flux at a reference depth
    with the rate of change of heat stored in the soil layer above that
    depth.

    .. math::

        G_{0,M}(t) = -\\lambda \\frac{\\partial T}{\\partial z} \\bigg|_{z=0.2m}
                     + \\int_{z=0}^{0.2m} c_v \\frac{\\partial T}{\\partial t} dz

    Parameters
    ----------
    temp_profile : numpy.ndarray
        A 2D array of soil temperatures (°C or K) with shape
        `(n_depths, n_times)`.
    depths : Sequence[float]
        A sequence of measurement depths in meters (positive downward),
        corresponding to the rows of `temp_profile`.
    times : Sequence[float]
        A sequence of time stamps in seconds (e.g., Unix timestamps),
        corresponding to the columns of `temp_profile`.
    cv : float
        Volumetric heat capacity of the soil (J m⁻³ K⁻¹). Assumed
        to be constant with depth.
    thermal_conductivity : float
        Soil thermal conductivity λ (W m⁻¹ K⁻¹). Assumed to be
        constant with depth.
    gradient_depth : float, optional
        The depth (m) at which the vertical temperature gradient is
        evaluated, by default 0.20 m.

    Returns
    -------
    numpy.ndarray
        A 1D array of the instantaneous ground‑heat flux *G₀,M* (W m⁻²)
        at the surface for each time step. Positive values indicate
        downward flux.

    Raises
    ------
    ValueError
        If `temp_profile` shape does not match the lengths of `depths`
        and `times`.
    """
    depths = np.asarray(depths, dtype=float)
    times = np.asarray(times, dtype=float)
    T = np.asarray(temp_profile, dtype=float)

    if T.shape != (depths.size, times.size):
        raise ValueError("temp_profile shape must be (n_depths, n_times)")

    # ------------ gradient term
    # Calculate spatial gradient (∂T/∂z) for each time step
    dT_dz = np.array([_central_gradient(T[:, i], depths) for i in range(T.shape[1])]).T

    # Interpolate ∂T/∂z to the specified gradient_depth for each time step
    grad_T_at_z = np.array([np.interp(gradient_depth, depths, dT_dz[:, i]) for i in range(T.shape[1])])

    # ------------ storage (calorimetry) term
    # Calculate temporal gradient (∂T/∂t) for each depth
    dT_dt = np.array([_central_gradient(T[i, :], times) for i in range(T.shape[0])])

    # Integrate storage term over depth using the trapezoidal rule
    storage = cv * getattr(np, "trapezoid", getattr(np, "trapz"))(dT_dt, depths, axis=0)

    # Combine terms to get surface heat flux
    return -thermal_conductivity * grad_T_at_z + storage


# ---------------------------------------------------------------------
#  Eq. (2) – Percentage‑of‑net‑radiation parameterisation
# ---------------------------------------------------------------------


def ground_heat_flux_pr(qs: np.ndarray, p: float) -> np.ndarray:
    """
    Estimate ground heat flux as a fixed fraction of net radiation
    (Liebethal & Foken 2006, Eq. 2).

    This is a simple empirical relationship where the ground heat flux
    is assumed to be a constant proportion of the net radiation at the
    surface.

    .. math:: G_{0,PR}(t) = -p \\cdot Q^*_s(t)

    Parameters
    ----------
    qs : numpy.ndarray
        Time series of net radiation (W m⁻²). Positive values are
        typically downward.
    p : float
        The fraction of net radiation that is partitioned into
        ground heat flux (dimensionless, typically 0 < p < 1).

    Returns
    -------
    numpy.ndarray
        Time series of the estimated ground heat flux *G₀,PR* (W m⁻²).
    """
    return -p * np.asarray(qs, dtype=float)


# ---------------------------------------------------------------------
#  Eq. (3) – Linear regression against net radiation (with lag)
# ---------------------------------------------------------------------


def ground_heat_flux_lr(
    qs: np.ndarray, a: float, b: float, lag_steps: int = 0
) -> np.ndarray:
    """
    Estimate ground heat flux using a linear regression against net
    radiation, with an optional time lag (Liebethal & Foken 2006, Eq. 3).

    .. math:: G_{0,LR}(t) = a \\cdot Q^*_s(t + \\Delta t_G) + b

    Parameters
    ----------
    qs : numpy.ndarray
        Time series of net radiation (W m⁻²).
    a : float
        The slope of the linear regression (dimensionless).
    b : float
        The intercept of the linear regression (W m⁻²).
    lag_steps : int, optional
        The integer time lag (number of array elements) to apply to the
        net radiation series. A positive value advances the series
        (e.g., `qs[t+lag]` is used for `G[t]`), by default 0.

    Returns
    -------
    numpy.ndarray
        Time series of the estimated ground heat flux *G₀,LR* (W m⁻²).
    """
    qs = np.asarray(qs, dtype=float)
    if lag_steps != 0:
        # A negative roll shift corresponds to a positive time lag
        qs = np.roll(qs, -lag_steps)
    return a * qs + b


# ---------------------------------------------------------------------
#  Eq. (5–6) – Universal net‑radiation parameters A, B
# ---------------------------------------------------------------------


def ur_coefficients(delta_ts: float | np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """
    Compute the parameters A and B for the universal net-radiation
    parameterization, based on the diurnal surface temperature amplitude
    (Liebethal & Foken 2006, Eq. 5 & 6).

    .. math::
        A = 0.0074 \\cdot \\Delta T_s + 0.088
        B = 1729 \\cdot \\Delta T_s + 65013

    Parameters
    ----------
    delta_ts : float or numpy.ndarray
        The diurnal amplitude of the surface temperature (K or °C).

    Returns
    -------
    Tuple[numpy.ndarray, numpy.ndarray]
        A tuple containing the computed parameters (A, B).
        - A is dimensionless.
        - B is in seconds.
    """
    delta_ts = np.asarray(delta_ts, dtype=float)
    A = 0.0074 * delta_ts + 0.088
    B = 1729.0 * delta_ts + 65013.0
    return A, B


# ---------------------------------------------------------------------
#  Eq. (4) – Universal net‑radiation parameterisation
# ---------------------------------------------------------------------


def ground_heat_flux_ur(
    qs: np.ndarray, times_sec: np.ndarray, delta_ts: float
) -> np.ndarray:
    """
    Estimate ground heat flux using the universal net-radiation
    parameterization by Santanello & Friedl (2003), as cited in
    Liebethal & Foken (2006, Eq. 4).

    This method modulates the fraction of net radiation partitioned to
    ground heat flux with a cosine function of the time of day.

    .. math::
        G_{0,UR}(t) = -A \\cdot \\cos\\left(\\frac{2\\pi (t + 10800)}{B}\\right)
                      \\cdot Q^*_s(t)

    Parameters
    ----------
    qs : numpy.ndarray
        Time series of net radiation (W m⁻²).
    times_sec : numpy.ndarray
        Time stamps in **seconds relative to solar noon**. Positive values
        indicate the afternoon.
    delta_ts : float
        The diurnal amplitude of the surface temperature (K or °C).

    Returns
    -------
    numpy.ndarray
        Time series of the estimated ground heat flux *G₀,UR* (W m⁻²).
    """
    A, B = ur_coefficients(delta_ts)
    phase = np.cos(2 * np.pi * (times_sec + 10_800.0) / B)
    return -A * phase * np.asarray(qs, dtype=float)


# ---------------------------------------------------------------------
#  Eq. (7–8) – Surface‑temperature amplitude from two depths
# ---------------------------------------------------------------------


def surface_temp_amplitude(
    delta_t1: float, delta_t2: float, z1: float, z2: float
) -> float:
    """
    Estimate the diurnal surface-temperature amplitude (ΔT_s) from
    temperature amplitudes measured at two different soil depths
    (Liebethal & Foken 2006, Eq. 8).

    This method assumes an exponential decay of the temperature wave
    amplitude with depth.

    .. math::
        \\Delta T_s = \\Delta T_1 + \\Delta T_2 \\cdot
                      \\exp\\left(\\frac{z_2}{z_2 - z_1}\\right)

    Parameters
    ----------
    delta_t1 : float
        Diurnal temperature amplitude (K or °C) measured at depth `z1`.
    delta_t2 : float
        Diurnal temperature amplitude (K or °C) measured at depth `z2`.
    z1 : float
        The shallower depth in meters (positive downward).
    z2 : float
        The deeper depth in meters (positive downward).

    Returns
    -------
    float
        The estimated diurnal surface-temperature amplitude ΔT_s (K or °C).

    Raises
    ------
    ValueError
        If `z2` is not greater than `z1`.
    """
    if z2 <= z1:
        raise ValueError("Depth z2 must be greater than z1.")
    exponent = z2 / (z2 - z1)
    return delta_t1 + delta_t2 * np.exp(exponent)


# ---------------------------------------------------------------------
#  Eq. (9–10) – Sensible‑heat function parameterisation
# ---------------------------------------------------------------------


def phi_from_soil_moisture(
    theta_0_10: float, a_phi: float = 9.62, b_phi: float = 0.402
) -> float:
    """
    Calculate the empirical parameter φ based on soil moisture content
    (Liebethal & Foken 2006, Eq. 10).

    .. math:: \\phi = a_\\phi \\cdot \\theta_{0-10} + b_\\phi

    Parameters
    ----------
    theta_0_10 : float
        Average volumetric soil moisture content in the top 10 cm
        (m³ m⁻³).
    a_phi : float, optional
        Empirical coefficient, by default 9.62.
    b_phi : float, optional
        Empirical coefficient, by default 0.402.

    Returns
    -------
    float
        The dimensionless parameter φ.
    """
    return a_phi * theta_0_10 + b_phi


def ground_heat_flux_sh(
    h: np.ndarray,
    phase_g0: Sequence[float],
    phase_h: Sequence[float],
    u_mean: float,
    phi: float,
    omega: float = 2 * np.pi / 86_400.0,
) -> np.ndarray:
    """
    Estimate ground heat flux from the sensible heat flux (H)
    (Liebethal & Foken 2006, Eq. 9).

    This method relates the ground heat flux to the sensible heat flux
    through a phase-shifted and scaled relationship.

    .. math::
        G_{0,SH}(t) = -\\frac{\\phi}{\\sqrt{\\bar{u}}}
                      \\frac{\\cos(\\omega t + \\varphi(G_0))}
                           {\\cos(\\omega t + \\varphi(H))} H(t)

    Parameters
    ----------
    h : numpy.ndarray
        Time series of sensible heat flux (W m⁻²).
    phase_g0 : Sequence[float]
        Phase lags of the ground heat flux, φ(G₀), in **radians**. Must
        have the same length as `h`.
    phase_h : Sequence[float]
        Phase lags of the sensible heat flux, φ(H), in **radians**. Must
        have the same length as `h`.
    u_mean : float
        Mean horizontal wind speed during the daytime (m s⁻¹).
    phi : float
        An empirical dimensionless parameter, often derived from soil
        moisture via `phi_from_soil_moisture`.
    omega : float, optional
        The diurnal angular frequency (rad s⁻¹), by default `2π/86400`.

    Returns
    -------
    numpy.ndarray
        Time series of the estimated ground heat flux *G₀,SH* (W m⁻²).

    Raises
    ------
    ValueError
        If the length of phase arrays does not match the length of `h`.
    """
    h = np.asarray(h, dtype=float)
    if len(phase_g0) != len(h) or len(phase_h) != len(h):
        raise ValueError("Phase arrays must match the length of h.")

    t_steps = np.arange(len(h))
    cos_g0 = np.cos(omega * t_steps + np.asarray(phase_g0))
    cos_h = np.cos(omega * t_steps + np.asarray(phase_h))

    ratio = cos_g0 / cos_h
    return -(phi / np.sqrt(u_mean)) * ratio * h


# ---------------------------------------------------------------------
#  Eq. (11) – Simple‑measurement (heat‑flux plate) method
# ---------------------------------------------------------------------


def ground_heat_flux_sm(
    gp: np.ndarray,
    t1: np.ndarray,
    delta_t: np.ndarray,
    cv: float,
    zp: float,
    dt_seconds: float,
) -> np.ndarray:
    """
    Estimate surface ground heat flux using the "simple measurement"
    parameterization, correcting a heat flux plate measurement with a
    storage term (Liebethal & Foken 2006, Eq. 11).

    .. math::
        G_{0,SM}(t) = G_p(t) + c_v z_p \\left( \\frac{dT_1}{dt} +
                      \\frac{1}{2} \\frac{d(\\Delta T)}{dt} \\right)

    Parameters
    ----------
    gp : numpy.ndarray
        Heat flux plate measurements at depth `zp` (W m⁻²).
    t1 : numpy.ndarray
        Soil temperature at 0.01 m depth (K or °C).
    delta_t : numpy.ndarray
        Temperature difference T(0.01 m) - T(zp) (K).
    cv : float
        Volumetric heat capacity of the soil (J m⁻³ K⁻¹).
    zp : float
        Depth of the heat flux plate (m, positive downward).
    dt_seconds : float
        The constant time step between consecutive samples (s).

    Returns
    -------
    numpy.ndarray
        Time series of the estimated ground heat flux *G₀,SM* (W m⁻²).
        The first element will be NaN due to the backward difference.
    """
    gp = np.asarray(gp, dtype=float)
    t1 = np.asarray(t1, dtype=float)
    delta_t = np.asarray(delta_t, dtype=float)

    # Use central differences for better accuracy if possible, but paper
    # implies a time-stepping scheme. Using backward difference for dT/dt.
    dT1_dt = np.full_like(t1, np.nan)
    dDeltaT_dt = np.full_like(delta_t, np.nan)

    dT1_dt[1:] = (t1[1:] - t1[:-1]) / dt_seconds
    dDeltaT_dt[1:] = (delta_t[1:] - delta_t[:-1]) / dt_seconds

    storage_term = cv * zp * (dT1_dt + 0.5 * dDeltaT_dt)
    return gp + storage_term


# ---------------------------------------------------------------------
#  Eq. (12–13) – Force‑restore method
# ---------------------------------------------------------------------


def active_layer_thickness(
    lambda_: float, cv: float, omega: float = 2 * np.pi / 86_400
) -> float:
    """
    Calculate the thickness of the active soil layer (δz) for the
    force-restore method (Liebethal & Foken 2006, Eq. 13).

    .. math:: \\delta_z = \\sqrt{\\frac{\\lambda}{2 c_v \\omega}}

    Parameters
    ----------
    lambda_ : float
        Soil thermal conductivity (W m⁻¹ K⁻¹).
    cv : float
        Volumetric heat capacity of the soil (J m⁻³ K⁻¹).
    omega : float, optional
        The diurnal angular frequency (rad s⁻¹), by default `2π/86400`.

    Returns
    -------
    float
        The thickness of the active soil layer δz (m).
    """
    return np.sqrt(lambda_ / (2 * cv * omega))


def ground_heat_flux_fr(
    tg: np.ndarray,
    tg_avg: float,
    cv: float,
    lambda_: float,
    delta_z: float | None = None,
    times: np.ndarray | None = None,
) -> np.ndarray:
    """
    Estimate ground heat flux using the two-layer force-restore method
    (Liebethal & Foken 2006, Eq. 12).

    .. math::
        G_{0,FR}(t) = -\\delta_z c_v \\frac{dT_g}{dt} +
                      \\sqrt{\\lambda \\omega c_v} \\left(
                      \\frac{1}{\\omega}\\frac{dT_g}{dt} + (T_g - \\bar{T_g})
                      \\right)

    Parameters
    ----------
    tg : numpy.ndarray
        Time series of the upper (surface) layer temperature, Tg(t) (K).
    tg_avg : float
        The long-term average or "restoring" temperature, T̄g (K).
    cv : float
        Volumetric heat capacity of the soil (J m⁻³ K⁻¹).
    lambda_ : float
        Soil thermal conductivity λ (W m⁻¹ K⁻¹).
    delta_z : float, optional
        Thickness of the active soil layer δz (m). If `None`, it is
        calculated internally using `active_layer_thickness`,
        by default None.
    times : numpy.ndarray, optional
        Time stamps in seconds corresponding to `tg`. Required for
        calculating the time derivative. If `None`, assumes a uniform
        time step of 1 second, by default None.

    Returns
    -------
    numpy.ndarray
        Time series of the estimated ground heat flux *G₀,FR* (W m⁻²).
    """
    tg = np.asarray(tg, dtype=float)
    if times is None:
        times = np.arange(tg.size, dtype=float)
    else:
        times = np.asarray(times, dtype=float)

    dt_tg = _central_gradient(tg, times)

    if delta_z is None:
        delta_z = active_layer_thickness(lambda_, cv)

    omega = 2 * np.pi / 86_400.0  # diurnal frequency
    term1 = -delta_z * cv * dt_tg
    term2 = np.sqrt(lambda_ * omega * cv) * (dt_tg / omega + (tg - tg_avg))
    return term1 + term2


# ---------------------------------------------------------------------
#  Convenience enumerations of all public callables
# ---------------------------------------------------------------------

__all__ = [
    "reference_ground_heat_flux",
    "ground_heat_flux_pr",
    "ground_heat_flux_lr",
    "ur_coefficients",
    "ground_heat_flux_ur",
    "surface_temp_amplitude",
    "phi_from_soil_moisture",
    "ground_heat_flux_sh",
    "ground_heat_flux_sm",
    "active_layer_thickness",
    "ground_heat_flux_fr",
]
