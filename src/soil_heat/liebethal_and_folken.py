"""liebethal_and_folken.py
==================================
A collection of Python functions that implement every numbered
equation from Liebethal & Foken (2006) *Evaluation of six
parameterization approaches for the ground heat flux*.

Each public function is named after the paper section and
equation number for easy cross‑referencing.  Helper utilities
for finite‐difference gradients and unit handling are provided
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
    """Central finite‑difference gradient with edge‐order 1.

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
    """Compute the reference ground‑heat flux *G₀,M* (Eq. 1).

    Equation
    --------
    G₀,M(t) = -λ ∂T/∂z |_(z=0.2 m) + ∫_{z=0}^{0.2 m} c_v ∂T/∂t dz

    Parameters
    ----------
    temp_profile : ndarray, shape (n_z, n_t)
        Soil temperatures (°C or K) at the depths specified by *depths* and
        time stamps *times*.
    depths : sequence of float, length *n_z*
        Measurement depths (m, **positive downward**).
    times : sequence of float, length *n_t*
        Epoch time in **seconds** (may be monotonic pandas DatetimeIndex
        converted via ``astype('int64')/1e9``).
    cv : float
        Volumetric heat capacity of the soil (J m⁻³ K⁻¹).
    thermal_conductivity : float
        Soil thermal conductivity λ (W m⁻¹ K⁻¹).
    gradient_depth : float, default 0.20
        Depth (m) at which the vertical gradient term is evaluated.

    Returns
    -------
    ndarray, shape (n_t,)
        Instantaneous ground‑heat flux *G₀,M* (W m⁻²). Positive = downward.
    """
    depths = np.asarray(depths, dtype=float)  # type: ignore
    times = np.asarray(times, dtype=float)  # type: ignore
    T = np.asarray(temp_profile, dtype=float)

    if T.shape != (depths.size, times.size):  # type: ignore
        raise ValueError("temp_profile shape must be (n_depths, n_times)")

    # ------------ gradient term
    dT_dz = _central_gradient(T, depths[:, None])  # type: ignore  shape (n_z, n_t)

    # Interpolate ∂T/∂z to the gradient_depth
    grad_T_at_z = np.interp(
        gradient_depth,
        depths,
        dT_dz,
        left=np.nan,
        right=np.nan,
    )

    # ------------ storage (calorimetry) term
    dT_dt = _central_gradient(T, times[None, :])  # type: ignore shape (n_z, n_t)

    # Integrate over depth using trapezoidal rule (axis=0 is depth)
    storage = cv * np.trapezoid(dT_dt, depths, axis=0)

    return -thermal_conductivity * grad_T_at_z + storage  # type: ignore


# ---------------------------------------------------------------------
#  Eq. (2) – Percentage‑of‑net‑radiation parameterisation
# ---------------------------------------------------------------------


def ground_heat_flux_pr(qs: np.ndarray, p: float) -> np.ndarray:
    """Ground heat flux using a fixed *p* fraction of net radiation (Eq. 2).

    G₀,PR(t) = ‑p · Q*ₛ(t)

    Parameters
    ----------
    qs : ndarray
        Net radiation time series (W m⁻²). Positive = downward.
    p : float
        Fraction of net radiation that becomes ground‑heat flux (0–1).

    Returns
    -------
    ndarray
        G₀,PR (W m⁻²).
    """
    return -p * np.asarray(qs, dtype=float)


# ---------------------------------------------------------------------
#  Eq. (3) – Linear regression against net radiation (with lag)
# ---------------------------------------------------------------------


def ground_heat_flux_lr(
    qs: np.ndarray, a: float, b: float, lag_steps: int = 0
) -> np.ndarray:
    """Linear net‑radiation parameterisation (Eq. 3).

    G₀,LR(t) = a·Q*ₛ(t+Δt_G) + b

    Parameters
    ----------
    qs : ndarray
        Net radiation (W m⁻²).
    a, b : float
        Regression coefficients.
    lag_steps : int, default 0
        Integer lag (number of samples) by which *qs* is advanced.

    Returns
    -------
    ndarray
        G₀,LR (W m⁻²).
    """
    qs = np.asarray(qs, dtype=float)
    if lag_steps != 0:
        qs = np.roll(qs, -lag_steps)
    return a * qs + b


# ---------------------------------------------------------------------
#  Eq. (5–6) – Universal net‑radiation parameters A, B
# ---------------------------------------------------------------------


def ur_coefficients(delta_ts: float | np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """Compute universal‑function parameters *A* and *B* (Eq. 5 & 6).

    Parameters
    ----------
    delta_ts : float or ndarray
        Diurnal amplitude of *surface* temperature (K).

    Returns
    -------
    A : ndarray
    B : ndarray (seconds)
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
    """Universal net‑radiation parameterisation (Eq. 4).

    Implements Santanello & Friedl (2003):
        G₀,UR(t) = -A · cos[2π (t + 10800) / B] · Q*ₛ(t)

    *t* is **seconds since solar noon** (positive in afternoon).

    Parameters
    ----------
    qs : ndarray
        Net radiation (W m⁻²).
    times_sec : ndarray
        Seconds relative to solar noon (s).
    delta_ts : float
        Diurnal surface‑temperature amplitude (K).

    Returns
    -------
    ndarray
        G₀,UR (W m⁻²).
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
    """Compute diurnal surface‑temperature amplitude ΔT_s (Eq. 8).

    Parameters
    ----------
    delta_t1, delta_t2 : float
        Diurnal temperature amplitudes (K) measured at depths *z1* and *z2*.
    z1, z2 : float
        Depths in meters (**positive downward**, with z2 > z1 > 0).

    Returns
    -------
    float
        Estimated ΔT_s (K).
    """
    if z2 <= z1:
        raise ValueError("Require z2 > z1")
    exponent = z2 / (z2 - z1)
    return delta_t1 + delta_t2 * np.exp(exponent)


# ---------------------------------------------------------------------
#  Eq. (9–10) – Sensible‑heat function parameterisation
# ---------------------------------------------------------------------


def phi_from_soil_moisture(
    theta_0_10: float, a_phi: float = 9.62, b_phi: float = 0.402
) -> float:
    """Soil‑moisture dependent φ (Eq. 10)."""
    return a_phi * theta_0_10 + b_phi


def ground_heat_flux_sh(
    h: np.ndarray,
    phase_g0: Sequence[float],
    phase_h: Sequence[float],
    u_mean: float,
    phi: float,
    omega: float = 2 * np.pi / 86_400.0,
) -> np.ndarray:
    """Ground‑heat flux from sensible heat flux H (Eq. 9).

    Parameters
    ----------
    h : ndarray
        Sensible heat flux time series (W m⁻²).
    phase_g0, phase_h : sequence of float
        Phase lags φ(G₀) and φ(H) in **radians**.
    u_mean : float
        Mean horizontal wind speed during daytime (m s⁻¹).
    phi : float
        Empirical parameter (dimensionless), see `phi_from_soil_moisture`.
    omega : float, default 2π/86400
        Diurnal angular frequency (s⁻¹).

    Returns
    -------
    ndarray
        G₀,SH (W m⁻²).
    """
    if len(phase_g0) != len(h) or len(phase_h) != len(h):
        raise ValueError("Phase arrays must match length of h")
    ratio = np.cos(omega * np.arange(len(h)) + phase_g0) / np.cos(
        omega * np.arange(len(h)) + phase_h
    )
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
    """Simple‑measurement parameterisation (Eq. 11).

    Parameters
    ----------
    gp : ndarray
        Heat‑flux plate measurement at depth *zp* (W m⁻²).
    t1 : ndarray
        Soil temperature at 0.01 m depth (K or °C).
    delta_t : ndarray
        Temperature difference T(0.01 m) – T(z_p) (K).
    cv : float
        Volumetric heat capacity (J m⁻³ K⁻¹).
    zp : float
        Plate depth (m, positive downward).
    dt_seconds : float
        Time step between consecutive samples (s).

    Returns
    -------
    ndarray
        G₀,SM (W m⁻²).
    """
    gp = np.asarray(gp, dtype=float)
    t1 = np.asarray(t1, dtype=float)
    delta_t = np.asarray(delta_t, dtype=float)

    # Finite differences using backward stencil to match t‑Δt in paper.
    dT1_dt = np.empty_like(t1)
    dDeltaT_dt = np.empty_like(delta_t)
    dT1_dt[1:] = (t1[1:] - t1[:-1]) / dt_seconds
    dT1_dt[0] = np.nan
    dDeltaT_dt[1:] = (delta_t[1:] - delta_t[:-1]) / dt_seconds
    dDeltaT_dt[0] = np.nan

    storage_term = cv * zp * (dT1_dt + 0.5 * dDeltaT_dt)
    return gp + storage_term


# ---------------------------------------------------------------------
#  Eq. (12–13) – Force‑restore method
# ---------------------------------------------------------------------


def active_layer_thickness(
    lambda_: float, cv: float, omega: float = 2 * np.pi / 86_400
) -> float:
    """Thickness δz of the active soil layer (Eq. 13)."""
    return np.sqrt(lambda_ / (2 * cv * omega))


def ground_heat_flux_fr(
    tg: np.ndarray,
    tg_avg: float,
    cv: float,
    lambda_: float,
    delta_z: float | None = None,
    times: np.ndarray | None = None,
) -> np.ndarray:
    """Force‑restore ground‑heat flux (Eq. 12).

    Implements the two‑layer force‑restore formulation with an optional
    diagnostic *δz* computed via Eq. 13 if not supplied.

    Parameters
    ----------
    tg : ndarray
        Temperature of the upper (surface) layer Tg(t).
    tg_avg : float
        Long‑term average or restoring temperature Tḡ (K).
    cv : float
        Volumetric heat capacity (J m⁻³ K⁻¹).
    lambda_ : float
        Soil thermal conductivity λ (W m⁻¹ K⁻¹).
    delta_z : float, optional
        Thickness of the active soil layer δz (m).  If *None*, computed
        from ``active_layer_thickness``.
    times : ndarray, optional
        Time stamps in seconds.  Required when *delta_z* is None or when
        irregular sampling; defaults to `np.arange(len(tg))` seconds.

    Returns
    -------
    ndarray
        G₀,FR (W m⁻²).
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
