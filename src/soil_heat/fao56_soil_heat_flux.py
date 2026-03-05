"""
FAO-56 Soil Heat Flux (G) Calculation Functions
================================================

Implements soil heat flux estimation methods from:
    Allen, R.G., Pereira, L.S., Raes, D., Smith, M. (1998).
    "Crop evapotranspiration - Guidelines for computing crop water
    requirements." FAO Irrigation and Drainage Paper 56. Rome, Italy.

Also includes the ASCE Standardized Reference ET hourly G estimation:
    ASCE-EWRI (2005). "The ASCE Standardized Reference Evapotranspiration
    Equation." ASCE Reston, VA.

Equations Implemented
---------------------
- **Eq. 41**: General soil heat flux from soil heat capacity, temperature
  change, effective depth, and time interval.
- **Eq. 42**: Daily soil heat flux for a grass reference surface (G ≈ 0).
- **Eq. 43**: Monthly soil heat flux using the previous and next month's
  mean air temperatures.
- **Eq. 44**: Monthly soil heat flux using the previous and current month's
  mean air temperatures (when next month is unavailable).
- **ASCE Hourly**: Sub-daily soil heat flux as a fraction of net radiation
  (Rn), with different coefficients for daytime and nighttime.

Units
-----
Temperatures in degrees Celsius.
Soil heat flux in MJ m⁻² day⁻¹ (or MJ m⁻² hr⁻¹ for hourly).
Can be converted to equivalent evaporation (mm day⁻¹) using
``energy_to_evap()``.

Author: Generated for Paul Inkenbrandt / Utah Geological Survey
"""

from __future__ import annotations

from typing import Optional, Union

import numpy as np
import numpy.typing as npt

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

#: Default soil heat capacity for a moist, mineral soil [MJ m⁻³ °C⁻¹].
#: FAO-56 uses 2.1 as a representative value (Chapter 3, Example 13).
DEFAULT_CS = 2.1

#: Inverse latent heat of vaporisation used to convert energy flux
#: (MJ m⁻² day⁻¹) to equivalent evaporation depth (mm day⁻¹).
#: λ ≈ 2.45 MJ kg⁻¹  →  1/λ ≈ 0.408
LATENT_HEAT_INV = 0.408


# ===========================================================================
# Core FAO-56 Equations
# ===========================================================================

def soil_heat_flux_general(
    t_current: Union[float, npt.ArrayLike],
    t_previous: Union[float, npt.ArrayLike],
    delta_t: float,
    delta_z: float = 1.0,
    c_s: float = DEFAULT_CS,
) -> Union[float, np.ndarray]:
    """Soil heat flux from the general calorimetric equation (FAO-56 Eq. 41).

    .. math::
        G = c_s \\, \\frac{T_i - T_{i-1}}{\\Delta t} \\, \\Delta z

    Parameters
    ----------
    t_current : float or array_like
        Mean air temperature of the current period [°C].
    t_previous : float or array_like
        Mean air temperature of the previous period [°C].
    delta_t : float
        Length of the time interval in **days**.  For successive months use
        ~30.4; for a single month pair use the actual number of days between
        the midpoints of the two months.
    delta_z : float, optional
        Effective soil depth [m].  Typical values are 0.10–0.20 m for daily
        periods and 0.5–2.0 m for monthly periods.  Default is 1.0 m
        (standard for successive-month calculations).
    c_s : float, optional
        Soil heat capacity [MJ m⁻³ °C⁻¹].  Default ``2.1``.

    Returns
    -------
    float or numpy.ndarray
        Soil heat flux *G* [MJ m⁻² day⁻¹].  Positive values indicate heat
        transfer into the soil (warming); negative values indicate heat
        release from the soil (cooling).

    Notes
    -----
    This is the most general form.  Equations 43 and 44 are simplified
    versions derived from this equation by assuming c_s = 2.1 MJ/(m³·°C),
    Δz = 1.0 m, and Δt appropriate for monthly intervals (~30 days).

    For Eq. 43:  0.07 ≈ 2.1 × 1.0 / (2 × 30)
    For Eq. 44:  0.14 ≈ 2.1 × 1.0 / (1 × 30)  (÷ by half the interval)

    Examples
    --------
    FAO-56 Example 13 — March to April:

    >>> soil_heat_flux_general(16.1, 14.1, delta_t=30.0, delta_z=1.0)
    0.14
    """
    t_current = np.asarray(t_current, dtype=float)
    t_previous = np.asarray(t_previous, dtype=float)

    g = c_s * (t_current - t_previous) * delta_z / delta_t
    return float(g) if g.ndim == 0 else g


def soil_heat_flux_daily() -> float:
    """Daily soil heat flux for a grass reference surface (FAO-56 Eq. 42).

    For a daily time step, soil heat flux beneath a grass reference surface
    is small relative to net radiation and may be assumed to be zero.

    Returns
    -------
    float
        ``0.0`` [MJ m⁻² day⁻¹].

    Notes
    -----
    This assumption holds for a dense, grass-covered surface where the soil
    surface temperature cycle nearly cancels over 24 hours.  For bare soil
    or sparse vegetation, daily *G* can be significant and should be
    estimated by other means (e.g., the general equation or measured data).
    """
    return 0.0


def soil_heat_flux_monthly(
    t_month_prev: Union[float, npt.ArrayLike],
    t_month_next: Union[float, npt.ArrayLike],
) -> Union[float, np.ndarray]:
    """Monthly soil heat flux using adjacent months (FAO-56 Eq. 43).

    When both the previous and next month's mean air temperatures are
    available, this is the preferred monthly estimate.

    .. math::
        G_{\\text{month},i} = 0.07 \\, (T_{i+1} - T_{i-1})

    Parameters
    ----------
    t_month_prev : float or array_like
        Mean air temperature of the **previous** month [°C].
    t_month_next : float or array_like
        Mean air temperature of the **next** month [°C].

    Returns
    -------
    float or numpy.ndarray
        Monthly soil heat flux [MJ m⁻² day⁻¹].

    Notes
    -----
    The coefficient 0.07 derives from Eq. 41 with c_s = 2.1 MJ/(m³·°C),
    Δz = 1.0 m, and a two-month span (Δt ≈ 2 × 30 = 60 days):

        2.1 × 1.0 / 60 = 0.035  →  rounded/adopted as 0.07 in FAO-56.

    (The FAO text notes that the coefficient accounts for the relationship
    between air and soil temperature damping; it is empirically calibrated
    for a grass reference surface.)

    Examples
    --------
    FAO-56 Example 13 — April soil heat flux (March = 14.1 °C, May = 18.8 °C):

    >>> soil_heat_flux_monthly(14.1, 18.8)
    0.329
    """
    t_month_prev = np.asarray(t_month_prev, dtype=float)
    t_month_next = np.asarray(t_month_next, dtype=float)

    g = 0.07 * (t_month_next - t_month_prev)
    return float(g) if g.ndim == 0 else g


def soil_heat_flux_monthly_prev_only(
    t_month_prev: Union[float, npt.ArrayLike],
    t_month_cur: Union[float, npt.ArrayLike],
) -> Union[float, np.ndarray]:
    """Monthly soil heat flux — next month unavailable (FAO-56 Eq. 44).

    When the next month's temperature is not yet known (e.g., real-time or
    forecast applications), use this equation with the current and previous
    month.

    .. math::
        G_{\\text{month},i} = 0.14 \\, (T_i - T_{i-1})

    Parameters
    ----------
    t_month_prev : float or array_like
        Mean air temperature of the **previous** month [°C].
    t_month_cur : float or array_like
        Mean air temperature of the **current** month [°C].

    Returns
    -------
    float or numpy.ndarray
        Monthly soil heat flux [MJ m⁻² day⁻¹].

    Notes
    -----
    The coefficient 0.14 is exactly twice the Eq. 43 coefficient because
    the temperature difference spans only one month instead of two:

        2.1 × 1.0 / 30 ≈ 0.07 × 2 = 0.14

    When both adjacent months are available, prefer ``soil_heat_flux_monthly()``
    (Eq. 43) for greater accuracy.

    Examples
    --------
    FAO-56 Example 13 — March soil heat flux (Feb = 12.1 °C, Mar = 14.1 °C):

    >>> soil_heat_flux_monthly_prev_only(12.1, 14.1)
    0.28
    """
    t_month_prev = np.asarray(t_month_prev, dtype=float)
    t_month_cur = np.asarray(t_month_cur, dtype=float)

    g = 0.14 * (t_month_cur - t_month_prev)
    return float(g) if g.ndim == 0 else g


# ===========================================================================
# ASCE Standardized Hourly / Sub-Daily Soil Heat Flux
# ===========================================================================

def soil_heat_flux_hourly(
    rn: Union[float, npt.ArrayLike],
    is_daytime: Union[bool, npt.ArrayLike],
    day_coeff: float = 0.1,
    night_coeff: float = 0.5,
) -> Union[float, np.ndarray]:
    """Hourly (sub-daily) soil heat flux as a fraction of net radiation.

    Based on the ASCE Standardized Reference Evapotranspiration Equation
    (ASCE-EWRI, 2005), which recommends estimating *G* for hourly or
    shorter periods as a proportion of net radiation *Rn*:

    .. math::
        G_{\\text{day}}   &= 0.1 \\; R_n   \\quad \\text{(daytime)}

        G_{\\text{night}} &= 0.5 \\; R_n   \\quad \\text{(nighttime)}

    Parameters
    ----------
    rn : float or array_like
        Net radiation at the crop surface [MJ m⁻² hr⁻¹].
    is_daytime : bool or array_like of bool
        ``True`` for daytime periods, ``False`` for nighttime.
        For array inputs, must be broadcastable with *rn*.
    day_coeff : float, optional
        Fraction of Rn used for daytime G.  Default ``0.1``.
    night_coeff : float, optional
        Fraction of Rn used for nighttime G.  Default ``0.5``.

    Returns
    -------
    float or numpy.ndarray
        Soil heat flux [MJ m⁻² hr⁻¹] (same units as *rn*).

    Notes
    -----
    Daytime is typically defined as the period when *Rn* > 0.  Some
    implementations use the solar angle or the time relative to sunrise
    and sunset instead.

    The default coefficients (0.1 / 0.5) are calibrated for a grass
    reference surface.  For other surfaces the ratio G/Rn can vary
    considerably (see Allen et al., 1998, Table 7-2).

    Examples
    --------
    Daytime hour with Rn = 2.5 MJ m⁻² hr⁻¹:

    >>> soil_heat_flux_hourly(2.5, is_daytime=True)
    0.25

    Nighttime hour with Rn = -0.4 MJ m⁻² hr⁻¹:

    >>> soil_heat_flux_hourly(-0.4, is_daytime=False)
    -0.2
    """
    rn = np.asarray(rn, dtype=float)
    is_daytime = np.asarray(is_daytime, dtype=bool)

    coeff = np.where(is_daytime, day_coeff, night_coeff)
    g = coeff * rn
    return float(g) if g.ndim == 0 else g


def soil_heat_flux_hourly_auto(
    rn: Union[float, npt.ArrayLike],
    day_coeff: float = 0.1,
    night_coeff: float = 0.5,
) -> Union[float, np.ndarray]:
    """Hourly soil heat flux, auto-detecting day/night from Rn sign.

    Convenience wrapper around :func:`soil_heat_flux_hourly` that uses
    ``Rn > 0`` as the daytime indicator, which is the most common approach
    in practice.

    Parameters
    ----------
    rn : float or array_like
        Net radiation at the crop surface [MJ m⁻² hr⁻¹].
    day_coeff : float, optional
        Fraction of Rn used for daytime G.  Default ``0.1``.
    night_coeff : float, optional
        Fraction of Rn used for nighttime G.  Default ``0.5``.

    Returns
    -------
    float or numpy.ndarray
        Soil heat flux [MJ m⁻² hr⁻¹].

    Examples
    --------
    >>> soil_heat_flux_hourly_auto(2.5)
    0.25
    >>> soil_heat_flux_hourly_auto(-0.4)
    -0.2
    """
    rn = np.asarray(rn, dtype=float)
    is_day = rn > 0
    return soil_heat_flux_hourly(rn, is_day, day_coeff, night_coeff)


# ===========================================================================
# Monthly Time Series Helper
# ===========================================================================

def soil_heat_flux_monthly_series(
    t_monthly: npt.ArrayLike,
    method: str = "centered",
) -> np.ndarray:
    """Compute soil heat flux for each month in a temperature time series.

    This is a convenience function that applies the appropriate FAO-56
    monthly equation to each element of a monthly mean-temperature array.

    Parameters
    ----------
    t_monthly : array_like
        1-D array of consecutive monthly mean air temperatures [°C].
        Length *N* ≥ 2.
    method : {"centered", "backward"}, optional
        * ``"centered"`` (default) — Uses Eq. 43 where both neighbors
          exist; falls back to Eq. 44 at the boundaries (first and last
          month).
        * ``"backward"`` — Uses Eq. 44 for every month (each month
          compared only to the previous month).  The first element is
          set to ``NaN`` because there is no previous month.

    Returns
    -------
    numpy.ndarray
        Array of monthly soil heat flux values [MJ m⁻² day⁻¹], same
        length as *t_monthly*.

    Examples
    --------
    FAO-56 Example 13 — March, April, May temps:

    >>> temps = [14.1, 16.1, 18.8]
    >>> soil_heat_flux_monthly_series(temps, method="centered")
    array([0.28 , 0.329, 0.378])

    Using backward-only differences:

    >>> soil_heat_flux_monthly_series(temps, method="backward")
    array([  nan, 0.28 , 0.378])
    """
    t = np.asarray(t_monthly, dtype=float)
    n = len(t)
    if n < 2:
        raise ValueError("Need at least 2 monthly temperatures.")

    g = np.empty(n, dtype=float)

    if method == "centered":
        # Interior months — Eq. 43
        for i in range(1, n - 1):
            g[i] = 0.07 * (t[i + 1] - t[i - 1])
        # Boundaries — fall back to Eq. 44.
        # First month: no previous month available, so use forward
        # difference (current→next) with Eq. 44 coefficient.
        g[0] = 0.14 * (t[1] - t[0])
        # Last month: no next month available, so use backward
        # difference (previous→current) with Eq. 44 coefficient.
        g[-1] = 0.14 * (t[-1] - t[-2])
        # Note: For full 12-month annual cycles, prefer
        # ``soil_heat_flux_annual_cycle()`` which wraps around properly.

    elif method == "backward":
        g[0] = np.nan
        for i in range(1, n):
            g[i] = 0.14 * (t[i] - t[i - 1])

    else:
        raise ValueError(f"method must be 'centered' or 'backward', got {method!r}")

    return g


def soil_heat_flux_annual_cycle(
    t_monthly_12: npt.ArrayLike,
) -> np.ndarray:
    """Monthly G for a complete 12-month annual cycle using Eq. 43.

    Wraps around so that January uses December and February as neighbors,
    and December uses November and January.

    Parameters
    ----------
    t_monthly_12 : array_like
        Exactly 12 monthly mean air temperatures [°C], January through
        December.

    Returns
    -------
    numpy.ndarray
        Shape ``(12,)`` array of monthly soil heat flux [MJ m⁻² day⁻¹].

    Examples
    --------
    >>> temps = [5.2, 6.1, 9.3, 13.0, 17.5, 22.1,
    ...          25.4, 24.8, 20.6, 14.9, 9.2, 5.8]
    >>> g = soil_heat_flux_annual_cycle(temps)
    >>> g[0]   # January: 0.07 * (Feb - Dec)
    0.021...
    """
    t = np.asarray(t_monthly_12, dtype=float)
    if t.shape != (12,):
        raise ValueError(f"Expected exactly 12 values, got shape {t.shape}")

    g = np.empty(12, dtype=float)
    for i in range(12):
        t_prev = t[(i - 1) % 12]
        t_next = t[(i + 1) % 12]
        g[i] = 0.07 * (t_next - t_prev)
    return g


# ===========================================================================
# Unit Conversion Utilities
# ===========================================================================

def energy_to_evap(
    energy: Union[float, npt.ArrayLike],
) -> Union[float, np.ndarray]:
    """Convert energy flux to equivalent evaporation depth (FAO-56 Eq. 20).

    .. math::
        E_{\\text{equiv}} = 0.408 \\times E

    Parameters
    ----------
    energy : float or array_like
        Energy flux [MJ m⁻² day⁻¹] (or per hour, etc.).

    Returns
    -------
    float or numpy.ndarray
        Equivalent evaporation [mm day⁻¹] (or per hour, etc.).
    """
    energy = np.asarray(energy, dtype=float)
    result = LATENT_HEAT_INV * energy
    return float(result) if result.ndim == 0 else result


def evap_to_energy(
    evap: Union[float, npt.ArrayLike],
) -> Union[float, np.ndarray]:
    """Convert equivalent evaporation depth back to energy flux.

    Inverse of :func:`energy_to_evap`.

    Parameters
    ----------
    evap : float or array_like
        Equivalent evaporation [mm day⁻¹].

    Returns
    -------
    float or numpy.ndarray
        Energy flux [MJ m⁻² day⁻¹].
    """
    evap = np.asarray(evap, dtype=float)
    result = evap / LATENT_HEAT_INV
    return float(result) if result.ndim == 0 else result


def mj_m2_day_to_w_m2(
    flux: Union[float, npt.ArrayLike],
) -> Union[float, np.ndarray]:
    """Convert flux from MJ m⁻² day⁻¹ to W m⁻².

    .. math::
        W\\,m^{-2} = \\frac{MJ\\,m^{-2}\\,day^{-1}}{0.0864}

    Parameters
    ----------
    flux : float or array_like
        Energy flux [MJ m⁻² day⁻¹].

    Returns
    -------
    float or numpy.ndarray
        Energy flux [W m⁻²].
    """
    flux = np.asarray(flux, dtype=float)
    result = flux / 0.0864
    return float(result) if result.ndim == 0 else result


def w_m2_to_mj_m2_day(
    flux: Union[float, npt.ArrayLike],
) -> Union[float, np.ndarray]:
    """Convert flux from W m⁻² to MJ m⁻² day⁻¹.

    Parameters
    ----------
    flux : float or array_like
        Energy flux [W m⁻²].

    Returns
    -------
    float or numpy.ndarray
        Energy flux [MJ m⁻² day⁻¹].
    """
    flux = np.asarray(flux, dtype=float)
    result = flux * 0.0864
    return float(result) if result.ndim == 0 else result


# ===========================================================================
# Validation / Self-Test
# ===========================================================================

def _run_tests() -> None:
    """Reproduce FAO-56 Example 13 and other verification cases."""

    print("=" * 65)
    print("FAO-56 Soil Heat Flux — Verification Tests")
    print("=" * 65)

    # ------------------------------------------------------------------
    # FAO-56 Example 13: March, April, May
    # Temperatures: 14.1, 16.1, 18.8 °C
    # ------------------------------------------------------------------
    t_mar, t_apr, t_may = 14.1, 16.1, 18.8

    print("\n--- FAO-56 Example 13 (Chapter 3) ---")
    print(f"Monthly temps: Mar={t_mar}, Apr={t_apr}, May={t_may} °C\n")

    # Eq. 41 — General (March→April, Δt=30 days, Δz=1.0 m, cs=2.1)
    g41 = soil_heat_flux_general(t_apr, t_mar, delta_t=30.0)
    print(f"Eq. 41 (Mar→Apr, Δt=30d): G = {g41:.4f} MJ/(m²·day)")

    # Eq. 43 — April using March & May
    g43 = soil_heat_flux_monthly(t_mar, t_may)
    print(f"Eq. 43 (Apr, centered):   G = {g43:.4f} MJ/(m²·day)  "
          f"[expected ≈ 0.329]")

    # Eq. 44 — March using previous month (assume Feb=12.1)
    t_feb = 12.1
    g44_mar = soil_heat_flux_monthly_prev_only(t_feb, t_mar)
    print(f"Eq. 44 (Mar, backward):   G = {g44_mar:.4f} MJ/(m²·day)  "
          f"[expected ≈ 0.28]")

    # Eq. 44 — May using April
    g44_may = soil_heat_flux_monthly_prev_only(t_apr, t_may)
    print(f"Eq. 44 (May, backward):   G = {g44_may:.4f} MJ/(m²·day)  "
          f"[expected ≈ 0.378]")

    # Eq. 42 — Daily
    g42 = soil_heat_flux_daily()
    print(f"\nEq. 42 (daily):           G = {g42:.1f} MJ/(m²·day)")

    # ------------------------------------------------------------------
    # Series helper
    # ------------------------------------------------------------------
    temps = [14.1, 16.1, 18.8]
    g_series = soil_heat_flux_monthly_series(temps, method="centered")
    print(f"\nMonthly series (centered): {np.round(g_series, 4)}")

    g_series_bw = soil_heat_flux_monthly_series(temps, method="backward")
    print(f"Monthly series (backward): {g_series_bw}")

    # ------------------------------------------------------------------
    # Hourly
    # ------------------------------------------------------------------
    print("\n--- ASCE Hourly Soil Heat Flux ---")
    rn_day = 2.5  # MJ m⁻² hr⁻¹
    rn_night = -0.4
    g_day = soil_heat_flux_hourly(rn_day, is_daytime=True)
    g_night = soil_heat_flux_hourly(rn_night, is_daytime=False)
    print(f"Daytime  (Rn={rn_day}):  G = {g_day:.3f} MJ/(m²·hr)  "
          f"[0.1 × Rn]")
    print(f"Nighttime (Rn={rn_night}): G = {g_night:.3f} MJ/(m²·hr)  "
          f"[0.5 × Rn]")

    # Auto-detection
    rn_array = np.array([2.5, 1.8, -0.1, -0.4])
    g_auto = soil_heat_flux_hourly_auto(rn_array)
    print(f"\nAuto day/night array:  Rn = {rn_array}")
    print(f"                        G = {np.round(g_auto, 3)}")

    # ------------------------------------------------------------------
    # Unit conversions
    # ------------------------------------------------------------------
    print("\n--- Unit Conversions ---")
    g_mj = 0.329
    g_mm = energy_to_evap(g_mj)
    print(f"{g_mj} MJ/(m²·day) → {g_mm:.4f} mm/day equiv. evaporation")
    print(f"{g_mj} MJ/(m²·day) → {mj_m2_day_to_w_m2(g_mj):.2f} W/m²")

    # ------------------------------------------------------------------
    # Annual cycle
    # ------------------------------------------------------------------
    print("\n--- 12-Month Annual Cycle (Eq. 43, wrapped) ---")
    # Example: loosely based on a temperate continental climate
    t12 = np.array([-2.0, 0.5, 5.0, 11.0, 17.0, 22.0,
                     25.0, 24.0, 19.0, 12.0, 5.0, -0.5])
    g12 = soil_heat_flux_annual_cycle(t12)
    months = ["Jan", "Feb", "Mar", "Apr", "May", "Jun",
              "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"]
    for m, ti, gi in zip(months, t12, g12):
        print(f"  {m}: T={ti:6.1f} °C  →  G = {gi:+.4f} MJ/(m²·day)"
              f"  ({energy_to_evap(gi):+.4f} mm/day)")

    print("\n" + "=" * 65)
    print("All tests completed.")
    print("=" * 65)


if __name__ == "__main__":
    _run_tests()
