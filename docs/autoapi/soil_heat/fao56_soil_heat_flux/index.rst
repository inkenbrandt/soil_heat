soil_heat.fao56_soil_heat_flux
==============================

.. py:module:: soil_heat.fao56_soil_heat_flux

.. autoapi-nested-parse::

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



Attributes
----------

.. autoapisummary::

   soil_heat.fao56_soil_heat_flux.DEFAULT_CS
   soil_heat.fao56_soil_heat_flux.LATENT_HEAT_INV


Functions
---------

.. autoapisummary::

   soil_heat.fao56_soil_heat_flux.soil_heat_flux_general
   soil_heat.fao56_soil_heat_flux.soil_heat_flux_daily
   soil_heat.fao56_soil_heat_flux.soil_heat_flux_monthly
   soil_heat.fao56_soil_heat_flux.soil_heat_flux_monthly_prev_only
   soil_heat.fao56_soil_heat_flux.soil_heat_flux_hourly
   soil_heat.fao56_soil_heat_flux.soil_heat_flux_hourly_auto
   soil_heat.fao56_soil_heat_flux.soil_heat_flux_monthly_series
   soil_heat.fao56_soil_heat_flux.soil_heat_flux_annual_cycle
   soil_heat.fao56_soil_heat_flux.energy_to_evap
   soil_heat.fao56_soil_heat_flux.evap_to_energy
   soil_heat.fao56_soil_heat_flux.mj_m2_day_to_w_m2
   soil_heat.fao56_soil_heat_flux.w_m2_to_mj_m2_day
   soil_heat.fao56_soil_heat_flux._run_tests


Module Contents
---------------

.. py:data:: DEFAULT_CS
   :value: 2.1


.. py:data:: LATENT_HEAT_INV
   :value: 0.408


.. py:function:: soil_heat_flux_general(t_current: Union[float, numpy.typing.ArrayLike], t_previous: Union[float, numpy.typing.ArrayLike], delta_t: float, delta_z: float = 1.0, c_s: float = DEFAULT_CS) -> Union[float, numpy.ndarray]

   Soil heat flux from the general calorimetric equation (FAO-56 Eq. 41).

   .. math::
       G = c_s \, \frac{T_i - T_{i-1}}{\Delta t} \, \Delta z

   :param t_current: Mean air temperature of the current period [°C].
   :type t_current: float or array_like
   :param t_previous: Mean air temperature of the previous period [°C].
   :type t_previous: float or array_like
   :param delta_t: Length of the time interval in **days**.  For successive months use
                   ~30.4; for a single month pair use the actual number of days between
                   the midpoints of the two months.
   :type delta_t: float
   :param delta_z: Effective soil depth [m].  Typical values are 0.10–0.20 m for daily
                   periods and 0.5–2.0 m for monthly periods.  Default is 1.0 m
                   (standard for successive-month calculations).
   :type delta_z: float, optional
   :param c_s: Soil heat capacity [MJ m⁻³ °C⁻¹].  Default ``2.1``.
   :type c_s: float, optional

   :returns: Soil heat flux *G* [MJ m⁻² day⁻¹].  Positive values indicate heat
             transfer into the soil (warming); negative values indicate heat
             release from the soil (cooling).
   :rtype: float or numpy.ndarray

   .. rubric:: Notes

   This is the most general form.  Equations 43 and 44 are simplified
   versions derived from this equation by assuming c_s = 2.1 MJ/(m³·°C),
   Δz = 1.0 m, and Δt appropriate for monthly intervals (~30 days).

   For Eq. 43:  0.07 ≈ 2.1 × 1.0 / (2 × 30)
   For Eq. 44:  0.14 ≈ 2.1 × 1.0 / (1 × 30)  (÷ by half the interval)

   .. rubric:: Examples

   FAO-56 Example 13 — March to April:

   >>> soil_heat_flux_general(16.1, 14.1, delta_t=30.0, delta_z=1.0)
   0.14


.. py:function:: soil_heat_flux_daily() -> float

   Daily soil heat flux for a grass reference surface (FAO-56 Eq. 42).

   For a daily time step, soil heat flux beneath a grass reference surface
   is small relative to net radiation and may be assumed to be zero.

   :returns: ``0.0`` [MJ m⁻² day⁻¹].
   :rtype: float

   .. rubric:: Notes

   This assumption holds for a dense, grass-covered surface where the soil
   surface temperature cycle nearly cancels over 24 hours.  For bare soil
   or sparse vegetation, daily *G* can be significant and should be
   estimated by other means (e.g., the general equation or measured data).


.. py:function:: soil_heat_flux_monthly(t_month_prev: Union[float, numpy.typing.ArrayLike], t_month_next: Union[float, numpy.typing.ArrayLike]) -> Union[float, numpy.ndarray]

   Monthly soil heat flux using adjacent months (FAO-56 Eq. 43).

   When both the previous and next month's mean air temperatures are
   available, this is the preferred monthly estimate.

   .. math::
       G_{\text{month},i} = 0.07 \, (T_{i+1} - T_{i-1})

   :param t_month_prev: Mean air temperature of the **previous** month [°C].
   :type t_month_prev: float or array_like
   :param t_month_next: Mean air temperature of the **next** month [°C].
   :type t_month_next: float or array_like

   :returns: Monthly soil heat flux [MJ m⁻² day⁻¹].
   :rtype: float or numpy.ndarray

   .. rubric:: Notes

   The coefficient 0.07 derives from Eq. 41 with c_s = 2.1 MJ/(m³·°C),
   Δz = 1.0 m, and a two-month span (Δt ≈ 2 × 30 = 60 days):

       2.1 × 1.0 / 60 = 0.035  →  rounded/adopted as 0.07 in FAO-56.

   (The FAO text notes that the coefficient accounts for the relationship
   between air and soil temperature damping; it is empirically calibrated
   for a grass reference surface.)

   .. rubric:: Examples

   FAO-56 Example 13 — April soil heat flux (March = 14.1 °C, May = 18.8 °C):

   >>> soil_heat_flux_monthly(14.1, 18.8)
   0.329


.. py:function:: soil_heat_flux_monthly_prev_only(t_month_prev: Union[float, numpy.typing.ArrayLike], t_month_cur: Union[float, numpy.typing.ArrayLike]) -> Union[float, numpy.ndarray]

   Monthly soil heat flux — next month unavailable (FAO-56 Eq. 44).

   When the next month's temperature is not yet known (e.g., real-time or
   forecast applications), use this equation with the current and previous
   month.

   .. math::
       G_{\text{month},i} = 0.14 \, (T_i - T_{i-1})

   :param t_month_prev: Mean air temperature of the **previous** month [°C].
   :type t_month_prev: float or array_like
   :param t_month_cur: Mean air temperature of the **current** month [°C].
   :type t_month_cur: float or array_like

   :returns: Monthly soil heat flux [MJ m⁻² day⁻¹].
   :rtype: float or numpy.ndarray

   .. rubric:: Notes

   The coefficient 0.14 is exactly twice the Eq. 43 coefficient because
   the temperature difference spans only one month instead of two:

       2.1 × 1.0 / 30 ≈ 0.07 × 2 = 0.14

   When both adjacent months are available, prefer ``soil_heat_flux_monthly()``
   (Eq. 43) for greater accuracy.

   .. rubric:: Examples

   FAO-56 Example 13 — March soil heat flux (Feb = 12.1 °C, Mar = 14.1 °C):

   >>> soil_heat_flux_monthly_prev_only(12.1, 14.1)
   0.28


.. py:function:: soil_heat_flux_hourly(rn: Union[float, numpy.typing.ArrayLike], is_daytime: Union[bool, numpy.typing.ArrayLike], day_coeff: float = 0.1, night_coeff: float = 0.5) -> Union[float, numpy.ndarray]

   Hourly (sub-daily) soil heat flux as a fraction of net radiation.

   Based on the ASCE Standardized Reference Evapotranspiration Equation
   (ASCE-EWRI, 2005), which recommends estimating *G* for hourly or
   shorter periods as a proportion of net radiation *Rn*:

   .. math::
       G_{\text{day}}   &= 0.1 \; R_n   \quad \text{(daytime)}

       G_{\text{night}} &= 0.5 \; R_n   \quad \text{(nighttime)}

   :param rn: Net radiation at the crop surface [MJ m⁻² hr⁻¹].
   :type rn: float or array_like
   :param is_daytime: ``True`` for daytime periods, ``False`` for nighttime.
                      For array inputs, must be broadcastable with *rn*.
   :type is_daytime: bool or array_like of bool
   :param day_coeff: Fraction of Rn used for daytime G.  Default ``0.1``.
   :type day_coeff: float, optional
   :param night_coeff: Fraction of Rn used for nighttime G.  Default ``0.5``.
   :type night_coeff: float, optional

   :returns: Soil heat flux [MJ m⁻² hr⁻¹] (same units as *rn*).
   :rtype: float or numpy.ndarray

   .. rubric:: Notes

   Daytime is typically defined as the period when *Rn* > 0.  Some
   implementations use the solar angle or the time relative to sunrise
   and sunset instead.

   The default coefficients (0.1 / 0.5) are calibrated for a grass
   reference surface.  For other surfaces the ratio G/Rn can vary
   considerably (see Allen et al., 1998, Table 7-2).

   .. rubric:: Examples

   Daytime hour with Rn = 2.5 MJ m⁻² hr⁻¹:

   >>> soil_heat_flux_hourly(2.5, is_daytime=True)
   0.25

   Nighttime hour with Rn = -0.4 MJ m⁻² hr⁻¹:

   >>> soil_heat_flux_hourly(-0.4, is_daytime=False)
   -0.2


.. py:function:: soil_heat_flux_hourly_auto(rn: Union[float, numpy.typing.ArrayLike], day_coeff: float = 0.1, night_coeff: float = 0.5) -> Union[float, numpy.ndarray]

   Hourly soil heat flux, auto-detecting day/night from Rn sign.

   Convenience wrapper around :func:`soil_heat_flux_hourly` that uses
   ``Rn > 0`` as the daytime indicator, which is the most common approach
   in practice.

   :param rn: Net radiation at the crop surface [MJ m⁻² hr⁻¹].
   :type rn: float or array_like
   :param day_coeff: Fraction of Rn used for daytime G.  Default ``0.1``.
   :type day_coeff: float, optional
   :param night_coeff: Fraction of Rn used for nighttime G.  Default ``0.5``.
   :type night_coeff: float, optional

   :returns: Soil heat flux [MJ m⁻² hr⁻¹].
   :rtype: float or numpy.ndarray

   .. rubric:: Examples

   >>> soil_heat_flux_hourly_auto(2.5)
   0.25
   >>> soil_heat_flux_hourly_auto(-0.4)
   -0.2


.. py:function:: soil_heat_flux_monthly_series(t_monthly: numpy.typing.ArrayLike, method: str = 'centered') -> numpy.ndarray

   Compute soil heat flux for each month in a temperature time series.

   This is a convenience function that applies the appropriate FAO-56
   monthly equation to each element of a monthly mean-temperature array.

   :param t_monthly: 1-D array of consecutive monthly mean air temperatures [°C].
                     Length *N* ≥ 2.
   :type t_monthly: array_like
   :param method:
                  * ``"centered"`` (default) — Uses Eq. 43 where both neighbors
                    exist; falls back to Eq. 44 at the boundaries (first and last
                    month).
                  * ``"backward"`` — Uses Eq. 44 for every month (each month
                    compared only to the previous month).  The first element is
                    set to ``NaN`` because there is no previous month.
   :type method: {"centered", "backward"}, optional

   :returns: Array of monthly soil heat flux values [MJ m⁻² day⁻¹], same
             length as *t_monthly*.
   :rtype: numpy.ndarray

   .. rubric:: Examples

   FAO-56 Example 13 — March, April, May temps:

   >>> temps = [14.1, 16.1, 18.8]
   >>> soil_heat_flux_monthly_series(temps, method="centered")
   array([0.28 , 0.329, 0.378])

   Using backward-only differences:

   >>> soil_heat_flux_monthly_series(temps, method="backward")
   array([  nan, 0.28 , 0.378])


.. py:function:: soil_heat_flux_annual_cycle(t_monthly_12: numpy.typing.ArrayLike) -> numpy.ndarray

   Monthly G for a complete 12-month annual cycle using Eq. 43.

   Wraps around so that January uses December and February as neighbors,
   and December uses November and January.

   :param t_monthly_12: Exactly 12 monthly mean air temperatures [°C], January through
                        December.
   :type t_monthly_12: array_like

   :returns: Shape ``(12,)`` array of monthly soil heat flux [MJ m⁻² day⁻¹].
   :rtype: numpy.ndarray

   .. rubric:: Examples

   >>> temps = [5.2, 6.1, 9.3, 13.0, 17.5, 22.1,
   ...          25.4, 24.8, 20.6, 14.9, 9.2, 5.8]
   >>> g = soil_heat_flux_annual_cycle(temps)
   >>> g[0]   # January: 0.07 * (Feb - Dec)
   0.021...


.. py:function:: energy_to_evap(energy: Union[float, numpy.typing.ArrayLike]) -> Union[float, numpy.ndarray]

   Convert energy flux to equivalent evaporation depth (FAO-56 Eq. 20).

   .. math::
       E_{\text{equiv}} = 0.408 \times E

   :param energy: Energy flux [MJ m⁻² day⁻¹] (or per hour, etc.).
   :type energy: float or array_like

   :returns: Equivalent evaporation [mm day⁻¹] (or per hour, etc.).
   :rtype: float or numpy.ndarray


.. py:function:: evap_to_energy(evap: Union[float, numpy.typing.ArrayLike]) -> Union[float, numpy.ndarray]

   Convert equivalent evaporation depth back to energy flux.

   Inverse of :func:`energy_to_evap`.

   :param evap: Equivalent evaporation [mm day⁻¹].
   :type evap: float or array_like

   :returns: Energy flux [MJ m⁻² day⁻¹].
   :rtype: float or numpy.ndarray


.. py:function:: mj_m2_day_to_w_m2(flux: Union[float, numpy.typing.ArrayLike]) -> Union[float, numpy.ndarray]

   Convert flux from MJ m⁻² day⁻¹ to W m⁻².

   .. math::
       W\,m^{-2} = \frac{MJ\,m^{-2}\,day^{-1}}{0.0864}

   :param flux: Energy flux [MJ m⁻² day⁻¹].
   :type flux: float or array_like

   :returns: Energy flux [W m⁻²].
   :rtype: float or numpy.ndarray


.. py:function:: w_m2_to_mj_m2_day(flux: Union[float, numpy.typing.ArrayLike]) -> Union[float, numpy.ndarray]

   Convert flux from W m⁻² to MJ m⁻² day⁻¹.

   :param flux: Energy flux [W m⁻²].
   :type flux: float or array_like

   :returns: Energy flux [MJ m⁻² day⁻¹].
   :rtype: float or numpy.ndarray


.. py:function:: _run_tests() -> None

   Reproduce FAO-56 Example 13 and other verification cases.
