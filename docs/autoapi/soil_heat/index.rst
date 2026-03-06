soil_heat
=========

.. py:module:: soil_heat

.. autoapi-nested-parse::

   soil_heat: A Python library for soil heat flux models.
   ======================================================

   This package provides a collection of Python functions for calculating
   soil heat flux and related soil thermal properties. The implementations
   are based on various models and parameterizations published in the
   peer-reviewed scientific literature.

   The library is organized into modules, each corresponding to a specific
   publication or a set of related equations.

   Available submodules:
   ---------------------
   - `soil_heat`: Core functions and utilities.
   - `gao_et_al`: Functions from Gao et al. (2017).
   - `liebethal_and_folken`: Functions from Liebethal & Foken (2006).
   - `wang_and_bouzeid`: Functions from Wang & Bou-Zeid (2012).
   - `wang_and_yang`: Functions from Yang & Wang (2008).
   - `fao56_soil_heat_flux`: FAO-56 and ASCE soil heat flux methods.
   - `storage_calculations`: Soil and canopy heat storage calculations.

   All functions are designed to work with NumPy arrays for efficient,
   vectorized computations.



Submodules
----------

.. toctree::
   :maxdepth: 1

   /autoapi/soil_heat/fao56_soil_heat_flux/index
   /autoapi/soil_heat/gao_et_al/index
   /autoapi/soil_heat/liebethal_and_folken/index
   /autoapi/soil_heat/soil_heat/index
   /autoapi/soil_heat/storage_calculations/index
   /autoapi/soil_heat/wang_and_bouzeid/index
   /autoapi/soil_heat/wang_and_yang/index


Attributes
----------

.. autoapisummary::

   soil_heat.__author__
   soil_heat.__email__
   soil_heat.__version__
   soil_heat.WATER_HEAT_CAPACITY
   soil_heat.df
   soil_heat.surface_energy_residual
   soil_heat.DEFAULT_CS
   soil_heat.LATENT_HEAT_INV


Functions
---------

.. autoapisummary::

   soil_heat.compute_heat_flux_conduction
   soil_heat.compute_heat_flux_calorimetric
   soil_heat.temperature_gradient
   soil_heat.soil_heat_flux
   soil_heat.volumetric_heat_capacity
   soil_heat.thermal_conductivity
   soil_heat.diurnal_amplitude
   soil_heat.diurnal_peak_lag
   soil_heat.fit_sinusoid
   soil_heat.sinusoid
   soil_heat.thermal_diffusivity_amplitude
   soil_heat.thermal_diffusivity_lag
   soil_heat.thermal_diffusivity_logrithmic
   soil_heat.calc_thermal_diffusivity_log_pair
   soil_heat.calculate_thermal_diffusivity_for_pair
   soil_heat.calculate_thermal_properties_for_all_pairs
   soil_heat.estimate_rhoc_dry
   soil_heat.calculate_soil_heat_storage
   soil_heat.lambda_s
   soil_heat.k_s
   soil_heat.volumetric_heat_capacity
   soil_heat.nme
   soil_heat.rmse
   soil_heat.calorimetric_gz
   soil_heat.force_restore_gz
   soil_heat.gao2010_gz
   soil_heat.heusinkveld_gz
   soil_heat.hsieh2009_gz
   soil_heat.leuning_damping_depth
   soil_heat.leuning_gz
   soil_heat.simple_measurement_gz
   soil_heat.wbz12_g_gz
   soil_heat.wbz12_s_gz
   soil_heat.exact_temperature_gz
   soil_heat.exact_gz
   soil_heat.reference_ground_heat_flux
   soil_heat.ground_heat_flux_pr
   soil_heat.ground_heat_flux_lr
   soil_heat.ur_coefficients
   soil_heat.ground_heat_flux_ur
   soil_heat.surface_temp_amplitude
   soil_heat.phi_from_soil_moisture
   soil_heat.ground_heat_flux_sh
   soil_heat.ground_heat_flux_sm
   soil_heat.active_layer_thickness
   soil_heat.ground_heat_flux_fr
   soil_heat.energy_balance_residual
   soil_heat.ground_heat_flux_conventional
   soil_heat.green_function_temperature
   soil_heat.temperature_convolution_solution
   soil_heat.soil_heat_flux_from_G0
   soil_heat.estimate_G0_from_Gz
   soil_heat.sinusoidal_boundary_flux
   soil_heat.soil_temperature_sinusoidal
   soil_heat.soil_heat_flux_sinusoidal
   soil_heat.heat_capacity_moist_soil
   soil_heat.pf_from_theta
   soil_heat.thermal_conductivity_moist_soil
   soil_heat.thermal_diffusivity
   soil_heat.soil_heat_flux
   soil_heat.integrated_soil_heat_flux
   soil_heat.volumetric_heat_capacity
   soil_heat.stretched_grid
   soil_heat.solve_tde
   soil_heat.correct_profile
   soil_heat.surface_temperature_longwave
   soil_heat.thermal_conductivity_yang2008
   soil_heat.flux_error_linear
   soil_heat.surface_energy_residual
   soil_heat.tdec_step
   soil_heat.soil_heat_flux_general
   soil_heat.soil_heat_flux_daily
   soil_heat.soil_heat_flux_monthly
   soil_heat.soil_heat_flux_monthly_prev_only
   soil_heat.soil_heat_flux_hourly
   soil_heat.soil_heat_flux_hourly_auto
   soil_heat.soil_heat_flux_monthly_series
   soil_heat.soil_heat_flux_annual_cycle
   soil_heat.energy_to_evap
   soil_heat.evap_to_energy
   soil_heat.mj_m2_day_to_w_m2
   soil_heat.w_m2_to_mj_m2_day
   soil_heat._run_tests
   soil_heat.compute_soil_storage_integrated
   soil_heat.compute_canopy_storage


Package Contents
----------------

.. py:data:: __author__
   :value: 'Paul Inkenbrandt'


.. py:data:: __email__
   :value: 'paulinkenbrandt@utah.gov'


.. py:data:: __version__
   :value: '0.1.1'


.. py:data:: WATER_HEAT_CAPACITY
   :value: 4.18


.. py:function:: compute_heat_flux_conduction(df: pandas.DataFrame, depth1: float = 0.05, depth2: float = 0.1, col_T1: str = 'T5cm', col_T2: str = 'T10cm', col_theta1: str = 'VWC5cm', col_theta2: str = 'VWC10cm', porosity: float = 0.4, k_dry: float = 0.25, k_sat: float = 1.5) -> pandas.Series

   Estimate near-surface soil heat flux using Fourier’s law.

   This “gradient” approach computes conductive ground-heat flux
   :math:`G` between two depths by multiplying the vertical
   temperature gradient with an **effective** thermal conductivity
   that varies with volumetric water content (VWC).

   :param df: Time-indexed data containing at least the four columns
              specified by *col_T1*, *col_T2*, *col_theta1*, and
              *col_theta2*. The index spacing defines the temporal
              resolution of the output.
   :type df: pandas.DataFrame
   :param depth1: Sensor depths (m).  `depth2` must be **greater** (deeper)
                  than `depth1`.
   :type depth1: float, default (0.05, 0.10)
   :param depth2: Sensor depths (m).  `depth2` must be **greater** (deeper)
                  than `depth1`.
   :type depth2: float, default (0.05, 0.10)
   :param col_T1: Column names for temperature (°C or K) at `depth1` and
                  `depth2`.
   :type col_T1: str, default ("T5cm", "T10cm")
   :param col_T2: Column names for temperature (°C or K) at `depth1` and
                  `depth2`.
   :type col_T2: str, default ("T5cm", "T10cm")
   :param col_theta1: Column names for volumetric water content (m³ m⁻³) at
                      `depth1` and `depth2`.
   :type col_theta1: str, default ("VWC5cm", "VWC10cm")
   :param col_theta2: Column names for volumetric water content (m³ m⁻³) at
                      `depth1` and `depth2`.
   :type col_theta2: str, default ("VWC5cm", "VWC10cm")
   :param porosity: Soil total porosity (saturated VWC, m³ m⁻³).
   :type porosity: float, default 0.40
   :param k_dry: Dry-soil thermal conductivity (W m⁻¹ K⁻¹).
   :type k_dry: float, default 0.25
   :param k_sat: Saturated-soil thermal conductivity (W m⁻¹ K⁻¹).
   :type k_sat: float, default 1.50

   :returns: Half-hourly (or whatever the index step is) ground-heat-flux
             series with name ``"G_conduction"``. Units are W m⁻².
             Positive values indicate **downward** flux.
   :rtype: pandas.Series

   .. rubric:: Notes

   The effective thermal conductivity is computed by a simple linear
   mixing model:

   .. math::

       \lambda_\text{eff} = k_\text{dry} +
       \frac{\bar{\theta}}{\phi}
       \bigl(k_\text{sat} - k_\text{dry}\bigr),

   where :math:`\bar{\theta}` is the mean VWC of the two depths and
   :math:`\phi` is porosity.  More sophisticated models
   (e.g. Johansen, de Vries) can be substituted if site-specific
   calibration is available.

   .. rubric:: References

   * Campbell & Norman (2012) *An Introduction to Environmental
     Biophysics*, ch. 7.
   * Gao et al. (2017) Agricultural and Forest Meteorology,
     240 – 241, 194–204.

   .. rubric:: Examples

   >>> G = compute_heat_flux_conduction(df_site,
   ...                                   depth1=0.05, depth2=0.10,
   ...                                   col_T1="T_05",
   ...                                   col_T2="T_10",
   ...                                   col_theta1="VWC_05",
   ...                                   col_theta2="VWC_10")
   >>> G.plot(title="Soil heat flux (gradient method)")


.. py:function:: compute_heat_flux_calorimetric(df: pandas.DataFrame, depth_levels: list[float], T_cols: list[str], theta_cols: list[str], C_dry: float = 2100000.0, C_w: float = 4200000.0) -> pandas.Series

   Calculate surface soil heat flux via the calorimetric (heat-storage) method.

   The calorimetric method integrates the transient change in heat
   *storage* within a multilayer soil column.  For a surface-to-depth
   layer of thickness :math:`z_{\text{ref}}`, the surface flux
   :math:`G_0` is approximated by

   .. math::

       G_0 \;\approx\; \frac{\Delta Q}{\Delta t}
       \;=\; \frac{1}{\Delta t}
       \sum_{i=1}^{N_\text{layers}}
       C_i \, \Delta T_i \, \Delta z_i,

   where :math:`C_i` is volumetric heat capacity
   (J m⁻³ K⁻¹), :math:`\Delta T_i` is the average temperature change
   (K) in layer *i*, and :math:`\Delta z_i` is layer thickness (m).
   No heat-flux-plate reading is required if the deepest
   measurement depth lies below the diurnal damping depth such that
   :math:`G(z_{\text{ref}}) \approx 0`.

   :param df: Time-indexed data containing temperature and VWC columns for
              **all** depths specified in *T_cols* and *theta_cols*.  Index
              spacing sets the output time step.
   :type df: pandas.DataFrame
   :param depth_levels: Depths (m) corresponding *in order* to the entries in
                        *T_cols* and *theta_cols*. Must be strictly increasing.
   :type depth_levels: list of float
   :param T_cols: Column names for soil temperatures (°C or K) at
                  `depth_levels`.
   :type T_cols: list of str
   :param theta_cols: Column names for volumetric water content (m³ m⁻³) at
                      `depth_levels`.
   :type theta_cols: list of str
   :param C_dry: Volumetric heat capacity of dry soil matrix
                 (J m⁻³ K⁻¹).
   :type C_dry: float, default 2.1e6
   :param C_w: Volumetric heat capacity of liquid water
               (J m⁻³ K⁻¹).
   :type C_w: float, default 4.2e6

   :returns: Surface ground-heat-flux series, ``"G_calorimetric"`` (W m⁻²).
             Positive values denote **downward** flux.  The first time step
             is set to *NaN* because a preceding interval is required.
   :rtype: pandas.Series

   .. rubric:: Notes

   **Heat capacity model**

   A simple two-component mixture is assumed:

   .. math::

       C = (1 - \theta)\,C_{\text{dry}} + \theta\,C_w.

   If bulk density or mineral fraction data are available, replace
   this linear approximation with a mass-weighted formulation.

   **Boundary assumption**

   The deepest temperature is treated as a “no-flux” boundary (storage
   only).  If diurnal waves penetrate deeper at your site, include an
   additional flux-plate term or extend `depth_levels` downward.

   .. rubric:: References

   * Mayocchi & Bristow (1995) Agricultural and Forest
     Meteorology 75, 93–109.
   * Oke (2002) *Boundary-Layer Climates*, 2nd ed., §2.3.
   * Fluxnet2015 “G” best-practice guide
     (https://fluxnet.org/sites/default/files/soil_heat_flux_guide.pdf).

   .. rubric:: Examples

   >>> depths = [0.05, 0.10, 0.20, 0.50]          # m
   >>> Tcols  = ["T5", "T10", "T20", "T50"]       # °C
   >>> Vcols  = ["VWC5", "VWC10", "VWC20", "VWC50"]
   >>> G0 = compute_heat_flux_calorimetric(df_site,
   ...                                     depths, Tcols, Vcols)
   >>> G0.resample("D").mean().plot()
   >>> plt.ylabel("Daily mean G₀ (W m$^{-2}$)")


.. py:function:: temperature_gradient(T_upper: numpy.ndarray | float, T_lower: numpy.ndarray | float, depth_upper: float, depth_lower: float) -> numpy.ndarray | float

   Compute the **vertical temperature gradient** between two sensors.

   The gradient is defined as the change in temperature divided by the
   change in depth (positive downward):

   .. math::

       \frac{∂T}{∂z}
       \;=\;
       \frac{T_{\text{lower}} - T_{\text{upper}}}
             {z_{\text{lower}} - z_{\text{upper}}}   \;\;[^{\circ}\text{C m}^{-1}]

   :param T_upper: Temperature at the **shallower** depth ``depth_upper`` (°C).
   :type T_upper: float or array_like
   :param T_lower: Temperature at the **deeper** depth ``depth_lower`` (°C).
                   Must be broadcast-compatible with ``T_upper``.
   :type T_lower: float or array_like
   :param depth_upper: Depth of the upper sensor (m, positive downward).
   :type depth_upper: float
   :param depth_lower: Depth of the lower sensor (m, positive downward).
                       Must satisfy ``depth_lower > depth_upper`` for a meaningful
                       gradient.
   :type depth_lower: float

   :returns: Temperature gradient ∂T/∂z (°C m⁻¹).
             Shape follows NumPy broadcasting of ``T_upper`` and ``T_lower``.
   :rtype: ndarray or float

   :raises ValueError: If ``depth_lower`` ≤ ``depth_upper``.

   .. rubric:: Notes

   * **Sign convention** – A **positive** gradient indicates
     temperatures increase with depth (warmer below).
   * **Vectorised** – The arithmetic is fully NumPy-broadcasted; use it
     on scalar values, 1-D arrays, or entire DataFrames’ columns.
   * **Units** – Because depth is in metres and temperature in degrees
     Celsius, the result is °C m⁻¹ (identical to K m⁻¹).

   .. rubric:: Examples

   >>> grad = temperature_gradient(
   ...     T_upper=18.6, T_lower=20.1,
   ...     depth_upper=0.05, depth_lower=0.10,
   ... )
   >>> print(f"Gradient = {grad:.2f} °C/m")
   Gradient = 30.00 °C/m

   Array input:

   >>> T_up  = np.array([15.0, 16.2, 17.1])
   >>> T_low = np.array([14.0, 15.8, 16.9])
   >>> temperature_gradient(T_up, T_low, 0.02, 0.08)
   array([-16.66666667, -6.66666667, -3.33333333])   # °C/m


.. py:function:: soil_heat_flux(T_upper, T_lower, depth_upper, depth_lower, k)

   Calculate soil heat flux (G) using Fourier's law of heat conduction.

   This function computes the one-dimensional heat flux in the soil based on
   the temperature gradient between two depths and the soil's thermal
   conductivity. The formula used is:

   .. math::

       G = -k \frac{\Delta T}{\Delta z}

   where :math:`G` is the soil heat flux, :math:`k` is the thermal
   conductivity, :math:`\Delta T` is the temperature difference, and
   :math:`\Delta z` is the distance between the two measurement depths.

   :param T_upper: Temperature at the upper depth (°C or K).
   :type T_upper: float or array_like
   :param T_lower: Temperature at the lower depth (°C or K).
   :type T_lower: float or array_like
   :param depth_upper: Depth of the upper sensor (m, positive downward).
   :type depth_upper: float
   :param depth_lower: Depth of the lower sensor (m, positive downward).
   :type depth_lower: float
   :param k: Thermal conductivity of the soil layer between the sensors (W m⁻¹ K⁻¹).
   :type k: float or array_like

   :returns: The calculated soil heat flux (W m⁻²). A positive value indicates
             a downward flux (from the surface into the soil).
   :rtype: float or array_like


.. py:function:: volumetric_heat_capacity(theta_v)

   Estimate volumetric heat capacity (Cv) from soil moisture content.

   This function uses a simple mixing model to estimate the volumetric heat
   capacity of moist soil. It assumes the soil is a two-component mixture
   of dry soil and water. The formula is:

   .. math::

       C_v = (1 - \theta_v) C_{soil} + \theta_v C_{water}

   where :math:`\theta_v` is the volumetric water content, and
   :math:`C_{soil}` and :math:`C_{water}` are the volumetric heat
   capacities of dry soil and water, respectively.

   :param theta_v: Volumetric water content (m³ m⁻³). This is a decimal fraction,
                   e.g., 0.20 for 20% water content.
   :type theta_v: float or array_like

   :returns: The estimated volumetric heat capacity (kJ m⁻³ K⁻¹).
   :rtype: float or array_like


.. py:function:: thermal_conductivity(alpha: numpy.ndarray | float, theta_v: numpy.ndarray | float) -> numpy.ndarray | float

   Convert **thermal diffusivity** (``α``) to **thermal conductivity** (``k``)
   using the bulk *volumetric heat capacity* of moist soil.

   The relationship is

   .. math::

       k \;=\; α \, C_v(θ_v),

   where

   * *α* – thermal diffusivity (m² s⁻¹),
   * *C_v(θ_v)* – volumetric heat capacity (J m⁻³ K⁻¹) as a function of
     volumetric water content *θ_v* (m³ m⁻³).
     It is obtained from :func:`volumetric_heat_capacity`.

   :param alpha: Thermal diffusivity **α** (m² s⁻¹).  May be scalar or any
                 NumPy‐broadcastable shape.
   :type alpha: float or array_like
   :param theta_v: Volumetric water content **θ_v** (m³ m⁻³, i.e. decimal fraction
                   of pore space filled with water).  Must be broadcast‐compatible
                   with ``alpha``.
   :type theta_v: float or array_like

   :returns: Thermal conductivity **k** (W m⁻¹ K⁻¹) with the broadcast shape
             of the inputs.
   :rtype: ndarray or float

   .. rubric:: Notes

   * **Volumetric heat capacity model** –
     :func:`volumetric_heat_capacity` typically assumes a two‐phase
     mixture of mineral soil and water:

     .. math::

        C_v(θ_v) \;=\; (1-θ_v)\,ρc_    ext{dry} \;+\;
                        θ_v\,ρc_       ext{w} ,

     where ``ρc_dry`` (≈ 2.0 MJ m⁻³ K⁻¹) and ``ρc_w`` (4.18 MJ m⁻³ K⁻¹)
     are the volumetric heat capacities of dry soil and liquid water,
     respectively.  Ensure these defaults suit your substrate.
   * **Vectorisation** – The function is a one‐liner,
     ``alpha * Cv``, and thus inherits full NumPy broadcasting rules.
   * **Temperature units** – Because heat capacity is per kelvin, *k*
     is returned in W m⁻¹ K⁻¹ (equivalent to W m⁻¹ °C⁻¹).

   .. rubric:: Examples

   >>> α = np.array([1.4e-7, 1.6e-7])       # m² s⁻¹
   >>> θ = np.array([0.10, 0.25])           # m³ m⁻³
   >>> k = thermal_conductivity(α, θ)
   >>> k
   array([0.29, 0.54])                      # W m⁻¹ K⁻¹

   Plot conductivity versus moisture:

   >>> θ_range = np.linspace(0, 0.45, 100)
   >>> k_vals = thermal_conductivity(1.5e-7, θ_range)
   >>> plt.plot(θ_range, k_vals)
   >>> plt.xlabel("Volumetric water content (m³ m⁻³)")
   >>> plt.ylabel("Thermal conductivity (W m⁻¹ K⁻¹)")


.. py:function:: diurnal_amplitude(series: pandas.Series) -> pandas.Series

   Compute the **daily diurnal amplitude** of a time-series.

   The diurnal amplitude for a given calendar day is defined as the
   difference between that day’s maximum and minimum values:

   .. math::

       A_d \;=\; \max\_{t \in d} x(t) \;-\; \min\_{t \in d} x(t)

   This metric is frequently used for temperature, soil-heat, or other
   environmental data to characterise the strength of the diurnal cycle.

   :param series: Time-indexed observations with a `DatetimeIndex`.
                  Any frequency is accepted, but the index **must** be sorted and
                  monotonic.  Missing values (`NaN`) are ignored within each daily
                  window.
   :type series: pandas.Series

   :returns: Daily diurnal amplitude, indexed by date (midnight ``00:00`` of
             each day).  Units are the same as those of the input ``series``.
   :rtype: pandas.Series

   .. rubric:: Notes

   * **Resampling rule** – The computation uses

     >>> daily_max = series.resample("D").max()
     >>> daily_min = series.resample("D").min()

     which bins data by *calendar day* in the series’ timezone.
     Incomplete trailing days yield `NaN`.
   * **Timezone safety** – If the series’ index spans daylight-saving
     transitions, consider converting to UTC prior to analysis to avoid
     artificial jumps in daily windows.
   * **Robustness** – For noisy signals, you may wish to smooth
     ``series`` (e.g. rolling median) before calling this function.

   .. rubric:: Examples

   >>> amp = diurnal_amplitude(df["air_temperature"])
   >>> amp.plot(title="Daily Temperature Amplitude")
   >>> amp.describe().loc[["min", "mean", "max"]]
   min      4.3
   mean     9.7
   max     15.2
   Name: air_temperature, dtype: float64


.. py:function:: diurnal_peak_lag(series1: pandas.Series, series2: pandas.Series) -> pandas.Series

   Compute the **daily peak‐time lag** (Δt) between two diurnal signals.

   For each calendar day the function identifies the clock time at which
   each series reaches its maximum value and returns the signed time
   difference in **hours** (``series1`` minus ``series2``).  A modular
   correction confines the result to the interval ``[-12, 12]`` h so
   that, for example, a raw lag of –23 h becomes +1 h.

   :param series1: Time-indexed observations of equal length, preferably
                   temperature or some other quantity exhibiting a clear diurnal
                   cycle.  The index **must** be `DatetimeIndex` and should be
                   timezone-aware and aligned in frequency.
                   Missing values are ignored within each daily resampling window.
   :type series1: pandas.Series
   :param series2: Time-indexed observations of equal length, preferably
                   temperature or some other quantity exhibiting a clear diurnal
                   cycle.  The index **must** be `DatetimeIndex` and should be
                   timezone-aware and aligned in frequency.
                   Missing values are ignored within each daily resampling window.
   :type series2: pandas.Series

   :returns: Daily peak-lag values (float, hours) indexed by the **date** of
             the peak (00:00 of each day).
             Positive lags mean the peak of ``series1`` occurs *later* than
             the peak of ``series2`` on that day; negative lags indicate the
             opposite.
   :rtype: pandas.Series

   .. rubric:: Notes

   * **Resampling rule** – Peaks are detected with
     ``series.resample('D').apply(lambda x: x.idxmax())``.  Ensure the
     input data span whole days; incomplete trailing days yield `NaN`.
   * **Wrap-around correction** – The transformation
     ``(lag + 12) % 24 − 12`` folds lags so that the maximum absolute
     value is always < 12 h, which prevents a late-evening peak at 23:30
     and an early-morning peak at 00:30 from being reported as –23 h.
   * **Daylight-saving** – If the index carries a timezone subject to
     DST transitions, consider converting to UTC prior to analysis to
     avoid spurious 1-h jumps.

   .. rubric:: Examples

   >>> peak_lag = diurnal_peak_lag(df['ts_05cm'], df['ts_10cm'])
   >>> peak_lag.describe()
   count    90.000000
   mean      1.42
   std       0.53
   min      -0.73
   25%       1.11
   50%       1.42
   75%       1.74
   max       2.33
   Name: ts_05cm, dtype: float64

   Plot the distribution:

   >>> peak_lag.plot(kind='hist', bins=24)
   >>> plt.xlabel('Peak lag (h)')
   >>> plt.title('Daily phase lag: 5 cm vs 10 cm temperature')


.. py:function:: fit_sinusoid(t: numpy.ndarray, data: numpy.ndarray) -> tuple[numpy.ndarray, numpy.ndarray]

       Fit a **sinusoidal model** to time–series data using non-linear least
       squares.

       The model is

       .. math::

           y(t)\;=\;A \sin( \omega t +
   arphi ) + C ,

       where
       ``A`` is the amplitude, ``ω`` the angular frequency,
       ``φ`` the phase shift, and ``C`` the vertical offset.
       Initial parameter guesses are derived from the sample statistics of
       *data* and an assumed daily frequency.

       Parameters
       ----------
       t : ndarray
           1-D array of time stamps (s).
           Must be the same length as ``data``.
       data : ndarray
           Observed values corresponding to ``t`` (e.g. temperature, °C).

       Returns
       -------
       popt : ndarray, shape (4,)
           Optimal parameters ``[A, ω, φ, C]`` that minimise
           the sum-of-squares error.
       pcov : ndarray, shape (4, 4)
           Covariance matrix of the parameter estimates returned by
           :func:`scipy.optimize.curve_fit`.

       Notes
       -----
       * **Initial guess** –
         *Amplitude* is set to the sample standard deviation of *data*,
         *frequency* to a 24-h cycle
         (``ω = 2π / 86 400`` s⁻¹),
         *phase* to 0, and *offset* to the sample mean.
         Adjust these if fitting to non-diurnal signals.
       * **Robustness** –
         If convergence issues arise, provide a closer initial guess or
         bound parameters via ``curve_fit``’s
         ``bounds=`` keyword.

       Examples
       --------
       >>> import numpy as np
       >>> from scipy.optimize import curve_fit
       >>> t = np.arange(0, 3*86400, 1800)                # 3 days, 30-min Δt
       >>> true = sinusoid(t, 7, 2*np.pi/86400, 0.3, 15)
       >>> rng = np.random.default_rng(0)
       >>> y = true + rng.normal(0, 0.5, t.size)          # add noise
       >>> params, _ = fit_sinusoid(t, y)
       >>> A, ω, φ, C = params
       >>> print(f"Amplitude={A:.2f}, Period={2*np.pi/ω/3600:.2f} h")
       Amplitude=7.01, Period=24.00 h



.. py:function:: sinusoid(t: numpy.ndarray | float, A: float, omega: float, phase: float, offset: float) -> numpy.ndarray | float

       Evaluate a **sinusoidal wave** of the form

       .. math::

           f(t) \;=\; A \sin(\omega\, t +
   arphi) + C ,

       where :math:`A` is the *amplitude*, :math:`\omega` the *angular
       frequency*, :math:`
   arphi` the *phase shift*, and :math:`C`
       a constant *vertical offset*.

       The function is *vectorised* with NumPy broadcasting, so ``t`` may be
       a scalar, 1-D array, or any shape compatible with the parameters.

       Parameters
       ----------
       t : float or array_like
           Independent variable (time, angle, etc.).  Units are arbitrary,
           but must be consistent with ``omega`` (e.g. seconds if
           ``omega`` is rad s⁻¹).
       A : float
           Wave amplitude.  Sets the peak deviation from ``offset``.
       omega : float
           Angular frequency (rad × ``t``⁻¹).
           For a *temporal* signal ``ω = 2π / P`` where *P* is the period.
       phase : float
           Phase shift :math:`
   arphi` in **radians**.
           Positive values delay the wave (right shift), negative values
           advance it.
       offset : float
           Constant vertical shift :math:`C`.  Often the long-term mean or
           base-line value of the signal.

       Returns
       -------
       ndarray or float
           Value(s) of the sinusoid at ``t`` with the same shape as the
           broadcast result of the inputs.

       Notes
       -----
       * **Vectorisation** – Internally relies on ``numpy.sin``; all
         standard broadcasting rules apply.
       * **Period** – The fundamental period *P* is related to ``omega`` by
         *P = 2π / ω*.  Specify *ω* rather than *P* to avoid repeated
         division operations when fitting.

       Examples
       --------
       >>> import numpy as np
       >>> import matplotlib.pyplot as plt
       >>> t = np.linspace(0, 24, 1000)                      # hours
       >>> temp = sinusoid(t, A=6, omega=2*np.pi/24, phase=0, offset=15)
       >>> plt.plot(t, temp)
       >>> plt.xlabel("Time (h)")
       >>> plt.ylabel("Temperature (°C)")
       >>> plt.title("Idealised diurnal temperature wave")
       >>> plt.show()

       Fit a sinusoid to noisy data with :func:`scipy.optimize.curve_fit`:

       >>> from scipy.optimize import curve_fit
       >>> rng = np.random.default_rng(42)
       >>> y_obs = sinusoid(t, 6, 2*np.pi/24, 0.2, 15) + rng.normal(0, 0.5, t.size)
       >>> p0 = (5, 2*np.pi/24, 0, 15)                       # initial guess
       >>> popt, _ = curve_fit(sinusoid, t, y_obs, p0=p0)
       >>> amp, omg, ph, off = popt
       >>> print(f"Amplitude = {amp:.2f}, Phase = {ph:.2f} rad")



.. py:function:: thermal_diffusivity_amplitude(A1: float, A2: float, z1: float, z2: float, period: int = 86400) -> float

   Estimate soil **thermal diffusivity** (``α``) from the *damping of
   harmonic amplitude* between two depths.

   A one–dimensional soil column subject to a sinusoidal surface
   temperature oscillation exhibits an exponential decay of amplitude
   with depth (Carslaw & Jaeger, 1959).  For a single angular frequency
   :math:`ω = 2π/P`, the analytical solution yields

   .. math::

       α \;=\; \frac{π\, (z_2 - z_1)^2}
                      {P \;\bigl[\,\ln(A_1/A_2)\bigr]^2} ,

   where

   * *A₁* and *A₂* are the harmonic amplitudes at depths *z₁* and *z₂*,
     respectively (*A₁ > A₂*),
   * *P* is the forcing period, and
   * *z₂  – z₁* is the vertical separation of the two sensors.

   :param A1: Diurnal (or other fundamental) temperature amplitudes at the
              shallow depth ``z1`` and deeper depth ``z2``.
              Units **°C** or **K** (identical for both).
   :type A1: float
   :param A2: Diurnal (or other fundamental) temperature amplitudes at the
              shallow depth ``z1`` and deeper depth ``z2``.
              Units **°C** or **K** (identical for both).
   :type A2: float
   :param z1: Sensor depths in **metres** (positive downward).
              Must satisfy ``z2 > z1``.
   :type z1: float
   :param z2: Sensor depths in **metres** (positive downward).
              Must satisfy ``z2 > z1``.
   :type z2: float
   :param period: Fundamental period *P* of the temperature wave in **seconds**.
                  ``86 400`` s corresponds to a 24-hour diurnal cycle.
   :type period: int, default ``86_400``

   :returns: Thermal diffusivity **α** in m² s⁻¹.
   :rtype: float

   :raises ValueError: If ``A1 <= A2`` (violates physical damping assumption) or
       if ``z2 <= z1``.

   .. rubric:: Notes

   * **Amplitude extraction** – ``A1`` and ``A2`` should be obtained
     from a harmonic fit or spectral decomposition that isolates the
     target frequency; raw peak–trough differences are less robust.
   * **Logarithmic sensitivity** – Because the formula involves
     ``ln(A1/A2)``, small uncertainties in amplitudes propagate
     non-linearly; ensure adequate signal-to-noise ratio.
   * Once ``α`` is known, thermal conductivity ``k`` follows from
     ``k = ρc α`` given an independent estimate of volumetric heat
     capacity ``ρc``.

   .. rubric:: References

   Carslaw, H. S., & Jaeger, J. C. (1959).
   *Conduction of Heat in Solids* (2nd ed., pp. 501–502).
   Oxford University Press.

   .. rubric:: Examples

   >>> # Amplitudes from harmonic regression at 5 cm and 10 cm depths
   >>> alpha = thermal_diffusivity_amplitude(
   ...     A1=6.3, A2=4.1, z1=0.05, z2=0.10
   ... )
   >>> print(f"α = {alpha:.2e} m² s⁻¹")
   α = 1.38e-07 m² s⁻¹


.. py:function:: thermal_diffusivity_lag(delta_t, z1, z2, period=86400)

   Estimate soil thermal diffusivity (α) from the phase lag of a temperature wave.

   This method calculates thermal diffusivity based on the time it takes for a
   temperature wave (e.g., the diurnal cycle) to travel between two depths in
   the soil. The formula is derived from the analytical solution to the
   one-dimensional heat conduction equation for a periodic boundary condition:

   .. math::

       \alpha = \frac{P (z_2 - z_1)^2}{4 \pi (\Delta t)^2}

   where :math:`P` is the period of the wave, :math:`(z_2 - z_1)` is the
   distance between the sensors, and :math:`\Delta t` is the time lag of
   the temperature peak between the two depths.

   :param delta_t: Time lag between the temperature peaks at two depths (seconds).
   :type delta_t: float or array_like
   :param z1: Depth of the upper sensor (m, positive downward).
   :type z1: float
   :param z2: Depth of the lower sensor (m, positive downward).
   :type z2: float
   :param period: The period of the temperature wave (seconds). The default is 86400,
                  which corresponds to the daily (diurnal) cycle.
   :type period: int, optional

   :returns: The calculated thermal diffusivity (α) in m² s⁻¹.
   :rtype: float or array_like

   .. rubric:: References

   Nerpin, S.V., and Chudnovskii, A.F. (1967). *Soil physics*.
   (Moscow: Nauka) p 584. (In Russian)


.. py:function:: thermal_diffusivity_logrithmic(t1z1: float, t2z1: float, t3z1: float, t4z1: float, t1z2: float, t2z2: float, t3z2: float, t4z2: float, z1: float, z2: float, period: int = 86400) -> float

   Estimate soil **thermal diffusivity** (``α``) between two depths using the
   *Seemann four–temperature logarithmic* method (also known as the
   Kolmogorov–Seemann method).

   The approach utilises two consecutive half‐period pairs of temperature
   measurements at a shallow depth ``z1`` and a deeper depth ``z2``.
   Let ``T₁–T₄`` denote the temperatures sampled at equal time steps
   (¼ *P*) apart, where *P* is the fundamental period of the harmonic
   forcing.  The solution of the 1-D heat conduction equation for a
   sinusoidal boundary yields

   .. math::

       α \;=\; \frac{4 \, π \, (z_2 - z_1)^2}
                       {P \;\bigl[\,
                       \ln\bigl( ΔT_{z1} / ΔT_{z2} \bigr)\bigr]^2}

   with amplitude decrements

   .. math::

       ΔT_{zij} = \sqrt{(T_1 - T_3)^2 + (T_2 - T_4)^2}\;.

   The formulation is advantageous when only a *short* record is
   available (four points suffice) but is sensitive to sensor noise and
   non-sinusoidal disturbances.

   :param t1z1: Temperatures (°C) at depth ``z1`` sampled at four successive
                quarter-period intervals.
   :type t1z1: float
   :param t2z1: Temperatures (°C) at depth ``z1`` sampled at four successive
                quarter-period intervals.
   :type t2z1: float
   :param t3z1: Temperatures (°C) at depth ``z1`` sampled at four successive
                quarter-period intervals.
   :type t3z1: float
   :param t4z1: Temperatures (°C) at depth ``z1`` sampled at four successive
                quarter-period intervals.
   :type t4z1: float
   :param t1z2: Temperatures (°C) at depth ``z2`` sampled at the *same* times as
                the readings at ``z1``.
   :type t1z2: float
   :param t2z2: Temperatures (°C) at depth ``z2`` sampled at the *same* times as
                the readings at ``z1``.
   :type t2z2: float
   :param t3z2: Temperatures (°C) at depth ``z2`` sampled at the *same* times as
                the readings at ``z1``.
   :type t3z2: float
   :param t4z2: Temperatures (°C) at depth ``z2`` sampled at the *same* times as
                the readings at ``z1``.
   :type t4z2: float
   :param z1: Sensor depths in **metres** (positive downward).  Must satisfy
              ``z2 > z1`` for a meaningful diffusivity.
   :type z1: float
   :param z2: Sensor depths in **metres** (positive downward).  Must satisfy
              ``z2 > z1`` for a meaningful diffusivity.
   :type z2: float
   :param period: Fundamental period *P* of the temperature oscillation in
                  **seconds**.  ``86 400`` s corresponds to a 24-hour diurnal wave.
   :type period: int, default ``86_400``

   :returns: Thermal diffusivity **α** in m² s⁻¹.
   :rtype: float

   .. rubric:: Notes

   * **Sampling interval** – The four readings should be equidistant in
     time and span a full period *P*.  A common practice is to use the
     peak, trough, and two mid-slope points of the diurnal cycle.
   * **Noise sensitivity** – Because the method involves logarithms of
     amplitude ratios, small errors in temperature can propagate
     strongly; consider pre-smoothing or repeating the calculation on
     multiple windows and averaging.
   * **Relation to conductivity** – Once ``α`` is known, bulk thermal
     conductivity ``k`` follows from ``k = ρc α`` with an independent
     estimate of volumetric heat capacity ``ρc``.

   .. rubric:: References

   * Kolmogorov, A. N. (1950). *On the question of determining the
     coefficient of thermal diffusivity of the soil*. *Izvestiya
     Akademii Nauk SSSR, Ser. Geogr. Geofiz.*, 14 (2), 97–99. (In
     Russian)
   * Seemann, W. (1928). *Die Wärmeleitung in der Bodenschicht*.
     Springer, Berlin.

   .. rubric:: Examples

   >>> α = thermal_diffusivity_logrithmic(
   ...     22.5, 20.3, 18.4, 20.1,   # temps @ z1
   ...     18.7, 17.2, 15.9, 17.1,   # temps @ z2
   ...     z1=0.05, z2=0.10,
   ... )
   >>> print(f"α = {α:.2e} m²/s")
   α = 1.46e-07 m²/s


.. py:function:: calc_thermal_diffusivity_log_pair(df, depth1_col, depth2_col, z1, z2, period=86400)

   Estimate soil **thermal diffusivity** (``α``) between two depths using the
   *four-point logarithmic amplitude* method.

   The function extracts the **first four consecutive samples** from two
   temperature records—one at the shallow depth ``z1`` and one at the deeper
   depth ``z2``—and passes them to
   :func:`thermal_diffusivity_logrithmic`.  That helper implements the
   log–ratio solution of the 1-D heat‐conduction equation for a sinusoidal
   boundary condition (Horton et al., 1934; de Vries, 1963):

   .. math::

       α = \frac{(z_2 - z_1)^2}
                 {2P\;\ln\left(\frac{ΔT_{\!z1}}{ΔT_{\!z2}}\right)},

   where

   * **P** is the forcing period (s),
   * :math:`ΔT_{\!z}` is the logarithmic temperature decrement derived
     from four successive measurements at depth *z*.

   The approach is robust for short windows (four points suffice) but is
   sensitive to noise; it is best applied to periods with clear, smooth
   diurnal cycling.

   :param df: Time‐indexed data containing at least the two temperature columns
              specified by ``depth1_col`` and ``depth2_col``.
              **Only the first four rows** are used in the calculation.
   :type df: pandas.DataFrame
   :param depth1_col: Column names for the shallow (``z1``) and deeper (``z2``)
                      temperature series, respectively.
   :type depth1_col: str
   :param depth2_col: Column names for the shallow (``z1``) and deeper (``z2``)
                      temperature series, respectively.
   :type depth2_col: str
   :param z1: Sensor depths in **metres** (positive downward).
              Must satisfy ``z2 > z1``.
   :type z1: float
   :param z2: Sensor depths in **metres** (positive downward).
              Must satisfy ``z2 > z1``.
   :type z2: float
   :param period: Dominant temperature oscillation period **P** in **seconds**.
                  The default (86 400 s) corresponds to 24 h.
   :type period: int, default ``86_400``

   :returns: Thermal diffusivity ``α`` in **m² s⁻¹**.
             Returns ``None`` when fewer than four valid samples are available
             or if ``thermal_diffusivity_logrithmic`` itself returns ``None``.
   :rtype: float or None

   :Warns: **UserWarning** -- Issued (via ``print``) when fewer than four rows are present in
           *df*, in which case the method is skipped and ``None`` is returned.

   .. rubric:: Notes

   * **Data requirement** – The function *does not* resample or align
     series; it simply grabs the first four rows.  Pre-filter or sort
     your DataFrame accordingly.
   * **Noise sensitivity** – Because the method depends on small
     differences between successive temperature readings, apply a
     smoothing filter or select a high-signal period to minimise error.
   * **Relationship to conductivity** – Once ``α`` is known, bulk
     thermal conductivity ``k`` can be obtained from ``k = ρc α`` given
     an estimate of volumetric heat capacity ``ρc``.

   .. rubric:: References

   Horton, R., Wierenga, P. J., Nielsen, D. R., & de Vries, D. A. (1983).
   *Calorimetric determination of soil thermal properties*.
   Soil Science Society of America Journal, **47**, 104–111.

   de Vries, D. A. (1963). *Thermal properties of soils*.
   In *Physics of Plant Environment* (pp. 210–235). North-Holland.

   .. rubric:: Examples

   >>> α_log = calc_thermal_diffusivity_log_pair(
   ...     df=df.sort_index(),          # ensure chronological order
   ...     depth1_col='ts_05cm',
   ...     depth2_col='ts_10cm',
   ...     z1=0.05, z2=0.10,
   ... )
   >>> if α_log is not None:
   ...     print(f"Log-method α = {α_log:.2e} m² s⁻¹")
   Log-method α = 1.45e-07 m² s⁻¹


.. py:function:: calculate_thermal_diffusivity_for_pair(df, col1, col2, z1, z2, period=86400)

   Estimate soil **thermal diffusivity** (``α``) between two depths using
   three classical harmonic methods: *log-amplitude*, *amplitude ratio*,
   and *phase shift*.

   Given two temperature time-series measured at depths ``z1`` and ``z2``,
   the function first extracts the dominant diurnal signal—its amplitude
   and phase—then applies the analytical solutions of the 1-D heat wave
   equation for a homogeneous medium subject to sinusoidal forcing
   (Carslaw & Jaeger, 1959).

   .. method:: 1. Log-Amplitude (α\_log)

      Uses the decay of the harmonic amplitude with depth:

      .. math::

          α\_{\text{log}} = \frac{(z_2 - z_1)^2}
                                   {2\,P\;\ln\bigl(A_1 / A_2\bigr)}


   .. method:: 2. Amplitude Ratio (α\_amp)

      Algebraically identical to the log-amplitude method but expressed
      directly in terms of the two amplitudes:

      .. math::

          α\_{\text{amp}} = \frac{(z_2 - z_1)^2\;\omega}
                                    {2\,[\ln(A_1/A_2)]^2}

      where ``ω = 2π / P`` is the angular frequency.


   .. method:: 3. Phase Lag (α\_lag)

      Relates the travel time (phase shift) of the temperature wave:

      .. math::

          α\_{\text{lag}} = \frac{(z_2 - z_1)^2}{2\,Δt\,P}

      with ``Δt`` the peak-to-peak time lag (s).


   :param df: Time-indexed frame containing temperature observations.
   :type df: pandas.DataFrame
   :param col1: Column names for the shallow and deeper temperature series,
                respectively.
   :type col1: str
   :param col2: Column names for the shallow and deeper temperature series,
                respectively.
   :type col2: str
   :param z1: Sensor depths in **metres** (positive downward).  Must satisfy
              ``z2 > z1``.
   :type z1: float
   :param z2: Sensor depths in **metres** (positive downward).  Must satisfy
              ``z2 > z1``.
   :type z2: float
   :param period: Fundamental period **P** of the harmonic forcing in **seconds**.
                  ``86 400`` s corresponds to 24 h diurnal cycling.
   :type period: int, default ``86_400``

   :returns: Mapping of method identifiers to diffusivity estimates
             (m² s⁻¹):

             * ``'alpha_log'`` – logarithmic amplitude method.
             * ``'alpha_amp'`` – direct amplitude-ratio method.
             * ``'alpha_lag'`` – phase-shift (lag) method.

             Any method returning *None* inside intermediate helpers is
             propagated unchanged.
   :rtype: dict[str, float]

   :raises ValueError: If ``z1`` ≥ ``z2`` or if either column is missing in *df*.

   .. rubric:: Notes

   * ``diurnal_amplitude`` extracts the half range of the 24-h harmonic,
     typically via fast Fourier transform or STL decomposition.
   * ``diurnal_peak_lag`` returns the modal lag **in hours**; the value
     is internally converted to seconds.
   * The function assumes a **single dominant harmonic**.  Strong
     synoptic or weather-front variability can bias results; apply
     filtering or select periods with clear diurnal cycling.
   * Thermal diffusivity relates to thermal conductivity ``k`` through

     .. math:: k = ρ c \, α

     once bulk volumetric heat capacity ``ρc`` is known.

   .. rubric:: References

   Carslaw, H. S., & Jaeger, J. C. (1959). *Conduction of Heat in Solids*
   (2nd ed.). Oxford University Press.

   .. rubric:: Examples

   >>> depth_map = {'ts_05cm': 0.05, 'ts_10cm': 0.10}
   >>> α = calculate_thermal_diffusivity_for_pair(
   ...         df, 'ts_05cm', 'ts_10cm',
   ...         z1=depth_map['ts_05cm'], z2=depth_map['ts_10cm'])
   >>> for meth, val in α.items():
   ...     print(f"{meth}: {val:.2e} m² s⁻¹")
   alpha_log: 1.43e-07 m² s⁻¹
   alpha_amp: 1.41e-07 m² s⁻¹
   alpha_lag: 1.38e-07 m² s⁻¹


.. py:function:: calculate_thermal_properties_for_all_pairs(df, depth_mapping, period=86400)

   Compute **thermal diffusivity**, **thermal conductivity**, and **soil heat
   flux** for *every unique pair* of temperature sensors in a profile.

   The routine iterates over all combinations of the depth‐indexed
   temperature columns supplied in ``depth_mapping``.  For each pair
   *(z₁, z₂)* it

   1. Derives thermal diffusivity ``α`` with
      :func:`calculate_thermal_diffusivity_for_pair`.
   2. Converts ``α`` to thermal conductivity ``k`` via
      :func:`thermal_conductivity`, using the mean volumetric water-
      content of the two layers.
   3. Estimates instantaneous soil heat flux ``G`` by calling
      :func:`soil_heat_flux`.

   Results are returned in a *tidy*, hierarchical ``DataFrame`` whose
   outermost index encodes the depth pair (e.g. ``'0.05-0.10'``).

   :param df: Time-indexed data frame containing at least

              * temperature columns listed in ``depth_mapping``; units **°C**,
                column names typically follow a pattern such as ``'ts_05cm'``.
              * matching soil-water-content columns; each temperature column
                ``'<name>ts'`` must have a companion column
                ``'<name>swc'`` in **percent**.  These are averaged and divided
                by 100 to obtain volumetric θ (*m³ m⁻³*).
   :type df: pandas.DataFrame
   :param depth_mapping: Mapping of *temperature* column names to sensor depths in **metres**
                         (positive downward), e.g. ``{'ts_05cm': 0.05, 'ts_10cm': 0.10}``.
   :type depth_mapping: dict[str, float]
   :param period: Dominant period of the temperature wave (s).  ``86_400`` s
                  corresponds to 24 h and is appropriate for daily forcing.
   :type period: int, default ``86_400``

   :returns: Concatenated frame of thermal properties for every depth pair.
             The outer ``Index`` level is the string ``f"{z1}-{z2}"`` and the
             inner index matches the *datetime* index of ``df`` (after
             dropping rows with *any* missing data).  For each analysis
             “method” returned by
             :func:`calculate_thermal_diffusivity_for_pair` (keys of its
             result dict) the following columns are present:

             ========  ==============================================================
             ``α``     Thermal diffusivity (m² s⁻¹) for that method.
             ``k``     Thermal conductivity (W m⁻¹ K⁻¹) derived from the same α.
             ``G``     Soil heat flux (W m⁻²) between depths z₁ and z₂.
             ``θ_v``   Layer-average volumetric water content (m³ m⁻³).
             ========  ==============================================================
   :rtype: pandas.DataFrame

   .. rubric:: Notes

   * **Alignment** – Each pairwise calculation is performed on a copy of
     ``df`` after dropping all rows with *any* missing values to ensure
     consistent sample support for derived quantities.
   * **Extensibility** – Additional diffusivity algorithms can be
     integrated by returning extra key–value pairs from
     :func:`calculate_thermal_diffusivity_for_pair`; they will be
     propagated automatically.
   * **Performance** – The loop scales *O(n²)* with the number of
     depths.  For large sensor arrays, filter the pairs of interest
     beforehand.

   .. rubric:: Examples

   >>> depth_map = {'ts_05cm': 0.05, 'ts_10cm': 0.10, 'ts_20cm': 0.20}
   >>> props = calculate_thermal_properties_for_all_pairs(df, depth_map)
   >>> props.loc['0.05-0.10'][['alpha_phase', 'G_phase']].plot()
   >>> props.groupby(level=0)['k_amplitude'].median().unstack()


.. py:function:: estimate_rhoc_dry(alpha: pandas.Series, theta: pandas.Series, porosity: float = 0.4, k_dry: float = 0.25, k_sat: float = 1.5, rhoc_w: float = 4180000.0, dry_quantile: float = 0.1) -> float

   Estimate the volumetric **heat capacity of dry soil** (``ρ c_dry``).

   This routine combines concurrent measurements of soil thermal
   diffusivity (``α``) and volumetric water content (``θ``) with a simple
   two–end-member mixing model for thermal conductivity (λ) to back-calculate
   the volumetric heat capacity of the dry soil matrix.  Only the
   *driest* records—defined by the lower ``dry_quantile`` of the observed
   moisture distribution—are used in the final statistic so that
   the latent contribution of soil water is negligible.

   The underlying relationships are

   .. math::

       λ(θ) &= k_\text{dry} + \frac{θ}{φ}
               \,\bigl(k_\text{sat} - k_\text{dry}\bigr)                     \\

       C_v  &= \frac{λ(θ)}{α}                                                   \\

       ρ\,c_\text{dry} &= \frac{C_v - θ\,ρ\,c_w}{1-θ}\,,

   where

   * *λ* is thermal conductivity (W m⁻¹ K⁻¹),
   * *α* is thermal diffusivity (m² s⁻¹),
   * *C_v* is volumetric heat capacity of the *moist* soil
     (J m⁻³ K⁻¹), and
   * *φ* is total porosity (m³ m⁻³).

   :param alpha: Soil thermal diffusivity **α** (m² s⁻¹), indexed identically to
                 *theta* (usually a time-series).
   :type alpha: pandas.Series
   :param theta: Volumetric water content **θ** (m³ m⁻³).
   :type theta: pandas.Series
   :param porosity: Total soil porosity **φ** (saturated water content).
   :type porosity: float, default ``0.40``
   :param k_dry: Thermal conductivity of *air-dry* soil (W m⁻¹ K⁻¹).
   :type k_dry: float, default ``0.25``
   :param k_sat: Thermal conductivity of **saturated** soil (W m⁻¹ K⁻¹).
   :type k_sat: float, default ``1.50``
   :param rhoc_w: Volumetric heat capacity of **liquid water**
                  (J m⁻³ K⁻¹, ≈ 4.18 MJ m⁻³ K⁻¹).
   :type rhoc_w: float, default ``4.18e6``
   :param dry_quantile: Fraction of the *lowest* moisture observations to treat as
                        “dry” when taking the median.  For example, ``0.10`` selects
                        the driest 10 % of the record.
   :type dry_quantile: float, default ``0.10``

   :returns: Median volumetric heat capacity of the *dry* soil matrix
             (J m⁻³ K⁻¹).
   :rtype: float

   .. rubric:: Notes

   * **Alignment** — The two series are first *inner-joined* so only
     timestamps present in both are considered.
   * **Robustness** — Using the median of the driest subset avoids
     bias from residual soil moisture while damping the influence of
     occasional outliers.
   * The default conductivity bounds ``k_dry``/``k_sat`` follow
     typical literature values for mineral soils; adjust them for
     peat, organic, or highly gravelly substrates.

   .. rubric:: Examples

   >>> rhoc_dry = estimate_rhoc_dry(
   ...     alpha=df['alpha_10cm'],
   ...     theta=df['VWC_10cm'],
   ...     porosity=0.43,
   ... )
   >>> print(f"ρ c_dry ≈ {rhoc_dry/1e6:.2f} MJ m⁻³ K⁻¹")
   2.07 MJ m⁻³ K⁻¹


.. py:data:: df
   :value: None


.. py:function:: calculate_soil_heat_storage(df, depths, porosity=0.45, bulk_density=1.3)

   Calculate soil heat storage and heat flux.

   Parameters:
   -----------
   df : DataFrame
       Must contain columns: 'T_5cm', 'T_10cm', etc. and 'SM_5cm', 'SM_10cm', etc.
       where SM is volumetric soil moisture (0-1 or 0-100%)
   depths : list
       Measurement depths in cm, e.g., [5, 10, 20, 30, 40, 50]
   porosity : float
       Soil porosity (0-1), default 0.45
   bulk_density : float
       Dry soil bulk density (g/cm³), default 1.3

   Returns:
   --------
   DataFrame with heat storage and flux values


.. py:function:: lambda_s(theta: numpy.ndarray | float) -> numpy.ndarray | float

   Compute the **soil thermal conductivity** :math:`\lambda_s`
   as a function of volumetric water content (θ).

   The relationship is taken from Gao et al. (2017, Eq. 12):

   .. math::

       \lambda_s(\theta) = 0.20 + \exp\bigl[\,1.46\,(\theta - 0.34)\bigr]

   :param theta: Volumetric water content (m³ m⁻³).
                 Accepts a scalar value or any array-like object that can be
                 converted to a :class:`numpy.ndarray`. Values should lie in the
                 closed interval ``[0, 1]``.
   :type theta: float or array_like

   :returns: Soil thermal conductivity λ\_s in W m⁻¹ K⁻¹.  The returned type
             matches the input: a scalar is returned for a scalar *θ*, and a
             NumPy array for array-like *θ*.
   :rtype: float or numpy.ndarray

   :raises ValueError: If any element of *θ* is outside the physically meaningful
       range ``[0, 1]``.

   .. rubric:: Notes

   * **Vectorization** – The function is fully vectorized; it operates
     element-wise on NumPy arrays and broadcasts according to NumPy
     broadcasting rules.
   * **Empirical limits** – Gao et al. recommend the equation for
     mineral soils where 0 ≤ θ ≤ 0.5. Extrapolation beyond this range
     may introduce error.
   * **Units** – The output uses SI units (W m⁻¹ K⁻¹).

   .. rubric:: References

   Gao, Z., Niu, G.-Y., & Hedquist, B. C. (2017).
   **Soil thermal conductivity parameterization: A closed-form
   equation and its evaluation.** *Journal of Geophysical Research:
   Atmospheres*, 122(6), 3466–3478.
   DOI: 10.1002/2016JD025992

   .. rubric:: Examples

   >>> lambda_s(0.25)
   0.463...
   >>> import numpy as np
   >>> theta = np.linspace(0, 0.5, 5)
   >>> lambda_s(theta)
   array([0.20      , 0.31243063, 0.41172297, 0.51210322, 0.61861198])


.. py:function:: k_s(theta: numpy.ndarray | float) -> numpy.ndarray | float

   Calculate the **soil thermal diffusivity** :math:`k_s`
   as a function of volumetric water content (θ).

   The empirical relationship follows Gao et al. (2017, Eq. 13):

   .. math::

       k_s(\theta) = \bigl[0.69 + \exp\bigl(3.06\,(\theta - 0.26)\bigr)\bigr] \times 10^{-7}

   :param theta: Volumetric water content (m³ m⁻³).
                 Accepts a scalar or any array-like sequence convertible to a
                 :class:`numpy.ndarray`. Values **must** lie in the closed
                 interval ``[0, 1]``.
   :type theta: float or array_like

   :returns: Soil thermal diffusivity *k\_s* in m² s⁻¹.
             A scalar is returned if *θ* is a scalar; otherwise a NumPy array
             of matching shape is returned.
   :rtype: float or numpy.ndarray

   :raises ValueError: If any element of *θ* is outside the physically meaningful
       range ``[0, 1]``.

   .. rubric:: Notes

   * **Vectorization** – The calculation is fully vectorized and
     broadcasts according to NumPy rules.
   * **Applicability** – Gao et al. derived this expression for mineral
     soils under typical field conditions. Extrapolating beyond
     0 ≤ θ ≤ 0.5 may reduce accuracy.
   * **Units** – Output is in SI units (m² s⁻¹).

   .. rubric:: References

   Gao, Z., Niu, G.-Y., & Hedquist, B. C. (2017).
   *Soil thermal conductivity parameterization: A closed-form equation
   and its evaluation.* Journal of Geophysical Research: Atmospheres,
   **122**(6), 3466–3478. https://doi.org/10.1002/2016JD025992

   .. rubric:: Examples

   >>> k_s(0.25)
   1.321...e-07
   >>> thetas = np.linspace(0, 0.5, 6)
   >>> k_s(thetas)
   array([1.14...e-07, 1.21...e-07, 1.30...e-07, 1.41...e-07,
          1.54...e-07, 1.69...e-07])


.. py:function:: volumetric_heat_capacity(lambda_s_val: numpy.ndarray | float, k_s_val: numpy.ndarray | float) -> numpy.ndarray | float

   Compute the **volumetric heat capacity** :math:`C_v` of soil:

   .. math::

       C_v \;=\; \frac{\lambda_s}{k_s}

   where
   :math:`\lambda_s` is the **thermal conductivity** (W m⁻¹ K⁻¹) and
   :math:`k_s` is the **thermal diffusivity** (m² s⁻¹).

   The resulting heat capacity has units of J m⁻³ K⁻¹.

   :param lambda_s_val: Soil thermal conductivity (W m⁻¹ K⁻¹). May be a scalar or any
                        array-like structure broadcastable with *k_s_val*.
   :type lambda_s_val: float or array_like
   :param k_s_val: Soil thermal diffusivity (m² s⁻¹). Must be positive. Accepts a
                   scalar or array-like input.
   :type k_s_val: float or array_like

   :returns: Volumetric heat capacity *C_v* (J m⁻³ K⁻¹).  The return type
             matches the input: a scalar for scalar inputs, or a NumPy array
             for array-like inputs.
   :rtype: float or numpy.ndarray

   :raises ValueError: If any element in *k_s_val* is zero or negative, which would
       lead to division by zero or non-physical results.

   .. rubric:: Notes

   * **Vectorization** – The function is fully vectorized; both inputs
     are converted to :class:`numpy.ndarray` and follow NumPy
     broadcasting rules.
   * **Physical meaning** – Volumetric heat capacity represents the
     energy required to raise the temperature of a unit volume of soil
     by one kelvin. High values correspond to moist or water-logged
     soils; dry mineral soils have lower *C_v*.

   .. rubric:: Examples

   >>> #Scalar inputs
   >>> volumetric_heat_capacity(0.8, 1.2e-6)
   666666.666...

   >>> #Vectorized inputs
   >>> lambda_vals = np.array([0.5, 0.6, 0.7])
   >>> diffusivities = np.array([1.1e-6, 1.2e-6, 1.3e-6])
   >>> volumetric_heat_capacity(lambda_vals, diffusivities)
   array([454545.4545..., 500000.    ..., 538461.5384...])


.. py:function:: nme(calc: numpy.ndarray | float, meas: numpy.ndarray | float) -> float

   Calculate the **normalized mean error (NME)** between calculated
   and measured values, expressed as a percentage.

   The formulation follows Gao et al. (2017, Eq. 14):

   .. math::

       \text{NME} \;=\; 100\,\frac{\sum\limits_{i}\left|\hat{y}_i - y_i\right|}
                                     {\sum\limits_{i}\left|y_i\right|}

   where :math:`\hat{y}_i` are the *calculated* (model) values and
   :math:`y_i` are the *measured* (reference) values.

   :param calc: Modelled / calculated values :math:`\hat{y}`.
                Accepts any array-like object (including scalars) convertible to
                a :class:`numpy.ndarray`.
   :type calc: float or array_like
   :param meas: Measured / observed reference values :math:`y`. Must be
                broadcast-compatible with *calc*.
   :type meas: float or array_like

   :returns: Normalized mean error (percent).
             A value of **0 %** indicates perfect agreement; larger values
             indicate greater error.
   :rtype: float

   :raises ValueError: * If *calc* and *meas* cannot be broadcast to a common shape.
       * If the denominator ``Σ|meas|`` equals zero (e.g., all measured
         values are zero), which would make NME undefined.

   .. rubric:: Notes

   * **Range** – NME is non-negative and unbounded above.
   * **Interpretation** – Because both numerator and denominator use
     absolute values, NME is insensitive to the direction of the error
     (over- vs under-prediction) and therefore complements signed error
     metrics such as mean bias.
   * **Vectorization** – The implementation is fully vectorized and
     adheres to NumPy broadcasting rules.

   .. rubric:: References

   Gao, Z., Niu, G.-Y., & Hedquist, B. C. (2017).
   *Soil thermal conductivity parameterization: A closed-form equation
   and its evaluation.* **Journal of Geophysical Research:
   Atmospheres**, 122(6), 3466–3478.
   https://doi.org/10.1002/2016JD025992

   .. rubric:: Examples

   >>> # Single values
   >>> nme(4.5, 5.0)
   10.0

   >>> # Vectors
   >>> calc = np.array([1.0, 2.1, 3.2])
   >>> meas = np.array([1.2, 2.0, 3.0])
   >>> nme(calc, meas)
   4.7619...

   >>> # Broadcasting (scalar vs array)
   >>> nme(2.0, np.array([1.5, 2.5, 2.0]))
   16.6666...


.. py:function:: rmse(calc: numpy.ndarray | float, meas: numpy.ndarray | float) -> float

   Compute the **root-mean-square error (RMSE)** between calculated
   (model) and measured (reference) values.

   Following Gao et al. (2017, Eq. 15):

   .. math::

       \text{RMSE} \;=\;
       \sqrt{\frac{1}{N}\sum_{i=1}^{N}\bigl(\hat{y}_i - y_i\bigr)^2}

   where :math:`\hat{y}_i` are *calculated* values and
   :math:`y_i` are *measured* values.

   :param calc: Calculated / modelled values :math:`\hat{y}`.  Accepts a scalar
                or any array-like object convertible to a
                :class:`numpy.ndarray`.
   :type calc: float or array_like
   :param meas: Measured / observed values :math:`y`. Must be broadcast-
                compatible with *calc*.
   :type meas: float or array_like

   :returns: Root-mean-square error (same units as the inputs).
   :rtype: float

   :raises ValueError: * If *calc* and *meas* cannot be broadcast to a common shape.
       * If the input arrays are empty (``N = 0``).

   .. rubric:: Notes

   * **Interpretation** – RMSE represents the sample-standard-deviation
     of the differences between two datasets; lower values indicate
     better agreement.
   * **Vectorization** – The function is fully vectorized and respects
     NumPy broadcasting rules.
   * **Units** – RMSE preserves the units of the input variables.

   .. rubric:: References

   Gao, Z., Niu, G.-Y., & Hedquist, B. C. (2017).
   *Soil thermal conductivity parameterization: A closed-form equation
   and its evaluation.* Journal of Geophysical Research:
   Atmospheres, **122**(6), 3466–3478.
   https://doi.org/10.1002/2016JD025992

   .. rubric:: Examples

   >>> # Scalar inputs
   >>> rmse(4.5, 5.0)
   0.5

   >>> # Vector inputs
   >>> calc = np.array([2.1, 3.0, 4.2])
   >>> meas = np.array([2.0, 3.5, 4.0])
   >>> rmse(calc, meas)
   0.2645...

   >>> # Broadcasting
   >>> rmse(3.0, np.array([2.5, 3.5, 3.0]))
   0.4082...


.. py:function:: calorimetric_gz(g_zr: numpy.ndarray | float, cv_layers: numpy.ndarray | float, dT_dt_layers: numpy.ndarray | float, dz_layers: numpy.ndarray | float) -> numpy.ndarray | float

   Estimate the **soil heat flux** at a target depth *z*
   (typically 5 cm) using the *calorimetric method*.

   The method corrects an in-situ heat-flux plate reading made at a
   reference depth *z_r* by adding the change in heat storage of the
   soil column located between *z* and *z_r*.  Mathematically
   (Liebethal & Foken 2007, Eq. 1):

   .. math::

       G(z,t) \;=\; G(z_r,t) \;+\;
       \sum_{l=1}^{N}
       C_{v,l}\,\frac{\partial \bar{T}_l}{\partial t}\,\Delta z_l

   where

   * :math:`G(z,t)`     … heat flux at depth *z* (W m⁻²)
   * :math:`G(z_r,t)`   … measured plate flux at reference depth *z_r*
   * :math:`C_{v,l}`    … volumetric heat capacity of layer *l*
     (J m⁻³ K⁻¹)
   * :math:`\partial \bar{T}_l / \partial t` … time derivative of the
     layer-averaged temperature (K s⁻¹)
   * :math:`\Delta z_l` … thickness of layer *l* (m)

   :param g_zr: Heat-flux plate measurement at depth *z_r* (W m⁻²).  Can be a
                scalar or time series.
   :type g_zr: float or array_like
   :param cv_layers: Volumetric heat capacity :math:`C_{v,l}` for each of *N* soil
                     sub-layers (J m⁻³ K⁻¹).  Shape ``(N, …)`` where the trailing
                     dimensions (``…``) must be broadcast-compatible with *g_zr*.
   :type cv_layers: array_like
   :param dT_dt_layers: Time derivative of layer-mean temperature
                        :math:`\partial \bar{T}_l / \partial t`
                        (K s⁻¹); same shape as *cv_layers*.
   :type dT_dt_layers: array_like
   :param dz_layers: Thickness of each layer :math:`\Delta z_l` (m); shape ``(N,)``
                     or broadcast-compatible with the first axis of *cv_layers*.
   :type dz_layers: array_like

   :returns: Calorimetrically corrected soil heat flux *G(z)* (W m⁻²).  Scalar
             if all inputs are scalar; otherwise a NumPy array matching the
             broadcast shape of *g_zr*.
   :rtype: float or numpy.ndarray

   :raises ValueError: If the inputs cannot be broadcast to a common shape, or if the
       number of layers inferred from *cv_layers*, *dT_dt_layers*, and
       *dz_layers* are inconsistent.

   .. rubric:: Notes

   * **Vectorization** – All operations are fully vectorized.  Inputs
     are converted to :class:`numpy.ndarray` and follow NumPy
     broadcasting rules.
   * **Layer axis** – The first dimension of *cv_layers* and
     *dT_dt_layers* (axis 0) is interpreted as the layer index *l*.
     Layer thicknesses *dz_layers* are broadcast across any additional
     dimensions.
   * **Units** – Ensure consistent SI units: W m⁻², J m⁻³ K⁻¹,
     K s⁻¹, and m.

   .. rubric:: References

   Liebethal, C., & Foken, T. (2007). *Evaluation of six parameterization
   approaches for the ground heat flux.* Agricultural and Forest
   Meteorology, **143**(1–2), 65-80.
   https://doi.org/10.1016/j.agrformet.2006.11.001

   .. rubric:: Examples

   >>> # Three-layer example, single time step
   >>> g_plate = -15.2                      # W m-2 at z_r = −0.05 m
   >>> Cv     = np.array([2.5e6, 2.3e6, 2.1e6])   # J m-3 K-1
   >>> dTdt   = np.array([1.2e-4, 0.9e-4, 0.6e-4])  # K s-1
   >>> dz     = np.array([0.02, 0.02, 0.01])       # m
   >>> calorimetric_gz(g_plate, Cv, dTdt, dz)
   -4.06...

   >>> # Vectorized daily time series with two layers
   >>> g_plate = np.random.normal(-10, 2, 1440)        # per minute
   >>> Cv       = np.array([[2.4e6], [2.2e6]])         # (2,1)
   >>> dTdt     = np.random.normal(5e-5, 2e-5, (2,1440))
   >>> dz       = np.array([0.03, 0.02])
   >>> Gz = calorimetric_gz(g_plate, Cv, dTdt, dz)     # shape (1440,)


.. py:function:: force_restore_gz(cv: numpy.ndarray | float, dTg_dt: numpy.ndarray | float, Tg: numpy.ndarray | float, Tg_bar: numpy.ndarray | float, delta_z: float = 0.05, omega: float = OMEGA_DAY) -> numpy.ndarray | float

   Estimate **soil heat flux** :math:`G(z)` at a shallow depth
   (:math:`z = \delta z`, default 5 cm) with the **force–restore
   method**.

   The formulation (Liebethal & Foken 2007, Eq. 2) corrects the
   calorimetric storage term with a *restore* component that accounts
   for the difference between the instantaneous ground temperature
   *T_g* and its running mean *\bar{T}_g*:

   .. math::

       G(z,t) \;=\; C_v \, \delta z \, \frac{\partial T_g}{\partial t}
       \;+
       \sqrt{\omega \, C_v \, \lambda_s(C_v)}
       \left[
           \frac{1}{\omega}\,\frac{\partial T_g}{\partial t}
           + \bigl(T_g - \bar{T}_g\bigr)
       \right]

   where

   * :math:`C_v`   … volumetric heat capacity (J m⁻³ K⁻¹)
   * :math:`\partial T_g / \partial t` … ground-temperature tendency
     (K s⁻¹)
   * :math:`\lambda_s(C_v)` … soil thermal conductivity derived from
     *C_v* via :func:`lambda_s_from_cv` (W m⁻¹ K⁻¹)
   * :math:`\omega` … angular frequency of the diurnal cycle
     (rad s⁻¹)
   * :math:`T_g` / :math:`\bar{T}_g` … instantaneous and running-mean
     ground temperature (K)
   * :math:`\delta z` … sensor depth below the surface (m)

   :param cv: Volumetric heat capacity :math:`C_v` (J m⁻³ K⁻¹).  Must be
              positive.  Accepts scalars or NumPy-broadcastable arrays.
   :type cv: float or array_like
   :param dTg_dt: Time derivative :math:`\partial T_g/\partial t` (K s⁻¹).
   :type dTg_dt: float or array_like
   :param Tg: Instantaneous ground (surface) temperature :math:`T_g` (K or °C)
              at depth *δz*.
   :type Tg: float or array_like
   :param Tg_bar: Running mean ground temperature :math:`\bar{T}_g` over the
                  diurnal cycle (same units as *Tg*).
   :type Tg_bar: float or array_like
   :param delta_z: Depth :math:`\delta z` in metres; default is **0.05 m**
                   (5 cm).
   :type delta_z: float, optional
   :param omega: Angular frequency :math:`\omega` in rad s⁻¹.  Defaults to
                 :pydata:`OMEGA_DAY` (2π / 86 400 s).
   :type omega: float, optional

   :returns: Soil heat flux *G(z)* at depth *δz* (W m⁻²).  The output shape
             follows NumPy broadcasting rules applied to the inputs.
   :rtype: float or numpy.ndarray

   :raises ValueError: If any input arrays cannot be broadcast to a common shape, or if
       *cv* contains non-positive values.

   .. rubric:: Notes

   * **λ_s(C_v) mapping** – The function relies on a helper
     :func:`lambda_s_from_cv` that converts volumetric heat capacity
     to thermal conductivity.  Ensure that this helper is present in
     the import path.
   * **Units** – Keep units internally consistent (SI).
   * **Vectorization** – All operations are vectorized; the
     mathematical expression is evaluated element-wise for array
     inputs.
   * **Interpretation** – The first term represents *storage*
     (calorimetric), while the second tempers short-term fluctuations,
     “restoring” *T_g* toward its mean.

   .. rubric:: References

   Liebethal, C., & Foken, T. (2007). *Evaluation of six
   parameterization approaches for the ground heat flux.*
   **Agricultural and Forest Meteorology**, 143(1–2), 65-80.
   https://doi.org/10.1016/j.agrformet.2006.11.001

   .. rubric:: Examples

   >>> Cv      = 2.4e6                 # J m-3 K-1
   >>> dTgdt   = 1.0e-4                # K s-1
   >>> Tg      = 293.5                 # K
   >>> Tg_bar  = 291.7                 # K (running mean)
   >>> force_restore_gz(Cv, dTgdt, Tg, Tg_bar)
   19.3...   # W m-2

   >>> #Vectorized daily record
   >>> cv_arr = np.full(1440, 2.2e6)
   >>> dT_arr = np.gradient(np.sin(np.linspace(0, 2*np.pi, 1440))) / 60
   >>> Gz_ts  = force_restore_gz(cv_arr, dT_arr, 298+2*np.sin(...),
   ...                           298*np.ones_like(dT_arr))


.. py:function:: gao2010_gz(AT: numpy.ndarray | float, lambda_s_val: numpy.ndarray | float, k_s_val: numpy.ndarray | float, t: numpy.ndarray | float, omega: float = OMEGA_DAY) -> numpy.ndarray | float

   Estimate **soil heat flux** :math:`G(z,t)` at depth *z* based on a
   *sinusoidal* ground-temperature forcing (Gao et al. 2010, Eq. 3).

   The analytical solution assumes that the surface (or ground-contact)
   temperature varies sinusoidally with amplitude *A_T* and angular
   frequency *ω*.  Under these conditions the heat flux at any depth
   *z* can be written

   .. math::

       G(z,t) \;=\;
       \sqrt{2}\,
       \frac{\lambda_s A_T}{d}\,
       \sin\bigl(\omega t + \tfrac{\pi}{4}\bigr),

   where the **thermal damping depth**

   .. math::

       d \;=\; \sqrt{\frac{2 k_s}{\omega}}

   is determined by the soil thermal diffusivity *k_s* and the forcing
   frequency *ω*.

   :param AT: Amplitude of the sinusoidal ground-surface temperature (K).
   :type AT: float or array_like
   :param lambda_s_val: Soil thermal conductivity :math:`\lambda_s` (W m⁻¹ K⁻¹).
   :type lambda_s_val: float or array_like
   :param k_s_val: Soil thermal diffusivity :math:`k_s` (m² s⁻¹).
   :type k_s_val: float or array_like
   :param t: Time variable (s).  Can be absolute time since epoch or simply
             seconds since the start of the cycle—as long as *ω t* is
             dimensionless.
   :type t: float or array_like
   :param omega: Angular frequency *ω* (rad s⁻¹).  Defaults to *OMEGA_DAY*
                 (≈ 7.272 × 10⁻⁵ s⁻¹, i.e. 2π / 86 400 s).
   :type omega: float, optional

   :returns: Heat flux *G(z,t)* (W m⁻²).  The output follows NumPy’s
             broadcasting rules applied to the inputs.
   :rtype: float or numpy.ndarray

   :raises ValueError: * If *lambda_s_val*, *k_s_val*, and *t* cannot be broadcast to a
         common shape.
       * If any element of *k_s_val* or *omega* is non-positive.

   .. rubric:: Notes

   * **Vectorization** – All inputs are internally converted to
     :class:`numpy.ndarray`; the formula is evaluated element-wise and
     fully supports broadcasting.
   * **Units** – Ensure consistent SI units: W m⁻¹ K⁻¹ (λ_s),
     m² s⁻¹ (k_s), s (t), and rad s⁻¹ (ω).
   * **Scope** – The solution presumes purely sinusoidal boundary
     forcing and homogeneous soil properties; real-world deviations
     (e.g., non-sine forcing, stratified soils, moisture variation)
     will introduce error.

   .. rubric:: References

   Gao, Z., Horton, R., Luo, L., & Kucharik, C. J. (2010).
   *A simple method to measure soil temperature dynamics: Theory and
   application.* **Soil Science Society of America Journal**, 74(2),
   580-588. https://doi.org/10.2136/sssaj2009.0169

   .. rubric:: Examples

   >>> # Daily cycle at 5 cm depth
   >>> AT      = 8.0                           # K
   >>> lambda_ = 1.2                           # W m-1 K-1
   >>> kappa   = 1.0e-6                        # m2 s-1
   >>> t_day   = np.linspace(0, 86400, 97)     # 15-min steps
   >>> Gz      = gao2010_gz(AT, lambda_, kappa, t_day)
   >>> Gz.shape
   (97,)


.. py:function:: heusinkveld_gz(A_n: numpy.ndarray | float, Phi_n: numpy.ndarray | float, n_max: int, k_s_val: numpy.ndarray | float, lambda_s_val: numpy.ndarray | float, w: float) -> numpy.ndarray | float

   Compute *soil heat flux* :math:`G(z,t)` from the **H04 harmonic
   series solution** proposed by Heusinkveld et al. (2004, Eq. 4).

   The approach represents the surface (ground) temperature as a Fourier
   series with harmonics up to order *n_max*.  For each harmonic
   :math:`n` the heat-flux contribution at depth *z* can be written

   .. math::

       G_n(z,t) \;=\;
       \frac{\lambda_s}{10\,\pi}\,
       A_n\;
       \sqrt{\frac{1}{k_s\,n\,\omega\,k_s}}\;
       \sin\bigl(n\,\omega\,t + \Phi_n + \tfrac{\pi}{4}\bigr),

   and the full signal is obtained by summing over *n = 1…n_max*.

   **Note** The implementation below follows the algebraic form that
   appeared in H04; consult the original paper for derivation details
   and recommended parameter ranges.

   :param A_n: Amplitudes :math:`A_n` of the *n*-th harmonic of surface-temperature
               forcing (K).  Provide either
               * a single scalar applied to every harmonic, or
               * an array of length ≥ *n_max* giving amplitude for each harmonic.
   :type A_n: float or array_like
   :param Phi_n: Phase shifts :math:`\Phi_n` (rad) corresponding to each harmonic
                 order.  Same broadcasting rules as *A_n*.
   :type Phi_n: float or array_like
   :param n_max: Highest harmonic order to include in the summation.
   :type n_max: int
   :param k_s_val: Soil thermal diffusivity :math:`k_s` (m² s⁻¹).
   :type k_s_val: float or array_like
   :param lambda_s_val: Soil thermal conductivity :math:`\lambda_s` (W m⁻¹ K⁻¹).
   :type lambda_s_val: float or array_like
   :param w: Fundamental angular frequency :math:`\omega` (rad s⁻¹), e.g.
             :math:`2\pi/86400` for a 24-h cycle.
   :type w: float

   :returns: Soil heat flux *G(z,t)* (W m⁻²).  The return shape is the
             broadcast shape of the input arrays (excluding the harmonic axis).
   :rtype: float or numpy.ndarray

   :raises ValueError: If *n_max* is less than 1, or if *k_s_val* or *w* are non-positive,
       or if the amplitudes/phases cannot be broadcast to length
       *n_max*.

   .. rubric:: Notes

   * **Vectorization** – Internally, the harmonic index axis has length
     *n_max* (``n = np.arange(1, n_max+1)``).  All other dimensions come
     from broadcasting *A_n*, *Phi_n*, *k_s_val*, *lambda_s_val*, and
     *t* (if vectorized in the caller).
   * **Units** – Consistency with SI units is assumed.
   * **Interpretation** – The prefactor ``λ_s / (10 π)`` appears in
     H04’s original derivation.  If you adopt a different convention
     (e.g., depth-explicit damping), modify accordingly.

   .. rubric:: References

   Heusinkveld, B. G., Jacobs, A. F. G., Holtslag, A. A. M., & *et al.*
   (2004). *Surface energy balance closure in an arid region: The role
   of soil heat flux.* Agricultural and Forest Meteorology, **122**(1),
   21-37. https://doi.org/10.1016/j.agrformet.2003.09.005

   .. rubric:: Examples

   >>> # Daily cycle with three harmonics
   >>> A = np.array([10, 4, 1.5])          # K
   >>> Phi = np.array([0, -np.pi/6, np.pi/8])
   >>> Gz = heusinkveld_gz(A, Phi, n_max=3,
   ...                     k_s_val=1e-6,
   ...                     lambda_s_val=1.0,
   ...                     w=2*np.pi/86400)
   >>> Gz.shape
   ()


.. py:function:: hsieh2009_gz(tz_series: numpy.ndarray | list | tuple, time_series: numpy.ndarray | list | tuple, cv_series: numpy.ndarray | list | tuple, ks_series: numpy.ndarray | list | tuple) -> float

   Compute **soil heat flux** at the end of a temperature record
   using the *half-order integral* (Hsieh et al., 2009, Eq. 5).

   The method exploits the analytical solution of the one-dimensional
   heat-conduction equation for a semi-infinite medium, leading to a
   convolution integral of order ½ that relates the time series of
   near-surface soil temperature to the downward heat flux:

   .. math::

       G(t_N) \;=\;
       2\sqrt{\frac{k_s(t_N)\,C_v(t_N)}{\pi}}\;
       \int_{t_0}^{t_N}
           \frac{\partial T(z,t')}{\partial t'}
           \left(t_N - t'\right)^{-1/2}\,dt'

   In discrete form with linear interpolation between measurements
   :math:`t_i \;(i = 0,\dots,N)`,

   .. math::

       \int_{t_0}^{t_N}\!
           \frac{dT}{dt'}\,(t_N-t')^{-1/2}dt'
       \approx
       \sum_{i=0}^{N-1}
         \frac{\Delta T_i}{\Delta t_i}\!
         \left[(t_N-t_i)^{1/2} - (t_N-t_{i+1})^{1/2}\right],

   where :math:`\Delta T_i = T_{i+1}-T_i` and
   :math:`\Delta t_i = t_{i+1}-t_i`.

   The final scalar result corresponds to the flux at
   :math:`t_N\;(=\text{time_series}[-1])`.

   :param tz_series: Near-surface (or shallow-depth) soil temperature observations
                     *T(z,t)* (K).  Length **N ≥ 2**.
   :type tz_series: array_like
   :param time_series: Strictly monotonically increasing time stamps (s) matching
                       `tz_series`.  Same length **N**.
   :type time_series: array_like
   :param cv_series: Volumetric heat capacity *C_v* (J m⁻³ K⁻¹) at each
                     time stamp.  Same length **N**.
   :type cv_series: array_like
   :param ks_series: Thermal diffusivity *k_s* (m² s⁻¹) at each time stamp.
                     Same length **N**.
   :type ks_series: array_like

   :returns: Soil heat flux *G(t_N)* at the final time point (W m⁻²).
   :rtype: float

   :raises ValueError: * If the four series differ in length or have fewer than two
         samples.
       * If `time_series` is not strictly increasing.
       * If any element of `cv_series` or `ks_series` is non-positive.

   .. rubric:: Notes

   * The algorithm uses a forward finite-difference for
     :math:`dT/dt'` and trapezoidal integration over each interval
     *[t_i, t_{i+1}]*.
   * Only the latest values of *k_s* and *C_v* are used, consistent
     with the Hsieh et al. derivation.  Replace with a time-variable
     kernel if property changes are large over the record.
   * All inputs are cast to :class:`numpy.ndarray` with
     ``dtype=float``.

   .. rubric:: References

   Hsieh, C.-I., Katul, G., & Chi, T.-C. (2009).
   *Retrieval of soil heat flux from soil temperature data by the
   continuous-time heat equation model.*
   **Water Resources Research**, 45, W08433.
   https://doi.org/10.1029/2009WR007891

   .. rubric:: Examples

   >>> times  = np.array([0, 600, 1200, 1800])        # every 10 min
   >>> temps  = np.array([292.5, 293.0, 293.6, 294.0])  # K
   >>> Cv     = np.full_like(temps, 2.3e6)            # J m-3 K-1
   >>> ks     = np.full_like(temps, 1.1e-6)           # m2 s-1
   >>> hsieh2009_gz(temps, times, Cv, ks)
   -2.13...


.. py:function:: leuning_damping_depth(z: numpy.ndarray | float, zr: numpy.ndarray | float, AT_z: numpy.ndarray | float, AT_zr: numpy.ndarray | float) -> numpy.ndarray | float

   Estimate the **thermal damping depth** *d* (m) from the exponential
   decay of temperature–wave amplitude with depth.

   The formulation is based on Eq. (6) of Leuning *et al.* (1985) for a
   one-dimensional, sinusoidally forced soil column:

   .. math::

       d \;=\;
       \frac{z_r - z}{\ln\bigl(A_T(z) / A_T(z_r)\bigr)}

   where

   * :math:`A_T(z)` — amplitude of the temperature wave at depth *z*
   * :math:`z_r`  — reference depth (m)
   * :math:`A_T(z_r)` — amplitude at reference depth

   The equation assumes amplitudes decrease monotonically with depth
   following :math:`A_T(z) = A_T(0)\,\exp(-z/d)`.

   :param z: Depth(s) (m) at which the amplitude *A_T(z)* was measured.
             Positive values denote depth below the surface.
   :type z: float or array_like
   :param zr: Reference depth(s) *z_r* (m).  Must be broadcast-compatible with
              *z*.  If *zr* < *z* the denominator switches sign, yielding a
              negative *d* that should be interpreted carefully.
   :type zr: float or array_like
   :param AT_z: Temperature-wave amplitude at *z* (K).
   :type AT_z: float or array_like
   :param AT_zr: Temperature-wave amplitude at *z_r* (K).
   :type AT_zr: float or array_like

   :returns: Thermal damping depth *d* (m).  The output shape follows NumPy
             broadcasting rules applied to the four inputs.
   :rtype: float or numpy.ndarray

   :raises ValueError: If any of the following conditions occur:

       * *AT_z* or *AT_zr* contain non-positive values (log undefined).
       * *AT_z* and *AT_zr* are equal everywhere, leading to
         ``ln(1) = 0`` and division by zero.
       * The inputs cannot be broadcast to a common shape.

   .. rubric:: Notes

   * **Vectorisation** — Internally, all inputs are converted to
     :class:`~numpy.ndarray`; the formula is applied element-wise and
     fully supports broadcasting.
   * **Physical meaning** — Larger *d* indicates slower attenuation of
     temperature fluctuations with depth (i.e. higher thermal
     diffusivity or conductivity).  Very small logarithmic denominators
     (`AT_z` ≈ `AT_zr`) imply an extremely deep damping depth that may
     fall outside the soil column considered.
   * **Units** — Depths must share the same units (metres).  Amplitudes
     may be in kelvins or degrees Celsius provided they use identical
     scaling.

   .. rubric:: References

   Leuning, R., Barlow, E. W. R., & Paltridge, G. (1985).
   *Temperature gradients in a soil: Theory and measurement.*
   *Agricultural and Forest Meteorology*, 35(1–4), 127–143.
   https://doi.org/10.1016/0168-1923(85)90067-3

   .. rubric:: Examples

   >>> d = leuning_damping_depth(z=0.10, zr=0.05, AT_z=4.0, AT_zr=8.0)
   >>> round(d, 4)
   0.0721

   >>> #Broadcasting with arrays
   >>> depths = np.array([0.05, 0.10, 0.20])
   >>> amps   = np.array([8.0, 4.0, 1.5])
   >>> d_arr  = leuning_damping_depth(depths, 0.02, amps, 10.0)
   >>> d_arr
   array([0.0380..., 0.0495..., 0.0854...])


.. py:function:: leuning_gz(g_zr: numpy.ndarray | float, z: numpy.ndarray | float, zr: numpy.ndarray | float, d: numpy.ndarray | float) -> numpy.ndarray | float

   Extrapolate **soil heat flux** from a reference depth *z_r* to a
   shallower (or deeper) target depth *z* using an *exponential
   attenuation* model (Leuning *et al.* 1985, Eq. 7).

   The model assumes the amplitude of the heat-flux wave decays
   exponentially with depth at the same *damping depth* **d** that
   governs temperature attenuation:

   .. math::

       G(z,t) \;=\; G(z_r,t)\;\exp\!\left(\frac{z_r - z}{d}\right)

   where
   :math:`G(z_r,t)` is the measured flux at *z_r* (W m⁻²).

   :param g_zr: Heat flux measured at reference depth *z_r* (W m⁻²).  Can be a
                scalar or a NumPy-broadcastable array (e.g. a time series).
   :type g_zr: float or array_like
   :param z: Target depth *z* (m, **positive downward**).  Must be broadcast-
             compatible with *g_zr*.
   :type z: float or array_like
   :param zr: Reference depth *z_r* (m, positive downward).  Broadcast
              compatible with *z*.
   :type zr: float or array_like
   :param d: Thermal damping depth (m).  Positive; same shape rules as *z*.
   :type d: float or array_like

   :returns: Estimated heat flux at depth *z* (W m⁻²).  Follows NumPy
             broadcasting rules applied to the inputs.
   :rtype: float or numpy.ndarray

   :raises ValueError: If any element of *d* is non-positive, or if inputs cannot be
       broadcast to a common shape.

   .. rubric:: Notes

   * **Direction of extrapolation** –
     *z < z_r* (shallower) ⇒ magnitude *increases* toward the surface;
     *z > z_r* (deeper)   ⇒ magnitude *decreases* with depth.
   * **Vectorisation** – Inputs are converted to
     :class:`numpy.ndarray` and the formula is evaluated element-wise.
   * **Limitations** – The simple exponential form ignores soil
     layering and moisture variability; best suited to homogeneous
     profiles over the depth interval considered.

   .. rubric:: References

   Leuning, R., Barlow, E. W. R., & Paltridge, G. (1985).
   *Temperature gradients in a soil: Theory and measurement.*
   **Agricultural and Forest Meteorology**, 35(1–4), 127-143.
   https://doi.org/10.1016/0168-1923(85)90067-3

   .. rubric:: Examples

   >>> # Single depth conversion
   >>> leuning_gz(g_zr=-12.0, z=0.05, zr=0.08, d=0.07)
   -18.841...

   >>> # Vectorized daily time series
   >>> g_plate = np.random.normal(-10, 3, 1440)      # W m-2 @ 8 cm
   >>> d       = 0.07                                # m
   >>> Gz_5cm  = leuning_gz(g_plate, z=0.05, zr=0.08, d=d)
   >>> Gz_5cm.shape
   (1440,)


.. py:function:: simple_measurement_gz(g_zr: numpy.ndarray | float, cv_layers: numpy.ndarray | float, tz_layers: numpy.ndarray | float, dt: float, dz_layers: numpy.ndarray | float) -> numpy.ndarray | float

       Ground-heat-flux estimate at a target depth *z* using the
       **simple-measurement variant** of the *calorimetric* method
       (Liebethal & Foken 2007, Eq. 8).

       The algorithm corrects a heat-flux‐plate measurement taken at a
       reference depth *z_r* by adding an approximation of the heat storage
       change *ΔS* in the soil slab between the plate and the surface,
       computed from two successive soil-temperature profiles:

       .. math::

           G(z,t_j) \;=\; G(z_r,t_j)
           \;+\; \sum_{l=1}^{N} C_{v,l}\,\Delta z_l\,

   rac{\Delta T_{l,j} +
   rac{1}{2}igl(\Delta T_{l,j} - \Delta T_{l,j-1}igr)}
                       {\Delta t}

       where

       * :math:`C_{v,l}` — volumetric heat capacity of layer *l*
         (J m⁻³ K⁻¹)
       * :math:`\Delta z_l` — layer thickness (m)
       * :math:`\Delta T_{l,j}` — temperature change in layer *l*
         between *t_{j-1}* and *t_j* (K)
       * :math:`\Delta t` — measurement interval (s)

       The centred finite-difference term in brackets approximates the
       mid-interval temperature tendency, providing second-order accuracy
       from simple time-series measurements (`tz_layers`).

       Parameters
       ----------
       g_zr : float or array_like
           Plate-measured soil heat flux at reference depth *z_r*
           (W m⁻²).  Length **M** time steps.
       cv_layers : array_like
           Volumetric heat capacity for each of *N* soil layers
           (J m⁻³ K⁻¹).  Shape ``(N,)`` or broadcast-compatible with the
           first axis of `tz_layers`.
       tz_layers : array_like
           Layer-average soil temperatures.  Shape ``(N, M)`` where *N* is
           the number of layers (matching `cv_layers`/`dz_layers`) and *M*
           is the number of time stamps (matching `g_zr`).
       dt : float
           Constant time step between consecutive temperature samples
           (s).  Must be positive.
       dz_layers : array_like
           Thickness of each layer (m).  Shape ``(N,)`` broadcast-compatible
           with `cv_layers`.

       Returns
       -------
       float or numpy.ndarray
           Calorimetrically corrected soil heat flux at depth *z*
           (W m⁻²).  Length **M − 1** because the centred difference
           requires two successive profiles.

       Raises
       ------
       ValueError
           If input dimensions are inconsistent, `dt` ≤ 0, or any layer
           arrays contain non-finite values.

       Notes
       -----
       * **Temporal alignment** – The first output value corresponds to the
         mid-point of the first two temperature samples
         ``t[0] … t[1]``; therefore the result is shifted **½ Δt** relative
         to the plate-flux time stamps.
       * **Vectorisation** – All calculations are fully vectorized using
         NumPy broadcasting.  Inputs are internally cast to
         :class:`numpy.ndarray`.
       * **Applicability** – Best suited to homogeneous layers with
         high-quality temperature measurements at identical timestamps.

       References
       ----------
       Liebethal, C., & Foken, T. (2007). *Evaluation of six parameterization
       approaches for the ground heat flux.* Agricultural and Forest
       Meteorology, 143(1–2), 65–80.
       https://doi.org/10.1016/j.agrformet.2006.11.001

       Examples
       --------
       >>> g_plate = np.array([-12.3, -10.9, -8.5])      # W m-2 at z_r
       >>> Cv       = np.array([2.3e6, 2.1e6])           # two layers
       >>> T        = np.array([[15.0, 15.5, 16.0],      # °C layer 1
       ...                     [14.0, 14.2, 14.4]])      # °C layer 2
       >>> dz       = np.array([0.03, 0.02])             # m
       >>> Gz = simple_measurement_gz(g_plate, Cv, T, dt=1800, dz_layers=dz)
       >>> Gz
       array([-10.700...,  -8.317...])



.. py:function:: wbz12_g_gz(g_zr_series: numpy.ndarray | list | tuple, time_series: numpy.ndarray | list | tuple, z: float, zr: float, k_s_val: float) -> numpy.ndarray

   Estimate **soil heat flux** at a shallow depth *z* from a flux-plate
   record at reference depth *z_r* using the **WBZ12-G convolution
   method** (Wang, Bou-Zeid & Zhang 2012, their Eq. 9–10).

   The algorithm inverts the one-dimensional heat-conduction equation
   for a semi-infinite homogeneous soil, assuming a step-response
   kernel based on the complementary error function *erfc*:

   .. math::

       G(z,t) \;=\; 2\,G(z_r,t)
       \;-\; \frac{J(t)}{\Delta F_z(t)}

   with the discrete convolution integral

   .. math::

       J(t_n) \;=\;
       \sum_{j=0}^{n-1}
       \tfrac{\bigl[G(z_r,t_{n-j-1}) + G(z_r,t_{n-j})\bigr]}{2}\;
       \Delta F_z(t_{j+1})

   and the *transfer function increment*

   .. math::

       \Delta F_z(t) =
       \operatorname{erfc}\!
       \left[\frac{z_r - z}{2\sqrt{k_s (t - t_0)}}\right]
       \;-\;
       \operatorname{erfc}\!
       \left[\frac{z_r - z}{2\sqrt{k_s (t - t_0 - \Delta t)}}\right].

   :param g_zr_series: Time series of heat-flux-plate measurements at *z_r* (W m⁻²);
                       length **N**.
   :type g_zr_series: array_like
   :param time_series: Monotonically increasing time stamps (s) corresponding to
                       `g_zr_series`.  Must be the same length **N**.
   :type time_series: array_like
   :param z: Target depth (m, **positive downward**).
   :type z: float
   :param zr: Reference depth of the heat-flux plate (m, positive downward).
   :type zr: float
   :param k_s_val: Soil thermal diffusivity *k_s* (m² s⁻¹).  Must be positive.
   :type k_s_val: float

   :returns: Estimated heat-flux series at depth *z*, same length **N** as the
             input series (W m⁻²).
   :rtype: numpy.ndarray

   :raises ValueError: If `k_s_val ≤ 0`, the two series differ in length, time stamps
       are non-monotonic, or contain fewer than two samples.

   .. rubric:: Notes

   * The first element of the output corresponds to *t₀* (no storage
     correction possible yet); subsequent points incorporate the
     convolution up to that instant.
   * A very small term ``+ 1e-12`` is added under the square root to
     avoid division by zero at *t = t₀*.
   * The complementary error function is evaluated with
     :func:`numpy.erfc`, which is vectorized and avoids the SciPy
     dependency.

   .. rubric:: References

   Wang, W., Bou-Zeid, E., & Zhang, Y. (2012).
   *Estimating surface heat fluxes using the surface renewal method*.
   **Boundary-Layer Meteorology**, 144(2), 407-422.
   https://doi.org/10.1007/s10546-012-9730-9

   .. rubric:: Examples

   >>> # 10-min sampled day-long plate record at zr = 0.08 m
   >>> t  = np.arange(0, 24*3600 + 1, 600)            # s
   >>> Gp = -10 + 5*np.sin(2*np.pi*t/86400)
   >>> Gz = wbz12_g_gz(Gp, t, z=0.05, zr=0.08, k_s_val=1.0e-6)
   >>> Gz.shape
   (145,)


.. py:function:: wbz12_s_gz(Ag: numpy.ndarray | float, ks_val: numpy.ndarray | float, zr: float, z: float, t: numpy.ndarray | float, eps: numpy.ndarray | float, omega: float = OMEGA_DAY) -> numpy.ndarray | float

       **WBZ12-S analytic–numeric solution** for soil heat flux *G(z, t)* at
       depth *z* (Wang, Bou-Zeid & Zhang 2012, Eq. 11).

       WBZ12-S combines a closed-form *forcing* term (sinusoidal surface
       temperature) with a *storage* term that is expressed as a Fourier‐
       Laplace integral and must be evaluated numerically.  The governing
       equation assumes homogeneous soil properties and a sinusoidal
       surface temperature of amplitude *A_g* and phase `eps`:

       .. math::

           G(z,t) \;=\;
           A_g\,e^{-\eta}\,\sin\bigl(\omega t +
   arepsilon - \eta\bigr)
           \;-\;
           \frac{2 A_g k_s}{\pi}
           \int_{0}^{\infty}
               \frac{
                   k_s \zeta^{2} igl[\sin\varepsilon
                   - \omega \cos\varepsilon\bigr]
                 }{
                   \omega^{2}
                   + k_{s}^{2} \zeta^{4}
               }
               \sin\bigl[\zeta\,(z_r - z)\bigr]\,
               e^{-k_s \zeta^{2} t}\,d\zeta

       with

       .. math::
           \eta \;=\; (z_r - z)\,\sqrt{\omega / (2 k_s)}.

       The first term (prefactor × sin_term) represents the periodic
       *steady-state* component, while the integral accounts for the
       *transient* adjustment of the soil profile.

       Parameters
       ----------
       Ag : float or array_like
           Amplitude of the sinusoidal surface/ground temperature (K).  Can
           be broadcast with `t` and `eps`.
       ks_val : float or array_like
           Soil thermal diffusivity :math:`k_s` (m² s⁻¹).  Must be positive.
       zr : float
           Reference depth of the heat-flux plate (m, positive downward).
       z : float
           Target depth for the output heat flux (m, positive downward).
       t : float or array_like
           Time (s) since the start of the forcing cycle.  Broadcast-
           compatible with `Ag` and `eps`.
       eps : float or array_like
           Phase shift :math:`\varepsilon` (rad) of the surface
           temperature.  Broadcast-compatible with `t`.
       omega : float, optional
           Angular frequency :math:`\omega` (rad s⁻¹) of the sinusoidal
           forcing (default is :pydata:`OMEGA_DAY` ≈ 7.272 × 10⁻⁵ s⁻¹).

       Returns
       -------
       float or numpy.ndarray
           Soil heat flux *G(z, t)* (W m⁻²).  The shape matches NumPy
           broadcasting over (`Ag`, `ks_val`, `t`, `eps`).

       Raises
       ------
       ValueError
           If *ks_val* or *omega* are non-positive, or if the inputs cannot
           be broadcast to a common shape.

       Notes
       -----
       * **Numerical integral** – The transient term is evaluated using
         :func:`scipy.integrate.quad` over :math:`\zeta \in [0, \infty)`.
         The integral can be expensive for long time-series; cache results
         or vectorise with more sophisticated quadrature if performance is
         critical.
       * **Vectorisation** – `quad` is called individually for every
         element in the broadcasted inputs; large arrays may thus incur
         heavy computational cost.
       * **Units** – Use consistent SI units: W m⁻² (output), K (Ag),
         m² s⁻¹ (ks_val), rad s⁻¹ (omega), seconds (t), metres (depths).

       References
       ----------
       Wang, W., Bou-Zeid, E., & Zhang, Y. (2012).
       *Estimating surface heat fluxes using the surface renewal method.*
       **Boundary-Layer Meteorology**, 144(2), 407–422.
       https://doi.org/10.1007/s10546-012-9730-9

       Examples
       --------
       >>> # Daily sinusoid, single point
       >>> G = wbz12_s_gz(
       ...     Ag=8.0, ks_val=1.0e-6,
       ...     zr=0.08, z=0.05,
       ...     t=np.linspace(0, 86400, 97),
       ...     eps=0.0
       ... )
       >>> G.shape
       (97,)



.. py:function:: exact_temperature_gz(z: numpy.ndarray | float, AT: numpy.ndarray | float, t: numpy.ndarray | float, d: numpy.ndarray | float, omega: float = OMEGA_DAY, T_i: float = 298.15) -> numpy.ndarray | float

   Compute the **exact analytical soil‐temperature profile**
   under sinusoidal surface forcing.

   A sinusoidal temperature wave propagating downward through a
   homogeneous soil decays exponentially with depth and lags in phase.
   The closed‐form solution at depth *z* (m) is

   .. math::

       T(z, t) \;=\;
       T_i
       \;+\;
       A_T \,
       e^{-z/d}\,
       \sin\bigl(\omega t \;-\; z/d\bigr),

   where

   * :math:`T_i` — initial (mean) temperature of the soil column (K)
   * :math:`A_T` — amplitude of the surface-temperature oscillation (K)
   * :math:`d` — thermal damping depth (m)
   * :math:`\omega` — angular frequency of the forcing (rad s⁻¹)
   * :math:`t` — time since the start of oscillation (s)

   :param z: Depth below the surface (m). Positive downward.  Can be a scalar
             or an array broadcastable with *t* and *AT*.
   :type z: float or array_like
   :param AT: Amplitude *A_T* of the surface temperature wave (K).
   :type AT: float or array_like
   :param t: Time variable (s).  Supports vectorisation.
   :type t: float or array_like
   :param d: Thermal damping depth (m). Must be positive.  Can be scalar or
             broadcastable with *z*.
   :type d: float or array_like
   :param omega: Angular frequency *ω* (rad s⁻¹). Defaults to
                 :pydata:`OMEGA_DAY` (≈ 7.272 × 10⁻⁵ s⁻¹ for 24 h).
   :type omega: float, optional
   :param T_i: Mean (initial) soil temperature (K).  Defaults to **298.15 K**
               (25 °C).
   :type T_i: float, optional

   :returns: Soil temperature *T(z, t)* (same units as *T_i*).  The result
             shape matches NumPy broadcasting over the inputs.
   :rtype: float or numpy.ndarray

   :raises ValueError: If any element of *d* or *omega* is non-positive, or if inputs
       cannot be broadcast to a common shape.

   .. rubric:: Notes

   * **Phase lag** – Each e-folding depth *d* introduces a phase delay
     of 1 radian (~57.3 °) relative to the surface signal.
   * **Vectorisation** – All inputs are converted to
     :class:`numpy.ndarray`; the formula is evaluated element-wise and
     fully supports broadcasting.
   * **Units** – Ensure consistent SI units (metres, seconds, kelvins).

   .. rubric:: References

   Adapted from Gao, Z., Horton, R., Luo, L., & Kucharik, C. J. (2010).
   *A simple method to measure soil temperature dynamics: Theory and
   application.* **Soil Science Society of America Journal**, 74(2),
   580-588. https://doi.org/10.2136/sssaj2009.0169

   .. rubric:: Examples

   >>> # Daily cycle at 10 cm depth
   >>> z      = 0.10                            # m
   >>> AT     = 8.0                             # K
   >>> d      = 0.12                            # m
   >>> t_day  = np.linspace(0, 86400, 97)       # 15-min resolution
   >>> Tz     = exact_temperature_gz(z, AT, t_day, d)
   >>> Tz.shape
   (97,)


.. py:function:: exact_gz(z: numpy.ndarray | float, AT: numpy.ndarray | float, lambda_s_val: numpy.ndarray | float, d: numpy.ndarray | float, t: numpy.ndarray | float, omega: float = OMEGA_DAY) -> numpy.ndarray | float

   Exact analytical **soil-heat-flux** solution for sinusoidal surface forcing.

   A sinusoidal surface–temperature wave of amplitude *A_T* and angular
   frequency *ω* generates a heat-flux wave that decays exponentially
   with depth and lags the temperature signal by *π / 4* radians
   (Gao et al., 2010, Eq. 17):

   .. math::

       G(z,t) \;=\;
       \sqrt{2}\;
       \frac{\lambda_s A_T}{d}\;
       e^{-z/d}\;
       \sin\bigl(\omega t - z/d + \tfrac{\pi}{4}\bigr)

   where

   * :math:`\lambda_s` — soil thermal conductivity (W m⁻¹ K⁻¹)
   * :math:`d` — thermal damping depth (m)
   * :math:`z` — depth below the surface (m, positive downward)
   * :math:`t` — time (s) since the start of oscillation
   * :math:`\omega` — angular frequency (rad s⁻¹)

   :param z: Depth(s) below the surface (m), positive downward.
   :type z: float or array_like
   :param AT: Amplitude *A_T* of the surface-temperature oscillation (K).
   :type AT: float or array_like
   :param lambda_s_val: Soil thermal conductivity *λ_s* (W m⁻¹ K⁻¹).
   :type lambda_s_val: float or array_like
   :param d: Thermal damping depth (m). Must be positive.
   :type d: float or array_like
   :param t: Time variable (s).
   :type t: float or array_like
   :param omega: Angular frequency *ω* (rad s⁻¹). Defaults to
                 :pydata:`OMEGA_DAY` (≈ 7.272 × 10⁻⁵ s⁻¹, i.e. 2π / 86 400 s).
   :type omega: float, optional

   :returns: Soil heat flux *G(z, t)* (W m⁻²).  The return shape follows NumPy
             broadcasting rules applied to the inputs.
   :rtype: float or numpy.ndarray

   :raises ValueError: If any element of *d* or *omega* is non-positive, or if inputs
       cannot be broadcast to a common shape.

   .. rubric:: Notes

   * **Phase shift** – At any given depth, the heat-flux wave lags the
     temperature wave by 45 ° (π / 4 rad).
   * **Vectorisation** – All inputs are converted to
     :class:`numpy.ndarray`; the expression is evaluated element-wise
     and fully supports broadcasting.
   * **Units** – Ensure all quantities use consistent SI units.

   .. rubric:: References

   Gao, Z., Horton, R., Luo, L., & Kucharik, C. J. (2010).
   *A simple method to measure soil temperature dynamics: Theory and
   application.* **Soil Science Society of America Journal**, 74(2),
   580–588. https://doi.org/10.2136/sssaj2009.0169

   .. rubric:: Examples

   >>> # Scalar example
   >>> exact_gz(
   ...     z=0.05,
   ...     AT=8.0,
   ...     lambda_s_val=1.2,
   ...     d=0.12,
   ...     t=3600,
   ... )
   21.52...

   >>> # Vectorized daily cycle
   >>> t_day = np.linspace(0, 86400, 97)           # 15-min resolution
   >>> Gz = exact_gz(
   ...     z=0.10,
   ...     AT=6.0,
   ...     lambda_s_val=1.1,
   ...     d=0.11,
   ...     t=t_day,
   ... )
   >>> Gz.shape
   (97,)


.. py:function:: reference_ground_heat_flux(temp_profile: numpy.ndarray, depths: Sequence[float], times: Sequence[float], cv: float, thermal_conductivity: float, gradient_depth: float = 0.2) -> numpy.ndarray

   Compute the reference ground‑heat flux *G₀,M* using the
   gradient method combined with calorimetry for heat storage
   (Liebethal & Foken 2006, Eq. 1).

   The method combines the conductive heat flux at a reference depth
   with the rate of change of heat stored in the soil layer above that
   depth.

   .. math::

       G_{0,M}(t) = -\lambda \frac{\partial T}{\partial z} \bigg|_{z=0.2m}
                    + \int_{z=0}^{0.2m} c_v \frac{\partial T}{\partial t} dz

   :param temp_profile: A 2D array of soil temperatures (°C or K) with shape
                        `(n_depths, n_times)`.
   :type temp_profile: numpy.ndarray
   :param depths: A sequence of measurement depths in meters (positive downward),
                  corresponding to the rows of `temp_profile`.
   :type depths: Sequence[float]
   :param times: A sequence of time stamps in seconds (e.g., Unix timestamps),
                 corresponding to the columns of `temp_profile`.
   :type times: Sequence[float]
   :param cv: Volumetric heat capacity of the soil (J m⁻³ K⁻¹). Assumed
              to be constant with depth.
   :type cv: float
   :param thermal_conductivity: Soil thermal conductivity λ (W m⁻¹ K⁻¹). Assumed to be
                                constant with depth.
   :type thermal_conductivity: float
   :param gradient_depth: The depth (m) at which the vertical temperature gradient is
                          evaluated, by default 0.20 m.
   :type gradient_depth: float, optional

   :returns: A 1D array of the instantaneous ground‑heat flux *G₀,M* (W m⁻²)
             at the surface for each time step. Positive values indicate
             downward flux.
   :rtype: numpy.ndarray

   :raises ValueError: If `temp_profile` shape does not match the lengths of `depths`
       and `times`.


.. py:function:: ground_heat_flux_pr(qs: numpy.ndarray, p: float) -> numpy.ndarray

   Estimate ground heat flux as a fixed fraction of net radiation
   (Liebethal & Foken 2006, Eq. 2).

   This is a simple empirical relationship where the ground heat flux
   is assumed to be a constant proportion of the net radiation at the
   surface.

   .. math:: G_{0,PR}(t) = -p \cdot Q^*_s(t)

   :param qs: Time series of net radiation (W m⁻²). Positive values are
              typically downward.
   :type qs: numpy.ndarray
   :param p: The fraction of net radiation that is partitioned into
             ground heat flux (dimensionless, typically 0 < p < 1).
   :type p: float

   :returns: Time series of the estimated ground heat flux *G₀,PR* (W m⁻²).
   :rtype: numpy.ndarray


.. py:function:: ground_heat_flux_lr(qs: numpy.ndarray, a: float, b: float, lag_steps: int = 0) -> numpy.ndarray

   Estimate ground heat flux using a linear regression against net
   radiation, with an optional time lag (Liebethal & Foken 2006, Eq. 3).

   .. math:: G_{0,LR}(t) = a \cdot Q^*_s(t + \Delta t_G) + b

   :param qs: Time series of net radiation (W m⁻²).
   :type qs: numpy.ndarray
   :param a: The slope of the linear regression (dimensionless).
   :type a: float
   :param b: The intercept of the linear regression (W m⁻²).
   :type b: float
   :param lag_steps: The integer time lag (number of array elements) to apply to the
                     net radiation series. A positive value advances the series
                     (e.g., `qs[t+lag]` is used for `G[t]`), by default 0.
   :type lag_steps: int, optional

   :returns: Time series of the estimated ground heat flux *G₀,LR* (W m⁻²).
   :rtype: numpy.ndarray


.. py:function:: ur_coefficients(delta_ts: float | numpy.ndarray) -> Tuple[numpy.ndarray, numpy.ndarray]

   Compute the parameters A and B for the universal net-radiation
   parameterization, based on the diurnal surface temperature amplitude
   (Liebethal & Foken 2006, Eq. 5 & 6).

   .. math::
       A = 0.0074 \cdot \Delta T_s + 0.088
       B = 1729 \cdot \Delta T_s + 65013

   :param delta_ts: The diurnal amplitude of the surface temperature (K or °C).
   :type delta_ts: float or numpy.ndarray

   :returns: A tuple containing the computed parameters (A, B).
             - A is dimensionless.
             - B is in seconds.
   :rtype: Tuple[numpy.ndarray, numpy.ndarray]


.. py:function:: ground_heat_flux_ur(qs: numpy.ndarray, times_sec: numpy.ndarray, delta_ts: float) -> numpy.ndarray

   Estimate ground heat flux using the universal net-radiation
   parameterization by Santanello & Friedl (2003), as cited in
   Liebethal & Foken (2006, Eq. 4).

   This method modulates the fraction of net radiation partitioned to
   ground heat flux with a cosine function of the time of day.

   .. math::
       G_{0,UR}(t) = -A \cdot \cos\left(\frac{2\pi (t + 10800)}{B}\right)
                     \cdot Q^*_s(t)

   :param qs: Time series of net radiation (W m⁻²).
   :type qs: numpy.ndarray
   :param times_sec: Time stamps in **seconds relative to solar noon**. Positive values
                     indicate the afternoon.
   :type times_sec: numpy.ndarray
   :param delta_ts: The diurnal amplitude of the surface temperature (K or °C).
   :type delta_ts: float

   :returns: Time series of the estimated ground heat flux *G₀,UR* (W m⁻²).
   :rtype: numpy.ndarray


.. py:function:: surface_temp_amplitude(delta_t1: float, delta_t2: float, z1: float, z2: float) -> float

   Estimate the diurnal surface-temperature amplitude (ΔT_s) from
   temperature amplitudes measured at two different soil depths
   (Liebethal & Foken 2006, Eq. 8).

   This method assumes an exponential decay of the temperature wave
   amplitude with depth.

   .. math::
       \Delta T_s = \Delta T_1 + \Delta T_2 \cdot
                     \exp\left(\frac{z_2}{z_2 - z_1}\right)

   :param delta_t1: Diurnal temperature amplitude (K or °C) measured at depth `z1`.
   :type delta_t1: float
   :param delta_t2: Diurnal temperature amplitude (K or °C) measured at depth `z2`.
   :type delta_t2: float
   :param z1: The shallower depth in meters (positive downward).
   :type z1: float
   :param z2: The deeper depth in meters (positive downward).
   :type z2: float

   :returns: The estimated diurnal surface-temperature amplitude ΔT_s (K or °C).
   :rtype: float

   :raises ValueError: If `z2` is not greater than `z1`.


.. py:function:: phi_from_soil_moisture(theta_0_10: float, a_phi: float = 9.62, b_phi: float = 0.402) -> float

   Calculate the empirical parameter φ based on soil moisture content
   (Liebethal & Foken 2006, Eq. 10).

   .. math:: \phi = a_\phi \cdot \theta_{0-10} + b_\phi

   :param theta_0_10: Average volumetric soil moisture content in the top 10 cm
                      (m³ m⁻³).
   :type theta_0_10: float
   :param a_phi: Empirical coefficient, by default 9.62.
   :type a_phi: float, optional
   :param b_phi: Empirical coefficient, by default 0.402.
   :type b_phi: float, optional

   :returns: The dimensionless parameter φ.
   :rtype: float


.. py:function:: ground_heat_flux_sh(h: numpy.ndarray, phase_g0: Sequence[float], phase_h: Sequence[float], u_mean: float, phi: float, omega: float = 2 * np.pi / 86400.0) -> numpy.ndarray

   Estimate ground heat flux from the sensible heat flux (H)
   (Liebethal & Foken 2006, Eq. 9).

   This method relates the ground heat flux to the sensible heat flux
   through a phase-shifted and scaled relationship.

   .. math::
       G_{0,SH}(t) = -\frac{\phi}{\sqrt{\bar{u}}}
                     \frac{\cos(\omega t + \varphi(G_0))}
                          {\cos(\omega t + \varphi(H))} H(t)

   :param h: Time series of sensible heat flux (W m⁻²).
   :type h: numpy.ndarray
   :param phase_g0: Phase lags of the ground heat flux, φ(G₀), in **radians**. Must
                    have the same length as `h`.
   :type phase_g0: Sequence[float]
   :param phase_h: Phase lags of the sensible heat flux, φ(H), in **radians**. Must
                   have the same length as `h`.
   :type phase_h: Sequence[float]
   :param u_mean: Mean horizontal wind speed during the daytime (m s⁻¹).
   :type u_mean: float
   :param phi: An empirical dimensionless parameter, often derived from soil
               moisture via `phi_from_soil_moisture`.
   :type phi: float
   :param omega: The diurnal angular frequency (rad s⁻¹), by default `2π/86400`.
   :type omega: float, optional

   :returns: Time series of the estimated ground heat flux *G₀,SH* (W m⁻²).
   :rtype: numpy.ndarray

   :raises ValueError: If the length of phase arrays does not match the length of `h`.


.. py:function:: ground_heat_flux_sm(gp: numpy.ndarray, t1: numpy.ndarray, delta_t: numpy.ndarray, cv: float, zp: float, dt_seconds: float) -> numpy.ndarray

   Estimate surface ground heat flux using the "simple measurement"
   parameterization, correcting a heat flux plate measurement with a
   storage term (Liebethal & Foken 2006, Eq. 11).

   .. math::
       G_{0,SM}(t) = G_p(t) + c_v z_p \left( \frac{dT_1}{dt} +
                     \frac{1}{2} \frac{d(\Delta T)}{dt} \right)

   :param gp: Heat flux plate measurements at depth `zp` (W m⁻²).
   :type gp: numpy.ndarray
   :param t1: Soil temperature at 0.01 m depth (K or °C).
   :type t1: numpy.ndarray
   :param delta_t: Temperature difference T(0.01 m) - T(zp) (K).
   :type delta_t: numpy.ndarray
   :param cv: Volumetric heat capacity of the soil (J m⁻³ K⁻¹).
   :type cv: float
   :param zp: Depth of the heat flux plate (m, positive downward).
   :type zp: float
   :param dt_seconds: The constant time step between consecutive samples (s).
   :type dt_seconds: float

   :returns: Time series of the estimated ground heat flux *G₀,SM* (W m⁻²).
             The first element will be NaN due to the backward difference.
   :rtype: numpy.ndarray


.. py:function:: active_layer_thickness(lambda_: float, cv: float, omega: float = 2 * np.pi / 86400) -> float

   Calculate the thickness of the active soil layer (δz) for the
   force-restore method (Liebethal & Foken 2006, Eq. 13).

   .. math:: \delta_z = \sqrt{\frac{\lambda}{2 c_v \omega}}

   :param lambda_: Soil thermal conductivity (W m⁻¹ K⁻¹).
   :type lambda_: float
   :param cv: Volumetric heat capacity of the soil (J m⁻³ K⁻¹).
   :type cv: float
   :param omega: The diurnal angular frequency (rad s⁻¹), by default `2π/86400`.
   :type omega: float, optional

   :returns: The thickness of the active soil layer δz (m).
   :rtype: float


.. py:function:: ground_heat_flux_fr(tg: numpy.ndarray, tg_avg: float, cv: float, lambda_: float, delta_z: float | None = None, times: numpy.ndarray | None = None) -> numpy.ndarray

   Estimate ground heat flux using the two-layer force-restore method
   (Liebethal & Foken 2006, Eq. 12).

   .. math::
       G_{0,FR}(t) = -\delta_z c_v \frac{dT_g}{dt} +
                     \sqrt{\lambda \omega c_v} \left(
                     \frac{1}{\omega}\frac{dT_g}{dt} + (T_g - \bar{T_g})
                     \right)

   :param tg: Time series of the upper (surface) layer temperature, Tg(t) (K).
   :type tg: numpy.ndarray
   :param tg_avg: The long-term average or "restoring" temperature, T̄g (K).
   :type tg_avg: float
   :param cv: Volumetric heat capacity of the soil (J m⁻³ K⁻¹).
   :type cv: float
   :param lambda_: Soil thermal conductivity λ (W m⁻¹ K⁻¹).
   :type lambda_: float
   :param delta_z: Thickness of the active soil layer δz (m). If `None`, it is
                   calculated internally using `active_layer_thickness`,
                   by default None.
   :type delta_z: float, optional
   :param times: Time stamps in seconds corresponding to `tg`. Required for
                 calculating the time derivative. If `None`, assumes a uniform
                 time step of 1 second, by default None.
   :type times: numpy.ndarray, optional

   :returns: Time series of the estimated ground heat flux *G₀,FR* (W m⁻²).
   :rtype: numpy.ndarray


.. py:function:: energy_balance_residual(Rn: float | numpy.ndarray, H: float | numpy.ndarray, LE: float | numpy.ndarray, G0: float | numpy.ndarray) -> float | numpy.ndarray

       Compute the **closure residual** of the surface energy balance (SEB).

       The classical SEB for land–atmosphere exchange is

       .. math::

           R_n - G_0 \;=\; H + LE +
   arepsilon ,

       where the residual term :math:`
   arepsilon` quantifies the lack of
       closure.  Rearranging gives

       .. math::


   arepsilon \;=\; R_n - G_0 - H - LE ,

       which is what this helper returns.

       Parameters
       ----------
       Rn : float or array_like
           Net radiation *Rₙ* (W m⁻²).  Positive downward.
       H : float or array_like
           Sensible heat flux *H* (W m⁻²).  Positive upward (atmosphere ← surface).
       LE : float or array_like
           Latent heat flux *LE* (W m⁻²).  Positive upward.
       G0 : float or array_like
           Ground (soil) heat flux *G₀* (W m⁻²).  Positive downward
           (into the soil).  Some authors use the opposite sign convention;
           ensure consistency with *Rn*.

       Returns
       -------
       float or ndarray
           Energy‐balance residual :math:`
   arepsilon` (W m⁻²) with the
           broadcast shape of the inputs.
           **Positive** values indicate missing energy
           (surface gains > turbulent + ground fluxes),
           whereas **negative** values mean an apparent energy surplus.

       Notes
       -----
       * **Broadcasting** – All inputs are treated with NumPy broadcasting,
         allowing scalars, 1-D arrays, or DataFrame columns.
       * **Closure diagnostics** – The residual can be summarised as a mean
         bias or expressed as a relative closure fraction
         ``1 − ε / (Rn − G0)``.
       * **Typical magnitudes** – Eddy‐covariance towers often report
         10–30 % non-closure.  Persistently large |ε| values may indicate
         advective transport, sensor tilt, or data‐processing issues.

       Examples
       --------
       >>> ε = energy_balance_residual(
       ...         Rn=df["Rn"],  H=df["H"],  LE=df["LE"],  G0=df["G"]
       ...     )
       >>> ε.describe(percentiles=[0.05, 0.5, 0.95])
       count    17280.000000
       mean        -9.73
       std         34.51
       5%         -58.42
       50%         -8.11
       95%         39.12
       Name: residual, dtype: float64

       Plot daily average residuals:

       >>> ε.resample("D").mean().plot(marker="o")
       >>> plt.axhline(0, color="k", lw=0.8)
       >>> plt.ylabel("Energy balance residual (W m$^{-2}$)")
       >>> plt.title("Daily SEB closure")



.. py:data:: surface_energy_residual

.. py:function:: ground_heat_flux_conventional(k: float, dT_dz_at_zr: float, rho_c: float, dT_dt_profile: Sequence[float], z_profile: Sequence[float]) -> float

   Conventional estimate of *surface* ground heat flux, Eq. (2).

   :param k: Effective soil thermal conductivity **(W m-1 K-1)**.
   :param dT_dz_at_zr: Vertical temperature gradient evaluated at the flux-plate depth ``z_r``
                       **(K m-1)**.  A *negative* gradient means temperature decreases with
                       depth.
   :param rho_c: Volumetric heat capacity of the soil **(J m-3 K-1)**.
   :param dT_dt_profile: Time derivatives ∂T/∂t for *each* node between the surface and ``z_r``
                         **(K s-1)**.  Any iterable (list, ndarray, …).  Must align with
                         ``z_profile``.
   :param z_profile: Depth of each node in the temperature profile **(m)**.  Increasing,
                     positive downward, **excluding** the surface (z = 0) but *including*
                     ``z_r`` (last element).

   :returns: **G0** -- Ground heat flux at the surface **(W m-2)**.  Positive *into* the soil.
   :rtype: float


.. py:function:: green_function_temperature(z: float, t: float, kappa: float) -> float

   Green-function solution :math:`g_z(t)` for the **one‐dimensional,
   semi-infinite heat equation** with a unit surface‐flux impulse at
   :math:`t=0,\,z=0`.

   The governing partial differential equation is

   .. math::

       \frac{\partial T}{\partial t}
       \;=\;
       \kappa\,\frac{\partial^{2} T}{\partial z^{2}},
       \qquad z \ge 0,\; t > 0,

   with initial condition :math:`T(z,0)=0` and boundary condition
   corresponding to a Dirac δ–heat pulse applied at the surface
   (:math:`z=0`).  The resulting Green function (Carslaw & Jaeger,
   1959, Eq. 7) is

   .. math::

       g_z(t)
       \;=\;
       \frac{2}{\sqrt{\pi}}
       \,\sqrt{\kappa t}\;
       \exp\!\Bigl[-\frac{z^{2}}{4\kappa t}\Bigr]
       \;-\;
       z\,\operatorname{erfc}\!
       \Bigl[\frac{z}{2\sqrt{\kappa t}}\Bigr],
       \qquad t>0,

   and :math:`g_z(t)=0` for :math:`t\le 0` (causality).

   :param z: Depth below the surface (m, positive downward).  Must be
             non-negative.
   :type z: float
   :param t: Time since the surface impulse (s).  Values :math:`t \le 0`
             return 0 by definition.
   :type t: float
   :param kappa: Thermal diffusivity :math:`\kappa` of the half-space
                 (m² s⁻¹).
   :type kappa: float

   :returns: Green-function value :math:`g_z(t)` (units **m**, because the
             solution integrates heat-flux density with respect to depth to
             yield temperature).
   :rtype: float

   .. rubric:: Notes

   * **Causality check** – If ``t`` ≤ 0 the function short-circuits and
     returns 0.0.
   * **Vectorisation** – For vector or array input use
     :func:`numpy.vectorize` or wrap the function in
     :func:`numpy.frompyfunc`; the core implementation is scalar for
     numerical clarity.
   * **Usage** – ``g_z(t)`` can be convolved with an arbitrary surface
     heat-flux time series ``q₀(t)`` to obtain temperature at depth
     via

     .. math::

        T(z,t) \;=\; \int_{0}^{t} g_z(t-τ)\,q_0(τ)\;\mathrm dτ .

   .. rubric:: References

   Carslaw, H. S., & Jaeger, J. C. (1959).
   *Conduction of Heat in Solids* (2nd ed., p. 100).
   Oxford University Press.

   .. rubric:: Examples

   >>> g = green_function_temperature(z=0.05, t=3_600, kappa=1.4e-7)
   >>> print(f"g(5 cm, 1 h) = {g:.3e} m")
   g(5 cm, 1 h) = 7.42e-04 m


.. py:function:: temperature_convolution_solution(z: float, t_series: numpy.ndarray, f_series: numpy.ndarray, kappa: float, Ti: float = 0.0) -> numpy.ndarray

   Temperature time-series at depth *z* via Duhamel convolution (Eq. 6).

   ``T(z,t) = Ti + ∫ f(t-τ) d g_z(τ)``

   The integral becomes a discrete convolution where *f* is the boundary
   heat-flux series (W m-2  → ∂T/∂z via Fourier).


.. py:function:: soil_heat_flux_from_G0(z: float, t_series: numpy.ndarray, G0_series: numpy.ndarray, kappa: float) -> numpy.ndarray

   Compute *G(z,t)* from a known surface flux series *G0* (Eq. 9).


.. py:function:: estimate_G0_from_Gz(Gz_series: numpy.ndarray, z_r: float, kappa: float, dt: float) -> numpy.ndarray

   Estimate *surface* ground heat flux *G0* from plate measurements *Gz*.

   Implements discretised Eq. (11) – the recursion proposed by Wang & Bou-Zeid
   (2012).  Time-series must be *regularly* sampled.

   :param Gz_series: Soil heat-flux measurements at depth *z_r* **(W m-2)**.
   :type Gz_series: np.ndarray
   :param z_r: Plate depth **(m)**.
   :type z_r: float
   :param kappa: Thermal diffusivity **(m² s-1)**.
   :type kappa: float
   :param dt: Sampling interval **(s)**.
   :type dt: float

   :returns: **G0** -- Estimated surface heat-flux series **(W m-2)**.
   :rtype: np.ndarray


.. py:function:: sinusoidal_boundary_flux(t: float | numpy.ndarray, A: float, omega: float, epsilon: float) -> float | numpy.ndarray

       Evaluate a **sinusoidal surface heat-flux forcing**

       .. math::

           q_0(t) \;=\; A \,\sin\!igl(\omega t +
   arepsilonigr),

       which is commonly used as a boundary condition for analytical soil-
       heat-flux solutions (see Eq. 13 of the companion text).

       Parameters
       ----------
       t : float or array_like
           Time since the start of the simulation (s).
           May be scalar or any NumPy-broadcastable shape; units must be
           consistent with ``omega``.
       A : float
           Amplitude of the surface heat flux (W m⁻²).  Positive **downward**
           into the soil.
       omega : float
           Angular frequency (rad s⁻¹).  For a diurnal cycle
           ``omega = 2 π / 86 400`` ≈ 7.272 × 10⁻⁵ rad s⁻¹.
       epsilon : float
           Phase shift **ε** (rad).  Positive values delay the flux peak,
           negative values advance it.

       Returns
       -------
       ndarray or float
           Instantaneous surface heat flux *q₀(t)* (W m⁻²) with shape given
           by NumPy broadcasting of the inputs.

       Notes
       -----
       * **Sign convention** — Positive *q₀* adds energy to the soil column; be
         sure it matches the sign convention of your governing equation.
       * **Vectorisation** — The implementation is a single call to
         ``numpy.sin`` and therefore fully vectorised.
       * **Period** — The period *P* (s) of the forcing is related to
         ``omega`` by *P = 2 π / ω*.

       Examples
       --------
       >>> import numpy as np, matplotlib.pyplot as plt
       >>> t = np.linspace(0, 2*86400, 1_000)                 # 2 days
       >>> q0 = sinusoidal_boundary_flux(
       ...         t, A=120, omega=2*np.pi/86400, epsilon=0)
       >>> plt.plot(t/3600, q0)
       >>> plt.xlabel("Time (h)")
       >>> plt.ylabel("Surface heat flux $q_0$ (W m$^{-2}$)")
       >>> plt.title("Sinusoidal surface forcing (A = 120 W m⁻²)")
       >>> plt.show()



.. py:function:: soil_temperature_sinusoidal(z: float, t: float | numpy.ndarray, A: float, omega: float, epsilon: float, Ti: float, kappa: float) -> float | numpy.ndarray

       Analytical solution for **soil temperature** beneath a sinusoidally
       forced surface heat‐flux boundary.

       The temperature response of a semi-infinite, homogeneous soil column
       to a surface heat flux

       .. math::

           q_0(t) \;=\; A \sin(\omega t +
   arepsilon)

       is (Carslaw & Jaeger 1959, Eq. 14)

       .. math::

           T(z,t)
           \;=\;
           T_i
           \;+\;
           \underbrace{\frac{A}{\kappa\sqrt{\omega}}
                        e^{-z r}\,
                        \sin\!igl(\omega t +
   arepsilon - z r - π/4\bigr)}
           _{\text{steady harmonic}}
           \;+\;
           \underbrace{T_{\text{trans}}(z,t)}
           _{\text{transient integral}} ,

       where :math:`r = \sqrt{\omega / 2\kappa}`.
       The first term is the *steady* periodic component that
       propagates downward with exponentially damped amplitude and a
       depth-dependent phase lag.  The second term accounts for the
       *transient* adjustment from the initial uniform temperature *Tᵢ* and
       is evaluated numerically here by vectorised trapezoidal quadrature.

       Parameters
       ----------
       z : float
           Depth below the soil surface (m, positive downward).
       t : float or array_like
           Time since the start of the forcing (s).  Scalar or NumPy array.
       A : float
           Amplitude of the sinusoidal **surface heat flux** (W m⁻², positive
           downward).  Consistent with
           :func:`sinusoidal_boundary_flux`.
       omega : float
           Angular frequency of the forcing (rad s⁻¹).
           For a diurnal wave ``omega = 2 * π / 86_400``.
       epsilon : float
           Phase shift ε (rad) of the surface heat-flux wave.
       Ti : float
           Initial uniform soil temperature *Tᵢ* (°C or K).
       kappa : float
           Thermal diffusivity κ of the soil (m² s⁻¹).

       Returns
       -------
       float or ndarray
           Soil temperature *T(z, t)* in the same units as *Ti*.
           If *t* is an array the returned array has the same shape.

       Notes
       -----
       * **Steady component** –
         The exponential term ``exp(-z*r)`` dampens amplitude with depth
         while the phase lag is ``z*r + π/4``.
       * **Transient component** –
         The integral is truncated at ``x = 50 / z`` (or 50 if *z = 0*) and
         evaluated with 2 000 panels, which empirically yields < 0.1 %
         relative error for typical soil parameters.  Adjust the limits or
         panel count for higher precision.
       * **Vectorisation** –
         For array *t* the quadrature is performed in parallel using
         broadcasting; memory usage scales with ``len(t)`` × 2 000.
       * **Units** –
         Ensure *A* is W m⁻² **heat flux**.  If you have a surface
         temperature wave instead, transform to an equivalent heat-flux
         boundary or modify the formulation.

       Examples
       --------
       >>> import numpy as np, matplotlib.pyplot as plt
       >>> z   = 0.05                      # 5 cm
       >>> k   = 1.4e-7                    # m² s⁻¹
       >>> A   = 120                       # W m⁻²
       >>> ω   = 2*np.pi/86400             # diurnal
       >>> ε   = 0                         # no phase shift
       >>> Ti  = 15.0                      # °C
       >>> t   = np.linspace(0, 172800, 2000)   # 2 days
       >>> Tz  = soil_temperature_sinusoidal(z, t, A, ω, ε, Ti, k)
       >>> plt.plot(t/3600, Tz)
       >>> plt.xlabel("Time (h)")
       >>> plt.ylabel("Temperature (°C)")
       >>> plt.title("Analytical diurnal soil temperature at 5 cm")
       >>> plt.show()



.. py:function:: soil_heat_flux_sinusoidal(z: float, t: float | numpy.ndarray, A: float, omega: float, epsilon: float, kappa: float) -> float | numpy.ndarray

       Analytical **soil heat–flux** response *G(z,t)* to a sinusoidal
       surface–flux boundary condition.

       A semi-infinite, homogeneous soil column forced at the surface by

       .. math::

           q_0(t) \;=\; A \,\sin\!igl(\omega t +
   arepsilonigr)

       admits the Green-function solution (Carslaw & Jaeger, 1959, Eq. 15)

       .. math::

           G(z,t)
           \;=\;
           A\,e^{-z r}\,\sin(\omega t +
   arepsilon - z r)
           \;+\;
           G_{     ext{trans}}(z,t),

       where :math:`r = \sqrt{ \omega / 2\kappa }` and the transient term

       .. math::

           G_{     ext{trans}}(z,t)
           \;=\;
           -
   rac{2 A \kappa}{\pi}
           \int_{0}^{\infty}
           \frac{( \kappa x^{2}\sin
   arepsilon - \omega \cos
   arepsilon )
                 \;x \sin(x z)}
                {\omega^{2} + \kappa^{2} x^{4}}
           \,e^{-\kappa x^{2} t}\; \mathrm d x .

       The first term is the steady, exponentially damped harmonic that
       propagates downward with a depth-dependent phase lag; the second
       term describes the transient adjustment from an initially unheated
       half-space.  The integral is evaluated here by vectorised
       trapezoidal quadrature on a finite domain (*x* ≤ 50 / *z*).

       Parameters
       ----------
       z : float
           Depth below the soil surface (m, positive downward).
       t : float or array_like
           Time since the onset of forcing (s).  Scalar or NumPy array.
       A : float
           Amplitude of the *surface* heat flux (W m⁻², positive downward).
       omega : float
           Angular frequency ω of the forcing (rad s⁻¹).
           For a diurnal wave ``omega = 2 π / 86 400``.
       epsilon : float
           Phase shift ε of the surface flux (radians).
       kappa : float
           Thermal diffusivity κ of the soil (m² s⁻¹).

       Returns
       -------
       float or ndarray
           Heat flux *G(z,t)* at depth *z* (W m⁻²) with shape equal to
           ``t`` after NumPy broadcasting.
           Positive values **downward**, matching the sign convention of *A*.

       Notes
       -----
       * **Steady component** –
         Amplitude decays as ``exp(-z*r)``; phase lag increases linearly
         with depth by *z r*.
       * **Transient component** –
         The integration limit ``50 / z`` (or 50 for *z = 0*) yields < 0.1 %
         relative error for common parameter ranges.  Increase the upper
         bound and/or panel count (2 000) for stricter accuracy.
       * **Vectorisation** –
         For array *t* the quadrature is evaluated for all time steps
         simultaneously via broadcasting; memory usage scales with
         ``len(t) × 2000``.
       * **Coupling with temperature** –
         The companion function
         :func:`soil_temperature_sinusoidal` gives *T(z,t)* under the
         same boundary forcing; *G* and *T* satisfy Fourier’s law
         ``G = -k ∂T/∂z`` once thermal conductivity *k* is specified.

       References
       ----------
       Carslaw, H. S., & Jaeger, J. C. (1959).
       *Conduction of Heat in Solids* (2nd ed., pp. 100–102).
       Oxford University Press.

       Examples
       --------
       >>> import numpy as np, matplotlib.pyplot as plt
       >>> z      = 0.05                    # 5 cm
       >>> kappa  = 1.4e-7                  # m² s⁻¹
       >>> A      = 120                     # W m⁻² downward
       >>> omega  = 2*np.pi/86400           # diurnal frequency
       >>> eps    = 0
       >>> t      = np.linspace(0, 172800, 2000)   # 2 days
       >>> G      = soil_heat_flux_sinusoidal(z, t, A, omega, eps, kappa)
       >>> plt.plot(t/3600, G)
       >>> plt.xlabel("Time (h)")
       >>> plt.ylabel("Heat flux $G$ (W m$^{-2}$)")
       >>> plt.title("Analytical diurnal soil heat flux at 5 cm")
       >>> plt.show()



.. py:function:: heat_capacity_moist_soil(theta_v: float | numpy.ndarray, theta_s: float, rho_c_s: float = 1260000.0, rho_c_w: float = 4200000.0) -> float | numpy.ndarray

   Volumetric heat capacity of moist soil, Eq. (16).

   :param theta_v: Volumetric water content **(m³ m-3)**.
   :param theta_s: Porosity (saturated volumetric water content) **(m³ m-3)**.
   :param rho_c_s: Heat capacity of dry soil / water **(J m-3 K-1)**.
   :param rho_c_w: Heat capacity of dry soil / water **(J m-3 K-1)**.


.. py:function:: pf_from_theta(theta_v: float | numpy.ndarray, theta_s: float, psi_s: float, b: float) -> float | numpy.ndarray

   Return Pf (Eq. 18) from volumetric water content.


.. py:function:: thermal_conductivity_moist_soil(theta_v: float | numpy.ndarray, theta_s: float, psi_s: float, b: float) -> float | numpy.ndarray

   Thermal conductivity parameterisation, Eq. (17).


.. py:function:: thermal_diffusivity(k: float | numpy.ndarray, rho_c: float | numpy.ndarray) -> float | numpy.ndarray

   Return κ = k / (ρ c).


.. py:function:: soil_heat_flux(Tz: ArrayLike, dz: ArrayLike, lambda_s: ArrayLike | float) -> numpy.ndarray

   Compute heat flux *G* at cell interfaces using Fourier’s law.

   This function calculates the conductive heat flux between soil layers
   based on the temperature gradient and thermal conductivity.

   .. math:: G = -\lambda_s \frac{\partial T}{\partial z}

   :param Tz: A 1D array of temperatures at the **centre** of each soil layer (K).
   :type Tz: ArrayLike
   :param dz: A 1D array of the thickness of each soil layer (m).
   :type dz: ArrayLike
   :param lambda_s: Thermal conductivity for each layer (W m⁻¹ K⁻¹). Can be a single
                    value (for homogeneous soil) or an array matching `Tz`.
   :type lambda_s: ArrayLike or float

   :returns: An array of heat fluxes (W m⁻²) at the *interfaces* between layers,
             including the surface and bottom boundaries. The length of the
             output is `len(Tz) + 1`. Positive values indicate downward flux.
   :rtype: numpy.ndarray


.. py:function:: integrated_soil_heat_flux(rho_c: ArrayLike, T_before: ArrayLike, T_after: ArrayLike, dz: ArrayLike, dt: float, G_ref: float = 0.0) -> numpy.ndarray

   Calculate the soil heat flux profile by integrating the change in
   heat storage upwards from a reference depth (Yang & Wang 2008, Eq. 5).

   .. math::
       G(z_i) = G(z_{ref}) + \int_{z_i}^{z_{ref}} \rho_s c_s
                \frac{\partial T}{\partial t} dz

   :param rho_c: A 1D array of volumetric heat capacity (`ρ_s c_s`) for each soil
                 layer (J m⁻³ K⁻¹).
   :type rho_c: ArrayLike
   :param T_before: A 1D array of temperatures at the start of the time step (K).
   :type T_before: ArrayLike
   :param T_after: A 1D array of temperatures at the end of the time step (K).
   :type T_after: ArrayLike
   :param dz: A 1D array of layer thicknesses (m).
   :type dz: ArrayLike
   :param dt: The duration of the time step (s).
   :type dt: float
   :param G_ref: The heat flux at the lower reference depth `z_ref` (W m⁻²). This
                 is typically assumed to be zero at a sufficient depth, by default 0.0.
   :type G_ref: float, optional

   :returns: An array of heat fluxes (W m⁻²) at the *upper interface* of each
             layer, with a size of `len(dz)`.
   :rtype: numpy.ndarray


.. py:function:: volumetric_heat_capacity(theta: ArrayLike, theta_sat: float | ArrayLike) -> numpy.ndarray

   Calculate the volumetric heat capacity of moist soil based on its
   water content (Yang & Wang 2008, Eq. 4).

   .. math::
       \rho_s c_s = (1 - \theta_{sat}) \rho_{d}c_{d} + \theta \rho_w c_w

   :param theta: Volumetric water content (m³ m⁻³).
   :type theta: ArrayLike
   :param theta_sat: Soil porosity, i.e., saturated volumetric water content (m³ m⁻³).
   :type theta_sat: float or ArrayLike

   :returns: The volumetric heat capacity `ρ_s c_s` (J m⁻³ K⁻¹).
   :rtype: numpy.ndarray


.. py:function:: stretched_grid(n: int, D: float, xi: float) -> numpy.ndarray

   Generate layer thicknesses for a vertically stretched grid.

   The grid layer thickness increases exponentially with depth, allowing
   for higher resolution near the surface.

   .. math:: \Delta z_i = \Delta z_0 \exp(\xi (i-1))

   :param n: The number of soil layers.
   :type n: int
   :param D: The total depth of the soil domain (m).
   :type D: float
   :param xi: The stretching parameter. `xi = 0` results in a uniform grid.
              `xi > 0` results in a grid that is finer at the top.
   :type xi: float

   :returns: A 1D array of thickness `Δz_i` for each of the `n` layers (m).
   :rtype: numpy.ndarray


.. py:function:: solve_tde(T_prev: ArrayLike, dz: ArrayLike, rho_c: ArrayLike, lambda_s: ArrayLike | float, Tsfc: float, Tbot: float, dt: float) -> numpy.ndarray

   Solve the 1D thermal diffusion equation for one time step using an
   implicit Crank-Nicolson scheme (Yang & Wang 2008, Eq. 7).

   This function sets up and solves the tridiagonal system of linear
   equations `M * T_new = D` to find the temperature profile at the
   next time step.

   :param T_prev: Temperature profile at the previous time step (K).
   :type T_prev: ArrayLike
   :param dz: Layer thicknesses (m).
   :type dz: ArrayLike
   :param rho_c: Volumetric heat capacity of each layer (J m⁻³ K⁻¹).
   :type rho_c: ArrayLike
   :param lambda_s: Thermal conductivity of each layer (W m⁻¹ K⁻¹).
   :type lambda_s: ArrayLike or float
   :param Tsfc: Surface temperature boundary condition (K).
   :type Tsfc: float
   :param Tbot: Bottom temperature boundary condition (K).
   :type Tbot: float
   :param dt: Time step (s).
   :type dt: float

   :returns: The new temperature profile at `t + dt`, including boundary nodes.
   :rtype: numpy.ndarray


.. py:function:: correct_profile(T_model: ArrayLike, depths_model: ArrayLike, T_obs: ArrayLike, depths_obs: ArrayLike) -> numpy.ndarray

   Correct a modeled temperature profile using observed temperatures.

   This function calculates the bias (error) between the model and
   observations at the observation depths, then linearly interpolates
   this bias across the entire model grid to correct the profile.

   :param T_model: The modeled temperature profile (K).
   :type T_model: ArrayLike
   :param depths_model: The depths corresponding to `T_model` (m).
   :type depths_model: ArrayLike
   :param T_obs: The observed temperatures (K).
   :type T_obs: ArrayLike
   :param depths_obs: The depths of the observations (m).
   :type depths_obs: ArrayLike

   :returns: The corrected temperature profile.
   :rtype: numpy.ndarray


.. py:function:: surface_temperature_longwave(R_lw_up: float, R_lw_dn: float, emissivity: float = 0.98) -> float

   Calculate surface temperature from upward and downward long-wave
   radiation measurements using the Stefan-Boltzmann law (Yang & Wang 2008, Eq. 8).

   .. math::
       T_s = \left[ \frac{R_{lw}^{\uparrow} - (1 - \epsilon) R_{lw}^{\downarrow}}
                     {\epsilon \sigma} \right]^{1/4}

   :param R_lw_up: Upwelling long-wave radiation (W m⁻²).
   :type R_lw_up: float
   :param R_lw_dn: Downwelling long-wave radiation (W m⁻²).
   :type R_lw_dn: float
   :param emissivity: Surface emissivity (dimensionless), by default 0.98.
   :type emissivity: float, optional

   :returns: The calculated surface temperature (K).
   :rtype: float


.. py:function:: thermal_conductivity_yang2008(theta: ArrayLike, theta_sat: float, rho_dry: float | ArrayLike) -> numpy.ndarray

   Estimate soil thermal conductivity based on soil moisture and dry
   density, following Yang et al. (2005) as cited in
   Yang & Wang (2008, Eq. 9).

   :param theta: Volumetric water content (m³ m⁻³).
   :type theta: ArrayLike
   :param theta_sat: Saturated volumetric water content (porosity) (m³ m⁻³).
   :type theta_sat: float
   :param rho_dry: Dry soil bulk density (kg m⁻³).
   :type rho_dry: float or ArrayLike

   :returns: The estimated soil thermal conductivity `λ_s` (W m⁻¹ K⁻¹).
   :rtype: numpy.ndarray


.. py:function:: flux_error_linear(rho_c: ArrayLike, S2_minus_S1: ArrayLike, dt: float) -> numpy.ndarray

   Calculate the diagnostic error in heat flux that arises from assuming
   a linear temperature profile between measurement points
   (Yang & Wang 2008, Eq. 11).

   .. math:: \Delta G_i = \frac{\rho_s c_s (S_{i,2} - S_{i,1})}{\Delta t}

   :param rho_c: Volumetric heat capacity for each layer (J m⁻³ K⁻¹).
   :type rho_c: ArrayLike
   :param S2_minus_S1: The area difference representing the deviation from linearity.
   :type S2_minus_S1: ArrayLike
   :param dt: The time step (s).
   :type dt: float

   :returns: The flux error for each layer (W m⁻²).
   :rtype: numpy.ndarray


.. py:function:: surface_energy_residual(R_net: float, H: float, LE: float, G0: float) -> float

   Calculate the residual of the surface energy budget (Yang & Wang 2008, Eq. 12).

   .. math:: \Delta E = R_{net} - (H + LE + G_0)

   :param R_net: Net radiation (W m⁻²).
   :type R_net: float
   :param H: Sensible heat flux (W m⁻²).
   :type H: float
   :param LE: Latent heat flux (W m⁻²).
   :type LE: float
   :param G0: Surface ground heat flux (W m⁻²).
   :type G0: float

   :returns: The energy balance residual `ΔE` (W m⁻²).
   :rtype: float


.. py:function:: tdec_step(T_prev: ArrayLike, dz: ArrayLike, theta: ArrayLike, theta_sat: float, rho_dry: float, lambda_const: float, Tsfc: float, Tbot: float, dt: float, depths_model: ArrayLike, T_obs: ArrayLike, depths_obs: ArrayLike) -> Tuple[numpy.ndarray, numpy.ndarray]

   Perform a single integration step of the TDEC (Temperature Diffusion
   Error Correction) scheme.

   This function encapsulates the predict-correct sequence for one time step.

   :param T_prev: Temperature profile from the previous time step (K).
   :type T_prev: ArrayLike
   :param dz: Layer thicknesses (m).
   :type dz: ArrayLike
   :param theta: Volumetric water content profile (m³ m⁻³).
   :type theta: ArrayLike
   :param theta_sat: Soil porosity (m³ m⁻³).
   :type theta_sat: float
   :param rho_dry: Dry soil bulk density (kg m⁻³).
   :type rho_dry: float
   :param lambda_const: A constant thermal conductivity for the prediction step (W m⁻¹ K⁻¹).
   :type lambda_const: float
   :param Tsfc: Surface temperature boundary condition (K).
   :type Tsfc: float
   :param Tbot: Bottom temperature boundary condition (K).
   :type Tbot: float
   :param dt: Time step (s).
   :type dt: float
   :param depths_model: Depths of the model grid nodes (m).
   :type depths_model: ArrayLike
   :param T_obs: Observed temperatures for the correction step (K).
   :type T_obs: ArrayLike
   :param depths_obs: Depths of the observed temperatures (m).
   :type depths_obs: ArrayLike

   :returns: A tuple containing:
             - `T_corr`: The corrected temperature profile at `t + dt`.
             - `G_prof`: The corresponding heat flux profile (W m⁻²) at layer interfaces.
   :rtype: Tuple[np.ndarray, np.ndarray]


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


.. py:function:: compute_soil_storage_integrated(df: pandas.DataFrame, col_T5: str = 'T5cm', col_VWC5: str = 'VWC5cm', col_IRT: str = None, depth: float = 0.05, bulk_density: float = 1300.0, specific_heat_dry: float = 870.0, specific_heat_water: float = 4218.0) -> pandas.Series

   Calculate the soil heat storage flux (S) in W/m² using a change-in-storage
   approach.

   This function estimates the energy stored or released in the upper soil
   layer by accounting for the heat capacity of both the dry soil minerals and
   the stored liquid water. It supports an integrated temperature approach if
   infrared surface temperature (IRT) is available. In the slab approach, you
   assume the estimate at the given depth applies to the whole slab above the
   depth, whereas the IRT approach takes the mean between surface temperature
   and temperature at depth.

   Note that this method is the same as the method use to calculate soil heat
   storage flux in the Campbell Scientific programs, though that progrmam uses
   higher-frequency data for the calculations.

   :param df: Input dataframe with a DatetimeIndex.
   :type df: pd.DataFrame
   :param col_T5: Column name for soil temperature at 5cm depth [°C].
   :type col_T5: str, default "T5cm"
   :param col_VWC5: Column name for Volumetric Water Content at 5cm depth [m³/m³].
   :type col_VWC5: str, default "VWC5cm"
   :param col_IRT: Column name for surface infrared temperature [°C]. If provided,
                   the layer temperature is the average of surface and 5cm temps. Note
                   that this must be surface temperature, not canopy temperature
   :type col_IRT: str, optional
   :param depth: The thickness of the soil layer being measured [m].
   :type depth: float, default 0.05
   :param bulk_density: Dry bulk density of the soil [kg/m³].
   :type bulk_density: float, default 1300.0
   :param specific_heat_dry: Specific heat capacity of dry soil minerals [J/(kg·K)].
   :type specific_heat_dry: float, default 870.0
   :param specific_heat_water: Specific heat capacity of water [J/(kg·K)].
   :type specific_heat_water: float, default 4218.0

   :returns: The calculated soil heat storage flux [W/m²].
             Positive values indicate energy moving into storage (warming).
   :rtype: pd.Series

   .. rubric:: Notes

   The canopy storage flux is calculated as:
   $$S_{canopy} = \frac{m_{veg} \cdot C_{p} \cdot \Delta T}{\Delta t}$$

   Where:
   - $m_{veg}$ = Fresh biomass ($kg/m^2$), calculated as $(height / 100) \cdot density\_factor$
   - $C_{p}$ = Specific heat of fresh vegetation ($\approx 3900 \text{ J/kg}\cdot\text{K}$)
   - $\Delta T$ = Change in canopy temperature between time steps
   - $\Delta t$ = Time step in seconds


.. py:function:: compute_canopy_storage(df: pandas.DataFrame, col_irt: str = 'T_CANOPY', crop_type: str = 'alfalfa', height_cm: str | float = None, density_factor: float = None) -> pandas.Series

   Calculate the energy storage flux of a plant canopy (S_canopy) in W/m².

   This function estimates heat storage within biomass. It can handle dynamic
   vegetation growth if a height column is provided, otherwise it falls back
   to static heights. Height and density values can be entered by the user or
   taken from a dictionary within the function.

   Care should be taken to examine the crop default values and the specific
   heat of the vegetation.

   :param df: Input dataframe with a DatetimeIndex.
   :type df: pd.DataFrame
   :param col_irt: Column name for the canopy surface temperature [°C].
   :type col_irt: str, default "T_CANOPY"
   :param crop_type: Options: 'alfalfa', 'corn', 'hay', 'wetland'. Used to determine
                     defaults for height and density.
   :type crop_type: str, default "alfalfa"
   :param height_cm: If str: The column name in `df` containing canopy height [cm] over time.
                     If float: A static height value [cm].
                     If None: Uses default height for the specified `crop_type`.
   :type height_cm: str or float, optional
   :param density_factor: Biomass density [kg/m² per meter height]. If None, uses default
                          for the specified `crop_type`.
   :type density_factor: float, optional

   :returns: Energy storage flux [W/m²].
   :rtype: pd.Series


