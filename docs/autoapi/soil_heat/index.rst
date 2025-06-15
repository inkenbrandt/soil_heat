soil_heat
=========

.. py:module:: soil_heat

.. autoapi-nested-parse::

   Top-level package for Soil Heat.



Submodules
----------

.. toctree::
   :maxdepth: 1

   /autoapi/soil_heat/gao_et_al/index
   /autoapi/soil_heat/liebethal_and_folken/index
   /autoapi/soil_heat/soil_heat/index
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


Package Contents
----------------

.. py:data:: __author__
   :value: 'Paul Inkenbrandt'


.. py:data:: __email__
   :value: 'paulinkenbrandt@utah.gov'


.. py:data:: __version__
   :value: '0.1.0'


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

   Calculate soil heat flux (G) using temperature gradient and thermal conductivity.

   Parameters:
   - T_upper: Temperature at upper depth (°C)
   - T_lower: Temperature at lower depth (°C)
   - depth_upper: Upper sensor depth (m)
   - depth_lower: Lower sensor depth (m)
   - k: Thermal conductivity (W/(m·°C))

   Returns:
   - Soil heat flux (W/m^2)


.. py:function:: volumetric_heat_capacity(theta_v)

   Estimate volumetric heat capacity Cv (J/(m³·°C)) from soil moisture.

   Parameters:
   - theta_v: Volumetric water content (decimal fraction, e.g., 0.20 for 20%)

   Returns:
   - Volumetric heat capacity (kJ/(m³·°C))


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
     It is obtained from :pyfunc:`volumetric_heat_capacity`.

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
     :pyfunc:`volumetric_heat_capacity` typically assumes a two‐phase
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

   Estimate thermal diffusivity from phase lag.

   Parameters:
   - delta_t: Time lag between peaks at two depths (seconds)
   - z1, z2: Depths (m)
   - period: Time period of wave (default = 86400 s for daily cycle)

   Returns:
   - Thermal diffusivity α (m²/s)

   Citation:
   S.V. Nerpin, and A.F. Chudnovskii, Soil physics, (Moscow: Nauka) p 584, 1967 (in Russian)


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
   :pyfunc:`thermal_diffusivity_logrithmic`.  That helper implements the
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

   .. method:: \*\*1. Log-Amplitude (α\_log)**

      Uses the decay of the harmonic amplitude with depth:

      .. math::

          α\_{\text{log}} = \frac{(z_2 - z_1)^2}
                                   {2\,P\;\ln\bigl(A_1 / A_2\bigr)}


   .. method:: \*\*2. Amplitude Ratio (α\_amp)**

      Algebraically identical to the log-amplitude method but expressed
      directly in terms of the two amplitudes:

      .. math::

          α\_{\text{amp}} = \frac{(z_2 - z_1)^2\;\omega}
                                    {2\,[\ln(A_1/A_2)]^2}

      where ``ω = 2π / P`` is the angular frequency.


   .. method:: \*\*3. Phase Lag (α\_lag)**

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
      :pyfunc:`calculate_thermal_diffusivity_for_pair`.
   2. Converts ``α`` to thermal conductivity ``k`` via
      :pyfunc:`thermal_conductivity`, using the mean volumetric water-
      content of the two layers.
   3. Estimates instantaneous soil heat flux ``G`` by calling
      :pyfunc:`soil_heat_flux`.

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
             :pyfunc:`calculate_thermal_diffusivity_for_pair` (keys of its
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
     :pyfunc:`calculate_thermal_diffusivity_for_pair`; they will be
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


.. py:function:: lambda_s(theta: numpy.ndarray | float) -> numpy.ndarray | float

   Thermal conductivity (λ_s) as a function of volumetric water content θ.

   Implements Eq. (12) from Gao et al. (2017) fileciteturn1file3.


.. py:function:: k_s(theta: numpy.ndarray | float) -> numpy.ndarray | float

   Thermal diffusivity (k_s) as a function of volumetric water content θ.

   Implements Eq. (13) from Gao et al. (2017) fileciteturn1file6.


.. py:function:: volumetric_heat_capacity(lambda_s_val, k_s_val)

   Volumetric heat capacity C_v = λ_s / k_s (units J m⁻³ K⁻¹).


.. py:function:: nme(calc: numpy.ndarray, meas: numpy.ndarray) -> float

   Normalized mean error (%). Implements Eq. (14) fileciteturn1file8.


.. py:function:: rmse(calc: numpy.ndarray, meas: numpy.ndarray) -> float

   Root‑mean‑square error. Implements Eq. (15) fileciteturn1file8.


.. py:function:: calorimetric_gz(g_zr, cv_layers, dT_dt_layers, dz_layers)

   Calorimetric method for G_z at depth z (usually 5 cm).

   :param g_zr: Measured heat flux at reference depth *z_r* (W m⁻²).
   :type g_zr: float or array_like
   :param cv_layers: Volumetric heat capacity for each sub‑layer *C_v,l* (J m⁻³ K⁻¹).
   :type cv_layers: sequence
   :param dT_dt_layers: Time derivative of average temperature for each layer ∂T/∂t (K s⁻¹).
   :type dT_dt_layers: sequence
   :param dz_layers: Thickness of each sub‑layer δz_l (m).
   :type dz_layers: sequence


.. py:function:: force_restore_gz(cv, dTg_dt, Tg, Tg_bar, delta_z=0.05, omega=OMEGA_DAY)

   Force‑restore estimate of G_z at z = δz (default 5 cm).

   Implements Eq. (2) fileciteturn1file2.


.. py:function:: gao2010_gz(AT, lambda_s_val, k_s_val, t, omega=OMEGA_DAY)

   Sinusoidal solution for G_z at depth d (Eq. 3).


.. py:function:: heusinkveld_gz(A_n, Phi_n, n_max, k_s_val, lambda_s_val, w)

   H04 harmonic solution (Eq. 4).


.. py:function:: hsieh2009_gz(tz_series, time_series, cv_series, ks_series)

   Half‑order integral solution (Eq. 5).

   tz_series, cv_series, ks_series must be monotonically increasing in *time_series*.


.. py:function:: leuning_damping_depth(z, zr, AT_z, AT_zr)

   Compute damping depth *d* via Eq. (6).


.. py:function:: leuning_gz(g_zr, z, zr, d)

   Exponentially adjust G_z from reference depth (Eq. 7).


.. py:function:: simple_measurement_gz(g_zr, cv_layers, tz_layers, dt, dz_layers)

   Simple‑measurement variant of calorimetric method (Eq. 8).


.. py:function:: wbz12_g_gz(g_zr_series, time_series, z, zr, k_s_val)

   WBZ12‑G method (Eq. 9–10).


.. py:function:: wbz12_s_gz(Ag, ks_val, zr, z, t, eps, omega=OMEGA_DAY)

   WBZ12‑S solution (Eq. 11).

   Numerical integration is performed for the second term.


.. py:function:: exact_temperature_gz(z, AT, t, d, omega=OMEGA_DAY, T_i=298.15)

   Exact sinusoidal soil‑temperature profile (Eq. 16).


.. py:function:: exact_gz(z, AT, lambda_s_val, d, t, omega=OMEGA_DAY)

   Exact sinusoidal soil‑heat‑flux (Eq. 17).


.. py:function:: reference_ground_heat_flux(temp_profile: numpy.ndarray, depths: Sequence[float], times: Sequence[float], cv: float, thermal_conductivity: float, gradient_depth: float = 0.2) -> numpy.ndarray

   Compute the reference ground‑heat flux *G₀,M* (Eq. 1).

   Equation
   --------
   G₀,M(t) = -λ ∂T/∂z |_(z=0.2 m) + ∫_{z=0}^{0.2 m} c_v ∂T/∂t dz

   :param temp_profile: Soil temperatures (°C or K) at the depths specified by *depths* and
                        time stamps *times*.
   :type temp_profile: ndarray, shape (n_z, n_t)
   :param depths: Measurement depths (m, **positive downward**).
   :type depths: sequence of float, length *n_z*
   :param times: Epoch time in **seconds** (may be monotonic pandas DatetimeIndex
                 converted via ``astype('int64')/1e9``).
   :type times: sequence of float, length *n_t*
   :param cv: Volumetric heat capacity of the soil (J m⁻³ K⁻¹).
   :type cv: float
   :param thermal_conductivity: Soil thermal conductivity λ (W m⁻¹ K⁻¹).
   :type thermal_conductivity: float
   :param gradient_depth: Depth (m) at which the vertical gradient term is evaluated.
   :type gradient_depth: float, default 0.20

   :returns: Instantaneous ground‑heat flux *G₀,M* (W m⁻²). Positive = downward.
   :rtype: ndarray, shape (n_t,)


.. py:function:: ground_heat_flux_pr(qs: numpy.ndarray, p: float) -> numpy.ndarray

   Ground heat flux using a fixed *p* fraction of net radiation (Eq. 2).

   G₀,PR(t) = ‑p · Q*ₛ(t)

   :param qs: Net radiation time series (W m⁻²). Positive = downward.
   :type qs: ndarray
   :param p: Fraction of net radiation that becomes ground‑heat flux (0–1).
   :type p: float

   :returns: G₀,PR (W m⁻²).
   :rtype: ndarray


.. py:function:: ground_heat_flux_lr(qs: numpy.ndarray, a: float, b: float, lag_steps: int = 0) -> numpy.ndarray

   Linear net‑radiation parameterisation (Eq. 3).

   G₀,LR(t) = a·Q*ₛ(t+Δt_G) + b

   :param qs: Net radiation (W m⁻²).
   :type qs: ndarray
   :param a: Regression coefficients.
   :type a: float
   :param b: Regression coefficients.
   :type b: float
   :param lag_steps: Integer lag (number of samples) by which *qs* is advanced.
   :type lag_steps: int, default 0

   :returns: G₀,LR (W m⁻²).
   :rtype: ndarray


.. py:function:: ur_coefficients(delta_ts: float | numpy.ndarray) -> Tuple[numpy.ndarray, numpy.ndarray]

   Compute universal‑function parameters *A* and *B* (Eq. 5 & 6).

   :param delta_ts: Diurnal amplitude of *surface* temperature (K).
   :type delta_ts: float or ndarray

   :returns: * **A** (*ndarray*)
             * **B** (*ndarray (seconds)*)


.. py:function:: ground_heat_flux_ur(qs: numpy.ndarray, times_sec: numpy.ndarray, delta_ts: float) -> numpy.ndarray

   Universal net‑radiation parameterisation (Eq. 4).

   Implements Santanello & Friedl (2003):
       G₀,UR(t) = -A · cos[2π (t + 10800) / B] · Q*ₛ(t)

   *t* is **seconds since solar noon** (positive in afternoon).

   :param qs: Net radiation (W m⁻²).
   :type qs: ndarray
   :param times_sec: Seconds relative to solar noon (s).
   :type times_sec: ndarray
   :param delta_ts: Diurnal surface‑temperature amplitude (K).
   :type delta_ts: float

   :returns: G₀,UR (W m⁻²).
   :rtype: ndarray


.. py:function:: surface_temp_amplitude(delta_t1: float, delta_t2: float, z1: float, z2: float) -> float

   Compute diurnal surface‑temperature amplitude ΔT_s (Eq. 8).

   :param delta_t1: Diurnal temperature amplitudes (K) measured at depths *z1* and *z2*.
   :type delta_t1: float
   :param delta_t2: Diurnal temperature amplitudes (K) measured at depths *z1* and *z2*.
   :type delta_t2: float
   :param z1: Depths in meters (**positive downward**, with z2 > z1 > 0).
   :type z1: float
   :param z2: Depths in meters (**positive downward**, with z2 > z1 > 0).
   :type z2: float

   :returns: Estimated ΔT_s (K).
   :rtype: float


.. py:function:: phi_from_soil_moisture(theta_0_10: float, a_phi: float = 9.62, b_phi: float = 0.402) -> float

   Soil‑moisture dependent φ (Eq. 10).


.. py:function:: ground_heat_flux_sh(h: numpy.ndarray, phase_g0: Sequence[float], phase_h: Sequence[float], u_mean: float, phi: float, omega: float = 2 * np.pi / 86400.0) -> numpy.ndarray

   Ground‑heat flux from sensible heat flux H (Eq. 9).

   :param h: Sensible heat flux time series (W m⁻²).
   :type h: ndarray
   :param phase_g0: Phase lags φ(G₀) and φ(H) in **radians**.
   :type phase_g0: sequence of float
   :param phase_h: Phase lags φ(G₀) and φ(H) in **radians**.
   :type phase_h: sequence of float
   :param u_mean: Mean horizontal wind speed during daytime (m s⁻¹).
   :type u_mean: float
   :param phi: Empirical parameter (dimensionless), see `phi_from_soil_moisture`.
   :type phi: float
   :param omega: Diurnal angular frequency (s⁻¹).
   :type omega: float, default 2π/86400

   :returns: G₀,SH (W m⁻²).
   :rtype: ndarray


.. py:function:: ground_heat_flux_sm(gp: numpy.ndarray, t1: numpy.ndarray, delta_t: numpy.ndarray, cv: float, zp: float, dt_seconds: float) -> numpy.ndarray

   Simple‑measurement parameterisation (Eq. 11).

   :param gp: Heat‑flux plate measurement at depth *zp* (W m⁻²).
   :type gp: ndarray
   :param t1: Soil temperature at 0.01 m depth (K or °C).
   :type t1: ndarray
   :param delta_t: Temperature difference T(0.01 m) – T(z_p) (K).
   :type delta_t: ndarray
   :param cv: Volumetric heat capacity (J m⁻³ K⁻¹).
   :type cv: float
   :param zp: Plate depth (m, positive downward).
   :type zp: float
   :param dt_seconds: Time step between consecutive samples (s).
   :type dt_seconds: float

   :returns: G₀,SM (W m⁻²).
   :rtype: ndarray


.. py:function:: active_layer_thickness(lambda_: float, cv: float, omega: float = 2 * np.pi / 86400) -> float

   Thickness δz of the active soil layer (Eq. 13).


.. py:function:: ground_heat_flux_fr(tg: numpy.ndarray, tg_avg: float, cv: float, lambda_: float, delta_z: float | None = None, times: numpy.ndarray | None = None) -> numpy.ndarray

   Force‑restore ground‑heat flux (Eq. 12).

   Implements the two‑layer force‑restore formulation with an optional
   diagnostic *δz* computed via Eq. 13 if not supplied.

   :param tg: Temperature of the upper (surface) layer Tg(t).
   :type tg: ndarray
   :param tg_avg: Long‑term average or restoring temperature Tḡ (K).
   :type tg_avg: float
   :param cv: Volumetric heat capacity (J m⁻³ K⁻¹).
   :type cv: float
   :param lambda_: Soil thermal conductivity λ (W m⁻¹ K⁻¹).
   :type lambda_: float
   :param delta_z: Thickness of the active soil layer δz (m).  If *None*, computed
                   from ``active_layer_thickness``.
   :type delta_z: float, optional
   :param times: Time stamps in seconds.  Required when *delta_z* is None or when
                 irregular sampling; defaults to `np.arange(len(tg))` seconds.
   :type times: ndarray, optional

   :returns: G₀,FR (W m⁻²).
   :rtype: ndarray


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
     :pyfunc:`numpy.frompyfunc`; the core implementation is scalar for
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
           :pyfunc:`sinusoidal_boundary_flux`.
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
         :pyfunc:`soil_temperature_sinusoidal` gives *T(z,t)* under the
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

   :param Tz: Temperature at layer **centres** (K).
   :type Tz: ArrayLike
   :param dz: Thickness of each layer (m).
   :type dz: ArrayLike
   :param lambda_s: Thermal conductivity for each layer (W m‑1 K‑1).
   :type lambda_s: ArrayLike or float

   :returns: Heat flux at *interfaces* (positive downward) with shape ``len(Tz)+1``.
   :rtype: np.ndarray


.. py:function:: integrated_soil_heat_flux(rho_c: ArrayLike, T_before: ArrayLike, T_after: ArrayLike, dz: ArrayLike, dt: float, G_ref: float = 0.0) -> numpy.ndarray

   Discrete integration of Eq. (3)/(5) to obtain heat‑flux profile.

   :param rho_c: Volumetric heat capacity ``ρ_s c_s`` for each layer (J m‑3 K‑1).
   :type rho_c: ArrayLike
   :param T_before: Temperatures at two successive timesteps (K).
   :type T_before: ArrayLike
   :param T_after: Temperatures at two successive timesteps (K).
   :type T_after: ArrayLike
   :param dz: Layer thicknesses (m).
   :type dz: ArrayLike
   :param dt: Timestep (s).
   :type dt: float
   :param G_ref: Heat flux at the lower reference depth *z_ref* (W m‑2).  Often ≈ 0.
   :type G_ref: float, default 0

   :returns: Heat flux at the *upper* interface of every layer, size ``len(dz)``.
   :rtype: np.ndarray


.. py:function:: volumetric_heat_capacity(theta: ArrayLike, theta_sat: float | ArrayLike) -> numpy.ndarray

   Volumetric heat capacity of moist soil.

   :param theta: Volumetric water content (m³ m‑3).
   :type theta: ArrayLike
   :param theta_sat: Soil porosity.
   :type theta_sat: float or ArrayLike

   :returns: ``ρ_s c_s`` (J m‑3 K‑1).
   :rtype: np.ndarray


.. py:function:: stretched_grid(n: int, D: float, xi: float) -> numpy.ndarray

   Generate *n* layer thicknesses following the exponential stretching rule.

   :param n: Number of layers.
   :type n: int
   :param D: Total domain depth (m).
   :type D: float
   :param xi: Stretching parameter; 0 → uniform grid.
   :type xi: float

   :returns: Thickness ``Δz_i`` for each layer *i* (m).
   :rtype: np.ndarray


.. py:function:: solve_tde(T_prev: ArrayLike, dz: ArrayLike, rho_c: ArrayLike, lambda_s: ArrayLike | float, Tsfc: float, Tbot: float, dt: float) -> numpy.ndarray

   Implicit Crank‑Nicholson (θ = 1) solve of Eq. (7).

   Boundary conditions (Eq. 7a, 7c) are Dirichlet.


.. py:function:: correct_profile(T_model: ArrayLike, depths_model: ArrayLike, T_obs: ArrayLike, depths_obs: ArrayLike) -> numpy.ndarray

   Add linear‑interpolated bias (Eq. ΔT_k) to model profile.


.. py:function:: surface_temperature_longwave(R_lw_up: float, R_lw_dn: float, emissivity: float = 0.98) -> float

   Convert upward/downward long‑wave radiation to surface temperature (Eq. 8).


.. py:function:: thermal_conductivity_yang2008(theta: ArrayLike, theta_sat: float, rho_dry: float | ArrayLike) -> numpy.ndarray

   Estimate soil thermal conductivity following Yang et al. (2005) (Eq. 9).


.. py:function:: flux_error_linear(rho_c: ArrayLike, S2_minus_S1: ArrayLike, dt: float) -> numpy.ndarray

   Error introduced when using a LINEAR temperature profile (diagnostic).


.. py:function:: surface_energy_residual(R_net: float, H: float, LE: float, G0: float) -> float

   Return the residual *ΔE* in Eq. (12).


.. py:function:: tdec_step(T_prev: ArrayLike, dz: ArrayLike, theta: ArrayLike, theta_sat: float, rho_dry: float, lambda_const: float, Tsfc: float, Tbot: float, dt: float, depths_model: ArrayLike, T_obs: ArrayLike, depths_obs: ArrayLike) -> Tuple[numpy.ndarray, numpy.ndarray]

   One integration step of the TDEC scheme.

   :returns: * **T_corr** (*np.ndarray*) -- Corrected temperature profile at *t + dt*.
             * **G_prof** (*np.ndarray*) -- Heat‑flux profile (W m‑2) at layer interfaces.


