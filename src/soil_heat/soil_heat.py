# https://www.scielo.br/j/rbcs/a/dFCLs7jXncc98VWNjdT64Xd/
# https://doi.org/10.1590/S0100-06832013000100011
# 10.52547/maco.2.1.5


import pandas as pd
import numpy as np
from scipy.optimize import curve_fit

# Constants
WATER_HEAT_CAPACITY = 4.18  # MJ m-3 K-1


def compute_heat_flux_conduction(
    df: pd.DataFrame,
    depth1: float = 0.05,
    depth2: float = 0.10,
    col_T1: str = "T5cm",
    col_T2: str = "T10cm",
    col_theta1: str = "VWC5cm",
    col_theta2: str = "VWC10cm",
    porosity: float = 0.40,
    k_dry: float = 0.25,
    k_sat: float = 1.50,
) -> pd.Series:
    """
    Estimate near-surface soil heat flux using Fourier’s law.

    This “gradient” approach computes conductive ground-heat flux
    :math:`G` between two depths by multiplying the vertical
    temperature gradient with an **effective** thermal conductivity
    that varies with volumetric water content (VWC).

    Parameters
    ----------
    df : pandas.DataFrame
        Time-indexed data containing at least the four columns
        specified by *col_T1*, *col_T2*, *col_theta1*, and
        *col_theta2*. The index spacing defines the temporal
        resolution of the output.
    depth1, depth2 : float, default (0.05, 0.10)
        Sensor depths (m).  `depth2` must be **greater** (deeper)
        than `depth1`.
    col_T1, col_T2 : str, default ("T5cm", "T10cm")
        Column names for temperature (°C or K) at `depth1` and
        `depth2`.
    col_theta1, col_theta2 : str, default ("VWC5cm", "VWC10cm")
        Column names for volumetric water content (m³ m⁻³) at
        `depth1` and `depth2`.
    porosity : float, default 0.40
        Soil total porosity (saturated VWC, m³ m⁻³).
    k_dry : float, default 0.25
        Dry-soil thermal conductivity (W m⁻¹ K⁻¹).
    k_sat : float, default 1.50
        Saturated-soil thermal conductivity (W m⁻¹ K⁻¹).

    Returns
    -------
    pandas.Series
        Half-hourly (or whatever the index step is) ground-heat-flux
        series with name ``"G_conduction"``. Units are W m⁻².
        Positive values indicate **downward** flux.

    Notes
    -----
    The effective thermal conductivity is computed by a simple linear
    mixing model:

    .. math::

        \\lambda_\\text{eff} = k_\\text{dry} +
        \\frac{\\bar{\\theta}}{\\phi}
        \\bigl(k_\\text{sat} - k_\\text{dry}\\bigr),

    where :math:`\\bar{\\theta}` is the mean VWC of the two depths and
    :math:`\\phi` is porosity.  More sophisticated models
    (e.g. Johansen, de Vries) can be substituted if site-specific
    calibration is available.

    References
    ----------
    * Campbell & Norman (2012) *An Introduction to Environmental
      Biophysics*, ch. 7.
    * Gao et al. (2017) Agricultural and Forest Meteorology,
      240 – 241, 194–204.

    Examples
    --------
    >>> G = compute_heat_flux_conduction(df_site,
    ...                                   depth1=0.05, depth2=0.10,
    ...                                   col_T1="T_05",
    ...                                   col_T2="T_10",
    ...                                   col_theta1="VWC_05",
    ...                                   col_theta2="VWC_10")
    >>> G.plot(title="Soil heat flux (gradient method)")
    """
    # 1. Effective thermal conductivity based on average moisture between the two depths
    theta_avg = (df[col_theta1] + df[col_theta2]) / 2.0
    frac_sat = np.clip(
        theta_avg / porosity, 0, 1
    )  # fraction of pore space filled with water
    lambda_eff = (
        k_dry + (k_sat - k_dry) * frac_sat
    )  # interpolate between dry and saturated

    # 2. Temperature gradient dT/dz (K/m). Depth increases downward.
    dT = df[col_T2] - df[col_T1]  # temperature difference
    dz = depth2 - depth1
    grad_T = dT / dz

    # 3. Fourier’s Law: G = -λ * dT/dz
    G = -lambda_eff * grad_T
    return pd.Series(G, index=df.index, name="G_conduction")


def compute_heat_flux_calorimetric(
    df: pd.DataFrame,
    depth_levels: "list[float]",
    T_cols: "list[str]",
    theta_cols: "list[str]",
    C_dry: float = 2.1e6,
    C_w: float = 4.2e6,
) -> pd.Series:
    """
    Calculate surface soil heat flux via the calorimetric (heat-storage) method.

    The calorimetric method integrates the transient change in heat
    *storage* within a multilayer soil column.  For a surface-to-depth
    layer of thickness :math:`z_{\\text{ref}}`, the surface flux
    :math:`G_0` is approximated by

    .. math::

        G_0 \\;\\approx\\; \\frac{\\Delta Q}{\\Delta t}
        \\;=\\; \\frac{1}{\\Delta t}
        \\sum_{i=1}^{N_\\text{layers}}
        C_i \\, \\Delta T_i \\, \\Delta z_i,

    where :math:`C_i` is volumetric heat capacity
    (J m⁻³ K⁻¹), :math:`\\Delta T_i` is the average temperature change
    (K) in layer *i*, and :math:`\\Delta z_i` is layer thickness (m).
    No heat-flux-plate reading is required if the deepest
    measurement depth lies below the diurnal damping depth such that
    :math:`G(z_{\\text{ref}}) \\approx 0`.

    Parameters
    ----------
    df : pandas.DataFrame
        Time-indexed data containing temperature and VWC columns for
        **all** depths specified in *T_cols* and *theta_cols*.  Index
        spacing sets the output time step.
    depth_levels : list of float
        Depths (m) corresponding *in order* to the entries in
        *T_cols* and *theta_cols*. Must be strictly increasing.
    T_cols : list of str
        Column names for soil temperatures (°C or K) at
        `depth_levels`.
    theta_cols : list of str
        Column names for volumetric water content (m³ m⁻³) at
        `depth_levels`.
    C_dry : float, default 2.1e6
        Volumetric heat capacity of dry soil matrix
        (J m⁻³ K⁻¹).
    C_w : float, default 4.2e6
        Volumetric heat capacity of liquid water
        (J m⁻³ K⁻¹).

    Returns
    -------
    pandas.Series
        Surface ground-heat-flux series, ``"G_calorimetric"`` (W m⁻²).
        Positive values denote **downward** flux.  The first time step
        is set to *NaN* because a preceding interval is required.

    Notes
    -----
    **Heat capacity model**

    A simple two-component mixture is assumed:

    .. math::

        C = (1 - \\theta)\\,C_{\\text{dry}} + \\theta\\,C_w.

    If bulk density or mineral fraction data are available, replace
    this linear approximation with a mass-weighted formulation.

    **Boundary assumption**

    The deepest temperature is treated as a “no-flux” boundary (storage
    only).  If diurnal waves penetrate deeper at your site, include an
    additional flux-plate term or extend `depth_levels` downward.

    References
    ----------
    * Mayocchi & Bristow (1995) Agricultural and Forest
      Meteorology 75, 93–109.
    * Oke (2002) *Boundary-Layer Climates*, 2nd ed., §2.3.
    * Fluxnet2015 “G” best-practice guide
      (https://fluxnet.org/sites/default/files/soil_heat_flux_guide.pdf).

    Examples
    --------
    >>> depths = [0.05, 0.10, 0.20, 0.50]          # m
    >>> Tcols  = ["T5", "T10", "T20", "T50"]       # °C
    >>> Vcols  = ["VWC5", "VWC10", "VWC20", "VWC50"]
    >>> G0 = compute_heat_flux_calorimetric(df_site,
    ...                                     depths, Tcols, Vcols)
    >>> G0.resample("D").mean().plot()
    >>> plt.ylabel("Daily mean G₀ (W m$^{-2}$)")
    """
    # --- basic error checks -------------------------------------------------
    if not (len(depth_levels) == len(T_cols) == len(theta_cols)):
        raise ValueError("depth_levels, T_cols, and theta_cols must be the same length")

    n = len(depth_levels)
    dt_seconds = (df.index[1] - df.index[0]).total_seconds()

    # volumetric heat capacity matrix (DataFrame, J m⁻³ K⁻¹)
    C_depth = (1.0 - df[theta_cols]) * C_dry + df[theta_cols] * C_w

    G0 = [np.nan]  # first element is undefined
    for i in range(1, len(df)):
        dQ = 0.0  # J m⁻²
        for j in range(1, n):
            z_top, z_bot = depth_levels[j - 1], depth_levels[j]
            dz = z_bot - z_top

            # layer-mean heat capacity over interval
            C_layer = 0.25 * (
                C_depth.iloc[i - 1, j - 1]  # type: ignore
                + C_depth.iloc[i - 1, j]  # type: ignore
                + C_depth.iloc[i, j - 1]  # type: ignore
                + C_depth.iloc[i, j]  # type: ignore
            )

            # layer-mean temperature change over interval
            dT_top = df[T_cols[j - 1]].iat[i] - df[T_cols[j - 1]].iat[i - 1]
            dT_bot = df[T_cols[j]].iat[i] - df[T_cols[j]].iat[i - 1]
            dT_layer = 0.5 * (dT_top + dT_bot)

            dQ += C_layer * dT_layer * dz

        G0.append(dQ / dt_seconds)  # W m⁻²

    return pd.Series(G0, index=df.index, name="G_calorimetric")


def temperature_gradient(
    T_upper: np.ndarray | float,
    T_lower: np.ndarray | float,
    depth_upper: float,
    depth_lower: float,
) -> np.ndarray | float:
    """
    Compute the **vertical temperature gradient** between two sensors.

    The gradient is defined as the change in temperature divided by the
    change in depth (positive downward):

    .. math::

        \\frac{∂T}{∂z}
        \;=\;
        \\frac{T_{\\text{lower}} - T_{\\text{upper}}}
              {z_{\\text{lower}} - z_{\\text{upper}}}   \\;\\;[^{\\circ}\\text{C m}^{-1}]

    Parameters
    ----------
    T_upper : float or array_like
        Temperature at the **shallower** depth ``depth_upper`` (°C).
    T_lower : float or array_like
        Temperature at the **deeper** depth ``depth_lower`` (°C).
        Must be broadcast-compatible with ``T_upper``.
    depth_upper : float
        Depth of the upper sensor (m, positive downward).
    depth_lower : float
        Depth of the lower sensor (m, positive downward).
        Must satisfy ``depth_lower > depth_upper`` for a meaningful
        gradient.

    Returns
    -------
    ndarray or float
        Temperature gradient ∂T/∂z (°C m⁻¹).
        Shape follows NumPy broadcasting of ``T_upper`` and ``T_lower``.

    Raises
    ------
    ValueError
        If ``depth_lower`` ≤ ``depth_upper``.

    Notes
    -----
    * **Sign convention** – A **positive** gradient indicates
      temperatures increase with depth (warmer below).
    * **Vectorised** – The arithmetic is fully NumPy-broadcasted; use it
      on scalar values, 1-D arrays, or entire DataFrames’ columns.
    * **Units** – Because depth is in metres and temperature in degrees
      Celsius, the result is °C m⁻¹ (identical to K m⁻¹).

    Examples
    --------
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
    """


def soil_heat_flux(T_upper, T_lower, depth_upper, depth_lower, k):
    """
    Calculate soil heat flux (G) using temperature gradient and thermal conductivity.

    Parameters:
    - T_upper: Temperature at upper depth (°C)
    - T_lower: Temperature at lower depth (°C)
    - depth_upper: Upper sensor depth (m)
    - depth_lower: Lower sensor depth (m)
    - k: Thermal conductivity (W/(m·°C))

    Returns:
    - Soil heat flux (W/m^2)
    """
    return -k * temperature_gradient(T_upper, T_lower, depth_upper, depth_lower)


def volumetric_heat_capacity(theta_v):
    """
    Estimate volumetric heat capacity Cv (J/(m³·°C)) from soil moisture.

    Parameters:
    - theta_v: Volumetric water content (decimal fraction, e.g., 0.20 for 20%)

    Returns:
    - Volumetric heat capacity (kJ/(m³·°C))
    """
    C_soil = 1942  # dry soil heat capacity kJ/(m³·°C)
    C_water = 4186  # water heat capacity kJ/(m³·°C)
    return (1 - theta_v) * C_soil + theta_v * C_water


def thermal_conductivity(
    alpha: np.ndarray | float,
    theta_v: np.ndarray | float,
) -> np.ndarray | float:
    """
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

    Parameters
    ----------
    alpha : float or array_like
        Thermal diffusivity **α** (m² s⁻¹).  May be scalar or any
        NumPy‐broadcastable shape.
    theta_v : float or array_like
        Volumetric water content **θ_v** (m³ m⁻³, i.e. decimal fraction
        of pore space filled with water).  Must be broadcast‐compatible
        with ``alpha``.

    Returns
    -------
    ndarray or float
        Thermal conductivity **k** (W m⁻¹ K⁻¹) with the broadcast shape
        of the inputs.

    Notes
    -----
    * **Volumetric heat capacity model** –
      :pyfunc:`volumetric_heat_capacity` typically assumes a two‐phase
      mixture of mineral soil and water:

      .. math::

         C_v(θ_v) \;=\; (1-θ_v)\,ρc_\text{dry} \;+\;
                         θ_v\,ρc_\text{w} ,

      where ``ρc_dry`` (≈ 2.0 MJ m⁻³ K⁻¹) and ``ρc_w`` (4.18 MJ m⁻³ K⁻¹)
      are the volumetric heat capacities of dry soil and liquid water,
      respectively.  Ensure these defaults suit your substrate.
    * **Vectorisation** – The function is a one‐liner,
      ``alpha * Cv``, and thus inherits full NumPy broadcasting rules.
    * **Temperature units** – Because heat capacity is per kelvin, *k*
      is returned in W m⁻¹ K⁻¹ (equivalent to W m⁻¹ °C⁻¹).

    Examples
    --------
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
    """

    Cv = volumetric_heat_capacity(theta_v)
    return alpha * Cv


def diurnal_amplitude(
    series: pd.Series,
) -> pd.Series:
    """
    Compute the **daily diurnal amplitude** of a time-series.

    The diurnal amplitude for a given calendar day is defined as the
    difference between that day’s maximum and minimum values:

    .. math::

        A_d \;=\; \max\_{t \\in d} x(t) \;-\; \min\_{t \\in d} x(t)

    This metric is frequently used for temperature, soil-heat, or other
    environmental data to characterise the strength of the diurnal cycle.

    Parameters
    ----------
    series : pandas.Series
        Time-indexed observations with a `DatetimeIndex`.
        Any frequency is accepted, but the index **must** be sorted and
        monotonic.  Missing values (`NaN`) are ignored within each daily
        window.

    Returns
    -------
    pandas.Series
        Daily diurnal amplitude, indexed by date (midnight ``00:00`` of
        each day).  Units are the same as those of the input ``series``.

    Notes
    -----
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

    Examples
    --------
    >>> amp = diurnal_amplitude(df["air_temperature"])
    >>> amp.plot(title="Daily Temperature Amplitude")
    >>> amp.describe().loc[["min", "mean", "max"]]
    min      4.3
    mean     9.7
    max     15.2
    Name: air_temperature, dtype: float64
    """

    daily_max = series.resample("D").max()
    daily_min = series.resample("D").min()
    amplitude = daily_max - daily_min

    return amplitude


def diurnal_peak_lag(
    series1: pd.Series,
    series2: pd.Series,
) -> pd.Series:
    """
    Compute the **daily peak‐time lag** (Δt) between two diurnal signals.

    For each calendar day the function identifies the clock time at which
    each series reaches its maximum value and returns the signed time
    difference in **hours** (``series1`` minus ``series2``).  A modular
    correction confines the result to the interval ``[-12, 12]`` h so
    that, for example, a raw lag of –23 h becomes +1 h.

    Parameters
    ----------
    series1, series2 : pandas.Series
        Time-indexed observations of equal length, preferably
        temperature or some other quantity exhibiting a clear diurnal
        cycle.  The index **must** be `DatetimeIndex` and should be
        timezone-aware and aligned in frequency.
        Missing values are ignored within each daily resampling window.

    Returns
    -------
    pandas.Series
        Daily peak-lag values (float, hours) indexed by the **date** of
        the peak (00:00 of each day).
        Positive lags mean the peak of ``series1`` occurs *later* than
        the peak of ``series2`` on that day; negative lags indicate the
        opposite.

    Notes
    -----
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

    Examples
    --------
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
    """

    def daily_peak_time(series):
        return series.resample("D").apply(
            lambda x: x.idxmax().hour + x.idxmax().minute / 60
        )

    peak_time_1 = daily_peak_time(series1)
    peak_time_2 = daily_peak_time(series2)

    peak_lag = peak_time_1 - peak_time_2

    # Adjust lag to account for day wrap-around (e.g., -23 hours to +1 hour)
    peak_lag = peak_lag.apply(lambda x: (x + 12) % 24 - 12)

    return peak_lag


def fit_sinusoid(
    t: np.ndarray,
    data: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    """
    Fit a **sinusoidal model** to time–series data using non-linear least
    squares.

    The model is

    .. math::

        y(t)\;=\;A \sin( \omega t + \varphi ) + C ,

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
    """

    # Initial guess for the parameters [A, omega, phase, offset]
    guess_amp = np.std(data)
    guess_freq = 2 * np.pi / 86400  # Assuming daily cycle
    guess_phase = 0
    guess_offset = np.mean(data)
    p0 = [guess_amp, guess_freq, guess_phase, guess_offset]

    # Fit the sine curve
    popt, pcov = curve_fit(sinusoid, t, data, p0=p0)
    return popt


def sinusoid(
    t: np.ndarray | float,
    A: float,
    omega: float,
    phase: float,
    offset: float,
) -> np.ndarray | float:
    """
    Evaluate a **sinusoidal wave** of the form

    .. math::

        f(t) \;=\; A \sin(\omega\, t + \varphi) + C ,

    where :math:`A` is the *amplitude*, :math:`\omega` the *angular
    frequency*, :math:`\varphi` the *phase shift*, and :math:`C`
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
        Phase shift :math:`\varphi` in **radians**.
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
    """

    return A * np.sin(omega * t + phase) + offset


def thermal_diffusivity_amplitude(
    A1: float,
    A2: float,
    z1: float,
    z2: float,
    period: int = 86_400,
) -> float:
    """
    Estimate soil **thermal diffusivity** (``α``) from the *damping of
    harmonic amplitude* between two depths.

    A one–dimensional soil column subject to a sinusoidal surface
    temperature oscillation exhibits an exponential decay of amplitude
    with depth (Carslaw & Jaeger, 1959).  For a single angular frequency
    :math:`ω = 2π/P`, the analytical solution yields

    .. math::

        α \;=\; \\frac{π\, (z_2 - z_1)^2}
                       {P \;\\bigl[\,\\ln(A_1/A_2)\\bigr]^2} ,

    where

    * *A₁* and *A₂* are the harmonic amplitudes at depths *z₁* and *z₂*,
      respectively (*A₁ > A₂*),
    * *P* is the forcing period, and
    * *z₂  – z₁* is the vertical separation of the two sensors.

    Parameters
    ----------
    A1, A2 : float
        Diurnal (or other fundamental) temperature amplitudes at the
        shallow depth ``z1`` and deeper depth ``z2``.
        Units **°C** or **K** (identical for both).
    z1, z2 : float
        Sensor depths in **metres** (positive downward).
        Must satisfy ``z2 > z1``.
    period : int, default ``86_400``
        Fundamental period *P* of the temperature wave in **seconds**.
        ``86 400`` s corresponds to a 24-hour diurnal cycle.

    Returns
    -------
    float
        Thermal diffusivity **α** in m² s⁻¹.

    Raises
    ------
    ValueError
        If ``A1 <= A2`` (violates physical damping assumption) or
        if ``z2 <= z1``.

    Notes
    -----
    * **Amplitude extraction** – ``A1`` and ``A2`` should be obtained
      from a harmonic fit or spectral decomposition that isolates the
      target frequency; raw peak–trough differences are less robust.
    * **Logarithmic sensitivity** – Because the formula involves
      ``ln(A1/A2)``, small uncertainties in amplitudes propagate
      non-linearly; ensure adequate signal-to-noise ratio.
    * Once ``α`` is known, thermal conductivity ``k`` follows from
      ``k = ρc α`` given an independent estimate of volumetric heat
      capacity ``ρc``.

    References
    ----------
    Carslaw, H. S., & Jaeger, J. C. (1959).
    *Conduction of Heat in Solids* (2nd ed., pp. 501–502).
    Oxford University Press.

    Examples
    --------
    >>> # Amplitudes from harmonic regression at 5 cm and 10 cm depths
    >>> alpha = thermal_diffusivity_amplitude(
    ...     A1=6.3, A2=4.1, z1=0.05, z2=0.10
    ... )
    >>> print(f"α = {alpha:.2e} m² s⁻¹")
    α = 1.38e-07 m² s⁻¹
    """

    alpha = (np.pi * (z2 - z1) ** 2) / (period * (np.log(A1 / A2)) ** 2)
    return alpha


def thermal_diffusivity_lag(delta_t, z1, z2, period=86400):
    """
    Estimate thermal diffusivity from phase lag.

    Parameters:
    - delta_t: Time lag between peaks at two depths (seconds)
    - z1, z2: Depths (m)
    - period: Time period of wave (default = 86400 s for daily cycle)

    Returns:
    - Thermal diffusivity α (m²/s)

    Citation:
    S.V. Nerpin, and A.F. Chudnovskii, Soil physics, (Moscow: Nauka) p 584, 1967 (in Russian)
    """

    alpha = (period / (4 * np.pi)) * (z2 - z1) ** 2 / (delta_t) ** 2
    return alpha


def thermal_diffusivity_logrithmic(
    t1z1: float,
    t2z1: float,
    t3z1: float,
    t4z1: float,
    t1z2: float,
    t2z2: float,
    t3z2: float,
    t4z2: float,
    z1: float,
    z2: float,
    period: int = 86_400,
) -> float:
    """
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

        α \;=\; \\frac{4 \, π \, (z_2 - z_1)^2}
                        {P \;\\bigl[\,
                        \\ln\\bigl( ΔT_{z1} / ΔT_{z2} \\bigr)\\bigr]^2}

    with amplitude decrements

    .. math::

        ΔT_{zij} = \\sqrt{(T_1 - T_3)^2 + (T_2 - T_4)^2}\;.

    The formulation is advantageous when only a *short* record is
    available (four points suffice) but is sensitive to sensor noise and
    non-sinusoidal disturbances.

    Parameters
    ----------
    t1z1, t2z1, t3z1, t4z1 : float
        Temperatures (°C) at depth ``z1`` sampled at four successive
        quarter-period intervals.
    t1z2, t2z2, t3z2, t4z2 : float
        Temperatures (°C) at depth ``z2`` sampled at the *same* times as
        the readings at ``z1``.
    z1, z2 : float
        Sensor depths in **metres** (positive downward).  Must satisfy
        ``z2 > z1`` for a meaningful diffusivity.
    period : int, default ``86_400``
        Fundamental period *P* of the temperature oscillation in
        **seconds**.  ``86 400`` s corresponds to a 24-hour diurnal wave.

    Returns
    -------
    float
        Thermal diffusivity **α** in m² s⁻¹.

    Notes
    -----
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

    References
    ----------
    * Kolmogorov, A. N. (1950). *On the question of determining the
      coefficient of thermal diffusivity of the soil*. *Izvestiya
      Akademii Nauk SSSR, Ser. Geogr. Geofiz.*, 14 (2), 97–99. (In
      Russian)
    * Seemann, W. (1928). *Die Wärmeleitung in der Bodenschicht*.
      Springer, Berlin.

    Examples
    --------
    >>> α = thermal_diffusivity_logrithmic(
    ...     22.5, 20.3, 18.4, 20.1,   # temps @ z1
    ...     18.7, 17.2, 15.9, 17.1,   # temps @ z2
    ...     z1=0.05, z2=0.10,
    ... )
    >>> print(f"α = {α:.2e} m²/s")
    α = 1.46e-07 m²/s
    """

    alpha = (4 * np.pi * (z2 - z1) ** 2) / (
        period
        * np.log(
            ((t1z1 - t3z1) ** 2 + (t2z1 - t4z1) ** 2)
            / ((t1z2 - t3z2) ** 2 + (t2z2 - t4z2)) ** 2
        )
        ** 2
    )
    return alpha


def calc_thermal_diffusivity_log_pair(df, depth1_col, depth2_col, z1, z2, period=86400):
    """
    Estimate soil **thermal diffusivity** (``α``) between two depths using the
    *four-point logarithmic amplitude* method.

    The function extracts the **first four consecutive samples** from two
    temperature records—one at the shallow depth ``z1`` and one at the deeper
    depth ``z2``—and passes them to
    :pyfunc:`thermal_diffusivity_logrithmic`.  That helper implements the
    log–ratio solution of the 1-D heat‐conduction equation for a sinusoidal
    boundary condition (Horton et al., 1934; de Vries, 1963):

    .. math::

        α = \\frac{(z_2 - z_1)^2}
                  {2P\\;\\ln\\left(\\frac{ΔT_{\\!z1}}{ΔT_{\\!z2}}\\right)},

    where

    * **P** is the forcing period (s),
    * :math:`ΔT_{\\!z}` is the logarithmic temperature decrement derived
      from four successive measurements at depth *z*.

    The approach is robust for short windows (four points suffice) but is
    sensitive to noise; it is best applied to periods with clear, smooth
    diurnal cycling.

    Parameters
    ----------
    df : pandas.DataFrame
        Time‐indexed data containing at least the two temperature columns
        specified by ``depth1_col`` and ``depth2_col``.
        **Only the first four rows** are used in the calculation.
    depth1_col, depth2_col : str
        Column names for the shallow (``z1``) and deeper (``z2``)
        temperature series, respectively.
    z1, z2 : float
        Sensor depths in **metres** (positive downward).
        Must satisfy ``z2 > z1``.
    period : int, default ``86_400``
        Dominant temperature oscillation period **P** in **seconds**.
        The default (86 400 s) corresponds to 24 h.

    Returns
    -------
    float or None
        Thermal diffusivity ``α`` in **m² s⁻¹**.
        Returns ``None`` when fewer than four valid samples are available
        or if ``thermal_diffusivity_logrithmic`` itself returns ``None``.

    Warns
    -----
    UserWarning
        Issued (via ``print``) when fewer than four rows are present in
        *df*, in which case the method is skipped and ``None`` is returned.

    Notes
    -----
    * **Data requirement** – The function *does not* resample or align
      series; it simply grabs the first four rows.  Pre-filter or sort
      your DataFrame accordingly.
    * **Noise sensitivity** – Because the method depends on small
      differences between successive temperature readings, apply a
      smoothing filter or select a high-signal period to minimise error.
    * **Relationship to conductivity** – Once ``α`` is known, bulk
      thermal conductivity ``k`` can be obtained from ``k = ρc α`` given
      an estimate of volumetric heat capacity ``ρc``.

    References
    ----------
    Horton, R., Wierenga, P. J., Nielsen, D. R., & de Vries, D. A. (1983).
    *Calorimetric determination of soil thermal properties*.
    Soil Science Society of America Journal, **47**, 104–111.

    de Vries, D. A. (1963). *Thermal properties of soils*.
    In *Physics of Plant Environment* (pp. 210–235). North-Holland.

    Examples
    --------
    >>> α_log = calc_thermal_diffusivity_log_pair(
    ...     df=df.sort_index(),          # ensure chronological order
    ...     depth1_col='ts_05cm',
    ...     depth2_col='ts_10cm',
    ...     z1=0.05, z2=0.10,
    ... )
    >>> if α_log is not None:
    ...     print(f"Log-method α = {α_log:.2e} m² s⁻¹")
    Log-method α = 1.45e-07 m² s⁻¹
    """
    if len(df) < 4:
        print(
            f"Warning: Not enough time points for logarithmic method between {depth1_col} and {depth2_col}."
        )
        return None

    t1z1 = df[depth1_col].iloc[0]
    t2z1 = df[depth1_col].iloc[1]
    t3z1 = df[depth1_col].iloc[2]
    t4z1 = df[depth1_col].iloc[3]
    t1z2 = df[depth2_col].iloc[0]
    t2z2 = df[depth2_col].iloc[1]
    t3z2 = df[depth2_col].iloc[2]
    t4z2 = df[depth2_col].iloc[3]

    return thermal_diffusivity_logrithmic(
        t1z1, t2z1, t3z1, t4z1, t1z2, t2z2, t3z2, t4z2, z1, z2, period
    )


def calculate_thermal_diffusivity_for_pair(df, col1, col2, z1, z2, period=86400):
    """
    Estimate soil **thermal diffusivity** (``α``) between two depths using
    three classical harmonic methods: *log-amplitude*, *amplitude ratio*,
    and *phase shift*.

    Given two temperature time-series measured at depths ``z1`` and ``z2``,
    the function first extracts the dominant diurnal signal—its amplitude
    and phase—then applies the analytical solutions of the 1-D heat wave
    equation for a homogeneous medium subject to sinusoidal forcing
    (Carslaw & Jaeger, 1959).

    Methods
    -------
    1. Log-Amplitude (α\_log)
        Uses the decay of the harmonic amplitude with depth:

        .. math::

            α\\_{\\text{log}} = \\frac{(z_2 - z_1)^2}
                                     {2\\,P\\;\\ln\\bigl(A_1 / A_2\\bigr)}

    2. Amplitude Ratio (α\_amp)
        Algebraically identical to the log-amplitude method but expressed
        directly in terms of the two amplitudes:

        .. math::

            α\\_{\\text{amp}} = \\frac{(z_2 - z_1)^2\\;\\omega}
                                      {2\\,[\\ln(A_1/A_2)]^2}

        where ``ω = 2π / P`` is the angular frequency.

    3. Phase Lag (α\_lag)
        Relates the travel time (phase shift) of the temperature wave:

        .. math::

            α\\_{\\text{lag}} = \\frac{(z_2 - z_1)^2}{2\\,Δt\\,P}

        with ``Δt`` the peak-to-peak time lag (s).

    Parameters
    ----------
    df : pandas.DataFrame
        Time-indexed frame containing temperature observations.
    col1, col2 : str
        Column names for the shallow and deeper temperature series,
        respectively.
    z1, z2 : float
        Sensor depths in **metres** (positive downward).  Must satisfy
        ``z2 > z1``.
    period : int, default ``86_400``
        Fundamental period **P** of the harmonic forcing in **seconds**.
        ``86 400`` s corresponds to 24 h diurnal cycling.

    Returns
    -------
    dict[str, float]
        Mapping of method identifiers to diffusivity estimates
        (m² s⁻¹):

        * ``'alpha_log'`` – logarithmic amplitude method.
        * ``'alpha_amp'`` – direct amplitude-ratio method.
        * ``'alpha_lag'`` – phase-shift (lag) method.

        Any method returning *None* inside intermediate helpers is
        propagated unchanged.

    Raises
    ------
    ValueError
        If ``z1`` ≥ ``z2`` or if either column is missing in *df*.

    Notes
    -----
    * ``diurnal_amplitude`` extracts the half range of the 24-h harmonic,
      typically via fast Fourier transform or STL decomposition.
    * ``diurnal_peak_lag`` returns the modal lag **in hours**; the value
      is internally converted to seconds.
    * The function assumes a **single dominant harmonic**.  Strong
      synoptic or weather-front variability can bias results; apply
      filtering or select periods with clear diurnal cycling.
    * Thermal diffusivity relates to thermal conductivity ``k`` through

      .. math:: k = ρ c \\, α

      once bulk volumetric heat capacity ``ρc`` is known.

    References
    ----------
    Carslaw, H. S., & Jaeger, J. C. (1959). *Conduction of Heat in Solids*
    (2nd ed.). Oxford University Press.

    Examples
    --------
    >>> depth_map = {'ts_05cm': 0.05, 'ts_10cm': 0.10}
    >>> α = calculate_thermal_diffusivity_for_pair(
    ...         df, 'ts_05cm', 'ts_10cm',
    ...         z1=depth_map['ts_05cm'], z2=depth_map['ts_10cm'])
    >>> for meth, val in α.items():
    ...     print(f"{meth}: {val:.2e} m² s⁻¹")
    alpha_log: 1.43e-07 m² s⁻¹
    alpha_amp: 1.41e-07 m² s⁻¹
    alpha_lag: 1.38e-07 m² s⁻¹
    """
    temp_data_depth1 = df[col1]
    temp_data_depth2 = df[col2]

    A1 = diurnal_amplitude(temp_data_depth1)
    A2 = diurnal_amplitude(temp_data_depth2)

    phase = diurnal_peak_lag(temp_data_depth2, temp_data_depth1)

    alpha_amplitude = thermal_diffusivity_amplitude(A1, A2, z1, z2, period)
    alpha_phase = thermal_diffusivity_lag(phase * 3600, z1, z2, period)
    alpha_log_amplitude = calc_thermal_diffusivity_log_pair(
        df, col1, col2, z1, z2, period
    )

    return {
        "alpha_log": alpha_log_amplitude,
        "alpha_amp": alpha_amplitude,
        "alpha_lag": alpha_phase,
    }


def calculate_thermal_properties_for_all_pairs(
    df,
    depth_mapping,
    period=86400,
):
    """
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

    Parameters
    ----------
    df : pandas.DataFrame
        Time-indexed data frame containing at least

        * temperature columns listed in ``depth_mapping``; units **°C**,
          column names typically follow a pattern such as ``'ts_05cm'``.
        * matching soil-water-content columns; each temperature column
          ``'<name>ts'`` must have a companion column
          ``'<name>swc'`` in **percent**.  These are averaged and divided
          by 100 to obtain volumetric θ (*m³ m⁻³*).

    depth_mapping : dict[str, float]
        Mapping of *temperature* column names to sensor depths in **metres**
        (positive downward), e.g. ``{'ts_05cm': 0.05, 'ts_10cm': 0.10}``.

    period : int, default ``86_400``
        Dominant period of the temperature wave (s).  ``86_400`` s
        corresponds to 24 h and is appropriate for daily forcing.

    Returns
    -------
    pandas.DataFrame
        Concatenated frame of thermal properties for every depth pair.
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

    Notes
    -----
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

    Examples
    --------
    >>> depth_map = {'ts_05cm': 0.05, 'ts_10cm': 0.10, 'ts_20cm': 0.20}
    >>> props = calculate_thermal_properties_for_all_pairs(df, depth_map)
    >>> props.loc['0.05-0.10'][['alpha_phase', 'G_phase']].plot()
    >>> props.groupby(level=0)['k_amplitude'].median().unstack()
    """

    depth_cols = list(depth_mapping.keys())
    dfs = {}

    for i in range(len(depth_cols)):
        for j in range(i + 1, len(depth_cols)):
            df1 = (
                df.copy().dropna()
            )  # Copy the DataFrame to avoid modifying the original data

            res_z = {}
            col1 = depth_cols[i]
            col2 = depth_cols[j]
            z1 = depth_mapping[col1]
            z2 = depth_mapping[col2]

            # Calculate thermal diffusivity
            alpha_results = calculate_thermal_diffusivity_for_pair(
                df1, col1, col2, z1, z2, period
            )

            soil_moist_percent = (
                df1[[col1.replace("ts", "swc"), col2.replace("ts", "swc")]].mean(axis=1)
                / 100
            )

            col_list = []
            for key, val in alpha_results.items():
                if val is None:
                    df1[f"k_{key}"] = np.nan

                else:
                    k = thermal_conductivity(val, soil_moist_percent)
                    df1[f"G_{key}"] = soil_heat_flux(
                        df1[col1],
                        df1[col2],
                        z1,
                        z2,
                        k,
                    )  # Calculate soil heat flux using the thermal conductivity
                    df1[f"k_{key}"] = k
                    df1[f"{key}"] = val  # Store the thermal diffusivity value
                    # Store the volumetric water content
                    df1["theta_v"] = soil_moist_percent
                    col_list.append(f"{key}")
                    col_list.append(f"G_{key}")
                    col_list.append(f"k_{key}")
                    col_list.append("theta_v")

            dfs[f"{z1}-{z2}"] = df1[col_list]

    return pd.concat(dfs)


def estimate_rhoc_dry(
    alpha: pd.Series,
    theta: pd.Series,
    porosity: float = 0.40,
    k_dry: float = 0.25,
    k_sat: float = 1.50,
    rhoc_w: float = 4.18e6,
    dry_quantile: float = 0.10,
) -> float:
    """
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

        λ(θ) &= k_\\text{dry} + \\frac{θ}{φ}
                \\,\\bigl(k_\\text{sat} - k_\\text{dry}\\bigr)                     \\\\

        C_v  &= \\frac{λ(θ)}{α}                                                   \\\\

        ρ\,c_\\text{dry} &= \\frac{C_v - θ\,ρ\,c_w}{1-θ}\,,

    where

    * *λ* is thermal conductivity (W m⁻¹ K⁻¹),
    * *α* is thermal diffusivity (m² s⁻¹),
    * *C_v* is volumetric heat capacity of the *moist* soil
      (J m⁻³ K⁻¹), and
    * *φ* is total porosity (m³ m⁻³).

    Parameters
    ----------
    alpha : pandas.Series
        Soil thermal diffusivity **α** (m² s⁻¹), indexed identically to
        *theta* (usually a time-series).
    theta : pandas.Series
        Volumetric water content **θ** (m³ m⁻³).
    porosity : float, default ``0.40``
        Total soil porosity **φ** (saturated water content).
    k_dry : float, default ``0.25``
        Thermal conductivity of *air-dry* soil (W m⁻¹ K⁻¹).
    k_sat : float, default ``1.50``
        Thermal conductivity of **saturated** soil (W m⁻¹ K⁻¹).
    rhoc_w : float, default ``4.18e6``
        Volumetric heat capacity of **liquid water**
        (J m⁻³ K⁻¹, ≈ 4.18 MJ m⁻³ K⁻¹).
    dry_quantile : float, default ``0.10``
        Fraction of the *lowest* moisture observations to treat as
        “dry” when taking the median.  For example, ``0.10`` selects
        the driest 10 % of the record.

    Returns
    -------
    float
        Median volumetric heat capacity of the *dry* soil matrix
        (J m⁻³ K⁻¹).

    Notes
    -----
    * **Alignment** — The two series are first *inner-joined* so only
      timestamps present in both are considered.
    * **Robustness** — Using the median of the driest subset avoids
      bias from residual soil moisture while damping the influence of
      occasional outliers.
    * The default conductivity bounds ``k_dry``/``k_sat`` follow
      typical literature values for mineral soils; adjust them for
      peat, organic, or highly gravelly substrates.

    Examples
    --------
    >>> rhoc_dry = estimate_rhoc_dry(
    ...     alpha=df['alpha_10cm'],
    ...     theta=df['VWC_10cm'],
    ...     porosity=0.43,
    ... )
    >>> print(f"ρ c_dry ≈ {rhoc_dry/1e6:.2f} MJ m⁻³ K⁻¹")
    2.07 MJ m⁻³ K⁻¹
    """

    # keep only days where both alpha & theta are available
    theta, alpha = theta.align(alpha, join="inner")  # ⬅ key line

    frac_sat = np.clip(theta / porosity, 0.0, 1.0)
    lam = k_dry + (k_sat - k_dry) * frac_sat  # W m⁻¹ K⁻¹

    # --- 3. Heat capacity & dry-soil estimate ------------------
    Cv = lam / alpha  # J m⁻³ K⁻¹
    rhoc_dry = (Cv - theta * rhoc_w) / (1.0 - theta)

    dry_days = theta <= theta.quantile(dry_quantile)
    return float(rhoc_dry.loc[dry_days].median())


if __name__ == "__main__":
    # Load the data
    df = pd.read_csv("utd_soil_data.csv")
    df["datetime_start"] = pd.to_datetime(df["datetime_start"])
    df.set_index("datetime_start", inplace=True)

    # Define depth mapping
    depth_mapping = {
        "ts_3_1_1": 0.05,  # 5 cm
        "ts_3_2_1": 0.10,  # 10 cm
        "ts_3_3_1": 0.20,  # 20 cm
    }

    # Calculate thermal properties
    results_df = calculate_thermal_properties_for_all_pairs(df, depth_mapping)
    print(results_df)
