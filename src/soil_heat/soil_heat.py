# https://www.scielo.br/j/rbcs/a/dFCLs7jXncc98VWNjdT64Xd/
# https://doi.org/10.1590/S0100-06832013000100011
# 10.52547/maco.2.1.5


import pandas as pd
import numpy as np
from scipy.optimize import curve_fit

# Constants
WATER_HEAT_CAPACITY = 4.18  # MJ m-3 K-1

import numpy as np


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


def temperature_gradient(T_upper, T_lower, depth_upper, depth_lower):
    """
    Calculate the temperature gradient (°C/m) between two depths.

    Parameters:
    - T_upper: Temperature at upper depth (°C)
    - T_lower: Temperature at lower depth (°C)
    - depth_upper: Upper sensor depth (m)
    - depth_lower: Lower sensor depth (m)

    Returns:
    - Temperature gradient (°C/m)
    """
    return (T_lower - T_upper) / (depth_lower - depth_upper)


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


def thermal_conductivity(alpha, theta_v):
    """
    Calculate thermal conductivity (k) from diffusivity and moisture.

    Parameters:
    - alpha: Thermal diffusivity (m²/s)
    - theta_v: Volumetric water content (decimal fraction)

    Returns:
    - Thermal conductivity k (W/(m·°C))
    """
    Cv = volumetric_heat_capacity(theta_v)
    return alpha * Cv


def diurnal_amplitude(series: pd.Series) -> pd.Series:
    """
    Calculate the amplitude of diurnal fluctuations (daily max - daily min) from a pandas series.

    Parameters:
    series (pd.Series): Pandas Series with a datetime index.

    Returns:
    pd.Series: Daily amplitude values.
    """
    daily_max = series.resample("D").max()
    daily_min = series.resample("D").min()
    amplitude = daily_max - daily_min

    return amplitude


def diurnal_peak_lag(
    series1: pd.Series, series2: pd.Series
) -> pd.Series | pd.DataFrame:
    """
    Calculate the daily lag (offset in time) between peaks of two diurnal pandas series.

    Parameters:
    series1 (pd.Series): First pandas Series with datetime index.
    series2 (pd.Series): Second pandas Series with datetime index.

    Returns:
    pd.Series: Daily peak lag values in hours.
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


def fit_sinusoid(t, data):
    """
    Fits a sinusoidal curve to the temperature data.
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


def sinusoid(t, A, omega, phase, offset):
    """
    Sinusoidal function.
    """
    return A * np.sin(omega * t + phase) + offset


def thermal_diffusivity_amplitude(A1, A2, z1, z2, period=86400):
    """
    Estimate thermal diffusivity from amplitude damping.

    Parameters:
    - A1, A2: Amplitudes at shallow (z1) and deep (z2) depths
    - z1, z2: Depths (m)
    - period: Time period of wave (default = 86400 s for daily cycle)

    Returns:
    - Thermal diffusivity α (m²/s)

    Citation:
    H.J. Carslaw, and J.C. Jaeger, Conduction of heat in solids (2nd edition), Oxford University Press, NewYork, 510 p, 1959
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
    t1z1, t2z1, t3z1, t4z1, t1z2, t2z2, t3z2, t4z2, z1, z2, period=86400
):
    """
    Estimate thermal diffusivity from Seemann's method.

    Parameters:
    - t11,t12,t13,t14: Temperatures at depth z1 at 4 different times
    - t21,t22,t23,t24: Temperatures at depth z2 at 4 different times
    - z1, z2: Depths (m)

    Returns:
    - Thermal diffusivity α (m²/s)

    Citation:
      A. N. Kolmogorov, On the question of determining the coefficient of thermal diffusivity of the soil, Izv.Acad. Sci. USSR. Geogr.Geophys., 2(14), 97–99, 1950 (in Russian).
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
    """Calculates thermal diffusivity for a pair of depths using the logarithmic method."""
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
    Calculates thermal diffusivity using the log-amplitude, amplitude, and phase methods.

    Parameters:
        df (pd.DataFrame): DataFrame with temperature columns.
        col1 (str): Name of the first depth column.
        col2 (str): Name of the second depth column.
        z1 (float): Depth of the first sensor (m).
        z2 (float): Depth of the second sensor (m).
        period (int): Period of the temperature oscillation (seconds).

    Returns:
        dict: A dictionary containing thermal diffusivity calculated by each method.
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
    Calculates thermal diffusivity, volumetric heat capacity, and soil heat flux for all depth pairs.

    Parameters:
        df (pd.DataFrame): Time-indexed DataFrame with temperature columns for different depths.
        depth_mapping (dict): A dictionary where keys are column names and values are their corresponding depths (in meters).
        period (int): The period of the temperature oscillation (in seconds, default is 24 hours).
        dry_density (float): Dry density of soil (Mg/m^3).
        moisture_content (float): Volumetric moisture content of soil (m^3/m^3).

    Returns:
        pd.DataFrame: A DataFrame containing the calculated thermal properties for each depth pair.
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
