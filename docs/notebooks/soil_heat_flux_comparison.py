import marimo

__generated_with = "0.18.4"
app = marimo.App()


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## Import relevant libraries
    """)
    return


@app.cell
def _():
    import marimo as mo
    import plotly.express as px
    import plotly.graph_objects as go
    import pandas as pd
    import numpy as np
    import matplotlib.pyplot as plt
    from pathlib import Path
    import sys
    return go, mo, np, pd, plt, px, sys


@app.cell
def _(sys):
    sys.path.append('../../src')
    import soil_heat     
    return (soil_heat,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## Load SoilVUE Data
    """)
    return


@app.cell
def _(pd):
    # ---------- SoilVUE data -----------------------------------
    soilvue_path_1 = "21031_Statistics_AmeriFlux37.dat"
    soilvue_path_2 = "21031_Statistics_AmeriFlux38.dat"

    # Combine the two SoilVUE files into one DataFrame
    df_sv_1 = pd.read_csv(
        soilvue_path_1,
        na_values=['NAN'],
        engine='python',
    )
    df_sv_2 = pd.read_csv(
        soilvue_path_2,
        na_values=['NAN'],
        engine='python',
    )

    # Use the Python CSV engine for robust parsing (mixed quoting, many columns)
    df_sv = pd.concat([df_sv_1, df_sv_2], ignore_index=True)

    # Back‑fill the one timestamp so that every row carries a datetime label
    df_sv['TIMESTAMP_START'] = df_sv['TIMESTAMP_START'].ffill()

    # Parse to pandas datetime (format YYYYMMDDhhmm)
    df_sv['datetime'] = pd.to_datetime(df_sv['TIMESTAMP_START'], format='%Y%m%d%H%M')

    # For this demo we keep only the channels needed for the gradient method
    df_sv = df_sv.astype({'T_1_1_1': float,
                          'T_1_2_1': float,
                            'T_1_3_1': float,
                            'T_1_4_1': float,
                            'T_1_5_1': float,
                            'T_1_6_1': float,
                            'T_1_7_1': float,
                          'VWC_1_1_1': float,
                          'VWC_1_2_1': float,
                            'VWC_1_3_1': float,
                            'VWC_1_4_1': float,
                            'VWC_1_5_1': float,
                            'VWC_1_6_1': float,
                            'VWC_1_7_1': float,
                            'T_CANOPY': float, 
                            'NETRAD': float,
                            'SW_IN': float,
                            'SW_OUT': float,
                            'LW_IN': float,
                            'LW_OUT': float,
                          })
    return (df_sv,)


@app.cell
def _(df_sv):
    #df_sv['datetime'] = pd.date_range(start='2024-07-11 05:30', periods=589, freq='30min')
    # Package expects the DataFrame columns to be named explicitly
    df_grad =df_sv.rename(columns={
        'T_CANOPY':'T00cm',
        'T_1_1_1':'T05cm',
        'T_1_2_1':'T10cm',
        'T_1_3_1':'T20cm',
        'T_1_4_1':'T30cm',
        'T_1_5_1':'T40cm',
        'T_1_6_1':'T50cm',
        'T_1_7_1':'T60cm', 
        'VWC_1_1_1':'VWC05cm',
        'VWC_1_2_1':'VWC10cm',
        'VWC_1_3_1':'VWC20cm',
        'VWC_1_4_1':'VWC30cm',
        'VWC_1_5_1':'VWC40cm',
        'VWC_1_6_1':'VWC50cm',
        'VWC_1_7_1':'VWC60cm',
    }).set_index('datetime')
    return (df_grad,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    Include data from the SI-111, which is technically measuring Temperature of the canopy top, but is T0cm for the timeframe being examined.
    """)
    return


@app.cell
def _(df_grad, mo, px):
    _plot = px.line(
        df_grad.reset_index(), x="datetime", y="T00cm"
    )

    plot = mo.ui.plotly(_plot)
    return


@app.cell
def _(df_grad):
    df_g = df_grad.reset_index()
    return (df_g,)


@app.cell
def _(df_g, mo):
    # Get all columns except the date column
    data_columns = [col for col in df_g.columns if "cm" in col and "T" in col]

    # Create a multiselect widget for column selection
    column_selector = mo.ui.multiselect(
        options=data_columns,
        value=data_columns[:2],  # Default to first 2 columns
        label="Select columns to plot:"
    )
    column_selector
    return (column_selector,)


@app.cell
def _(column_selector, df_g, go, mo):
    # Create the plotly figure reactively based on selected columns
    selected_cols = column_selector.value

    fig = go.Figure()

    for col in selected_cols:
        fig.add_trace(go.Scatter(
            x=df_g['datetime'],
            y=df_g[col],
            mode='lines',
            name=col
        ))

    fig.update_layout(
        title='Timeseries Data',
        xaxis_title='Date',
        yaxis_title='Temperature',
        hovermode='x unified',
        height=500
    )

    mo.ui.plotly(fig)
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## Initial Estimate of G0
    """)
    return


@app.cell
def _(df_grad, soil_heat):
    # ---------- Parameters ------------------------------------------------------
    DEPTH_TOP  = 0.05   # m  (sensor 1)
    DEPTH_LOW  = 0.10   # m  (sensor 2)

    # ---- Call helper – returns a pd.Series (W m‑2, positive downward)
    G0_est = soil_heat.compute_heat_flux_conduction(
        df_grad,
        depth1=DEPTH_TOP,
        depth2=DEPTH_LOW,
        col_T1='T05cm',
        col_T2='T10cm',
        col_theta1='VWC05cm',
        col_theta2='VWC10cm',
        porosity = 0.5,
        k_dry = 0.45,
        k_sat = 2
    )

    G0_est_2 = soil_heat.compute_heat_flux_calorimetric(df_grad,depth_levels=[0.05,0.1,0.2,0.3],T_cols=['T05cm', 'T10cm','T20cm','T30cm'],theta_cols=['VWC05cm', 'VWC10cm', 'VWC20cm', 'VWC30cm'], )
    return G0_est, G0_est_2


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## Calculate amplitude and lag of data
    """)
    return


@app.cell
def _(df_grad, pd, soil_heat):
    amp_values = pd.concat([soil_heat.diurnal_amplitude(df_grad['T00cm']),
                            soil_heat.diurnal_amplitude(df_grad['T05cm']),
                            soil_heat.diurnal_amplitude(df_grad['T10cm']),
                            soil_heat.diurnal_amplitude(df_grad['T20cm']),
                            soil_heat.diurnal_amplitude(df_grad['T30cm']),
    ],axis=1).round(2).rename(columns={
        'T00cm':'T00cm_amp',
        'T05cm':'T05cm_amp',
        'T10cm':'T10cm_amp',
        'T20cm':'T20cm_amp',
        'T30cm':'T30cm_amp',
    })

    amp_values
    return (amp_values,)


@app.cell
def _(amp_values, np, pd, soil_heat):
    diff_00_10 = soil_heat.thermal_diffusivity_amplitude(amp_values['T00cm_amp'], amp_values['T10cm_amp'],z1=0.00,z2=0.10)#.median()#.plot(label='T00cm to T10cm')
    diff_00_05 = soil_heat.thermal_diffusivity_amplitude(amp_values['T00cm_amp'], amp_values['T05cm_amp'],z1=0.00,z2=0.05)#.median()#.plot(label='T00cm to T05cm')
    diff_05_10 = soil_heat.thermal_diffusivity_amplitude(amp_values['T05cm_amp'], amp_values['T10cm_amp'],z1=0.05,z2=0.10)#.median()#.plot(label='T05cm to T10cm')
    diff = pd.concat([diff_00_10, diff_00_05, diff_05_10], axis=1)#.plot(title='Thermal diffusivity (amplitude method)', ylabel='m² s⁻¹')
    diff['median'] = diff.median(axis=1)
    #soil_heat.thermal_diffusivity_amplitude(amp_values['T05cm_amp'], amp_values['T20cm_amp'],z1=0.00,z2=0.20).plot(label='T05cm to T20cm')
    #soil_heat.thermal_diffusivity_amplitude(amp_values['T05cm_amp'], amp_values['T30cm_amp'],z1=0.05,z2=0.30).plot()
    #soil_heat.thermal_diffusivity_amplitude(amp_values['T10cm_amp'], amp_values['T20cm_amp'],z1=0.10,z2=0.20).plot()
    print(f"Thermal diffusivity (T00cm to T10cm): {diff_00_10.median():.7f} m² s⁻¹")
    print(f"Thermal diffusivity (T00cm to T05cm): {diff_00_05.median():.7f} m² s⁻¹")
    print(f"Thermal diffusivity (T05cm to T10cm): {diff_05_10.median():.7f} m² s⁻¹")
    print(f"Mean Thermal diffusivity (T00cm to T10cm): {np.mean([diff_00_10.median(),diff_00_05.median(),diff_05_10.median()]):.7f} m² s⁻¹")
    return (diff,)


@app.cell
def _(df_grad, diff, plt):
    diff['median'].plot(title='Thermal diffusivity (amplitude method)', ylabel='m² s⁻¹')
    plt.twinx()
    df_grad['VWC05cm'].resample('D').mean().plot(color='red')
    return


@app.cell
def _(np, pd):
    from soil_heat import diurnal_amplitude

    def estimate_rhoc_dry(alpha: pd.Series, theta: pd.Series, porosity: float=0.4, k_dry: float=0.25, k_sat: float=1.5, rhoc_w: float=4180000.0, dry_quantile: float=0.1) -> float:
        """Return volumetric heat capacity of *dry* soil (J m⁻³ K⁻¹)."""
        theta, alpha = theta.align(alpha, join='inner')
        frac_sat = np.clip(theta / porosity, 0.0, 1.0)
        lam = k_dry + (k_sat - k_dry) * frac_sat
        Cv = lam / alpha
        rhoc_dry = (Cv - theta * rhoc_w) / (1.0 - theta)
        dry_days = theta <= theta.quantile(dry_quantile)
        return float(rhoc_dry.loc[dry_days].median())  # keep only days where both alpha & theta are available  # ⬅ key line  # W m⁻¹ K⁻¹  # --- 3. Heat capacity & dry-soil estimate ------------------  # J m⁻³ K⁻¹
    return (estimate_rhoc_dry,)


@app.cell
def _(df_grad, diff, estimate_rhoc_dry):
    alpha_daily = diff['median']
    theta_daily = df_grad['VWC05cm'].resample('D').mean()

    estimate_rhoc_dry(alpha_daily, theta_daily, porosity=0.40, k_dry=0.25, k_sat=1.50, rhoc_w=4.18e6, dry_quantile=0.10)*1e-6
    return


@app.cell
def _(df_grad, soil_heat):
    G0 = soil_heat.compute_heat_flux_calorimetric(
        df_grad,
        depth_levels = [0.05, 0.10, 0.20, 0.30],
        T_cols = ['T05cm', 'T10cm', 'T20cm', 'T30cm'],
        theta_cols = ['VWC05cm', 'VWC10cm', 'VWC20cm', 'VWC30cm'],
        C_dry = 1.72e6,
        C_w = 4.2e6,)
    G0.plot()
    return (G0,)


@app.cell
def _(df_grad, pd, soil_heat):
    lag_values = pd.concat([soil_heat.diurnal_peak_lag(df_grad['T05cm'], df_grad['T10cm']),
                           soil_heat.diurnal_peak_lag(df_grad['T05cm'], df_grad['T20cm']),
                           soil_heat.diurnal_peak_lag(df_grad['T05cm'], df_grad['T30cm']),],axis=1)
    lag_values
    return


@app.cell
def _(G0_est, G0_est_2):
    G0_est.plot(label='G0_est_grad', color='red')
    G0_est_2.plot(label='G0_est_calorimetric', color='blue')
    return


@app.cell
def _(pd):
    # ---------- Eddy‑covariance flux data ---------------------------------------
    # ---------- SoilVUE data ----------------------------------------------------
    flux_path_1 = "21314_Flux_AmeriFluxFormat_21.dat"
    flux_path_2 = "21314_Flux_AmeriFluxFormat_22.dat"

    df_flux_1 = pd.read_csv(flux_path_1, na_values=['NAN'])
    df_flux_2 = pd.read_csv(flux_path_2, na_values=['NAN'])
    df_flux = pd.concat([df_flux_1, df_flux_2], ignore_index=True)
    # Parse timestamp and keep only ground‑heat flux variable ‘G’
    df_flux['datetime'] = pd.to_datetime(df_flux['TIMESTAMP_START'].astype(str), format='%Y%m%d%H%M.0')

    df_flux = df_flux[['datetime', 'G']].dropna(subset=['G']).astype({'G': float})

    # Sub‑set to 11 July 2024 (matches the SoilVUE example file)
    #mask_0711 = df_flux['datetime'].dt.strftime('%Y%m%d') == '20240711'
    #df_flux_0711 = df_flux.loc[mask_0711].set_index('datetime')
    #df_flux_0711.head()
    df_flux
    return (df_flux,)


@app.cell
def _(G0, df_grad):
    df_grad['NETRAD'].plot(label='NETRAD', color='orange')
    G0.plot(label='G0_est_grad', color='red')
    return


@app.cell
def _(G0, df_grad, pd):
    df_NETRAD_G0 = pd.concat([df_grad['NETRAD'], G0], axis=1).dropna()

    x = df_NETRAD_G0['NETRAD']
    y = df_NETRAD_G0['G_calorimetric']

    import statsmodels.api as sm

    X = sm.add_constant(x)      # adds intercept term
    model = sm.OLS(y, X).fit()  # fit the model
    print(model.summary())
    return


@app.cell
def _(G0, df_flux):
    df_flux.plot(x='datetime', y='G', label='G (eddy‑covariance)', color='blue')
    G0.plot(label='G0_est_grad', color='red')
    return


@app.cell
def _(np, pd):
    def calculate_soil_heat_storage(df, depths, porosity=0.45, bulk_density=1.3):
        """
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
        """
    
        # Constants
        Cw = 4.18e6  # J m-3 K-1, volumetric heat capacity of water
        Cs = 2.0e6   # J m-3 K-1, volumetric heat capacity of soil minerals
        Ca = 1.2e3   # J m-3 K-1, volumetric heat capacity of air (often negligible)
    
        # Define layer boundaries (midpoints between sensors)
        # First layer: 0 to midpoint between surface and first sensor
        # Last layer: from midpoint to some assumed bottom depth
    
        layer_boundaries = [0]  # Start at surface
        for i in range(len(depths)-1):
            midpoint = (depths[i] + depths[i+1]) / 2
            layer_boundaries.append(midpoint)
        # Assume bottom boundary extends another half-interval
        last_interval = depths[-1] - depths[-2]
        layer_boundaries.append(depths[-1] + last_interval/2)
    
        # Convert to meters
        layer_boundaries = np.array(layer_boundaries) / 100
        depths_m = np.array(depths) / 100
    
        # Calculate layer thicknesses
        layer_thicknesses = np.diff(layer_boundaries)
    
        print(f"Layer boundaries (m): {layer_boundaries}")
        print(f"Layer thicknesses (m): {layer_thicknesses}")
        print(f"Total depth (m): {layer_boundaries[-1]}")
    
        # Initialize result dataframe
        result = pd.DataFrame(index=df.index)
    
        # Calculate volumetric heat capacity for each layer
        for i, depth in enumerate(depths):
            # Get temperature and soil moisture
            T_col = f'T_{depth}cm'
            SM_col = f'SM_{depth}cm'
        
            # Soil moisture (convert to 0-1 if given as percentage)
            theta_w = df[SM_col].copy()
            if theta_w.max() > 1:
                theta_w = theta_w / 100
        
            # Volume fractions
            theta_s = bulk_density / 2.65  # Assuming particle density of 2.65 g/cm³
            theta_a = porosity - theta_w
            theta_a = theta_a.clip(lower=0)  # Can't be negative
        
            # Volumetric heat capacity (J m-3 K-1)
            Cv = theta_w * Cw + theta_s * Cs + theta_a * Ca
        
            # Heat content for this layer (J m-2)
            dz = layer_thicknesses[i]
            heat_content = Cv * dz * df[T_col]
        
            result[f'Cv_{depth}cm'] = Cv
            result[f'heat_content_{depth}cm'] = heat_content
    
        # Total heat storage (sum across all layers)
        heat_cols = [c for c in result.columns if c.startswith('heat_content_')]
        result['total_heat_storage'] = result[heat_cols].sum(axis=1)
    
        # Calculate heat storage change (soil heat flux, G)
        # This is dS/dt - typically calculated as change from some reference time
        result['heat_storage_change'] = result['total_heat_storage'].diff()  # J m-2 per timestep
    
        # Convert to W m-2 (assuming hourly data)
        dt_seconds = (df.index[1] - df.index[0]).total_seconds()
        result['G'] = result['heat_storage_change'] / dt_seconds  # W m-2
    
        return result
    return


if __name__ == "__main__":
    app.run()
