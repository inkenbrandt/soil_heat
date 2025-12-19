import pandas as pd
import numpy as np
import logging

def compute_soil_storage_integrated(
    df: pd.DataFrame,
    col_T5: str = "T5cm",
    col_VWC5: str = "VWC5cm",
    col_IRT: str = None, # must be SURFACE temp, not vegetation temperature
    depth: float = 0.05,
    bulk_density: float = 1300.0, # value in Campbell
    specific_heat_dry: float = 870.0, # value in Campbell
    specific_heat_water: float = 4218.0 # value in Campbell
) -> pd.Series:
    r"""
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

    Parameters
    ----------
    df : pd.DataFrame
        Input dataframe with a DatetimeIndex.
    col_T5 : str, default "T5cm"
        Column name for soil temperature at 5cm depth [°C].
    col_VWC5 : str, default "VWC5cm"
        Column name for Volumetric Water Content at 5cm depth [m³/m³].
    col_IRT : str, optional
        Column name for surface infrared temperature [°C]. If provided, 
        the layer temperature is the average of surface and 5cm temps. Note 
        that this must be surface temperature, not canopy temperature
    depth : float, default 0.05
        The thickness of the soil layer being measured [m].
    bulk_density : float, default 1300.0
        Dry bulk density of the soil [kg/m³].
    specific_heat_dry : float, default 870.0
        Specific heat capacity of dry soil minerals [J/(kg·K)].
    specific_heat_water : float, default 4218.0
        Specific heat capacity of water [J/(kg·K)].

    Returns
    -------
    pd.Series
        The calculated soil heat storage flux [W/m²]. 
        Positive values indicate energy moving into storage (warming).

    Notes
    -----
    The canopy storage flux is calculated as:
    $$S_{canopy} = \frac{m_{veg} \cdot C_{p} \cdot \Delta T}{\Delta t}$$
    
    Where:
    - $m_{veg}$ = Fresh biomass ($kg/m^2$), calculated as $(height / 100) \cdot density\_factor$
    - $C_{p}$ = Specific heat of fresh vegetation ($\approx 3900 \text{ J/kg}\cdot\text{K}$)
    - $\Delta T$ = Change in canopy temperature between time steps
    - $\Delta t$ = Time step in seconds
    """

    rho_w = 1000.0
    dt = df.index.to_series().diff().dt.total_seconds()

    # 1. Determine the representative temperature for the 0-5cm layer
    if col_IRT and col_IRT in df.columns:
        # Integrated approach: Average of surface and 5cm depth
        T_layer = (df[col_IRT] + df[col_T5]) / 2.0
        method_label = "IRT-Adjusted"
    else:
        # Standard approach: 5cm sensor represents the whole slab
        T_layer = df[col_T5]
        method_label = "Slab-Only"

    T_layer_prev = T_layer.shift(1)
    VWC_curr = df[col_VWC5]
    VWC_prev = df[col_VWC5].shift(1)

    # 2. Dry Soil Energy Change (J/m^3)
    # (T_curr - T_prev) * Cds * rho_b
    dry_storage = (T_layer - T_layer_prev) * specific_heat_dry * bulk_density

    # 3. Water Energy Change (J/m^3)
    # (T_curr * VWC_curr - T_prev * VWC_prev) * rho_w * Cw
    water_storage = (T_layer * VWC_curr - T_layer_prev * VWC_prev) * rho_w * specific_heat_water

    # 4. Final Flux Conversion (W/m^2)
    storage = (dry_storage + water_storage) * depth / dt

    # Repl
    storage = storage.replace([np.inf, -np.inf], np.nan)
    
    return pd.Series(storage, index=df.index, name=f"S_{method_label}")


def compute_canopy_storage(
    df: pd.DataFrame, 
    col_irt: str = "T_CANOPY", 
    crop_type: str = "alfalfa",
    height_cm: str | float = None,
    density_factor: float = None
) -> pd.Series:
    """
    Calculate the energy storage flux of a plant canopy (S_canopy) in W/m².

    This function estimates heat storage within biomass. It can handle dynamic 
    vegetation growth if a height column is provided, otherwise it falls back 
    to static heights. Height and density values can be entered by the user or 
    taken from a dictionary within the function.

    Care should be taken to examine the crop default values and the specific 
    heat of the vegetation. 

    Parameters
    ----------
    df : pd.DataFrame
        Input dataframe with a DatetimeIndex.
    col_irt : str, default "T_CANOPY"
        Column name for the canopy surface temperature [°C].
    crop_type : str, default "alfalfa"
        Options: 'alfalfa', 'corn', 'hay', 'wetland'. Used to determine 
        defaults for height and density.
    height_cm : str or float, optional
        If str: The column name in `df` containing canopy height [cm] over time.
        If float: A static height value [cm].
        If None: Uses default height for the specified `crop_type`.
    density_factor : float, optional
        Biomass density [kg/m² per meter height]. If None, uses default 
        for the specified `crop_type`.

    Returns
    -------
    pd.Series
        Energy storage flux [W/m²].
    """
    
    # 1. Define crop-specific defaults
    crop_defaults = {
        "alfalfa": {"height": 30.0, "density": 5.5},
        "corn":    {"height": 200.0, "density": 4.5},
        "hay":     {"height": 30.0, "density": 3.5},
        "wetland": {"height": 120.0, "density": 7.0}
    }

    crop_key = crop_type.lower()
    if crop_key not in crop_defaults:
        logging.warning(
            f"Crop type '{crop_type}' not found in defaults. "
            f"Falling back to 'alfalfa' settings (H: 40cm, D: 5.5). "
            f"Available crops: {list(crop_defaults.keys())}"
        )
        defaults = crop_defaults["alfalfa"]
    else:
        defaults = crop_defaults[crop_key]
    defaults = crop_defaults.get(crop_type.lower(), crop_defaults["alfalfa"])

    # 2. Resolve Height (Column vs. Scalar vs. Default)
    if height_cm is None:
        h_val = defaults["height"]
    elif isinstance(height_cm, str) and height_cm in df.columns:
        h_val = df[height_cm]
    else:
        h_val = height_cm  # Assume it's a numeric scalar

    # 3. Resolve Density
    final_d = density_factor if density_factor is not None else defaults["density"]

    # 4. Physical Constants
    cp_veg = 3900 # Specific heat of fresh vegetation (J/kg·K)
    
    # 5. Calculate dynamic fresh biomass (kg/m2)
    # This results in a Series if h_val is a Series, or a scalar if h_val is a float
    m_veg = (h_val / 100.0) * final_d
    
    # 6. Calculate flux
    dt = df.index.to_series().diff().dt.total_seconds()
    dT = df[col_irt].diff()
    
    # S = (m * Cp * dT) / dt
    s_canopy = (m_veg * cp_veg * dT) / dt
    
    return s_canopy.replace([np.inf, -np.inf], 0).fillna(0)

