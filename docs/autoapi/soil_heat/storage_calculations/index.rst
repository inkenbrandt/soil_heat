soil_heat.storage_calculations
==============================

.. py:module:: soil_heat.storage_calculations


Functions
---------

.. autoapisummary::

   soil_heat.storage_calculations.compute_soil_storage_integrated
   soil_heat.storage_calculations.compute_canopy_storage


Module Contents
---------------

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
