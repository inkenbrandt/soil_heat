=====
Usage
=====

To use Soil Heat in a project::

    import soil_heat

FAO-56 Soil Heat Flux
---------------------

The `fao56_soil_heat_flux` module implements methods from the FAO-56 paper and the ASCE Standardized Reference ET equation.

Daily Soil Heat Flux
^^^^^^^^^^^^^^^^^^^^

For daily time steps, soil heat flux is often assumed to be zero for a grass reference surface:

.. code-block:: python

    from soil_heat import soil_heat_flux_daily

    g_daily = soil_heat_flux_daily()
    print(f"Daily G: {g_daily} MJ m-2 day-1")

Monthly Soil Heat Flux
^^^^^^^^^^^^^^^^^^^^^^

When both the previous and next month's mean air temperatures are available:

.. code-block:: python

    from soil_heat import soil_heat_flux_monthly

    # Temps for March and May
    t_mar, t_may = 14.1, 18.8
    g_apr = soil_heat_flux_monthly(t_mar, t_may)
    print(f"April G: {g_apr:.3f} MJ m-2 day-1")

Hourly Soil Heat Flux
^^^^^^^^^^^^^^^^^^^^^

For hourly or sub-daily periods, G is estimated as a fraction of net radiation (Rn):

.. code-block:: python

    from soil_heat import soil_heat_flux_hourly

    rn = 2.5  # MJ m-2 hr-1
    g_hourly = soil_heat_flux_hourly(rn, is_daytime=True)
    print(f"Hourly G: {g_hourly} MJ m-2 hr-1")

Storage Calculations
--------------------

The `storage_calculations` module provides methods to estimate energy storage in the soil and canopy.

Soil Heat Storage
^^^^^^^^^^^^^^^^^

Using the change-in-storage approach for a soil layer:

.. code-block:: python

    import pandas as pd
    import numpy as np
    from soil_heat import compute_soil_storage_integrated

    # Create sample data
    index = pd.date_range("2023-01-01", periods=2, freq="15T")
    df = pd.DataFrame({
        "T5cm": [10.0, 10.5],
        "VWC5cm": [0.25, 0.26]
    }, index=index)

    # Compute storage flux (W/m2)
    s_soil = compute_soil_storage_integrated(df, col_T5="T5cm", col_VWC5="VWC5cm")
    print(s_soil)

Canopy Heat Storage
^^^^^^^^^^^^^^^^^^^

Estimating heat storage within biomass:

.. code-block:: python

    from soil_heat import compute_canopy_storage

    # Sample data with canopy temperature
    df["T_CANOPY"] = [15.0, 16.0]

    s_canopy = compute_canopy_storage(df, col_irt="T_CANOPY", crop_type="alfalfa")
    print(s_canopy)

Soil Thermal Properties and Advanced Flux (Johansen)
----------------------------------------------------

The `johansen` module provides sophisticated models for soil thermal properties and heat flux.

Thermal Properties
^^^^^^^^^^^^^^^^^^

Estimating volumetric heat capacity and thermal conductivity from water content:

.. code-block:: python

    import numpy as np
    from soil_heat import volumetric_heat_capacity, thermal_conductivity_johansen

    theta = 0.25  # m3/m3
    cv = volumetric_heat_capacity(theta)
    lam = thermal_conductivity_johansen(theta)

    print(f"Heat Capacity: {cv:.2e} J m-3 K-1")
    print(f"Conductivity: {lam:.3f} W m-1 K-1")

Surface Heat Flux from Gradient and Storage
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The `method_gradient_plus_storage` function estimates surface ground heat flux (G0) by combining a conductive flux at a reference depth with calorimetric storage in the layers above:

.. code-block:: python

    import pandas as pd
    from soil_heat import method_gradient_plus_storage, DEPTHS_CM

    # Prepare sample temperature and water content dataframes
    index = pd.date_range("2023-01-01", periods=10, freq="h")
    ts = pd.DataFrame({d: np.linspace(10, 15, 10) for d in DEPTHS_CM}, index=index)
    swc = pd.DataFrame({d: [0.2]*10 for d in DEPTHS_CM}, index=index)

    # Estimate surface G0 using 10 cm as the reference depth
    g0 = method_gradient_plus_storage(ts, swc, ref_depth_idx=2)
    print(g0.head())
