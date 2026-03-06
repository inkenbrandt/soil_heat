import numpy as np
import pandas as pd
import pytest
from soil_heat import (
    soil_heat_flux_general,
    soil_heat_flux_monthly,
    soil_heat_flux_hourly,
    compute_soil_storage_integrated,
    compute_canopy_storage
)

def test_fao56_general():
    # FAO-56 Example 13 — March to April
    g = soil_heat_flux_general(16.1, 14.1, delta_t=30.0, delta_z=1.0)
    assert np.isclose(g, 0.14)

def test_fao56_monthly():
    # FAO-56 Example 13 — April soil heat flux (March = 14.1, May = 18.8)
    g = soil_heat_flux_monthly(14.1, 18.8)
    assert np.isclose(g, 0.329)

def test_fao56_hourly():
    # Daytime hour with Rn = 2.5 MJ m-2 hr-1
    g = soil_heat_flux_hourly(2.5, is_daytime=True)
    assert np.isclose(g, 0.25)

def test_compute_soil_storage_integrated():
    index = pd.date_range("2023-01-01", periods=2, freq="15min")
    df = pd.DataFrame({
        "T5cm": [10.0, 10.5],
        "VWC5cm": [0.25, 0.25]
    }, index=index)
    s_soil = compute_soil_storage_integrated(df, col_T5="T5cm", col_VWC5="VWC5cm")
    assert isinstance(s_soil, pd.Series)
    assert not s_soil.isna().all()

def test_compute_canopy_storage():
    index = pd.date_range("2023-01-01", periods=2, freq="15min")
    df = pd.DataFrame({
        "T_CANOPY": [15.0, 16.0]
    }, index=index)
    s_canopy = compute_canopy_storage(df, col_irt="T_CANOPY", crop_type="alfalfa")
    assert isinstance(s_canopy, pd.Series)
    assert s_canopy.iloc[1] > 0
