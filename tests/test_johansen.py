import numpy as np
import pandas as pd
import pytest
from soil_heat.johansen import (
    volumetric_heat_capacity,
    thermal_conductivity_johansen,
    thermal_conductivity_lu2007,
    gradient_surface,
    compute_storage_above_plate,
    gradient_plus_storage,
    DEPTHS_CM
)

def test_volumetric_heat_capacity():
    # Test scalar input
    theta = 0.2
    cv = volumetric_heat_capacity(theta)
    assert isinstance(cv, float)
    assert cv > 0

    # Test array input
    theta_arr = np.array([0.1, 0.2, 0.3])
    cv_arr = volumetric_heat_capacity(theta_arr)
    assert isinstance(cv_arr, np.ndarray)
    assert len(cv_arr) == 3
    assert np.all(cv_arr > 0)

def test_thermal_conductivity_johansen():
    # Test scalar input
    theta = 0.2
    lam = thermal_conductivity_johansen(theta)
    assert isinstance(lam, float)
    assert lam > 0

    # Test array input
    theta_arr = np.array([0.1, 0.2, 0.3])
    lam_arr = thermal_conductivity_johansen(theta_arr)
    assert isinstance(lam_arr, np.ndarray)
    assert len(lam_arr) == 3
    assert np.all(lam_arr > 0)

def test_thermal_conductivity_lu2007():
    # Test scalar input
    theta = 0.2
    lam = thermal_conductivity_lu2007(theta)
    assert isinstance(lam, float)
    assert lam > 0

    # Test array input
    theta_arr = np.array([0.1, 0.2, 0.3])
    lam_arr = thermal_conductivity_lu2007(theta_arr)
    assert isinstance(lam_arr, np.ndarray)
    assert len(lam_arr) == 3
    assert np.all(lam_arr > 0)

def test_gradient_surface():
    index = pd.date_range("2023-01-01", periods=5, freq="h")
    ts = pd.DataFrame({5.0: [10.0, 11.0, 12.0, 11.0, 10.0]}, index=index)
    swc = pd.DataFrame({5.0: [0.2, 0.2, 0.2, 0.2, 0.2]}, index=index)
    df = pd.DataFrame({"T_CANOPY_1_1_1": [15.0, 16.0, 17.0, 16.0, 15.0]}, index=index)

    g_surf = gradient_surface(ts, swc, df)
    assert isinstance(g_surf, pd.Series)
    assert len(g_surf) == 5
    assert g_surf.name == "G_gradient_surface"

def test_compute_storage_above_plate():
    index = pd.date_range("2023-01-01", periods=5, freq="h")
    ts = pd.DataFrame({5.0: [10.0, 11.0, 12.0, 11.0, 10.0]}, index=index)
    swc = pd.DataFrame({5.0: [0.2, 0.2, 0.2, 0.2, 0.2]}, index=index)

    s_plate = compute_storage_above_plate(ts, swc)
    assert isinstance(s_plate, pd.Series)
    assert len(s_plate) == 5
    assert s_plate.name == "S_above_plate"
    assert pd.isna(s_plate.iloc[0])

def test_gradient_plus_storage():
    index = pd.date_range("2023-01-01", periods=5, freq="h")
    ts_data = {d: np.linspace(10, 15, 5) for d in DEPTHS_CM}
    ts = pd.DataFrame(ts_data, index=index)
    swc_data = {d: [0.2]*5 for d in DEPTHS_CM}
    swc = pd.DataFrame(swc_data, index=index)

    g_stor = gradient_plus_storage(ts, swc, ref_depth_idx=2)
    assert isinstance(g_stor, pd.Series)
    assert len(g_stor) == 5
    assert g_stor.name == "G_grad+stor_20cm"

    # Test with Lu2007 model
    g_stor_lu = gradient_plus_storage(ts, swc, ref_depth_idx=2, lam_model="lu2007")
    assert isinstance(g_stor_lu, pd.Series)
