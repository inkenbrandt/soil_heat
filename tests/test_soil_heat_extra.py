import numpy as np
import pandas as pd
import pytest
from soil_heat import soil_heat

@pytest.fixture
def sample_dataframe():
    """Create a sample dataframe for testing."""
    dates = pd.to_datetime(
        [
            "2023-01-01 00:00",
            "2023-01-01 00:30",
            "2023-01-01 01:00",
            "2023-01-01 01:30",
        ]
    )
    data = {
        "T5cm": [10.0, 10.5, 11.0, 11.5],
        "T10cm": [9.5, 10.0, 10.5, 11.0],
        "VWC5cm": [0.2, 0.21, 0.22, 0.23],
        "VWC10cm": [0.25, 0.26, 0.27, 0.28],
    }
    df = pd.DataFrame(data, index=dates)
    return df

def test_temperature_gradient():
    """Test the temperature_gradient function."""
    # Test with scalar inputs
    grad = soil_heat.temperature_gradient(
        T_upper=18.6, T_lower=20.1, depth_upper=0.05, depth_lower=0.10
    )
    assert np.isclose(grad, 30.0)

    # Test with array inputs
    T_up = np.array([15.0, 16.2, 17.1])
    T_low = np.array([14.0, 15.8, 16.9])
    grad_arr = soil_heat.temperature_gradient(T_up, T_low, 0.02, 0.08)
    expected = np.array([-16.66666667, -6.66666667, -3.33333333])
    assert np.allclose(grad_arr, expected)

    # Test with invalid depth inputs
    with pytest.raises(ValueError):
        soil_heat.temperature_gradient(18.6, 20.1, 0.10, 0.05)
    with pytest.raises(ValueError):
        soil_heat.temperature_gradient(18.6, 20.1, 0.10, 0.10)

def test_soil_heat_flux():
    """Test the soil_heat_flux function."""
    # Test with scalar inputs
    flux = soil_heat.soil_heat_flux(
        T_upper=18.6, T_lower=20.1, depth_upper=0.05, depth_lower=0.10, k=0.25
    )
    assert np.isclose(flux, -7.5)

    # Test with array inputs
    T_up = np.array([15.0, 16.2, 17.1])
    T_low = np.array([14.0, 15.8, 16.9])
    k = np.array([0.2, 0.3, 0.4])
    flux_arr = soil_heat.soil_heat_flux(T_up, T_low, 0.02, 0.08, k)
    expected = np.array([3.33333333, 2.0, 1.33333333])
    assert np.allclose(flux_arr, expected)

def test_volumetric_heat_capacity():
    """Test the volumetric_heat_capacity function."""
    # Test with scalar input
    cv = soil_heat.volumetric_heat_capacity(0.2)
    assert np.isclose(cv, 2390.8)

    # Test with array input
    theta_v = np.array([0.1, 0.2, 0.3])
    cv_arr = soil_heat.volumetric_heat_capacity(theta_v)
    expected = np.array([2166.4, 2390.8, 2615.2])
    assert np.allclose(cv_arr, expected)

def test_thermal_conductivity():
    """Test the thermal_conductivity function."""
    # Test with scalar inputs
    k = soil_heat.thermal_conductivity(1.5e-7, 0.2)
    assert np.isclose(k, 0.00035862)

    # Test with array inputs
    alpha = np.array([1.4e-7, 1.6e-7])
    theta = np.array([0.10, 0.25])
    k_arr = soil_heat.thermal_conductivity(alpha, theta)
    expected = np.array([0.000303296, 0.000410048])
    assert np.allclose(k_arr, expected, atol=1e-5)

def test_sinusoid():
    """Test the sinusoid function."""
    t = np.linspace(0, 24, 1000)
    temp = soil_heat.sinusoid(t, A=6, omega=2 * np.pi / 24, phase=0, offset=15)
    assert np.isclose(temp[0], 15.0)
    assert np.isclose(temp[250], 21.0)
    assert np.isclose(temp[500], 15.0, atol=0.02)
    assert np.isclose(temp[750], 9.0, atol=0.02)

def test_compute_heat_flux_conduction(sample_dataframe):
    """Test the compute_heat_flux_conduction function."""
    G = soil_heat.compute_heat_flux_conduction(sample_dataframe)
    assert isinstance(G, pd.Series)
    assert G.name == "G_conduction"
    assert G.shape == (4,)
    expected = pd.Series(
        [9.53125, 9.84375, 10.15625, 10.46875],
        index=sample_dataframe.index,
        name="G_conduction",
    )
    pd.testing.assert_series_equal(G, expected)

def test_compute_heat_flux_calorimetric(sample_dataframe):
    depths = [0.05, 0.10]
    Tcols = ["T5cm", "T10cm"]
    Vcols = ["VWC5cm", "VWC10cm"]
    G0 = soil_heat.compute_heat_flux_calorimetric(sample_dataframe, depths, Tcols, Vcols)
    assert np.isnan(G0.iloc[0])
    assert np.isclose(G0.iloc[1], 35.875, atol=1e-3)

def test_diurnal_amplitude(sample_dataframe):
    amp = soil_heat.diurnal_amplitude(sample_dataframe["T5cm"])
    assert amp.iloc[0] == 1.5

def test_diurnal_peak_lag(sample_dataframe):
    lag = soil_heat.diurnal_peak_lag(sample_dataframe["T5cm"], sample_dataframe["T10cm"])
    assert lag.iloc[0] == 0.0

def test_fit_sinusoid():
    t = np.linspace(0, 86400, 100)
    y = 10 * np.sin(2 * np.pi / 86400 * t + 0.5) + 20
    popt = soil_heat.fit_sinusoid(t, y)
    assert np.isclose(popt[0], 10, rtol=1e-2)
    assert np.isclose(popt[3], 20, rtol=1e-2)

def test_thermal_diffusivity_amplitude():
    alpha = soil_heat.thermal_diffusivity_amplitude(10, 5, 0.05, 0.10)
    assert np.isclose(alpha, 1.892e-7, atol=1e-10)

def test_thermal_diffusivity_lag():
    alpha = soil_heat.thermal_diffusivity_lag(3600, 0.05, 0.10)
    assert np.isclose(alpha, 1.326e-6, atol=1e-9)

def test_thermal_diffusivity_logrithmic():
    alpha = soil_heat.thermal_diffusivity_logrithmic(
        20, 22, 20, 18, 18, 19, 18, 17, 0.05, 0.10
    )
    assert np.isclose(alpha, 1.892e-7, atol=1e-10)

def test_calculate_thermal_properties_for_all_pairs(sample_dataframe):
    df = sample_dataframe.copy()
    df["T5swc"] = 20.0
    df["T10swc"] = 25.0
    df = df.rename(columns={"T5cm": "T5ts", "T10cm": "T10ts"})
    mapping = {"T5ts": 0.05, "T10ts": 0.10}
    res = soil_heat.calculate_thermal_properties_for_all_pairs(df, mapping)
    assert not res.empty

def test_estimate_rhoc_dry():
    alpha = pd.Series([1e-7, 1.1e-7], index=[0, 1])
    theta = pd.Series([0.1, 0.12], index=[0, 1])
    rhoc = soil_heat.estimate_rhoc_dry(alpha, theta)
    assert rhoc > 0

def test_calculate_soil_heat_storage(sample_dataframe):
    df = sample_dataframe.copy()
    df = df.rename(columns={"T5cm": "T_5cm", "VWC5cm": "SM_5cm", "T10cm": "T_10cm", "VWC10cm": "SM_10cm"})
    res = soil_heat.calculate_soil_heat_storage(df, [5, 10])
    assert "G" in res
