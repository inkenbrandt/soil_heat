import numpy as np
import pytest

from soil_heat import gao_et_al


def test_lambda_s():
    """Test the lambda_s function."""
    # Test with a scalar value
    assert np.isclose(gao_et_al.lambda_s(0.25), 1.076866957449353)

    # Test with a numpy array
    theta = np.linspace(0, 0.5, 5)
    expected = 0.20 + np.exp(1.46 * (theta - 0.34))
    assert np.allclose(gao_et_al.lambda_s(theta), expected)

    # Test with invalid input
    with pytest.raises(ValueError):
        gao_et_al.lambda_s(-0.1)
    with pytest.raises(ValueError):
        gao_et_al.lambda_s(1.1)


def test_k_s():
    """Test the k_s function."""
    # Test with a scalar value
    assert np.isclose(gao_et_al.k_s(0.25), 1.6598e-07, atol=1e-9)

    # Test with a numpy array
    thetas = np.linspace(0, 0.5, 6)
    expected = (0.69 + np.exp(3.06 * (thetas - 0.26))) * 1e-7
    assert np.allclose(gao_et_al.k_s(thetas), expected)

    # Test with invalid input
    with pytest.raises(ValueError):
        gao_et_al.k_s(-0.1)
    with pytest.raises(ValueError):
        gao_et_al.k_s(1.1)


def test_volumetric_heat_capacity():
    """Test the volumetric_heat_capacity function."""
    # Test with scalar inputs
    assert np.isclose(
        gao_et_al.volumetric_heat_capacity(0.8, 1.2e-6), 666666.666
    )

    # Test with vectorized inputs
    lambda_vals = np.array([0.5, 0.6, 0.7])
    diffusivities = np.array([1.1e-6, 1.2e-6, 1.3e-6])
    expected = np.array([454545.4545, 500000.0, 538461.5384])
    assert np.allclose(
        gao_et_al.volumetric_heat_capacity(lambda_vals, diffusivities),
        expected,
        atol=1e-2,
    )

    # Test with non-positive k_s_val
    with pytest.raises(ValueError):
        gao_et_al.volumetric_heat_capacity(0.8, 0)
    with pytest.raises(ValueError):
        gao_et_al.volumetric_heat_capacity(0.8, -1.2e-6)


def test_nme():
    """Test the nme function."""
    # Test with single values
    assert np.isclose(gao_et_al.nme(4.5, 5.0), 10.0)

    # Test with vectors
    calc = np.array([1.0, 2.1, 3.2])
    meas = np.array([1.2, 2.0, 3.0])
    assert np.isclose(gao_et_al.nme(calc, meas), 8.0645, atol=1e-4)

    # Test with broadcasting
    assert np.isclose(gao_et_al.nme(2.0, np.array([1.5, 2.5, 2.0])), 16.6666, atol=1e-4)

    # Test with zero denominator
    with pytest.raises(ValueError):
        gao_et_al.nme(np.array([1, 2, 3]), np.array([0, 0, 0]))

    # Test with incompatible shapes
    with pytest.raises(ValueError):
        gao_et_al.nme(np.array([1, 2, 3]), np.array([1, 2]))


def test_rmse():
    """Test the rmse function."""
    # Test with scalar inputs
    assert np.isclose(gao_et_al.rmse(4.5, 5.0), 0.5)

    # Test with vector inputs
    calc = np.array([2.1, 3.0, 4.2])
    meas = np.array([2.0, 3.5, 4.0])
    assert np.isclose(gao_et_al.rmse(calc, meas), 0.3162, atol=1e-4)

    # Test with broadcasting
    assert np.isclose(gao_et_al.rmse(3.0, np.array([2.5, 3.5, 3.0])), 0.4082, atol=1e-4)

    # Test with empty inputs
    with pytest.raises(ValueError):
        gao_et_al.rmse(np.array([]), np.array([]))


def test_calorimetric_gz():
    """Test the calorimetric_gz function."""
    # Three-layer example, single time step
    g_plate = -15.2
    Cv = np.array([2.5e6, 2.3e6, 2.1e6])
    dTdt = np.array([1.2e-4, 0.9e-4, 0.6e-4])
    dz = np.array([0.02, 0.02, 0.01])
    expected = np.array([-0.2, -4.85, -8.9])
    assert np.allclose(gao_et_al.calorimetric_gz(g_plate, Cv, dTdt, dz), expected)

    # Vectorized daily time series with two layers
    g_plate = np.random.normal(-10, 2, 1440)
    Cv = np.array([[2.4e6], [2.2e6]])
    dTdt = np.random.normal(5e-5, 2e-5, (2, 1440))
    dz = np.array([0.03, 0.02])
    Gz = gao_et_al.calorimetric_gz(g_plate, Cv, dTdt, dz)
    assert Gz.shape == (1440,)

    # Test with incompatible shapes by making dz have a different number of layers
    with pytest.raises(ValueError):
        gao_et_al.calorimetric_gz(g_plate, Cv, dTdt, np.array([0.1, 0.2, 0.3]))


def test_force_restore_gz():
    cv = 2e6
    dTg_dt = 1e-4
    Tg = 300
    Tg_bar = 295
    g = gao_et_al.force_restore_gz(cv, dTg_dt, Tg, Tg_bar)
    assert np.isclose(g, 83.305, atol=1e-3)


def test_lambda_s_from_cv():
    assert gao_et_al.lambda_s_from_cv(2.2e6) == 1.0
    with pytest.raises(ValueError):
        gao_et_al.lambda_s_from_cv(0)


def test_gao2010_gz():
    g = gao_et_al.gao2010_gz(8.0, 1.0, 1e-6, 3600)
    assert np.isclose(g, 59.081, atol=1e-3)


def test_heusinkveld_gz():
    g = gao_et_al.heusinkveld_gz([10], [0], 1, 1e-6, 1.0, 2 * np.pi / 86400)
    assert np.isclose(g[0], 26395728.157, atol=1e-3)


def test_hsieh2009_gz():
    t = [0, 3600]
    T = [290, 291]
    cv = [2e6, 2e6]
    ks = [1e-6, 1e-6]
    g = gao_et_al.hsieh2009_gz(T, t, cv, ks)
    assert np.isclose(g, 0.0265, atol=1e-3)


def test_leuning_damping_depth():
    d = gao_et_al.leuning_damping_depth(0.1, 0.05, 4, 8)
    assert np.isclose(d, 0.07213, atol=1e-5)


def test_leuning_gz():
    g = gao_et_al.leuning_gz(10, 0.05, 0.1, 0.1)
    assert np.isclose(g, 16.487, atol=1e-3)


def test_simple_measurement_gz():
    g_zr = [10, 12, 11]
    T = [[15, 16, 15.5], [14, 14.2, 14.1]]
    cv = [2e6, 2e6]
    dz = [0.03, 0.02]
    g = gao_et_al.simple_measurement_gz(g_zr, cv, T, 3600, dz)
    assert len(g) == 1
    assert np.isclose(g[0], 7.277, atol=1e-3)


def test_wbz12_g_gz():
    g_zr = [10, 11, 12]
    t = [0, 600, 1200]
    g = gao_et_al.wbz12_g_gz(g_zr, t, 0.05, 0.08, 1e-6)
    assert len(g) == 3
    assert np.isclose(g[2], -5.357, atol=1e-3)


def test_wbz12_s_gz():
    g = gao_et_al.wbz12_s_gz(8.0, 1e-6, 0.08, 0.05, 3600, 0.1)
    assert np.isclose(g, 1.193, atol=1e-3)


def test_exact_temperature_gz():
    T = gao_et_al.exact_temperature_gz(0.05, 8, 3600, 0.1)
    assert np.isclose(T, 297.005, atol=1e-3)


def test_exact_gz():
    g = gao_et_al.exact_gz(0.05, 8, 1.0, 0.1, 3600)
    assert np.isclose(g, 35.703, atol=1e-3)
