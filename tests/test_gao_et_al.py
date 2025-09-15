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
