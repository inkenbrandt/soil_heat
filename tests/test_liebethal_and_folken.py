import numpy as np
import pytest
from soil_heat import liebethal_and_folken

def test_central_gradient():
    y = np.array([1, 2, 4])
    x = np.array([0, 1, 2])
    # grad at 0: (2-1)/(1-0) = 1.0
    # grad at 1: (4-1)/(2-0) = 1.5
    # grad at 2: (4-2)/(2-1) = 2.0
    grad = liebethal_and_folken._central_gradient(y, x)
    assert np.allclose(grad, [1.0, 1.5, 2.0])

    with pytest.raises(ValueError):
        liebethal_and_folken._central_gradient(y, np.array([0, 1]))

def test_pad_nan_like():
    arr = np.array([1, 2])
    nan_arr = liebethal_and_folken._pad_nan_like(arr)
    assert nan_arr.shape == arr.shape
    assert np.all(np.isnan(nan_arr))

def test_reference_ground_heat_flux():
    depths = [0.0, 0.1, 0.2]
    times = [0, 3600]
    temp_profile = np.array([[10, 11], [12, 12.5], [14, 14]])
    cv = 2.4e6
    k = 1.0
    g0 = liebethal_and_folken.reference_ground_heat_flux(temp_profile, depths, times, cv, k)
    assert g0.shape == (2,)

def test_ground_heat_flux_pr():
    qs = np.array([100, 200])
    p = 0.1
    g0 = liebethal_and_folken.ground_heat_flux_pr(qs, p)
    assert np.allclose(g0, [-10, -20])

def test_ground_heat_flux_lr():
    qs = np.array([100, 200, 150])
    a = 0.5
    b = 5.0
    g0 = liebethal_and_folken.ground_heat_flux_lr(qs, a, b, lag_steps=1)
    # lag 1: qs becomes [200, 150, 100]
    # g0: [200*0.5+5, 150*0.5+5, 100*0.5+5] = [105, 80, 55]
    assert np.allclose(g0, [105, 80, 55])

def test_ur_coefficients():
    A, B = liebethal_and_folken.ur_coefficients(10.0)
    assert np.isclose(A, 0.0074 * 10.0 + 0.088)
    assert np.isclose(B, 1729.0 * 10.0 + 65013.0)

def test_ground_heat_flux_ur():
    qs = np.array([100, 200])
    times = np.array([0, 3600])
    g0 = liebethal_and_folken.ground_heat_flux_ur(qs, times, 10.0)
    assert g0.shape == (2,)

def test_surface_temp_amplitude():
    amp = liebethal_and_folken.surface_temp_amplitude(5.0, 2.0, 0.05, 0.10)
    assert amp > 5.0
    with pytest.raises(ValueError):
        liebethal_and_folken.surface_temp_amplitude(5.0, 2.0, 0.10, 0.05)

def test_phi_from_soil_moisture():
    phi = liebethal_and_folken.phi_from_soil_moisture(0.25)
    assert phi == 9.62 * 0.25 + 0.402

def test_ground_heat_flux_sh():
    h = np.array([50, 60])
    phase_g0 = [0, 0.1]
    phase_h = [0.1, 0.2]
    u_mean = 2.0
    phi = 2.0
    g0 = liebethal_and_folken.ground_heat_flux_sh(h, phase_g0, phase_h, u_mean, phi)
    assert g0.shape == (2,)

    with pytest.raises(ValueError):
        liebethal_and_folken.ground_heat_flux_sh(h, [0], [0], u_mean, phi)

def test_ground_heat_flux_sm():
    gp = np.array([10, 12])
    t1 = np.array([15, 15.5])
    delta_t = np.array([1, 1.1])
    cv = 2.4e6
    zp = 0.08
    dt = 3600
    g0 = liebethal_and_folken.ground_heat_flux_sm(gp, t1, delta_t, cv, zp, dt)
    assert np.isnan(g0[0])
    assert not np.isnan(g0[1])

def test_active_layer_thickness():
    dz = liebethal_and_folken.active_layer_thickness(1.0, 2.4e6)
    assert dz > 0

def test_ground_heat_flux_fr():
    tg = np.array([15, 16, 15.5])
    tg_avg = 15.2
    cv = 2.4e6
    k = 1.0
    g0 = liebethal_and_folken.ground_heat_flux_fr(tg, tg_avg, cv, k)
    assert g0.shape == (3,)
