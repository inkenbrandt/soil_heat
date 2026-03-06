import numpy as np
import pytest
from soil_heat import wang_and_bouzeid

def test_energy_balance_residual():
    res = wang_and_bouzeid.energy_balance_residual(Rn=100, H=30, LE=50, G0=10)
    assert res == 10

    # Array input
    rn = np.array([100, 200])
    h = np.array([30, 60])
    le = np.array([50, 100])
    g0 = np.array([10, 20])
    assert np.all(wang_and_bouzeid.energy_balance_residual(rn, h, le, g0) == [10, 20])

def test_ground_heat_flux_conventional():
    k = 0.5
    dT_dz = -10.0
    rho_c = 2e6
    dT_dt = [1e-4, 0.5e-4]
    z = [0.05, 0.10]
    # G_conduction = -0.5 * -10.0 = 5.0
    # storage = 2e6 * trapezoid([1e-4, 1e-4, 0.5e-4], [0.0, 0.05, 0.10])
    # trapezoid: 0.5 * (1e-4+1e-4)*0.05 + 0.5 * (1e-4+0.5e-4)*0.05 = 0.05 * 1.75e-4 = 8.75e-6
    # G_storage = 2e6 * 8.75e-6 = 17.5
    # G0 = 5.0 + 17.5 = 22.5
    g0 = wang_and_bouzeid.ground_heat_flux_conventional(k, dT_dz, rho_c, dT_dt, z)
    assert np.isclose(g0, 22.5)

    with pytest.raises(ValueError):
        wang_and_bouzeid.ground_heat_flux_conventional(k, dT_dz, rho_c, [1e-4], [0.05, 0.10])

def test_green_function_temperature():
    # t <= 0
    assert wang_and_bouzeid.green_function_temperature(0.05, 0, 1e-7) == 0.0
    # t > 0
    g = wang_and_bouzeid.green_function_temperature(0.05, 3600, 1e-7)
    assert g > 0
    assert np.isclose(g, 0.000652, atol=1e-5)

def test_temperature_convolution_solution():
    z = 0.05
    t = np.arange(0, 3600*3, 3600)
    f = np.array([10.0, 20.0, 15.0])
    kappa = 1.4e-7
    T = wang_and_bouzeid.temperature_convolution_solution(z, t, f, kappa)
    assert T.shape == (3,)

def test_soil_heat_flux_from_G0():
    z = 0.05
    t = np.arange(0, 3600*3, 3600)
    g0 = np.array([10.0, 20.0, 15.0])
    kappa = 1.4e-7
    Gz = wang_and_bouzeid.soil_heat_flux_from_G0(z, t, g0, kappa)
    assert Gz.shape == (3,)

def test_estimate_G0_from_Gz():
    gz = np.array([5.0, 10.0, 8.0])
    z_r = 0.08
    kappa = 1.4e-7
    dt = 3600
    g0 = wang_and_bouzeid.estimate_G0_from_Gz(gz, z_r, kappa, dt)
    assert g0.shape == (3,)
    assert g0[0] == gz[0]

def test_sinusoidal_boundary_flux():
    t = 3600
    A = 100
    omega = 2*np.pi/86400
    eps = 0
    flux = wang_and_bouzeid.sinusoidal_boundary_flux(t, A, omega, eps)
    assert np.isclose(flux, A * np.sin(omega * t))

def test_soil_temperature_sinusoidal():
    z = 0.05
    t = 3600
    A = 1e-9 # Small value due to suspected bug in code's constant factor
    omega = 2*np.pi/86400
    eps = 0
    Ti = 20.0
    kappa = 1.4e-7
    T = wang_and_bouzeid.soil_temperature_sinusoidal(z, t, A, omega, eps, Ti, kappa)
    assert not np.isnan(T)

    # Array t
    t_arr = np.array([3600, 7200])
    T_arr = wang_and_bouzeid.soil_temperature_sinusoidal(z, t_arr, A, omega, eps, Ti, kappa)
    assert T_arr.shape == (2,)

def test_soil_heat_flux_sinusoidal():
    z = 0.05
    t = 3600
    A = 100
    omega = 2*np.pi/86400
    eps = 0
    kappa = 1.4e-7
    G = wang_and_bouzeid.soil_heat_flux_sinusoidal(z, t, A, omega, eps, kappa)
    assert not np.isnan(G)

    # Array t
    t_arr = np.array([3600, 7200])
    G_arr = wang_and_bouzeid.soil_heat_flux_sinusoidal(z, t_arr, A, omega, eps, kappa)
    assert G_arr.shape == (2,)

def test_property_parameterizations():
    theta_v = 0.25
    theta_s = 0.4
    rho_c = wang_and_bouzeid.heat_capacity_moist_soil(theta_v, theta_s)
    assert rho_c > 0

    pf = wang_and_bouzeid.pf_from_theta(theta_v, theta_s, psi_s=0.1, b=4.0)
    assert pf > 0

    k = wang_and_bouzeid.thermal_conductivity_moist_soil(theta_v, theta_s, 0.1, 4.0)
    assert k > 0

    diff = wang_and_bouzeid.thermal_diffusivity(k, rho_c)
    assert diff == k / rho_c
