import numpy as np
import pytest
from soil_heat import wang_and_yang

def test_soil_heat_flux():
    Tz = [10, 12, 14]
    dz = [0.05, 0.05, 0.05]
    lam = 1.0
    G = wang_and_yang.soil_heat_flux(Tz, dz, lam)
    # n_layers = 3, output length should be 4
    assert len(G) == 4
    assert np.allclose(G[1:3], [-1.0 * (12-10)/0.05, -1.0 * (14-12)/0.05])

def test_integrated_soil_heat_flux():
    rho_c = [2e6, 2e6]
    T_before = [10, 11]
    T_after = [10.5, 11.2]
    dz = [0.05, 0.05]
    dt = 3600
    G = wang_and_yang.integrated_soil_heat_flux(rho_c, T_before, T_after, dz, dt)
    assert len(G) == 2

def test_volumetric_heat_capacity():
    theta = 0.25
    theta_sat = 0.4
    cv = wang_and_yang.volumetric_heat_capacity(theta, theta_sat)
    assert cv == (1.0 - 0.4) * 2.1e6 + 4.2e6 * 0.25

def test_stretched_grid():
    dz = wang_and_yang.stretched_grid(5, 1.0, 0.1)
    assert len(dz) == 5
    assert np.isclose(np.sum(dz), 1.0)

    dz_unif = wang_and_yang.stretched_grid(5, 1.0, 0)
    assert np.allclose(dz_unif, 0.2)

def test_tridiagonal_coeffs():
    dz = [0.05, 0.05, 0.05]
    rho_c = [2e6, 2e6, 2e6]
    lam = 1.0
    dt = 3600
    A, B, C = wang_and_yang.tridiagonal_coeffs(dz, rho_c, lam, dt)
    assert len(A) == 2
    assert len(B) == 3
    assert len(C) == 2

def test_solve_tde():
    T_prev = [290, 290, 290]
    dz = [0.05, 0.05, 0.05]
    rho_c = [2e6, 2e6, 2e6]
    lam = 1.0
    Tsfc = 295
    Tbot = 285
    dt = 3600
    T_new = wang_and_yang.solve_tde(T_prev, dz, rho_c, lam, Tsfc, Tbot, dt)
    assert len(T_new) == 5

def test_correct_profile():
    T_model = np.array([10, 12, 14])
    depths_model = np.array([0.05, 0.10, 0.20])
    T_obs = np.array([11, 13.5])
    depths_obs = np.array([0.05, 0.15])
    T_corr = wang_and_yang.correct_profile(T_model, depths_model, T_obs, depths_obs)
    assert T_corr.shape == T_model.shape

def test_surface_temperature_longwave():
    T = wang_and_yang.surface_temperature_longwave(400, 300)
    assert T > 0

def test_thermal_conductivity_yang2008():
    k = wang_and_yang.thermal_conductivity_yang2008(0.2, 0.4, 1300)
    assert k > 0

def test_flux_error_linear():
    err = wang_and_yang.flux_error_linear(2e6, 0.1, 3600)
    assert err == 2e6 * 0.1 / 3600

def test_surface_energy_residual():
    res = wang_and_yang.surface_energy_residual(100, 30, 50, 10)
    assert res == 10

def test_tdec_step():
    T_prev = [290, 290, 290]
    dz = [0.05, 0.05, 0.05]
    theta = [0.2, 0.2, 0.2]
    theta_sat = 0.4
    rho_dry = 1300
    lam = 1.0
    Tsfc = 295
    Tbot = 285
    dt = 3600
    depths_model = [0, 0.05, 0.10, 0.15, 0.20] # Including boundaries
    T_obs = [294, 289]
    depths_obs = [0.05, 0.15]

    T_corr, G_prof = wang_and_yang.tdec_step(
        T_prev, dz, theta, theta_sat, rho_dry, lam, Tsfc, Tbot, dt,
        depths_model, T_obs, depths_obs
    )
    assert len(T_corr) == 5
    assert len(G_prof) == 3
