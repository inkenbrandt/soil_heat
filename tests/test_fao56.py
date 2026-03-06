import numpy as np
import pytest
from soil_heat import fao56_soil_heat_flux

def test_soil_heat_flux_general():
    # FAO-56 Example 13 — March to April
    g = fao56_soil_heat_flux.soil_heat_flux_general(16.1, 14.1, delta_t=30.0, delta_z=1.0)
    assert np.isclose(g, 0.14)

    # Array input
    t_curr = np.array([16.1, 18.8])
    t_prev = np.array([14.1, 16.1])
    g_arr = fao56_soil_heat_flux.soil_heat_flux_general(t_curr, t_prev, delta_t=30.0)
    assert np.allclose(g_arr, [0.14, 0.189], atol=1e-3)

def test_soil_heat_flux_daily():
    assert fao56_soil_heat_flux.soil_heat_flux_daily() == 0.0

def test_soil_heat_flux_monthly():
    # FAO-56 Example 13 — April soil heat flux (March = 14.1, May = 18.8)
    g = fao56_soil_heat_flux.soil_heat_flux_monthly(14.1, 18.8)
    assert np.isclose(g, 0.329)

    # Array input
    t_prev = np.array([14.1, 16.1])
    t_next = np.array([18.8, 20.0])
    g_arr = fao56_soil_heat_flux.soil_heat_flux_monthly(t_prev, t_next)
    assert np.allclose(g_arr, [0.329, 0.273])

def test_soil_heat_flux_monthly_prev_only():
    # FAO-56 Example 13 — March soil heat flux (Feb = 12.1, Mar = 14.1)
    g = fao56_soil_heat_flux.soil_heat_flux_monthly_prev_only(12.1, 14.1)
    assert np.isclose(g, 0.28)

def test_soil_heat_flux_hourly():
    # Daytime
    assert np.isclose(fao56_soil_heat_flux.soil_heat_flux_hourly(2.5, True), 0.25)
    # Nighttime
    assert np.isclose(fao56_soil_heat_flux.soil_heat_flux_hourly(-0.4, False), -0.2)
    # Array
    rn = np.array([2.5, -0.4])
    day = np.array([True, False])
    assert np.allclose(fao56_soil_heat_flux.soil_heat_flux_hourly(rn, day), [0.25, -0.2])

def test_soil_heat_flux_hourly_auto():
    assert np.isclose(fao56_soil_heat_flux.soil_heat_flux_hourly_auto(2.5), 0.25)
    assert np.isclose(fao56_soil_heat_flux.soil_heat_flux_hourly_auto(-0.4), -0.2)

def test_soil_heat_flux_monthly_series():
    temps = [14.1, 16.1, 18.8]
    # Centered
    g_centered = fao56_soil_heat_flux.soil_heat_flux_monthly_series(temps, method="centered")
    expected_centered = [0.14 * (16.1 - 14.1), 0.07 * (18.8 - 14.1), 0.14 * (18.8 - 16.1)]
    assert np.allclose(g_centered, expected_centered)

    # Backward
    g_backward = fao56_soil_heat_flux.soil_heat_flux_monthly_series(temps, method="backward")
    assert np.isnan(g_backward[0])
    assert np.isclose(g_backward[1], 0.14 * (16.1 - 14.1))

    with pytest.raises(ValueError):
        fao56_soil_heat_flux.soil_heat_flux_monthly_series([14.1], method="centered")
    with pytest.raises(ValueError):
        fao56_soil_heat_flux.soil_heat_flux_monthly_series(temps, method="invalid")

def test_soil_heat_flux_annual_cycle():
    temps = [5.2, 6.1, 9.3, 13.0, 17.5, 22.1, 25.4, 24.8, 20.6, 14.9, 9.2, 5.8]
    g = fao56_soil_heat_flux.soil_heat_flux_annual_cycle(temps)
    assert len(g) == 12
    # Jan: 0.07 * (Feb - Dec)
    assert np.isclose(g[0], 0.07 * (6.1 - 5.8))
    # Dec: 0.07 * (Jan - Nov)
    assert np.isclose(g[11], 0.07 * (5.2 - 9.2))

    with pytest.raises(ValueError):
        fao56_soil_heat_flux.soil_heat_flux_annual_cycle(temps[:11])

def test_conversions():
    assert np.isclose(fao56_soil_heat_flux.energy_to_evap(0.329), 0.329 * 0.408)
    assert np.isclose(fao56_soil_heat_flux.evap_to_energy(0.134232), 0.134232 / 0.408)
    assert np.isclose(fao56_soil_heat_flux.mj_m2_day_to_w_m2(0.329), 0.329 / 0.0864)
    assert np.isclose(fao56_soil_heat_flux.w_m2_to_mj_m2_day(3.80787), 3.80787 * 0.0864, atol=1e-5)
