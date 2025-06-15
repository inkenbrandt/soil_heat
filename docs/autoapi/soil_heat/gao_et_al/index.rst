soil_heat.gao_et_al
===================

.. py:module:: soil_heat.gao_et_al


Functions
---------

.. autoapisummary::

   soil_heat.gao_et_al.lambda_s
   soil_heat.gao_et_al.k_s
   soil_heat.gao_et_al.volumetric_heat_capacity
   soil_heat.gao_et_al.nme
   soil_heat.gao_et_al.rmse
   soil_heat.gao_et_al.calorimetric_gz
   soil_heat.gao_et_al.force_restore_gz
   soil_heat.gao_et_al.gao2010_gz
   soil_heat.gao_et_al.heusinkveld_gz
   soil_heat.gao_et_al.hsieh2009_gz
   soil_heat.gao_et_al.leuning_damping_depth
   soil_heat.gao_et_al.leuning_gz
   soil_heat.gao_et_al.simple_measurement_gz
   soil_heat.gao_et_al.wbz12_g_gz
   soil_heat.gao_et_al.wbz12_s_gz
   soil_heat.gao_et_al.exact_temperature_gz
   soil_heat.gao_et_al.exact_gz


Module Contents
---------------

.. py:function:: lambda_s(theta: numpy.ndarray | float) -> numpy.ndarray | float

   Thermal conductivity (λ_s) as a function of volumetric water content θ.

   Implements Eq. (12) from Gao et al. (2017) fileciteturn1file3.


.. py:function:: k_s(theta: numpy.ndarray | float) -> numpy.ndarray | float

   Thermal diffusivity (k_s) as a function of volumetric water content θ.

   Implements Eq. (13) from Gao et al. (2017) fileciteturn1file6.


.. py:function:: volumetric_heat_capacity(lambda_s_val, k_s_val)

   Volumetric heat capacity C_v = λ_s / k_s (units J m⁻³ K⁻¹).


.. py:function:: nme(calc: numpy.ndarray, meas: numpy.ndarray) -> float

   Normalized mean error (%). Implements Eq. (14) fileciteturn1file8.


.. py:function:: rmse(calc: numpy.ndarray, meas: numpy.ndarray) -> float

   Root‑mean‑square error. Implements Eq. (15) fileciteturn1file8.


.. py:function:: calorimetric_gz(g_zr, cv_layers, dT_dt_layers, dz_layers)

   Calorimetric method for G_z at depth z (usually 5 cm).

   :param g_zr: Measured heat flux at reference depth *z_r* (W m⁻²).
   :type g_zr: float or array_like
   :param cv_layers: Volumetric heat capacity for each sub‑layer *C_v,l* (J m⁻³ K⁻¹).
   :type cv_layers: sequence
   :param dT_dt_layers: Time derivative of average temperature for each layer ∂T/∂t (K s⁻¹).
   :type dT_dt_layers: sequence
   :param dz_layers: Thickness of each sub‑layer δz_l (m).
   :type dz_layers: sequence


.. py:function:: force_restore_gz(cv, dTg_dt, Tg, Tg_bar, delta_z=0.05, omega=OMEGA_DAY)

   Force‑restore estimate of G_z at z = δz (default 5 cm).

   Implements Eq. (2) fileciteturn1file2.


.. py:function:: gao2010_gz(AT, lambda_s_val, k_s_val, t, omega=OMEGA_DAY)

   Sinusoidal solution for G_z at depth d (Eq. 3).


.. py:function:: heusinkveld_gz(A_n, Phi_n, n_max, k_s_val, lambda_s_val, w)

   H04 harmonic solution (Eq. 4).


.. py:function:: hsieh2009_gz(tz_series, time_series, cv_series, ks_series)

   Half‑order integral solution (Eq. 5).

   tz_series, cv_series, ks_series must be monotonically increasing in *time_series*.


.. py:function:: leuning_damping_depth(z, zr, AT_z, AT_zr)

   Compute damping depth *d* via Eq. (6).


.. py:function:: leuning_gz(g_zr, z, zr, d)

   Exponentially adjust G_z from reference depth (Eq. 7).


.. py:function:: simple_measurement_gz(g_zr, cv_layers, tz_layers, dt, dz_layers)

   Simple‑measurement variant of calorimetric method (Eq. 8).


.. py:function:: wbz12_g_gz(g_zr_series, time_series, z, zr, k_s_val)

   WBZ12‑G method (Eq. 9–10).


.. py:function:: wbz12_s_gz(Ag, ks_val, zr, z, t, eps, omega=OMEGA_DAY)

   WBZ12‑S solution (Eq. 11).

   Numerical integration is performed for the second term.


.. py:function:: exact_temperature_gz(z, AT, t, d, omega=OMEGA_DAY, T_i=298.15)

   Exact sinusoidal soil‑temperature profile (Eq. 16).


.. py:function:: exact_gz(z, AT, lambda_s_val, d, t, omega=OMEGA_DAY)

   Exact sinusoidal soil‑heat‑flux (Eq. 17).


