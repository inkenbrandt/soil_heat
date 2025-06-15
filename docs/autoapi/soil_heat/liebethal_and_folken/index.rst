soil_heat.liebethal_and_folken
==============================

.. py:module:: soil_heat.liebethal_and_folken

.. autoapi-nested-parse::

   liebethal_and_folken.py
   ==================================
   A collection of Python functions that implement every numbered
   equation from Liebethal & Foken (2006) *Evaluation of six
   parameterization approaches for the ground heat flux*.

   Each public function is named after the paper section and
   equation number for easy cross‑referencing.  Helper utilities
   for finite‐difference gradients and unit handling are provided
   at the end of the module.

   .. rubric:: References

   Liebethal, C., & Foken, T. (2006). Evaluation of six parameterization
   approaches for the ground heat flux. *Theoretical and Applied Climatology*.
   DOI:10.1007/s00704‑005‑0234‑0



Functions
---------

.. autoapisummary::

   soil_heat.liebethal_and_folken.reference_ground_heat_flux
   soil_heat.liebethal_and_folken.ground_heat_flux_pr
   soil_heat.liebethal_and_folken.ground_heat_flux_lr
   soil_heat.liebethal_and_folken.ur_coefficients
   soil_heat.liebethal_and_folken.ground_heat_flux_ur
   soil_heat.liebethal_and_folken.surface_temp_amplitude
   soil_heat.liebethal_and_folken.phi_from_soil_moisture
   soil_heat.liebethal_and_folken.ground_heat_flux_sh
   soil_heat.liebethal_and_folken.ground_heat_flux_sm
   soil_heat.liebethal_and_folken.active_layer_thickness
   soil_heat.liebethal_and_folken.ground_heat_flux_fr


Module Contents
---------------

.. py:function:: reference_ground_heat_flux(temp_profile: numpy.ndarray, depths: Sequence[float], times: Sequence[float], cv: float, thermal_conductivity: float, gradient_depth: float = 0.2) -> numpy.ndarray

   Compute the reference ground‑heat flux *G₀,M* (Eq. 1).

   Equation
   --------
   G₀,M(t) = -λ ∂T/∂z |_(z=0.2 m) + ∫_{z=0}^{0.2 m} c_v ∂T/∂t dz

   :param temp_profile: Soil temperatures (°C or K) at the depths specified by *depths* and
                        time stamps *times*.
   :type temp_profile: ndarray, shape (n_z, n_t)
   :param depths: Measurement depths (m, **positive downward**).
   :type depths: sequence of float, length *n_z*
   :param times: Epoch time in **seconds** (may be monotonic pandas DatetimeIndex
                 converted via ``astype('int64')/1e9``).
   :type times: sequence of float, length *n_t*
   :param cv: Volumetric heat capacity of the soil (J m⁻³ K⁻¹).
   :type cv: float
   :param thermal_conductivity: Soil thermal conductivity λ (W m⁻¹ K⁻¹).
   :type thermal_conductivity: float
   :param gradient_depth: Depth (m) at which the vertical gradient term is evaluated.
   :type gradient_depth: float, default 0.20

   :returns: Instantaneous ground‑heat flux *G₀,M* (W m⁻²). Positive = downward.
   :rtype: ndarray, shape (n_t,)


.. py:function:: ground_heat_flux_pr(qs: numpy.ndarray, p: float) -> numpy.ndarray

   Ground heat flux using a fixed *p* fraction of net radiation (Eq. 2).

   G₀,PR(t) = ‑p · Q*ₛ(t)

   :param qs: Net radiation time series (W m⁻²). Positive = downward.
   :type qs: ndarray
   :param p: Fraction of net radiation that becomes ground‑heat flux (0–1).
   :type p: float

   :returns: G₀,PR (W m⁻²).
   :rtype: ndarray


.. py:function:: ground_heat_flux_lr(qs: numpy.ndarray, a: float, b: float, lag_steps: int = 0) -> numpy.ndarray

   Linear net‑radiation parameterisation (Eq. 3).

   G₀,LR(t) = a·Q*ₛ(t+Δt_G) + b

   :param qs: Net radiation (W m⁻²).
   :type qs: ndarray
   :param a: Regression coefficients.
   :type a: float
   :param b: Regression coefficients.
   :type b: float
   :param lag_steps: Integer lag (number of samples) by which *qs* is advanced.
   :type lag_steps: int, default 0

   :returns: G₀,LR (W m⁻²).
   :rtype: ndarray


.. py:function:: ur_coefficients(delta_ts: float | numpy.ndarray) -> Tuple[numpy.ndarray, numpy.ndarray]

   Compute universal‑function parameters *A* and *B* (Eq. 5 & 6).

   :param delta_ts: Diurnal amplitude of *surface* temperature (K).
   :type delta_ts: float or ndarray

   :returns: * **A** (*ndarray*)
             * **B** (*ndarray (seconds)*)


.. py:function:: ground_heat_flux_ur(qs: numpy.ndarray, times_sec: numpy.ndarray, delta_ts: float) -> numpy.ndarray

   Universal net‑radiation parameterisation (Eq. 4).

   Implements Santanello & Friedl (2003):
       G₀,UR(t) = -A · cos[2π (t + 10800) / B] · Q*ₛ(t)

   *t* is **seconds since solar noon** (positive in afternoon).

   :param qs: Net radiation (W m⁻²).
   :type qs: ndarray
   :param times_sec: Seconds relative to solar noon (s).
   :type times_sec: ndarray
   :param delta_ts: Diurnal surface‑temperature amplitude (K).
   :type delta_ts: float

   :returns: G₀,UR (W m⁻²).
   :rtype: ndarray


.. py:function:: surface_temp_amplitude(delta_t1: float, delta_t2: float, z1: float, z2: float) -> float

   Compute diurnal surface‑temperature amplitude ΔT_s (Eq. 8).

   :param delta_t1: Diurnal temperature amplitudes (K) measured at depths *z1* and *z2*.
   :type delta_t1: float
   :param delta_t2: Diurnal temperature amplitudes (K) measured at depths *z1* and *z2*.
   :type delta_t2: float
   :param z1: Depths in meters (**positive downward**, with z2 > z1 > 0).
   :type z1: float
   :param z2: Depths in meters (**positive downward**, with z2 > z1 > 0).
   :type z2: float

   :returns: Estimated ΔT_s (K).
   :rtype: float


.. py:function:: phi_from_soil_moisture(theta_0_10: float, a_phi: float = 9.62, b_phi: float = 0.402) -> float

   Soil‑moisture dependent φ (Eq. 10).


.. py:function:: ground_heat_flux_sh(h: numpy.ndarray, phase_g0: Sequence[float], phase_h: Sequence[float], u_mean: float, phi: float, omega: float = 2 * np.pi / 86400.0) -> numpy.ndarray

   Ground‑heat flux from sensible heat flux H (Eq. 9).

   :param h: Sensible heat flux time series (W m⁻²).
   :type h: ndarray
   :param phase_g0: Phase lags φ(G₀) and φ(H) in **radians**.
   :type phase_g0: sequence of float
   :param phase_h: Phase lags φ(G₀) and φ(H) in **radians**.
   :type phase_h: sequence of float
   :param u_mean: Mean horizontal wind speed during daytime (m s⁻¹).
   :type u_mean: float
   :param phi: Empirical parameter (dimensionless), see `phi_from_soil_moisture`.
   :type phi: float
   :param omega: Diurnal angular frequency (s⁻¹).
   :type omega: float, default 2π/86400

   :returns: G₀,SH (W m⁻²).
   :rtype: ndarray


.. py:function:: ground_heat_flux_sm(gp: numpy.ndarray, t1: numpy.ndarray, delta_t: numpy.ndarray, cv: float, zp: float, dt_seconds: float) -> numpy.ndarray

   Simple‑measurement parameterisation (Eq. 11).

   :param gp: Heat‑flux plate measurement at depth *zp* (W m⁻²).
   :type gp: ndarray
   :param t1: Soil temperature at 0.01 m depth (K or °C).
   :type t1: ndarray
   :param delta_t: Temperature difference T(0.01 m) – T(z_p) (K).
   :type delta_t: ndarray
   :param cv: Volumetric heat capacity (J m⁻³ K⁻¹).
   :type cv: float
   :param zp: Plate depth (m, positive downward).
   :type zp: float
   :param dt_seconds: Time step between consecutive samples (s).
   :type dt_seconds: float

   :returns: G₀,SM (W m⁻²).
   :rtype: ndarray


.. py:function:: active_layer_thickness(lambda_: float, cv: float, omega: float = 2 * np.pi / 86400) -> float

   Thickness δz of the active soil layer (Eq. 13).


.. py:function:: ground_heat_flux_fr(tg: numpy.ndarray, tg_avg: float, cv: float, lambda_: float, delta_z: float | None = None, times: numpy.ndarray | None = None) -> numpy.ndarray

   Force‑restore ground‑heat flux (Eq. 12).

   Implements the two‑layer force‑restore formulation with an optional
   diagnostic *δz* computed via Eq. 13 if not supplied.

   :param tg: Temperature of the upper (surface) layer Tg(t).
   :type tg: ndarray
   :param tg_avg: Long‑term average or restoring temperature Tḡ (K).
   :type tg_avg: float
   :param cv: Volumetric heat capacity (J m⁻³ K⁻¹).
   :type cv: float
   :param lambda_: Soil thermal conductivity λ (W m⁻¹ K⁻¹).
   :type lambda_: float
   :param delta_z: Thickness of the active soil layer δz (m).  If *None*, computed
                   from ``active_layer_thickness``.
   :type delta_z: float, optional
   :param times: Time stamps in seconds.  Required when *delta_z* is None or when
                 irregular sampling; defaults to `np.arange(len(tg))` seconds.
   :type times: ndarray, optional

   :returns: G₀,FR (W m⁻²).
   :rtype: ndarray


