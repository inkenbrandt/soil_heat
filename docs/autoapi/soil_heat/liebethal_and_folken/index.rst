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
   for finite‑difference gradients and unit handling are provided
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

   Compute the reference ground‑heat flux *G₀,M* using the
   gradient method combined with calorimetry for heat storage
   (Liebethal & Foken 2006, Eq. 1).

   The method combines the conductive heat flux at a reference depth
   with the rate of change of heat stored in the soil layer above that
   depth.

   .. math::

       G_{0,M}(t) = -\lambda \frac{\partial T}{\partial z} \bigg|_{z=0.2m}
                    + \int_{z=0}^{0.2m} c_v \frac{\partial T}{\partial t} dz

   :param temp_profile: A 2D array of soil temperatures (°C or K) with shape
                        `(n_depths, n_times)`.
   :type temp_profile: numpy.ndarray
   :param depths: A sequence of measurement depths in meters (positive downward),
                  corresponding to the rows of `temp_profile`.
   :type depths: Sequence[float]
   :param times: A sequence of time stamps in seconds (e.g., Unix timestamps),
                 corresponding to the columns of `temp_profile`.
   :type times: Sequence[float]
   :param cv: Volumetric heat capacity of the soil (J m⁻³ K⁻¹). Assumed
              to be constant with depth.
   :type cv: float
   :param thermal_conductivity: Soil thermal conductivity λ (W m⁻¹ K⁻¹). Assumed to be
                                constant with depth.
   :type thermal_conductivity: float
   :param gradient_depth: The depth (m) at which the vertical temperature gradient is
                          evaluated, by default 0.20 m.
   :type gradient_depth: float, optional

   :returns: A 1D array of the instantaneous ground‑heat flux *G₀,M* (W m⁻²)
             at the surface for each time step. Positive values indicate
             downward flux.
   :rtype: numpy.ndarray

   :raises ValueError: If `temp_profile` shape does not match the lengths of `depths`
       and `times`.


.. py:function:: ground_heat_flux_pr(qs: numpy.ndarray, p: float) -> numpy.ndarray

   Estimate ground heat flux as a fixed fraction of net radiation
   (Liebethal & Foken 2006, Eq. 2).

   This is a simple empirical relationship where the ground heat flux
   is assumed to be a constant proportion of the net radiation at the
   surface.

   .. math:: G_{0,PR}(t) = -p \cdot Q^*_s(t)

   :param qs: Time series of net radiation (W m⁻²). Positive values are
              typically downward.
   :type qs: numpy.ndarray
   :param p: The fraction of net radiation that is partitioned into
             ground heat flux (dimensionless, typically 0 < p < 1).
   :type p: float

   :returns: Time series of the estimated ground heat flux *G₀,PR* (W m⁻²).
   :rtype: numpy.ndarray


.. py:function:: ground_heat_flux_lr(qs: numpy.ndarray, a: float, b: float, lag_steps: int = 0) -> numpy.ndarray

   Estimate ground heat flux using a linear regression against net
   radiation, with an optional time lag (Liebethal & Foken 2006, Eq. 3).

   .. math:: G_{0,LR}(t) = a \cdot Q^*_s(t + \Delta t_G) + b

   :param qs: Time series of net radiation (W m⁻²).
   :type qs: numpy.ndarray
   :param a: The slope of the linear regression (dimensionless).
   :type a: float
   :param b: The intercept of the linear regression (W m⁻²).
   :type b: float
   :param lag_steps: The integer time lag (number of array elements) to apply to the
                     net radiation series. A positive value advances the series
                     (e.g., `qs[t+lag]` is used for `G[t]`), by default 0.
   :type lag_steps: int, optional

   :returns: Time series of the estimated ground heat flux *G₀,LR* (W m⁻²).
   :rtype: numpy.ndarray


.. py:function:: ur_coefficients(delta_ts: float | numpy.ndarray) -> Tuple[numpy.ndarray, numpy.ndarray]

   Compute the parameters A and B for the universal net-radiation
   parameterization, based on the diurnal surface temperature amplitude
   (Liebethal & Foken 2006, Eq. 5 & 6).

   .. math::
       A = 0.0074 \cdot \Delta T_s + 0.088
       B = 1729 \cdot \Delta T_s + 65013

   :param delta_ts: The diurnal amplitude of the surface temperature (K or °C).
   :type delta_ts: float or numpy.ndarray

   :returns: A tuple containing the computed parameters (A, B).
             - A is dimensionless.
             - B is in seconds.
   :rtype: Tuple[numpy.ndarray, numpy.ndarray]


.. py:function:: ground_heat_flux_ur(qs: numpy.ndarray, times_sec: numpy.ndarray, delta_ts: float) -> numpy.ndarray

   Estimate ground heat flux using the universal net-radiation
   parameterization by Santanello & Friedl (2003), as cited in
   Liebethal & Foken (2006, Eq. 4).

   This method modulates the fraction of net radiation partitioned to
   ground heat flux with a cosine function of the time of day.

   .. math::
       G_{0,UR}(t) = -A \cdot \cos\left(\frac{2\pi (t + 10800)}{B}\right)
                     \cdot Q^*_s(t)

   :param qs: Time series of net radiation (W m⁻²).
   :type qs: numpy.ndarray
   :param times_sec: Time stamps in **seconds relative to solar noon**. Positive values
                     indicate the afternoon.
   :type times_sec: numpy.ndarray
   :param delta_ts: The diurnal amplitude of the surface temperature (K or °C).
   :type delta_ts: float

   :returns: Time series of the estimated ground heat flux *G₀,UR* (W m⁻²).
   :rtype: numpy.ndarray


.. py:function:: surface_temp_amplitude(delta_t1: float, delta_t2: float, z1: float, z2: float) -> float

   Estimate the diurnal surface-temperature amplitude (ΔT_s) from
   temperature amplitudes measured at two different soil depths
   (Liebethal & Foken 2006, Eq. 8).

   This method assumes an exponential decay of the temperature wave
   amplitude with depth.

   .. math::
       \Delta T_s = \Delta T_1 + \Delta T_2 \cdot
                     \exp\left(\frac{z_2}{z_2 - z_1}\right)

   :param delta_t1: Diurnal temperature amplitude (K or °C) measured at depth `z1`.
   :type delta_t1: float
   :param delta_t2: Diurnal temperature amplitude (K or °C) measured at depth `z2`.
   :type delta_t2: float
   :param z1: The shallower depth in meters (positive downward).
   :type z1: float
   :param z2: The deeper depth in meters (positive downward).
   :type z2: float

   :returns: The estimated diurnal surface-temperature amplitude ΔT_s (K or °C).
   :rtype: float

   :raises ValueError: If `z2` is not greater than `z1`.


.. py:function:: phi_from_soil_moisture(theta_0_10: float, a_phi: float = 9.62, b_phi: float = 0.402) -> float

   Calculate the empirical parameter φ based on soil moisture content
   (Liebethal & Foken 2006, Eq. 10).

   .. math:: \phi = a_\phi \cdot \theta_{0-10} + b_\phi

   :param theta_0_10: Average volumetric soil moisture content in the top 10 cm
                      (m³ m⁻³).
   :type theta_0_10: float
   :param a_phi: Empirical coefficient, by default 9.62.
   :type a_phi: float, optional
   :param b_phi: Empirical coefficient, by default 0.402.
   :type b_phi: float, optional

   :returns: The dimensionless parameter φ.
   :rtype: float


.. py:function:: ground_heat_flux_sh(h: numpy.ndarray, phase_g0: Sequence[float], phase_h: Sequence[float], u_mean: float, phi: float, omega: float = 2 * np.pi / 86400.0) -> numpy.ndarray

   Estimate ground heat flux from the sensible heat flux (H)
   (Liebethal & Foken 2006, Eq. 9).

   This method relates the ground heat flux to the sensible heat flux
   through a phase-shifted and scaled relationship.

   .. math::
       G_{0,SH}(t) = -\frac{\phi}{\sqrt{\bar{u}}}
                     \frac{\cos(\omega t + \varphi(G_0))}
                          {\cos(\omega t + \varphi(H))} H(t)

   :param h: Time series of sensible heat flux (W m⁻²).
   :type h: numpy.ndarray
   :param phase_g0: Phase lags of the ground heat flux, φ(G₀), in **radians**. Must
                    have the same length as `h`.
   :type phase_g0: Sequence[float]
   :param phase_h: Phase lags of the sensible heat flux, φ(H), in **radians**. Must
                   have the same length as `h`.
   :type phase_h: Sequence[float]
   :param u_mean: Mean horizontal wind speed during the daytime (m s⁻¹).
   :type u_mean: float
   :param phi: An empirical dimensionless parameter, often derived from soil
               moisture via `phi_from_soil_moisture`.
   :type phi: float
   :param omega: The diurnal angular frequency (rad s⁻¹), by default `2π/86400`.
   :type omega: float, optional

   :returns: Time series of the estimated ground heat flux *G₀,SH* (W m⁻²).
   :rtype: numpy.ndarray

   :raises ValueError: If the length of phase arrays does not match the length of `h`.


.. py:function:: ground_heat_flux_sm(gp: numpy.ndarray, t1: numpy.ndarray, delta_t: numpy.ndarray, cv: float, zp: float, dt_seconds: float) -> numpy.ndarray

   Estimate surface ground heat flux using the "simple measurement"
   parameterization, correcting a heat flux plate measurement with a
   storage term (Liebethal & Foken 2006, Eq. 11).

   .. math::
       G_{0,SM}(t) = G_p(t) + c_v z_p \left( \frac{dT_1}{dt} +
                     \frac{1}{2} \frac{d(\Delta T)}{dt} \right)

   :param gp: Heat flux plate measurements at depth `zp` (W m⁻²).
   :type gp: numpy.ndarray
   :param t1: Soil temperature at 0.01 m depth (K or °C).
   :type t1: numpy.ndarray
   :param delta_t: Temperature difference T(0.01 m) - T(zp) (K).
   :type delta_t: numpy.ndarray
   :param cv: Volumetric heat capacity of the soil (J m⁻³ K⁻¹).
   :type cv: float
   :param zp: Depth of the heat flux plate (m, positive downward).
   :type zp: float
   :param dt_seconds: The constant time step between consecutive samples (s).
   :type dt_seconds: float

   :returns: Time series of the estimated ground heat flux *G₀,SM* (W m⁻²).
             The first element will be NaN due to the backward difference.
   :rtype: numpy.ndarray


.. py:function:: active_layer_thickness(lambda_: float, cv: float, omega: float = 2 * np.pi / 86400) -> float

   Calculate the thickness of the active soil layer (δz) for the
   force-restore method (Liebethal & Foken 2006, Eq. 13).

   .. math:: \delta_z = \sqrt{\frac{\lambda}{2 c_v \omega}}

   :param lambda_: Soil thermal conductivity (W m⁻¹ K⁻¹).
   :type lambda_: float
   :param cv: Volumetric heat capacity of the soil (J m⁻³ K⁻¹).
   :type cv: float
   :param omega: The diurnal angular frequency (rad s⁻¹), by default `2π/86400`.
   :type omega: float, optional

   :returns: The thickness of the active soil layer δz (m).
   :rtype: float


.. py:function:: ground_heat_flux_fr(tg: numpy.ndarray, tg_avg: float, cv: float, lambda_: float, delta_z: float | None = None, times: numpy.ndarray | None = None) -> numpy.ndarray

   Estimate ground heat flux using the two-layer force-restore method
   (Liebethal & Foken 2006, Eq. 12).

   .. math::
       G_{0,FR}(t) = -\delta_z c_v \frac{dT_g}{dt} +
                     \sqrt{\lambda \omega c_v} \left(
                     \frac{1}{\omega}\frac{dT_g}{dt} + (T_g - \bar{T_g})
                     \right)

   :param tg: Time series of the upper (surface) layer temperature, Tg(t) (K).
   :type tg: numpy.ndarray
   :param tg_avg: The long-term average or "restoring" temperature, T̄g (K).
   :type tg_avg: float
   :param cv: Volumetric heat capacity of the soil (J m⁻³ K⁻¹).
   :type cv: float
   :param lambda_: Soil thermal conductivity λ (W m⁻¹ K⁻¹).
   :type lambda_: float
   :param delta_z: Thickness of the active soil layer δz (m). If `None`, it is
                   calculated internally using `active_layer_thickness`,
                   by default None.
   :type delta_z: float, optional
   :param times: Time stamps in seconds corresponding to `tg`. Required for
                 calculating the time derivative. If `None`, assumes a uniform
                 time step of 1 second, by default None.
   :type times: numpy.ndarray, optional

   :returns: Time series of the estimated ground heat flux *G₀,FR* (W m⁻²).
   :rtype: numpy.ndarray


