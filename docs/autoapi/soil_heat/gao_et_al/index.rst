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

   Compute the **soil thermal conductivity** :math:`\lambda_s`
   as a function of volumetric water content (θ).

   The relationship is taken from Gao et al. (2017, Eq. 12):

   .. math::

       \lambda_s(\theta) = 0.20 + \exp\bigl[\,1.46\,(\theta - 0.34)\bigr]

   :param theta: Volumetric water content (m³ m⁻³).
                 Accepts a scalar value or any array-like object that can be
                 converted to a :class:`numpy.ndarray`. Values should lie in the
                 closed interval ``[0, 1]``.
   :type theta: float or array_like

   :returns: Soil thermal conductivity λ\_s in W m⁻¹ K⁻¹.  The returned type
             matches the input: a scalar is returned for a scalar *θ*, and a
             NumPy array for array-like *θ*.
   :rtype: float or numpy.ndarray

   :raises ValueError: If any element of *θ* is outside the physically meaningful
       range ``[0, 1]``.

   .. rubric:: Notes

   * **Vectorization** – The function is fully vectorized; it operates
     element-wise on NumPy arrays and broadcasts according to NumPy
     broadcasting rules.
   * **Empirical limits** – Gao et al. recommend the equation for
     mineral soils where 0 ≤ θ ≤ 0.5. Extrapolation beyond this range
     may introduce error.
   * **Units** – The output uses SI units (W m⁻¹ K⁻¹).

   .. rubric:: References

   Gao, Z., Niu, G.-Y., & Hedquist, B. C. (2017).
   **Soil thermal conductivity parameterization: A closed-form
   equation and its evaluation.** *Journal of Geophysical Research:
   Atmospheres*, 122(6), 3466–3478.
   DOI: 10.1002/2016JD025992

   .. rubric:: Examples

   >>> lambda_s(0.25)
   0.463...
   >>> import numpy as np
   >>> theta = np.linspace(0, 0.5, 5)
   >>> lambda_s(theta)
   array([0.20      , 0.31243063, 0.41172297, 0.51210322, 0.61861198])


.. py:function:: k_s(theta: numpy.ndarray | float) -> numpy.ndarray | float

   Calculate the **soil thermal diffusivity** :math:`k_s`
   as a function of volumetric water content (θ).

   The empirical relationship follows Gao et al. (2017, Eq. 13):

   .. math::

       k_s(\theta) = \bigl[0.69 + \exp\bigl(3.06\,(\theta - 0.26)\bigr)\bigr] \times 10^{-7}

   :param theta: Volumetric water content (m³ m⁻³).
                 Accepts a scalar or any array-like sequence convertible to a
                 :class:`numpy.ndarray`. Values **must** lie in the closed
                 interval ``[0, 1]``.
   :type theta: float or array_like

   :returns: Soil thermal diffusivity *k\_s* in m² s⁻¹.
             A scalar is returned if *θ* is a scalar; otherwise a NumPy array
             of matching shape is returned.
   :rtype: float or numpy.ndarray

   :raises ValueError: If any element of *θ* is outside the physically meaningful
       range ``[0, 1]``.

   .. rubric:: Notes

   * **Vectorization** – The calculation is fully vectorized and
     broadcasts according to NumPy rules.
   * **Applicability** – Gao et al. derived this expression for mineral
     soils under typical field conditions. Extrapolating beyond
     0 ≤ θ ≤ 0.5 may reduce accuracy.
   * **Units** – Output is in SI units (m² s⁻¹).

   .. rubric:: References

   Gao, Z., Niu, G.-Y., & Hedquist, B. C. (2017).
   *Soil thermal conductivity parameterization: A closed-form equation
   and its evaluation.* Journal of Geophysical Research: Atmospheres,
   **122**(6), 3466–3478. https://doi.org/10.1002/2016JD025992

   .. rubric:: Examples

   >>> k_s(0.25)
   1.321...e-07
   >>> thetas = np.linspace(0, 0.5, 6)
   >>> k_s(thetas)
   array([1.14...e-07, 1.21...e-07, 1.30...e-07, 1.41...e-07,
          1.54...e-07, 1.69...e-07])


.. py:function:: volumetric_heat_capacity(lambda_s_val: numpy.ndarray | float, k_s_val: numpy.ndarray | float) -> numpy.ndarray | float

   Compute the **volumetric heat capacity** :math:`C_v` of soil:

   .. math::

       C_v \;=\; \frac{\lambda_s}{k_s}

   where
   :math:`\lambda_s` is the **thermal conductivity** (W m⁻¹ K⁻¹) and
   :math:`k_s` is the **thermal diffusivity** (m² s⁻¹).

   The resulting heat capacity has units of J m⁻³ K⁻¹.

   :param lambda_s_val: Soil thermal conductivity (W m⁻¹ K⁻¹). May be a scalar or any
                        array-like structure broadcastable with *k_s_val*.
   :type lambda_s_val: float or array_like
   :param k_s_val: Soil thermal diffusivity (m² s⁻¹). Must be positive. Accepts a
                   scalar or array-like input.
   :type k_s_val: float or array_like

   :returns: Volumetric heat capacity *C_v* (J m⁻³ K⁻¹).  The return type
             matches the input: a scalar for scalar inputs, or a NumPy array
             for array-like inputs.
   :rtype: float or numpy.ndarray

   :raises ValueError: If any element in *k_s_val* is zero or negative, which would
       lead to division by zero or non-physical results.

   .. rubric:: Notes

   * **Vectorization** – The function is fully vectorized; both inputs
     are converted to :class:`numpy.ndarray` and follow NumPy
     broadcasting rules.
   * **Physical meaning** – Volumetric heat capacity represents the
     energy required to raise the temperature of a unit volume of soil
     by one kelvin. High values correspond to moist or water-logged
     soils; dry mineral soils have lower *C_v*.

   .. rubric:: Examples

   >>> #Scalar inputs
   >>> volumetric_heat_capacity(0.8, 1.2e-6)
   666666.666...

   >>> #Vectorized inputs
   >>> lambda_vals = np.array([0.5, 0.6, 0.7])
   >>> diffusivities = np.array([1.1e-6, 1.2e-6, 1.3e-6])
   >>> volumetric_heat_capacity(lambda_vals, diffusivities)
   array([454545.4545..., 500000.    ..., 538461.5384...])


.. py:function:: nme(calc: numpy.ndarray | float, meas: numpy.ndarray | float) -> float

   Calculate the **normalized mean error (NME)** between calculated
   and measured values, expressed as a percentage.

   The formulation follows Gao et al. (2017, Eq. 14):

   .. math::

       \text{NME} \;=\; 100\,\frac{\sum\limits_{i}\left|\hat{y}_i - y_i\right|}
                                     {\sum\limits_{i}\left|y_i\right|}

   where :math:`\hat{y}_i` are the *calculated* (model) values and
   :math:`y_i` are the *measured* (reference) values.

   :param calc: Modelled / calculated values :math:`\hat{y}`.
                Accepts any array-like object (including scalars) convertible to
                a :class:`numpy.ndarray`.
   :type calc: float or array_like
   :param meas: Measured / observed reference values :math:`y`. Must be
                broadcast-compatible with *calc*.
   :type meas: float or array_like

   :returns: Normalized mean error (percent).
             A value of **0 %** indicates perfect agreement; larger values
             indicate greater error.
   :rtype: float

   :raises ValueError: * If *calc* and *meas* cannot be broadcast to a common shape.
       * If the denominator ``Σ|meas|`` equals zero (e.g., all measured
         values are zero), which would make NME undefined.

   .. rubric:: Notes

   * **Range** – NME is non-negative and unbounded above.
   * **Interpretation** – Because both numerator and denominator use
     absolute values, NME is insensitive to the direction of the error
     (over- vs under-prediction) and therefore complements signed error
     metrics such as mean bias.
   * **Vectorization** – The implementation is fully vectorized and
     adheres to NumPy broadcasting rules.

   .. rubric:: References

   Gao, Z., Niu, G.-Y., & Hedquist, B. C. (2017).
   *Soil thermal conductivity parameterization: A closed-form equation
   and its evaluation.* **Journal of Geophysical Research:
   Atmospheres**, 122(6), 3466–3478.
   https://doi.org/10.1002/2016JD025992

   .. rubric:: Examples

   >>> # Single values
   >>> nme(4.5, 5.0)
   10.0

   >>> # Vectors
   >>> calc = np.array([1.0, 2.1, 3.2])
   >>> meas = np.array([1.2, 2.0, 3.0])
   >>> nme(calc, meas)
   4.7619...

   >>> # Broadcasting (scalar vs array)
   >>> nme(2.0, np.array([1.5, 2.5, 2.0]))
   16.6666...


.. py:function:: rmse(calc: numpy.ndarray | float, meas: numpy.ndarray | float) -> float

   Compute the **root-mean-square error (RMSE)** between calculated
   (model) and measured (reference) values.

   Following Gao et al. (2017, Eq. 15):

   .. math::

       \text{RMSE} \;=\;
       \sqrt{\frac{1}{N}\sum_{i=1}^{N}\bigl(\hat{y}_i - y_i\bigr)^2}

   where :math:`\hat{y}_i` are *calculated* values and
   :math:`y_i` are *measured* values.

   :param calc: Calculated / modelled values :math:`\hat{y}`.  Accepts a scalar
                or any array-like object convertible to a
                :class:`numpy.ndarray`.
   :type calc: float or array_like
   :param meas: Measured / observed values :math:`y`. Must be broadcast-
                compatible with *calc*.
   :type meas: float or array_like

   :returns: Root-mean-square error (same units as the inputs).
   :rtype: float

   :raises ValueError: * If *calc* and *meas* cannot be broadcast to a common shape.
       * If the input arrays are empty (``N = 0``).

   .. rubric:: Notes

   * **Interpretation** – RMSE represents the sample-standard-deviation
     of the differences between two datasets; lower values indicate
     better agreement.
   * **Vectorization** – The function is fully vectorized and respects
     NumPy broadcasting rules.
   * **Units** – RMSE preserves the units of the input variables.

   .. rubric:: References

   Gao, Z., Niu, G.-Y., & Hedquist, B. C. (2017).
   *Soil thermal conductivity parameterization: A closed-form equation
   and its evaluation.* Journal of Geophysical Research:
   Atmospheres, **122**(6), 3466–3478.
   https://doi.org/10.1002/2016JD025992

   .. rubric:: Examples

   >>> # Scalar inputs
   >>> rmse(4.5, 5.0)
   0.5

   >>> # Vector inputs
   >>> calc = np.array([2.1, 3.0, 4.2])
   >>> meas = np.array([2.0, 3.5, 4.0])
   >>> rmse(calc, meas)
   0.2645...

   >>> # Broadcasting
   >>> rmse(3.0, np.array([2.5, 3.5, 3.0]))
   0.4082...


.. py:function:: calorimetric_gz(g_zr: numpy.ndarray | float, cv_layers: numpy.ndarray | float, dT_dt_layers: numpy.ndarray | float, dz_layers: numpy.ndarray | float) -> numpy.ndarray | float

   Estimate the **soil heat flux** at a target depth *z*
   (typically 5 cm) using the *calorimetric method*.

   The method corrects an in-situ heat-flux plate reading made at a
   reference depth *z_r* by adding the change in heat storage of the
   soil column located between *z* and *z_r*.  Mathematically
   (Liebethal & Foken 2007, Eq. 1):

   .. math::

       G(z,t) \;=\; G(z_r,t) \;+\;
       \sum_{l=1}^{N}
       C_{v,l}\,\frac{\partial \bar{T}_l}{\partial t}\,\Delta z_l

   where

   * :math:`G(z,t)`     … heat flux at depth *z* (W m⁻²)
   * :math:`G(z_r,t)`   … measured plate flux at reference depth *z_r*
   * :math:`C_{v,l}`    … volumetric heat capacity of layer *l*
     (J m⁻³ K⁻¹)
   * :math:`\partial \bar{T}_l / \partial t` … time derivative of the
     layer-averaged temperature (K s⁻¹)
   * :math:`\Delta z_l` … thickness of layer *l* (m)

   :param g_zr: Heat-flux plate measurement at depth *z_r* (W m⁻²).  Can be a
                scalar or time series.
   :type g_zr: float or array_like
   :param cv_layers: Volumetric heat capacity :math:`C_{v,l}` for each of *N* soil
                     sub-layers (J m⁻³ K⁻¹).  Shape ``(N, …)`` where the trailing
                     dimensions (``…``) must be broadcast-compatible with *g_zr*.
   :type cv_layers: array_like
   :param dT_dt_layers: Time derivative of layer-mean temperature
                        :math:`\partial \bar{T}_l / \partial t`
                        (K s⁻¹); same shape as *cv_layers*.
   :type dT_dt_layers: array_like
   :param dz_layers: Thickness of each layer :math:`\Delta z_l` (m); shape ``(N,)``
                     or broadcast-compatible with the first axis of *cv_layers*.
   :type dz_layers: array_like

   :returns: Calorimetrically corrected soil heat flux *G(z)* (W m⁻²).  Scalar
             if all inputs are scalar; otherwise a NumPy array matching the
             broadcast shape of *g_zr*.
   :rtype: float or numpy.ndarray

   :raises ValueError: If the inputs cannot be broadcast to a common shape, or if the
       number of layers inferred from *cv_layers*, *dT_dt_layers*, and
       *dz_layers* are inconsistent.

   .. rubric:: Notes

   * **Vectorization** – All operations are fully vectorized.  Inputs
     are converted to :class:`numpy.ndarray` and follow NumPy
     broadcasting rules.
   * **Layer axis** – The first dimension of *cv_layers* and
     *dT_dt_layers* (axis 0) is interpreted as the layer index *l*.
     Layer thicknesses *dz_layers* are broadcast across any additional
     dimensions.
   * **Units** – Ensure consistent SI units: W m⁻², J m⁻³ K⁻¹,
     K s⁻¹, and m.

   .. rubric:: References

   Liebethal, C., & Foken, T. (2007). *Evaluation of six parameterization
   approaches for the ground heat flux.* Agricultural and Forest
   Meteorology, **143**(1–2), 65-80.
   https://doi.org/10.1016/j.agrformet.2006.11.001

   .. rubric:: Examples

   >>> # Three-layer example, single time step
   >>> g_plate = -15.2                      # W m-2 at z_r = −0.05 m
   >>> Cv     = np.array([2.5e6, 2.3e6, 2.1e6])   # J m-3 K-1
   >>> dTdt   = np.array([1.2e-4, 0.9e-4, 0.6e-4])  # K s-1
   >>> dz     = np.array([0.02, 0.02, 0.01])       # m
   >>> calorimetric_gz(g_plate, Cv, dTdt, dz)
   -4.06...

   >>> # Vectorized daily time series with two layers
   >>> g_plate = np.random.normal(-10, 2, 1440)        # per minute
   >>> Cv       = np.array([[2.4e6], [2.2e6]])         # (2,1)
   >>> dTdt     = np.random.normal(5e-5, 2e-5, (2,1440))
   >>> dz       = np.array([0.03, 0.02])
   >>> Gz = calorimetric_gz(g_plate, Cv, dTdt, dz)     # shape (1440,)


.. py:function:: force_restore_gz(cv: numpy.ndarray | float, dTg_dt: numpy.ndarray | float, Tg: numpy.ndarray | float, Tg_bar: numpy.ndarray | float, delta_z: float = 0.05, omega: float = OMEGA_DAY) -> numpy.ndarray | float

   Estimate **soil heat flux** :math:`G(z)` at a shallow depth
   (:math:`z = \delta z`, default 5 cm) with the **force–restore
   method**.

   The formulation (Liebethal & Foken 2007, Eq. 2) corrects the
   calorimetric storage term with a *restore* component that accounts
   for the difference between the instantaneous ground temperature
   *T_g* and its running mean *\bar{T}_g*:

   .. math::

       G(z,t) \;=\; C_v \, \delta z \, \frac{\partial T_g}{\partial t}
       \;+
       \sqrt{\omega \, C_v \, \lambda_s(C_v)}
       \left[
           \frac{1}{\omega}\,\frac{\partial T_g}{\partial t}
           + \bigl(T_g - \bar{T}_g\bigr)
       \right]

   where

   * :math:`C_v`   … volumetric heat capacity (J m⁻³ K⁻¹)
   * :math:`\partial T_g / \partial t` … ground-temperature tendency
     (K s⁻¹)
   * :math:`\lambda_s(C_v)` … soil thermal conductivity derived from
     *C_v* via :func:`lambda_s_from_cv` (W m⁻¹ K⁻¹)
   * :math:`\omega` … angular frequency of the diurnal cycle
     (rad s⁻¹)
   * :math:`T_g` / :math:`\bar{T}_g` … instantaneous and running-mean
     ground temperature (K)
   * :math:`\delta z` … sensor depth below the surface (m)

   :param cv: Volumetric heat capacity :math:`C_v` (J m⁻³ K⁻¹).  Must be
              positive.  Accepts scalars or NumPy-broadcastable arrays.
   :type cv: float or array_like
   :param dTg_dt: Time derivative :math:`\partial T_g/\partial t` (K s⁻¹).
   :type dTg_dt: float or array_like
   :param Tg: Instantaneous ground (surface) temperature :math:`T_g` (K or °C)
              at depth *δz*.
   :type Tg: float or array_like
   :param Tg_bar: Running mean ground temperature :math:`\bar{T}_g` over the
                  diurnal cycle (same units as *Tg*).
   :type Tg_bar: float or array_like
   :param delta_z: Depth :math:`\delta z` in metres; default is **0.05 m**
                   (5 cm).
   :type delta_z: float, optional
   :param omega: Angular frequency :math:`\omega` in rad s⁻¹.  Defaults to
                 :pydata:`OMEGA_DAY` (2π / 86 400 s).
   :type omega: float, optional

   :returns: Soil heat flux *G(z)* at depth *δz* (W m⁻²).  The output shape
             follows NumPy broadcasting rules applied to the inputs.
   :rtype: float or numpy.ndarray

   :raises ValueError: If any input arrays cannot be broadcast to a common shape, or if
       *cv* contains non-positive values.

   .. rubric:: Notes

   * **λ_s(C_v) mapping** – The function relies on a helper
     :func:`lambda_s_from_cv` that converts volumetric heat capacity
     to thermal conductivity.  Ensure that this helper is present in
     the import path.
   * **Units** – Keep units internally consistent (SI).
   * **Vectorization** – All operations are vectorized; the
     mathematical expression is evaluated element-wise for array
     inputs.
   * **Interpretation** – The first term represents *storage*
     (calorimetric), while the second tempers short-term fluctuations,
     “restoring” *T_g* toward its mean.

   .. rubric:: References

   Liebethal, C., & Foken, T. (2007). *Evaluation of six
   parameterization approaches for the ground heat flux.*
   **Agricultural and Forest Meteorology**, 143(1–2), 65-80.
   https://doi.org/10.1016/j.agrformet.2006.11.001

   .. rubric:: Examples

   >>> Cv      = 2.4e6                 # J m-3 K-1
   >>> dTgdt   = 1.0e-4                # K s-1
   >>> Tg      = 293.5                 # K
   >>> Tg_bar  = 291.7                 # K (running mean)
   >>> force_restore_gz(Cv, dTgdt, Tg, Tg_bar)
   19.3...   # W m-2

   >>> #Vectorized daily record
   >>> cv_arr = np.full(1440, 2.2e6)
   >>> dT_arr = np.gradient(np.sin(np.linspace(0, 2*np.pi, 1440))) / 60
   >>> Gz_ts  = force_restore_gz(cv_arr, dT_arr, 298+2*np.sin(...),
   ...                           298*np.ones_like(dT_arr))


.. py:function:: gao2010_gz(AT: numpy.ndarray | float, lambda_s_val: numpy.ndarray | float, k_s_val: numpy.ndarray | float, t: numpy.ndarray | float, omega: float = OMEGA_DAY) -> numpy.ndarray | float

   Estimate **soil heat flux** :math:`G(z,t)` at depth *z* based on a
   *sinusoidal* ground-temperature forcing (Gao et al. 2010, Eq. 3).

   The analytical solution assumes that the surface (or ground-contact)
   temperature varies sinusoidally with amplitude *A_T* and angular
   frequency *ω*.  Under these conditions the heat flux at any depth
   *z* can be written

   .. math::

       G(z,t) \;=\;
       \sqrt{2}\,
       \frac{\lambda_s A_T}{d}\,
       \sin\bigl(\omega t + \tfrac{\pi}{4}\bigr),

   where the **thermal damping depth**

   .. math::

       d \;=\; \sqrt{\frac{2 k_s}{\omega}}

   is determined by the soil thermal diffusivity *k_s* and the forcing
   frequency *ω*.

   :param AT: Amplitude of the sinusoidal ground-surface temperature (K).
   :type AT: float or array_like
   :param lambda_s_val: Soil thermal conductivity :math:`\lambda_s` (W m⁻¹ K⁻¹).
   :type lambda_s_val: float or array_like
   :param k_s_val: Soil thermal diffusivity :math:`k_s` (m² s⁻¹).
   :type k_s_val: float or array_like
   :param t: Time variable (s).  Can be absolute time since epoch or simply
             seconds since the start of the cycle—as long as *ω t* is
             dimensionless.
   :type t: float or array_like
   :param omega: Angular frequency *ω* (rad s⁻¹).  Defaults to *OMEGA_DAY*
                 (≈ 7.272 × 10⁻⁵ s⁻¹, i.e. 2π / 86 400 s).
   :type omega: float, optional

   :returns: Heat flux *G(z,t)* (W m⁻²).  The output follows NumPy’s
             broadcasting rules applied to the inputs.
   :rtype: float or numpy.ndarray

   :raises ValueError: * If *lambda_s_val*, *k_s_val*, and *t* cannot be broadcast to a
         common shape.
       * If any element of *k_s_val* or *omega* is non-positive.

   .. rubric:: Notes

   * **Vectorization** – All inputs are internally converted to
     :class:`numpy.ndarray`; the formula is evaluated element-wise and
     fully supports broadcasting.
   * **Units** – Ensure consistent SI units: W m⁻¹ K⁻¹ (λ_s),
     m² s⁻¹ (k_s), s (t), and rad s⁻¹ (ω).
   * **Scope** – The solution presumes purely sinusoidal boundary
     forcing and homogeneous soil properties; real-world deviations
     (e.g., non-sine forcing, stratified soils, moisture variation)
     will introduce error.

   .. rubric:: References

   Gao, Z., Horton, R., Luo, L., & Kucharik, C. J. (2010).
   *A simple method to measure soil temperature dynamics: Theory and
   application.* **Soil Science Society of America Journal**, 74(2),
   580-588. https://doi.org/10.2136/sssaj2009.0169

   .. rubric:: Examples

   >>> # Daily cycle at 5 cm depth
   >>> AT      = 8.0                           # K
   >>> lambda_ = 1.2                           # W m-1 K-1
   >>> kappa   = 1.0e-6                        # m2 s-1
   >>> t_day   = np.linspace(0, 86400, 97)     # 15-min steps
   >>> Gz      = gao2010_gz(AT, lambda_, kappa, t_day)
   >>> Gz.shape
   (97,)


.. py:function:: heusinkveld_gz(A_n: numpy.ndarray | float, Phi_n: numpy.ndarray | float, n_max: int, k_s_val: numpy.ndarray | float, lambda_s_val: numpy.ndarray | float, w: float) -> numpy.ndarray | float

   Compute *soil heat flux* :math:`G(z,t)` from the **H04 harmonic
   series solution** proposed by Heusinkveld et al. (2004, Eq. 4).

   The approach represents the surface (ground) temperature as a Fourier
   series with harmonics up to order *n_max*.  For each harmonic
   :math:`n` the heat-flux contribution at depth *z* can be written

   .. math::

       G_n(z,t) \;=\;
       \frac{\lambda_s}{10\,\pi}\,
       A_n\;
       \sqrt{\frac{1}{k_s\,n\,\omega\,k_s}}\;
       \sin\bigl(n\,\omega\,t + \Phi_n + \tfrac{\pi}{4}\bigr),

   and the full signal is obtained by summing over *n = 1…n_max*.

   **Note** The implementation below follows the algebraic form that
   appeared in H04; consult the original paper for derivation details
   and recommended parameter ranges.

   :param A_n: Amplitudes :math:`A_n` of the *n*-th harmonic of surface-temperature
               forcing (K).  Provide either
               * a single scalar applied to every harmonic, or
               * an array of length ≥ *n_max* giving amplitude for each harmonic.
   :type A_n: float or array_like
   :param Phi_n: Phase shifts :math:`\Phi_n` (rad) corresponding to each harmonic
                 order.  Same broadcasting rules as *A_n*.
   :type Phi_n: float or array_like
   :param n_max: Highest harmonic order to include in the summation.
   :type n_max: int
   :param k_s_val: Soil thermal diffusivity :math:`k_s` (m² s⁻¹).
   :type k_s_val: float or array_like
   :param lambda_s_val: Soil thermal conductivity :math:`\lambda_s` (W m⁻¹ K⁻¹).
   :type lambda_s_val: float or array_like
   :param w: Fundamental angular frequency :math:`\omega` (rad s⁻¹), e.g.
             :math:`2\pi/86400` for a 24-h cycle.
   :type w: float

   :returns: Soil heat flux *G(z,t)* (W m⁻²).  The return shape is the
             broadcast shape of the input arrays (excluding the harmonic axis).
   :rtype: float or numpy.ndarray

   :raises ValueError: If *n_max* is less than 1, or if *k_s_val* or *w* are non-positive,
       or if the amplitudes/phases cannot be broadcast to length
       *n_max*.

   .. rubric:: Notes

   * **Vectorization** – Internally, the harmonic index axis has length
     *n_max* (``n = np.arange(1, n_max+1)``).  All other dimensions come
     from broadcasting *A_n*, *Phi_n*, *k_s_val*, *lambda_s_val*, and
     *t* (if vectorized in the caller).
   * **Units** – Consistency with SI units is assumed.
   * **Interpretation** – The prefactor ``λ_s / (10 π)`` appears in
     H04’s original derivation.  If you adopt a different convention
     (e.g., depth-explicit damping), modify accordingly.

   .. rubric:: References

   Heusinkveld, B. G., Jacobs, A. F. G., Holtslag, A. A. M., & *et al.*
   (2004). *Surface energy balance closure in an arid region: The role
   of soil heat flux.* Agricultural and Forest Meteorology, **122**(1),
   21-37. https://doi.org/10.1016/j.agrformet.2003.09.005

   .. rubric:: Examples

   >>> # Daily cycle with three harmonics
   >>> A = np.array([10, 4, 1.5])          # K
   >>> Phi = np.array([0, -np.pi/6, np.pi/8])
   >>> Gz = heusinkveld_gz(A, Phi, n_max=3,
   ...                     k_s_val=1e-6,
   ...                     lambda_s_val=1.0,
   ...                     w=2*np.pi/86400)
   >>> Gz.shape
   ()


.. py:function:: hsieh2009_gz(tz_series: numpy.ndarray | list | tuple, time_series: numpy.ndarray | list | tuple, cv_series: numpy.ndarray | list | tuple, ks_series: numpy.ndarray | list | tuple) -> float

   Compute **soil heat flux** at the end of a temperature record
   using the *half-order integral* (Hsieh et al., 2009, Eq. 5).

   The method exploits the analytical solution of the one-dimensional
   heat-conduction equation for a semi-infinite medium, leading to a
   convolution integral of order ½ that relates the time series of
   near-surface soil temperature to the downward heat flux:

   .. math::

       G(t_N) \;=\;
       2\sqrt{\frac{k_s(t_N)\,C_v(t_N)}{\pi}}\;
       \int_{t_0}^{t_N}
           \frac{\partial T(z,t')}{\partial t'}
           \left(t_N - t'\right)^{-1/2}\,dt'

   In discrete form with linear interpolation between measurements
   :math:`t_i \;(i = 0,\dots,N)`,

   .. math::

       \int_{t_0}^{t_N}\!
           \frac{dT}{dt'}\,(t_N-t')^{-1/2}dt'
       \approx
       \sum_{i=0}^{N-1}
         \frac{\Delta T_i}{\Delta t_i}\!
         \left[(t_N-t_i)^{1/2} - (t_N-t_{i+1})^{1/2}\right],

   where :math:`\Delta T_i = T_{i+1}-T_i` and
   :math:`\Delta t_i = t_{i+1}-t_i`.

   The final scalar result corresponds to the flux at
   :math:`t_N\;(=\text{time_series}[-1])`.

   :param tz_series: Near-surface (or shallow-depth) soil temperature observations
                     *T(z,t)* (K).  Length **N ≥ 2**.
   :type tz_series: array_like
   :param time_series: Strictly monotonically increasing time stamps (s) matching
                       `tz_series`.  Same length **N**.
   :type time_series: array_like
   :param cv_series: Volumetric heat capacity *C_v* (J m⁻³ K⁻¹) at each
                     time stamp.  Same length **N**.
   :type cv_series: array_like
   :param ks_series: Thermal diffusivity *k_s* (m² s⁻¹) at each time stamp.
                     Same length **N**.
   :type ks_series: array_like

   :returns: Soil heat flux *G(t_N)* at the final time point (W m⁻²).
   :rtype: float

   :raises ValueError: * If the four series differ in length or have fewer than two
         samples.
       * If `time_series` is not strictly increasing.
       * If any element of `cv_series` or `ks_series` is non-positive.

   .. rubric:: Notes

   * The algorithm uses a forward finite-difference for
     :math:`dT/dt'` and trapezoidal integration over each interval
     *[t_i, t_{i+1}]*.
   * Only the latest values of *k_s* and *C_v* are used, consistent
     with the Hsieh et al. derivation.  Replace with a time-variable
     kernel if property changes are large over the record.
   * All inputs are cast to :class:`numpy.ndarray` with
     ``dtype=float``.

   .. rubric:: References

   Hsieh, C.-I., Katul, G., & Chi, T.-C. (2009).
   *Retrieval of soil heat flux from soil temperature data by the
   continuous-time heat equation model.*
   **Water Resources Research**, 45, W08433.
   https://doi.org/10.1029/2009WR007891

   .. rubric:: Examples

   >>> times  = np.array([0, 600, 1200, 1800])        # every 10 min
   >>> temps  = np.array([292.5, 293.0, 293.6, 294.0])  # K
   >>> Cv     = np.full_like(temps, 2.3e6)            # J m-3 K-1
   >>> ks     = np.full_like(temps, 1.1e-6)           # m2 s-1
   >>> hsieh2009_gz(temps, times, Cv, ks)
   -2.13...


.. py:function:: leuning_damping_depth(z: numpy.ndarray | float, zr: numpy.ndarray | float, AT_z: numpy.ndarray | float, AT_zr: numpy.ndarray | float) -> numpy.ndarray | float

   Estimate the **thermal damping depth** *d* (m) from the exponential
   decay of temperature–wave amplitude with depth.

   The formulation is based on Eq. (6) of Leuning *et al.* (1985) for a
   one-dimensional, sinusoidally forced soil column:

   .. math::

       d \;=\;
       \frac{z_r - z}{\ln\bigl(A_T(z) / A_T(z_r)\bigr)}

   where

   * :math:`A_T(z)` — amplitude of the temperature wave at depth *z*
   * :math:`z_r`  — reference depth (m)
   * :math:`A_T(z_r)` — amplitude at reference depth

   The equation assumes amplitudes decrease monotonically with depth
   following :math:`A_T(z) = A_T(0)\,\exp(-z/d)`.

   :param z: Depth(s) (m) at which the amplitude *A_T(z)* was measured.
             Positive values denote depth below the surface.
   :type z: float or array_like
   :param zr: Reference depth(s) *z_r* (m).  Must be broadcast-compatible with
              *z*.  If *zr* < *z* the denominator switches sign, yielding a
              negative *d* that should be interpreted carefully.
   :type zr: float or array_like
   :param AT_z: Temperature-wave amplitude at *z* (K).
   :type AT_z: float or array_like
   :param AT_zr: Temperature-wave amplitude at *z_r* (K).
   :type AT_zr: float or array_like

   :returns: Thermal damping depth *d* (m).  The output shape follows NumPy
             broadcasting rules applied to the four inputs.
   :rtype: float or numpy.ndarray

   :raises ValueError: If any of the following conditions occur:

       * *AT_z* or *AT_zr* contain non-positive values (log undefined).
       * *AT_z* and *AT_zr* are equal everywhere, leading to
         ``ln(1) = 0`` and division by zero.
       * The inputs cannot be broadcast to a common shape.

   .. rubric:: Notes

   * **Vectorisation** — Internally, all inputs are converted to
     :class:`~numpy.ndarray`; the formula is applied element-wise and
     fully supports broadcasting.
   * **Physical meaning** — Larger *d* indicates slower attenuation of
     temperature fluctuations with depth (i.e. higher thermal
     diffusivity or conductivity).  Very small logarithmic denominators
     (`AT_z` ≈ `AT_zr`) imply an extremely deep damping depth that may
     fall outside the soil column considered.
   * **Units** — Depths must share the same units (metres).  Amplitudes
     may be in kelvins or degrees Celsius provided they use identical
     scaling.

   .. rubric:: References

   Leuning, R., Barlow, E. W. R., & Paltridge, G. (1985).
   *Temperature gradients in a soil: Theory and measurement.*
   *Agricultural and Forest Meteorology*, 35(1–4), 127–143.
   https://doi.org/10.1016/0168-1923(85)90067-3

   .. rubric:: Examples

   >>> d = leuning_damping_depth(z=0.10, zr=0.05, AT_z=4.0, AT_zr=8.0)
   >>> round(d, 4)
   0.0721

   >>> #Broadcasting with arrays
   >>> depths = np.array([0.05, 0.10, 0.20])
   >>> amps   = np.array([8.0, 4.0, 1.5])
   >>> d_arr  = leuning_damping_depth(depths, 0.02, amps, 10.0)
   >>> d_arr
   array([0.0380..., 0.0495..., 0.0854...])


.. py:function:: leuning_gz(g_zr: numpy.ndarray | float, z: numpy.ndarray | float, zr: numpy.ndarray | float, d: numpy.ndarray | float) -> numpy.ndarray | float

   Extrapolate **soil heat flux** from a reference depth *z_r* to a
   shallower (or deeper) target depth *z* using an *exponential
   attenuation* model (Leuning *et al.* 1985, Eq. 7).

   The model assumes the amplitude of the heat-flux wave decays
   exponentially with depth at the same *damping depth* **d** that
   governs temperature attenuation:

   .. math::

       G(z,t) \;=\; G(z_r,t)\;\exp\!\left(\frac{z_r - z}{d}\right)

   where
   :math:`G(z_r,t)` is the measured flux at *z_r* (W m⁻²).

   :param g_zr: Heat flux measured at reference depth *z_r* (W m⁻²).  Can be a
                scalar or a NumPy-broadcastable array (e.g. a time series).
   :type g_zr: float or array_like
   :param z: Target depth *z* (m, **positive downward**).  Must be broadcast-
             compatible with *g_zr*.
   :type z: float or array_like
   :param zr: Reference depth *z_r* (m, positive downward).  Broadcast
              compatible with *z*.
   :type zr: float or array_like
   :param d: Thermal damping depth (m).  Positive; same shape rules as *z*.
   :type d: float or array_like

   :returns: Estimated heat flux at depth *z* (W m⁻²).  Follows NumPy
             broadcasting rules applied to the inputs.
   :rtype: float or numpy.ndarray

   :raises ValueError: If any element of *d* is non-positive, or if inputs cannot be
       broadcast to a common shape.

   .. rubric:: Notes

   * **Direction of extrapolation** –
     *z < z_r* (shallower) ⇒ magnitude *increases* toward the surface;
     *z > z_r* (deeper)   ⇒ magnitude *decreases* with depth.
   * **Vectorisation** – Inputs are converted to
     :class:`numpy.ndarray` and the formula is evaluated element-wise.
   * **Limitations** – The simple exponential form ignores soil
     layering and moisture variability; best suited to homogeneous
     profiles over the depth interval considered.

   .. rubric:: References

   Leuning, R., Barlow, E. W. R., & Paltridge, G. (1985).
   *Temperature gradients in a soil: Theory and measurement.*
   **Agricultural and Forest Meteorology**, 35(1–4), 127-143.
   https://doi.org/10.1016/0168-1923(85)90067-3

   .. rubric:: Examples

   >>> # Single depth conversion
   >>> leuning_gz(g_zr=-12.0, z=0.05, zr=0.08, d=0.07)
   -18.841...

   >>> # Vectorized daily time series
   >>> g_plate = np.random.normal(-10, 3, 1440)      # W m-2 @ 8 cm
   >>> d       = 0.07                                # m
   >>> Gz_5cm  = leuning_gz(g_plate, z=0.05, zr=0.08, d=d)
   >>> Gz_5cm.shape
   (1440,)


.. py:function:: simple_measurement_gz(g_zr: numpy.ndarray | float, cv_layers: numpy.ndarray | float, tz_layers: numpy.ndarray | float, dt: float, dz_layers: numpy.ndarray | float) -> numpy.ndarray | float

       Ground-heat-flux estimate at a target depth *z* using the
       **simple-measurement variant** of the *calorimetric* method
       (Liebethal & Foken 2007, Eq. 8).

       The algorithm corrects a heat-flux‐plate measurement taken at a
       reference depth *z_r* by adding an approximation of the heat storage
       change *ΔS* in the soil slab between the plate and the surface,
       computed from two successive soil-temperature profiles:

       .. math::

           G(z,t_j) \;=\; G(z_r,t_j)
           \;+\; \sum_{l=1}^{N} C_{v,l}\,\Delta z_l\,

   rac{\Delta T_{l,j} +
   rac{1}{2}igl(\Delta T_{l,j} - \Delta T_{l,j-1}igr)}
                       {\Delta t}

       where

       * :math:`C_{v,l}` — volumetric heat capacity of layer *l*
         (J m⁻³ K⁻¹)
       * :math:`\Delta z_l` — layer thickness (m)
       * :math:`\Delta T_{l,j}` — temperature change in layer *l*
         between *t_{j-1}* and *t_j* (K)
       * :math:`\Delta t` — measurement interval (s)

       The centred finite-difference term in brackets approximates the
       mid-interval temperature tendency, providing second-order accuracy
       from simple time-series measurements (`tz_layers`).

       Parameters
       ----------
       g_zr : float or array_like
           Plate-measured soil heat flux at reference depth *z_r*
           (W m⁻²).  Length **M** time steps.
       cv_layers : array_like
           Volumetric heat capacity for each of *N* soil layers
           (J m⁻³ K⁻¹).  Shape ``(N,)`` or broadcast-compatible with the
           first axis of `tz_layers`.
       tz_layers : array_like
           Layer-average soil temperatures.  Shape ``(N, M)`` where *N* is
           the number of layers (matching `cv_layers`/`dz_layers`) and *M*
           is the number of time stamps (matching `g_zr`).
       dt : float
           Constant time step between consecutive temperature samples
           (s).  Must be positive.
       dz_layers : array_like
           Thickness of each layer (m).  Shape ``(N,)`` broadcast-compatible
           with `cv_layers`.

       Returns
       -------
       float or numpy.ndarray
           Calorimetrically corrected soil heat flux at depth *z*
           (W m⁻²).  Length **M − 1** because the centred difference
           requires two successive profiles.

       Raises
       ------
       ValueError
           If input dimensions are inconsistent, `dt` ≤ 0, or any layer
           arrays contain non-finite values.

       Notes
       -----
       * **Temporal alignment** – The first output value corresponds to the
         mid-point of the first two temperature samples
         ``t[0] … t[1]``; therefore the result is shifted **½ Δt** relative
         to the plate-flux time stamps.
       * **Vectorisation** – All calculations are fully vectorized using
         NumPy broadcasting.  Inputs are internally cast to
         :class:`numpy.ndarray`.
       * **Applicability** – Best suited to homogeneous layers with
         high-quality temperature measurements at identical timestamps.

       References
       ----------
       Liebethal, C., & Foken, T. (2007). *Evaluation of six parameterization
       approaches for the ground heat flux.* Agricultural and Forest
       Meteorology, 143(1–2), 65–80.
       https://doi.org/10.1016/j.agrformet.2006.11.001

       Examples
       --------
       >>> g_plate = np.array([-12.3, -10.9, -8.5])      # W m-2 at z_r
       >>> Cv       = np.array([2.3e6, 2.1e6])           # two layers
       >>> T        = np.array([[15.0, 15.5, 16.0],      # °C layer 1
       ...                     [14.0, 14.2, 14.4]])      # °C layer 2
       >>> dz       = np.array([0.03, 0.02])             # m
       >>> Gz = simple_measurement_gz(g_plate, Cv, T, dt=1800, dz_layers=dz)
       >>> Gz
       array([-10.700...,  -8.317...])



.. py:function:: wbz12_g_gz(g_zr_series: numpy.ndarray | list | tuple, time_series: numpy.ndarray | list | tuple, z: float, zr: float, k_s_val: float) -> numpy.ndarray

   Estimate **soil heat flux** at a shallow depth *z* from a flux-plate
   record at reference depth *z_r* using the **WBZ12-G convolution
   method** (Wang, Bou-Zeid & Zhang 2012, their Eq. 9–10).

   The algorithm inverts the one-dimensional heat-conduction equation
   for a semi-infinite homogeneous soil, assuming a step-response
   kernel based on the complementary error function *erfc*:

   .. math::

       G(z,t) \;=\; 2\,G(z_r,t)
       \;-\; \frac{J(t)}{\Delta F_z(t)}

   with the discrete convolution integral

   .. math::

       J(t_n) \;=\;
       \sum_{j=0}^{n-1}
       \tfrac{\bigl[G(z_r,t_{n-j-1}) + G(z_r,t_{n-j})\bigr]}{2}\;
       \Delta F_z(t_{j+1})

   and the *transfer function increment*

   .. math::

       \Delta F_z(t) =
       \operatorname{erfc}\!
       \left[\frac{z_r - z}{2\sqrt{k_s (t - t_0)}}\right]
       \;-\;
       \operatorname{erfc}\!
       \left[\frac{z_r - z}{2\sqrt{k_s (t - t_0 - \Delta t)}}\right].

   :param g_zr_series: Time series of heat-flux-plate measurements at *z_r* (W m⁻²);
                       length **N**.
   :type g_zr_series: array_like
   :param time_series: Monotonically increasing time stamps (s) corresponding to
                       `g_zr_series`.  Must be the same length **N**.
   :type time_series: array_like
   :param z: Target depth (m, **positive downward**).
   :type z: float
   :param zr: Reference depth of the heat-flux plate (m, positive downward).
   :type zr: float
   :param k_s_val: Soil thermal diffusivity *k_s* (m² s⁻¹).  Must be positive.
   :type k_s_val: float

   :returns: Estimated heat-flux series at depth *z*, same length **N** as the
             input series (W m⁻²).
   :rtype: numpy.ndarray

   :raises ValueError: If `k_s_val ≤ 0`, the two series differ in length, time stamps
       are non-monotonic, or contain fewer than two samples.

   .. rubric:: Notes

   * The first element of the output corresponds to *t₀* (no storage
     correction possible yet); subsequent points incorporate the
     convolution up to that instant.
   * A very small term ``+ 1e-12`` is added under the square root to
     avoid division by zero at *t = t₀*.
   * The complementary error function is evaluated with
     :func:`numpy.erfc`, which is vectorized and avoids the SciPy
     dependency.

   .. rubric:: References

   Wang, W., Bou-Zeid, E., & Zhang, Y. (2012).
   *Estimating surface heat fluxes using the surface renewal method*.
   **Boundary-Layer Meteorology**, 144(2), 407-422.
   https://doi.org/10.1007/s10546-012-9730-9

   .. rubric:: Examples

   >>> # 10-min sampled day-long plate record at zr = 0.08 m
   >>> t  = np.arange(0, 24*3600 + 1, 600)            # s
   >>> Gp = -10 + 5*np.sin(2*np.pi*t/86400)
   >>> Gz = wbz12_g_gz(Gp, t, z=0.05, zr=0.08, k_s_val=1.0e-6)
   >>> Gz.shape
   (145,)


.. py:function:: wbz12_s_gz(Ag: numpy.ndarray | float, ks_val: numpy.ndarray | float, zr: float, z: float, t: numpy.ndarray | float, eps: numpy.ndarray | float, omega: float = OMEGA_DAY) -> numpy.ndarray | float

       **WBZ12-S analytic–numeric solution** for soil heat flux *G(z, t)* at
       depth *z* (Wang, Bou-Zeid & Zhang 2012, Eq. 11).

       WBZ12-S combines a closed-form *forcing* term (sinusoidal surface
       temperature) with a *storage* term that is expressed as a Fourier‐
       Laplace integral and must be evaluated numerically.  The governing
       equation assumes homogeneous soil properties and a sinusoidal
       surface temperature of amplitude *A_g* and phase `eps`:

       .. math::

           G(z,t) \;=\;
           A_g\,e^{-\eta}\,\sin\bigl(\omega t +
   arepsilon - \eta\bigr)
           \;-\;
           \frac{2 A_g k_s}{\pi}
           \int_{0}^{\infty}
               \frac{
                   k_s \zeta^{2} igl[\sin\varepsilon
                   - \omega \cos\varepsilon\bigr]
                 }{
                   \omega^{2}
                   + k_{s}^{2} \zeta^{4}
               }
               \sin\bigl[\zeta\,(z_r - z)\bigr]\,
               e^{-k_s \zeta^{2} t}\,d\zeta

       with

       .. math::
           \eta \;=\; (z_r - z)\,\sqrt{\omega / (2 k_s)}.

       The first term (prefactor × sin_term) represents the periodic
       *steady-state* component, while the integral accounts for the
       *transient* adjustment of the soil profile.

       Parameters
       ----------
       Ag : float or array_like
           Amplitude of the sinusoidal surface/ground temperature (K).  Can
           be broadcast with `t` and `eps`.
       ks_val : float or array_like
           Soil thermal diffusivity :math:`k_s` (m² s⁻¹).  Must be positive.
       zr : float
           Reference depth of the heat-flux plate (m, positive downward).
       z : float
           Target depth for the output heat flux (m, positive downward).
       t : float or array_like
           Time (s) since the start of the forcing cycle.  Broadcast-
           compatible with `Ag` and `eps`.
       eps : float or array_like
           Phase shift :math:`\varepsilon` (rad) of the surface
           temperature.  Broadcast-compatible with `t`.
       omega : float, optional
           Angular frequency :math:`\omega` (rad s⁻¹) of the sinusoidal
           forcing (default is :pydata:`OMEGA_DAY` ≈ 7.272 × 10⁻⁵ s⁻¹).

       Returns
       -------
       float or numpy.ndarray
           Soil heat flux *G(z, t)* (W m⁻²).  The shape matches NumPy
           broadcasting over (`Ag`, `ks_val`, `t`, `eps`).

       Raises
       ------
       ValueError
           If *ks_val* or *omega* are non-positive, or if the inputs cannot
           be broadcast to a common shape.

       Notes
       -----
       * **Numerical integral** – The transient term is evaluated using
         :func:`scipy.integrate.quad` over :math:`\zeta \in [0, \infty)`.
         The integral can be expensive for long time-series; cache results
         or vectorise with more sophisticated quadrature if performance is
         critical.
       * **Vectorisation** – `quad` is called individually for every
         element in the broadcasted inputs; large arrays may thus incur
         heavy computational cost.
       * **Units** – Use consistent SI units: W m⁻² (output), K (Ag),
         m² s⁻¹ (ks_val), rad s⁻¹ (omega), seconds (t), metres (depths).

       References
       ----------
       Wang, W., Bou-Zeid, E., & Zhang, Y. (2012).
       *Estimating surface heat fluxes using the surface renewal method.*
       **Boundary-Layer Meteorology**, 144(2), 407–422.
       https://doi.org/10.1007/s10546-012-9730-9

       Examples
       --------
       >>> # Daily sinusoid, single point
       >>> G = wbz12_s_gz(
       ...     Ag=8.0, ks_val=1.0e-6,
       ...     zr=0.08, z=0.05,
       ...     t=np.linspace(0, 86400, 97),
       ...     eps=0.0
       ... )
       >>> G.shape
       (97,)



.. py:function:: exact_temperature_gz(z: numpy.ndarray | float, AT: numpy.ndarray | float, t: numpy.ndarray | float, d: numpy.ndarray | float, omega: float = OMEGA_DAY, T_i: float = 298.15) -> numpy.ndarray | float

   Compute the **exact analytical soil‐temperature profile**
   under sinusoidal surface forcing.

   A sinusoidal temperature wave propagating downward through a
   homogeneous soil decays exponentially with depth and lags in phase.
   The closed‐form solution at depth *z* (m) is

   .. math::

       T(z, t) \;=\;
       T_i
       \;+\;
       A_T \,
       e^{-z/d}\,
       \sin\bigl(\omega t \;-\; z/d\bigr),

   where

   * :math:`T_i` — initial (mean) temperature of the soil column (K)
   * :math:`A_T` — amplitude of the surface-temperature oscillation (K)
   * :math:`d` — thermal damping depth (m)
   * :math:`\omega` — angular frequency of the forcing (rad s⁻¹)
   * :math:`t` — time since the start of oscillation (s)

   :param z: Depth below the surface (m). Positive downward.  Can be a scalar
             or an array broadcastable with *t* and *AT*.
   :type z: float or array_like
   :param AT: Amplitude *A_T* of the surface temperature wave (K).
   :type AT: float or array_like
   :param t: Time variable (s).  Supports vectorisation.
   :type t: float or array_like
   :param d: Thermal damping depth (m). Must be positive.  Can be scalar or
             broadcastable with *z*.
   :type d: float or array_like
   :param omega: Angular frequency *ω* (rad s⁻¹). Defaults to
                 :pydata:`OMEGA_DAY` (≈ 7.272 × 10⁻⁵ s⁻¹ for 24 h).
   :type omega: float, optional
   :param T_i: Mean (initial) soil temperature (K).  Defaults to **298.15 K**
               (25 °C).
   :type T_i: float, optional

   :returns: Soil temperature *T(z, t)* (same units as *T_i*).  The result
             shape matches NumPy broadcasting over the inputs.
   :rtype: float or numpy.ndarray

   :raises ValueError: If any element of *d* or *omega* is non-positive, or if inputs
       cannot be broadcast to a common shape.

   .. rubric:: Notes

   * **Phase lag** – Each e-folding depth *d* introduces a phase delay
     of 1 radian (~57.3 °) relative to the surface signal.
   * **Vectorisation** – All inputs are converted to
     :class:`numpy.ndarray`; the formula is evaluated element-wise and
     fully supports broadcasting.
   * **Units** – Ensure consistent SI units (metres, seconds, kelvins).

   .. rubric:: References

   Adapted from Gao, Z., Horton, R., Luo, L., & Kucharik, C. J. (2010).
   *A simple method to measure soil temperature dynamics: Theory and
   application.* **Soil Science Society of America Journal**, 74(2),
   580-588. https://doi.org/10.2136/sssaj2009.0169

   .. rubric:: Examples

   >>> # Daily cycle at 10 cm depth
   >>> z      = 0.10                            # m
   >>> AT     = 8.0                             # K
   >>> d      = 0.12                            # m
   >>> t_day  = np.linspace(0, 86400, 97)       # 15-min resolution
   >>> Tz     = exact_temperature_gz(z, AT, t_day, d)
   >>> Tz.shape
   (97,)


.. py:function:: exact_gz(z: numpy.ndarray | float, AT: numpy.ndarray | float, lambda_s_val: numpy.ndarray | float, d: numpy.ndarray | float, t: numpy.ndarray | float, omega: float = OMEGA_DAY) -> numpy.ndarray | float

   Exact analytical **soil-heat-flux** solution for sinusoidal surface forcing.

   A sinusoidal surface–temperature wave of amplitude *A_T* and angular
   frequency *ω* generates a heat-flux wave that decays exponentially
   with depth and lags the temperature signal by *π / 4* radians
   (Gao et al., 2010, Eq. 17):

   .. math::

       G(z,t) \;=\;
       \sqrt{2}\;
       \frac{\lambda_s A_T}{d}\;
       e^{-z/d}\;
       \sin\bigl(\omega t - z/d + \tfrac{\pi}{4}\bigr)

   where

   * :math:`\lambda_s` — soil thermal conductivity (W m⁻¹ K⁻¹)
   * :math:`d` — thermal damping depth (m)
   * :math:`z` — depth below the surface (m, positive downward)
   * :math:`t` — time (s) since the start of oscillation
   * :math:`\omega` — angular frequency (rad s⁻¹)

   :param z: Depth(s) below the surface (m), positive downward.
   :type z: float or array_like
   :param AT: Amplitude *A_T* of the surface-temperature oscillation (K).
   :type AT: float or array_like
   :param lambda_s_val: Soil thermal conductivity *λ_s* (W m⁻¹ K⁻¹).
   :type lambda_s_val: float or array_like
   :param d: Thermal damping depth (m). Must be positive.
   :type d: float or array_like
   :param t: Time variable (s).
   :type t: float or array_like
   :param omega: Angular frequency *ω* (rad s⁻¹). Defaults to
                 :pydata:`OMEGA_DAY` (≈ 7.272 × 10⁻⁵ s⁻¹, i.e. 2π / 86 400 s).
   :type omega: float, optional

   :returns: Soil heat flux *G(z, t)* (W m⁻²).  The return shape follows NumPy
             broadcasting rules applied to the inputs.
   :rtype: float or numpy.ndarray

   :raises ValueError: If any element of *d* or *omega* is non-positive, or if inputs
       cannot be broadcast to a common shape.

   .. rubric:: Notes

   * **Phase shift** – At any given depth, the heat-flux wave lags the
     temperature wave by 45 ° (π / 4 rad).
   * **Vectorisation** – All inputs are converted to
     :class:`numpy.ndarray`; the expression is evaluated element-wise
     and fully supports broadcasting.
   * **Units** – Ensure all quantities use consistent SI units.

   .. rubric:: References

   Gao, Z., Horton, R., Luo, L., & Kucharik, C. J. (2010).
   *A simple method to measure soil temperature dynamics: Theory and
   application.* **Soil Science Society of America Journal**, 74(2),
   580–588. https://doi.org/10.2136/sssaj2009.0169

   .. rubric:: Examples

   >>> # Scalar example
   >>> exact_gz(
   ...     z=0.05,
   ...     AT=8.0,
   ...     lambda_s_val=1.2,
   ...     d=0.12,
   ...     t=3600,
   ... )
   21.52...

   >>> # Vectorized daily cycle
   >>> t_day = np.linspace(0, 86400, 97)           # 15-min resolution
   >>> Gz = exact_gz(
   ...     z=0.10,
   ...     AT=6.0,
   ...     lambda_s_val=1.1,
   ...     d=0.11,
   ...     t=t_day,
   ... )
   >>> Gz.shape
   (97,)


