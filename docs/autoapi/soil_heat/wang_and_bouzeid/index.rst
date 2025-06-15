soil_heat.wang_and_bouzeid
==========================

.. py:module:: soil_heat.wang_and_bouzeid

.. autoapi-nested-parse::

   soil_ground_heat_flux.py
   ====================================================
   Python utilities implementing the equations from:

       Wang, Z.-H., & Bou-Zeid, E. (2012). *A novel approach for the estimation of
       soil ground heat flux*. *Agricultural and Forest Meteorology*, 154-155,
       214-221.

   All equations appearing in that paper are reproduced as vectorised Python
   functions, together with a few convenience helpers.  The library is fully
   NumPy-aware and can be used on scalars or on time-series arrays.

   Units
   -----
   All functions assume SI units **throughout**:

   * depth ``z``            ― m  (positive downward)
   * time ``t``             ― s
   * temperature ``T``      ― K  (or °C provided the 0-offset is handled
     consistently)
   * heat flux ``G, H, LE`` ― W m-2
   * radiative flux ``Rn``  ― W m-2
   * soil properties
       * thermal conductivity ``k`` ― W m-1 K-1
       * heat capacity ``rho_c``    ― J m-3 K-1
       * thermal diffusivity ``kappa = k / rho_c`` ― m² s-1

   Dependencies
   ------------
   >>> pip install numpy scipy

   .. rubric:: Example

   >>> import numpy as np, wang_and_bouzeid as sghf
   >>> # 30-minute series (dt = 1800 s) of flux-plate measurements at z = 0.08 m
   >>> Gz = np.loadtxt('Gz_8cm.txt')
   >>> G0 = sghf.estimate_G0_from_Gz(Gz, z_r=0.08, kappa=0.7e-6, dt=1800)



Attributes
----------

.. autoapisummary::

   soil_heat.wang_and_bouzeid.surface_energy_residual


Functions
---------

.. autoapisummary::

   soil_heat.wang_and_bouzeid.energy_balance_residual
   soil_heat.wang_and_bouzeid.ground_heat_flux_conventional
   soil_heat.wang_and_bouzeid.green_function_temperature
   soil_heat.wang_and_bouzeid.temperature_convolution_solution
   soil_heat.wang_and_bouzeid.soil_heat_flux_from_G0
   soil_heat.wang_and_bouzeid.estimate_G0_from_Gz
   soil_heat.wang_and_bouzeid.sinusoidal_boundary_flux
   soil_heat.wang_and_bouzeid.soil_temperature_sinusoidal
   soil_heat.wang_and_bouzeid.soil_heat_flux_sinusoidal
   soil_heat.wang_and_bouzeid.heat_capacity_moist_soil
   soil_heat.wang_and_bouzeid.pf_from_theta
   soil_heat.wang_and_bouzeid.thermal_conductivity_moist_soil
   soil_heat.wang_and_bouzeid.thermal_diffusivity


Module Contents
---------------

.. py:function:: energy_balance_residual(Rn: float | numpy.ndarray, H: float | numpy.ndarray, LE: float | numpy.ndarray, G0: float | numpy.ndarray) -> float | numpy.ndarray

       Compute the **closure residual** of the surface energy balance (SEB).

       The classical SEB for land–atmosphere exchange is

       .. math::

           R_n - G_0 \;=\; H + LE +
   arepsilon ,

       where the residual term :math:`
   arepsilon` quantifies the lack of
       closure.  Rearranging gives

       .. math::


   arepsilon \;=\; R_n - G_0 - H - LE ,

       which is what this helper returns.

       Parameters
       ----------
       Rn : float or array_like
           Net radiation *Rₙ* (W m⁻²).  Positive downward.
       H : float or array_like
           Sensible heat flux *H* (W m⁻²).  Positive upward (atmosphere ← surface).
       LE : float or array_like
           Latent heat flux *LE* (W m⁻²).  Positive upward.
       G0 : float or array_like
           Ground (soil) heat flux *G₀* (W m⁻²).  Positive downward
           (into the soil).  Some authors use the opposite sign convention;
           ensure consistency with *Rn*.

       Returns
       -------
       float or ndarray
           Energy‐balance residual :math:`
   arepsilon` (W m⁻²) with the
           broadcast shape of the inputs.
           **Positive** values indicate missing energy
           (surface gains > turbulent + ground fluxes),
           whereas **negative** values mean an apparent energy surplus.

       Notes
       -----
       * **Broadcasting** – All inputs are treated with NumPy broadcasting,
         allowing scalars, 1-D arrays, or DataFrame columns.
       * **Closure diagnostics** – The residual can be summarised as a mean
         bias or expressed as a relative closure fraction
         ``1 − ε / (Rn − G0)``.
       * **Typical magnitudes** – Eddy‐covariance towers often report
         10–30 % non-closure.  Persistently large |ε| values may indicate
         advective transport, sensor tilt, or data‐processing issues.

       Examples
       --------
       >>> ε = energy_balance_residual(
       ...         Rn=df["Rn"],  H=df["H"],  LE=df["LE"],  G0=df["G"]
       ...     )
       >>> ε.describe(percentiles=[0.05, 0.5, 0.95])
       count    17280.000000
       mean        -9.73
       std         34.51
       5%         -58.42
       50%         -8.11
       95%         39.12
       Name: residual, dtype: float64

       Plot daily average residuals:

       >>> ε.resample("D").mean().plot(marker="o")
       >>> plt.axhline(0, color="k", lw=0.8)
       >>> plt.ylabel("Energy balance residual (W m$^{-2}$)")
       >>> plt.title("Daily SEB closure")



.. py:data:: surface_energy_residual

.. py:function:: ground_heat_flux_conventional(k: float, dT_dz_at_zr: float, rho_c: float, dT_dt_profile: Sequence[float], z_profile: Sequence[float]) -> float

   Conventional estimate of *surface* ground heat flux, Eq. (2).

   :param k: Effective soil thermal conductivity **(W m-1 K-1)**.
   :param dT_dz_at_zr: Vertical temperature gradient evaluated at the flux-plate depth ``z_r``
                       **(K m-1)**.  A *negative* gradient means temperature decreases with
                       depth.
   :param rho_c: Volumetric heat capacity of the soil **(J m-3 K-1)**.
   :param dT_dt_profile: Time derivatives ∂T/∂t for *each* node between the surface and ``z_r``
                         **(K s-1)**.  Any iterable (list, ndarray, …).  Must align with
                         ``z_profile``.
   :param z_profile: Depth of each node in the temperature profile **(m)**.  Increasing,
                     positive downward, **excluding** the surface (z = 0) but *including*
                     ``z_r`` (last element).

   :returns: **G0** -- Ground heat flux at the surface **(W m-2)**.  Positive *into* the soil.
   :rtype: float


.. py:function:: green_function_temperature(z: float, t: float, kappa: float) -> float

   Green-function solution :math:`g_z(t)` for the **one‐dimensional,
   semi-infinite heat equation** with a unit surface‐flux impulse at
   :math:`t=0,\,z=0`.

   The governing partial differential equation is

   .. math::

       \frac{\partial T}{\partial t}
       \;=\;
       \kappa\,\frac{\partial^{2} T}{\partial z^{2}},
       \qquad z \ge 0,\; t > 0,

   with initial condition :math:`T(z,0)=0` and boundary condition
   corresponding to a Dirac δ–heat pulse applied at the surface
   (:math:`z=0`).  The resulting Green function (Carslaw & Jaeger,
   1959, Eq. 7) is

   .. math::

       g_z(t)
       \;=\;
       \frac{2}{\sqrt{\pi}}
       \,\sqrt{\kappa t}\;
       \exp\!\Bigl[-\frac{z^{2}}{4\kappa t}\Bigr]
       \;-\;
       z\,\operatorname{erfc}\!
       \Bigl[\frac{z}{2\sqrt{\kappa t}}\Bigr],
       \qquad t>0,

   and :math:`g_z(t)=0` for :math:`t\le 0` (causality).

   :param z: Depth below the surface (m, positive downward).  Must be
             non-negative.
   :type z: float
   :param t: Time since the surface impulse (s).  Values :math:`t \le 0`
             return 0 by definition.
   :type t: float
   :param kappa: Thermal diffusivity :math:`\kappa` of the half-space
                 (m² s⁻¹).
   :type kappa: float

   :returns: Green-function value :math:`g_z(t)` (units **m**, because the
             solution integrates heat-flux density with respect to depth to
             yield temperature).
   :rtype: float

   .. rubric:: Notes

   * **Causality check** – If ``t`` ≤ 0 the function short-circuits and
     returns 0.0.
   * **Vectorisation** – For vector or array input use
     :func:`numpy.vectorize` or wrap the function in
     :pyfunc:`numpy.frompyfunc`; the core implementation is scalar for
     numerical clarity.
   * **Usage** – ``g_z(t)`` can be convolved with an arbitrary surface
     heat-flux time series ``q₀(t)`` to obtain temperature at depth
     via

     .. math::

        T(z,t) \;=\; \int_{0}^{t} g_z(t-τ)\,q_0(τ)\;\mathrm dτ .

   .. rubric:: References

   Carslaw, H. S., & Jaeger, J. C. (1959).
   *Conduction of Heat in Solids* (2nd ed., p. 100).
   Oxford University Press.

   .. rubric:: Examples

   >>> g = green_function_temperature(z=0.05, t=3_600, kappa=1.4e-7)
   >>> print(f"g(5 cm, 1 h) = {g:.3e} m")
   g(5 cm, 1 h) = 7.42e-04 m


.. py:function:: temperature_convolution_solution(z: float, t_series: numpy.ndarray, f_series: numpy.ndarray, kappa: float, Ti: float = 0.0) -> numpy.ndarray

   Temperature time-series at depth *z* via Duhamel convolution (Eq. 6).

   ``T(z,t) = Ti + ∫ f(t-τ) d g_z(τ)``

   The integral becomes a discrete convolution where *f* is the boundary
   heat-flux series (W m-2  → ∂T/∂z via Fourier).


.. py:function:: soil_heat_flux_from_G0(z: float, t_series: numpy.ndarray, G0_series: numpy.ndarray, kappa: float) -> numpy.ndarray

   Compute *G(z,t)* from a known surface flux series *G0* (Eq. 9).


.. py:function:: estimate_G0_from_Gz(Gz_series: numpy.ndarray, z_r: float, kappa: float, dt: float) -> numpy.ndarray

   Estimate *surface* ground heat flux *G0* from plate measurements *Gz*.

   Implements discretised Eq. (11) – the recursion proposed by Wang & Bou-Zeid
   (2012).  Time-series must be *regularly* sampled.

   :param Gz_series: Soil heat-flux measurements at depth *z_r* **(W m-2)**.
   :type Gz_series: np.ndarray
   :param z_r: Plate depth **(m)**.
   :type z_r: float
   :param kappa: Thermal diffusivity **(m² s-1)**.
   :type kappa: float
   :param dt: Sampling interval **(s)**.
   :type dt: float

   :returns: **G0** -- Estimated surface heat-flux series **(W m-2)**.
   :rtype: np.ndarray


.. py:function:: sinusoidal_boundary_flux(t: float | numpy.ndarray, A: float, omega: float, epsilon: float) -> float | numpy.ndarray

       Evaluate a **sinusoidal surface heat-flux forcing**

       .. math::

           q_0(t) \;=\; A \,\sin\!igl(\omega t +
   arepsilonigr),

       which is commonly used as a boundary condition for analytical soil-
       heat-flux solutions (see Eq. 13 of the companion text).

       Parameters
       ----------
       t : float or array_like
           Time since the start of the simulation (s).
           May be scalar or any NumPy-broadcastable shape; units must be
           consistent with ``omega``.
       A : float
           Amplitude of the surface heat flux (W m⁻²).  Positive **downward**
           into the soil.
       omega : float
           Angular frequency (rad s⁻¹).  For a diurnal cycle
           ``omega = 2 π / 86 400`` ≈ 7.272 × 10⁻⁵ rad s⁻¹.
       epsilon : float
           Phase shift **ε** (rad).  Positive values delay the flux peak,
           negative values advance it.

       Returns
       -------
       ndarray or float
           Instantaneous surface heat flux *q₀(t)* (W m⁻²) with shape given
           by NumPy broadcasting of the inputs.

       Notes
       -----
       * **Sign convention** — Positive *q₀* adds energy to the soil column; be
         sure it matches the sign convention of your governing equation.
       * **Vectorisation** — The implementation is a single call to
         ``numpy.sin`` and therefore fully vectorised.
       * **Period** — The period *P* (s) of the forcing is related to
         ``omega`` by *P = 2 π / ω*.

       Examples
       --------
       >>> import numpy as np, matplotlib.pyplot as plt
       >>> t = np.linspace(0, 2*86400, 1_000)                 # 2 days
       >>> q0 = sinusoidal_boundary_flux(
       ...         t, A=120, omega=2*np.pi/86400, epsilon=0)
       >>> plt.plot(t/3600, q0)
       >>> plt.xlabel("Time (h)")
       >>> plt.ylabel("Surface heat flux $q_0$ (W m$^{-2}$)")
       >>> plt.title("Sinusoidal surface forcing (A = 120 W m⁻²)")
       >>> plt.show()



.. py:function:: soil_temperature_sinusoidal(z: float, t: float | numpy.ndarray, A: float, omega: float, epsilon: float, Ti: float, kappa: float) -> float | numpy.ndarray

       Analytical solution for **soil temperature** beneath a sinusoidally
       forced surface heat‐flux boundary.

       The temperature response of a semi-infinite, homogeneous soil column
       to a surface heat flux

       .. math::

           q_0(t) \;=\; A \sin(\omega t +
   arepsilon)

       is (Carslaw & Jaeger 1959, Eq. 14)

       .. math::

           T(z,t)
           \;=\;
           T_i
           \;+\;
           \underbrace{\frac{A}{\kappa\sqrt{\omega}}
                        e^{-z r}\,
                        \sin\!igl(\omega t +
   arepsilon - z r - π/4\bigr)}
           _{\text{steady harmonic}}
           \;+\;
           \underbrace{T_{\text{trans}}(z,t)}
           _{\text{transient integral}} ,

       where :math:`r = \sqrt{\omega / 2\kappa}`.
       The first term is the *steady* periodic component that
       propagates downward with exponentially damped amplitude and a
       depth-dependent phase lag.  The second term accounts for the
       *transient* adjustment from the initial uniform temperature *Tᵢ* and
       is evaluated numerically here by vectorised trapezoidal quadrature.

       Parameters
       ----------
       z : float
           Depth below the soil surface (m, positive downward).
       t : float or array_like
           Time since the start of the forcing (s).  Scalar or NumPy array.
       A : float
           Amplitude of the sinusoidal **surface heat flux** (W m⁻², positive
           downward).  Consistent with
           :pyfunc:`sinusoidal_boundary_flux`.
       omega : float
           Angular frequency of the forcing (rad s⁻¹).
           For a diurnal wave ``omega = 2 * π / 86_400``.
       epsilon : float
           Phase shift ε (rad) of the surface heat-flux wave.
       Ti : float
           Initial uniform soil temperature *Tᵢ* (°C or K).
       kappa : float
           Thermal diffusivity κ of the soil (m² s⁻¹).

       Returns
       -------
       float or ndarray
           Soil temperature *T(z, t)* in the same units as *Ti*.
           If *t* is an array the returned array has the same shape.

       Notes
       -----
       * **Steady component** –
         The exponential term ``exp(-z*r)`` dampens amplitude with depth
         while the phase lag is ``z*r + π/4``.
       * **Transient component** –
         The integral is truncated at ``x = 50 / z`` (or 50 if *z = 0*) and
         evaluated with 2 000 panels, which empirically yields < 0.1 %
         relative error for typical soil parameters.  Adjust the limits or
         panel count for higher precision.
       * **Vectorisation** –
         For array *t* the quadrature is performed in parallel using
         broadcasting; memory usage scales with ``len(t)`` × 2 000.
       * **Units** –
         Ensure *A* is W m⁻² **heat flux**.  If you have a surface
         temperature wave instead, transform to an equivalent heat-flux
         boundary or modify the formulation.

       Examples
       --------
       >>> import numpy as np, matplotlib.pyplot as plt
       >>> z   = 0.05                      # 5 cm
       >>> k   = 1.4e-7                    # m² s⁻¹
       >>> A   = 120                       # W m⁻²
       >>> ω   = 2*np.pi/86400             # diurnal
       >>> ε   = 0                         # no phase shift
       >>> Ti  = 15.0                      # °C
       >>> t   = np.linspace(0, 172800, 2000)   # 2 days
       >>> Tz  = soil_temperature_sinusoidal(z, t, A, ω, ε, Ti, k)
       >>> plt.plot(t/3600, Tz)
       >>> plt.xlabel("Time (h)")
       >>> plt.ylabel("Temperature (°C)")
       >>> plt.title("Analytical diurnal soil temperature at 5 cm")
       >>> plt.show()



.. py:function:: soil_heat_flux_sinusoidal(z: float, t: float | numpy.ndarray, A: float, omega: float, epsilon: float, kappa: float) -> float | numpy.ndarray

       Analytical **soil heat–flux** response *G(z,t)* to a sinusoidal
       surface–flux boundary condition.

       A semi-infinite, homogeneous soil column forced at the surface by

       .. math::

           q_0(t) \;=\; A \,\sin\!igl(\omega t +
   arepsilonigr)

       admits the Green-function solution (Carslaw & Jaeger, 1959, Eq. 15)

       .. math::

           G(z,t)
           \;=\;
           A\,e^{-z r}\,\sin(\omega t +
   arepsilon - z r)
           \;+\;
           G_{     ext{trans}}(z,t),

       where :math:`r = \sqrt{ \omega / 2\kappa }` and the transient term

       .. math::

           G_{     ext{trans}}(z,t)
           \;=\;
           -
   rac{2 A \kappa}{\pi}
           \int_{0}^{\infty}
           \frac{( \kappa x^{2}\sin
   arepsilon - \omega \cos
   arepsilon )
                 \;x \sin(x z)}
                {\omega^{2} + \kappa^{2} x^{4}}
           \,e^{-\kappa x^{2} t}\; \mathrm d x .

       The first term is the steady, exponentially damped harmonic that
       propagates downward with a depth-dependent phase lag; the second
       term describes the transient adjustment from an initially unheated
       half-space.  The integral is evaluated here by vectorised
       trapezoidal quadrature on a finite domain (*x* ≤ 50 / *z*).

       Parameters
       ----------
       z : float
           Depth below the soil surface (m, positive downward).
       t : float or array_like
           Time since the onset of forcing (s).  Scalar or NumPy array.
       A : float
           Amplitude of the *surface* heat flux (W m⁻², positive downward).
       omega : float
           Angular frequency ω of the forcing (rad s⁻¹).
           For a diurnal wave ``omega = 2 π / 86 400``.
       epsilon : float
           Phase shift ε of the surface flux (radians).
       kappa : float
           Thermal diffusivity κ of the soil (m² s⁻¹).

       Returns
       -------
       float or ndarray
           Heat flux *G(z,t)* at depth *z* (W m⁻²) with shape equal to
           ``t`` after NumPy broadcasting.
           Positive values **downward**, matching the sign convention of *A*.

       Notes
       -----
       * **Steady component** –
         Amplitude decays as ``exp(-z*r)``; phase lag increases linearly
         with depth by *z r*.
       * **Transient component** –
         The integration limit ``50 / z`` (or 50 for *z = 0*) yields < 0.1 %
         relative error for common parameter ranges.  Increase the upper
         bound and/or panel count (2 000) for stricter accuracy.
       * **Vectorisation** –
         For array *t* the quadrature is evaluated for all time steps
         simultaneously via broadcasting; memory usage scales with
         ``len(t) × 2000``.
       * **Coupling with temperature** –
         The companion function
         :pyfunc:`soil_temperature_sinusoidal` gives *T(z,t)* under the
         same boundary forcing; *G* and *T* satisfy Fourier’s law
         ``G = -k ∂T/∂z`` once thermal conductivity *k* is specified.

       References
       ----------
       Carslaw, H. S., & Jaeger, J. C. (1959).
       *Conduction of Heat in Solids* (2nd ed., pp. 100–102).
       Oxford University Press.

       Examples
       --------
       >>> import numpy as np, matplotlib.pyplot as plt
       >>> z      = 0.05                    # 5 cm
       >>> kappa  = 1.4e-7                  # m² s⁻¹
       >>> A      = 120                     # W m⁻² downward
       >>> omega  = 2*np.pi/86400           # diurnal frequency
       >>> eps    = 0
       >>> t      = np.linspace(0, 172800, 2000)   # 2 days
       >>> G      = soil_heat_flux_sinusoidal(z, t, A, omega, eps, kappa)
       >>> plt.plot(t/3600, G)
       >>> plt.xlabel("Time (h)")
       >>> plt.ylabel("Heat flux $G$ (W m$^{-2}$)")
       >>> plt.title("Analytical diurnal soil heat flux at 5 cm")
       >>> plt.show()



.. py:function:: heat_capacity_moist_soil(theta_v: float | numpy.ndarray, theta_s: float, rho_c_s: float = 1260000.0, rho_c_w: float = 4200000.0) -> float | numpy.ndarray

   Volumetric heat capacity of moist soil, Eq. (16).

   :param theta_v: Volumetric water content **(m³ m-3)**.
   :param theta_s: Porosity (saturated volumetric water content) **(m³ m-3)**.
   :param rho_c_s: Heat capacity of dry soil / water **(J m-3 K-1)**.
   :param rho_c_w: Heat capacity of dry soil / water **(J m-3 K-1)**.


.. py:function:: pf_from_theta(theta_v: float | numpy.ndarray, theta_s: float, psi_s: float, b: float) -> float | numpy.ndarray

   Return Pf (Eq. 18) from volumetric water content.


.. py:function:: thermal_conductivity_moist_soil(theta_v: float | numpy.ndarray, theta_s: float, psi_s: float, b: float) -> float | numpy.ndarray

   Thermal conductivity parameterisation, Eq. (17).


.. py:function:: thermal_diffusivity(k: float | numpy.ndarray, rho_c: float | numpy.ndarray) -> float | numpy.ndarray

   Return κ = k / (ρ c).


