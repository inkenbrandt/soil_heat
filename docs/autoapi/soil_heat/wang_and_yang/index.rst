soil_heat.wang_and_yang
=======================

.. py:module:: soil_heat.wang_and_yang

.. autoapi-nested-parse::

   Yang & Wang (2008) – Equation Set Implementation
   =================================================
   Python translation of every numbered equation in:

   > Yang, K., & Wang, J. (2008). *A temperature prediction‑correction method for
   > estimating surface soil heat flux from soil temperature and moisture data.*
   > *Science in China Series D: Earth Sciences, 51*(5), 721‑729.
   > https://doi.org/10.1007/s11430‑008‑0036‑1

   The paper introduces a *Temperature‑Diffusion plus Error‑Correction* (TDEC)
   approach for estimating surface soil heat flux.  This module provides a direct
   one‑to‑one mapping from each equation in the paper (Eqs. 1–12) to a Python
   function.  Helper utilities required by those equations—matrix creators, grid
   stretching, tridiagonal solvers, etc.—are also included.

   All functions accept **NumPy arrays** or scalars where applicable and are fully
   type‑annotated.  Docstrings use the **NumPy docstring standard** so they can be
   rendered by *Sphinx‑napoleon*.



Functions
---------

.. autoapisummary::

   soil_heat.wang_and_yang.soil_heat_flux
   soil_heat.wang_and_yang.integrated_soil_heat_flux
   soil_heat.wang_and_yang.volumetric_heat_capacity
   soil_heat.wang_and_yang.stretched_grid
   soil_heat.wang_and_yang.solve_tde
   soil_heat.wang_and_yang.correct_profile
   soil_heat.wang_and_yang.surface_temperature_longwave
   soil_heat.wang_and_yang.thermal_conductivity_yang2008
   soil_heat.wang_and_yang.flux_error_linear
   soil_heat.wang_and_yang.surface_energy_residual
   soil_heat.wang_and_yang.tdec_step


Module Contents
---------------

.. py:function:: soil_heat_flux(Tz: ArrayLike, dz: ArrayLike, lambda_s: ArrayLike | float) -> numpy.ndarray

   Compute heat flux *G* at cell interfaces using Fourier’s law.

   :param Tz: Temperature at layer **centres** (K).
   :type Tz: ArrayLike
   :param dz: Thickness of each layer (m).
   :type dz: ArrayLike
   :param lambda_s: Thermal conductivity for each layer (W m‑1 K‑1).
   :type lambda_s: ArrayLike or float

   :returns: Heat flux at *interfaces* (positive downward) with shape ``len(Tz)+1``.
   :rtype: np.ndarray


.. py:function:: integrated_soil_heat_flux(rho_c: ArrayLike, T_before: ArrayLike, T_after: ArrayLike, dz: ArrayLike, dt: float, G_ref: float = 0.0) -> numpy.ndarray

   Discrete integration of Eq. (3)/(5) to obtain heat‑flux profile.

   :param rho_c: Volumetric heat capacity ``ρ_s c_s`` for each layer (J m‑3 K‑1).
   :type rho_c: ArrayLike
   :param T_before: Temperatures at two successive timesteps (K).
   :type T_before: ArrayLike
   :param T_after: Temperatures at two successive timesteps (K).
   :type T_after: ArrayLike
   :param dz: Layer thicknesses (m).
   :type dz: ArrayLike
   :param dt: Timestep (s).
   :type dt: float
   :param G_ref: Heat flux at the lower reference depth *z_ref* (W m‑2).  Often ≈ 0.
   :type G_ref: float, default 0

   :returns: Heat flux at the *upper* interface of every layer, size ``len(dz)``.
   :rtype: np.ndarray


.. py:function:: volumetric_heat_capacity(theta: ArrayLike, theta_sat: float | ArrayLike) -> numpy.ndarray

   Volumetric heat capacity of moist soil.

   :param theta: Volumetric water content (m³ m‑3).
   :type theta: ArrayLike
   :param theta_sat: Soil porosity.
   :type theta_sat: float or ArrayLike

   :returns: ``ρ_s c_s`` (J m‑3 K‑1).
   :rtype: np.ndarray


.. py:function:: stretched_grid(n: int, D: float, xi: float) -> numpy.ndarray

   Generate *n* layer thicknesses following the exponential stretching rule.

   :param n: Number of layers.
   :type n: int
   :param D: Total domain depth (m).
   :type D: float
   :param xi: Stretching parameter; 0 → uniform grid.
   :type xi: float

   :returns: Thickness ``Δz_i`` for each layer *i* (m).
   :rtype: np.ndarray


.. py:function:: solve_tde(T_prev: ArrayLike, dz: ArrayLike, rho_c: ArrayLike, lambda_s: ArrayLike | float, Tsfc: float, Tbot: float, dt: float) -> numpy.ndarray

   Implicit Crank‑Nicholson (θ = 1) solve of Eq. (7).

   Boundary conditions (Eq. 7a, 7c) are Dirichlet.


.. py:function:: correct_profile(T_model: ArrayLike, depths_model: ArrayLike, T_obs: ArrayLike, depths_obs: ArrayLike) -> numpy.ndarray

   Add linear‑interpolated bias (Eq. ΔT_k) to model profile.


.. py:function:: surface_temperature_longwave(R_lw_up: float, R_lw_dn: float, emissivity: float = 0.98) -> float

   Convert upward/downward long‑wave radiation to surface temperature (Eq. 8).


.. py:function:: thermal_conductivity_yang2008(theta: ArrayLike, theta_sat: float, rho_dry: float | ArrayLike) -> numpy.ndarray

   Estimate soil thermal conductivity following Yang et al. (2005) (Eq. 9).


.. py:function:: flux_error_linear(rho_c: ArrayLike, S2_minus_S1: ArrayLike, dt: float) -> numpy.ndarray

   Error introduced when using a LINEAR temperature profile (diagnostic).


.. py:function:: surface_energy_residual(R_net: float, H: float, LE: float, G0: float) -> float

   Return the residual *ΔE* in Eq. (12).


.. py:function:: tdec_step(T_prev: ArrayLike, dz: ArrayLike, theta: ArrayLike, theta_sat: float, rho_dry: float, lambda_const: float, Tsfc: float, Tbot: float, dt: float, depths_model: ArrayLike, T_obs: ArrayLike, depths_obs: ArrayLike) -> Tuple[numpy.ndarray, numpy.ndarray]

   One integration step of the TDEC scheme.

   :returns: * **T_corr** (*np.ndarray*) -- Corrected temperature profile at *t + dt*.
             * **G_prof** (*np.ndarray*) -- Heat‑flux profile (W m‑2) at layer interfaces.


