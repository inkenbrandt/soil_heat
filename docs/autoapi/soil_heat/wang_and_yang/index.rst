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

   This function calculates the conductive heat flux between soil layers
   based on the temperature gradient and thermal conductivity.

   .. math:: G = -\lambda_s \frac{\partial T}{\partial z}

   :param Tz: A 1D array of temperatures at the **centre** of each soil layer (K).
   :type Tz: ArrayLike
   :param dz: A 1D array of the thickness of each soil layer (m).
   :type dz: ArrayLike
   :param lambda_s: Thermal conductivity for each layer (W m⁻¹ K⁻¹). Can be a single
                    value (for homogeneous soil) or an array matching `Tz`.
   :type lambda_s: ArrayLike or float

   :returns: An array of heat fluxes (W m⁻²) at the *interfaces* between layers,
             including the surface and bottom boundaries. The length of the
             output is `len(Tz) + 1`. Positive values indicate downward flux.
   :rtype: numpy.ndarray


.. py:function:: integrated_soil_heat_flux(rho_c: ArrayLike, T_before: ArrayLike, T_after: ArrayLike, dz: ArrayLike, dt: float, G_ref: float = 0.0) -> numpy.ndarray

   Calculate the soil heat flux profile by integrating the change in
   heat storage upwards from a reference depth (Yang & Wang 2008, Eq. 5).

   .. math::
       G(z_i) = G(z_{ref}) + \int_{z_i}^{z_{ref}} \rho_s c_s
                \frac{\partial T}{\partial t} dz

   :param rho_c: A 1D array of volumetric heat capacity (`ρ_s c_s`) for each soil
                 layer (J m⁻³ K⁻¹).
   :type rho_c: ArrayLike
   :param T_before: A 1D array of temperatures at the start of the time step (K).
   :type T_before: ArrayLike
   :param T_after: A 1D array of temperatures at the end of the time step (K).
   :type T_after: ArrayLike
   :param dz: A 1D array of layer thicknesses (m).
   :type dz: ArrayLike
   :param dt: The duration of the time step (s).
   :type dt: float
   :param G_ref: The heat flux at the lower reference depth `z_ref` (W m⁻²). This
                 is typically assumed to be zero at a sufficient depth, by default 0.0.
   :type G_ref: float, optional

   :returns: An array of heat fluxes (W m⁻²) at the *upper interface* of each
             layer, with a size of `len(dz)`.
   :rtype: numpy.ndarray


.. py:function:: volumetric_heat_capacity(theta: ArrayLike, theta_sat: float | ArrayLike) -> numpy.ndarray

   Calculate the volumetric heat capacity of moist soil based on its
   water content (Yang & Wang 2008, Eq. 4).

   .. math::
       \rho_s c_s = (1 - \theta_{sat}) \rho_{d}c_{d} + \theta \rho_w c_w

   :param theta: Volumetric water content (m³ m⁻³).
   :type theta: ArrayLike
   :param theta_sat: Soil porosity, i.e., saturated volumetric water content (m³ m⁻³).
   :type theta_sat: float or ArrayLike

   :returns: The volumetric heat capacity `ρ_s c_s` (J m⁻³ K⁻¹).
   :rtype: numpy.ndarray


.. py:function:: stretched_grid(n: int, D: float, xi: float) -> numpy.ndarray

   Generate layer thicknesses for a vertically stretched grid.

   The grid layer thickness increases exponentially with depth, allowing
   for higher resolution near the surface.

   .. math:: \Delta z_i = \Delta z_0 \exp(\xi (i-1))

   :param n: The number of soil layers.
   :type n: int
   :param D: The total depth of the soil domain (m).
   :type D: float
   :param xi: The stretching parameter. `xi = 0` results in a uniform grid.
              `xi > 0` results in a grid that is finer at the top.
   :type xi: float

   :returns: A 1D array of thickness `Δz_i` for each of the `n` layers (m).
   :rtype: numpy.ndarray


.. py:function:: solve_tde(T_prev: ArrayLike, dz: ArrayLike, rho_c: ArrayLike, lambda_s: ArrayLike | float, Tsfc: float, Tbot: float, dt: float) -> numpy.ndarray

   Solve the 1D thermal diffusion equation for one time step using an
   implicit Crank-Nicolson scheme (Yang & Wang 2008, Eq. 7).

   This function sets up and solves the tridiagonal system of linear
   equations `M * T_new = D` to find the temperature profile at the
   next time step.

   :param T_prev: Temperature profile at the previous time step (K).
   :type T_prev: ArrayLike
   :param dz: Layer thicknesses (m).
   :type dz: ArrayLike
   :param rho_c: Volumetric heat capacity of each layer (J m⁻³ K⁻¹).
   :type rho_c: ArrayLike
   :param lambda_s: Thermal conductivity of each layer (W m⁻¹ K⁻¹).
   :type lambda_s: ArrayLike or float
   :param Tsfc: Surface temperature boundary condition (K).
   :type Tsfc: float
   :param Tbot: Bottom temperature boundary condition (K).
   :type Tbot: float
   :param dt: Time step (s).
   :type dt: float

   :returns: The new temperature profile at `t + dt`, including boundary nodes.
   :rtype: numpy.ndarray


.. py:function:: correct_profile(T_model: ArrayLike, depths_model: ArrayLike, T_obs: ArrayLike, depths_obs: ArrayLike) -> numpy.ndarray

   Correct a modeled temperature profile using observed temperatures.

   This function calculates the bias (error) between the model and
   observations at the observation depths, then linearly interpolates
   this bias across the entire model grid to correct the profile.

   :param T_model: The modeled temperature profile (K).
   :type T_model: ArrayLike
   :param depths_model: The depths corresponding to `T_model` (m).
   :type depths_model: ArrayLike
   :param T_obs: The observed temperatures (K).
   :type T_obs: ArrayLike
   :param depths_obs: The depths of the observations (m).
   :type depths_obs: ArrayLike

   :returns: The corrected temperature profile.
   :rtype: numpy.ndarray


.. py:function:: surface_temperature_longwave(R_lw_up: float, R_lw_dn: float, emissivity: float = 0.98) -> float

   Calculate surface temperature from upward and downward long-wave
   radiation measurements using the Stefan-Boltzmann law (Yang & Wang 2008, Eq. 8).

   .. math::
       T_s = \left[ \frac{R_{lw}^{\uparrow} - (1 - \epsilon) R_{lw}^{\downarrow}}
                     {\epsilon \sigma} \right]^{1/4}

   :param R_lw_up: Upwelling long-wave radiation (W m⁻²).
   :type R_lw_up: float
   :param R_lw_dn: Downwelling long-wave radiation (W m⁻²).
   :type R_lw_dn: float
   :param emissivity: Surface emissivity (dimensionless), by default 0.98.
   :type emissivity: float, optional

   :returns: The calculated surface temperature (K).
   :rtype: float


.. py:function:: thermal_conductivity_yang2008(theta: ArrayLike, theta_sat: float, rho_dry: float | ArrayLike) -> numpy.ndarray

   Estimate soil thermal conductivity based on soil moisture and dry
   density, following Yang et al. (2005) as cited in
   Yang & Wang (2008, Eq. 9).

   :param theta: Volumetric water content (m³ m⁻³).
   :type theta: ArrayLike
   :param theta_sat: Saturated volumetric water content (porosity) (m³ m⁻³).
   :type theta_sat: float
   :param rho_dry: Dry soil bulk density (kg m⁻³).
   :type rho_dry: float or ArrayLike

   :returns: The estimated soil thermal conductivity `λ_s` (W m⁻¹ K⁻¹).
   :rtype: numpy.ndarray


.. py:function:: flux_error_linear(rho_c: ArrayLike, S2_minus_S1: ArrayLike, dt: float) -> numpy.ndarray

   Calculate the diagnostic error in heat flux that arises from assuming
   a linear temperature profile between measurement points
   (Yang & Wang 2008, Eq. 11).

   .. math:: \Delta G_i = \frac{\rho_s c_s (S_{i,2} - S_{i,1})}{\Delta t}

   :param rho_c: Volumetric heat capacity for each layer (J m⁻³ K⁻¹).
   :type rho_c: ArrayLike
   :param S2_minus_S1: The area difference representing the deviation from linearity.
   :type S2_minus_S1: ArrayLike
   :param dt: The time step (s).
   :type dt: float

   :returns: The flux error for each layer (W m⁻²).
   :rtype: numpy.ndarray


.. py:function:: surface_energy_residual(R_net: float, H: float, LE: float, G0: float) -> float

   Calculate the residual of the surface energy budget (Yang & Wang 2008, Eq. 12).

   .. math:: \Delta E = R_{net} - (H + LE + G_0)

   :param R_net: Net radiation (W m⁻²).
   :type R_net: float
   :param H: Sensible heat flux (W m⁻²).
   :type H: float
   :param LE: Latent heat flux (W m⁻²).
   :type LE: float
   :param G0: Surface ground heat flux (W m⁻²).
   :type G0: float

   :returns: The energy balance residual `ΔE` (W m⁻²).
   :rtype: float


.. py:function:: tdec_step(T_prev: ArrayLike, dz: ArrayLike, theta: ArrayLike, theta_sat: float, rho_dry: float, lambda_const: float, Tsfc: float, Tbot: float, dt: float, depths_model: ArrayLike, T_obs: ArrayLike, depths_obs: ArrayLike) -> Tuple[numpy.ndarray, numpy.ndarray]

   Perform a single integration step of the TDEC (Temperature Diffusion
   Error Correction) scheme.

   This function encapsulates the predict-correct sequence for one time step.

   :param T_prev: Temperature profile from the previous time step (K).
   :type T_prev: ArrayLike
   :param dz: Layer thicknesses (m).
   :type dz: ArrayLike
   :param theta: Volumetric water content profile (m³ m⁻³).
   :type theta: ArrayLike
   :param theta_sat: Soil porosity (m³ m⁻³).
   :type theta_sat: float
   :param rho_dry: Dry soil bulk density (kg m⁻³).
   :type rho_dry: float
   :param lambda_const: A constant thermal conductivity for the prediction step (W m⁻¹ K⁻¹).
   :type lambda_const: float
   :param Tsfc: Surface temperature boundary condition (K).
   :type Tsfc: float
   :param Tbot: Bottom temperature boundary condition (K).
   :type Tbot: float
   :param dt: Time step (s).
   :type dt: float
   :param depths_model: Depths of the model grid nodes (m).
   :type depths_model: ArrayLike
   :param T_obs: Observed temperatures for the correction step (K).
   :type T_obs: ArrayLike
   :param depths_obs: Depths of the observed temperatures (m).
   :type depths_obs: ArrayLike

   :returns: A tuple containing:
             - `T_corr`: The corrected temperature profile at `t + dt`.
             - `G_prof`: The corresponding heat flux profile (W m⁻²) at layer interfaces.
   :rtype: Tuple[np.ndarray, np.ndarray]


