Soil Thermal Conductivity and Heat-Flux Estimation
==================================================

.. contents::
   :depth: 2
   :local:

Introduction
------------

This document explains how to estimate soil thermal diffusivity
(:math:`\alpha`), volumetric heat capacity (:math:`C_v`),
thermal conductivity (:math:`k`), and soil heat flux (:math:`G`)
from a network of temperature–moisture sensors installed at nine depths
(5 cm → 90 cm, sampled every 30 min).

Core concepts
~~~~~~~~~~~~~

* **Amplitude damping** – Daily temperature‐wave amplitude declines with depth.
* **Phase (time-lag)** – The peak temperature at deeper layers lags behind the
  surface.
* **Fourier’s law** – Governs vertical conductive heat flow:

  .. math::

     G = -\,k \,\frac{\partial T}{\partial z}

  *:math:`G` positive downward.*

* **Link between parameters** –

  .. math::

     k = \alpha \, C_v

Measurement-based estimation
----------------------------

Daily angular frequency
~~~~~~~~~~~~~~~~~~~~~~~

.. math::

   \omega = \frac{2\pi}{P}

where :math:`P = 86\,400\ \text{s}` (24 h).

Thermal diffusivity from amplitude
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. math::

   \alpha_\mathrm{amp}
   = \frac{\omega\,(z_2 - z_1)^2}
          {2\,\ln\!\bigl(A_1 / A_2\bigr)}

Thermal diffusivity from phase lag
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. math::

   \alpha_\mathrm{lag}
   = \frac{\omega\,(z_2 - z_1)^2}
          {2\,\Delta t^{\,2}}

Volumetric heat capacity
~~~~~~~~~~~~~~~~~~~~~~~~

.. math::

   C_v = (1 - \theta_v)\,C_\text{soil} + \theta_v\,C_\text{water}

Typical constants:

* :math:`C_\text{soil} = 1.9\times10^{6}\ \text{J m}^{-3}\,{}^\circ\text{C}^{-1}`
* :math:`C_\text{water} = 4.18\times10^{6}\ \text{J m}^{-3}\,{}^\circ\text{C}^{-1}`

Python implementation
---------------------

.. code-block:: python

   import numpy as np
   import pandas as pd

   SEC_PER_DAY = 86_400
   OMEGA = 2 * np.pi / SEC_PER_DAY


   # ────────────────────────────────
   # 1.  Thermal-diffusivity helpers
   # ────────────────────────────────
   def thermal_diffusivity_amplitude(A1, A2, z1, z2, period=SEC_PER_DAY):
       """
       Diffusivity from amplitude damping (m² s⁻¹).
       """
       omega = 2 * np.pi / period
       return omega * (z2 - z1) ** 2 / (2 * np.log(A1 / A2))


   def thermal_diffusivity_lag(delta_t, z1, z2, period=SEC_PER_DAY):
       """
       Diffusivity from phase lag (m² s⁻¹).
       """
       omega = 2 * np.pi / period
       return omega * (z2 - z1) ** 2 / (2 * delta_t ** 2)


   # ───────────────────────────────
   # 2.  Volumetric heat-capacity Cv
   # ───────────────────────────────
   def volumetric_heat_capacity(theta_v,
                                c_soil=1.9e6,
                                c_water=4.18e6):
       """
       Cv from volumetric water content θ_v (fraction).
       """
       return (1 - theta_v) * c_soil + theta_v * c_water


   # ───────────────────────────────
   # 3.  Thermal conductivity   k
   # ───────────────────────────────
   def thermal_conductivity(alpha, theta_v):
       """
       k (W m⁻¹ °C⁻¹) from α and θ_v.
       """
       Cv = volumetric_heat_capacity(theta_v)
       return alpha * Cv


   # ───────────────────────────────
   # 4.  Temperature gradient  ∂T/∂z
   # ───────────────────────────────
   def temperature_gradient(T_upper, T_lower, z1, z2):
       """
       ∂T/∂z (°C m⁻¹) between two depths.
       """
       return (T_lower - T_upper) / (z2 - z1)


   # ───────────────────────────────
   # 5.  Soil heat flux          G
   # ───────────────────────────────
   def soil_heat_flux(T_upper, T_lower, z1, z2, k):
       """
       G (W m⁻²).  Positive downward.
       """
       grad = temperature_gradient(T_upper, T_lower, z1, z2)
       return -k * grad


Example workflow
----------------

.. code-block:: python

   # Example inputs
   A1, A2 = 6.0, 1.5          # °C
   z1, z2 = 0.05, 0.30        # m
   dt_lag = 3 * 3600          # s
   theta_v = 0.25             # 25 % moisture
   T_5cm, T_30cm = 18.0, 17.0 # °C snapshots

   alpha1 = thermal_diffusivity_amplitude(A1, A2, z1, z2)
   alpha2 = thermal_diffusivity_lag(dt_lag, z1, z2)
   alpha = (alpha1 + alpha2) / 2               # m² s⁻¹

   k_est = thermal_conductivity(alpha, theta_v)  # W m⁻¹ °C⁻¹
   G = soil_heat_flux(T_5cm, T_30cm, z1, z2, k_est)

   print(f"α = {alpha:.2e} m² s⁻¹")
   print(f"k = {k_est:.2f} W m⁻¹ °C⁻¹")
   print(f"G = {G:.1f} W m⁻²  (positive = downward)")

Typical conductivity ranges
---------------------------

+---------------------+--------------------------------------+
| Soil condition      | Approx. :math:`k` (W m⁻¹ °C⁻¹)       |
+=====================+======================================+
| Dry sandy soil      | 0.25 – 0.40                          |
+---------------------+--------------------------------------+
| Moist sandy soil    | 0.80 – 1.50                          |
+---------------------+--------------------------------------+
| Moist loamy soil    | 0.60 – 1.30                          |
+---------------------+--------------------------------------+
| Moist clay soil     | 0.80 – 1.60                          |
+---------------------+--------------------------------------+


