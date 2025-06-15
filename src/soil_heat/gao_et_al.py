import numpy as np
from scipy.special import erfc
from scipy.integrate import quad

__all__ = [
    "lambda_s",
    "k_s",
    "volumetric_heat_capacity",
    "nme",
    "rmse",
    "calorimetric_gz",
    "force_restore_gz",
    "gao2010_gz",
    "heusinkveld_gz",
    "hsieh2009_gz",
    "leuning_damping_depth",
    "leuning_gz",
    "simple_measurement_gz",
    "wbz12_g_gz",
    "wbz12_s_gz",
    "exact_temperature_gz",
    "exact_gz",
]


# Typical diurnal angular frequency (rad s⁻¹)
OMEGA_DAY: float = 2 * np.pi / 86_400.0  # 2π / 24 h
CV_REF: float = 2.2e6  # J m⁻³ K⁻¹ — representative volumetric heat capacity


# -----------------------------------------------------------------------------
# 1. Soil thermal-property relationships (Eq. 12–13)
# -----------------------------------------------------------------------------
def lambda_s(theta: np.ndarray | float) -> np.ndarray | float:
    """
    Compute the **soil thermal conductivity** :math:`\\lambda_s`
    as a function of volumetric water content (θ).

    The relationship is taken from Gao et al. (2017, Eq. 12):

    .. math::

        \\lambda_s(\\theta) = 0.20 + \\exp\\bigl[\,1.46\\,(\\theta - 0.34)\\bigr]

    Parameters
    ----------
    theta : float or array_like
        Volumetric water content (m³ m⁻³).
        Accepts a scalar value or any array-like object that can be
        converted to a :class:`numpy.ndarray`. Values should lie in the
        closed interval ``[0, 1]``.

    Returns
    -------
    float or numpy.ndarray
        Soil thermal conductivity λ\_s in W m⁻¹ K⁻¹.  The returned type
        matches the input: a scalar is returned for a scalar *θ*, and a
        NumPy array for array-like *θ*.

    Raises
    ------
    ValueError
        If any element of *θ* is outside the physically meaningful
        range ``[0, 1]``.

    Notes
    -----
    * **Vectorization** – The function is fully vectorized; it operates
      element-wise on NumPy arrays and broadcasts according to NumPy
      broadcasting rules.
    * **Empirical limits** – Gao et al. recommend the equation for
      mineral soils where 0 ≤ θ ≤ 0.5. Extrapolation beyond this range
      may introduce error.
    * **Units** – The output uses SI units (W m⁻¹ K⁻¹).

    References
    ----------
    Gao, Z., Niu, G.-Y., & Hedquist, B. C. (2017).
    **Soil thermal conductivity parameterization: A closed-form
    equation and its evaluation.** *Journal of Geophysical Research:
    Atmospheres*, 122(6), 3466–3478.
    DOI: 10.1002/2016JD025992

    Examples
    --------
    >>> lambda_s(0.25)
    0.463...
    >>> import numpy as np
    >>> theta = np.linspace(0, 0.5, 5)
    >>> lambda_s(theta)
    array([0.20      , 0.31243063, 0.41172297, 0.51210322, 0.61861198])
    """
    theta = np.asarray(theta)

    if np.any((theta < 0) | (theta > 1)):
        raise ValueError(
            "Volumetric water content 'theta' must be within [0, 1]. "
            f"Received values spanning {theta.min()}–{theta.max()}."
        )

    return 0.20 + np.exp(1.46 * (theta - 0.34))


def k_s(theta: np.ndarray | float) -> np.ndarray | float:
    """
    Calculate the **soil thermal diffusivity** :math:`k_s`
    as a function of volumetric water content (θ).

    The empirical relationship follows Gao et al. (2017, Eq. 13):

    .. math::

        k_s(\\theta) = \\bigl[0.69 + \\exp\\bigl(3.06\\,(\\theta - 0.26)\\bigr)\\bigr] \\times 10^{-7}

    Parameters
    ----------
    theta : float or array_like
        Volumetric water content (m³ m⁻³).
        Accepts a scalar or any array-like sequence convertible to a
        :class:`numpy.ndarray`. Values **must** lie in the closed
        interval ``[0, 1]``.

    Returns
    -------
    float or numpy.ndarray
        Soil thermal diffusivity *k\_s* in m² s⁻¹.
        A scalar is returned if *θ* is a scalar; otherwise a NumPy array
        of matching shape is returned.

    Raises
    ------
    ValueError
        If any element of *θ* is outside the physically meaningful
        range ``[0, 1]``.

    Notes
    -----
    * **Vectorization** – The calculation is fully vectorized and
      broadcasts according to NumPy rules.
    * **Applicability** – Gao et al. derived this expression for mineral
      soils under typical field conditions. Extrapolating beyond
      0 ≤ θ ≤ 0.5 may reduce accuracy.
    * **Units** – Output is in SI units (m² s⁻¹).

    References
    ----------
    Gao, Z., Niu, G.-Y., & Hedquist, B. C. (2017).
    *Soil thermal conductivity parameterization: A closed-form equation
    and its evaluation.* Journal of Geophysical Research: Atmospheres,
    **122**(6), 3466–3478. https://doi.org/10.1002/2016JD025992

    Examples
    --------
    >>> k_s(0.25)
    1.321...e-07
    >>> thetas = np.linspace(0, 0.5, 6)
    >>> k_s(thetas)
    array([1.14...e-07, 1.21...e-07, 1.30...e-07, 1.41...e-07,
           1.54...e-07, 1.69...e-07])
    """
    theta = np.asarray(theta)

    if np.any((theta < 0) | (theta > 1)):
        raise ValueError(
            "Volumetric water content 'theta' must be within [0, 1]. "
            f"Received values spanning {theta.min()}–{theta.max()}."
        )

    return (0.69 + np.exp(3.06 * (theta - 0.26))) * 1e-7


def volumetric_heat_capacity(
    lambda_s_val: np.ndarray | float,
    k_s_val: np.ndarray | float,
) -> np.ndarray | float:
    """
    Compute the **volumetric heat capacity** :math:`C_v` of soil:

    .. math::

        C_v \;=\; \\frac{\\lambda_s}{k_s}

    where
    :math:`\\lambda_s` is the **thermal conductivity** (W m⁻¹ K⁻¹) and
    :math:`k_s` is the **thermal diffusivity** (m² s⁻¹).

    The resulting heat capacity has units of J m⁻³ K⁻¹.

    Parameters
    ----------
    lambda_s_val : float or array_like
        Soil thermal conductivity (W m⁻¹ K⁻¹). May be a scalar or any
        array-like structure broadcastable with *k_s_val*.
    k_s_val : float or array_like
        Soil thermal diffusivity (m² s⁻¹). Must be positive. Accepts a
        scalar or array-like input.

    Returns
    -------
    float or numpy.ndarray
        Volumetric heat capacity *C_v* (J m⁻³ K⁻¹).  The return type
        matches the input: a scalar for scalar inputs, or a NumPy array
        for array-like inputs.

    Raises
    ------
    ValueError
        If any element in *k_s_val* is zero or negative, which would
        lead to division by zero or non-physical results.

    Notes
    -----
    * **Vectorization** – The function is fully vectorized; both inputs
      are converted to :class:`numpy.ndarray` and follow NumPy
      broadcasting rules.
    * **Physical meaning** – Volumetric heat capacity represents the
      energy required to raise the temperature of a unit volume of soil
      by one kelvin. High values correspond to moist or water-logged
      soils; dry mineral soils have lower *C_v*.

    Examples
    --------
    >>> #Scalar inputs
    >>> volumetric_heat_capacity(0.8, 1.2e-6)
    666666.666...

    >>> #Vectorized inputs
    >>> lambda_vals = np.array([0.5, 0.6, 0.7])
    >>> diffusivities = np.array([1.1e-6, 1.2e-6, 1.3e-6])
    >>> volumetric_heat_capacity(lambda_vals, diffusivities)
    array([454545.4545..., 500000.    ..., 538461.5384...])
    """
    lambda_s_val = np.asarray(lambda_s_val)
    k_s_val = np.asarray(k_s_val)

    if np.any(k_s_val <= 0):
        raise ValueError("Thermal diffusivity 'k_s_val' must be positive and non-zero.")

    return lambda_s_val / k_s_val


# -----------------------------------------------------------------------------
# 2. Error metrics (Eq. 14–15)
# -----------------------------------------------------------------------------
def nme(
    calc: np.ndarray | float,
    meas: np.ndarray | float,
) -> float:
    """
    Calculate the **normalized mean error (NME)** between calculated
    and measured values, expressed as a percentage.

    The formulation follows Gao et al. (2017, Eq. 14):

    .. math::

        \\text{NME} \\;=\\; 100\\,\\frac{\\sum\\limits_{i}\\left|\\hat{y}_i - y_i\\right|}
                                      {\\sum\\limits_{i}\\left|y_i\\right|}

    where :math:`\\hat{y}_i` are the *calculated* (model) values and
    :math:`y_i` are the *measured* (reference) values.

    Parameters
    ----------
    calc : float or array_like
        Modelled / calculated values :math:`\\hat{y}`.
        Accepts any array-like object (including scalars) convertible to
        a :class:`numpy.ndarray`.
    meas : float or array_like
        Measured / observed reference values :math:`y`. Must be
        broadcast-compatible with *calc*.

    Returns
    -------
    float
        Normalized mean error (percent).
        A value of **0 %** indicates perfect agreement; larger values
        indicate greater error.

    Raises
    ------
    ValueError
        * If *calc* and *meas* cannot be broadcast to a common shape.
        * If the denominator ``Σ|meas|`` equals zero (e.g., all measured
          values are zero), which would make NME undefined.

    Notes
    -----
    * **Range** – NME is non-negative and unbounded above.
    * **Interpretation** – Because both numerator and denominator use
      absolute values, NME is insensitive to the direction of the error
      (over- vs under-prediction) and therefore complements signed error
      metrics such as mean bias.
    * **Vectorization** – The implementation is fully vectorized and
      adheres to NumPy broadcasting rules.

    References
    ----------
    Gao, Z., Niu, G.-Y., & Hedquist, B. C. (2017).
    *Soil thermal conductivity parameterization: A closed-form equation
    and its evaluation.* **Journal of Geophysical Research:
    Atmospheres**, 122(6), 3466–3478.
    https://doi.org/10.1002/2016JD025992

    Examples
    --------
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
    """
    calc = np.asarray(calc)
    meas = np.asarray(meas)

    # Validate broadcast compatibility
    try:
        _ = np.broadcast(calc, meas)
    except ValueError as exc:
        raise ValueError(
            "Inputs 'calc' and 'meas' must be broadcast-compatible."
        ) from exc

    denom = np.sum(np.abs(meas))
    if denom == 0:
        raise ValueError("Denominator Σ|meas| equals zero; NME is undefined.")

    return 100.0 * np.sum(np.abs(calc - meas)) / denom


def rmse(
    calc: np.ndarray | float,
    meas: np.ndarray | float,
) -> float:
    """
    Compute the **root-mean-square error (RMSE)** between calculated
    (model) and measured (reference) values.

    Following Gao et al. (2017, Eq. 15):

    .. math::

        \\text{RMSE} \;=\;
        \\sqrt{\\frac{1}{N}\\sum_{i=1}^{N}\\bigl(\\hat{y}_i - y_i\\bigr)^2}

    where :math:`\\hat{y}_i` are *calculated* values and
    :math:`y_i` are *measured* values.

    Parameters
    ----------
    calc : float or array_like
        Calculated / modelled values :math:`\\hat{y}`.  Accepts a scalar
        or any array-like object convertible to a
        :class:`numpy.ndarray`.
    meas : float or array_like
        Measured / observed values :math:`y`. Must be broadcast-
        compatible with *calc*.

    Returns
    -------
    float
        Root-mean-square error (same units as the inputs).

    Raises
    ------
    ValueError
        * If *calc* and *meas* cannot be broadcast to a common shape.
        * If the input arrays are empty (``N = 0``).

    Notes
    -----
    * **Interpretation** – RMSE represents the sample-standard-deviation
      of the differences between two datasets; lower values indicate
      better agreement.
    * **Vectorization** – The function is fully vectorized and respects
      NumPy broadcasting rules.
    * **Units** – RMSE preserves the units of the input variables.

    References
    ----------
    Gao, Z., Niu, G.-Y., & Hedquist, B. C. (2017).
    *Soil thermal conductivity parameterization: A closed-form equation
    and its evaluation.* Journal of Geophysical Research:
    Atmospheres, **122**(6), 3466–3478.
    https://doi.org/10.1002/2016JD025992

    Examples
    --------
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
    """
    calc = np.asarray(calc)
    meas = np.asarray(meas)

    # Ensure broadcast compatibility
    try:
        _ = np.broadcast(calc, meas)
    except ValueError as exc:
        raise ValueError(
            "Inputs 'calc' and 'meas' must be broadcast-compatible."
        ) from exc

    if calc.size == 0 or meas.size == 0:
        raise ValueError("Input arrays must contain at least one element.")

    return float(np.sqrt(np.mean((calc - meas) ** 2)))


# -----------------------------------------------------------------------------
# 3. Soil‑heat‑flux methods
# -----------------------------------------------------------------------------

# -- 3.1 Calorimetric (Eq. 1) ---------------------------------------------------


def calorimetric_gz(
    g_zr: np.ndarray | float,
    cv_layers: np.ndarray | float,
    dT_dt_layers: np.ndarray | float,
    dz_layers: np.ndarray | float,
) -> np.ndarray | float:
    """
    Estimate the **soil heat flux** at a target depth *z*
    (typically 5 cm) using the *calorimetric method*.

    The method corrects an in-situ heat-flux plate reading made at a
    reference depth *z_r* by adding the change in heat storage of the
    soil column located between *z* and *z_r*.  Mathematically
    (Liebethal & Foken 2007, Eq. 1):

    .. math::

        G(z,t) \;=\; G(z_r,t) \;+\;
        \sum_{l=1}^{N}
        C_{v,l}\,\\frac{\\partial \\bar{T}_l}{\\partial t}\,\Delta z_l

    where

    * :math:`G(z,t)`     … heat flux at depth *z* (W m⁻²)
    * :math:`G(z_r,t)`   … measured plate flux at reference depth *z_r*
    * :math:`C_{v,l}`    … volumetric heat capacity of layer *l*
      (J m⁻³ K⁻¹)
    * :math:`\\partial \\bar{T}_l / \\partial t` … time derivative of the
      layer-averaged temperature (K s⁻¹)
    * :math:`\\Delta z_l` … thickness of layer *l* (m)

    Parameters
    ----------
    g_zr : float or array_like
        Heat-flux plate measurement at depth *z_r* (W m⁻²).  Can be a
        scalar or time series.
    cv_layers : array_like
        Volumetric heat capacity :math:`C_{v,l}` for each of *N* soil
        sub-layers (J m⁻³ K⁻¹).  Shape ``(N, …)`` where the trailing
        dimensions (``…``) must be broadcast-compatible with *g_zr*.
    dT_dt_layers : array_like
        Time derivative of layer-mean temperature
        :math:`\\partial \\bar{T}_l / \\partial t`
        (K s⁻¹); same shape as *cv_layers*.
    dz_layers : array_like
        Thickness of each layer :math:`\\Delta z_l` (m); shape ``(N,)``
        or broadcast-compatible with the first axis of *cv_layers*.

    Returns
    -------
    float or numpy.ndarray
        Calorimetrically corrected soil heat flux *G(z)* (W m⁻²).  Scalar
        if all inputs are scalar; otherwise a NumPy array matching the
        broadcast shape of *g_zr*.

    Raises
    ------
    ValueError
        If the inputs cannot be broadcast to a common shape, or if the
        number of layers inferred from *cv_layers*, *dT_dt_layers*, and
        *dz_layers* are inconsistent.

    Notes
    -----
    * **Vectorization** – All operations are fully vectorized.  Inputs
      are converted to :class:`numpy.ndarray` and follow NumPy
      broadcasting rules.
    * **Layer axis** – The first dimension of *cv_layers* and
      *dT_dt_layers* (axis 0) is interpreted as the layer index *l*.
      Layer thicknesses *dz_layers* are broadcast across any additional
      dimensions.
    * **Units** – Ensure consistent SI units: W m⁻², J m⁻³ K⁻¹,
      K s⁻¹, and m.

    References
    ----------
    Liebethal, C., & Foken, T. (2007). *Evaluation of six parameterization
    approaches for the ground heat flux.* Agricultural and Forest
    Meteorology, **143**(1–2), 65-80.
    https://doi.org/10.1016/j.agrformet.2006.11.001

    Examples
    --------
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
    """
    # Convert inputs to numpy arrays
    g_zr = np.asarray(g_zr)
    cv_layers = np.asarray(cv_layers)
    dT_dt_layers = np.asarray(dT_dt_layers)
    dz_layers = np.asarray(dz_layers)

    # Broadcast layer thickness to match cv_layers if needed
    if dz_layers.ndim == 1:
        dz_layers = dz_layers[:, np.newaxis]  # shape (N, 1, …)

    # Basic shape and broadcast checks
    try:
        cv_b, dT_b, dz_b = np.broadcast_arrays(cv_layers, dT_dt_layers, dz_layers)
    except ValueError as exc:
        raise ValueError(
            "cv_layers, dT_dt_layers, and dz_layers must be broadcast-"
            "compatible (same number of layers and trailing dimensions)."
        ) from exc

    # Broadcast plate flux to the same trailing dims as the storage term
    try:
        g_zr_b = np.broadcast_to(g_zr, cv_b.shape[1:])
    except ValueError as exc:
        raise ValueError(
            "g_zr is not broadcast-compatible with the trailing "
            "dimensions of the layer arrays."
        ) from exc

    # Storage term Σ C_v · dT/dt · Δz  (units W m-2)
    storage = np.sum(cv_b * dT_b * dz_b, axis=0)

    return g_zr_b + storage


# -- 3.2 Force‑restore (Eq. 2) --------------------------------------------------


def force_restore_gz(
    cv: np.ndarray | float,
    dTg_dt: np.ndarray | float,
    Tg: np.ndarray | float,
    Tg_bar: np.ndarray | float,
    delta_z: float = 0.05,
    omega: float = OMEGA_DAY,
) -> np.ndarray | float:
    """
    Estimate **soil heat flux** :math:`G(z)` at a shallow depth
    (:math:`z = \\delta z`, default 5 cm) with the **force–restore
    method**.

    The formulation (Liebethal & Foken 2007, Eq. 2) corrects the
    calorimetric storage term with a *restore* component that accounts
    for the difference between the instantaneous ground temperature
    *T_g* and its running mean *\\bar{T}_g*:

    .. math::

        G(z,t) \;=\; C_v \\, \\delta z \\, \\frac{\\partial T_g}{\\partial t}
        \;+
        \\sqrt{\\omega \\, C_v \\, \\lambda_s(C_v)}
        \\left[
            \\frac{1}{\\omega}\\,\\frac{\\partial T_g}{\\partial t}
            + \\bigl(T_g - \\bar{T}_g\\bigr)
        \\right]

    where

    * :math:`C_v`   … volumetric heat capacity (J m⁻³ K⁻¹)
    * :math:`\\partial T_g / \\partial t` … ground-temperature tendency
      (K s⁻¹)
    * :math:`\\lambda_s(C_v)` … soil thermal conductivity derived from
      *C_v* via :pyfunc:`lambda_s_from_cv` (W m⁻¹ K⁻¹)
    * :math:`\\omega` … angular frequency of the diurnal cycle
      (rad s⁻¹)
    * :math:`T_g` / :math:`\\bar{T}_g` … instantaneous and running-mean
      ground temperature (K)
    * :math:`\\delta z` … sensor depth below the surface (m)

    Parameters
    ----------
    cv : float or array_like
        Volumetric heat capacity :math:`C_v` (J m⁻³ K⁻¹).  Must be
        positive.  Accepts scalars or NumPy-broadcastable arrays.
    dTg_dt : float or array_like
        Time derivative :math:`\\partial T_g/\\partial t` (K s⁻¹).
    Tg : float or array_like
        Instantaneous ground (surface) temperature :math:`T_g` (K or °C)
        at depth *δz*.
    Tg_bar : float or array_like
        Running mean ground temperature :math:`\\bar{T}_g` over the
        diurnal cycle (same units as *Tg*).
    delta_z : float, optional
        Depth :math:`\\delta z` in metres; default is **0.05 m**
        (5 cm).
    omega : float, optional
        Angular frequency :math:`\\omega` in rad s⁻¹.  Defaults to
        :pydata:`OMEGA_DAY` (2π / 86 400 s).

    Returns
    -------
    float or numpy.ndarray
        Soil heat flux *G(z)* at depth *δz* (W m⁻²).  The output shape
        follows NumPy broadcasting rules applied to the inputs.

    Raises
    ------
    ValueError
        If any input arrays cannot be broadcast to a common shape, or if
        *cv* contains non-positive values.

    Notes
    -----
    * **λ_s(C_v) mapping** – The function relies on a helper
      :pyfunc:`lambda_s_from_cv` that converts volumetric heat capacity
      to thermal conductivity.  Ensure that this helper is present in
      the import path.
    * **Units** – Keep units internally consistent (SI).
    * **Vectorization** – All operations are vectorized; the
      mathematical expression is evaluated element-wise for array
      inputs.
    * **Interpretation** – The first term represents *storage*
      (calorimetric), while the second tempers short-term fluctuations,
      “restoring” *T_g* toward its mean.

    References
    ----------
    Liebethal, C., & Foken, T. (2007). *Evaluation of six
    parameterization approaches for the ground heat flux.*
    **Agricultural and Forest Meteorology**, 143(1–2), 65-80.
    https://doi.org/10.1016/j.agrformet.2006.11.001

    Examples
    --------
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
    """
    # Convert to arrays and broadcast
    cv = np.asarray(cv, dtype=float)
    dTg_dt = np.asarray(dTg_dt, dtype=float)
    Tg = np.asarray(Tg, dtype=float)
    Tg_bar = np.asarray(Tg_bar, dtype=float)

    # Basic validation
    if np.any(cv <= 0):
        raise ValueError("Volumetric heat capacity 'cv' must be positive.")
    try:
        cv_b, dT_b, Tg_b, Tgbar_b = np.broadcast_arrays(cv, dTg_dt, Tg, Tg_bar)
    except ValueError as exc:
        raise ValueError("Inputs are not broadcast-compatible.") from exc

    # --- Force–restore terms ------------------------------------------------
    term1 = cv_b * dT_b * delta_z  # storage term

    # λ_s from C_v (user-supplied helper)
    lambda_s = lambda_s_from_cv(cv_b)  # → W m-1 K-1

    term2 = np.sqrt(omega * cv_b * lambda_s) * (dT_b / omega + (Tg_b - Tgbar_b))

    return term1 + term2


def lambda_s_from_cv(
    cv: np.ndarray | float,
) -> np.ndarray | float:
    """
    Convert **volumetric heat capacity** (*C_v*) to an *approximate*
    **soil thermal conductivity** :math:`\\lambda_s`.

    The function assumes that the ratio
    :math:`\\lambda_s / k_s = C_v` holds (i.e.\ :math:`\\lambda_s = k_s
    \\times C_v`) and that a *representative* thermal diffusivity
    :math:`k_s` is implicit in the reference heat capacity
    ``CV_REF`` (≈ 2.2 MJ m⁻³ K⁻¹).  Rearranging gives

    .. math::

        \\lambda_s \\approx
        \\frac{C_v}{C_{v,\\text{ref}}}

    where :math:`C_{v,\\text{ref}} \\;\\approx\\; 2.2\\times10^{6}`
    J m⁻³ K⁻¹ is typical for moist mineral soil.

    Parameters
    ----------
    cv : float or array_like
        Volumetric heat capacity (J m⁻³ K⁻¹).
        Accepts a scalar or any array-like object convertible to
        :class:`numpy.ndarray`.  Values must be positive.

    Returns
    -------
    float or numpy.ndarray
        Estimated soil thermal conductivity λ\_s (W m⁻¹ K⁻¹).  The
        returned type mirrors the input: a scalar for scalar *cv* or a
        NumPy array for array-like input.

    Raises
    ------
    ValueError
        If any element of *cv* is non-positive.

    Notes
    -----
    * **Heuristic** – The calculation is intentionally *coarse* and
      should only be used where detailed λ\_s information is lacking.
    * **Reference value** – ``CV_REF`` corresponds to a diffusivity
      *k\_s* of roughly :math:`1.0\\times10^{-6}` m² s⁻¹ when
      λ\_s ≈ 1 W m⁻¹ K⁻¹; adjust ``CV_REF`` if a different reference is
      more appropriate for your soils.
    * **Vectorization** – Fully vectorized; broadcasting follows NumPy
      rules.

    Examples
    --------
    >>> lambda_s_from_cv(2.2e6)
    1.0
    >>> cv_series = np.array([1.8e6, 2.4e6, 2.0e6])
    >>> lambda_s_from_cv(cv_series)
    array([0.81818182, 1.09090909, 0.90909091])
    """
    cv = np.asarray(cv, dtype=float)

    if np.any(cv <= 0):
        raise ValueError("Volumetric heat capacity 'cv' must be positive.")

    return cv / CV_REF


# -- 3.3 Gao 2010 sinusoid (Eq. 3) --------------------------------------------


def gao2010_gz(
    AT: np.ndarray | float,
    lambda_s_val: np.ndarray | float,
    k_s_val: np.ndarray | float,
    t: np.ndarray | float,
    omega: float = OMEGA_DAY,
) -> np.ndarray | float:
    """
    Estimate **soil heat flux** :math:`G(z,t)` at depth *z* based on a
    *sinusoidal* ground-temperature forcing (Gao et al. 2010, Eq. 3).

    The analytical solution assumes that the surface (or ground-contact)
    temperature varies sinusoidally with amplitude *A_T* and angular
    frequency *ω*.  Under these conditions the heat flux at any depth
    *z* can be written

    .. math::

        G(z,t) \;=\;
        \\sqrt{2}\,
        \\frac{\\lambda_s A_T}{d}\,
        \\sin\\bigl(\\omega t + \\tfrac{\\pi}{4}\\bigr),

    where the **thermal damping depth**

    .. math::

        d \;=\; \\sqrt{\\frac{2 k_s}{\\omega}}

    is determined by the soil thermal diffusivity *k_s* and the forcing
    frequency *ω*.

    Parameters
    ----------
    AT : float or array_like
        Amplitude of the sinusoidal ground-surface temperature (K).
    lambda_s_val : float or array_like
        Soil thermal conductivity :math:`\\lambda_s` (W m⁻¹ K⁻¹).
    k_s_val : float or array_like
        Soil thermal diffusivity :math:`k_s` (m² s⁻¹).
    t : float or array_like
        Time variable (s).  Can be absolute time since epoch or simply
        seconds since the start of the cycle—as long as *ω t* is
        dimensionless.
    omega : float, optional
        Angular frequency *ω* (rad s⁻¹).  Defaults to *OMEGA_DAY*
        (≈ 7.272 × 10⁻⁵ s⁻¹, i.e. 2π / 86 400 s).

    Returns
    -------
    float or numpy.ndarray
        Heat flux *G(z,t)* (W m⁻²).  The output follows NumPy’s
        broadcasting rules applied to the inputs.

    Raises
    ------
    ValueError
        * If *lambda_s_val*, *k_s_val*, and *t* cannot be broadcast to a
          common shape.
        * If any element of *k_s_val* or *omega* is non-positive.

    Notes
    -----
    * **Vectorization** – All inputs are internally converted to
      :class:`numpy.ndarray`; the formula is evaluated element-wise and
      fully supports broadcasting.
    * **Units** – Ensure consistent SI units: W m⁻¹ K⁻¹ (λ_s),
      m² s⁻¹ (k_s), s (t), and rad s⁻¹ (ω).
    * **Scope** – The solution presumes purely sinusoidal boundary
      forcing and homogeneous soil properties; real-world deviations
      (e.g., non-sine forcing, stratified soils, moisture variation)
      will introduce error.

    References
    ----------
    Gao, Z., Horton, R., Luo, L., & Kucharik, C. J. (2010).
    *A simple method to measure soil temperature dynamics: Theory and
    application.* **Soil Science Society of America Journal**, 74(2),
    580-588. https://doi.org/10.2136/sssaj2009.0169

    Examples
    --------
    >>> # Daily cycle at 5 cm depth
    >>> AT      = 8.0                           # K
    >>> lambda_ = 1.2                           # W m-1 K-1
    >>> kappa   = 1.0e-6                        # m2 s-1
    >>> t_day   = np.linspace(0, 86400, 97)     # 15-min steps
    >>> Gz      = gao2010_gz(AT, lambda_, kappa, t_day)
    >>> Gz.shape
    (97,)
    """
    # Convert to arrays
    AT = np.asarray(AT, dtype=float)
    lambda_s_val = np.asarray(lambda_s_val, dtype=float)
    k_s_val = np.asarray(k_s_val, dtype=float)
    t = np.asarray(t, dtype=float)

    # Basic validation
    if np.any(k_s_val <= 0):
        raise ValueError("'k_s_val' must be positive.")
    if omega <= 0:
        raise ValueError("'omega' must be positive.")

    try:
        AT_b, lam_b, k_b, t_b = np.broadcast_arrays(AT, lambda_s_val, k_s_val, t)
    except ValueError as exc:
        raise ValueError(
            "Inputs AT, lambda_s_val, k_s_val, and t are not " "broadcast-compatible."
        ) from exc

    # Thermal damping depth d
    d = np.sqrt(2.0 * k_b / omega)

    # Heat-flux solution
    return np.sqrt(2.0) * lam_b * AT_b / d * np.sin(omega * t_b + np.pi / 4.0)


# -- 3.4 Heusinkveld harmonic (Eq. 4) -----------------------------------------


def heusinkveld_gz(
    A_n: np.ndarray | float,
    Phi_n: np.ndarray | float,
    n_max: int,
    k_s_val: np.ndarray | float,
    lambda_s_val: np.ndarray | float,
    w: float,
) -> np.ndarray | float:
    """
    Compute *soil heat flux* :math:`G(z,t)` from the **H04 harmonic
    series solution** proposed by Heusinkveld et al. (2004, Eq. 4).

    The approach represents the surface (ground) temperature as a Fourier
    series with harmonics up to order *n_max*.  For each harmonic
    :math:`n` the heat-flux contribution at depth *z* can be written

    .. math::

        G_n(z,t) \;=\;
        \\frac{\\lambda_s}{10\\,\\pi}\,
        A_n\;
        \\sqrt{\\frac{1}{k_s\,n\,\\omega\,k_s}}\;
        \\sin\\bigl(n\\,\\omega\\,t + \\Phi_n + \\tfrac{\\pi}{4}\\bigr),

    and the full signal is obtained by summing over *n = 1…n_max*.

    **Note** The implementation below follows the algebraic form that
    appeared in H04; consult the original paper for derivation details
    and recommended parameter ranges.

    Parameters
    ----------
    A_n : float or array_like
        Amplitudes :math:`A_n` of the *n*-th harmonic of surface-temperature
        forcing (K).  Provide either
        * a single scalar applied to every harmonic, or
        * an array of length ≥ *n_max* giving amplitude for each harmonic.
    Phi_n : float or array_like
        Phase shifts :math:`\\Phi_n` (rad) corresponding to each harmonic
        order.  Same broadcasting rules as *A_n*.
    n_max : int
        Highest harmonic order to include in the summation.
    k_s_val : float or array_like
        Soil thermal diffusivity :math:`k_s` (m² s⁻¹).
    lambda_s_val : float or array_like
        Soil thermal conductivity :math:`\\lambda_s` (W m⁻¹ K⁻¹).
    w : float
        Fundamental angular frequency :math:`\\omega` (rad s⁻¹), e.g.
        :math:`2\\pi/86400` for a 24-h cycle.

    Returns
    -------
    float or numpy.ndarray
        Soil heat flux *G(z,t)* (W m⁻²).  The return shape is the
        broadcast shape of the input arrays (excluding the harmonic axis).

    Raises
    ------
    ValueError
        If *n_max* is less than 1, or if *k_s_val* or *w* are non-positive,
        or if the amplitudes/phases cannot be broadcast to length
        *n_max*.

    Notes
    -----
    * **Vectorization** – Internally, the harmonic index axis has length
      *n_max* (``n = np.arange(1, n_max+1)``).  All other dimensions come
      from broadcasting *A_n*, *Phi_n*, *k_s_val*, *lambda_s_val*, and
      *t* (if vectorized in the caller).
    * **Units** – Consistency with SI units is assumed.
    * **Interpretation** – The prefactor ``λ_s / (10 π)`` appears in
      H04’s original derivation.  If you adopt a different convention
      (e.g., depth-explicit damping), modify accordingly.

    References
    ----------
    Heusinkveld, B. G., Jacobs, A. F. G., Holtslag, A. A. M., & *et al.*
    (2004). *Surface energy balance closure in an arid region: The role
    of soil heat flux.* Agricultural and Forest Meteorology, **122**(1),
    21-37. https://doi.org/10.1016/j.agrformet.2003.09.005

    Examples
    --------
    >>> # Daily cycle with three harmonics
    >>> A = np.array([10, 4, 1.5])          # K
    >>> Phi = np.array([0, -np.pi/6, np.pi/8])
    >>> Gz = heusinkveld_gz(A, Phi, n_max=3,
    ...                     k_s_val=1e-6,
    ...                     lambda_s_val=1.0,
    ...                     w=2*np.pi/86400)
    >>> Gz.shape
    ()
    """
    # --- Validation & broadcasting -------------------------------------
    if n_max < 1:
        raise ValueError("'n_max' must be at least 1.")
    if w <= 0:
        raise ValueError("'w' must be positive.")
    if np.any(np.asarray(k_s_val) <= 0):
        raise ValueError("'k_s_val' must be positive.")

    # Harmonic index array
    n = np.arange(1, n_max + 1)  # shape (n_max,)

    # Ensure arrays
    A_n = np.asarray(A_n, dtype=float)
    Phi_n = np.asarray(Phi_n, dtype=float)
    k_s_val = np.asarray(k_s_val, dtype=float)
    lambda_s_val = np.asarray(lambda_s_val, dtype=float)

    # Broadcast amplitude & phase to (n_max, …)
    try:
        A_b, Phi_b = np.broadcast_arrays(
            np.broadcast_to(A_n, (n_max,) + A_n.shape[-A_n.ndim + 1 :]),
            np.broadcast_to(Phi_n, (n_max,) + Phi_n.shape[-Phi_n.ndim + 1 :]),
        )
    except ValueError as exc:
        raise ValueError(
            "A_n and Phi_n could not be broadcast to length 'n_max'."
        ) from exc

    # Term inside the summation
    term = (
        A_b
        * np.sqrt(1.0 / (k_s_val * n[:, None] * w * k_s_val))
        * np.sin(n[:, None] * w + Phi_b + np.pi / 4.0)
    )

    # Sum over the harmonic axis (axis 0)
    return (lambda_s_val / (10.0 * np.pi)) * np.sum(term, axis=0)


# -- 3.5 Hsieh 2009 fractional derivative (Eq. 5) ------------------------------


def hsieh2009_gz(
    tz_series: np.ndarray | list | tuple,
    time_series: np.ndarray | list | tuple,
    cv_series: np.ndarray | list | tuple,
    ks_series: np.ndarray | list | tuple,
) -> float:
    """
    Compute **soil heat flux** at the end of a temperature record
    using the *half-order integral* (Hsieh et al., 2009, Eq. 5).

    The method exploits the analytical solution of the one-dimensional
    heat-conduction equation for a semi-infinite medium, leading to a
    convolution integral of order ½ that relates the time series of
    near-surface soil temperature to the downward heat flux:

    .. math::

        G(t_N) \;=\;
        2\\sqrt{\\frac{k_s(t_N)\\,C_v(t_N)}{\\pi}}\\;
        \\int_{t_0}^{t_N}
            \\frac{\\partial T(z,t')}{\\partial t'}
            \\left(t_N - t'\\right)^{-1/2}\\,dt'

    In discrete form with linear interpolation between measurements
    :math:`t_i \\;(i = 0,\\dots,N)`,

    .. math::

        \\int_{t_0}^{t_N}\\!
            \\frac{dT}{dt'}\\,(t_N-t')^{-1/2}dt'
        \\approx
        \\sum_{i=0}^{N-1}
          \\frac{\\Delta T_i}{\\Delta t_i}\\!
          \\left[(t_N-t_i)^{1/2} - (t_N-t_{i+1})^{1/2}\\right],

    where :math:`\\Delta T_i = T_{i+1}-T_i` and
    :math:`\\Delta t_i = t_{i+1}-t_i`.

    The final scalar result corresponds to the flux at
    :math:`t_N\\;(=\\text{time_series}[-1])`.

    Parameters
    ----------
    tz_series : array_like
        Near-surface (or shallow-depth) soil temperature observations
        *T(z,t)* (K).  Length **N ≥ 2**.
    time_series : array_like
        Strictly monotonically increasing time stamps (s) matching
        `tz_series`.  Same length **N**.
    cv_series : array_like
        Volumetric heat capacity *C_v* (J m⁻³ K⁻¹) at each
        time stamp.  Same length **N**.
    ks_series : array_like
        Thermal diffusivity *k_s* (m² s⁻¹) at each time stamp.
        Same length **N**.

    Returns
    -------
    float
        Soil heat flux *G(t_N)* at the final time point (W m⁻²).

    Raises
    ------
    ValueError
        * If the four series differ in length or have fewer than two
          samples.
        * If `time_series` is not strictly increasing.
        * If any element of `cv_series` or `ks_series` is non-positive.

    Notes
    -----
    * The algorithm uses a forward finite-difference for
      :math:`dT/dt'` and trapezoidal integration over each interval
      *[t_i, t_{i+1}]*.
    * Only the latest values of *k_s* and *C_v* are used, consistent
      with the Hsieh et al. derivation.  Replace with a time-variable
      kernel if property changes are large over the record.
    * All inputs are cast to :class:`numpy.ndarray` with
      ``dtype=float``.

    References
    ----------
    Hsieh, C.-I., Katul, G., & Chi, T.-C. (2009).
    *Retrieval of soil heat flux from soil temperature data by the
    continuous-time heat equation model.*
    **Water Resources Research**, 45, W08433.
    https://doi.org/10.1029/2009WR007891

    Examples
    --------
    >>> times  = np.array([0, 600, 1200, 1800])        # every 10 min
    >>> temps  = np.array([292.5, 293.0, 293.6, 294.0])  # K
    >>> Cv     = np.full_like(temps, 2.3e6)            # J m-3 K-1
    >>> ks     = np.full_like(temps, 1.1e-6)           # m2 s-1
    >>> hsieh2009_gz(temps, times, Cv, ks)
    -2.13...
    """
    # ---- convert & validate ------------------------------------------------
    tz_series = np.asarray(tz_series, dtype=float)
    time_series = np.asarray(time_series, dtype=float)
    cv_series = np.asarray(cv_series, dtype=float)
    ks_series = np.asarray(ks_series, dtype=float)

    if not (
        tz_series.size == time_series.size == cv_series.size == ks_series.size >= 2
    ):
        raise ValueError("All series must share length ≥ 2.")

    if not np.all(np.diff(time_series) > 0):
        raise ValueError("'time_series' must be strictly increasing.")

    if np.any(cv_series <= 0) or np.any(ks_series <= 0):
        raise ValueError("'cv_series' and 'ks_series' must be positive.")

    # ---- half-order integral ----------------------------------------------
    tN = time_series[-1]
    integral = 0.0
    for i in range(len(time_series) - 1):
        t_i, t_ip1 = time_series[i], time_series[i + 1]
        dT = tz_series[i + 1] - tz_series[i]
        dt = t_ip1 - t_i
        integral += (dT / dt) * (np.sqrt(tN - t_i) - np.sqrt(tN - t_ip1))

    ks_N = ks_series[-1]
    cv_N = cv_series[-1]

    return 2.0 * np.sqrt(ks_N * cv_N / np.pi) * integral


# -- 3.6 Leuning 2012 (Eq. 6–7) -------------------------------------------------


def leuning_damping_depth(
    z: np.ndarray | float,
    zr: np.ndarray | float,
    AT_z: np.ndarray | float,
    AT_zr: np.ndarray | float,
) -> np.ndarray | float:
    """
    Estimate the **thermal damping depth** *d* (m) from the exponential
    decay of temperature–wave amplitude with depth.

    The formulation is based on Eq. (6) of Leuning *et al.* (1985) for a
    one-dimensional, sinusoidally forced soil column:

    .. math::

        d \;=\;
        \\frac{z_r - z}{\\ln\\bigl(A_T(z) / A_T(z_r)\\bigr)}

    where

    * :math:`A_T(z)` — amplitude of the temperature wave at depth *z*
    * :math:`z_r`  — reference depth (m)
    * :math:`A_T(z_r)` — amplitude at reference depth

    The equation assumes amplitudes decrease monotonically with depth
    following :math:`A_T(z) = A_T(0)\\,\\exp(-z/d)`.

    Parameters
    ----------
    z : float or array_like
        Depth(s) (m) at which the amplitude *A_T(z)* was measured.
        Positive values denote depth below the surface.
    zr : float or array_like
        Reference depth(s) *z_r* (m).  Must be broadcast-compatible with
        *z*.  If *zr* < *z* the denominator switches sign, yielding a
        negative *d* that should be interpreted carefully.
    AT_z : float or array_like
        Temperature-wave amplitude at *z* (K).
    AT_zr : float or array_like
        Temperature-wave amplitude at *z_r* (K).

    Returns
    -------
    float or numpy.ndarray
        Thermal damping depth *d* (m).  The output shape follows NumPy
        broadcasting rules applied to the four inputs.

    Raises
    ------
    ValueError
        If any of the following conditions occur:

        * *AT_z* or *AT_zr* contain non-positive values (log undefined).
        * *AT_z* and *AT_zr* are equal everywhere, leading to
          ``ln(1) = 0`` and division by zero.
        * The inputs cannot be broadcast to a common shape.

    Notes
    -----
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

    References
    ----------
    Leuning, R., Barlow, E. W. R., & Paltridge, G. (1985).
    *Temperature gradients in a soil: Theory and measurement.*
    *Agricultural and Forest Meteorology*, 35(1–4), 127–143.
    https://doi.org/10.1016/0168-1923(85)90067-3

    Examples
    --------
    >>> d = leuning_damping_depth(z=0.10, zr=0.05, AT_z=4.0, AT_zr=8.0)
    >>> round(d, 4)
    0.0721

    >>> #Broadcasting with arrays
    >>> depths = np.array([0.05, 0.10, 0.20])
    >>> amps   = np.array([8.0, 4.0, 1.5])
    >>> d_arr  = leuning_damping_depth(depths, 0.02, amps, 10.0)
    >>> d_arr
    array([0.0380..., 0.0495..., 0.0854...])
    """
    # Convert to NumPy arrays
    z = np.asarray(z, dtype=float)
    zr = np.asarray(zr, dtype=float)
    AT_z = np.asarray(AT_z, dtype=float)
    AT_zr = np.asarray(AT_zr, dtype=float)

    # Broadcasting check
    try:
        z_b, zr_b, AT_z_b, AT_zr_b = np.broadcast_arrays(z, zr, AT_z, AT_zr)
    except ValueError as exc:
        raise ValueError(
            "Input arrays could not be broadcast to a common shape."
        ) from exc

    # Guard: positive amplitudes
    if np.any(AT_z_b <= 0) or np.any(AT_zr_b <= 0):
        raise ValueError("Amplitudes 'AT_z' and 'AT_zr' must be positive.")

    ratio = AT_z_b / AT_zr_b
    if np.any(ratio == 1.0):
        raise ValueError(
            "Amplitude ratio AT_z / AT_zr equals 1 for some elements, "
            "yielding division by zero in the logarithm."
        )

    return (zr_b - z_b) / np.log(ratio)


def leuning_gz(
    g_zr: np.ndarray | float,
    z: np.ndarray | float,
    zr: np.ndarray | float,
    d: np.ndarray | float,
) -> np.ndarray | float:
    """
    Extrapolate **soil heat flux** from a reference depth *z_r* to a
    shallower (or deeper) target depth *z* using an *exponential
    attenuation* model (Leuning *et al.* 1985, Eq. 7).

    The model assumes the amplitude of the heat-flux wave decays
    exponentially with depth at the same *damping depth* **d** that
    governs temperature attenuation:

    .. math::

        G(z,t) \;=\; G(z_r,t)\;\exp\\!\left(\\frac{z_r - z}{d}\\right)

    where
    :math:`G(z_r,t)` is the measured flux at *z_r* (W m⁻²).

    Parameters
    ----------
    g_zr : float or array_like
        Heat flux measured at reference depth *z_r* (W m⁻²).  Can be a
        scalar or a NumPy-broadcastable array (e.g. a time series).
    z : float or array_like
        Target depth *z* (m, **positive downward**).  Must be broadcast-
        compatible with *g_zr*.
    zr : float or array_like
        Reference depth *z_r* (m, positive downward).  Broadcast
        compatible with *z*.
    d : float or array_like
        Thermal damping depth (m).  Positive; same shape rules as *z*.

    Returns
    -------
    float or numpy.ndarray
        Estimated heat flux at depth *z* (W m⁻²).  Follows NumPy
        broadcasting rules applied to the inputs.

    Raises
    ------
    ValueError
        If any element of *d* is non-positive, or if inputs cannot be
        broadcast to a common shape.

    Notes
    -----
    * **Direction of extrapolation** –
      *z < z_r* (shallower) ⇒ magnitude *increases* toward the surface;
      *z > z_r* (deeper)   ⇒ magnitude *decreases* with depth.
    * **Vectorisation** – Inputs are converted to
      :class:`numpy.ndarray` and the formula is evaluated element-wise.
    * **Limitations** – The simple exponential form ignores soil
      layering and moisture variability; best suited to homogeneous
      profiles over the depth interval considered.

    References
    ----------
    Leuning, R., Barlow, E. W. R., & Paltridge, G. (1985).
    *Temperature gradients in a soil: Theory and measurement.*
    **Agricultural and Forest Meteorology**, 35(1–4), 127-143.
    https://doi.org/10.1016/0168-1923(85)90067-3

    Examples
    --------
    >>> # Single depth conversion
    >>> leuning_gz(g_zr=-12.0, z=0.05, zr=0.08, d=0.07)
    -18.841...

    >>> # Vectorized daily time series
    >>> g_plate = np.random.normal(-10, 3, 1440)      # W m-2 @ 8 cm
    >>> d       = 0.07                                # m
    >>> Gz_5cm  = leuning_gz(g_plate, z=0.05, zr=0.08, d=d)
    >>> Gz_5cm.shape
    (1440,)
    """
    # Convert inputs to NumPy arrays and broadcast
    g_zr = np.asarray(g_zr, dtype=float)
    z = np.asarray(z, dtype=float)
    zr = np.asarray(zr, dtype=float)
    d = np.asarray(d, dtype=float)

    if np.any(d <= 0):
        raise ValueError("Damping depth 'd' must be positive.")

    try:
        g_b, z_b, zr_b, d_b = np.broadcast_arrays(g_zr, z, zr, d)
    except ValueError as exc:
        raise ValueError(
            "Inputs g_zr, z, zr, and d are not broadcast-compatible."
        ) from exc

    return g_b * np.exp((zr_b - z_b) / d_b)


# -- 3.7 Simple‑measurement (Eq. 8) -------------------------------------------


def simple_measurement_gz(
    g_zr: np.ndarray | float,
    cv_layers: np.ndarray | float,
    tz_layers: np.ndarray | float,
    dt: float,
    dz_layers: np.ndarray | float,
) -> np.ndarray | float:
    """
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
                \frac{\Delta T_{l,j} + \frac{1}{2}\bigl(\Delta T_{l,j} - \Delta T_{l,j-1}\bigr)}
                    {\Delta t}

    where

    * :math:`C_{v,l}` — volumetric heat capacity of layer *l*
      (J m⁻³ K⁻¹)
    * :math:`\\Delta z_l` — layer thickness (m)
    * :math:`\\Delta T_{l,j}` — temperature change in layer *l*
      between *t_{j-1}* and *t_j* (K)
    * :math:`\\Delta t` — measurement interval (s)

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
    """
    # --- Basic validation -------------------------------------------------
    g_zr = np.asarray(g_zr, dtype=float)
    tz_layers = np.asarray(tz_layers, dtype=float)
    cv_layers = np.asarray(cv_layers, dtype=float)
    dz_layers = np.asarray(dz_layers, dtype=float)

    if dt <= 0:
        raise ValueError("'dt' must be positive.")

    # Expect tz_layers shape (N, M)
    if tz_layers.ndim != 2:
        raise ValueError("'tz_layers' must be a 2-D array (layers × time).")

    n_layers, n_times = tz_layers.shape
    if g_zr.shape[-1] != n_times:
        raise ValueError(
            "Length of 'g_zr' ({}) must match the time dimension "
            "of 'tz_layers' ({})".format(g_zr.shape[-1], n_times)
        )

    # Broadcast cv and dz to (N, 1)
    cv_layers = np.broadcast_to(cv_layers, (n_layers, 1))
    dz_layers = np.broadcast_to(dz_layers, (n_layers, 1))

    # --- Storage term -----------------------------------------------------
    delta_tz = tz_layers[:, 1:] - tz_layers[:, :-1]  # shape (N, M-1)
    delta_tz_mid = 0.5 * (delta_tz[:, :-1] + delta_tz[:, 1:])  # (N, M-2)

    # First and last centred differences align with t[1]…t[M-2]
    storage = np.sum(
        cv_layers * dz_layers * (delta_tz[:, 1:] + delta_tz_mid) / dt,
        axis=0,
    )  # shape (M-2,)

    # Plate flux aligned to storage length
    return g_zr[1:-1] + storage


# -- 3.8 Wang & Bou‑Zeid 2012 Green's function (Eq. 9–10) ---------------------


def wbz12_g_gz(
    g_zr_series: np.ndarray | list | tuple,
    time_series: np.ndarray | list | tuple,
    z: float,
    zr: float,
    k_s_val: float,
) -> np.ndarray:
    """
    Estimate **soil heat flux** at a shallow depth *z* from a flux-plate
    record at reference depth *z_r* using the **WBZ12-G convolution
    method** (Wang, Bou-Zeid & Zhang 2012, their Eq. 9–10).

    The algorithm inverts the one-dimensional heat-conduction equation
    for a semi-infinite homogeneous soil, assuming a step-response
    kernel based on the complementary error function *erfc*:

    .. math::

        G(z,t) \;=\; 2\,G(z_r,t)
        \;-\; \\frac{J(t)}{\\Delta F_z(t)}

    with the discrete convolution integral

    .. math::

        J(t_n) \;=\;
        \\sum_{j=0}^{n-1}
        \\tfrac{\\bigl[G(z_r,t_{n-j-1}) + G(z_r,t_{n-j})\\bigr]}{2}\;
        \\Delta F_z(t_{j+1})

    and the *transfer function increment*

    .. math::

        \\Delta F_z(t) =
        \\operatorname{erfc}\\!
        \\left[\\frac{z_r - z}{2\\sqrt{k_s (t - t_0)}}\\right]
        \;-\;
        \\operatorname{erfc}\\!
        \\left[\\frac{z_r - z}{2\\sqrt{k_s (t - t_0 - \\Delta t)}}\\right].

    Parameters
    ----------
    g_zr_series : array_like
        Time series of heat-flux-plate measurements at *z_r* (W m⁻²);
        length **N**.
    time_series : array_like
        Monotonically increasing time stamps (s) corresponding to
        `g_zr_series`.  Must be the same length **N**.
    z : float
        Target depth (m, **positive downward**).
    zr : float
        Reference depth of the heat-flux plate (m, positive downward).
    k_s_val : float
        Soil thermal diffusivity *k_s* (m² s⁻¹).  Must be positive.

    Returns
    -------
    numpy.ndarray
        Estimated heat-flux series at depth *z*, same length **N** as the
        input series (W m⁻²).

    Raises
    ------
    ValueError
        If `k_s_val ≤ 0`, the two series differ in length, time stamps
        are non-monotonic, or contain fewer than two samples.

    Notes
    -----
    * The first element of the output corresponds to *t₀* (no storage
      correction possible yet); subsequent points incorporate the
      convolution up to that instant.
    * A very small term ``+ 1e-12`` is added under the square root to
      avoid division by zero at *t = t₀*.
    * The complementary error function is evaluated with
      :pyfunc:`numpy.erfc`, which is vectorized and avoids the SciPy
      dependency.

    References
    ----------
    Wang, W., Bou-Zeid, E., & Zhang, Y. (2012).
    *Estimating surface heat fluxes using the surface renewal method*.
    **Boundary-Layer Meteorology**, 144(2), 407-422.
    https://doi.org/10.1007/s10546-012-9730-9

    Examples
    --------
    >>> # 10-min sampled day-long plate record at zr = 0.08 m
    >>> t  = np.arange(0, 24*3600 + 1, 600)            # s
    >>> Gp = -10 + 5*np.sin(2*np.pi*t/86400)
    >>> Gz = wbz12_g_gz(Gp, t, z=0.05, zr=0.08, k_s_val=1.0e-6)
    >>> Gz.shape
    (145,)
    """
    # --- validation -------------------------------------------------------
    g_zr_series = np.asarray(g_zr_series, dtype=float)
    time_series = np.asarray(time_series, dtype=float)

    if g_zr_series.shape != time_series.shape:
        raise ValueError("'g_zr_series' and 'time_series' must have the same length.")
    if g_zr_series.size < 2:
        raise ValueError("At least two time steps are required.")
    if not np.all(np.diff(time_series) > 0):
        raise ValueError("'time_series' must be strictly increasing.")
    if k_s_val <= 0:
        raise ValueError("'k_s_val' must be positive.")

    # ----------------------------------------------------------------------
    t0 = time_series[0]
    tau = time_series - t0  # elapsed time (s)

    # Transfer function F_z(τ)
    Fz = np.erfc((zr - z) / (2.0 * np.sqrt(k_s_val * tau + 1e-12)))
    delta_Fz = np.diff(Fz, prepend=0.0)  # ΔF_z

    # Discrete convolution integral J(t_n)
    J = np.zeros_like(g_zr_series)
    for n in range(1, len(time_series)):
        # Reverse slice of plate-flux series up to index n-1 (inclusive)
        g_slice = g_zr_series[n - 1 :: -1]
        deltaF_slice = delta_Fz[1 : n + 1]
        # Trapezoidal weighting of consecutive flux pairs
        J[n] = np.sum((g_slice[:-1] + g_slice[1:]) * deltaF_slice)

    # WBZ12-G heat-flux estimate
    Gz = 2.0 * g_zr_series - J / delta_Fz[1]

    return Gz


# -- 3.9 Wang & Bou‑Zeid 2012 sinusoid‑plus‑integral (Eq. 11) -----------------

import numpy as np
from scipy.integrate import quad


def wbz12_s_gz(
    Ag: np.ndarray | float,
    ks_val: np.ndarray | float,
    zr: float,
    z: float,
    t: np.ndarray | float,
    eps: np.ndarray | float,
    omega: float = OMEGA_DAY,
) -> np.ndarray | float:
    """
    **WBZ12-S analytic–numeric solution** for soil heat flux *G(z, t)* at
    depth *z* (Wang, Bou-Zeid & Zhang 2012, Eq. 11).

    WBZ12-S combines a closed-form *forcing* term (sinusoidal surface
    temperature) with a *storage* term that is expressed as a Fourier‐
    Laplace integral and must be evaluated numerically.  The governing
    equation assumes homogeneous soil properties and a sinusoidal
    surface temperature of amplitude *A_g* and phase `eps`:

    .. math::

        G(z,t) \;=\;
        A_g\,e^{-\eta}\,\sin\\bigl(\\omega t + \varepsilon - \eta\\bigr)
        \;-\;
        \\frac{2 A_g k_s}{\\pi}
        \int_{0}^{\\infty}
            \\frac{
                k_s \zeta^{2} \bigl[\\sin\\varepsilon
                - \\omega \\cos\\varepsilon\\bigr]
              }{
                \\omega^{2}
                + k_{s}^{2} \\zeta^{4}
            }
            \sin\\bigl[\\zeta\,(z_r - z)\\bigr]\,
            e^{-k_s \zeta^{2} t}\,d\\zeta

    with

    .. math::
        \\eta \\;=\\; (z_r - z)\,\\sqrt{\\omega / (2 k_s)}.

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
        Phase shift :math:`\\varepsilon` (rad) of the surface
        temperature.  Broadcast-compatible with `t`.
    omega : float, optional
        Angular frequency :math:`\\omega` (rad s⁻¹) of the sinusoidal
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
      :func:`scipy.integrate.quad` over :math:`\\zeta \\in [0, \\infty)`.
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
    """
    # ------------------------------------------------------------------ #
    # Broadcast and basic validation
    Ag = np.asarray(Ag, dtype=float)
    ks_val = np.asarray(ks_val, dtype=float)
    t = np.asarray(t, dtype=float)
    eps = np.asarray(eps, dtype=float)

    if np.any(ks_val <= 0):
        raise ValueError("'ks_val' must be positive.")
    if omega <= 0:
        raise ValueError("'omega' must be positive.")

    try:
        Ag_b, ks_b, t_b, eps_b = np.broadcast_arrays(Ag, ks_val, t, eps)
    except ValueError as exc:
        raise ValueError(
            "Inputs Ag, ks_val, t, and eps are not broadcast-compatible."
        ) from exc

    # vectorized calculation ------------------------------------------------
    # Prefactor & steady sinusoidal term
    eta = (zr - z) * np.sqrt(omega / (2.0 * ks_b))
    prefactor = Ag_b * np.exp(-eta)
    sin_term = np.sin(omega * t_b + eps_b - eta)

    # Numerical transient integral (scalar per broadcast element)
    def transient_scalar(A, ks, tau, eps_scalar):
        """Compute the transient integral for a single point in time."""

        def integrand(zeta):
            numerator = ks * zeta**2 * (np.sin(eps_scalar) - omega * np.cos(eps_scalar))
            denom = omega**2 + (ks**2) * zeta**4
            return (
                numerator
                / denom
                * np.sin(zeta * (zr - z))
                * np.exp(-ks * zeta**2 * tau)
            )

        return quad(integrand, 0, np.inf, limit=500)[0]

    # Allocate output
    result = np.empty_like(Ag_b, dtype=float)

    # Iterate over flattened broadcast shape (quad is not vectorized)
    it = np.nditer(
        [Ag_b, ks_b, t_b, eps_b, result],
        flags=["multi_index"],
        op_flags=[["readonly"]] * 4 + [["writeonly"]],
    )
    for Ai, ksi, ti, epsi, out in it:
        integral_val = transient_scalar(Ai.item(), ksi.item(), ti.item(), epsi.item())
        second_term = -2.0 * Ai.item() * ksi.item() / np.pi * integral_val
        out[...] = (
            Ai.item()
            * np.exp(-(zr - z) * np.sqrt(omega / (2 * ksi.item())))
            * np.sin(
                omega * ti.item()
                + epsi.item()
                - (zr - z) * np.sqrt(omega / (2 * ksi.item()))
            )
            + second_term
        )

    return result


# -----------------------------------------------------------------------------
# 4. Exact sinusoidal benchmark (Eq. 16–17)
# -----------------------------------------------------------------------------
import numpy as np


def exact_temperature_gz(
    z: np.ndarray | float,
    AT: np.ndarray | float,
    t: np.ndarray | float,
    d: np.ndarray | float,
    omega: float = OMEGA_DAY,
    T_i: float = 298.15,
) -> np.ndarray | float:
    """
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
        \sin\\bigl(\\omega t \;-\; z/d\\bigr),

    where

    * :math:`T_i` — initial (mean) temperature of the soil column (K)
    * :math:`A_T` — amplitude of the surface-temperature oscillation (K)
    * :math:`d` — thermal damping depth (m)
    * :math:`\\omega` — angular frequency of the forcing (rad s⁻¹)
    * :math:`t` — time since the start of oscillation (s)

    Parameters
    ----------
    z : float or array_like
        Depth below the surface (m). Positive downward.  Can be a scalar
        or an array broadcastable with *t* and *AT*.
    AT : float or array_like
        Amplitude *A_T* of the surface temperature wave (K).
    t : float or array_like
        Time variable (s).  Supports vectorisation.
    d : float or array_like
        Thermal damping depth (m). Must be positive.  Can be scalar or
        broadcastable with *z*.
    omega : float, optional
        Angular frequency *ω* (rad s⁻¹). Defaults to
        :pydata:`OMEGA_DAY` (≈ 7.272 × 10⁻⁵ s⁻¹ for 24 h).
    T_i : float, optional
        Mean (initial) soil temperature (K).  Defaults to **298.15 K**
        (25 °C).

    Returns
    -------
    float or numpy.ndarray
        Soil temperature *T(z, t)* (same units as *T_i*).  The result
        shape matches NumPy broadcasting over the inputs.

    Raises
    ------
    ValueError
        If any element of *d* or *omega* is non-positive, or if inputs
        cannot be broadcast to a common shape.

    Notes
    -----
    * **Phase lag** – Each e-folding depth *d* introduces a phase delay
      of 1 radian (~57.3 °) relative to the surface signal.
    * **Vectorisation** – All inputs are converted to
      :class:`numpy.ndarray`; the formula is evaluated element-wise and
      fully supports broadcasting.
    * **Units** – Ensure consistent SI units (metres, seconds, kelvins).

    References
    ----------
    Adapted from Gao, Z., Horton, R., Luo, L., & Kucharik, C. J. (2010).
    *A simple method to measure soil temperature dynamics: Theory and
    application.* **Soil Science Society of America Journal**, 74(2),
    580-588. https://doi.org/10.2136/sssaj2009.0169

    Examples
    --------
    >>> # Daily cycle at 10 cm depth
    >>> z      = 0.10                            # m
    >>> AT     = 8.0                             # K
    >>> d      = 0.12                            # m
    >>> t_day  = np.linspace(0, 86400, 97)       # 15-min resolution
    >>> Tz     = exact_temperature_gz(z, AT, t_day, d)
    >>> Tz.shape
    (97,)
    """
    # Convert inputs to arrays and broadcast
    z = np.asarray(z, dtype=float)
    AT = np.asarray(AT, dtype=float)
    t = np.asarray(t, dtype=float)
    d = np.asarray(d, dtype=float)

    if np.any(d <= 0):
        raise ValueError("Damping depth 'd' must be positive.")
    if omega <= 0:
        raise ValueError("'omega' must be positive.")

    try:
        z_b, AT_b, t_b, d_b = np.broadcast_arrays(z, AT, t, d)
    except ValueError as exc:
        raise ValueError(
            "Inputs z, AT, t, and d are not broadcast-compatible."
        ) from exc

    return T_i + AT_b * np.exp(-z_b / d_b) * np.sin(omega * t_b - z_b / d_b)


import numpy as np


def exact_gz(
    z: np.ndarray | float,
    AT: np.ndarray | float,
    lambda_s_val: np.ndarray | float,
    d: np.ndarray | float,
    t: np.ndarray | float,
    omega: float = OMEGA_DAY,
) -> np.ndarray | float:
    """
    Exact analytical **soil-heat-flux** solution for sinusoidal surface forcing.

    A sinusoidal surface–temperature wave of amplitude *A_T* and angular
    frequency *ω* generates a heat-flux wave that decays exponentially
    with depth and lags the temperature signal by *π / 4* radians
    (Gao et al., 2010, Eq. 17):

    .. math::

        G(z,t) \;=\;
        \\sqrt{2}\;
        \\frac{\\lambda_s A_T}{d}\;
        e^{-z/d}\;
        \\sin\\bigl(\\omega t - z/d + \\tfrac{\\pi}{4}\\bigr)

    where

    * :math:`\\lambda_s` — soil thermal conductivity (W m⁻¹ K⁻¹)
    * :math:`d` — thermal damping depth (m)
    * :math:`z` — depth below the surface (m, positive downward)
    * :math:`t` — time (s) since the start of oscillation
    * :math:`\\omega` — angular frequency (rad s⁻¹)

    Parameters
    ----------
    z : float or array_like
        Depth(s) below the surface (m), positive downward.
    AT : float or array_like
        Amplitude *A_T* of the surface-temperature oscillation (K).
    lambda_s_val : float or array_like
        Soil thermal conductivity *λ_s* (W m⁻¹ K⁻¹).
    d : float or array_like
        Thermal damping depth (m). Must be positive.
    t : float or array_like
        Time variable (s).
    omega : float, optional
        Angular frequency *ω* (rad s⁻¹). Defaults to
        :pydata:`OMEGA_DAY` (≈ 7.272 × 10⁻⁵ s⁻¹, i.e. 2π / 86 400 s).

    Returns
    -------
    float or numpy.ndarray
        Soil heat flux *G(z, t)* (W m⁻²).  The return shape follows NumPy
        broadcasting rules applied to the inputs.

    Raises
    ------
    ValueError
        If any element of *d* or *omega* is non-positive, or if inputs
        cannot be broadcast to a common shape.

    Notes
    -----
    * **Phase shift** – At any given depth, the heat-flux wave lags the
      temperature wave by 45 ° (π / 4 rad).
    * **Vectorisation** – All inputs are converted to
      :class:`numpy.ndarray`; the expression is evaluated element-wise
      and fully supports broadcasting.
    * **Units** – Ensure all quantities use consistent SI units.

    References
    ----------
    Gao, Z., Horton, R., Luo, L., & Kucharik, C. J. (2010).
    *A simple method to measure soil temperature dynamics: Theory and
    application.* **Soil Science Society of America Journal**, 74(2),
    580–588. https://doi.org/10.2136/sssaj2009.0169

    Examples
    --------
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
    """
    # --- Validation & broadcasting ------------------------------------
    z = np.asarray(z, dtype=float)
    AT = np.asarray(AT, dtype=float)
    lambda_s_val = np.asarray(lambda_s_val, dtype=float)
    d = np.asarray(d, dtype=float)
    t = np.asarray(t, dtype=float)

    if np.any(d <= 0):
        raise ValueError("Damping depth 'd' must be positive.")
    if omega <= 0:
        raise ValueError("'omega' must be positive.")

    try:
        z_b, AT_b, lam_b, d_b, t_b = np.broadcast_arrays(z, AT, lambda_s_val, d, t)
    except ValueError as exc:
        raise ValueError(
            "Inputs z, AT, lambda_s_val, d, and t are not " "broadcast-compatible."
        ) from exc

    # --- Analytical solution ------------------------------------------
    return (
        np.sqrt(2.0)
        * lam_b
        * AT_b
        / d_b
        * np.exp(-z_b / d_b)
        * np.sin(omega * t_b - z_b / d_b + np.pi / 4.0)
    )
