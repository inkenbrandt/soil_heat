import numpy as np
import pandas as pd

# ============================================================
# CONFIGURATION
# ============================================================

# SoilVUE depths
DEPTHS_CM = np.array([5, 10, 20, 30, 40, 50, 60], dtype=float)
DEPTHS_M = DEPTHS_CM / 100.0

# Layer boundaries: midpoints between sensors, with surface and bottom extension
LAYER_BOUNDS_CM = np.array([0, 7.5, 15, 25, 35, 45, 55, 65], dtype=float)
LAYER_THICK_M = np.diff(LAYER_BOUNDS_CM) / 100.0

# Thermal constants
CV_W = 4.186e6  # water volumetric heat capacity (J/(m³·K))
CV_S = 2.0e6  # mineral solids (J/(m³·K))
CV_A = 1250.0  # air (J/(m³·K))
CV_ORG = 2.5e6  # organic matter (J/(m³·K))
LAMBDA_W = 0.594  # water thermal conductivity (W/(m·K))
LAMBDA_S = 2.5  # quartz (W/(m·K))

# Default soil assumptions
POROSITY = 0.42
F_QUARTZ = 0.5
F_ORGANIC = 0.02
DT_SEC = 3600.0  # hourly data

# Assumed heat flux plate depth (cm) for storage comparison
PLATE_DEPTH_CM = 8.0


# ============================================================
# THERMAL PROPERTIES
# ============================================================


def volumetric_heat_capacity(theta, porosity=POROSITY, f_org=F_ORGANIC):
    """
    Estimate soil volumetric heat capacity using the de Vries (1963) linear
    mixing model.

    Computes Cv as a volume-fraction-weighted sum of the heat capacities of
    the four soil constituents — mineral solids, organic matter, liquid water,
    and air:

        Cv(θ) = f_min·Cv_s + f_org·Cv_org + θ·Cv_w + θ_a·Cv_a

    where the air-filled porosity θ_a = max(porosity - θ, 0).

    Parameters
    ----------
    theta : array-like
        Volumetric water content [m³ m⁻³]. Accepts scalars, lists, or NumPy
        arrays. Values above porosity are implicitly handled by flooring
        air-filled porosity at zero.
    porosity : float, optional
        Total soil porosity [m³ m⁻³]. Used to derive both the mineral solid
        fraction (f_min = 1 - porosity - f_org) and the air-filled pore
        fraction. Defaults to the module-level constant ``POROSITY``.
    f_org : float, optional
        Volumetric fraction of organic matter in the solid phase [m³ m⁻³].
        Defaults to the module-level constant ``F_ORGANIC``.

    Returns
    -------
    np.ndarray
        Volumetric heat capacity Cv [J m⁻³ K⁻¹] with the same shape as
        ``theta``. Typical values range from ~1.0×10⁶ J m⁻³ K⁻¹ (dry) to
        ~3.0×10⁶ J m⁻³ K⁻¹ (saturated mineral soil).

    Notes
    -----
    Component heat capacities are taken from module-level constants:

    - ``CV_S``   — mineral solids  (~2.0×10⁶ J m⁻³ K⁻¹)
    - ``CV_ORG`` — organic matter  (~2.5×10⁶ J m⁻³ K⁻¹)
    - ``CV_W``   — liquid water    (~4.18×10⁶ J m⁻³ K⁻¹)
    - ``CV_A``   — air             (~1.2×10³ J m⁻³ K⁻¹)

    The air contribution is negligible at most water contents but is retained
    for physical completeness. The model assumes no phase change and is not
    valid for frozen soils where ice content would require a separate term.

    References
    ----------
    de Vries, D. A. (1963). Thermal properties of soils. In W. R. van Wijk
    (Ed.), *Physics of Plant Environment* (pp. 210–235). North-Holland,
    Amsterdam.
    """
    f_min = 1.0 - porosity - f_org
    theta_a = np.maximum(porosity - theta, 0)
    return f_min * CV_S + f_org * CV_ORG + theta * CV_W + theta_a * CV_A


def thermal_conductivity_johansen(theta, porosity=POROSITY, f_quartz=F_QUARTZ):
    """
    Estimate soil thermal conductivity using the Johansen (1975) model with the
    Kersten number formulation from Côté & Konrad (2005).

    Computes λ by interpolating between dry and saturated end-member conductivities
    using the Kersten number (Ke) as a saturation-dependent weighting function:

        λ(θ) = λ_dry + Ke(Sr) · (λ_sat - λ_dry)

    End members and saturation are derived from soil texture (quartz fraction),
    porosity, and standard mineral/water conductivity constants.

    Parameters
    ----------
    theta : array-like
        Volumetric water content [m³ m⁻³]. Accepts scalars, lists, or NumPy
        arrays. Clipped internally to the range [0, porosity] via the degree
        of saturation calculation.
    porosity : float, optional
        Total soil porosity [m³ m⁻³], i.e., the volumetric fraction of pore
        space. Defaults to the module-level constant ``POROSITY``.
    f_quartz : float, optional
        Quartz fraction of the solid phase [dimensionless, 0–1]. Higher quartz
        content increases the solid-phase conductivity λ_s and therefore λ_sat.
        Defaults to the module-level constant ``F_QUARTZ``.

    Returns
    -------
    np.ndarray
        Thermal conductivity λ [W m⁻¹ K⁻¹] with the same shape as ``theta``.

    Notes
    -----
    The three derived end-member quantities are:

    - **λ_dry**: dry soil conductivity estimated from bulk density using the
      empirical formula of Johansen (1975), assuming a particle density of
      2700 kg m⁻³.
    - **λ_s**: solid-phase conductivity computed as a geometric mean weighted
      by quartz fraction, using module constants ``LAMBDA_S`` (mineral) and
      2.0 W m⁻¹ K⁻¹ for non-quartz minerals.
    - **λ_sat**: saturated conductivity computed as a geometric mean of the
      solid and water (``LAMBDA_W``) phase conductivities weighted by their
      volume fractions (1 - porosity and porosity, respectively).

    The Kersten number follows the log-linear relationship for unfrozen mineral
    soils from Côté & Konrad (2005):

        Ke = 0.7 · log10(Sr) + 1.0   for Sr > 0.05
        Ke = 0.0                       for Sr ≤ 0.05

    Ke is clipped to [0, 1] and Sr is floored at 0.01 inside the log to
    prevent numerical issues at very low saturation. This formulation is valid
    for unfrozen, mineral soils; it is not appropriate for frozen conditions or
    highly organic soils.

    References
    ----------
    Johansen, O. (1975). *Thermal conductivity of soils*. PhD thesis,
    University of Trondheim, Norway. (CRREL Draft Translation 637, 1977).

    Côté, J., & Konrad, J.-M. (2005). A generalized thermal conductivity model
    for soils and construction materials. *Canadian Geotechnical Journal*,
    42(2), 443–458. https://doi.org/10.1139/t04-106
    """
    rho_s = 2700.0  # particle density
    lambda_dry = (0.135 * rho_s * (1 - porosity) + 64.7) / (
        rho_s - 0.947 * rho_s * (1 - porosity)
    )
    lambda_s = LAMBDA_S**f_quartz * 2.0 ** (1 - f_quartz)
    lambda_sat = lambda_s ** (1 - porosity) * LAMBDA_W**porosity

    Sr = np.clip(theta / porosity, 0, 1)
    Ke = np.where(Sr > 0.05, 0.7 * np.log10(np.maximum(Sr, 0.01)) + 1.0, 0.0)
    Ke = np.clip(Ke, 0, 1)

    return lambda_dry + Ke * (lambda_sat - lambda_dry)


def thermal_conductivity_lu2007(theta):
    """
    Estimate soil thermal conductivity using the Lu et al. (2007) empirical model.

    Computes λ as a function of volumetric water content using a nonlinear
    empirical equation fitted across a wide range of mineral soil textures:

        λ(θ) = a + b·θ - (a - d) · exp(-(c·θ)^e)

    with coefficients a=0.56, b=0.51, c=0.96, d=0.27, e=0.87.

    Parameters
    ----------
    theta : array-like
        Volumetric water content [m³ m⁻³]. Accepts scalars, lists, or NumPy
        arrays. Values are expected in the range [0, porosity]; behavior outside
        this range is not physically meaningful.

    Returns
    -------
    np.ndarray
        Thermal conductivity λ [W m⁻¹ K⁻¹] with the same shape as ``theta``.

    Notes
    -----
    The model was calibrated on 36 soils spanning sandy to clayey textures and
    dry to saturated conditions, with reported RMSE of ~0.1 W m⁻¹ K⁻¹. It does
    not account for soil texture, organic matter content, or temperature
    explicitly, so accuracy may degrade for organic-rich or highly heterogeneous
    soils.

    References
    ----------
    Lu, S., Ren, T., Gong, Y., & Horton, R. (2007). An improved model for
    predicting soil thermal conductivity from water content at room temperature.
    *Soil Science Society of America Journal*, 71(1), 8–14.
    https://doi.org/10.2136/sssaj2006.0041
    """
    a, b, c, d, e = 0.56, 0.51, 0.96, 0.27, 0.87
    return a + b * theta - (a - d) * np.exp(-((c * theta) ** e))


def gradient_surface(ts, swc, df):
    """
    Estimate surface ground heat flux using a temperature gradient from the canopy
    surface down to the 5 cm soil sensor (Method 2).

    Applies Fourier's law over the 0 → 5 cm layer, treating the canopy air
    temperature as a proxy for the soil surface temperature:

        G0 = -λ(θ_5cm) · (T_5cm - T_surface) / dz

    where dz = 0.05 m.

    Parameters
    ----------
    ts : pd.DataFrame
        Soil temperature time series [°C], with column keys as sensor depths in cm.
        Must contain a column for the 5.0 cm depth.
    swc : pd.DataFrame
        Volumetric soil water content time series [m³ m⁻³], with column keys matching
        those of ``ts``. The 5.0 cm column is used to compute thermal conductivity.
    df : pd.DataFrame
        Main eddy covariance / met data DataFrame. Must contain the column
        ``"T_CANOPY_1_1_1"`` [°C], used as a proxy for the soil surface temperature.

    Returns
    -------
    pd.Series
        Estimated surface ground heat flux [W m⁻²], indexed to match ``ts``.
        Named ``"G_gradient_surface"``. Positive values indicate downward flux
        (heat flowing into the soil).

    Notes
    -----
    Using canopy air temperature as a surface proxy introduces uncertainty under
    conditions where the canopy-to-soil surface temperature decouples significantly
    (e.g., high radiative loading, dense litter layers, or sparse canopy cover).
    Thermal conductivity is estimated at the 5 cm depth via the Johansen (1975)
    model, which may underestimate conductivity in very dry or coarse-textured soils.

    References
    ----------
    Johansen, O. (1975). *Thermal conductivity of soils*. PhD thesis, University
    of Trondheim, Norway.
    """
    T_surface = df["T_CANOPY_1_1_1"]
    T_5cm = ts[5.0]
    theta_surf = swc[5.0]
    lam = thermal_conductivity_johansen(theta_surf.values)
    G = -lam * (T_5cm.values - T_surface.values) / 0.05
    return pd.Series(G, index=ts.index, name="G_gradient_surface")


# ============================================================
# SOIL HEAT FLUX METHODS
# ============================================================
def compute_storage_above_plate(ts, swc, plate_depth_cm=PLATE_DEPTH_CM, dt=DT_SEC):
    """
    Compute soil heat storage in the layer from the surface to the heat flux plate depth.

    Integrates the calorimetric storage term (Cv * dT/dt * dz) over the 0 → plate_depth
    layer using the 5 cm temperature and soil water content sensors as representative
    values for that layer. The result can be added to the plate-measured soil heat flux
    (SG) to estimate total ground heat flux at the surface (G0 = SG + S_above_plate).

    Since the plate is typically installed at 8 cm — between the 5 cm and 10 cm sensors —
    the 5 cm sensor is used as a pragmatic representative of the above-plate layer rather
    than performing a full depth interpolation.

    Parameters
    ----------
    ts : pd.DataFrame
        Soil temperature time series [°C], with column keys as sensor depths in cm
        (e.g., 5.0, 10.0). Must contain a column for the 5.0 cm depth.
    swc : pd.DataFrame
        Volumetric soil water content time series [m³ m⁻³], with column keys matching
        those of ``ts``. Must contain a column for the 5.0 cm depth.
    plate_depth_cm : float, optional
        Depth of the soil heat flux plate [cm]. Defines the thickness of the storage
        layer (dz = plate_depth_cm / 100). Defaults to ``PLATE_DEPTH_CM``.
    dt : float, optional
        Timestep duration [s] used to compute dT/dt via finite differencing.
        Defaults to ``DT_SEC``.

    Returns
    -------
    pd.Series
        Heat storage rate in the 0 → plate_depth layer [W m⁻²], indexed to match
        ``ts``. Named ``"S_above_plate"``. Positive values indicate heat being stored
        (warming), negative values indicate heat being released (cooling). The first
        element will be NaN due to the forward-difference computation.

    Notes
    -----
    Storage is calculated as:

        S = Cv(θ) · (ΔT / Δt) · dz

    where Cv is the volumetric heat capacity [J m⁻³ K⁻¹] derived from the soil water
    content at 5 cm, ΔT is the temperature change over one timestep, and dz is the
    layer thickness in meters.
    """
    # Simple approach: use the 5 cm sensor to represent the 0→plate_depth layer
    # Layer thickness = plate_depth / 100 (m)
    dz = plate_depth_cm / 100.0

    dTdt_5cm = ts[5.0].diff() / dt
    cv_5cm = volumetric_heat_capacity(swc[5.0].values)

    S_plate = cv_5cm * dTdt_5cm.values * dz
    return pd.Series(S_plate, index=ts.index, name="S_above_plate")


def gradient_plus_storage(ts, swc, ref_depth_idx=2, dt=DT_SEC, lam_model="johansen"):
    """
    Estimate surface ground heat flux using a conductive gradient at a reference depth
    plus calorimetric storage integrated over all layers above it.

    Applies the standard soil heat flux correction:

        G0 = G(z_ref) + S(0 → z_ref)

    The conductive flux at the reference depth is estimated from Fourier's law using
    temperatures at the sensor depths bracketing ``z_ref`` (i.e., one layer above and
    one layer below). The storage term is computed calorimetrically for each discrete
    layer from the surface down to and including ``z_ref``, then summed.

    Parameters
    ----------
    ts : pd.DataFrame
        Soil temperature time series [°C], with column keys as sensor depths in cm
        matching ``DEPTHS_CM``. Must span at least the depths indexed by
        ``ref_depth_idx - 1`` through ``ref_depth_idx + 1``.
    swc : pd.DataFrame
        Volumetric soil water content time series [m³ m⁻³], with column keys matching
        those of ``ts``. Used both to compute thermal conductivity at the reference
        depth and volumetric heat capacity for all above-reference layers.
    ref_depth_idx : int, optional
        Index into ``DEPTHS_CM`` identifying the reference depth at which the
        conductive gradient is evaluated. Depths at ``ref_depth_idx - 1`` and
        ``ref_depth_idx + 1`` are used as the bracketing pair for the gradient.
        Defaults to 2, typically corresponding to the 10 cm sensor.
    dt : float, optional
        Timestep duration [s] used to compute dT/dt via finite differencing.
        Defaults to ``DT_SEC``.
    lam_model : {"johansen", "lu2007"}, optional
        Thermal conductivity model to apply at the reference depth.
        ``"johansen"`` uses the Johansen (1975) mixing model;
        ``"lu2007"`` uses the Lu et al. (2007) model.
        Defaults to ``"johansen"``.

    Returns
    -------
    pd.Series
        Estimated surface ground heat flux [W m⁻²], indexed to match ``ts``.
        Named ``"G_grad+stor_{z_ref:.0f}cm"`` where ``z_ref`` is the reference depth
        in cm (e.g., ``"G_grad+stor_10cm"``). Positive values indicate downward flux
        (heat flowing into the soil). The first element will be NaN due to the
        finite-difference storage calculation.

    Notes
    -----
    Conductive flux at the reference depth is computed as:

        G(z_ref) = -λ(θ_ref) · (T_below - T_above) / dz_grad

    where ``dz_grad`` is the distance in meters between the two bracketing sensor
    depths and λ is thermal conductivity [W m⁻¹ K⁻¹] estimated from VWC at the
    reference depth.

    Storage in each above-reference layer is:

        S_i = Cv(θ_i) · (ΔT_i / Δt) · dz_i

    and the total above-reference storage is:

        S(0 → z_ref) = Σ S_i  for i = 0 … ref_depth_idx

    Layer thicknesses ``dz_i`` are taken from the module-level ``LAYER_THICK_M``
    array. NaN values in individual layer storage terms are excluded from the
    column-wise sum via ``np.nansum``.

    References
    ----------
    Johansen, O. (1975). *Thermal conductivity of soils*. PhD thesis, University
    of Trondheim, Norway.

    Lu, S., Ren, T., Gong, Y., & Horton, R. (2007). An improved model for predicting
    soil thermal conductivity from water content at room temperature.
    *Soil Science Society of America Journal*, 71(1), 8–14.
    """
    idx_a = ref_depth_idx - 1
    idx_b = ref_depth_idx + 1
    d_a, d_b = DEPTHS_CM[idx_a], DEPTHS_CM[idx_b]
    dz_grad = (d_b - d_a) / 100.0

    T_a = ts[d_a].values
    T_b = ts[d_b].values
    theta_ref = swc[DEPTHS_CM[ref_depth_idx]].values

    if lam_model == "johansen":
        lam = thermal_conductivity_johansen(theta_ref)
    else:
        lam = thermal_conductivity_lu2007(theta_ref)

    G_ref = -lam * (T_b - T_a) / dz_grad

    # Storage above reference depth
    n_above = ref_depth_idx + 1
    dTdt = ts.diff() / dt
    cv = volumetric_heat_capacity(swc.values)
    S_above = np.nansum(
        cv[:, :n_above] * dTdt.values[:, :n_above] * LAYER_THICK_M[:n_above], axis=1
    )

    G_surface = G_ref + S_above
    return pd.Series(
        G_surface, index=ts.index, name=f"G_grad+stor_{DEPTHS_CM[ref_depth_idx]:.0f}cm"
    )
