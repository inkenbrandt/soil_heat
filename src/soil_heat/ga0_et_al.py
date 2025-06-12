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

OMEGA_DAY = 2 * np.pi / 86400.0  # rad/s for 24 h period

# -----------------------------------------------------------------------------
# 1. Soil thermal-property relationships (Eq. 12–13)
# -----------------------------------------------------------------------------


def lambda_s(theta: np.ndarray | float) -> np.ndarray | float:
    """Thermal conductivity (λ_s) as a function of volumetric water content θ.

    Implements Eq. (12) from Gao et al. (2017) fileciteturn1file3.
    """
    theta = np.asarray(theta)
    return 0.20 + np.exp(1.46 * (theta - 0.34))


def k_s(theta: np.ndarray | float) -> np.ndarray | float:
    """Thermal diffusivity (k_s) as a function of volumetric water content θ.

    Implements Eq. (13) from Gao et al. (2017) fileciteturn1file6.
    """
    theta = np.asarray(theta)
    return (0.69 + np.exp(3.06 * (theta - 0.26))) * 1e-7


def volumetric_heat_capacity(lambda_s_val, k_s_val):
    """Volumetric heat capacity C_v = λ_s / k_s (units J m⁻³ K⁻¹)."""
    return lambda_s_val / k_s_val


# -----------------------------------------------------------------------------
# 2. Error metrics (Eq. 14–15)
# -----------------------------------------------------------------------------


def nme(calc: np.ndarray, meas: np.ndarray) -> float:
    """Normalized mean error (%). Implements Eq. (14) fileciteturn1file8."""
    return 100.0 * np.sum(np.abs(calc - meas)) / np.sum(np.abs(meas))


def rmse(calc: np.ndarray, meas: np.ndarray) -> float:
    """Root‑mean‑square error. Implements Eq. (15) fileciteturn1file8."""
    return np.sqrt(np.mean((calc - meas) ** 2))


# -----------------------------------------------------------------------------
# 3. Soil‑heat‑flux methods
# -----------------------------------------------------------------------------

# -- 3.1 Calorimetric (Eq. 1) ---------------------------------------------------


def calorimetric_gz(g_zr, cv_layers, dT_dt_layers, dz_layers):
    """Calorimetric method for G_z at depth z (usually 5 cm).

    Parameters
    ----------
    g_zr : float or array_like
        Measured heat flux at reference depth *z_r* (W m⁻²).
    cv_layers : sequence
        Volumetric heat capacity for each sub‑layer *C_v,l* (J m⁻³ K⁻¹).
    dT_dt_layers : sequence
        Time derivative of average temperature for each layer ∂T/∂t (K s⁻¹).
    dz_layers : sequence
        Thickness of each sub‑layer δz_l (m).
    """
    cv_layers = np.asarray(cv_layers)
    dT_dt_layers = np.asarray(dT_dt_layers)
    dz_layers = np.asarray(dz_layers)
    storage = np.sum(cv_layers * dT_dt_layers * dz_layers, axis=0)
    return g_zr + storage  # Eq. (1) fileciteturn1file1


# -- 3.2 Force‑restore (Eq. 2) --------------------------------------------------


def force_restore_gz(cv, dTg_dt, Tg, Tg_bar, delta_z=0.05, omega=OMEGA_DAY):
    """Force‑restore estimate of G_z at z = δz (default 5 cm).

    Implements Eq. (2) fileciteturn1file2.
    """
    term1 = cv * dTg_dt * delta_z
    term2 = np.sqrt(omega * cv * lambda_s_from_cv(cv)) * (
        1 / omega * dTg_dt + (Tg - Tg_bar)
    )
    return term1 + term2


def lambda_s_from_cv(cv):
    """Helper: rough λ_s from C_v and standard soil k_s derived via λ_s/k_s=C_v."""
    # guard: avoid division by zero
    return (
        cv / 2.2e6
    )  # assume typical k_s = λ_s / C_v; here placeholder (~2.2 MJ m⁻³ K⁻¹)


# -- 3.3 Gao 2010 sinusoid (Eq. 3) --------------------------------------------


def gao2010_gz(AT, lambda_s_val, k_s_val, t, omega=OMEGA_DAY):
    """Sinusoidal solution for G_z at depth d (Eq. 3)."""
    d = np.sqrt(2 * k_s_val / omega)
    return np.sqrt(2) * lambda_s_val * AT / d * np.sin(omega * t + np.pi / 4)


# -- 3.4 Heusinkveld harmonic (Eq. 4) -----------------------------------------


def heusinkveld_gz(A_n, Phi_n, n_max, k_s_val, lambda_s_val, w):
    """H04 harmonic solution (Eq. 4)."""
    n = np.arange(1, n_max + 1)
    A_n = np.asarray(A_n)
    Phi_n = np.asarray(Phi_n)
    term = (
        A_n
        * np.sqrt(1 / (k_s_val * n * w * k_s_val))
        * np.sin(n * w + Phi_n + np.pi / 4)
    )
    return (lambda_s_val / (10 * np.pi)) * np.sum(term, axis=0)


# -- 3.5 Hsieh 2009 fractional derivative (Eq. 5) ------------------------------


def hsieh2009_gz(tz_series, time_series, cv_series, ks_series):
    """Half‑order integral solution (Eq. 5).

    tz_series, cv_series, ks_series must be monotonically increasing in *time_series*.
    """
    tz_series = np.asarray(tz_series)
    time_series = np.asarray(time_series)
    cv_series = np.asarray(cv_series)
    ks_series = np.asarray(ks_series)

    t = time_series[-1]
    sum_int = 0.0
    for i in range(len(time_series) - 1):
        si, sip1 = time_series[i], time_series[i + 1]
        delta_tz = tz_series[i + 1] - tz_series[i]
        denom = sip1 - si
        term = delta_tz / denom * (np.sqrt(t - si) - np.sqrt(t - sip1))
        sum_int += term
    ks_t = ks_series[-1]
    cv_t = cv_series[-1]
    return 2 * np.sqrt(ks_t * cv_t / np.pi) * sum_int


# -- 3.6 Leuning 2012 (Eq. 6–7) -------------------------------------------------


def leuning_damping_depth(z, zr, AT_z, AT_zr):
    """Compute damping depth *d* via Eq. (6)."""
    return (zr - z) / np.log(AT_z / AT_zr)


def leuning_gz(g_zr, z, zr, d):
    """Exponentially adjust G_z from reference depth (Eq. 7)."""
    return g_zr * np.exp((zr - z) / d)


# -- 3.7 Simple‑measurement (Eq. 8) -------------------------------------------


def simple_measurement_gz(g_zr, cv_layers, tz_layers, dt, dz_layers):
    """Simple‑measurement variant of calorimetric method (Eq. 8)."""
    tz_layers = np.asarray(tz_layers)
    cv_layers = np.asarray(cv_layers)
    dz_layers = np.asarray(dz_layers)

    delta_tz = tz_layers[:, 1:] - tz_layers[:, :-1]
    delta_tz_mid = 0.5 * (delta_tz[:, :-1] + delta_tz[:, 1:])
    storage = np.sum(
        cv_layers
        * dz_layers[:, None]
        * (tz_layers[:, 1:] - tz_layers[:, :-1] + delta_tz_mid)
        / dt,
        axis=0,
    )
    return g_zr[1:] + storage


# -- 3.8 Wang & Bou‑Zeid 2012 Green's function (Eq. 9–10) ---------------------


def wbz12_g_gz(g_zr_series, time_series, z, zr, k_s_val):
    """WBZ12‑G method (Eq. 9–10)."""
    g_zr_series = np.asarray(g_zr_series)
    t = time_series

    # Pre‑compute ΔF_z
    tau = t - t[0]
    Fz = erfc((zr - z) / (2 * np.sqrt(k_s_val * tau + 1e-12)))
    delta_Fz = np.diff(Fz, prepend=0)

    # Discretised convolution integral J
    J = np.zeros_like(g_zr_series)
    for n in range(1, len(t)):
        g_slice = g_zr_series[n - 1 :: -1]  # reverse order up to n-1
        deltaF_slice = delta_Fz[1 : n + 1]
        J[n] = np.sum(
            (g_slice[:-1] + g_slice[1:]) * deltaF_slice
        )  # trapezoidal approx multiplied by ΔF

    Gz = 2 * g_zr_series - J / delta_Fz[1]
    return Gz


# -- 3.9 Wang & Bou‑Zeid 2012 sinusoid‑plus‑integral (Eq. 11) -----------------


def wbz12_s_gz(Ag, ks_val, zr, z, t, eps, omega=OMEGA_DAY):
    """WBZ12‑S solution (Eq. 11).

    Numerical integration is performed for the second term.
    """
    prefactor = Ag * np.exp(-(zr - z) * np.sqrt(omega / (2 * ks_val)))
    sin_term = np.sin(omega * t + eps - (zr - z) * np.sqrt(omega / (2 * ks_val)))

    def integrand(zeta):
        numerator = ks_val * zeta**2 * (np.sin(eps) - omega * np.cos(eps))
        denom = omega**2 + (ks_val**2) * zeta**4
        return (
            numerator / denom * np.sin(zeta * (zr - z)) * np.exp(-ks_val * zeta**2 * t)
        )

    integral_val = quad(integrand, 0, np.inf, limit=500)[0]
    second_term = -2 * Ag * ks_val / np.pi * integral_val
    return prefactor * sin_term + second_term


# -----------------------------------------------------------------------------
# 4. Exact sinusoidal benchmark (Eq. 16–17)
# -----------------------------------------------------------------------------


def exact_temperature_gz(z, AT, t, d, omega=OMEGA_DAY, T_i=298.15):
    """Exact sinusoidal soil‑temperature profile (Eq. 16)."""
    return T_i + AT * np.exp(-z / d) * np.sin(omega * t - z / d)


def exact_gz(z, AT, lambda_s_val, d, t, omega=OMEGA_DAY):
    """Exact sinusoidal soil‑heat‑flux (Eq. 17)."""
    return (
        np.sqrt(2)
        * lambda_s_val
        * AT
        / d
        * np.exp(-z / d)
        * np.sin(omega * t - z / d + np.pi / 4)
    )
