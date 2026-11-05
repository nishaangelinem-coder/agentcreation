"""SIW design equations (all lengths in metres unless stated)."""
import numpy as np
from .constants import c0, eta0, eps0


def effective_width(a, d, p):
    """Cassivi / Xu-Wu equivalent rectangular-waveguide width of an SIW with
    via-row spacing a, via diameter d and pitch p.
    a_eff = a - d^2/(0.95 p)            (Cassivi et al. 2002)
    """
    return a - d ** 2 / (0.95 * p)


def effective_width_xu(a, d, p):
    """Xu & Wu (2005) refined expression."""
    return a - 1.08 * d ** 2 / p + 0.1 * d ** 2 / a


def cutoff_frequency(a_eff, er):
    return c0 / (2 * a_eff * np.sqrt(er))


def beta(f, a_eff, er):
    """TE10 phase constant (rad/m); 0 below cutoff."""
    k = 2 * np.pi * f / c0 * np.sqrt(er)
    kc = np.pi / a_eff
    return np.sqrt(np.maximum(k ** 2 - kc ** 2, 0.0))


def guided_wavelength(f, a_eff, er):
    return 2 * np.pi / beta(f, a_eff, er)


def wave_impedance_TE(f, a_eff, er):
    k = 2 * np.pi * f / c0 * np.sqrt(er)
    return (eta0 / np.sqrt(er)) * k / beta(f, a_eff, er)


def dielectric_attenuation(f, a_eff, er, tand):
    """Dielectric attenuation (Np/m) of TE10: alpha_d = k^2 tan(d) / (2 beta)."""
    k = 2 * np.pi * f / c0 * np.sqrt(er)
    return k ** 2 * tand / (2 * beta(f, a_eff, er))


def conductor_attenuation(f, a_eff, h, er, sigma_c=5.8e7):
    """Conductor attenuation (Np/m) of TE10 in an equivalent rectangular
    waveguide (Pozar), width a_eff, height h."""
    Rs = np.sqrt(2 * np.pi * f * 4e-7 * np.pi / (2 * sigma_c))
    k = 2 * np.pi * f / c0 * np.sqrt(er)
    b = beta(f, a_eff, er)
    eta = eta0 / np.sqrt(er)
    return Rs / (a_eff ** 3 * h * b * k * eta) * (2 * h * np.pi ** 2 + a_eff ** 3 * k ** 2)


def leakage_check(d, p, lam_c):
    """Design-rule checks (Deslandes & Wu 2006, Xu & Wu 2005)."""
    return dict(pitch_le_2d=p <= 2 * d, d_lt_lamc_over_5=d < lam_c / 5, p_lt_lamc_over_2=p < lam_c / 2)
