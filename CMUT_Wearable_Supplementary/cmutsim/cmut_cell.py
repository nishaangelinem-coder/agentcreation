"""Lumped electromechanical model of a clamped circular CMUT cell and the
Butterworth-Van Dyke (BVD) model of an array element."""
import numpy as np
from scipy.optimize import brentq
from . import params as P


def flexural_rigidity(t=P.PLATE_T, E=P.E_PLATE, nu=P.NU_PLATE):
    return E * t**3 / (12 * (1 - nu**2))


def spring_constant(a=P.CELL_RADIUS, t=P.PLATE_T, E=P.E_PLATE, nu=P.NU_PLATE):
    """Piston-equivalent stiffness of a clamped plate, k = 64*pi*D/a^2."""
    return 64 * np.pi * flexural_rigidity(t, E, nu) / a**2


def effective_gap(g0=P.GAP, t_ins=P.T_INSULATOR, eps_ins=P.EPS_INS):
    return g0 + t_ins / eps_ins


def collapse_voltage(a=P.CELL_RADIUS, t=P.PLATE_T, g0=P.GAP, **kw):
    k = spring_constant(a, t)
    g = effective_gap(g0)
    A = np.pi * a**2
    return np.sqrt(8 * k * g**3 / (27 * P.EPS0 * A))


def clamped_capacitance(a=P.CELL_RADIUS, g0=P.GAP):
    """Vacuum-gap capacitance eps0*A/g0 (the definition used in Table I)."""
    return P.EPS0 * np.pi * a**2 / g0


def static_deflection(V, a=P.CELL_RADIUS, t=P.PLATE_T, g0=P.GAP):
    """Stable piston deflection x (m) under DC bias V (< Vc)."""
    k = spring_constant(a, t)
    g = effective_gap(g0)
    A = np.pi * a**2
    f = lambda x: k * x - P.EPS0 * A * V**2 / (2 * (g - x)**2)
    if V <= 0:
        return 0.0
    return brentq(f, 0.0, g / 3 * 0.99999)


def resonance_vacuum(a=P.CELL_RADIUS, t=P.PLATE_T):
    """First axisymmetric mode of a clamped plate (Hz)."""
    D = flexural_rigidity(t)
    return 10.21 / (2 * np.pi * a**2) * np.sqrt(D / (P.RHO_PLATE * t))


def resonance_immersed(a=P.CELL_RADIUS, t=P.PLATE_T, V=0.0, g0=P.GAP,
                       rho_f=P.RHO0):
    """Fluid-loaded resonance with added-mass factor and spring softening."""
    beta = 0.6689 * (rho_f / P.RHO_PLATE) * (a / t)
    f_vac = resonance_vacuum(a, t)
    if V > 0:
        k = spring_constant(a, t)
        g = effective_gap(g0)
        x = static_deflection(V, a, t, g0)
        k_soft = k - P.EPS0 * np.pi * a**2 * V**2 / (g - x)**3
        f_vac = f_vac * np.sqrt(max(k_soft, 0) / k)
    return f_vac / np.sqrt(1 + beta)


def radius_for_immersed_resonance(f_target, t=P.PLATE_T, V_frac=0.8):
    """Plate radius that places the biased immersed resonance at f_target."""
    def fn(a):
        Vc = collapse_voltage(a, t)
        return resonance_immersed(a, t, V=V_frac * Vc) - f_target
    return brentq(fn, 5e-6, 200e-6)


def cells_per_element(width=P.WIDTH, height=P.ELEV_HEIGHT, pitch=P.CELL_PITCH):
    return int(width // pitch) * int(height // pitch)


def bvd_element(f, n_cells=None, f0=P.F0, Q=P.Q_LOADED, k2=0.5, a=P.CELL_RADIUS,
                g0=P.GAP):
    """Electrical input impedance and normalised transmit response of one
    element (BVD model). k2 = effective electromechanical coupling."""
    if n_cells is None:
        n_cells = cells_per_element()
    C0 = n_cells * clamped_capacitance(a, g0)
    Cm = C0 * k2 / (1 - k2)
    Lm = 1 / ((2 * np.pi * f0)**2 * Cm)
    Rm = 2 * np.pi * f0 * Lm / Q
    w = 2 * np.pi * f
    Zm = Rm + 1j * w * Lm + 1 / (1j * w * Cm)
    Z0 = 1 / (1j * w * C0)
    Zin = Z0 * Zm / (Z0 + Zm)
    Htx = np.abs(1 / Zm)            # motional current ~ velocity ~ pressure
    Htx = Htx / Htx.max()
    return Zin, Htx, dict(C0=C0, Cm=Cm, Lm=Lm, Rm=Rm, n_cells=n_cells)


def fractional_bandwidth(f, H, level_db=-6):
    Hdb = 20 * np.log10(H / H.max())
    idx = np.where(Hdb >= level_db)[0]
    f_lo, f_hi = f[idx[0]], f[idx[-1]]
    fc = 0.5 * (f_lo + f_hi)
    return (f_hi - f_lo) / fc, fc


def surface_pressure_per_volt(V_dc, a=P.CELL_RADIUS, g0=P.GAP,
                              pitch=P.CELL_PITCH):
    """Radiation-dominated (fluid-loaded) small-signal surface pressure per
    volt of AC drive: p/Vac = eps0*Vdc/(g_b^2) * fill_factor."""
    x = static_deflection(V_dc, a, g0=g0)
    g_b = effective_gap(g0) - x
    fill = np.pi * a**2 / pitch**2
    return P.EPS0 * V_dc / g_b**2 * fill
