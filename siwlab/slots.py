"""Standing-wave longitudinal-slot array synthesis (Elliott's first-order
procedure) using full-wave single-slot characterisation data."""
import numpy as np


def stegen_conductance(x, a, b, lam0, lamg):
    """Stevenson/Stegen normalised resonant conductance of a longitudinal shunt
    slot of offset x in a rectangular waveguide a x b (air-filled form)."""
    return 2.09 * (a / b) * (lamg / lam0) * np.cos(np.pi * lam0 / (2 * lamg)) ** 2 * np.sin(np.pi * x / a) ** 2


def fit_offset_law(xo, g, a):
    """Fit g(x) = K sin^2(pi x / a) in the least-squares sense; returns K."""
    s2 = np.sin(np.pi * np.asarray(xo) / a) ** 2
    return float(np.dot(s2, g) / np.dot(s2, s2))


def offset_for_conductance(g, K, a):
    return a / np.pi * np.arcsin(np.sqrt(np.clip(g / K, 0, 1)))


def resonant_length_from_sweep(L, f_res, f_target):
    """Linear fit of resonant frequency versus slot length; returns L at f_target."""
    L = np.asarray(L); f_res = np.asarray(f_res)
    # f_res ~ c / L  -> fit 1/f = p L + q
    p, q = np.polyfit(L, 1.0 / f_res, 1)
    return float((1.0 / f_target - q) / p)


def aperture_distribution(N, kind='uniform', sll_db=20.0):
    if kind == 'uniform':
        return np.ones(N)
    if kind == 'taylor':
        from scipy.signal.windows import taylor
        return taylor(N, nbar=2, sll=sll_db)
    if kind == 'chebyshev':
        from scipy.signal.windows import chebwin
        return chebwin(N, at=sll_db)
    raise ValueError(kind)


def design_standing_wave_array(N, K, a, L_res, lamg, kind='uniform', sll_db=20.0, y_short=0.0):
    """Return slot list [(yc, xoff, L, w)] (mm) for a resonant array fed from -y,
    shorted at y_short, slots at lamg/2 spacing with alternating offsets.
    Conductances g_n ∝ V_n^2 normalised so that sum g_n = 1 (input match)."""
    V = aperture_distribution(N, kind, sll_db)
    g = V ** 2 / np.sum(V ** 2)
    xo = offset_for_conductance(g, K, a)
    slots = []
    for n in range(N):
        yc = y_short - lamg / 4 - n * lamg / 2
        slots.append((yc, (-1) ** n * xo[n], L_res, None))
    return slots, g, xo


def array_factor(positions, weights, k, theta, phi=np.pi / 2):
    """Array factor for elements at positions (N,2) in the x-y plane."""
    pos = np.asarray(positions)
    u = np.sin(theta) * np.cos(phi); v = np.sin(theta) * np.sin(phi)
    ph = k * (np.outer(u, pos[:, 0]) + np.outer(v, pos[:, 1]))
    return np.exp(1j * ph) @ np.asarray(weights)
