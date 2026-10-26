"""Electrochemical read-outs of a screen-printed carbon electrode (SPCE) on the phantom.

The electrochemical sub-study uses phantoms doped with an equimolar
[Fe(CN)6]3-/4- redox probe. The biofilm partially blocks the electrode
(fractional coverage ``theta``) which lowers the effective rate constant and the
double-layer capacitance; plasma-derived oxidants convert ferro- to ferricyanide,
shifting the open-circuit potential (Nernst). A Prussian-blue-modified SPCE
operated at -0.05 V gives an amperometric H2O2 read-out.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .constants import F, R, T

f_RT = F / (R * T)


@dataclass
class SPCE:
    area: float = 7.07e-6          # 3 mm carbon disc, m^2
    k0: float = 2.0e-5             # standard rate constant on bare carbon, m/s
    alpha: float = 0.5
    E0: float = 0.22               # formal potential vs Ag/AgCl (pseudo-ref), V
    D_O: float = 7.2e-10 * 0.9     # ferricyanide in hydrogel, m^2/s
    D_R: float = 6.7e-10 * 0.9
    c_total: float = 5.0           # [Fe(CN)6]3- + [Fe(CN)6]4-, mol/m^3
    C_dl_area: float = 0.15        # F/m^2 (15 uF/cm^2)
    n_dl: float = 0.90
    Ru: float = 180.0              # uncompensated phantom resistance, ohm
    pb_sensitivity: float = 0.55   # Prussian-blue H2O2 sensitivity, A M^-1 cm^-2


def effective_params(e: SPCE, theta):
    return e.k0 * (1 - theta), e.C_dl_area * e.area * (1 - 0.7 * theta)


def randles_spectrum(f, e: SPCE, theta, x_ox, Ru=None):
    """Faradaic EIS at open circuit for a probe with ferri fraction ``x_ox``."""
    w = 2 * np.pi * f
    Ru = e.Ru if Ru is None else Ru
    k_eff, Cdl = effective_params(e, theta)
    cO = e.c_total * x_ox
    cR = e.c_total * (1 - x_ox)
    i0 = F * e.area * k_eff * cO ** (1 - e.alpha) * cR ** e.alpha
    Rct = R * T / (F * i0)
    sw = R * T / (F ** 2 * e.area * np.sqrt(2)) * (1 / (np.sqrt(e.D_O) * cO) + 1 / (np.sqrt(e.D_R) * cR))
    zw = sw * (1 - 1j) / np.sqrt(w)
    Q = Cdl * (2 * np.pi * 1e3) ** (1 - e.n_dl)          # CPE with C_dl at 1 kHz
    z = Ru + 1.0 / (Q * (1j * w) ** e.n_dl + 1.0 / (Rct + zw))
    return z, {"Rct_true": Rct, "sigma_w_true": sw, "Cdl_true": Cdl}


def ocp(e: SPCE, x_ox):
    return e.E0 + np.log(x_ox / (1 - x_ox)) / f_RT


def cyclic_voltammogram(e: SPCE, theta, x_ox, E_start=-0.2, E_vertex=0.6, scan=0.05,
                        dx=2.5e-6, n_x=500, n_cycles=2):
    """Explicit finite-difference CV with Butler-Volmer kinetics and iR drop.

    ``n_cycles`` scans are simulated from the start potential and the last one is
    returned (quasi-steady voltammogram, as reported experimentally).
    """
    k_eff, Cdl = effective_params(e, theta)
    D = max(e.D_O, e.D_R)
    dt = 0.4 * dx ** 2 / D
    t_half = (E_vertex - E_start) / scan
    n_cyc = int(2 * t_half / dt)
    nt = n_cycles * n_cyc
    t = np.arange(nt) * dt
    tc = t % (2 * t_half)
    E_app = np.where(tc < t_half, E_start + scan * tc, E_vertex - scan * (tc - t_half))
    cO = np.full(n_x, e.c_total * x_ox)
    cR = np.full(n_x, e.c_total * (1 - x_ox))
    lO = e.D_O * dt / dx ** 2
    lR = e.D_R * dt / dx ** 2
    cur = np.zeros(nt)
    i_prev = 0.0
    dEdt = np.where(tc < t_half, scan, -scan)
    for k in range(nt):
        Eeff = E_app[k] - i_prev * e.Ru
        eta = Eeff - e.E0
        kf = k_eff * np.exp(-e.alpha * f_RT * eta)          # reduction O + e -> R
        kb = k_eff * np.exp((1 - e.alpha) * f_RT * eta)     # oxidation
        # surface flux balance with first interior node (Euler half-cell)
        cO1, cR1 = cO[1], cR[1]
        a = e.D_O / dx
        b = e.D_R / dx
        # solve: a(cO1-cO0) = kf cO0 - kb cR0 ;  b(cR1-cR0) = -(kf cO0 - kb cR0)
        det = (a + kf) * (b + kb) - kf * kb
        cO0 = (a * cO1 * (b + kb) + kb * b * cR1) / det
        cR0 = (b * cR1 * (a + kf) + kf * a * cO1) / det
        cO[0], cR[0] = cO0, cR0
        flux = kf * cO0 - kb * cR0                           # net reduction flux, mol m^-2 s^-1
        i_far = -F * e.area * flux                           # anodic positive
        i_cap = Cdl * dEdt[k]
        i_prev = i_far + i_cap
        cur[k] = i_prev
        cO[1:-1] += lO * (cO[2:] - 2 * cO[1:-1] + cO[:-2])
        cR[1:-1] += lR * (cR[2:] - 2 * cR[1:-1] + cR[:-2])
    last = slice(nt - n_cyc, nt)
    stride = max(1, n_cyc // 2000)
    return E_app[last][::stride], cur[last][::stride]


def cv_metrics(E, i):
    half = np.argmax(E)
    Ea = E[:half]; ia = i[:half]
    Ec = E[half:]; ic = i[half:]
    ka = np.argmax(ia); kc = np.argmin(ic)
    return {"Epa": Ea[ka], "ipa": ia[ka], "Epc": Ec[kc], "ipc": ic[kc],
            "dEp": Ea[ka] - Ec[kc], "E_half": 0.5 * (Ea[ka] + Ec[kc])}


def amperometric_h2o2(e: SPCE, c_h2o2_mM, rng, noise_nA=4.0, baseline_nA=15.0):
    """Steady-state cathodic current (nA) of the PB-modified SPCE."""
    area_cm2 = e.area * 1e4
    i = e.pb_sensitivity * area_cm2 * (np.asarray(c_h2o2_mM) * 1e-3) * 1e9
    return i + baseline_nA + noise_nA * rng.standard_normal(np.shape(c_h2o2_mM))
