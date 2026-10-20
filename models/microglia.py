#!/usr/bin/env python3
"""Minimal virtual microglia cell — receptor -> NF-kB -> {TNF-alpha, IL-6}.

A small deterministic ODE model of the canonical microglial pro-inflammatory
axis. A single lumped `receptor_activation` input (0-1), aggregated per polymer
from the Stage-2 ML signal (and Stage-3 docking when available), drives NF-kB
activation, which in turn drives TNF-alpha and IL-6 production. The steady-state
cytokine output is collapsed into a single neuroinflammation index.

====================================================================
ALL KINETIC PARAMETERS BELOW ARE ILLUSTRATIVE (ORDER-OF-MAGNITUDE),
NOT FITTED TO EXPERIMENTAL DATA. They encode qualitative biology
(saturable activation, first-order decay, NF-kB-gated transcription)
sufficient for a RELATIVE ranking across polymers. They must be
refit to microglial time-course data before any absolute claim.
====================================================================

Topology:
    d[NFkB]/dt   = k_on * u * (NFkB_tot - NFkB) - k_off * NFkB
    d[TNF]/dt    = s_tnf * hill(NFkB) - d_tnf * TNF
    d[IL6]/dt    = s_il6 * hill(NFkB) - d_il6 * IL6
where u = receptor_activation in [0,1], hill(x)=x^n/(K^n+x^n).
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.integrate import solve_ivp

# ----- ILLUSTRATIVE parameters (see header warning) -----
PARAMS = {
    "NFkB_tot": 1.0,   # total NF-kB pool (a.u.)          [illustrative]
    "k_on": 2.0,       # activation rate (per u per time)  [illustrative]
    "k_off": 0.8,      # NF-kB deactivation rate           [illustrative]
    "s_tnf": 1.5,      # max TNF-alpha synthesis            [illustrative]
    "d_tnf": 0.9,      # TNF-alpha decay                    [illustrative]
    "s_il6": 1.2,      # max IL-6 synthesis                 [illustrative]
    "d_il6": 0.7,      # IL-6 decay                         [illustrative]
    "K_hill": 0.4,     # NF-kB half-activation for transcription [illustrative]
    "n_hill": 2.0,     # Hill coefficient                   [illustrative]
}


@dataclass
class MicrogliaOutput:
    nfkb_ss: float
    tnf_ss: float
    il6_ss: float
    cytokine_sum: float  # TNF + IL6 at steady state


def _hill(x: float, K: float, n: float) -> float:
    xn = max(x, 0.0) ** n
    return xn / (K ** n + xn)


def _rhs(t, y, u, p):
    nfkb, tnf, il6 = y
    dnfkb = p["k_on"] * u * (p["NFkB_tot"] - nfkb) - p["k_off"] * nfkb
    h = _hill(nfkb, p["K_hill"], p["n_hill"])
    dtnf = p["s_tnf"] * h - p["d_tnf"] * tnf
    dil6 = p["s_il6"] * h - p["d_il6"] * il6
    return [dnfkb, dtnf, dil6]


def simulate(receptor_activation: float, params: dict | None = None,
             t_end: float = 100.0) -> MicrogliaOutput:
    """Integrate to (numerical) steady state for a given receptor drive u in [0,1]."""
    p = dict(PARAMS)
    if params:
        p.update(params)
    u = float(np.clip(receptor_activation, 0.0, 1.0))
    sol = solve_ivp(_rhs, (0.0, t_end), [0.0, 0.0, 0.0], args=(u, p),
                    method="LSODA", rtol=1e-8, atol=1e-10, dense_output=False)
    nfkb, tnf, il6 = sol.y[:, -1]
    return MicrogliaOutput(nfkb_ss=float(nfkb), tnf_ss=float(tnf),
                           il6_ss=float(il6), cytokine_sum=float(tnf + il6))


def max_cytokine_sum(params: dict | None = None) -> float:
    """Steady-state cytokine sum at full drive (u=1), for 0-1 normalisation."""
    return simulate(1.0, params=params).cytokine_sum


if __name__ == "__main__":
    mx = max_cytokine_sum()
    print("Illustrative microglia model self-test (u -> normalised cytokine sum):")
    for u in [0.0, 0.25, 0.5, 0.75, 1.0]:
        out = simulate(u)
        print(f"  u={u:.2f}  NFkB*={out.nfkb_ss:.3f}  TNF*={out.tnf_ss:.3f}  "
              f"IL6*={out.il6_ss:.3f}  index={out.cytokine_sum / mx:.3f}")
