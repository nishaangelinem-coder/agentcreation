"""Complex permittivity of the wound-phantom constituents.

* hydrogel / wound-fluid simulant: Debye water + dc conductivity
* bacterium: single-shell sphere (cytoplasm + membrane), Pauly-Schwan / Asami
* biofilm: Maxwell-Garnett mixture of intact and permeabilised cells in an EPS matrix

All functions take angular frequency ``w`` (rad/s) and return complex relative
permittivity eps* = eps' - j sigma/(w eps0) (engineering sign convention, e^{jwt}).
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .constants import EPS0


def debye_medium(w, sigma, eps_s=78.4, eps_inf=5.2, tau=8.27e-12):
    return eps_inf + (eps_s - eps_inf) / (1 + 1j * w * tau) - 1j * sigma / (w * EPS0)


@dataclass
class Bacterium:
    radius: float = 0.45e-6        # S. aureus coccus radius, m
    d_mem: float = 6e-9            # envelope (membrane) thickness, m
    eps_mem: float = 6.0
    sigma_mem: float = 1e-7        # intact membrane, S/m
    eps_cyt: float = 60.0
    sigma_cyt: float = 0.55        # cytoplasm, S/m

    def permittivity(self, w):
        e_c = self.eps_cyt - 1j * self.sigma_cyt / (w * EPS0)
        e_m = self.eps_mem - 1j * self.sigma_mem / (w * EPS0)
        v = (1 - self.d_mem / self.radius) ** 3
        return e_m * (2 * (1 - v) * e_m + (1 + 2 * v) * e_c) / ((2 + v) * e_m + (1 - v) * e_c)

    def permeabilised(self, sigma_medium, shrink=0.9):
        """Envelope-damaged cell: leaky membrane, cytoplasm equilibrating with the medium."""
        return Bacterium(radius=self.radius * shrink, d_mem=self.d_mem, eps_mem=self.eps_mem,
                         sigma_mem=2e-3, eps_cyt=self.eps_cyt,
                         sigma_cyt=0.5 * (self.sigma_cyt + sigma_medium))


def maxwell_garnett(e_host, inclusions):
    """Generalised Maxwell-Garnett: ``inclusions`` = [(volume_fraction, eps*), ...]."""
    s = sum(p * (e - e_host) / (e + 2 * e_host) for p, e in inclusions)
    return e_host * (1 + 2 * s) / (1 - s)


def biofilm_permittivity(w, sigma_eps, phi_cells, frac_perm, cell: Bacterium, shrink=0.9):
    """EPS matrix (conductivity ``sigma_eps``) loaded with live + permeabilised cells."""
    host = debye_medium(w, sigma_eps, eps_s=74.0)
    dead = cell.permeabilised(sigma_eps, shrink)
    phi_live = phi_cells * (1 - frac_perm)
    phi_dead = phi_cells * frac_perm * shrink ** 3
    return maxwell_garnett(host, [(phi_live, cell.permittivity(w)),
                                  (phi_dead, dead.permittivity(w))])
