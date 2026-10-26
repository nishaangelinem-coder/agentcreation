"""Forward model of a wideband (10 Hz - 100 MHz) coplanar-electrode measurement.

A flexible interdigitated electrode (IDE) is pressed onto the biofilm surface. The
IDE sensitivity decays exponentially with distance from the electrode plane
(penetration length ``lam``), so the biofilm and the underlying hydrogel act as
parallel admittance paths weighted by ``w_bio = 1 - exp(-d_bio/lam)``.

    Z(w) = Z_EP(w) + 1 / (Y_sample(w) + j w C_stray)
    Y_sample = j w eps0 [w_bio eps*_bio + (1 - w_bio) eps*_gel] / K
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .constants import EPS0


def frequency_grid(f_min=10.0, f_max=100e6, per_decade=20):
    n = int(round(np.log10(f_max / f_min) * per_decade)) + 1
    return np.logspace(np.log10(f_min), np.log10(f_max), n)


@dataclass
class IDESensor:
    K: float = 600.0              # cell constant, 1/m (from KCl calibration)
    lam: float = 45e-6            # field penetration length, m
    Q_ep: float = 4.0e-7          # electrode-polarisation CPE magnitude, S s^n
    n_ep: float = 0.82
    C_stray: float = 1.5e-12      # residual stray capacitance after open/short compensation, F

    def weight_bio(self, d_bio):
        return 1.0 - np.exp(-d_bio / self.lam)

    def impedance(self, w, eps_bio, eps_gel, d_bio):
        wb = self.weight_bio(d_bio)
        eps_mix = wb * eps_bio + (1 - wb) * eps_gel
        y = 1j * w * EPS0 * eps_mix / self.K + 1j * w * self.C_stray
        z_ep = 1.0 / (self.Q_ep * (1j * w) ** self.n_ep)
        return z_ep + 1.0 / y


def add_noise(z, rng, rel=2.5e-3, abs_ohm=0.02):
    """Instrument noise: proportional (|Z| accuracy) + additive, complex Gaussian."""
    sd = rel * np.abs(z) + abs_ohm
    return z + sd * (rng.standard_normal(z.size) + 1j * rng.standard_normal(z.size)) / np.sqrt(2)


def apparent_permittivity(w, z, K):
    """Convert measured impedance to apparent complex permittivity (eps', sigma)."""
    eps = K / (1j * w * EPS0 * z)
    return eps.real, -eps.imag * w * EPS0
