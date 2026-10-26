"""Delivery of plasma-generated long-lived RONS into a layered wound phantom.

The phantom is a 1-D column: a bacterial biofilm (thickness ``d_bio``) resting on a
hydrogel wound bed (thickness ``d_gel``). During treatment the plasma jet imposes a
surface flux of H2O2, NO2-, NO3- and H+ (from HNO2/HNO3). Species diffuse into the
column; H2O2 is consumed in the biofilm by bacterial catalase (first order in viable
density). After treatment the column is held for ``t_hold`` seconds (no flux) before
the electrical measurements, mimicking the post-treatment incubation step.

Concentrations are in mol/m^3 (= mM).
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from scipy.linalg import solve_banded


@dataclass
class PlasmaSource:
    """Surface fluxes of an Ar plasma jet at fixed distance (mol m^-2 s^-1)."""

    flux_h2o2: float = 2.0e-6
    flux_no2: float = 1.1e-6
    flux_no3: float = 0.7e-6
    evaporation: float = 4.0e-7     # water-height loss under the Ar gas flow, m/s
    plasma_on: bool = True

    def scaled(self, factor: float) -> "PlasmaSource":
        return PlasmaSource(self.flux_h2o2 * factor, self.flux_no2 * factor,
                            self.flux_no3 * factor, self.evaporation, self.plasma_on)


@dataclass
class Column:
    d_bio: float = 100e-6           # biofilm thickness, m
    d_gel: float = 3.0e-3           # hydrogel thickness, m
    n_bio: int = 40
    n_gel: int = 120
    # diffusion coefficients in water (m^2/s); hindrance factors per layer
    D_w: dict = field(default_factory=lambda: {"H2O2": 1.43e-9, "NO2": 1.91e-9,
                                               "NO3": 1.90e-9, "H": 9.31e-9})
    hind_bio: float = 0.45
    hind_gel: float = 0.92
    k_cat: float = 0.035            # catalase consumption of H2O2 at full viable density, 1/s

    def grid(self):
        zb = np.linspace(0, self.d_bio, self.n_bio, endpoint=False)
        zg = np.linspace(self.d_bio, self.d_bio + self.d_gel, self.n_gel)
        z = np.concatenate([zb, zg])
        in_bio = z < self.d_bio
        return z, in_bio


SPECIES = ("H2O2", "NO2", "NO3", "H")


def _diffuse_step(c, D_face, dz_c, dz_f, dt, surface_flux, sink):
    """One implicit (backward Euler) step of dc/dt = d/dz(D dc/dz) - sink*c.

    Finite-volume on a non-uniform grid; zero-flux at the bottom, prescribed influx
    at the top. Solved as a tridiagonal (banded) system.
    """
    n = c.size
    a = np.zeros(n)
    b = np.ones(n)
    cc = np.zeros(n)
    g = D_face / dz_f                    # conductance between node i and i+1 (n-1 faces)
    vol = dz_c
    b[:-1] += dt * g / vol[:-1]
    b[1:] += dt * g / vol[1:]
    cc[:-1] = -dt * g / vol[:-1]
    a[1:] = -dt * g / vol[1:]
    b += dt * sink
    rhs = c.copy()
    rhs[0] += dt * surface_flux / vol[0]
    ab = np.vstack([np.concatenate([[0.0], cc[:-1]]), b, np.concatenate([a[1:], [0.0]])])
    return solve_banded((1, 1), ab, rhs)


def simulate_delivery(t_treat: float, source: PlasmaSource, column: Column,
                      viable_density=None, t_hold: float = 600.0, dt: float = 1.0,
                      record_every: float = 5.0):
    """Integrate RONS transport.

    Returns a dict with grid ``z``, mask ``in_bio``, final profiles per species,
    cumulative exposure (time integral of concentration, mM s) per node and a
    time series of layer-averaged concentrations.
    """
    z, in_bio = column.grid()
    n = z.size
    edges = np.concatenate([[0.0], 0.5 * (z[1:] + z[:-1]), [z[-1]]])
    dz_c = np.diff(edges)
    dz_c[dz_c <= 0] = 1e-12
    dz_f = np.diff(z)
    hind = np.where(in_bio, column.hind_bio, column.hind_gel)
    hind_f = 2 * hind[1:] * hind[:-1] / (hind[1:] + hind[:-1])     # harmonic mean
    if viable_density is None:
        viable_density = np.ones(in_bio.sum())
    sink_h2o2 = np.zeros(n)
    sink_h2o2[in_bio] = column.k_cat * viable_density

    c = {s: np.zeros(n) for s in SPECIES}
    exposure = {s: np.zeros(n) for s in SPECIES}
    t_end = t_treat + t_hold
    times, series = [], {s: [] for s in SPECIES}
    t = 0.0
    next_rec = 0.0
    fluxes = {"H2O2": source.flux_h2o2, "NO2": source.flux_no2, "NO3": source.flux_no3,
              "H": source.flux_no2 + source.flux_no3}
    while t < t_end - 1e-9:
        step = min(dt, t_end - t)
        on = source.plasma_on and t < t_treat
        for s in SPECIES:
            D_f = column.D_w[s] * hind_f
            flux = fluxes[s] if on else 0.0
            sink = sink_h2o2 if s == "H2O2" else 0.0
            c[s] = _diffuse_step(c[s], D_f, dz_c, dz_f, step, flux, sink)
            exposure[s] += c[s] * step
        t += step
        if t >= next_rec - 1e-9:
            times.append(t)
            for s in SPECIES:
                series[s].append(np.average(c[s], weights=dz_c))
            next_rec += record_every
    return {"z": z, "in_bio": in_bio, "dz": dz_c, "conc": c, "exposure": exposure,
            "times": np.array(times), "series": {k: np.array(v) for k, v in series.items()}}


def buffered_ph(acid_added_mM, buffer_mM=10.0, pKa=7.2, pH0=7.4):
    """pH after adding strong acid to a phosphate-buffered phantom (vectorised).

    Solves  Cb*(h/(h+Ka) - h0/(h0+Ka)) + (h - h0) = Ca  for h by bisection in log10 h.
    """
    ca = np.atleast_1d(np.asarray(acid_added_mM, float)) * 1e-3      # mol/L
    cb = buffer_mM * 1e-3
    ka = 10 ** -pKa
    h0 = 10 ** -pH0
    lo = np.full_like(ca, -pH0)
    hi = np.full_like(ca, 0.0)

    def resid(logh):
        h = 10 ** logh
        return cb * (h / (h + ka) - h0 / (h0 + ka)) + (h - h0) - ca

    for _ in range(80):
        mid = 0.5 * (lo + hi)
        r = resid(mid)
        hi = np.where(r > 0, mid, hi)
        lo = np.where(r > 0, lo, mid)
    return -0.5 * (lo + hi)


def free_proton_mM(pH, pH0=7.4):
    return (10 ** -np.asarray(pH) - 10 ** -pH0) * 1e3


def hno2_fraction(pH, pKa=3.4):
    return 1.0 / (1.0 + 10 ** (np.asarray(pH) - pKa))
