"""Virtual experiment: replicate infected wound phantoms -> CAP -> measurements -> fits."""

from __future__ import annotations

from dataclasses import dataclass, replace

import numpy as np
import pandas as pd

from . import electrochem as ec
from .constants import LAMBDA
from .dielectric import Bacterium, biofilm_permittivity, debye_medium
from .fitting import fit_dielectric, fit_randles, lin_kk
from .impedance import IDESensor, add_noise, frequency_grid
from .microbiology import WeibullInactivation, biofilm_response
from .plasma_chemistry import Column, PlasmaSource, buffered_ph, free_proton_mM, simulate_delivery


@dataclass(frozen=True)
class Group:
    name: str
    t_treat: float
    plasma: bool = True
    infected: bool = True


GROUPS = (
    Group("Control", 0.0),
    Group("CAP-30", 30.0),
    Group("CAP-60", 60.0),
    Group("CAP-120", 120.0),
    Group("CAP-180", 180.0),
    Group("CAP-300", 300.0),
    Group("Gas-300", 300.0, plasma=False),
    Group("Sterile", 0.0, infected=False),
    Group("Sterile+CAP-300", 300.0, infected=False),
)

ETCH_RATE = 0.04e-6          # biofilm thinning under the plasma jet, m/s
THETA0 = 0.62                # SPCE blocking by an untreated biofilm
THETA_STERILE = 0.03         # conditioning-film blocking on sterile phantoms
GEL_SENSE_DEPTH = 200e-6     # depth of gel that dominates the IDE reading


def sample_replicate(rng, infected=True):
    """Draw phantom-to-phantom variability (one replicate)."""
    return {
        "log10_N0": rng.normal(8.0, 0.15) if infected else -np.inf,
        "d_bio": rng.normal(100e-6, 8e-6) if infected else 100e-6,
        "phi": float(np.clip(rng.normal(0.30, 0.02), 0.2, 0.4)) if infected else 0.0,
        "sigma_gel": rng.normal(1.00, 0.03),
        "flux_factor": float(np.exp(rng.normal(0, 0.10))),
        "delta_factor": float(np.exp(rng.normal(0, 0.12))),
        "Q_ep": rng.normal(4.0e-7, 0.2e-7),
        "n_ep": rng.normal(0.82, 0.01),
        "K_true": rng.normal(600.0, 6.0),
        "k0_factor": float(np.exp(rng.normal(0, 0.06))),
    }


def delta_sigma(c_no2, c_no3, pH):
    """Conductivity change (S/m) from plasma-derived ions in a phosphate buffer."""
    h = np.clip(free_proton_mM(pH), 0, None)
    bound = np.clip(c_no2 + c_no3 - h, 0, None)
    d_lambda_buffer = 2 * 57.0e-4 - 36.0e-4    # HPO4^2- -> H2PO4- on protonation
    return LAMBDA["NO2-"] * c_no2 + LAMBDA["NO3-"] * c_no3 + LAMBDA["H+"] * h - d_lambda_buffer * bound


def run_sample(group: Group, rep: dict, rng, freqs, eis_freqs, do_cv=True, noise=True):
    infected = group.infected
    col = Column(d_bio=rep["d_bio"])
    if not infected:
        col = replace(col, hind_bio=col.hind_gel, k_cat=0.0)
    src = PlasmaSource(plasma_on=group.plasma).scaled(rep["flux_factor"])
    src.plasma_on = group.plasma
    deliv = simulate_delivery(group.t_treat, src, col)
    acid = deliv["conc"]["H"]
    pH = buffered_ph(acid)
    m = deliv["in_bio"]
    dz = deliv["dz"]

    model = WeibullInactivation()
    model = replace(model, delta=model.delta * rep["delta_factor"])
    t_gas = group.t_treat
    if infected:
        bio = biofilm_response(deliv, pH, model, t_gas=t_gas, log10_N0=rep["log10_N0"])
    else:
        bio = {"viable_fraction": np.nan, "permeabilised_fraction": 0.0,
               "log_reduction": np.nan, "log10_cfu": np.nan}

    # --- ionic changes (top biofilm layer and sensed gel depth)
    zg = deliv["z"] - col.d_bio
    gsel = (~m) & (zg < GEL_SENSE_DEPTH)
    cB = {s: np.average(deliv["conc"][s][m], weights=dz[m]) for s in ("H2O2", "NO2", "NO3")}
    cG = {s: np.average(deliv["conc"][s][gsel], weights=dz[gsel]) for s in ("H2O2", "NO2", "NO3")}
    pH_bio = float(np.average(pH[m], weights=dz[m]))
    pH_gel = float(np.average(pH[gsel], weights=dz[gsel]))
    evap = 1.0 + src.evaporation * t_gas / (col.d_gel + col.d_bio)
    sig_gel = rep["sigma_gel"] * evap + delta_sigma(cG["NO2"], cG["NO3"], pH_gel)
    sig_eps = 0.85 * rep["sigma_gel"] * evap + delta_sigma(cB["NO2"], cB["NO3"], pH_bio)
    if not infected:
        sig_eps = sig_gel

    d_bio_eff = rep["d_bio"] - (ETCH_RATE * group.t_treat if (group.plasma and infected) else 0.0)

    # --- wideband impedance
    w = 2 * np.pi * freqs
    cell = Bacterium()
    eps_bio = biofilm_permittivity(w, sig_eps, rep["phi"], bio["permeabilised_fraction"], cell)
    eps_gel = debye_medium(w, sig_gel)
    sensor = IDESensor(K=rep["K_true"], Q_ep=rep["Q_ep"], n_ep=rep["n_ep"])
    z_true = sensor.impedance(w, eps_bio, eps_gel, d_bio_eff)
    z = add_noise(z_true, rng) if noise else z_true
    fit = fit_dielectric(freqs, z, K=600.0)
    M, mu, rr, ri = lin_kk(freqs, z)

    # --- electrochemistry (redox-probe phantom, same treatment)
    spce = ec.SPCE()
    spce = replace(spce, k0=spce.k0 * rep["k0_factor"], Ru=spce.Ru / (sig_gel / 1.0))
    if infected:
        theta = THETA0 * (d_bio_eff / 100e-6) * (1 - 0.30 * bio["permeabilised_fraction"])
    else:
        theta = THETA_STERILE
    theta = float(np.clip(theta, 0, 0.95))
    z_eis, truth = ec.randles_spectrum(eis_freqs, spce, theta, 0.5)
    z_eis = add_noise(z_eis, rng, rel=4e-3, abs_ohm=0.05) if noise else z_eis
    rf = fit_randles(eis_freqs, z_eis)
    cvm = {}
    E = I = None
    if do_cv:
        E, I = ec.cyclic_voltammogram(spce, theta, 0.5)
        if noise:
            I = I + 20e-9 * rng.standard_normal(I.size)
        cvm = ec.cv_metrics(E, I)
    i_h2o2 = float(ec.amperometric_h2o2(spce, cB["H2O2"], rng))
    # IrOx pH micro-sensor (-58 mV/pH) on the wound surface
    e_ph = 0.70 - 0.058 * pH_bio + 0.0015 * rng.standard_normal()
    pH_meas = (0.70 - e_ph) / 0.058

    row = {
        "group": group.name, "t_treat": group.t_treat, "plasma": group.plasma,
        "infected": infected,
        # ground truth
        "log10_N0": rep["log10_N0"], "log10_cfu": bio["log10_cfu"],
        "log_reduction": bio["log_reduction"], "perm_frac": bio["permeabilised_fraction"],
        "d_bio_eff_um": d_bio_eff * 1e6,
        "H2O2_bio_mM": cB["H2O2"], "NO2_bio_mM": cB["NO2"], "NO3_bio_mM": cB["NO3"],
        "H2O2_gel_mM": cG["H2O2"], "NO2_gel_mM": cG["NO2"], "NO3_gel_mM": cG["NO3"],
        "pH_bio": pH_bio, "pH_gel": pH_gel, "sigma_gel_true": sig_gel, "theta": theta,
        "Rct_true": truth["Rct_true"],
        # bioimpedance features
        **{k: v for k, v in fit.items() if k != "z_fit"},
        "absZ_1k": np.abs(z[np.argmin(abs(freqs - 1e3))]),
        "absZ_100k": np.abs(z[np.argmin(abs(freqs - 1e5))]),
        "absZ_10M": np.abs(z[np.argmin(abs(freqs - 1e7))]),
        "phase_5M": np.degrees(np.angle(z[np.argmin(abs(freqs - 5e6))])),
        "kk_M": M, "kk_mu": mu, "kk_rms": float(np.sqrt(np.mean(rr ** 2 + ri ** 2))),
        "kk_res_max": float(np.max(np.abs(np.r_[rr, ri]))),
        # electrochemical features
        "Rs_fit": rf["Rs"], "Q_dl_fit": rf["Q_dl"], "n_dl_fit": rf["n_dl"],
        "Rct_fit": rf["Rct"], "sigma_w_fit": rf["sigma_w"],
        **{f"cv_{k}": v for k, v in cvm.items()},
        "i_h2o2_nA": i_h2o2, "pH_meas": pH_meas,
    }
    traces = {"z": z, "z_fit": fit["z_fit"], "z_eis": z_eis, "z_eis_fit": rf["z_fit"],
              "E": E, "I": I, "deliv": deliv, "pH_profile": pH,
              "bio": bio, "kk": (rr, ri)}
    return row, traces


def run_experiment(n_rep=6, seed=2026, groups=GROUPS):
    rng = np.random.default_rng(seed)
    freqs = frequency_grid()
    eis_freqs = np.logspace(-1, 5, 61)
    rows, traces = [], {}
    for g in groups:
        for r in range(n_rep):
            rep = sample_replicate(rng, g.infected)
            row, tr = run_sample(g, rep, rng, freqs, eis_freqs, do_cv=True)
            row["replicate"] = r + 1
            rows.append(row)
            traces[(g.name, r)] = tr
    return pd.DataFrame(rows), traces, freqs, eis_freqs
