"""Run the complete in silico study: simulation -> fitting -> statistics -> figures.

Usage:  python scripts/run_all.py [--n-rep 6] [--seed 2026]
Outputs go to results/ (CSV/JSON/TeX tables) and figures/ (PDF + 600-dpi PNG).
"""

from __future__ import annotations

import argparse
import json
import sys
from dataclasses import replace
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats
from scipy.optimize import curve_fit

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from capwound import pipeline as pl                                      # noqa: E402
from capwound.dielectric import Bacterium, biofilm_permittivity, debye_medium  # noqa: E402
from capwound.fitting import fit_dielectric                              # noqa: E402
from capwound.impedance import IDESensor, apparent_permittivity          # noqa: E402
from capwound.plasma_chemistry import Column, PlasmaSource, buffered_ph, simulate_delivery  # noqa: E402
from capwound.stats import anova_tukey, linear_lod, nested_pls, spearman_table, vip_scores  # noqa: E402

RES = ROOT / "results"
FIG = ROOT / "figures"
RES.mkdir(exist_ok=True)
FIG.mkdir(exist_ok=True)

# ----------------------------------------------------------------- style
plt.rcParams.update({
    "font.family": "DejaVu Sans", "font.size": 7.5, "axes.titlesize": 8,
    "axes.labelsize": 7.5, "legend.fontsize": 6.5, "xtick.labelsize": 6.5,
    "ytick.labelsize": 6.5, "axes.linewidth": 0.6, "lines.linewidth": 1.4,
    "axes.spines.top": False, "axes.spines.right": False, "axes.grid": True,
    "grid.color": "#e4e3df", "grid.linewidth": 0.5, "savefig.dpi": 600,
    "axes.edgecolor": "#52514e", "xtick.color": "#52514e", "ytick.color": "#52514e",
    "axes.labelcolor": "#0b0b0b", "legend.frameon": False,
})
SEQ = {"CAP-30": "#86b6ef", "CAP-60": "#5598e7", "CAP-120": "#2a78d6",
       "CAP-180": "#1c5cab", "CAP-300": "#0d366b"}
COLORS = {"Control": "#8a8984", **SEQ, "Gas-300": "#eb6834",
          "Sterile": "#1baf7a", "Sterile+CAP-300": "#0f7a54"}
ORDER = [g.name for g in pl.GROUPS]
INFECTED = ["Control", "CAP-30", "CAP-60", "CAP-120", "CAP-180", "CAP-300", "Gas-300"]
TIMES = [0, 30, 60, 120, 180, 300]
W2, W1 = 7.0, 3.35            # Springer double / single column widths (in)


def save(fig, name):
    fig.savefig(FIG / f"{name}.pdf", bbox_inches="tight")
    fig.savefig(FIG / f"{name}.png", bbox_inches="tight")
    plt.close(fig)


def panel(ax, letter):
    ax.text(-0.26, 1.04, letter, transform=ax.transAxes, fontsize=9, fontweight="bold", va="bottom")


def msd(df, col, groups=ORDER):
    g = df.groupby("group")[col]
    return g.mean().reindex(groups), g.std().reindex(groups)


# ----------------------------------------------------------------- figures
def fig_schematic():
    fig, ax = plt.subplots(figsize=(W2, 2.3))
    ax.set_xlim(0, 10); ax.set_ylim(0, 4); ax.axis("off"); ax.grid(False)
    ax.add_patch(plt.Rectangle((0.6, 0.5), 4.2, 1.2, fc="#cde2fb", ec="#256abf", lw=0.8))
    ax.add_patch(plt.Rectangle((0.6, 1.7), 4.2, 0.35, fc="#f6d3a8", ec="#c98500", lw=0.8))
    rng = np.random.default_rng(1)
    xs = rng.uniform(0.7, 4.7, 160); ys = rng.uniform(1.73, 2.02, 160)
    ax.scatter(xs, ys, s=3, c="#8a5a00", lw=0)
    ax.text(2.7, 1.05, "Agarose / simulated wound fluid (3 mm)\n10 mM phosphate, σ = 1.0 S/m",
            ha="center", va="center", fontsize=6.5)
    ax.text(0.65, 2.42, "S. aureus biofilm\n(100 µm, φ = 0.30)", va="center", fontsize=6.3)
    ax.add_patch(plt.Polygon([[2.4, 3.8], [3.0, 3.8], [2.85, 2.9], [2.55, 2.9]], fc="#d9d8d3", ec="#52514e", lw=0.8))
    ax.add_patch(plt.Polygon([[2.58, 2.9], [2.82, 2.9], [2.95, 2.1], [2.45, 2.1]], fc="#b9a6f0", ec="none", alpha=0.8))
    ax.text(3.1, 3.45, "Ar plasma jet\n(H2O2, NO2-, NO3-, H+ fluxes)", fontsize=6.5, va="center")
    for i, (t, x) in enumerate([("RONS transport\n(1-D diffusion-reaction)", 6.6),
                                ("Weibull inactivation\n+ permeabilisation", 6.6),
                                ("IDE 10 Hz-100 MHz\nbioimpedance", 8.6),
                                ("SPCE: EIS, CV,\nH2O2 & pH sensors", 8.6)]):
        y = 3.2 if i % 2 == 0 else 2.0
        ax.add_patch(matplotlib.patches.FancyBboxPatch((x - 0.85, y - 0.4), 1.7, 0.8,
                     boxstyle="round,pad=0.03", fc="#fcfcfb", ec="#52514e", lw=0.7))
        ax.text(x, y, t, ha="center", va="center", fontsize=6.3)
    ax.annotate("", (5.75, 3.2), (3.2, 3.0), arrowprops=dict(arrowstyle="->", lw=0.7, color="#52514e"))
    ax.annotate("", (6.6, 2.4), (6.6, 2.8), arrowprops=dict(arrowstyle="->", lw=0.7))
    ax.annotate("", (7.75, 2.0), (7.45, 2.0), arrowprops=dict(arrowstyle="->", lw=0.7))
    ax.annotate("", (8.6, 2.8), (8.6, 2.4), arrowprops=dict(arrowstyle="<-", lw=0.7))
    ax.text(7.6, 0.75, "CNLS (Cole-Cole + CPE, Randles), Lin-KK,\nANOVA/Tukey, nested PLS regression",
            ha="center", fontsize=6.5, style="italic")
    save(fig, "fig1_schematic")


def fig_chemistry(df):
    fig, axs = plt.subplots(1, 3, figsize=(W2, 2.2))
    src, col = PlasmaSource(), Column()
    d = simulate_delivery(300, src, col, record_every=5)
    ax = axs[0]
    for s, lab, c in (("H2O2", "H$_2$O$_2$", "#2a78d6"), ("NO2", "NO$_2^-$", "#eb6834"), ("NO3", "NO$_3^-$", "#1baf7a")):
        ax.plot(d["times"] / 60, d["series"][s] * 1e3, label=lab, color=c)
    ax.axvline(5, color="#8a8984", lw=0.6, ls="--"); ax.text(5.15, ax.get_ylim()[1] * 0.9, "jet off", fontsize=6)
    ax.set_xlabel("Time (min)"); ax.set_ylabel("Column-mean conc. (µM)"); ax.legend(loc="lower right")
    panel(ax, "a")
    ax = axs[1]
    for t in (30, 120, 300):
        dd = simulate_delivery(t, src, col, t_hold=600)
        ax.plot(dd["z"] * 1e3, dd["conc"]["H2O2"] * 1e3, color=SEQ[f"CAP-{t}"], label=f"{t} s")
        dd0 = simulate_delivery(t, src, replace(col, hind_bio=col.hind_gel, k_cat=0.0), t_hold=600)
        ax.plot(dd0["z"] * 1e3, dd0["conc"]["H2O2"] * 1e3, color=SEQ[f"CAP-{t}"], ls=":", lw=1.0)
    ax.axvspan(0, 0.1, color="#f6d3a8", alpha=0.6, lw=0)
    ax.set_xlabel("Depth z (mm)"); ax.set_ylabel("H$_2$O$_2$ at measurement (µM)")
    ax.legend(title="CAP (solid: infected,\ndotted: sterile)", title_fontsize=6)
    panel(ax, "b")
    ax = axs[2]
    m, s = msd(df, "i_h2o2_nA", INFECTED[:6])
    t = np.array(TIMES)
    ax.errorbar(t, m, s, fmt="o-", color="#2a78d6", ms=4, capsize=2, label="infected")
    ms_, ss_ = msd(df, "i_h2o2_nA", ["Sterile", "Sterile+CAP-300"])
    ax.errorbar([0, 300], ms_ / 10, ss_ / 10, fmt="s--", color="#1baf7a", ms=4, capsize=2, label="sterile (/10)")
    ax.set_xlabel("CAP exposure (s)"); ax.set_ylabel("PB-SPCE current (nA)"); ax.legend()
    panel(ax, "c")
    fig.tight_layout(w_pad=1.5)
    save(fig, "fig2_chemistry")


def weibull_t(t, b, p):
    return (t / b) ** p


def fig_inactivation(df):
    fig, axs = plt.subplots(1, 2, figsize=(W2 * 0.72, 2.2))
    t = np.array(TIMES, float)
    sub = df[df.group.isin(INFECTED[:6])]
    ax = axs[0]
    for g, tt in zip(INFECTED[:6], TIMES):
        y = sub.loc[sub.group == g, "log_reduction"]
        ax.scatter(np.full(len(y), tt) + np.linspace(-4, 4, len(y)), y, s=9, color=COLORS[g], zorder=3, lw=0)
    popt, pcov = curve_fit(weibull_t, sub.t_treat[sub.t_treat > 0], sub.log_reduction[sub.t_treat > 0], p0=(100, 0.8))
    tt = np.linspace(0, 300, 200)
    ax.plot(tt, weibull_t(tt, *popt), color="#0b0b0b", lw=1)
    g = df[df.group == "Gas-300"].log_reduction
    ax.scatter(np.full(len(g), 310), g, s=9, color=COLORS["Gas-300"], marker="D", lw=0, label="Ar gas only")
    ax.set_xlabel("CAP exposure (s)"); ax.set_ylabel("log$_{10}$ reduction (CFU)")
    ax.text(10, weibull_t(300, *popt) * 0.92, f"Weibull fit: b = {popt[0]:.0f} s, p = {popt[1]:.2f}", fontsize=6.3)
    ax.legend(loc="lower right")
    panel(ax, "a")
    ax = axs[1]
    tr = TRACES[("CAP-120", 0)]["bio"]
    tr3 = TRACES[("CAP-300", 0)]["bio"]
    zb = TRACES[("CAP-120", 0)]["deliv"]["z"][TRACES[("CAP-120", 0)]["deliv"]["in_bio"]] * 1e6
    zb3 = TRACES[("CAP-300", 0)]["deliv"]["z"][TRACES[("CAP-300", 0)]["deliv"]["in_bio"]] * 1e6
    ax.semilogx(tr["viable_profile"], zb, color=SEQ["CAP-120"], label="120 s")
    ax.semilogx(tr3["viable_profile"], zb3, color=SEQ["CAP-300"], label="300 s")
    ax.invert_yaxis(); ax.set_xlabel("Viable fraction N/N$_0$"); ax.set_ylabel("Depth in biofilm (µm)")
    ax.legend()
    panel(ax, "b")
    fig.tight_layout(w_pad=1.5)
    save(fig, "fig3_inactivation")
    return popt, np.sqrt(np.diag(pcov))


def fig_spectra(df, freqs):
    w = 2 * np.pi * freqs
    groups = ["Control", "CAP-60", "CAP-120", "CAP-300", "Gas-300", "Sterile"]
    fig, axs = plt.subplots(2, 2, figsize=(W2, 4.4))
    for g in groups:
        zs = np.array([TRACES[(g, r)]["z"] for r in range(N_REP)])
        zm = zs.mean(0)
        ls = "--" if g in ("Gas-300",) else ("-." if g == "Sterile" else "-")
        axs[0, 0].loglog(freqs, np.abs(zm), color=COLORS[g], ls=ls, label=g)
        axs[0, 1].semilogx(freqs, -np.degrees(np.angle(zm)), color=COLORS[g], ls=ls)
        e, s = apparent_permittivity(w, zm, 600.0)
        # de-embed electrode polarisation using each replicate's fitted CPE
        zc = np.array([TRACES[(g, r)]["z"] - 1 / (df[(df.group == g) & (df.replicate == r + 1)].Q_ep.iloc[0]
                       * (1j * w) ** df[(df.group == g) & (df.replicate == r + 1)].n_ep.iloc[0])
                       for r in range(N_REP)]).mean(0)
        e, s = apparent_permittivity(w, zc, 600.0)
        sel = freqs > 1e5
        axs[1, 0].semilogx(freqs[sel], e[sel], color=COLORS[g], ls=ls)
        axs[1, 1].semilogx(freqs[sel], s[sel], color=COLORS[g], ls=ls)
    # inset Nyquist of control fit
    zc0 = TRACES[("Control", 0)]
    axs[0, 0].set_xlabel("Frequency (Hz)"); axs[0, 0].set_ylabel("|Z| (Ω)")
    axs[0, 0].legend(ncol=2, loc="upper right")
    axs[0, 1].set_xlabel("Frequency (Hz)"); axs[0, 1].set_ylabel("-Phase (deg)")
    ax_in = axs[0, 1].inset_axes([0.50, 0.42, 0.30, 0.5])
    hf = freqs > 2e5
    ax_in.plot(zc0["z"].real[hf], -zc0["z"].imag[hf], "o", ms=2, color="#8a8984", mew=0)
    ax_in.plot(zc0["z_fit"].real[hf], -zc0["z_fit"].imag[hf], color="#0b0b0b", lw=0.8)
    ax_in.set_title("Nyquist > 200 kHz (Control)", fontsize=5.5)
    ax_in.tick_params(labelsize=5)
    axs[1, 0].set_xlabel("Frequency (Hz)"); axs[1, 0].set_ylabel("Apparent ε′ (EP-corrected)")
    axs[1, 1].set_xlabel("Frequency (Hz)"); axs[1, 1].set_ylabel("Apparent conductivity (S/m)")
    for a, l in zip(axs.ravel(), "abcd"):
        panel(a, l)
    fig.tight_layout()
    save(fig, "fig4_spectra")


def bar_by_group(ax, df, col, groups, scale=1.0, ylabel=""):
    m, s = msd(df, col, groups)
    x = np.arange(len(groups))
    ax.bar(x, m * scale, yerr=s * scale, color=[COLORS[g] for g in groups], width=0.7,
           capsize=1.5, error_kw={"lw": 0.6}, ec="white", lw=0.6)
    for i, g in enumerate(groups):
        y = df.loc[df.group == g, col] * scale
        ax.scatter(np.full(len(y), i) + np.linspace(-0.18, 0.18, len(y)), y, s=3, c="#0b0b0b", lw=0, zorder=3)
    ax.set_xticks(x); ax.set_xticklabels([g.replace("Sterile+CAP-300", "Ster.+CAP").replace("Sterile", "Ster.") for g in groups],
                                         rotation=45, ha="right")
    ax.set_ylabel(ylabel); ax.grid(axis="x", visible=False)


def fig_parameters(df):
    fig, axs = plt.subplots(2, 3, figsize=(W2, 4.2))
    inf = INFECTED
    bar_by_group(axs[0, 0], df, "d_eps", inf, ylabel="Dispersion strength Δε")
    bar_by_group(axs[0, 1], df, "f_c", inf, 1e-6, ylabel="Characteristic freq. f$_c$ (MHz)")
    bar_by_group(axs[0, 2], df, "sigma", ORDER, ylabel="σ (S/m)")
    axs[0, 1].set_ylim(10, None)
    axs[0, 2].set_ylim(0.5, None)
    bar_by_group(axs[1, 0], df, "alpha", inf, ylabel="Cole–Cole α")
    bar_by_group(axs[1, 1], df, "absZ_100k", ORDER, ylabel="|Z| at 100 kHz (Ω)")
    axs[1, 1].set_ylim(500, None)
    ax = axs[1, 2]
    sub = df[df.group.isin(INFECTED)]
    for g in INFECTED:
        s = sub[sub.group == g]
        ax.scatter(s.perm_frac, s.d_eps, s=10, color=COLORS[g], lw=0, label=g)
    r = stats.linregress(sub.perm_frac, sub.d_eps)
    xx = np.linspace(0, 1, 10)
    ax.plot(xx, r.intercept + r.slope * xx, color="#0b0b0b", lw=0.8)
    ax.text(0.05, 0.08, f"R$^2$ = {r.rvalue**2:.3f}", transform=ax.transAxes, fontsize=6.5)
    ax.set_xlabel("Permeabilised fraction (model truth)"); ax.set_ylabel("Δε (fitted)")
    ax.legend(fontsize=5.2, ncol=2, loc="upper right")
    for a, l in zip(axs.ravel(), "abcdef"):
        panel(a, l)
    fig.tight_layout()
    save(fig, "fig5_parameters")


def fig_electrochem(df, eis_freqs):
    fig, axs = plt.subplots(2, 2, figsize=(W2, 4.4))
    ax = axs[0, 0]
    for g in ["Sterile", "Control", "CAP-60", "CAP-300"]:
        tr = TRACES[(g, 0)]
        ax.plot(tr["E"], tr["I"] * 1e6, color=COLORS[g], label=g)
    ax.set_xlabel("E vs Ag/AgCl (V)"); ax.set_ylabel("Current (µA)"); ax.legend(loc="upper left")
    ax = axs[0, 1]
    for g in ["Sterile", "Control", "CAP-60", "CAP-300"]:
        tr = TRACES[(g, 0)]
        ax.plot(tr["z_eis"].real, -tr["z_eis"].imag, "o", ms=2.2, color=COLORS[g], mew=0, label=g)
        ax.plot(tr["z_eis_fit"].real, -tr["z_eis_fit"].imag, color=COLORS[g], lw=0.8)
    ax.set_xlabel("Z′ (Ω)"); ax.set_ylabel("−Z″ (Ω)"); ax.set_aspect("equal", "datalim")
    bar_by_group(axs[1, 0], df, "Rct_fit", ORDER, 1e-3, ylabel="R$_{ct}$ (kΩ)")
    bar_by_group(axs[1, 1], df, "cv_dEp", ORDER, 1e3, ylabel="Peak separation ΔE$_p$ (mV)")
    axs[1, 1].set_ylim(100, None)
    for a, l in zip(axs.ravel(), "abcd"):
        panel(a, l)
    fig.tight_layout()
    save(fig, "fig6_electrochem")


def fig_ml(df, freqs, pls, fstat, vip):
    fig, axs = plt.subplots(1, 3, figsize=(W2, 2.3))
    ax = axs[0]
    ax.semilogx(freqs, fstat["F_mag"], color="#2a78d6", label="log|Z|")
    ax.semilogx(freqs, fstat["F_phase"], color="#eb6834", label="phase")
    ax.set_yscale("log"); ax.set_xlabel("Frequency (Hz)"); ax.set_ylabel("ANOVA F (infected groups)")
    ax.legend(); panel(ax, "a")
    ax = axs[1]
    y = PLS_Y
    for g in INFECTED:
        msk = (PLS_GROUPS == g)
        ax.scatter(y[msk], pls["pred"][msk], s=10, color=COLORS[g], lw=0, label=g)
    lim = [y.min() - 0.2, y.max() + 0.2]
    ax.plot(lim, lim, color="#8a8984", lw=0.7, ls="--")
    ax.set_xlabel("True log$_{10}$ CFU cm$^{-2}$"); ax.set_ylabel("Predicted (nested LOO-CV)")
    ax.text(0.05, 0.86, f"Q$^2$ = {pls['q2']:.3f}\nRMSECV = {pls['rmsecv']:.2f} log\n{pls['n_comp']} LV",
            transform=ax.transAxes, fontsize=6.3)
    ax.legend(fontsize=5, loc="lower right"); panel(ax, "b")
    ax = axs[2]
    nf = freqs.size
    ax.semilogx(freqs, vip[:nf], color="#2a78d6", label="log|Z|")
    ax.semilogx(freqs, vip[nf:], color="#eb6834", label="phase")
    ax.axhline(1, color="#8a8984", lw=0.6, ls="--")
    ax.set_xlabel("Frequency (Hz)"); ax.set_ylabel("PLS VIP score"); ax.legend(); panel(ax, "c")
    fig.tight_layout(w_pad=1.5)
    save(fig, "fig7_regression")


def fig_kk_sens(df, freqs, sens):
    fig, axs = plt.subplots(1, 2, figsize=(W2, 2.3))
    ax = axs[0]
    rr, ri = TRACES[("Control", 0)]["kk"]
    ax.semilogx(freqs, rr * 100, "o", ms=2, color="#2a78d6", mew=0, label="real")
    ax.semilogx(freqs, ri * 100, "s", ms=2, color="#eb6834", mew=0, label="imag")
    ax.axhline(0, color="#8a8984", lw=0.6)
    ax.set_xlabel("Frequency (Hz)"); ax.set_ylabel("Lin-KK residual (%)"); ax.legend()
    ax.text(0.03, 0.06, f"all 54 spectra: RMS = {df.kk_rms.mean()*100:.2f} ± {df.kk_rms.std()*100:.2f} %",
            transform=ax.transAxes, fontsize=6.3)
    panel(ax, "a")
    ax = axs[1]
    s = sens.sort_values("span")
    y = np.arange(len(s))
    ax.barh(y, s["minus"] - s["base"], left=s["base"], color="#86b6ef", height=0.6, label="-20 %")
    ax.barh(y, s["plus"] - s["base"], left=s["base"], color="#1c5cab", height=0.6, label="+20 %")
    ax.axvline(s["base"].iloc[0], color="#0b0b0b", lw=0.7)
    ax.set_yticks(y); ax.set_yticklabels(s["label"])
    ax.set_xlabel("Relative Δε loss, CAP-300 vs control (%)"); ax.legend(loc="lower right")
    ax.grid(axis="y", visible=False)
    panel(ax, "b")
    fig.tight_layout(w_pad=2)
    save(fig, "fig8_validation")


# ----------------------------------------------------------------- analyses
def frequency_anova(df, freqs):
    sub = [(g, r) for g in INFECTED for r in range(N_REP)]
    mag = np.array([np.log10(np.abs(TRACES[k]["z"])) for k in sub])
    ph = np.array([np.angle(TRACES[k]["z"]) for k in sub])
    lab = np.array([k[0] for k in sub])
    Fm = [stats.f_oneway(*[mag[lab == g, i] for g in INFECTED])[0] for i in range(freqs.size)]
    Fp = [stats.f_oneway(*[ph[lab == g, i] for g in INFECTED])[0] for i in range(freqs.size)]
    return {"F_mag": np.array(Fm), "F_phase": np.array(Fp)}, mag, ph, lab


def sensitivity(freqs):
    """One-at-a-time +/-20 % sensitivity of the CAP-300 vs control Delta-eps contrast."""
    w = 2 * np.pi * freqs

    def contrast(cell=Bacterium(), phi=0.30, d_bio=100e-6, lam=45e-6, sig=1.0, perm=None):
        out = []
        for pf in (0.0, perm if perm is not None else PERM300):
            eb = biofilm_permittivity(w, 0.85 * sig, phi, pf, cell)
            eg = debye_medium(w, sig)
            z = IDESensor(lam=lam).impedance(w, eb, eg, d_bio if pf == 0 else d_bio - 12e-6)
            out.append(fit_dielectric(freqs, z, 600.0)["d_eps"])
        return 100 * (1 - out[1] / out[0])

    base = contrast()
    b = Bacterium()
    cases = {
        "cell radius": lambda k: contrast(cell=replace(b, radius=b.radius * k)),
        "membrane ε (C$_m$)": lambda k: contrast(cell=replace(b, eps_mem=b.eps_mem * k)),
        "cytoplasm σ": lambda k: contrast(cell=replace(b, sigma_cyt=b.sigma_cyt * k)),
        "cell volume fraction": lambda k: contrast(phi=0.30 * k),
        "biofilm thickness": lambda k: contrast(d_bio=100e-6 * k),
        "IDE penetration length": lambda k: contrast(lam=45e-6 * k),
        "medium conductivity": lambda k: contrast(sig=1.0 * k),
        "permeabilised fraction": lambda k: contrast(perm=min(1.0, PERM300 * k)),
    }
    rows = [{"label": k, "base": base, "minus": f(0.8), "plus": f(1.2)} for k, f in cases.items()]
    s = pd.DataFrame(rows)
    s["span"] = (s.plus - s.minus).abs()
    return s


def write_tables(df, an, pls, wb, lod, sp):
    cols = [("log_reduction", "log$_{10}$ reduction", 1, 2), ("perm_frac", "Permeabilised fraction", 1, 2),
            ("H2O2_bio_mM", "H$_2$O$_2$ in biofilm ($\\mu$M)", 1e3, 1),
            ("pH_meas", "Surface pH (IrO$_x$)", 1, 2),
            ("sigma", "$\\sigma$ (S\\,m$^{-1}$)", 1, 3), ("d_eps", "$\\Delta\\varepsilon$", 1, 0),
            ("f_c", "$f_c$ (MHz)", 1e-6, 1), ("alpha", "$\\alpha$", 1, 3),
            ("Rct_fit", "$R_{ct}$ (k$\\Omega$)", 1e-3, 2), ("cv_dEp", "$\\Delta E_p$ (mV)", 1e3, 0),
            ("i_h2o2_nA", "$i_{\\mathrm{H_2O_2}}$ (nA)", 1, 0)]
    g = df.groupby("group")
    lines = []
    for c, lab, sc, nd in cols:
        cells = []
        for grp in ORDER:
            v = g.get_group(grp)[c] * sc
            if v.isna().all() or (c in ("d_eps", "f_c", "alpha") and not g.get_group(grp).infected.iloc[0]):
                cells.append("n.d.")
            else:
                cells.append(f"{v.mean():.{nd}f}$\\pm${v.std():.{nd}f}")
        lines.append(lab + " & " + " & ".join(cells) + " \\\\")
    (RES / "table_groups.tex").write_text("\n".join(lines) + "\n")
    alines = []
    for k, a in an.items():
        alines.append(f"{k} & {a['F']:.1f} & {a['p']:.1e} & {a['eta2']:.3f} \\\\")
    (RES / "table_anova.tex").write_text("\n".join(alines) + "\n")


def fig_graphical_abstract(df):
    """Single-panel graphical abstract: phantom + jet, channels, two headline results."""
    fig = plt.figure(figsize=(W2, 2.6))
    gs = fig.add_gridspec(1, 3, width_ratios=[1.25, 1, 1], wspace=0.45)
    ax = fig.add_subplot(gs[0]); ax.axis("off"); ax.set_xlim(0, 5); ax.set_ylim(0, 4)
    ax.add_patch(plt.Rectangle((0.2, 0.6), 4.6, 1.1, fc="#cde2fb", ec="#256abf", lw=0.8))
    ax.add_patch(plt.Rectangle((0.2, 1.7), 4.6, 0.35, fc="#f6d3a8", ec="#c98500", lw=0.8))
    rng = np.random.default_rng(1)
    ax.scatter(rng.uniform(0.3, 4.7, 140), rng.uniform(1.73, 2.02, 140), s=3, c="#8a5a00", lw=0)
    ax.add_patch(plt.Polygon([[2.2, 3.8], [2.8, 3.8], [2.65, 3.0], [2.35, 3.0]], fc="#d9d8d3", ec="#52514e", lw=0.8))
    ax.add_patch(plt.Polygon([[2.38, 3.0], [2.62, 3.0], [2.75, 2.1], [2.25, 2.1]], fc="#b9a6f0", ec="none"))
    ax.text(2.95, 3.4, "Ar CAP jet\nH$_2$O$_2$, NO$_2^-$, NO$_3^-$", fontsize=6.3, va="center")
    ax.text(2.5, 1.15, "Buffered agarose wound bed", ha="center", fontsize=6.3)
    ax.text(0.25, 2.35, "S. aureus biofilm", fontsize=6.3)
    ax.text(2.5, 0.2, "IDE 10 Hz–100 MHz  |  SPCE: EIS, CV, H$_2$O$_2$", ha="center", fontsize=6.3,
            fontweight="bold")
    ax.set_title("In silico infected wound phantom", fontsize=7.5)
    t = np.array(TIMES)
    ax2 = fig.add_subplot(gs[1])
    m, s = msd(df, "log_reduction", INFECTED[:6])
    ax2.errorbar(t, m, s, fmt="o-", color="#0b0b0b", ms=3.5, capsize=2)
    ax2.set_xlabel("CAP exposure (s)"); ax2.set_ylabel("log$_{10}$ reduction")
    ax2.set_title("Biofilm inactivation", fontsize=7.5)
    ax3 = fig.add_subplot(gs[2])
    m, s = msd(df, "f_c", INFECTED[:6])
    ax3.errorbar(t, m * 1e-6, s * 1e-6, fmt="o-", color="#2a78d6", ms=3.5, capsize=2)
    ax3.set_xlabel("CAP exposure (s)"); ax3.set_ylabel("β-dispersion f$_c$ (MHz)")
    ax3.set_title("Label-free impedance marker", fontsize=7.5)
    ax3.text(0.04, 0.95, "η² = 0.97\nρ = −0.95 vs CFU", transform=ax3.transAxes, va="top", fontsize=6)
    save(fig, "graphical_abstract")


def main():
    global TRACES, N_REP, PLS_Y, PLS_GROUPS, PERM300
    ap = argparse.ArgumentParser()
    ap.add_argument("--n-rep", type=int, default=6)
    ap.add_argument("--seed", type=int, default=2026)
    a = ap.parse_args()
    N_REP = a.n_rep
    df, TRACES, freqs, eis_freqs = pl.run_experiment(n_rep=N_REP, seed=a.seed)
    df.to_csv(RES / "samples.csv", index=False)
    summary = df.groupby("group", sort=False).agg(["mean", "std"])
    summary.to_csv(RES / "group_summary.csv")
    PERM300 = df[df.group == "CAP-300"].perm_frac.mean()

    inf = df[df.group.isin(INFECTED)]
    an = {}
    for col, lab in [("log_reduction", "log$_{10}$ reduction"), ("d_eps", "$\\Delta\\varepsilon$"),
                     ("f_c", "$f_c$"), ("sigma", "$\\sigma$"), ("absZ_100k", "$|Z|_{100\\,\\mathrm{kHz}}$"),
                     ("Rct_fit", "$R_{ct}$"), ("cv_dEp", "$\\Delta E_p$"), ("i_h2o2_nA", "$i_{\\mathrm{H_2O_2}}$")]:
        an[lab] = anova_tukey(inf, col, "group", INFECTED)
    pairs = {k: v["pairs"] for k, v in an.items()}
    pd.concat(pairs, names=["variable"]).to_csv(RES / "tukey_pairs.csv")
    ster = df[df.group.isin(["Sterile", "Sterile+CAP-300"])]
    t_ster = {c: stats.ttest_ind(ster[ster.group == "Sterile"][c], ster[ster.group == "Sterile+CAP-300"][c],
                                 equal_var=False) for c in ("sigma", "absZ_100k", "Rct_fit")}
    feats = ["d_eps", "f_c", "alpha", "sigma", "absZ_1k", "absZ_100k", "absZ_10M", "phase_5M",
             "Rct_fit", "cv_dEp", "cv_ipa", "i_h2o2_nA", "pH_meas"]
    sp = spearman_table(inf, feats, "log10_cfu")
    sp.to_csv(RES / "spearman.csv", index=False)

    fstat, mag, ph, lab = frequency_anova(df, freqs)
    X = np.hstack([mag, ph])
    PLS_Y = np.array([inf[(inf.group == g) & (inf.replicate == r + 1)].log10_cfu.iloc[0]
                      for g in INFECTED for r in range(N_REP)])
    PLS_GROUPS = lab
    pls = nested_pls(X, PLS_Y)
    vip = vip_scores(pls["model"][-1])
    # univariate calibration: Delta-eps vs log reduction
    lod = linear_lod(inf.log_reduction.to_numpy(), inf.d_eps.to_numpy())
    lod_rct = linear_lod(inf.log_reduction.to_numpy(), inf.Rct_fit.to_numpy())
    sens = sensitivity(freqs)
    sens.to_csv(RES / "sensitivity.csv", index=False)

    fig_schematic()
    fig_chemistry(df)
    wb, wb_se = fig_inactivation(df)
    fig_spectra(df, freqs)
    fig_parameters(df)
    fig_electrochem(df, eis_freqs)
    fig_ml(df, freqs, pls, fstat, vip)
    fig_kk_sens(df, freqs, sens)
    fig_graphical_abstract(df)
    write_tables(df, an, pls, wb, lod, sp)

    imax = int(np.argmax(fstat["F_mag"]))
    imaxp = int(np.argmax(fstat["F_phase"]))
    out = {
        "n_rep": N_REP, "seed": a.seed, "n_spectra": int(len(df)),
        "freq_points": int(freqs.size),
        "anova": {k: {kk: float(vv) for kk, vv in v.items() if kk != "pairs"} for k, v in an.items()},
        "sterile_welch": {k: {"t": float(v.statistic), "p": float(v.pvalue)} for k, v in t_ster.items()},
        "weibull_time": {"b_s": float(wb[0]), "p": float(wb[1]), "se_b": float(wb_se[0]), "se_p": float(wb_se[1])},
        "pls": {"q2": float(pls["q2"]), "rmsecv": float(pls["rmsecv"]), "n_comp": pls["n_comp"]},
        "lod_deps": {k: float(v) for k, v in lod.items()},
        "lod_rct": {k: float(v) for k, v in lod_rct.items()},
        "fmax_mag_Hz": float(freqs[imax]), "Fmax_mag": float(fstat["F_mag"][imax]),
        "fmax_phase_Hz": float(freqs[imaxp]), "Fmax_phase": float(fstat["F_phase"][imaxp]),
        "kk_rms_mean": float(df.kk_rms.mean()), "kk_rms_max": float(df.kk_rms.max()),
        "chi2_red_mean": float(df.chi2_red.mean()),
        "spearman": sp.set_index("feature").to_dict(orient="index"),
        "sensitivity": sens.to_dict(orient="records"),
        "perm300": float(PERM300),
    }
    (RES / "summary.json").write_text(json.dumps(out, indent=2, default=float))
    print(json.dumps({k: out[k] for k in ("anova", "sterile_welch", "weibull_time", "pls", "lod_deps",
                                          "lod_rct", "fmax_mag_Hz", "fmax_phase_Hz", "kk_rms_mean")},
                     indent=1, default=float))
    print(summary.loc[:, (["log_reduction", "perm_frac", "d_eps", "f_c", "sigma", "Rct_fit", "cv_dEp", "i_h2o2_nA",
                           "pH_meas", "H2O2_bio_mM", "absZ_100k"], "mean")].round(4).to_string())


if __name__ == "__main__":
    main()
