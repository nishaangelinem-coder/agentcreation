"""Draw the manuscript figures at print resolution (PDF vector + 600 dpi TIFF/PNG).

Fig. 1  Four-layer architecture of a closed-loop CAP system (schematic).
Fig. 2  In-silico dose error by plant-gain bin (reads ../supplementary/S1_simulation_results.json).
Fig. 3  Staged roadmap with metric gates.

Usage:  python make_figures.py
"""

import json
import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
RESULTS = os.path.join(HERE, "..", "supplementary", "S1_simulation_results.json")

INK = "#1a1a1a"
INK2 = "#555555"
GRID = "#e3e3e3"
GREY = "#8c8c8c"  # open-loop reference arm
BLUE = "#1f5fa8"  # closed loop, nominal sensor
ORANGE = "#d9822b"  # closed loop, degraded sensor
HILITE = "#fbe9d6"

plt.rcParams.update({
    "font.family": "DejaVu Sans",
    "font.size": 8,
    "axes.edgecolor": INK2,
    "axes.labelcolor": INK,
    "xtick.color": INK2,
    "ytick.color": INK2,
    "axes.linewidth": 0.6,
    "pdf.fonttype": 42,
})

SINGLE_COL = 3.54  # in (90 mm)
DOUBLE_COL = 7.48  # in (190 mm)


def save(fig, name):
    for ext, kw in (("pdf", {}), ("tiff", {"dpi": 600, "pil_kwargs": {"compression": "tiff_lzw"}}),
                    ("png", {"dpi": 600})):
        fig.savefig(os.path.join(HERE, f"{name}.{ext}"), bbox_inches="tight", **kw)
    plt.close(fig)


def box(ax, x, y, w, h, title, body, face="#ffffff", edge=INK2, lw=0.8):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.02,rounding_size=0.06",
                                facecolor=face, edgecolor=edge, linewidth=lw))
    ax.text(x + w / 2, y + h - 0.12, title, ha="center", va="top", fontsize=6.5, weight="bold", color=INK)
    ax.text(x + w / 2, y + h - 0.42, body, ha="center", va="top", fontsize=6.6, color=INK2, linespacing=1.35)


def arrow(ax, p0, p1, color=INK2, style="-|>", lw=0.9, ls="-", rad=0.0):
    ax.add_patch(FancyArrowPatch(p0, p1, arrowstyle=style, mutation_scale=9, color=color, lw=lw,
                                 linestyle=ls, connectionstyle=f"arc3,rad={rad}"))


def fig_architecture():
    fig, ax = plt.subplots(figsize=(DOUBLE_COL, 3.0))
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 4.25)
    ax.axis("off")
    w, h, y = 2.0, 1.55, 2.1
    xs = [0.2, 2.65, 5.1, 7.55]
    layers = [
        ("1  Embedded sensing", "Electrochemical RONS\n(H$_2$O$_2$, NO, ORP)\nOES · V/I probes · IR"),
        ("2  State estimation", "Signal fusion, drift\ncorrection, dose\nestimate $\\hat{D}(t)$"),
        ("3  Feedback control", "Constrained MPC,\nadaptive / learning\n+ safety supervisor"),
        ("4  Actuation", "Power, duty cycle,\ngas flow, stand-off;\nDBD / jet / patch"),
    ]
    maturity = ["Least mature\n(drift, fouling)", "Emerging", "Most mature\n(lab-demonstrated)", "Hardware exists;\nrarely co-integrated"]
    faces = [HILITE, "#ffffff", "#e6eef8", "#ffffff"]
    for x, (t, b), m, f in zip(xs, layers, maturity, faces):
        box(ax, x, y, w, h, t, b, face=f)
        ax.text(x + w / 2, y - 0.12, m, ha="center", va="top", fontsize=6.4, color=INK2, style="italic")
    for i in range(3):
        arrow(ax, (xs[i] + w + 0.03, y + h / 2), (xs[i + 1] - 0.03, y + h / 2))
    # Plant
    ax.add_patch(FancyBboxPatch((2.9, 0.15), 4.2, 0.75, boxstyle="round,pad=0.02,rounding_size=0.06",
                                facecolor="#f4f4f4", edgecolor=INK2, linewidth=0.8))
    ax.text(5.0, 0.53, "Plasma–target system: tissue, liquid, gas phase, heat",
            ha="center", va="center", fontsize=6.6, color=INK)
    arrow(ax, (xs[3] + w / 2, y - 0.62), (7.12, 0.55), rad=-0.25)
    arrow(ax, (2.88, 0.55), (xs[0] + w / 2, y - 0.62), rad=-0.25)
    ax.text(8.95, 0.7, "plasma dose", fontsize=6.4, color=INK2, ha="center")
    ax.text(1.05, 0.7, "measured signals", fontsize=6.4, color=INK2, ha="center")
    ax.text(5.0, 3.93, "Target dose $D_{ref}$ and safety limits ($T_{max}$, UV, leakage current)",
            ha="center", va="bottom", fontsize=7, color=INK)
    arrow(ax, (xs[2] + w / 2, 3.74), (xs[2] + w / 2, y + h + 0.03))
    save(fig, "Fig1_architecture")


def fig_dose_error():
    with open(RESULTS) as fh:
        res = json.load(fh)
    arms = [("open_loop", "Open loop", GREY, "o"),
            ("closed_nominal", "Closed loop, nominal sensor", BLUE, "s"),
            ("closed_degraded", "Closed loop, degraded sensor", ORANGE, "^")]
    bins = [b["bin"] for b in res["arms"]["open_loop"]["by_gain_bin"]]
    cap = 500.0
    fig, ax = plt.subplots(figsize=(DOUBLE_COL * 0.62, 3.0))
    offs = [-0.22, 0.0, 0.22]
    for (key, label, col, mk), off in zip(arms, offs):
        data = res["arms"][key]["by_gain_bin"]
        for i, b in enumerate(data):
            x = i + off
            lo, hi, med = b["p10_pct"], b["p90_pct"], b["median_pct"]
            ax.plot([x, x], [lo, min(hi, cap)], color=col, lw=1.4, solid_capstyle="round", zorder=2)
            ax.plot(x, min(med, cap), marker=mk, ms=5, color=col, mec="white", mew=0.6, zorder=3,
                    label=label if i == 0 else None, ls="none")
            if hi > cap:
                ax.annotate(f"p90 +{hi:.0f}%\n(clipped)", xy=(x, cap), xytext=(x - 0.55, cap - 60),
                            fontsize=6.2, color=INK2, ha="center", va="top",
                            arrowprops=dict(arrowstyle="-", color=INK2, lw=0.5))
                ax.plot(x, cap, marker=(3, 0, 0), ms=4, color=col, zorder=3)
    ax.axhline(0, color=INK, lw=0.8, zorder=1)
    ax.text(-0.45, 8, "target", fontsize=6.4, color=INK2, ha="left", va="bottom")
    ax.set_xticks(range(len(bins)))
    ax.set_xticklabels([b.replace("-", "–") for b in bins])
    ax.set_xlabel("Plant-gain multiplier K / K$_{nom}$")
    ax.set_ylabel("CEM43 dose error (%)")
    ax.set_ylim(-110, cap + 10)
    ax.set_xlim(-0.5, len(bins) - 0.5)
    ax.yaxis.grid(True, color=GRID, lw=0.5)
    ax.set_axisbelow(True)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    ax.legend(frameon=False, loc="lower left", bbox_to_anchor=(0, 1.0), ncol=3, fontsize=6.6,
              handletextpad=0.2, columnspacing=1.0)
    save(fig, "Fig2_dose_error")


def fig_roadmap():
    fig, ax = plt.subplots(figsize=(DOUBLE_COL, 2.2))
    ax.set_xlim(0, 10.6)
    ax.set_ylim(0, 2.6)
    ax.axis("off")
    stages = [
        ("Stage 1\nSensor-supervised", "Open-loop dosing +\nIR / current interlocks", "Feasible now", "#ffffff"),
        ("Stage 2\nThermal closed loop", "MPC on CEM43 +\noptical intensity", "Lab-demonstrated", "#e6eef8"),
        ("Stage 3\nSpecies feedback", "H$_2$O$_2$ / ORP sensors\ndrive the source", "Sensors only so far", HILITE),
        ("Stage 4\nMultimodal dose", "Fused chemical,\noptical, electrical", "Not yet shown", "#ffffff"),
    ]
    gates = ["Gate\nM3, M6", "Gate\nM1, M2, M4", "Gate\nM5"]
    w, h, y = 1.95, 1.75, 0.55
    xs = [0.05, 2.9, 5.75, 8.6]
    for x, (t, b, s, f) in zip(xs, stages):
        edge = ORANGE if f == HILITE else INK2
        ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.02,rounding_size=0.06",
                                    facecolor=f, edgecolor=edge, linewidth=1.1 if f == HILITE else 0.8))
        ax.text(x + w / 2, y + h - 0.1, t, ha="center", va="top", fontsize=6.4, weight="bold", color=INK)
        ax.text(x + w / 2, y + h - 0.78, b, ha="center", va="top", fontsize=6.6, color=INK2, linespacing=1.3)
        ax.text(x + w / 2, y - 0.1, s, ha="center", va="top", fontsize=6.4, color=INK2, style="italic")
    for i, g in enumerate(gates):
        x0, x1 = xs[i] + w + 0.04, xs[i + 1] - 0.04
        arrow(ax, (x0, y + h / 2), (x1, y + h / 2))
        ax.text((x0 + x1) / 2, y + h / 2 - 0.12, g, ha="center", va="top", fontsize=5.6, color=INK,
                bbox=dict(boxstyle="round,pad=0.2", fc="white", ec=INK2, lw=0.5))
    save(fig, "Fig3_roadmap")


if __name__ == "__main__":
    fig_architecture()
    fig_dose_error()
    fig_roadmap()
    print("figures written to", HERE)
