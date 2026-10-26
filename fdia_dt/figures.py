"""Publication figures (IEEE single-column width)."""
import pickle
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import networkx as nx
from sklearn.metrics import roc_curve

from .grid import Grid

PAL = {"BDD (chi-square)": "#52514e", "DT-Threshold": "#eda100", "RandomForest (raw)": "#1baf7a",
       "GCN (raw)": "#6250d6", "DT-MLP": "#eb6834", "DT-GCN (proposed)": "#2a78d6"}
LS = {"BDD (chi-square)": ":", "DT-Threshold": "-.", "RandomForest (raw)": "--",
      "GCN (raw)": (0, (5, 1, 1, 1)), "DT-MLP": (0, (3, 1)), "DT-GCN (proposed)": "-"}
plt.rcParams.update({"font.family": "serif", "font.size": 8, "axes.linewidth": 0.6,
                     "axes.grid": True, "grid.color": "#dddddd", "grid.linewidth": 0.5,
                     "axes.spines.top": False, "axes.spines.right": False,
                     "legend.frameon": False, "legend.fontsize": 6.5, "savefig.dpi": 300})
W = 3.45
CASES = [("ieee14", "IEEE 14-bus"), ("ieee118", "IEEE 118-bus")]
S = {c: pickle.load(open(f"results/store_{c}.pkl", "rb")) for c, _ in CASES}


def roc():
    fig, ax = plt.subplots(1, 2, figsize=(W * 2, 2.3))
    for a, (c, title) in zip(ax, CASES):
        s = S[c]
        yt = s["y"][s["te"]]
        for m, col in PAL.items():
            fpr, tpr, _ = roc_curve(yt, s["scores"][m])
            a.plot(fpr, tpr, color=col, ls=LS[m], lw=1.4 if "proposed" in m else 1.0, label=m)
        a.set_xscale("log"); a.set_xlim(1e-3, 1); a.set_ylim(0, 1.01)
        a.set_xlabel("False-positive rate (log)"); a.set_ylabel("True-positive rate"); a.set_title(title)
    ax[1].legend(loc="lower right")
    fig.tight_layout(); fig.savefig("figures/roc.pdf"); fig.savefig("figures/roc.png")


def recall_vs_mag():
    fig, ax = plt.subplots(1, 2, figsize=(W * 2, 2.2), sharey=True)
    bins = np.array([0.001, 0.004, 0.007, 0.010, 0.014, 0.020])
    rows = []
    for a, (c, title) in zip(ax, CASES):
        s = S[c]; te = s["te"]; kt = s["kind"][te]; mg = s["mag"][te]
        st = kt == 1
        cen = np.degrees((bins[:-1] + bins[1:]) / 2)
        for m, col in PAL.items():
            sc = s["scores"][m]
            if m == "BDD (chi-square)":
                pred = sc > s["thr"]
            elif m == "DT-Threshold":
                pred = sc > s["tau"]
            else:
                pred = sc > 0.5
            r = [pred[st & (mg >= lo) & (mg < hi)].mean() for lo, hi in zip(bins[:-1], bins[1:])]
            rows += [dict(case=c, method=m, bin_deg=round(x, 3), recall=v) for x, v in zip(cen, r)]
            a.plot(cen, r, color=col, ls=LS[m], marker="o", ms=3, lw=1.4 if "proposed" in m else 1.0, label=m)
        a.set_xlabel(r"Attack magnitude $|c_\theta|$ (deg)"); a.set_title(title)
    ax[0].set_ylabel("Detection rate (stealthy FDIA)")
    ax[1].legend(loc="lower right")
    fig.tight_layout(); fig.savefig("figures/recall_mag.pdf"); fig.savefig("figures/recall_mag.png")
    pd.DataFrame(rows).to_csv("results/recall_vs_magnitude.csv", index=False)


def sensitivity():
    fig, ax = plt.subplots(1, 2, figsize=(W * 2, 2.1))
    for c, title in CASES:
        try:
            df = pd.read_csv(f"results/sens_{c}.csv")
        except FileNotFoundError:
            continue
        g = df.groupby(["method", "level"])[["f1", "loc_f1"]].agg(["mean", "std"]).reset_index()
        for m in ["DT-Threshold", "DT-MLP", "DT-GCN (proposed)"]:
            q = g[g.method == m]
            for a, k in zip(ax, ["f1", "loc_f1"]):
                a.errorbar(q.level * 100, q[(k, "mean")], q[(k, "std")], color=PAL[m],
                           ls="-" if c == "ieee14" else "--", marker="o" if c == "ieee14" else "s",
                           ms=3, lw=1.0, capsize=2, label=f"{m}, {title}")
    ax[0].set_ylabel("Detection F1"); ax[1].set_ylabel("Localisation F1")
    for a in ax:
        a.set_xlabel(r"Twin forecast error $\sigma_f$ (%)")
    ax[1].legend(loc="lower left", fontsize=5.5)
    fig.tight_layout(); fig.savefig("figures/sensitivity.pdf"); fig.savefig("figures/sensitivity.png")


def timeline():
    s = S["ieee14"]; te = s["te"]
    k = s["kind"][te]
    sl = slice(0, 96)  # first test day
    t = np.arange(96) * 0.25
    fig, ax = plt.subplots(2, 1, figsize=(W, 2.6), sharex=True)
    J = s["J"][te][sl]
    ax[0].plot(t, J, color=PAL["BDD (chi-square)"], lw=0.9)
    ax[0].axhline(s["thr"], color="#e34948", lw=0.8, ls="--")
    ax[0].set_ylabel(r"BDD $J(\hat x)$"); ax[0].set_yscale("log")
    p = s["scores"]["DT-GCN (proposed)"][sl]
    ax[1].plot(t, p, color=PAL["DT-GCN (proposed)"], lw=0.9)
    ax[1].axhline(0.5, color="#e34948", lw=0.8, ls="--")
    ax[1].set_ylabel("DT-GCN score"); ax[1].set_xlabel("Time of day (h)")
    for a in ax:
        for i in np.where(k[sl] == 1)[0]:
            a.axvspan(t[i] - 0.125, t[i] + 0.125, color="#e87ba4", alpha=0.25, lw=0)
        for i in np.where(k[sl] == 2)[0]:
            a.axvspan(t[i] - 0.125, t[i] + 0.125, color="#eda100", alpha=0.35, lw=0)
    ax[0].text(0.01, 0.9, "shaded: stealthy (pink) / naive (amber) attacks", transform=ax[0].transAxes, fontsize=6)
    fig.tight_layout(); fig.savefig("figures/timeline.pdf"); fig.savefig("figures/timeline.png")


def localisation_example():
    g = Grid("ieee14"); s = S["ieee14"]; te = s["te"]
    k = s["kind"][te]; N = s["N"][te]; P = s["locs"]["DT-GCN (proposed)"]
    idx = [i for i in np.where(k == 1)[0] if N[i].sum() == 3][0]
    G = nx.Graph(); G.add_edges_from(zip(g.f, g.t))
    pos = nx.kamada_kawai_layout(G)
    fig, ax = plt.subplots(1, 2, figsize=(W, 1.8))
    for a, val, title in [(ax[0], N[idx], "Ground truth"), (ax[1], P[idx], "DT-GCN probability")]:
        nx.draw_networkx_edges(G, pos, ax=a, width=0.6, edge_color="#999999")
        sc = nx.draw_networkx_nodes(G, pos, ax=a, node_size=90, node_color=val, cmap="Blues",
                                    vmin=0, vmax=1, edgecolors="#333333", linewidths=0.5)
        nx.draw_networkx_labels(G, pos, {i: str(i + 1) for i in G}, ax=a, font_size=5)
        a.set_title(title, fontsize=7); a.axis("off")
    fig.colorbar(sc, ax=ax, shrink=0.8, pad=0.02)
    fig.savefig("figures/localisation.pdf", bbox_inches="tight"); fig.savefig("figures/localisation.png", bbox_inches="tight")


if __name__ == "__main__":
    roc(); recall_vs_mag(); timeline(); localisation_example(); sensitivity()
    print("figures done")
