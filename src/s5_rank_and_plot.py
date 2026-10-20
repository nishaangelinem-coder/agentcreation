#!/usr/bin/env python3
"""Stage 5 — Ranking, figures and summary.

Figures (PNG + SVG):
  fig1_polymer_index.*     per-polymer neuroinflammation index bar chart
  fig2_receptor_heatmap.*  endpoint(receptor-proxy) x polymer signal heatmap
  fig3_leachate_breakdown.* dominant-leachate contribution per polymer

Also writes results/summary.md, comparing the predicted ranking to
data/experimental_ranking.csv when the author has filled it in (skips the
comparison gracefully otherwise).
"""
from __future__ import annotations

import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
RESULTS = os.path.join(ROOT, "results")
FIG = os.path.join(RESULTS, "figures")
DATA = os.path.join(ROOT, "data")
POLYMERS = ["PE", "PP", "PS", "PVC", "PU"]
ENDPOINTS = ["NR-AhR", "SR-ARE", "SR-p53"]


def save(fig, stem):
    for ext in ("png", "svg"):
        fig.savefig(os.path.join(FIG, f"{stem}.{ext}"), bbox_inches="tight", dpi=150)
    plt.close(fig)


def fig_index(idx: pd.DataFrame):
    d = idx.sort_values("index", ascending=False)
    fig, ax = plt.subplots(figsize=(6, 4))
    bars = ax.bar(d["polymer"], d["index"], color="#3b6ea5", edgecolor="black")
    ax.set_ylabel("Predicted neuroinflammation index (0-1)")
    ax.set_xlabel("Polymer")
    ax.set_title("Predicted neuroinflammatory potential by polymer")
    ax.set_ylim(0, 1.05)
    for b, v in zip(bars, d["index"]):
        ax.text(b.get_x() + b.get_width() / 2, v + 0.02, f"{v:.2f}",
                ha="center", fontsize=9)
    ax.text(0.5, -0.22, "Computational prediction (illustrative ODE params); "
            "affinity != activation direction.", transform=ax.transAxes,
            ha="center", fontsize=7, style="italic", color="#555")
    save(fig, "fig1_polymer_index")


def fig_heatmap(scores: pd.DataFrame):
    mat = np.zeros((len(ENDPOINTS), len(POLYMERS)))
    for i, ep in enumerate(ENDPOINTS):
        for j, pol in enumerate(POLYMERS):
            sub = scores[(scores["endpoint"] == ep) & (scores["polymer"] == pol)]
            vals = sorted(sub["prob"].tolist(), reverse=True)[:3]
            mat[i, j] = sum(vals) / len(vals) if vals else 0.0
    fig, ax = plt.subplots(figsize=(6, 3.4))
    im = ax.imshow(mat, cmap="YlOrRd", vmin=0, vmax=max(0.5, mat.max()), aspect="auto")
    ax.set_xticks(range(len(POLYMERS)), POLYMERS)
    ax.set_yticks(range(len(ENDPOINTS)), ENDPOINTS)
    ax.set_title("Receptor-proxy signal (Tox21 top-3 mean) by polymer")
    for i in range(len(ENDPOINTS)):
        for j in range(len(POLYMERS)):
            ax.text(j, i, f"{mat[i, j]:.2f}", ha="center", va="center", fontsize=8)
    fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04, label="signal (0-1)")
    save(fig, "fig2_receptor_heatmap")


def fig_breakdown(idx: pd.DataFrame, scores: pd.DataFrame):
    fig, ax = plt.subplots(figsize=(7, 4))
    width = 0.25
    x = np.arange(len(POLYMERS))
    colors = {"NR-AhR": "#d1495b", "SR-ARE": "#edae49", "SR-p53": "#66a182"}
    for k, ep in enumerate(ENDPOINTS):
        vals = []
        for pol in POLYMERS:
            sub = scores[(scores["endpoint"] == ep) & (scores["polymer"] == pol)]
            top3 = sorted(sub["prob"].tolist(), reverse=True)[:3]
            vals.append(sum(top3) / len(top3) if top3 else 0.0)
        ax.bar(x + (k - 1) * width, vals, width, label=ep, color=colors[ep],
               edgecolor="black", linewidth=0.5)
    ax.set_xticks(x, POLYMERS)
    ax.set_ylabel("Top-3 mean Tox21 probability")
    ax.set_xlabel("Polymer")
    ax.set_title("Endpoint contribution breakdown by polymer")
    ax.legend(frameon=False, fontsize=8)
    save(fig, "fig3_leachate_breakdown")


def main() -> int:
    os.makedirs(FIG, exist_ok=True)
    idx = pd.read_csv(os.path.join(RESULTS, "neuroinflammation_index.csv"))
    scores = pd.read_csv(os.path.join(RESULTS, "ml_scores.csv"))
    metrics = pd.read_csv(os.path.join(RESULTS, "ml_metrics.csv"))

    fig_index(idx)
    fig_heatmap(scores)
    fig_breakdown(idx, scores)
    print(f"[S5] wrote 3 figures (png+svg) -> {FIG}")

    predicted = idx.sort_values("index", ascending=False)["polymer"].tolist()

    # Optional comparison to author-supplied experimental ranking.
    exp_path = os.path.join(DATA, "experimental_ranking.csv")
    comparison = ""
    if os.path.exists(exp_path):
        try:
            exp = pd.read_csv(exp_path)
            if {"polymer", "rank"}.issubset(exp.columns) and exp["rank"].notna().any():
                exp_order = exp.sort_values("rank")["polymer"].tolist()
                rho = _spearman(predicted, exp_order)
                comparison = (f"Experimental ranking (author-supplied): "
                              f"{' > '.join(exp_order)}\n\n"
                              f"Spearman rho (predicted vs experimental) = {rho:.2f}\n")
            else:
                comparison = "experimental_ranking.csv present but not yet filled in.\n"
        except Exception as exc:  # noqa: BLE001
            comparison = f"Could not parse experimental_ranking.csv: {exc}\n"
    else:
        comparison = ("No data/experimental_ranking.csv supplied yet — comparison "
                      "skipped. Author fills it from review screening to validate.\n")

    with open(os.path.join(RESULTS, "summary.md"), "w") as fh:
        fh.write("# Pipeline summary — Microplastics & Neuroinflammation\n\n")
        fh.write("## Predicted polymer ranking (high -> low neuroinflammatory potential)\n\n")
        fh.write("**" + " > ".join(predicted) + "**\n\n")
        fh.write(idx.to_markdown(index=False) + "\n\n")
        fh.write("## Stage-2 scaffold-split metrics (imbalance-aware)\n\n")
        fh.write(metrics.to_markdown(index=False) + "\n\n")
        fh.write("## Validation vs experiment\n\n")
        fh.write(comparison + "\n")
        fh.write("## Single most important caveat\n\n")
        fh.write("**Docking/ML affinity is binding STRENGTH, not activation "
                 "DIRECTION.** The index assumes the engaged receptors are driven "
                 "toward a pro-inflammatory state; agonist/antagonist direction must "
                 "be confirmed per leachate (here the Tox21 label supplies direction "
                 "for the AhR/stress axes). All virtual-cell kinetic parameters are "
                 "illustrative, not fitted. Stage-3 receptor docking (TLR4-MD2, "
                 "NLRP3) could not be run because structure databases were "
                 "unreachable in the execution sandbox (see docking_status.md); the "
                 "index therefore rests on the Tox21 ML signal. The ranking is a "
                 "falsifiable PREDICTION for wet-lab testing, not established fact.\n")
    print("[S5] wrote results/summary.md")
    print("[S5] predicted ranking: " + " > ".join(predicted))
    return 0


def _spearman(pred_order, exp_order):
    common = [p for p in pred_order if p in exp_order]
    pr = {p: i for i, p in enumerate(pred_order) if p in common}
    er = {p: i for i, p in enumerate(exp_order) if p in common}
    n = len(common)
    if n < 2:
        return float("nan")
    d2 = sum((pr[p] - er[p]) ** 2 for p in common)
    return 1 - 6 * d2 / (n * (n ** 2 - 1))


if __name__ == "__main__":
    sys.exit(main())
