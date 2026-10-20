#!/usr/bin/env python3
"""Stage 4 — Aggregate per polymer and run the virtual microglia cell.

Converts the Stage-2 ML signal (and Stage-3 docking when available) into a single
receptor-activation drive `u` per polymer, runs the illustrative microglia ODE
model, and reports a normalised neuroinflammation index per polymer.

Aggregation (documented, deterministic):
  * Per polymer, per endpoint: take the mean of that polymer's TOP-3 leachate
    probabilities (dominant-actor view; concentrations are unknown so abundance
    is not weighted).
  * Combine endpoints into u with mechanistic weights:
        NR-AhR 0.50  (primary neuroinflammatory receptor axis)
        SR-ARE 0.30  (oxidative-stress / Nrf2)
        SR-p53 0.20  (genotoxic stress)
  * u is already in [0,1] (probabilities), so it feeds the ODE drive directly.

Stage-3 docking: if results/docking.csv has rows, a normalised docking term would
be blended in here; when it is empty (structures unreachable) the stage records
used_docking=False and relies on the ML signal, matching docking_status.md.

Inputs : results/ml_scores.csv, results/docking.csv
Output : results/neuroinflammation_index.csv
         (polymer, index, receptor_activation_u, dominant_receptor,
          dominant_leachate, top_endpoint_scores, used_docking)
"""
from __future__ import annotations

import os
import sys

import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(ROOT, "models"))
RESULTS = os.path.join(ROOT, "results")

from microglia import max_cytokine_sum, simulate  # noqa: E402

ENDPOINT_WEIGHTS = {"NR-AhR": 0.50, "SR-ARE": 0.30, "SR-p53": 0.20}
TOP_K = 3
POLYMERS = ["PE", "PP", "PS", "PVC", "PU"]


def top_k_mean(series: pd.Series, k: int = TOP_K) -> float:
    vals = sorted(series.tolist(), reverse=True)[:k]
    return float(sum(vals) / len(vals)) if vals else 0.0


def main() -> int:
    scores = pd.read_csv(os.path.join(RESULTS, "ml_scores.csv"))
    dock_path = os.path.join(RESULTS, "docking.csv")
    used_docking = False
    if os.path.exists(dock_path):
        dock = pd.read_csv(dock_path)
        used_docking = len(dock) > 0
    print(f"[S4] docking signal available: {used_docking} "
          f"(ML-only fallback if False)")

    mx = max_cytokine_sum()
    rows = []
    for pol in POLYMERS:
        sub = scores[scores["polymer"] == pol]
        if sub.empty:
            print(f"[S4] WARNING: no leachate scores for {pol}; skipping")
            continue
        ep_scores = {}
        for ep, w in ENDPOINT_WEIGHTS.items():
            ep_sub = sub[sub["endpoint"] == ep]
            ep_scores[ep] = top_k_mean(ep_sub["prob"])
        u = sum(ENDPOINT_WEIGHTS[ep] * ep_scores[ep] for ep in ENDPOINT_WEIGHTS)

        # dominant receptor = highest weighted endpoint contribution
        contrib = {ep: ENDPOINT_WEIGHTS[ep] * ep_scores[ep] for ep in ENDPOINT_WEIGHTS}
        dominant_receptor = max(contrib, key=contrib.get)
        dom_sub = sub[sub["endpoint"] == dominant_receptor]
        dominant_leachate = (
            dom_sub.sort_values("prob", ascending=False)["chemical_name"].iloc[0]
            if not dom_sub.empty else "n/a"
        )

        out = simulate(u)
        index = out.cytokine_sum / mx
        rows.append({
            "polymer": pol,
            "index": round(index, 4),
            "receptor_activation_u": round(u, 4),
            "dominant_receptor": dominant_receptor,
            "dominant_leachate": dominant_leachate,
            "NR-AhR_top3": round(ep_scores["NR-AhR"], 4),
            "SR-ARE_top3": round(ep_scores["SR-ARE"], 4),
            "SR-p53_top3": round(ep_scores["SR-p53"], 4),
            "used_docking": used_docking,
        })
        print(f"[S4] {pol}: u={u:.3f} index={index:.3f} "
              f"dominant={dominant_receptor} ({dominant_leachate})")

    df = pd.DataFrame(rows).sort_values("index", ascending=False).reset_index(drop=True)
    out_path = os.path.join(RESULTS, "neuroinflammation_index.csv")
    df.to_csv(out_path, index=False)
    print(f"[S4] wrote {out_path}")
    print("[S4] predicted ranking (high -> low): "
          + " > ".join(df["polymer"].tolist()))
    return 0


if __name__ == "__main__":
    sys.exit(main())
