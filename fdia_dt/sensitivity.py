"""Sensitivity of DT-based detectors to digital-twin forecast fidelity."""
import pickle
import sys
import numpy as np
import torch
torch.set_num_threads(4)
import pandas as pd

from .grid import Grid
from .data import node_features, redo_twin
from .models import GCNDetector, MLPDetector, train_torch, predict_torch
from .run_experiments import split, det_metrics, loc_metrics, standardise

if __name__ == "__main__":
    case = sys.argv[1]
    levels = [0.01, 0.02, 0.04, 0.06, 0.08] if case == "ieee14" else [0.01, 0.02, 0.04, 0.08]
    seeds = [0, 1, 2] if case == "ieee14" else [0, 1]
    d0 = pickle.load(open(f"data/{case}.pkl", "rb"))
    g = Grid(case)
    T = len(d0["y"])
    tr, va, te = split(T)
    y, N, kind = d0["y"].astype(float), d0["ynode"].astype(float), d0["kind"]
    st = te[kind[te] == 1]
    rows = []
    for lv in levels:
        d = redo_twin(g, d0, lv, seed=100)
        X, _ = node_features(g, d)
        X = standardise(X, tr[kind[tr] == 0]).astype(np.float32)
        ns = np.abs(X).max(-1)
        gs = ns.max(-1)
        tau = np.quantile(gs[va][kind[va] == 0], 0.99)
        rows.append(dict(level=lv, method="DT-Threshold", seed=0,
                         **det_metrics(y[te], (gs[te] > tau).astype(int), gs[te]),
                         **loc_metrics(N[st], (ns[st] > tau).astype(float))))
        for name, mk in [("DT-MLP", lambda: MLPDetector(g.n, X.shape[-1])),
                         ("DT-GCN (proposed)", lambda: GCNDetector(g.adj, X.shape[-1]))]:
            for sd in seeds:
                m, _ = train_torch(mk(), X[tr], y[tr], N[tr], X[va], y[va], N[va], seed=sd)
                pg, pn = predict_torch(m, X[te])
                rows.append(dict(level=lv, method=name, seed=sd,
                                 **det_metrics(y[te], (pg > 0.5).astype(int), pg),
                                 **loc_metrics(N[st], (pn[kind[te] == 1] > 0.5).astype(float))))
        print(lv, "done", flush=True)
    pd.DataFrame(rows).to_csv(f"results/sens_{case}.csv", index=False)
