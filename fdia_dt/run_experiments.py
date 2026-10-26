"""Main experiments: detection/localisation benchmark, sensitivity studies, figures."""
import json
import pickle
import sys
import time
import numpy as np
import pandas as pd
from scipy.stats import chi2
from sklearn.ensemble import RandomForestClassifier
from sklearn.svm import SVC
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import roc_auc_score, roc_curve
import torch

from .grid import Grid
from .data import node_features, raw_node_features, redo_twin
from .models import GCNDetector, MLPDetector, train_torch, predict_torch

torch.set_num_threads(4)
SEEDS = [0, 1, 2, 3, 4]


def split(T):
    a, b = int(0.7 * T), int(0.85 * T)
    return np.arange(a), np.arange(a, b), np.arange(b, T)


def det_metrics(y, pred, score=None):
    tp = np.sum((pred == 1) & (y == 1)); fp = np.sum((pred == 1) & (y == 0))
    fn = np.sum((pred == 0) & (y == 1)); tn = np.sum((pred == 0) & (y == 0))
    p = tp / max(tp + fp, 1); r = tp / max(tp + fn, 1)
    out = dict(acc=(tp + tn) / len(y), prec=p, rec=r, f1=2 * p * r / max(p + r, 1e-12), fpr=fp / max(fp + tn, 1))
    out["auc"] = roc_auc_score(y, score) if score is not None else np.nan
    return out


def loc_metrics(N, P):
    tp = np.sum(P * N); fp = np.sum(P * (1 - N)); fn = np.sum((1 - P) * N)
    p = tp / max(tp + fp, 1); r = tp / max(tp + fn, 1)
    exact = np.mean(np.all(P == N, 1))
    return dict(loc_prec=p, loc_rec=r, loc_f1=2 * p * r / max(p + r, 1e-12), loc_exact=exact)


def standardise(X, idx_norm):
    mu = X[idx_norm].mean(0, keepdims=True)
    sd = X[idx_norm].std(0, keepdims=True) + 1e-6
    return (X - mu) / sd


def run_case(name, d, seeds=SEEDS, fast=False, store=None):
    g = Grid(name)
    T = len(d["y"])
    tr, va, te = split(T)
    y, N, kind, mag = d["y"].astype(float), d["ynode"].astype(float), d["kind"], d["mag"]
    Xdt, _ = node_features(g, d)
    Xraw = raw_node_features(g, d)
    norm_tr = tr[kind[tr] == 0]
    Xdt = standardise(Xdt, norm_tr).astype(np.float32)
    Xraw = standardise(Xraw, norm_tr).astype(np.float32)
    zs = StandardScaler().fit(d["z"][norm_tr]).transform(d["z"])
    rows, scores, locs = [], {}, {}
    stealth = te[kind[te] == 1]

    # 1) chi-square BDD
    thr = chi2.ppf(0.99, g.nz - (2 * g.n - 1))
    s = d["J"][te]
    rows.append(dict(method="BDD (chi-square)", seed=0, **det_metrics(y[te], (s > thr).astype(int), s)))
    scores["BDD (chi-square)"] = s

    # 2) DT residual threshold (no learning)
    node_s = np.abs(Xdt).max(-1)
    gs = node_s.max(-1)
    tau = np.quantile(gs[va][kind[va] == 0], 0.99)
    rows.append(dict(method="DT-Threshold", seed=0, **det_metrics(y[te], (gs[te] > tau).astype(int), gs[te]),
                     **loc_metrics(N[stealth], (node_s[stealth] > tau).astype(float))))
    scores["DT-Threshold"] = gs[te]

    timing = {}
    # 3) classical ML on raw telemetry
    for mname, mk in [("SVM-RBF (raw)", lambda s_: SVC(C=10, gamma="scale", probability=True, random_state=s_)),
                      ("RandomForest (raw)", lambda s_: RandomForestClassifier(300, n_jobs=4, random_state=s_))]:
        for sd in (seeds if "Forest" in mname else seeds[:1]):
            clf = mk(sd).fit(zs[tr], y[tr])
            t0 = time.time(); pr = clf.predict_proba(zs[te])[:, 1]; timing[mname] = (time.time() - t0) / len(te)
            rows.append(dict(method=mname, seed=sd, **det_metrics(y[te], (pr > 0.5).astype(int), pr)))
            scores.setdefault(mname, pr)

    # 4) deep models
    deep = [("MLP (raw)", lambda: MLPDetector(g.n, Xraw.shape[-1]), Xraw),
            ("GCN (raw)", lambda: GCNDetector(g.adj, Xraw.shape[-1]), Xraw),
            ("DT-MLP", lambda: MLPDetector(g.n, Xdt.shape[-1]), Xdt),
            ("DT-GCN (proposed)", lambda: GCNDetector(g.adj, Xdt.shape[-1]), Xdt)]
    hists = {}
    for mname, mk, X in deep:
        for sd in seeds:
            torch.manual_seed(sd)
            model, hist = train_torch(mk(), X[tr], y[tr], N[tr], X[va], y[va], N[va], seed=sd,
                                      epochs=40 if fast else 60)
            t0 = time.time(); pg, pn = predict_torch(model, X[te]); timing[mname] = (time.time() - t0) / len(te)
            Pst = (pn[kind[te] == 1] > 0.5).astype(float)
            rows.append(dict(method=mname, seed=sd, **det_metrics(y[te], (pg > 0.5).astype(int), pg),
                             **loc_metrics(N[stealth], Pst)))
            print(f"  {name} {mname} seed {sd} f1={rows[-1]['f1']:.4f}", flush=True)
            if sd == seeds[0]:
                scores[mname] = pg
                locs[mname] = pn
                hists[mname] = hist
    res = pd.DataFrame(rows)
    if store is not None:
        store.update(dict(scores=scores, locs=locs, te=te, y=y, N=N, kind=kind, mag=mag, timing=timing,
                          hists=hists, nodes=g.n, tau=tau, node_s=node_s, gs=gs, J=d["J"], thr=thr,
                          hours=d["hours"]))
    return res


def summarise(res):
    cols = [c for c in ["acc", "prec", "rec", "f1", "fpr", "auc", "loc_prec", "loc_rec", "loc_f1", "loc_exact"]
            if c in res]
    order = list(dict.fromkeys(res["method"]))
    m = res.groupby("method")[cols].mean().loc[order]
    s = res.groupby("method")[cols].std().fillna(0).loc[order]
    return m, s


if __name__ == "__main__":
    case = sys.argv[1]
    d = pickle.load(open(f"data/{case}.pkl", "rb"))
    store = {}
    res = run_case(case, d, store=store)
    res.to_csv(f"results/raw_{case}.csv", index=False)
    m, s = summarise(res)
    m.to_csv(f"results/mean_{case}.csv"); s.to_csv(f"results/std_{case}.csv")
    print(m.round(4).to_string())
    pickle.dump(store, open(f"results/store_{case}.pkl", "wb"))
