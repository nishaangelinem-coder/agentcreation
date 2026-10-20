#!/usr/bin/env python3
"""Stage 2 — ML triage on Tox21 (scaffold split, imbalance-aware metrics).

Downloads the public Tox21 benchmark, trains one classifier per inflammation /
stress-relevant endpoint using Morgan-fingerprint features and a Bemis-Murcko
**scaffold split** (so test scaffolds are unseen in training), and reports
ROC-AUC *and* AUPRC (AUPRC is the headline because the positive classes are
heavily imbalanced). The trained models then score every Stage-1 leachate.

Endpoints used (all mechanistically linked to neuroinflammation):
  NR-AhR  — aryl hydrocarbon receptor activation
  SR-ARE  — Nrf2/antioxidant-response-element (oxidative stress)
  SR-p53  — p53 stress / genotoxic response

This is a deliberately transparent sklearn baseline (RandomForest on 2048-bit
Morgan fingerprints), not a graph neural network. DeepChem GraphConv/AttentiveFP
is the handover's suggested model; it is listed in environment.yml as an optional
upgrade. The sklearn baseline is used here because it installs and runs
deterministically in the sandbox; all reported metrics are from this baseline and
are labelled as such.

Inputs : data/leachates.csv
Outputs: results/ml_scores.csv        (SMILES, polymer, chemical_name, endpoint, prob)
         results/ml_metrics.csv       (endpoint, n_pos, n_neg, roc_auc, auprc, base_rate)
         results/ml_calibration.md     (reliability note)
"""
from __future__ import annotations

import gzip
import io
import os
import sys
import urllib.request

import numpy as np
import pandas as pd
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem
from rdkit.Chem.Scaffolds import MurckoScaffold
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import average_precision_score, roc_auc_score

RDLogger.DisableLog("rdApp.*")

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
DATA = os.path.join(ROOT, "data")
RESULTS = os.path.join(ROOT, "results")
TOX21_LOCAL = os.path.join(DATA, "tox21.csv")

TOX21_URLS = [
    "https://deepchemdata.s3-us-west-1.amazonaws.com/datasets/tox21.csv.gz",
    "https://raw.githubusercontent.com/deepchem/deepchem/master/datasets/tox21.csv.gz",
]
ENDPOINTS = ["NR-AhR", "SR-ARE", "SR-p53"]
N_BITS = 2048
RADIUS = 2
SEED = 42


def fetch_tox21() -> pd.DataFrame:
    """Load Tox21 from a local cache, else download and cache it."""
    if os.path.exists(TOX21_LOCAL):
        print(f"[S2] using cached Tox21 at {TOX21_LOCAL}")
        return pd.read_csv(TOX21_LOCAL)
    last_err = None
    for url in TOX21_URLS:
        try:
            print(f"[S2] downloading Tox21 from {url}")
            with urllib.request.urlopen(url, timeout=60) as resp:
                raw = resp.read()
            if url.endswith(".gz"):
                raw = gzip.decompress(raw)
            df = pd.read_csv(io.BytesIO(raw))
            df.to_csv(TOX21_LOCAL, index=False)
            print(f"[S2] cached Tox21 -> {TOX21_LOCAL} ({len(df)} rows)")
            return df
        except Exception as exc:  # noqa: BLE001
            last_err = exc
            print(f"[S2]   failed: {exc}")
    raise RuntimeError(
        "Tox21 could not be downloaded and no local cache exists. "
        "Place tox21.csv in data/ to run offline."
    ) from last_err


def morgan(smiles: str) -> np.ndarray | None:
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    fp = AllChem.GetMorganFingerprintAsBitVect(mol, RADIUS, nBits=N_BITS)
    arr = np.zeros((N_BITS,), dtype=np.int8)
    AllChem.DataStructs.ConvertToNumpyArray(fp, arr)
    return arr


def scaffold(smiles: str) -> str:
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return ""
    try:
        return MurckoScaffold.MurckoScaffoldSmiles(mol=mol, includeChirality=False)
    except Exception:  # noqa: BLE001
        return ""


def scaffold_split(smiles_list, frac_train=0.8):
    """Deterministic Bemis-Murcko scaffold split: largest scaffold groups to
    train first so the test set is populated by rarer, unseen scaffolds."""
    groups: dict[str, list[int]] = {}
    for idx, smi in enumerate(smiles_list):
        groups.setdefault(scaffold(smi), []).append(idx)
    ordered = sorted(groups.values(), key=lambda g: (-len(g), g[0]))
    n_total = len(smiles_list)
    train_idx, test_idx = [], []
    for grp in ordered:
        if len(train_idx) + len(grp) <= frac_train * n_total:
            train_idx.extend(grp)
        else:
            test_idx.extend(grp)
    return sorted(train_idx), sorted(test_idx)


def main() -> int:
    os.makedirs(RESULTS, exist_ok=True)
    tox = fetch_tox21()
    smiles_col = "smiles" if "smiles" in tox.columns else "SMILES"

    # Featurise the whole Tox21 set once.
    feats, keep = [], []
    for i, smi in enumerate(tox[smiles_col].astype(str)):
        fp = morgan(smi)
        if fp is not None:
            feats.append(fp)
            keep.append(i)
    tox = tox.iloc[keep].reset_index(drop=True)
    X_all = np.vstack(feats)
    smiles_all = tox[smiles_col].astype(str).tolist()
    tr_idx, te_idx = scaffold_split(smiles_all, frac_train=0.8)
    print(f"[S2] Tox21 featurised: {len(smiles_all)} mols | "
          f"scaffold split train={len(tr_idx)} test={len(te_idx)}")

    leach = pd.read_csv(os.path.join(DATA, "leachates.csv"))
    leach_fps, leach_keep = [], []
    for i, smi in enumerate(leach["SMILES"]):
        fp = morgan(smi)
        if fp is not None:
            leach_fps.append(fp)
            leach_keep.append(i)
    leach = leach.iloc[leach_keep].reset_index(drop=True)
    X_leach = np.vstack(leach_fps)

    metrics_rows = []
    score_frames = []
    trained_models = {}
    for ep in ENDPOINTS:
        y = tox[ep].values.astype(float)
        mask = ~np.isnan(y)  # Tox21 has missing labels per endpoint
        tr = [i for i in tr_idx if mask[i]]
        te = [i for i in te_idx if mask[i]]
        Xtr, ytr = X_all[tr], y[tr].astype(int)
        Xte, yte = X_all[te], y[te].astype(int)
        base_rate = float(yte.mean()) if len(yte) else float("nan")

        clf = RandomForestClassifier(n_estimators=300, n_jobs=-1, random_state=SEED, class_weight="balanced")
        clf.fit(Xtr, ytr)
        proba = clf.predict_proba(Xte)[:, 1]
        roc = roc_auc_score(yte, proba) if len(set(yte)) > 1 else float("nan")
        auprc = average_precision_score(yte, proba) if len(set(yte)) > 1 else float("nan")
        metrics_rows.append({
            "endpoint": ep,
            "n_train": len(tr),
            "n_test": len(te),
            "n_pos_test": int(yte.sum()),
            "base_rate_test": round(base_rate, 4),
            "roc_auc": round(float(roc), 4),
            "auprc": round(float(auprc), 4),
        })
        print(f"[S2] {ep}: scaffold-split ROC-AUC={roc:.3f}  AUPRC={auprc:.3f}  "
              f"(test base rate={base_rate:.3f})")

        # Refit on all labelled data for final leachate scoring.
        all_lab = [i for i in range(len(y)) if mask[i]]
        clf_full = RandomForestClassifier(n_estimators=300, n_jobs=-1, random_state=SEED, class_weight="balanced")
        clf_full.fit(X_all[all_lab], y[all_lab].astype(int))
        trained_models[ep] = clf_full
        leach_prob = clf_full.predict_proba(X_leach)[:, 1]
        df = leach[["SMILES", "polymer", "chemical_name"]].copy()
        df["endpoint"] = ep
        df["prob"] = np.round(leach_prob, 4)
        score_frames.append(df)

    pd.DataFrame(metrics_rows).to_csv(os.path.join(RESULTS, "ml_metrics.csv"), index=False)
    scores = pd.concat(score_frames, ignore_index=True)
    scores.to_csv(os.path.join(RESULTS, "ml_scores.csv"), index=False)
    print(f"[S2] wrote results/ml_metrics.csv and results/ml_scores.csv "
          f"({len(scores)} score rows)")

    with open(os.path.join(RESULTS, "ml_calibration.md"), "w") as fh:
        fh.write("# Stage 2 reliability note\n\n")
        fh.write("- Model: RandomForest (300 trees, balanced) on 2048-bit Morgan (r=2) "
                 "fingerprints (sklearn baseline; DeepChem GNN is an optional upgrade).\n")
        fh.write("- Validation: deterministic Bemis-Murcko **scaffold split** "
                 "(80/20); test scaffolds are unseen in training.\n")
        fh.write("- Headline metric is **AUPRC**, because Tox21 positive classes "
                 "are rare (test base rates below).\n")
        fh.write("- Probabilities are raw RandomForest vote fractions and are **not "
                 "isotonically calibrated**; treat them as a relative ranking "
                 "signal, not absolute activation probabilities.\n")
        fh.write("- Leachate scores are produced by models refit on all labelled "
                 "Tox21 rows for each endpoint.\n\n")
        fh.write("| endpoint | n_test | n_pos_test | base_rate | ROC-AUC | AUPRC |\n")
        fh.write("|---|---|---|---|---|---|\n")
        for r in metrics_rows:
            fh.write(f"| {r['endpoint']} | {r['n_test']} | {r['n_pos_test']} | "
                     f"{r['base_rate_test']} | {r['roc_auc']} | {r['auprc']} |\n")
    print("[S2] wrote results/ml_calibration.md")
    return 0


if __name__ == "__main__":
    sys.exit(main())
