#!/usr/bin/env python
"""Phase 3 validation gate: synthetic tumours with known exposures drawn inside the panel territory.
Custom multinomial simulator (ground-truth exposures are required, which SigProfilerSimulator does not provide)."""
import json, numpy as np, pandas as pd
from common import *; from sigs import *
rng = np.random.default_rng(SEED)
SCEN = {"UV": ["SBS7a", "SBS7b"], "Clock": [], "APOBEC": ["SBS2", "SBS13"], "HRD": ["SBS3"]}
def draw_expo(scen):
    e = pd.Series(0.0, index=REF_SIGS); f = rng.uniform(0.6, 0.9) if SCEN[scen] else 0.0
    if SCEN[scen]:
        w = rng.dirichlet(np.ones(len(SCEN[scen])) * 3); e[SCEN[scen]] = f * w
    c1 = rng.uniform(0.1, 0.5); e["SBS1"] = (1 - f) * c1; e["SBS5"] = (1 - f) * (1 - c1)
    return e.values
rows = []; NS = [3, 5, 10, 20, 50, 100, 200, 500, 2000]
for N in NS:
    scen = np.repeat(list(SCEN), 50); T = np.array([draw_expo(s) for s in scen])          # 200 tumours
    prof = T @ PAN.values.T; prof /= prof.sum(1, keepdims=True)
    M = np.array([rng.multinomial(N, p) for p in prof])
    truthF = fam_sum(T).values
    for name, S in [("uncorrected_genome", GEN), ("exome_SPA", EXO), ("panel_opportunity", PAN)]:
        X = fit_nnls(M, S); Xn = X / N; F = fam_sum(Xn).values
        cs_sig = [cosine(a, b) for a, b in zip(T, Xn)]; cs_fam = [cosine(a, b) for a, b in zip(truthF, F)]
        dom = np.array(FAMILIES)[F.argmax(1)]
        thr = np.full(len(M), 0.80 if N < 4 else 0.60); call, bc, _ = cosine_classify(M, S, thr)
        called = call != "Indeterminate"
        for sc in list(SCEN) + ["ALL"]:
            m = (scen == sc) if sc != "ALL" else np.ones(len(scen), bool)
            rows.append(dict(N=N, correction=name, scenario=sc, cos_signature=np.mean(np.array(cs_sig)[m]), cos_family=np.mean(np.array(cs_fam)[m]),
                             nnls_dominant_acc=np.mean(dom[m] == scen[m]), cosine_call_rate=called[m].mean(),
                             cosine_acc_when_called=(call[m & called] == scen[m & called]).mean() if (m & called).any() else np.nan))
R = pd.DataFrame(rows); R.to_csv(f"{RES}/simulation_gate.tsv", sep="\t", index=False)
A = R[R.scenario == "ALL"].pivot(index="N", columns="correction", values="cos_family").round(3)
print(A); print(R[R.scenario == "ALL"].pivot(index="N", columns="correction", values="cos_signature").round(3))
print(R[(R.correction == "panel_opportunity")].pivot(index="N", columns="scenario", values="nnls_dominant_acc").round(2))
print(R[(R.correction == "panel_opportunity")].pivot(index="N", columns="scenario", values="cosine_acc_when_called").round(2))
print(R[(R.correction == "panel_opportunity")].pivot(index="N", columns="scenario", values="cosine_call_rate").round(2))
gate = {c: int(A.index[(A[c] >= 0.90).values].min()) if (A[c] >= 0.90).any() else None for c in A.columns}
json.dump({"min_N_family_cosine_ge_0.90": gate, "family_cosine_by_N": A.to_dict()}, open(f"{RES}/simulation_gate.json", "w"), indent=2); print(gate)
