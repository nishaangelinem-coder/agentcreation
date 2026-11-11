#!/usr/bin/env python
"""Phase 5.4: constrained de novo extraction (SigProfilerExtractor, 250 NMF replicates) on basket samples whose best
cosine similarity to every known reference process is < 0.75; novelty criteria: stability >= 0.8 and max cosine to COSMIC v3.5 < 0.6."""
import os, json, glob, numpy as np, pandas as pd
from common import *; from sigs import *
from SigProfilerExtractor import sigpro as sig
if __name__ == "__main__":
    b = pd.read_csv(f"{RES}/basket_sample_calls.tsv", sep="\t"); M = pd.read_parquet(f"{RES}/sbs96_all_samples.parquet")
    ids = b[(b.FIT_CLASS != "DROP") & (b.M3_best_cosine < 0.75)].SAMPLE_ID
    X = M.loc[ids].T; X.index.name = "MutationType"; os.makedirs(f"{RES}/denovo", exist_ok=True); X.to_csv(f"{RES}/denovo/input.txt", sep="\t")
    sig.sigProfilerExtractor("matrix", f"{RES}/denovo/out", f"{RES}/denovo/input.txt", reference_genome="GRCh37", minimum_signatures=1, maximum_signatures=4,
                             nmf_replicates=250, cpu=4, seeds=f"{RES}/denovo/seeds.txt" if os.path.exists(f"{RES}/denovo/seeds.txt") else "random", make_decomposition_plots=False, exome=True)
    st = pd.read_csv(f"{RES}/denovo/out/SBS96/All_solutions_stat.csv") if os.path.exists(f"{RES}/denovo/out/SBS96/All_solutions_stat.csv") else None
    res = {"n_samples": int(len(ids)), "n_mutations": int(X.values.sum())}
    if st is not None: res["solutions"] = st.to_dict("records")
    cand = []
    for f in sorted(glob.glob(f"{RES}/denovo/out/SBS96/All_Solutions/SBS96_*_Signatures/Signatures/SBS96_S*_Signatures.txt")):
        S = pd.read_csv(f, sep="\t", index_col=0).reindex(CHANNELS)
        for c in S.columns:
            cs_ = {k: cosine(S[c], EXO_ALL[k]) for k in EXO_ALL.columns}; best = max(cs_, key=cs_.get)
            cand.append({"solution": os.path.basename(f), "signature": c, "best_COSMIC": best, "cosine": round(cs_[best], 3)})
    res["signatures_vs_COSMIC"] = cand
    json.dump(res, open(f"{RES}/denovo_summary.json", "w"), indent=1, default=str); print(json.dumps(res, indent=1, default=str))
