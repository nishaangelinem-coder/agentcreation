#!/usr/bin/env python
"""Phase 4, Method 1: SigProfilerAssignment cosmic_fit (COSMIC v3.5, GRCh37, exome=True) on basket samples (>=2 SNVs),
pooled group spectra and pan-cohort positive controls."""
import os, sys, time, numpy as np, pandas as pd
from common import *
from SigProfilerAssignment import Analyzer
M = pd.read_parquet(f"{RES}/sbs96_all_samples.parquet"); co = pd.read_csv(f"{RES}/cohort.tsv", sep="\t")
cs = read_cbio(f"{STUDY}/data_clinical_sample.txt").set_index("SAMPLE_ID")
n = M.sum(1); rng = np.random.default_rng(SEED)
basket = co[co.N_SNV >= 2].SAMPLE_ID.tolist()
ctrl = {}
for lab, mask in {"CTRL_UV_cutaneous_melanoma": cs.ONCOTREE_CODE.isin(["SKCM"]), "CTRL_Tobacco_lung": cs.CANCER_TYPE.eq("Non-Small Cell Lung Cancer"),
                  "CTRL_APOBEC_bladder": cs.CANCER_TYPE.eq("Bladder Cancer"), "CTRL_MSI_instable": cs.MSI_TYPE.eq("Instable"),
                  "CTRL_MSS_colorectal": cs.CANCER_TYPE.eq("Colorectal Cancer") & cs.MSI_TYPE.eq("Stable")}.items():
    ids = [i for i in cs.index[mask] if n[i] >= 4]; ctrl[lab] = list(rng.choice(ids, min(300, len(ids)), replace=False))
pd.DataFrame([(k, i) for k, v in ctrl.items() for i in v], columns=["set", "SAMPLE_ID"]).to_csv(f"{RES}/control_sets.tsv", sep="\t", index=False)
ids = list(dict.fromkeys(basket + [i for v in ctrl.values() for i in v]))
X = M.loc[ids].T
for g in ["ASI", "THYM", "ACC"]: X[f"POOLED_{g}"] = M.loc[co[co.GROUP == g].SAMPLE_ID].sum(0)
X.index.name = "MutationType"; os.makedirs(f"{RES}/spa", exist_ok=True); X.to_csv(f"{RES}/spa/input_SBS96.txt", sep="\t")
t = time.time()
Analyzer.cosmic_fit(samples=f"{RES}/spa/input_SBS96.txt", output=f"{RES}/spa/out", input_type="matrix", context_type="96", genome_build="GRCh37",
                    cosmic_version=3.5, exome=True, make_plots=False, sample_reconstruction_plots=None, export_probabilities=True, cpu=4)
print("SPA done", X.shape, round(time.time() - t), "s")
