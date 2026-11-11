#!/usr/bin/env python
"""Phase 4: three-method dominant-process calls, consensus, concordance (kappa), pan-cohort positive controls, pooled spectra."""
import json, re, numpy as np, pandas as pd
from sklearn.metrics import cohen_kappa_score, roc_auc_score
from common import *; from sigs import *
rng = np.random.default_rng(SEED)
M = pd.read_parquet(f"{RES}/sbs96_all_samples.parquet"); n = M.sum(1)
cs = read_cbio(f"{STUDY}/data_clinical_sample.txt").set_index("SAMPLE_ID"); co = pd.read_csv(f"{RES}/cohort.tsv", sep="\t")
FAM35 = dict(FAMILY, SBS10c="POLE", SBS10d="POLE", SBS14="MMR", SBS20="MMR", SBS7a="UV")
def fam_of(name): return FAM35.get(name, "Other")
# ---------- Method 1: SPA
A = pd.read_csv(f"{RES}/spa/out/Assignment_Solution/Activities/Assignment_Solution_Activities.txt", sep="\t", index_col=0)
st = pd.read_csv(f"{RES}/spa/out/Assignment_Solution/Solution_Stats/Assignment_Solution_Samples_Stats.txt", sep="\t", index_col=0)
AF = A.T.groupby(fam_of).sum().T; A.to_csv(f"{RES}/spa_activities.tsv", sep="\t")
thr_of = lambda ids: np.where(n.reindex(ids).fillna(99).values >= 4, 0.60, 0.80)
m1 = pd.Series(np.where(st["Cosine Similarity"].reindex(AF.index).values >= thr_of(AF.index), AF.idxmax(1), "Indeterminate"), index=AF.index)
# ---------- Method 2: MSK-supplied COSMIC v2 contributions (mskcc/mutation-signatures) -- independent algorithm, substitutes DeepSig
V = pd.read_csv(f"{STUDY}/data_mutational_signatures_contribution_v2.txt", sep="\t", index_col=0).drop(columns=["NAME", "DESCRIPTION", "URL"]).T
V.columns = [int(re.search(r"Signature(\d+)", c).group(1)) for c in V.columns]
VF = V.T.groupby(lambda k: FAMILY_V2.get(k, "Other")).sum().T
m2 = VF.idxmax(1).where(VF.sum(1) > 0, "Indeterminate")
# ---------- Method 3: cosine classifier (panel-renormalised references), and NNLS exposures for all samples
ok = n[n >= 2].index
call, bc, C = cosine_classify(M.loc[ok].values, PAN, np.where(n[ok] >= 4, 0.60, 0.80)); m3 = pd.Series(call, index=ok)
X = fit_nnls(M.loc[ok].values, PAN); NF = fam_sum(X); NF.index = ok; NFfrac = NF.div(NF.sum(1), axis=0)
calls = pd.DataFrame({"M1_SPA": m1, "M2_MSKv2": m2, "M3_cosine": m3, "M3_best_cosine": pd.Series(bc, index=ok)}).reindex(cs.index)
calls["N_SNV"] = n
def consensus(r):
    v = [x for x in r if isinstance(x, str) and x != "Indeterminate"]
    if not v: return "NO_CALL"
    vc = pd.Series(v).value_counts()
    return vc.index[0] if vc.iloc[0] >= 2 else "AMBIGUOUS"
calls["CONSENSUS"] = calls[["M1_SPA", "M2_MSKv2", "M3_cosine"]].apply(consensus, axis=1)
calls.to_csv(f"{RES}/calls_all_samples.tsv.gz", sep="\t")
out = {}
# ---------- concordance
def kappa(a, b, idx, drop_indet):
    d = calls.loc[idx, [a, b]].dropna()
    if drop_indet: d = d[(d[a] != "Indeterminate") & (d[b] != "Indeterminate")]
    return {"n": int(len(d)), "kappa": float(cohen_kappa_score(d[a], d[b])) if len(d) > 5 and d.nunique().max() > 1 else None, "agreement": float((d[a] == d[b]).mean()) if len(d) else None}
basket = co[co.N_SNV >= 2].SAMPLE_ID; spa_ids = [i for i in m1.index if not i.startswith("POOLED")]
pairs = [("M1_SPA", "M2_MSKv2"), ("M1_SPA", "M3_cosine"), ("M2_MSKv2", "M3_cosine")]
out["kappa"] = {s: {f"{a}~{b}": {"called_only": kappa(a, b, idx, True), "incl_indeterminate": kappa(a, b, idx, False)} for a, b in pairs}
                for s, idx in {"basket_ge2": basket, "basket_ge4": co[co.N_SNV >= 4].SAMPLE_ID, "basket_plus_controls_SPA_set": spa_ids,
                               "pan_cohort_ge4_(M2~M3 only)": n[n >= 4].index, "pan_cohort_ge20_(M2~M3 only)": n[n >= 20].index}.items()}
km = [v["called_only"]["kappa"] for v in out["kappa"]["basket_ge2"].values()]; out["kappa_gate_basket_mean"] = float(np.mean(km)); out["kappa_gate_pass"] = bool(np.mean(km) >= 0.6)
# ---------- basket per-sample table
b = co.merge(calls.drop(columns="N_SNV"), left_on="SAMPLE_ID", right_index=True, how="left")
b = b.merge(NFfrac.add_prefix("NNLS_"), left_on="SAMPLE_ID", right_index=True, how="left")
skin = {"Skin", "Scalp", "Forehead", "Neck", "Parotid"}
b["ASI_SITE"] = np.where(b.GROUP != "ASI", "", np.where(b.PRIMARY_SITE.isin(skin), "Cutaneous/H&N", np.where(b.PRIMARY_SITE == "Breast", "Breast", "Visceral/soft tissue/other")))
b.to_csv(f"{RES}/basket_sample_calls.tsv", sep="\t", index=False)
an = b[b.FIT_CLASS != "DROP"]
out["consensus_by_group"] = pd.crosstab(an.GROUP, an.CONSENSUS).to_dict("index")
out["consensus_by_ASI_site"] = pd.crosstab(an[an.GROUP == "ASI"].ASI_SITE, an[an.GROUP == "ASI"].CONSENSUS).to_dict("index")
out["method_calls_basket"] = {m: an[m].value_counts(dropna=False).to_dict() for m in ["M1_SPA", "M2_MSKv2", "M3_cosine"]}
# ---------- pan-cohort positive controls (>=4 SNVs)
g4 = calls[calls.N_SNV >= 4].join(cs[["ONCOTREE_CODE", "CANCER_TYPE", "MSI_TYPE", "MSI_SCORE"]]).join(NFfrac.add_prefix("NNLS_"))
def sens_spec(pos, neg, fam):
    r = {}
    for m in ["M2_MSKv2", "M3_cosine", "CONSENSUS", "M1_SPA"]:
        p, q = g4.loc[pos, m].dropna(), g4.loc[neg, m].dropna()
        if m == "M1_SPA" and (len(p) < 20 or len(q) < 20): continue
        r[m] = {"n_pos": int(len(p)), "sens": float((p == fam).mean()), "n_neg": int(len(q)), "spec": float((q != fam).mean())}
    return r
skinlike = g4.CANCER_TYPE.isin(["Melanoma", "Skin Cancer, Non-Melanoma"])
out["controls"] = {
    "UV_in_cutaneous_melanoma_vs_non_skin": sens_spec(g4.ONCOTREE_CODE.eq("SKCM"), ~skinlike & ~g4.CANCER_TYPE.str.contains("Unknown Primary|Head and Neck|Sarcoma", na=False), "UV"),
    "MMR_in_MSI_instable_vs_stable": sens_spec(g4.MSI_TYPE.eq("Instable"), g4.MSI_TYPE.eq("Stable"), "MMR"),
    "APOBEC_in_bladder_vs_colorectal": sens_spec(g4.CANCER_TYPE.eq("Bladder Cancer"), g4.CANCER_TYPE.eq("Colorectal Cancer"), "APOBEC"),
    "Tobacco_in_NSCLC_vs_colorectal": sens_spec(g4.CANCER_TYPE.eq("Non-Small Cell Lung Cancer"), g4.CANCER_TYPE.eq("Colorectal Cancer"), "Tobacco")}
d = g4[g4.MSI_TYPE.isin(["Instable", "Stable"])]
out["AUC_MMRfraction_vs_MSI"] = {f">={k}": {"n": int((d.N_SNV >= k).sum()), "n_MSI": int(((d.N_SNV >= k) & d.MSI_TYPE.eq("Instable")).sum()),
                                 "auc": float(roc_auc_score(d[d.N_SNV >= k].MSI_TYPE.eq("Instable"), d[d.N_SNV >= k].NNLS_MMR))} for k in [4, 10, 20]}
d = g4[g4.ONCOTREE_CODE.eq("SKCM") | (~skinlike)]
out["AUC_UVfraction_SKCM_vs_rest"] = {f">={k}": float(roc_auc_score(d[d.N_SNV >= k].ONCOTREE_CODE.eq("SKCM"), d[d.N_SNV >= k].NNLS_UV)) for k in [4, 10, 20]}
# ---------- pooled (pseudo-bulk) spectra with patient bootstrap
pooled = {}; boots = {}
strata = {g: b[b.GROUP == g].SAMPLE_ID.values for g in ["ASI", "THYM", "ACC"]}
for s in ["Cutaneous/H&N", "Breast", "Visceral/soft tissue/other"]: strata[f"ASI: {s}"] = b[b.ASI_SITE == s].SAMPLE_ID.values
strata["THYM: thymoma"] = b[b.ONCOTREE_CODE == "THYM"].SAMPLE_ID.values; strata["THYM: thymic carcinoma"] = b[b.ONCOTREE_CODE == "THYC"].SAMPLE_ID.values
strata["ASI (excl. hypermutator)"] = b[(b.GROUP == "ASI") & ~b.HYPERMUTATOR].SAMPLE_ID.values
spec = {}
for k, ids in strata.items():
    m = M.loc[ids]; tot = m.sum(0).values; spec[k] = tot
    f = fam_sum(fit_nnls(tot[None, :], PAN)).iloc[0]; f = f / f.sum()
    bs = np.array([M.loc[rng.choice(ids, len(ids))].sum(0).values for _ in range(1000)])
    bf = fam_sum(fit_nnls(bs, PAN)); bf = bf.div(bf.sum(1), axis=0)
    rec = PAN.values @ fit_nnls(tot[None, :], PAN)[0]
    pooled[k] = {"n_patients": int(len(ids)), "n_snv": int(tot.sum()), "recon_cosine": cosine(rec, tot),
                 "families": {fam: {"frac": float(f[fam]), "lo": float(bf[fam].quantile(.025)), "hi": float(bf[fam].quantile(.975))} for fam in FAMILIES},
                 "SPA_pooled": AF.loc[f"POOLED_{k}"].div(AF.loc[f"POOLED_{k}"].sum()).round(3).to_dict() if f"POOLED_{k}" in AF.index else None,
                 "SPA_pooled_sigs": A.loc[f"POOLED_{k}"][A.loc[f"POOLED_{k}"] > 0].to_dict() if f"POOLED_{k}" in A.index else None}
    boots[k] = bf
pd.DataFrame(spec, index=CHANNELS).to_csv(f"{RES}/pooled_spectra.tsv", sep="\t")
pd.concat(boots).to_parquet(f"{RES}/pooled_bootstrap.parquet")
# bootstrap difference UV fraction cutaneous vs non-cutaneous ASI
du = boots["ASI: Cutaneous/H&N"].UV.values - boots["ASI: Visceral/soft tissue/other"].UV.values
out["UV_frac_diff_cutaneous_vs_visceral_ASI"] = {"diff": float(pooled["ASI: Cutaneous/H&N"]["families"]["UV"]["frac"] - pooled["ASI: Visceral/soft tissue/other"]["families"]["UV"]["frac"]),
                                                  "ci": [float(np.quantile(du, .025)), float(np.quantile(du, .975))]}
out["pooled"] = pooled
json.dump(out, open(f"{RES}/consensus_summary.json", "w"), indent=2)
print(json.dumps({k: v for k, v in out.items() if k != "pooled"}, indent=1))
for k, v in pooled.items(): print(k, v["n_patients"], v["n_snv"], round(v["recon_cosine"], 3), {f: (round(x["frac"], 2), round(x["lo"], 2), round(x["hi"], 2)) for f, x in v["families"].items() if x["hi"] > 0.03}, v["SPA_pooled_sigs"])
