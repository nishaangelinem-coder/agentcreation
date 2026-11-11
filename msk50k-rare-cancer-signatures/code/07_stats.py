#!/usr/bin/env python
"""Phase 5: landscape test, survival, gene-signature enrichment, biomarker checks (basket cohort)."""
import json, warnings, numpy as np, pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests
from lifelines import KaplanMeierFitter, CoxPHFitter
from lifelines.statistics import multivariate_logrank_test, logrank_test, proportional_hazard_test
from common import *; from sigs import *
warnings.filterwarnings("ignore"); rng = np.random.default_rng(SEED); out = {}
b = pd.read_csv(f"{RES}/basket_sample_calls.tsv", sep="\t"); an = b[b.FIT_CLASS != "DROP"].copy()
# ---- 1 landscape
an["CAT"] = an.CONSENSUS.where(an.CONSENSUS.isin(["UV", "Clock", "MMR", "AMBIGUOUS"]), "Other called")
T = pd.crosstab(an.GROUP, an.CAT); chi2, p, dof, exp = stats.chi2_contingency(T)
resid = (T - exp) / np.sqrt(exp * np.outer(1 - T.sum(1) / T.values.sum(), 1 - T.sum(0) / T.values.sum()))
g, c = an.GROUP.values, an.CAT.values; perm = np.mean([stats.chi2_contingency(pd.crosstab(g, rng.permutation(c)))[0] >= chi2 for _ in range(10000)])
out["landscape"] = {"table": T.to_dict("index"), "chi2": chi2, "dof": int(dof), "p_asymptotic": p, "p_permutation": float(max(perm, 1e-4)), "min_expected": float(exp.min()),
                    "adj_std_residuals": resid.round(2).to_dict("index")}
a = an[an.GROUP == "ASI"]; t2 = pd.crosstab(a.ASI_SITE == "Cutaneous/H&N", a.CONSENSUS == "UV"); o, pf = stats.fisher_exact(t2)
out["ASI_UV_by_site_fisher"] = {"table": t2.values.tolist(), "odds_ratio": o, "p": pf}
# ---- 2 survival
s = b.dropna(subset=["OS_MONTHS", "OS_STATUS"]).copy(); s["E"] = s.OS_STATUS.str.contains("DECEASED").astype(int); s["T"] = s.OS_MONTHS.clip(lower=0.01)
s["log10TMB"] = np.log10(pd.to_numeric(s.TMB_SCORE, errors="coerce").replace([np.inf, -np.inf], np.nan).clip(lower=0) + 1); s["AGE"] = pd.to_numeric(s.AGE_AT_DX, errors="coerce"); s["MSI"] = s.MSI_SCORE
surv = {"n": int(len(s)), "events": int(s.E.sum()), "by_group": {}}
for gname, d in s.groupby("GROUP"):
    km = KaplanMeierFitter().fit(d["T"], d.E); surv["by_group"][gname] = {"n": int(len(d)), "events": int(d.E.sum()), "median_OS_months": float(km.median_survival_time_),
        "median_followup": float(d["T"].median()), "OS_24m": float(km.predict(24))}
lr = multivariate_logrank_test(s["T"], s.GROUP, s.E); surv["logrank_groups_p"] = float(lr.p_value)
def cox(d, covs, label, strata=None):
    d = d[covs + ["T", "E"] + ([strata] if strata else [])].dropna(); 
    if d.E.sum() < 10: return {"label": label, "skipped": "fewer than 10 events", "n": int(len(d)), "events": int(d.E.sum())}
    m = CoxPHFitter(penalizer=0.0).fit(d, "T", "E", strata=strata); ph = proportional_hazard_test(m, d, time_transform="rank")
    r = m.summary[["exp(coef)", "exp(coef) lower 95%", "exp(coef) upper 95%", "p"]].round(4)
    return {"label": label, "n": int(len(d)), "events": int(d.E.sum()), "concordance": float(m.concordance_index_), "hr": r.to_dict("index"),
            "schoenfeld_p": {k: float(v) for k, v in zip(ph.summary.index, ph.summary.p)}, "ph_ok": bool((ph.summary.p > 0.05).all())}
s["G_THYM"] = (s.GROUP == "THYM").astype(int); s["G_ACC"] = (s.GROUP == "ACC").astype(int)
s["UV"] = (s.CONSENSUS == "UV").astype(int); s["CUT"] = (s.ASI_SITE == "Cutaneous/H&N").astype(int); s["MET"] = (s.SAMPLE_TYPE == "Metastasis").astype(int)
base = ["AGE", "log10TMB", "N_SNV", "MSI"]
surv["cox"] = [cox(s, base + ["G_THYM", "G_ACC"], "A: all basket patients, protocol covariates (ref ASI)"),
               cox(s, base + ["MET"], "A2: as A, stratified by group", strata="GROUP"),
               cox(s[(s.GROUP == "ASI") & (s.FIT_CLASS != "DROP")], ["UV", "AGE", "log10TMB", "MET"], "B: ASI analysable, consensus UV-dominant"),
               cox(s[(s.GROUP == "ASI") & (s.FIT_CLASS != "DROP")], ["UV", "CUT", "AGE", "log10TMB", "MET"], "B2: as B + cutaneous/H&N site"),
               cox(s[(s.GROUP == "ASI") & (s.FIT_CLASS == "FULL")], ["UV", "AGE", "log10TMB", "MET"], "B-sens: ASI FULL class only (SPARSE excluded)"),
               cox(s[(s.GROUP == "ASI")], ["CUT", "AGE", "log10TMB", "MET"], "C: ASI all patients, cutaneous/H&N site"),
               cox(s[(s.GROUP == "ACC")], ["AGE", "log10TMB", "MET"], "D: ACC all patients"),
               cox(s[(s.GROUP == "THYM")], ["AGE", "log10TMB", "MET"], "E: THYM all patients")]
d = s[(s.GROUP == "ASI") & (s.FIT_CLASS != "DROP")]; l = logrank_test(d[d.UV == 1]["T"], d[d.UV == 0]["T"], d[d.UV == 1].E, d[d.UV == 0].E)
surv["ASI_UV_logrank"] = {"n_UV": int(d.UV.sum()), "events_UV": int(d[d.UV == 1].E.sum()), "n_nonUV": int((1 - d.UV).sum()), "p": float(l.p_value),
                          "median_UV": float(KaplanMeierFitter().fit(d[d.UV == 1]["T"], d[d.UV == 1].E).median_survival_time_),
                          "median_nonUV": float(KaplanMeierFitter().fit(d[d.UV == 0]["T"], d[d.UV == 0].E).median_survival_time_)}
out["survival"] = surv; s.to_csv(f"{RES}/survival_table.tsv", sep="\t", index=False)
# ---- 3 gene-signature enrichment (per-mutation maximum-posterior family from panel-frame NNLS exposures)
snv = pd.read_parquet(f"{RES}/snv_annotated.parquet"); snv = snv[snv.Tumor_Sample_Barcode.isin(an.SAMPLE_ID)].copy()
M = pd.read_parquet(f"{RES}/sbs96_all_samples.parquet").loc[an.SAMPLE_ID]; X = pd.DataFrame(fit_nnls(M.values, PAN), index=M.index, columns=REF_SIGS)
ci = snv.channel.map(CH_IDX).values; post = X.loc[snv.Tumor_Sample_Barcode].values * PAN.values[ci]           # (n_mut, k)
pf_ = fam_sum(post); snv["FAM"] = np.where(pf_.sum(1) > 0, pf_.idxmax(1), "NA"); snv["GROUP"] = snv.Tumor_Sample_Barcode.map(an.set_index("SAMPLE_ID").GROUP)
snv.drop(columns=["tri"]).to_csv(f"{RES}/basket_mutation_attribution.tsv", sep="\t", index=False)
rows = []; gc = snv.Hugo_Symbol.value_counts()
for gene in gc[gc >= 3].index:
    for fam in FAMILIES:
        ing, isf = snv.Hugo_Symbol == gene, snv.FAM == fam; k = int((ing & isf).sum())
        if isf.sum() < 10: continue
        orr, pv = stats.fisher_exact([[k, int((ing & ~isf).sum())], [int((~ing & isf).sum()), int((~ing & ~isf).sum())]])
        rows.append(dict(gene=gene, family=fam, n_gene=int(ing.sum()), n_gene_family=k, odds_ratio=orr, p=pv))
E_ = pd.DataFrame(rows); E_["q"] = multipletests(E_.p, method="fdr_bh")[1]; E_ = E_.sort_values("p"); E_.to_csv(f"{RES}/gene_signature_enrichment.tsv", sep="\t", index=False)
out["enrichment"] = {"n_mutations": int(len(snv)), "n_tests": int(len(E_)), "n_sig_q05": int((E_.q < 0.05).sum()), "top": E_.head(8).round(4).to_dict("records"),
                     "mutations_by_family": snv.FAM.value_counts().to_dict()}
# most frequently mutated genes per group (all nonsynonymous incl. indels) for context
maf = pd.read_csv(f"{STUDY}/data_mutations.txt", sep="\t", usecols=["Hugo_Symbol", "Tumor_Sample_Barcode", "Variant_Classification", "HGVSp_Short"], low_memory=False)
maf = maf[maf.Tumor_Sample_Barcode.isin(b.SAMPLE_ID) & ~maf.Variant_Classification.isin(["Silent", "Intron", "3'UTR", "5'UTR", "3'Flank"])]
maf["GROUP"] = maf.Tumor_Sample_Barcode.map(b.set_index("SAMPLE_ID").GROUP)
out["top_genes"] = {g: {k: [int(v), round(100 * v / (b.GROUP == g).sum(), 1)] for k, v in d.drop_duplicates(["Tumor_Sample_Barcode", "Hugo_Symbol"]).Hugo_Symbol.value_counts().head(8).items()} for g, d in maf.groupby("GROUP")}
# gene vs consensus-UV within ASI (sample level)
asi = an[an.GROUP == "ASI"]; mm = maf[maf.Tumor_Sample_Barcode.isin(asi.SAMPLE_ID)].drop_duplicates(["Tumor_Sample_Barcode", "Hugo_Symbol"]); r2 = []
for gene, k in mm.Hugo_Symbol.value_counts().items():
    if k < 5: break
    has = asi.SAMPLE_ID.isin(mm[mm.Hugo_Symbol == gene].Tumor_Sample_Barcode); uv = asi.CONSENSUS == "UV"
    orr, pv = stats.fisher_exact(pd.crosstab(has, uv)); r2.append(dict(gene=gene, n_mut=int(has.sum()), mut_in_UV=int((has & uv).sum()), n_UV=int(uv.sum()), odds_ratio=orr, p=pv))
R2 = pd.DataFrame(r2); R2["q"] = multipletests(R2.p, method="fdr_bh")[1]; out["ASI_gene_vs_UV"] = R2.sort_values("p").head(8).round(4).to_dict("records")
# ---- 5 biomarker checks in the basket
mmr = an[an.CONSENSUS == "MMR"]; out["biomarker"] = {"consensus_MMR_n": int(len(mmr)), "consensus_MMR_MSI_TYPE": mmr.MSI_TYPE.value_counts().to_dict(),
    "consensus_MMR_MSI_SCORE_median": float(mmr.MSI_SCORE.median()), "consensus_MMR_N_SNV": mmr.N_SNV.tolist(), "basket_MSI_instable_n": int((b.MSI_TYPE == "Instable").sum()),
    "instable_cases": b[b.MSI_TYPE == "Instable"][["GROUP", "N_SNV", "MSI_SCORE", "CONSENSUS", "M1_SPA", "M2_MSKv2", "M3_cosine"]].to_dict("records"),
    "MSI_score_MMRcalled_vs_rest_MWU_p": float(stats.mannwhitneyu(mmr.MSI_SCORE.dropna(), an[an.CONSENSUS != "MMR"].MSI_SCORE.dropna()).pvalue),
    "hypermutator": b[b.HYPERMUTATOR][["GROUP", "PRIMARY_SITE", "N_SNV", "TMB_SCORE", "MSI_TYPE", "MSI_SCORE", "CONSENSUS", "M1_SPA", "M2_MSKv2", "M3_cosine"]].to_dict("records")}
json.dump(out, open(f"{RES}/stats_summary.json", "w"), indent=1, default=float); print(json.dumps(out, indent=1, default=float))
