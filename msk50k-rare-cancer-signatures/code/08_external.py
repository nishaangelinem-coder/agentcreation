#!/usr/bin/env python
"""Phase 5.6 external validation in whole-exome cohorts (Angiosarcoma Project, TCGA-ACC, TCGA-THYM) and
in-silico panel down-sampling (restrict WES calls to the IMPACT505 territory and re-call)."""
import bisect, json, warnings, numpy as np, pandas as pd, py2bit
from scipy import stats
from lifelines import CoxPHFitter
from common import *; from sigs import *
warnings.filterwarnings("ignore"); rng = np.random.default_rng(SEED); tb = py2bit.open(f"{REF}/hg19.2bit"); chroms = tb.chroms(); out = {}
bed = pd.read_csv(f"{REF}/IMPACT505.cds_flank10.hg19.bed", sep="\t", names=["c", "s", "e"]); idx = {c: list(zip(d.s, d.e)) for c, d in bed.groupby("c")}
def inside(c, p):
    L = idx.get(c, []); i = bisect.bisect_right(L, (p, 10**12)) - 1
    return i >= 0 and L[i][0] < p <= L[i][1]
def load(study):
    m = pd.read_csv(f"{RAW}/{study}/data_mutations.txt", sep="\t", comment="#", low_memory=False, usecols=["Hugo_Symbol", "Chromosome", "Start_Position", "Variant_Type", "Reference_Allele", "Tumor_Seq_Allele2", "Tumor_Sample_Barcode"])
    m = m[(m.Variant_Type == "SNP") & m.Reference_Allele.isin(list("ACGT")) & m.Tumor_Seq_Allele2.isin(list("ACGT"))].copy()
    m["chr"] = "chr" + m.Chromosome.astype(str).str.replace("chr", ""); m = m[m.chr.isin(chroms)]
    m["channel"] = [channel(tb.sequence(c, int(p) - 2, int(p) + 1), r, a) for c, p, r, a in zip(m.chr, m.Start_Position, m.Reference_Allele, m.Tumor_Seq_Allele2)]
    m = m.dropna(subset=["channel"]); m["in_panel"] = [inside(c, int(p)) for c, p in zip(m.chr, m.Start_Position)]
    W = pd.crosstab(m.Tumor_Sample_Barcode, m.channel).reindex(columns=CHANNELS, fill_value=0)
    P = pd.crosstab(m[m.in_panel].Tumor_Sample_Barcode, m[m.in_panel].channel).reindex(index=W.index, columns=CHANNELS, fill_value=0).fillna(0).astype(int)
    return m, W, P
def famfrac(Mx, S):
    F = fam_sum(fit_nnls(Mx.values, S)); F.index = Mx.index; return F.div(F.sum(1).replace(0, np.nan), axis=0)
def pooled(tot, ids_M, S):
    f = fam_sum(fit_nnls(tot[None, :], S)).iloc[0]; f = f / f.sum()
    bs = np.array([ids_M.values[rng.integers(0, len(ids_M), len(ids_M))].sum(0) for _ in range(500)]); bf = fam_sum(fit_nnls(bs, S)); bf = bf.div(bf.sum(1), axis=0)
    return {k: [round(float(f[k]), 3), round(float(bf[k].quantile(.025)), 3), round(float(bf[k].quantile(.975)), 3)] for k in FAMILIES if bf[k].quantile(.975) > 0.03}
def downsample_eval(W, P, label):
    """WES truth (dominant family, exome frame, >=50 SNVs) vs panel-territory re-call of the same tumour."""
    wn, pn = W.sum(1), P.sum(1); truth = famfrac(W, EXO).idxmax(1)[wn >= 50]
    ok = [i for i in truth.index if pn[i] >= 2]; r = {"n_WES_ge50": int(len(truth)), "median_WES_SNV": float(wn.median()), "median_panel_SNV": float(pn.median()),
                                                     "n_panel_ge2": int(len(ok)), "n_panel_ge4": int((pn[truth.index] >= 4).sum())}
    if ok:
        call, bc, _ = cosine_classify(P.loc[ok].values, PAN, np.where(pn[ok] >= 4, 0.60, 0.80)); call = pd.Series(call, index=ok)
        nn = famfrac(P.loc[ok], PAN).idxmax(1); called = call != "Indeterminate"
        r.update({"cosine_call_rate": float(called.mean()), "cosine_acc_when_called": float((call[called] == truth[ok][called]).mean()) if called.any() else None,
                  "nnls_dominant_acc": float((nn == truth[ok]).mean()), "truth_distribution": truth[ok].value_counts().to_dict()})
        for fam in ["UV"]:
            t = truth[ok] == fam
            if t.sum() >= 5: r[f"{fam}_sens_cosine"] = float((call[t] == fam).mean()); r[f"{fam}_spec_cosine"] = float((call[~t] != fam).mean()); r[f"{fam}_n_truth"] = int(t.sum())
    return r
# ---------- Angiosarcoma Project (provisional 2025 release; WES)
for study in ["angs_painter_2025", "angs_project_painter_2018"]:
    m, W, P = load(study); cs = read_cbio(f"{RAW}/{study}/data_clinical_sample.txt"); cp = read_cbio(f"{RAW}/{study}/data_clinical_patient.txt")
    d = cs.merge(cp, on="PATIENT_ID", how="left").set_index("SAMPLE_ID").reindex(W.index); r = {"n_samples_with_SNV": int(len(W)), "n_patients": int(d.PATIENT_ID.nunique())}
    loc = d.BX_LOCATION.astype(str).str.lower() if "BX_LOCATION" in d else pd.Series("", index=d.index)
    hn = loc.str.contains("scalp|face|head|neck|forehead|nose|cheek|ear|eyelid|temple")
    r["BX_LOCATION_top"] = loc.value_counts().head(12).to_dict()
    F = famfrac(W, EXO); dom = F.idxmax(1); n = W.sum(1); keep = n >= 20
    one = d[keep].assign(n=n[keep]).sort_values("n", ascending=False).drop_duplicates("PATIENT_ID").index       # one sample per patient
    r["n_patients_ge20SNV"] = int(len(one)); t = pd.crosstab(hn[one], dom[one] == "UV"); r["UV_dominant_by_headneck"] = {str(k): v for k, v in t.to_dict("index").items()}
    if t.shape == (2, 2): o, p = stats.fisher_exact(t); r["fisher_OR"] = float(o); r["fisher_p"] = float(p)
    r["UV_frac_median_headneck_vs_other"] = [float(F.UV[one][hn[one]].median()), float(F.UV[one][~hn[one]].median())]; r["MWU_p"] = float(stats.mannwhitneyu(F.UV[one][hn[one]], F.UV[one][~hn[one]]).pvalue)
    r["pooled_headneck"] = pooled(W.loc[one][hn[one]].sum(0).values, W.loc[one][hn[one]], EXO); r["pooled_other"] = pooled(W.loc[one][~hn[one]].sum(0).values, W.loc[one][~hn[one]], EXO)
    r["panel_downsampling"] = downsample_eval(W.loc[one], P.loc[one], study)
    # pooled panel-territory spectrum of the same tumours (what a panel would see)
    r["pooled_headneck_panel_territory"] = pooled(P.loc[one][hn[one]].sum(0).values, P.loc[one][hn[one]], PAN)
    pd.DataFrame({"n_wes": n, "n_panel": P.sum(1), "UV_frac": F.UV, "dominant": dom, "head_neck": hn}).to_csv(f"{RES}/external_{study}.tsv", sep="\t")
    out[study] = r
# ---------- TCGA ACC / THYM
for study, g in [("acc_tcga_pan_can_atlas_2018", "ACC"), ("thym_tcga_pan_can_atlas_2018", "THYM")]:
    m, W, P = load(study); cs = read_cbio(f"{RAW}/{study}/data_clinical_sample.txt"); cp = read_cbio(f"{RAW}/{study}/data_clinical_patient.txt")
    d = cs.merge(cp, on="PATIENT_ID", how="left").set_index("SAMPLE_ID").reindex(W.index); n = W.sum(1); F = famfrac(W, EXO)
    r = {"n": int(len(W)), "median_SNV": float(n.median()), "dominant_family_ge20": F.idxmax(1)[n >= 20].value_counts().to_dict(), "n_ge20": int((n >= 20).sum()),
         "pooled_WES": pooled(W.sum(0).values, W, EXO), "pooled_panel_territory": pooled(P.sum(0).values, P, PAN), "panel_territory_SNV_total": int(P.values.sum()),
         "panel_downsampling": downsample_eval(W, P, study)}
    d["E"] = d.OS_STATUS.astype(str).str.startswith("1").astype(int); d["T"] = pd.to_numeric(d.OS_MONTHS, errors="coerce").clip(lower=0.01); d["log10TMB"] = np.log10(pd.to_numeric(d.TMB_NONSYNONYMOUS, errors="coerce") + 1)
    d["AGE"] = pd.to_numeric(d.AGE, errors="coerce"); x = d[["T", "E", "AGE", "log10TMB"]].dropna()
    cm = CoxPHFitter().fit(x, "T", "E"); r["cox_AGE_log10TMB"] = {"n": int(len(x)), "events": int(x.E.sum()), "hr": cm.summary[["exp(coef)", "exp(coef) lower 95%", "exp(coef) upper 95%", "p"]].round(4).to_dict("index")}
    if study.startswith("thym"): r["histology"] = d.CANCER_TYPE_DETAILED.value_counts().to_dict()
    out[study] = r
json.dump(out, open(f"{RES}/external_validation.json", "w"), indent=1); print(json.dumps(out, indent=1))
