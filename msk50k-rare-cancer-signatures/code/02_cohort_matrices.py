#!/usr/bin/env python
"""Phase 2-3a: cohort definition/QC and 96-channel SBS matrices (all 54,331 samples; GRCh37)."""
import json, re, numpy as np, pandas as pd, py2bit
from common import *
cs = read_cbio(f"{STUDY}/data_clinical_sample.txt"); cp = read_cbio(f"{STUDY}/data_clinical_patient.txt")
maf = pd.read_csv(f"{STUDY}/data_mutations.txt", sep="\t", low_memory=False,
                  usecols=["Hugo_Symbol", "NCBI_Build", "Chromosome", "Start_Position", "Variant_Classification", "Variant_Type",
                           "Reference_Allele", "Tumor_Seq_Allele2", "Tumor_Sample_Barcode", "HGVSp_Short", "ASCN.CLONAL"])
assert set(maf.NCBI_Build.dropna().unique()) == {"GRCh37"}
snv = maf[(maf.Variant_Type == "SNP") & maf.Reference_Allele.isin(list("ACGT")) & maf.Tumor_Seq_Allele2.isin(list("ACGT"))].copy()
tb = py2bit.open(f"{REF}/hg19.2bit"); chroms = tb.chroms()
def tri(c, p):
    c = "chr" + str(c).replace("chr", "").replace("MT", "M")
    return tb.sequence(c, p - 2, p + 1) if c in chroms else ""
snv["tri"] = [tri(c, p) for c, p in zip(snv.Chromosome, snv.Start_Position)]
snv["channel"] = [channel(t, r, a) for t, r, a in zip(snv.tri, snv.Reference_Allele, snv.Tumor_Seq_Allele2)]
qc = {"n_maf_rows": int(len(maf)), "n_snv": int(len(snv)), "n_ref_mismatch": int(snv.channel.isna().sum())}
snv = snv.dropna(subset=["channel"])
snv.to_parquet(f"{RES}/snv_annotated.parquet")
M = pd.crosstab(snv.Tumor_Sample_Barcode, snv.channel).reindex(index=cs.SAMPLE_ID, columns=CHANNELS, fill_value=0).fillna(0).astype(int)
M.to_parquet(f"{RES}/sbs96_all_samples.parquet")
cs["N_SNV"] = M.sum(1).values
# ---- cohort
GROUPS = {"ASI": ["ANGS", "BA"], "THYM": ["THYM", "THYC"], "ACC": ["ACC"]}
KW = {"ASI": "angiosarcoma", "THYM": "thymoma|thymic carcinoma", "ACC": "adrenocortical carcinoma"}
cs["CANCER_TYPE_DETAILED"] = cs.CANCER_TYPE_DETAILED.fillna(""); cs["GROUP"] = None
for g, codes in GROUPS.items():
    m = cs.ONCOTREE_CODE.isin(codes) & cs.CANCER_TYPE_DETAILED.str.contains(KW[g], case=False, regex=True) & ~cs.CANCER_TYPE_DETAILED.str.contains("glomangio", case=False)
    cs.loc[m, "GROUP"] = g
co = cs[cs.GROUP.notna()].copy(); flow = {"all_samples": int(len(cs)), "all_patients": int(cs.PATIENT_ID.nunique())}
flow["basket_samples"] = co.GROUP.value_counts().to_dict(); flow["basket_patients"] = co.groupby("GROUP").PATIENT_ID.nunique().to_dict()
co["matched_ok"] = co.SOMATIC_STATUS.eq("Matched"); flow["excluded_unmatched"] = int((~co.matched_ok).sum()); co = co[co.matched_ok]
co["prim_rank"] = co.SAMPLE_TYPE.map({"Primary": 0}).fillna(1)
co["tnum"] = co.SAMPLE_ID.str.extract(r"-T(\d+)-")[0].astype(int)
co = co.sort_values(["PATIENT_ID", "GROUP", "prim_rank", "tnum"]).drop_duplicates(["PATIENT_ID", "GROUP"])
flow["after_dedup"] = co.GROUP.value_counts().to_dict()
co["FIT_CLASS"] = np.select([co.N_SNV >= 4, co.N_SNV >= 2], ["FULL", "SPARSE"], "DROP")
co["HYPERMUTATOR"] = co.N_SNV > 50
flow["fit_class"] = pd.crosstab(co.GROUP, co.FIT_CLASS).to_dict("index")
flow["hypermutators"] = co[co.HYPERMUTATOR].GROUP.value_counts().to_dict()
flow["analysable"] = co[co.FIT_CLASS != "DROP"].GROUP.value_counts().to_dict()
flow["gate_ge50_per_group"] = {g: bool(flow["analysable"].get(g, 0) >= 50) for g in GROUPS}
co = co.merge(cp, on="PATIENT_ID", how="left"); co["AGE_AT_DX"] = pd.to_numeric(co.AGE_AT_DX, errors="coerce")
co.to_csv(f"{RES}/cohort.tsv", sep="\t", index=False)
json.dump({"qc": qc, "flow": flow}, open(f"{RES}/cohort_flow.json", "w"), indent=2)
print(json.dumps({"qc": qc, "flow": flow}, indent=2))
print(co.groupby("GROUP").agg(n=("SAMPLE_ID", "size"), snv_median=("N_SNV", "median"), snv_max=("N_SNV", "max"), tmb=("TMB_SCORE", "median"),
      age=("AGE_AT_DX", "median"), deaths=("OS_STATUS", lambda s: s.astype(str).str.contains("DECEASED").sum()), os_na=("OS_MONTHS", lambda s: s.isna().sum())))
print(pd.crosstab(co.GROUP, [co.SAMPLE_TYPE])); print(pd.crosstab(co.GROUP, co.GENE_PANEL)); print(pd.crosstab(co.GROUP, co.SEX)); print(pd.crosstab(co.GROUP, co.MSI_TYPE))
