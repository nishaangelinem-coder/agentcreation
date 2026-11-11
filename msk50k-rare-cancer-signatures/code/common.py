"""Shared helpers: paths, 96-channel definitions, signature families."""
import os, itertools, numpy as np, pandas as pd
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RAW = os.path.join(ROOT, "data", "raw", "study"); REF = os.path.join(ROOT, "data", "ref")
RES = os.path.join(ROOT, "results"); FIG = os.path.join(ROOT, "figures")
STUDY = os.path.join(RAW, "msk_impact_50k_2026")
SEED = 42
SUBS = ["C>A", "C>G", "C>T", "T>A", "T>C", "T>G"]
CHANNELS = [f"{a}[{s}]{b}" for s in SUBS for a in "ACGT" for b in "ACGT"]   # COSMIC order
CH_IDX = {c: i for i, c in enumerate(CHANNELS)}
TRI_OF_CH = [c[0] + c[2] + c[6] for c in CHANNELS]
COMP = str.maketrans("ACGT", "TGCA")
def revcomp(s): return s.translate(COMP)[::-1]
def channel(tri, ref, alt):
    """tri = 3-mer on + strand centred on mutated base. Returns pyrimidine-frame channel or None."""
    tri = tri.upper()
    if len(tri) != 3 or set(tri) - set("ACGT") or tri[1] != ref or alt not in "ACGT" or alt == ref: return None
    if ref in "AG": tri, ref, alt = revcomp(tri), ref.translate(COMP), alt.translate(COMP)
    return f"{tri[0]}[{ref}>{alt}]{tri[2]}"
# family-level synonym map (COSMIC v3 SBS -> process family)
FAMILY = {"SBS1": "Clock", "SBS5": "Clock", "SBS2": "APOBEC", "SBS13": "APOBEC",
          "SBS7a": "UV", "SBS7b": "UV", "SBS7c": "UV", "SBS7d": "UV", "SBS38": "UV", "SBS4": "Tobacco",
          "SBS3": "HRD", "SBS10a": "POLE", "SBS10b": "POLE", "SBS6": "MMR", "SBS15": "MMR", "SBS21": "MMR",
          "SBS26": "MMR", "SBS44": "MMR", "SBS11": "TMZ"}
REF_SIGS = list(FAMILY)   # curated reference set used for panel-scale fitting
# COSMIC v2 (MSK-supplied contributions) -> family
FAMILY_V2 = {1: "Clock", 5: "Clock", 2: "APOBEC", 13: "APOBEC", 7: "UV", 4: "Tobacco", 3: "HRD", 10: "POLE",
             6: "MMR", 15: "MMR", 20: "MMR", 21: "MMR", 26: "MMR", 11: "TMZ"}
FAMILIES = ["Clock", "APOBEC", "UV", "Tobacco", "HRD", "POLE", "MMR", "TMZ"]
def read_cbio(path):
    return pd.read_csv(path, sep="\t", comment="#", low_memory=False)
def cosine(a, b):
    a = np.asarray(a, float); b = np.asarray(b, float)
    d = np.linalg.norm(a) * np.linalg.norm(b)
    return float(a @ b / d) if d > 0 else 0.0
