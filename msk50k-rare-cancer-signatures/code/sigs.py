"""Reference signatures (COSMIC v3.5, GRCh37) in genome, exome and panel frames + fitting primitives."""
import os, numpy as np, pandas as pd, SigProfilerAssignment
from scipy.optimize import nnls
from common import *
_D = os.path.join(os.path.dirname(SigProfilerAssignment.__file__), "data", "Reference_Signatures", "GRCh37")
def _load(fn): return pd.read_csv(os.path.join(_D, fn), sep="\t", index_col=0).reindex(CHANNELS)
GEN_ALL = _load("COSMIC_v3.5_SBS_GRCh37.txt"); EXO_ALL = _load("COSMIC_v3.5_SBS_GRCh37_exome.txt")
_P = pd.read_csv(f"{RES}/trinucleotide_opportunity_proportions.tsv", sep="\t", index_col=0)
RATIO = (_P.IMPACT505 / _P.genome).reindex(TRI_OF_CH).values           # panel/genome opportunity per channel
PAN_ALL = GEN_ALL.mul(RATIO, axis=0); PAN_ALL = PAN_ALL / PAN_ALL.sum()
GEN, EXO, PAN = GEN_ALL[REF_SIGS], EXO_ALL[REF_SIGS], PAN_ALL[REF_SIGS]
FAM_OF = np.array([FAMILY[s] for s in REF_SIGS])
def fam_sum(expo):
    """expo: (n, len(REF_SIGS)) -> DataFrame (n, families)"""
    return pd.DataFrame({f: expo[:, FAM_OF == f].sum(1) for f in FAMILIES})
def fit_nnls(M, S):
    """M: (n,96) counts; S: (96,k). Returns exposures in mutation units (n,k)."""
    A = np.asarray(S, float); out = np.zeros((len(M), A.shape[1]))
    for i, m in enumerate(np.asarray(M, float)):
        if m.sum() > 0: x, _ = nnls(A, m); out[i] = x * (m.sum() / x.sum()) if x.sum() > 0 else 0
    return out
# family reference profiles for the cosine classifier (members + equal-weight composites)
def family_refs(S):
    refs = {}
    for f in FAMILIES:
        mem = [s for s in REF_SIGS if FAMILY[s] == f]
        refs[f] = [S[s].values for s in mem] + ([S[mem].mean(1).values] if len(mem) > 1 else [])
    return refs
def cosine_classify(M, S, thr):
    """MESiCA-style dominant-process call: max cosine of the normalised 96-profile to family references.
    thr: per-sample threshold array. Returns (family or 'Indeterminate', best cosine, table of family cosines)."""
    refs = family_refs(S); M = np.asarray(M, float)
    nm = np.linalg.norm(M, axis=1, keepdims=True); nm[nm == 0] = 1; U = M / nm
    C = pd.DataFrame({f: np.max([U @ (r / np.linalg.norm(r)) for r in rs], axis=0) for f, rs in refs.items()})
    best = C.idxmax(1).values; bc = C.max(1).values
    call = np.where(bc >= thr, best, "Indeterminate")
    return call, bc, C
