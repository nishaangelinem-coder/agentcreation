"""Classical baselines: (i) target-decoy database search with a
hyperscore-like PSM score, (ii) spectrum-graph de novo sequencing (DP)."""
import numpy as np
from math import lgamma, log
from chem import AA_MASS, PROTON, H2O, peptide_mass, fragment_ladders

FRAG_TOL = 0.02


class DBSearch:
    def __init__(self, targets, prec_ppm=10.0):
        dec = [p[:-1][::-1] + p[-1] for p in targets]
        self.peps = list(targets) + dec
        self.is_decoy = np.r_[np.zeros(len(targets), bool), np.ones(len(dec), bool)]
        self.mass = np.array([peptide_mass(p) for p in self.peps])
        o = np.argsort(self.mass)
        self.peps = [self.peps[i] for i in o]; self.mass = self.mass[o]; self.is_decoy = self.is_decoy[o]
        self.ppm = prec_ppm

    @staticmethod
    def hyperscore(seq, mz, it, z):
        b, y = fragment_ladders(seq)
        sc, nb, ny = 0.0, 0, 0
        for ladder, kind in ((b, 'b'), (y, 'y')):
            th = ladder if z == 2 else np.r_[ladder, (ladder + PROTON) / 2]
            j = np.searchsorted(mz, th)
            j0 = np.clip(j - 1, 0, len(mz) - 1); j1 = np.clip(j, 0, len(mz) - 1)
            d0 = np.abs(mz[j0] - th); d1 = np.abs(mz[j1] - th)
            jj = np.where(d0 < d1, j0, j1); dd = np.minimum(d0, d1)
            hit = dd < FRAG_TOL
            sc += it[jj[hit]].sum()
            if kind == 'b': nb = hit.sum()
            else: ny = hit.sum()
        if sc <= 0: return 0.0
        return log(sc + 1e-9) + lgamma(nb + 1) + lgamma(ny + 1)

    def search(self, mz, it, pm, z):
        tol = pm * self.ppm * 1e-6
        lo, hi = np.searchsorted(self.mass, [pm - tol, pm + tol])
        best, bs, bd = None, -1.0, False
        for k in range(lo, hi):
            s = self.hyperscore(self.peps[k], mz, it, z)
            if s > bs: best, bs, bd = self.peps[k], s, self.is_decoy[k]
        return best, bs, bd


def fdr_filter(scores, decoy, q=0.01):
    """Target-decoy competition: return boolean accepted mask at FDR q."""
    scores = np.asarray(scores); decoy = np.asarray(decoy)
    o = np.argsort(-scores)
    d = np.cumsum(decoy[o]); t = np.cumsum(~decoy[o])
    fdr = (d + 1) / np.maximum(t, 1)
    qv = np.minimum.accumulate(fdr[::-1])[::-1]
    acc = np.zeros(len(scores), bool)
    acc[o] = (qv <= q) & (~decoy[o])
    return acc


# --------------------------------------------------------------- spectrum graph
_AA = [a for a in AA_MASS]
_single = [(AA_MASS[a], a) for a in _AA]
_pairs = [(AA_MASS[a] + AA_MASS[b], a + b) for i, a in enumerate(_AA) for b in _AA[i:]]
_EDGE = sorted(_single + _pairs)
_EM = np.array([e[0] for e in _EDGE]); _ES = [e[1] for e in _EDGE]


def spectrum_graph_denovo(mz, it, pm, z, tol=0.025):
    Mres = pm - H2O
    nodes = [(0.0, 0.0), (Mres, 0.0)]
    for m, s in zip(mz, it):
        nodes.append((m - PROTON, s))                 # b-ion interpretation
        nodes.append((Mres - (m - H2O - PROTON), s))  # y-ion complement
    nodes = [n for n in nodes if -0.01 <= n[0] <= Mres + 0.01]
    nodes.sort()
    # merge close nodes
    merged = []
    for m, s in nodes:
        if merged and m - merged[-1][0] < tol:
            pm_, ps = merged[-1]; merged[-1] = ((pm_ * ps + m * s) / (ps + s + 1e-12) if ps + s > 0 else pm_, ps + s)
        else:
            merged.append((m, s))
    M = np.array([n[0] for n in merged]); S = np.array([n[1] for n in merged])
    n = len(M)
    best = np.full(n, -np.inf); prev = np.full(n, -1); lab = [''] * n
    best[0] = 0.0
    for j in range(1, n):
        d = M[j] - M[:j]
        k = np.searchsorted(_EM, d)
        k0 = np.clip(k - 1, 0, len(_EM) - 1); k1 = np.clip(k, 0, len(_EM) - 1)
        e0 = np.abs(_EM[k0] - d); e1 = np.abs(_EM[k1] - d)
        ke = np.where(e0 < e1, k0, k1); err = np.minimum(e0, e1)
        ok = (err < tol) & np.isfinite(best[:j])
        if not ok.any(): continue
        pen = np.array([0.0 if len(_ES[q]) == 1 else 0.25 for q in ke])
        cand = np.where(ok, best[:j] + S[j] - pen, -np.inf)
        i = int(cand.argmax()); best[j] = cand[i]; prev[j] = i; lab[j] = _ES[ke[i]]
    end = int(np.argmin(np.abs(M - Mres)))
    if not np.isfinite(best[end]) or abs(M[end] - Mres) > 0.05:
        return '', 0.0
    seq = []; j = end
    while j > 0:
        seq.append(lab[j]); j = prev[j]
    seq = ''.join(seq[::-1])
    return seq, best[end] / max(len(seq), 1)
