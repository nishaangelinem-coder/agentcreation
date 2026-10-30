import numpy as np
from chem import AA_MASS


def aa_match(pred, true, cum_tol=0.5, aa_tol=0.1):
    """DeepNovo/Casanovo-style mass-aligned amino-acid matching."""
    i = j = 0; cp = ct = 0.0; m = 0
    while i < len(pred) and j < len(true):
        mp, mt = AA_MASS[pred[i]], AA_MASS[true[j]]
        if abs(cp - ct) < cum_tol:
            if abs(mp - mt) < aa_tol: m += 1
            cp += mp; ct += mt; i += 1; j += 1
        elif cp < ct:
            cp += mp; i += 1
        else:
            ct += mt; j += 1
    return m


def pep_match(pred, true):
    if not pred or len(pred) != len(true): return False
    return all(abs(AA_MASS[a] - AA_MASS[b]) < 0.1 for a, b in zip(pred, true))


def summarize(preds, trues, conf=None):
    nm = sum(aa_match(p, t) for p, t in zip(preds, trues))
    npred = sum(len(p) for p in preds); ntrue = sum(len(t) for t in trues)
    pm = np.array([pep_match(p, t) for p, t in zip(preds, trues)])
    out = dict(aa_precision=nm / max(npred, 1), aa_recall=nm / ntrue, peptide_recall=pm.mean())
    if conf is not None:
        o = np.argsort(-np.asarray(conf))
        cov = np.arange(1, len(o) + 1) / len(o)
        prec = np.cumsum(pm[o]) / np.arange(1, len(o) + 1)
        out['pr_auc'] = float(np.trapezoid(prec, cov))
        out['curve'] = (cov[::max(1, len(o) // 200)].tolist(), prec[::max(1, len(o) // 200)].tolist())
    return out, pm
