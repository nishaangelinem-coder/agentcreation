"""Bootstrap 95 % CIs for real-data de novo metrics (paired resampling of spectra)."""
import json, numpy as np, torch
from real_calib import load_massivekb, load_bsa
from real_denovo import evaluate, run_model, preprocess
from baselines import spectrum_graph_denovo
from metrics import aa_match, pep_match
from models import SpecNovo

mk = [r for r in load_massivekb() if 6 <= len(r['seq']) <= 14 and r['z'] in (2, 3)]
trues = [r['seq'] for r in mk]
preds = {}
for tag, ck in (('SpecNovo zero-shot', 'results/specnovo_ev.pt'), ('SpecNovo calibrated', 'results/specnovo_cal.pt')):
    m = SpecNovo(mass_aware=True, evidence=True); m.load_state_dict(torch.load(ck)); m.eval()
    preds[tag] = run_model(m, mk)
sg = [spectrum_graph_denovo(*preprocess(r), r['pmass'], r['z'], tol=0.025) for r in mk]
preds['SpectrumGraph-DP'] = ([s[0] for s in sg], [s[1] for s in sg])
rng = np.random.default_rng(0); n = len(mk); out = {}
per = {k: (np.array([pep_match(p, t) for p, t in zip(v[0], trues)]), np.array([aa_match(p, t) for p, t in zip(v[0], trues)]),
           np.array([len(p) for p in v[0]]), np.array(v[1])) for k, v in preds.items()}
L = np.array([len(t) for t in trues])
def prauc(pm, cf):
    o = np.argsort(-cf); prec = np.cumsum(pm[o]) / np.arange(1, len(o) + 1); return float(np.trapezoid(prec, np.arange(1, len(o) + 1) / len(o)))
for k, (pm, am, lp, cf) in per.items():
    bs = []
    for _ in range(2000):
        i = rng.integers(0, n, n); bs.append((pm[i].mean(), am[i].sum() / max(lp[i].sum(), 1), prauc(pm[i], cf[i])))
    bs = np.array(bs)
    out[k] = dict(peptide_recall=[float(pm.mean()), *np.percentile(bs[:, 0], [2.5, 97.5]).tolist()],
                  aa_precision=[float(am.sum() / lp.sum()), *np.percentile(bs[:, 1], [2.5, 97.5]).tolist()],
                  pr_auc=[prauc(pm, cf), *np.percentile(bs[:, 2], [2.5, 97.5]).tolist()])
d = per['SpecNovo calibrated'][0].astype(int) - per['SpecNovo zero-shot'][0].astype(int)
bd = [d[rng.integers(0, n, n)].mean() for _ in range(2000)]
out['delta_recall_calibrated_minus_zeroshot'] = [float(d.mean()), *np.percentile(bd, [2.5, 97.5]).tolist()]
out['n'] = n
json.dump(out, open('realdata/denovo_bootstrap.json', 'w'), indent=1); print(json.dumps(out, indent=1))
