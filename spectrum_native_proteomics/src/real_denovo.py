"""Real-data validation, step 3: de novo sequencing of real spectra (zero-shot and
after fine-tuning on calibrated simulation)."""
import sys, json, pickle, argparse, numpy as np, torch
from chem import ITOS, PROTON
from models import SpecNovo
from baselines import spectrum_graph_denovo
from metrics import summarize
from real_calib import load_massivekb, load_bsa

P = 80


def preprocess(r, max_peaks=P):
    mz, it = r['mz'], r['raw_it']
    keep = (mz > 50) & (mz < r['pmass']) & (np.abs(mz - r['pmz']) > 2.0)
    mz, it = mz[keep], it[keep]
    o = np.argsort(-it)[:max_peaks]; mz, it = mz[o], it[o]
    o = np.argsort(mz)
    return mz[o].astype(np.float32), np.sqrt(it[o] / it.max()).astype(np.float32)


def run_model(model, recs, beam=5):
    B = len(recs)
    mz = np.zeros((B, P), np.float32); it = np.zeros((B, P), np.float32); mask = np.ones((B, P), bool)
    for i, r in enumerate(recs):
        a, b = preprocess(r); mz[i, :len(a)] = a; it[i, :len(a)] = b; mask[i, :len(a)] = False
    pm = torch.tensor([r['pmass'] for r in recs], dtype=torch.float32)
    ch = torch.tensor([min(r['z'], 4) for r in recs])
    toks, conf = model.beam_search(torch.from_numpy(mz), torch.from_numpy(it), torch.from_numpy(mask), pm, ch, beam=beam)
    return [''.join(ITOS[x] for x in t if x >= 3) for t in toks.tolist()], conf.tolist()


def evaluate(ckpt, recs, tol_dp):
    out = {}
    trues = [r['seq'] for r in recs]
    if ckpt:
        m = SpecNovo(mass_aware=True, evidence=True); m.load_state_dict(torch.load(ckpt)); m.eval()
        pr, cf = run_model(m, recs); r, pm = summarize(pr, trues, cf); out['SpecNovo'] = r; out['_pred'] = pr
    sg = []
    for r in recs:
        a, b = preprocess(r); sg.append(spectrum_graph_denovo(a, b, r['pmass'], r['z'], tol=tol_dp))
    r2, _ = summarize([s[0] for s in sg], trues, [s[1] for s in sg]); out['SpectrumGraph-DP'] = r2
    return out


if __name__ == '__main__':
    ap = argparse.ArgumentParser(); ap.add_argument('--ckpt', default='results/specnovo_ev.pt'); ap.add_argument('--tag', default='zeroshot')
    a = ap.parse_args(); torch.set_num_threads(4)
    mk = [r for r in load_massivekb() if 6 <= len(r['seq']) <= 14 and r['z'] in (2, 3)]
    bsa, _ = load_bsa()
    bsa = [r for r in bsa if len(r['seq']) <= 14]
    res = dict(massivekb=evaluate(a.ckpt, mk, 0.025), bsa_cid=evaluate(a.ckpt, bsa, 0.5),
               n=dict(massivekb=len(mk), bsa=len(bsa), bsa_unique=len(set(r['seq'] for r in bsa))))
    for k in ('massivekb', 'bsa_cid'):
        print(k, {m: {kk: round(v, 3) for kk, v in d.items() if kk != 'curve'} for m, d in res[k].items() if not m.startswith('_')})
    if 'SpecNovo' in res['massivekb']:
        print(list(zip([r['seq'] for r in mk][:12], res['massivekb']['_pred'][:12])))
    json.dump(res, open(f'realdata/denovo_{a.tag}.json', 'w'), default=float)
