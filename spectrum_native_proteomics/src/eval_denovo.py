import sys, json, pickle, time, argparse
import numpy as np, torch
from chem import TOK, ITOS, EOS, PAD, fragment_ladders, PROTON, H2O, NH3
from models import SpecNovo
from baselines import DBSearch, fdr_filter, spectrum_graph_denovo
from metrics import summarize, pep_match
from simulator import build_ms2_set

ap = argparse.ArgumentParser()
ap.add_argument('--ckpt', default='results/specnovo_main.pt')
ap.add_argument('--tag', default='main')
ap.add_argument('--full', type=int, default=1)
ap.add_argument('--threads', type=int, default=4)
ap.add_argument('--mass_aware', type=int, default=1)
ap.add_argument('--evidence', type=int, default=1)
a = ap.parse_args()
torch.set_num_threads(a.threads)
D = pickle.load(open('data/ms2_corpus.pkl', 'rb'))
P = 80
model = SpecNovo(mass_aware=bool(a.mass_aware), evidence=bool(a.evidence)); model.load_state_dict(torch.load(a.ckpt)); model.eval()


def batch(items):
    B = len(items)
    mz = np.zeros((B, P), np.float32); it = np.zeros((B, P), np.float32); mask = np.ones((B, P), bool)
    for i, (m, s, *_ ) in enumerate(items):
        mz[i, :len(m)] = m; it[i, :len(m)] = s; mask[i, :len(m)] = False
    pm = np.array([x[2] for x in items], np.float32); ch = np.array([x[3] for x in items])
    return torch.from_numpy(mz), torch.from_numpy(it), torch.from_numpy(mask), torch.from_numpy(pm), torch.from_numpy(ch)


def run_model(items, bs=256, beam=5):
    preds, confs = [], []
    for i in range(0, len(items), bs):
        toks, conf = model.beam_search(*batch(items[i:i + bs]), beam=beam)
        for t in toks.tolist():
            preds.append(''.join(ITOS[x] for x in t if x >= 3))
        confs += conf.tolist()
    return preds, confs


tests = D['tests'] if a.full else {k: v for k, v in D['tests'].items() if k in ('inDB',)}
out = {}
db = DBSearch(D['db']) if a.full else None
for name, items in tests.items():
    t = time.time()
    trues = [x[4] for x in items]
    pr, cf = run_model(items); tm = time.time() - t
    r, pm = summarize(pr, trues, cf); r['secs_per_1k'] = tm / len(items) * 1000
    out.setdefault(name, {})['SpecNovo'] = r
    if a.full:
        sg = [spectrum_graph_denovo(*x[:4]) for x in items]
        r2, _ = summarize([s[0] for s in sg], trues, [s[1] for s in sg])
        out[name]['SpectrumGraph-DP'] = r2
        R = [db.search(*x[:4]) for x in items]
        acc = fdr_filter([q[1] for q in R], [q[2] for q in R])
        corr = np.array([q[0] is not None and pep_match(q[0], tt) for q, tt in zip(R, trues)])
        out[name]['DB-search'] = dict(peptide_recall=float((acc & corr).mean()), false_ids=float((acc & ~corr).mean()),
                                      accepted=float(acc.mean()))
        # SpecNovo at the score threshold giving 99% precision on inDB is computed later
        out[name]['_specnovo_conf'] = cf; out[name]['_specnovo_correct'] = pm.tolist()
        out[name]['_db_score'] = [q[1] for q in R]
    print(name, {k: {kk: (round(vv, 4) if isinstance(vv, float) else '') for kk, vv in v.items() if not kk.startswith('curve')}
                 for k, v in out[name].items() if not k.startswith('_')}, flush=True)
json.dump(out, open(f'results/denovo_eval_{a.tag}.json', 'w'))

if a.full:
    # -------- representation analysis: replicate spectra of 40 peptides
    rng = np.random.default_rng(7)
    peps = [D['tests']['inDB'][i][4] for i in rng.choice(3000, 40, replace=False)]
    reps = []
    for p in peps:
        reps += build_ms2_set([p] * 12, rng, max_peaks=80, p_obs=0.75, noise_peaks=60)
    with torch.no_grad():
        mz, it, mask, pm, ch = batch(reps)
        mem = model.encode(mz, it, mask, pm)
        emb = (mem * (~mask).unsqueeze(-1)).sum(1) / (~mask).sum(1, keepdim=True)
    # binned spectrum baseline (0.5 Da bins)
    binv = np.zeros((len(reps), 4000), np.float32)
    for i, (m, s, *_ ) in enumerate(reps):
        np.add.at(binv[i], np.clip((m / 0.5).astype(int), 0, 3999), s)
    labels = np.repeat(np.arange(40), 12)
    # -------- peak saliency on 300 test spectra
    items = D['tests']['inDB'][:300]
    sal_frag, int_frag, examples = [], [], []
    for x in items:
        mzt, itt, mk, pmt, cht = batch([x])
        itt = itt.clone().requires_grad_(True)
        seq = x[4]
        tok = torch.tensor([[TOK[c] for c in seq]])
        lg = model(mzt, itt, mk, pmt, cht, tok[:, :-1] if False else tok)
        tgt = torch.tensor([TOK[c] for c in seq] + [EOS])
        lp = torch.log_softmax(lg[0], -1)[torch.arange(len(tgt)), tgt].sum()
        lp.backward()
        g = (itt.grad[0] * itt[0]).abs().detach().numpy()[:len(x[0])]
        b, y = fragment_ladders(seq)
        th = np.r_[b, y, b - H2O, y - H2O, b - NH3, y - NH3, (b + PROTON) / 2, (y + PROTON) / 2]
        isf = np.array([np.min(np.abs(th - m)) < 0.02 for m in x[0]])
        sal_frag.append(g[isf].sum() / (g.sum() + 1e-12)); int_frag.append(x[1][isf].sum() / x[1].sum())
        if len(examples) < 3: examples.append(dict(mz=x[0], it=x[1], sal=g, isf=isf, seq=seq, z=x[3]))
    pickle.dump(dict(emb=emb.numpy(), binv=binv, labels=labels, peps=peps, sal_frag=sal_frag,
                     int_frag=int_frag, examples=examples), open(f'results/denovo_analysis_{a.tag}.pkl', 'wb'))
print('done')
