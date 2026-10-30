import argparse, json, time, pickle
import numpy as np, torch
from pheno import *

ap = argparse.ArgumentParser()
ap.add_argument('--mode', default='main')
ap.add_argument('--seed', type=int, default=0)
ap.add_argument('--frac', type=float, default=0.5)
ap.add_argument('--threads', type=int, default=1)
a = ap.parse_args()
torch.set_num_threads(a.threads); torch.manual_seed(a.seed)
rng = np.random.default_rng(100 + a.seed)
sim = CohortSimulator(rng, frac_unid=a.frac)
t0 = time.time()
Xu, _, _, _ = sim.sample(1200)                    # unlabeled corpus (all sites)
Xpool, ypool, bpool, _ = sim.sample(384)
Xte, yte, bte, _ = sim.sample(400)
mu, sd = Xu.mean(), Xu.std()
nz = lambda X: ((X - mu) / sd).astype(np.float32)
pre, mae_hist = pretrain_mae(nz(Xu), rng, steps=500 if a.mode == 'main' else 350)
print('pretrained', time.time() - t0, flush=True)
res = dict(mode=a.mode, seed=a.seed, frac=a.frac, mae_hist=mae_hist, rows=[])
ns = [24, 48, 96, 192, 384] if a.mode == 'main' else [192]
methods = ['ID-first', 'Pixel-LR', 'CNN', 'MS1Former-scratch', 'MS1Former-MAE'] if a.mode == 'main' \
    else ['ID-first', 'Pixel-LR', 'MS1Former-MAE']

def run_all(Xtr, ytr, Xev, yev, tag, keep=None):
    for meth in methods:
        tt = time.time()
        if meth == 'ID-first': p = id_first(sim, None, ytr, None, Xtr, Xev)
        elif meth == 'Pixel-LR': p = pixel_lr(nz(Xtr), ytr, nz(Xev))
        elif meth == 'CNN': m = train_clf(SmallCNN(), nz(Xtr), ytr, rng, lr=1e-3); p = predict(m, nz(Xev))
        elif meth == 'MS1Former-scratch': m = train_clf(MS1Former(), nz(Xtr), ytr, rng, lr=5e-4); p = predict(m, nz(Xev))
        else: m = finetune(pre, nz(Xtr), ytr, rng); p = predict(m, nz(Xev))
        r = evaluate(p, yev); r.update(method=meth, n=len(Xtr), setting=tag, secs=time.time() - tt)
        res['rows'].append(r); print(r, flush=True)
        if keep is not None:
            keep[meth] = p
            if meth == 'MS1Former-MAE': keep['model'] = m

for n in ns:
    idx = np.r_[np.where(ypool == 0)[0][:n // 2], np.where(ypool == 1)[0][:n // 2]]
    keep = {} if (a.mode == 'main' and n == 384) else None
    run_all(Xpool[idx], ypool[idx], Xte, yte, 'iid', keep)
if a.mode == 'main':
    # leave-one-batch-out: train on sites 0-2, test on unseen site 3
    tr = bpool != 3; te = bte == 3
    run_all(Xpool[tr], ypool[tr], Xte[te], yte[te], 'LOBO')
    if a.seed == 0:
        m = keep['model']; m.eval()
        X = torch.from_numpy(nz(Xte)).requires_grad_(True)
        logit = m(X)
        (logit[:, 1] - logit[:, 0]).sum().backward()
        sal = (X.grad * X).abs().detach().numpy().mean(0)
        pickle.dump(dict(sal=sal, Xte=Xte[:8], yte=yte, probs={k: v for k, v in keep.items() if k != 'model'},
                         sim=dict(mz=sim.mz, rt=sim.rt, z=sim.z, bio=sim.bio, in_lib=sim.in_lib, effect=sim.effect,
                                  batch_rt=sim.batch_rt), Xmean_case=Xte[yte == 1].mean(0), Xmean_ctrl=Xte[yte == 0].mean(0),
                         lobo_b=bte), open('results/pheno_seed0_extra.pkl', 'wb'))
json.dump(res, open(f'results/pheno_{a.mode}_f{a.frac}_s{a.seed}.json', 'w'))
print('total', time.time() - t0)
