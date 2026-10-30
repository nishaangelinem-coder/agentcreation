import argparse, json, time, pickle, warnings
import numpy as np, torch
warnings.filterwarnings('ignore')
from pheno2 import *

ap = argparse.ArgumentParser()
ap.add_argument('--mode', default='main')
ap.add_argument('--seed', type=int, default=0)
ap.add_argument('--fracs', default='0.5')
ap.add_argument('--threads', type=int, default=1)
a = ap.parse_args()
torch.set_num_threads(a.threads)
allres = []
for frac in [float(f) for f in a.fracs.split(',')]:
    torch.manual_seed(a.seed)
    rng = np.random.default_rng(1000 + a.seed)
    sim = CohortSimulator(rng, n_feat=300, frac_lib=0.7, n_bio=20, frac_unid=frac, log2fc=0.35)
    t0 = time.time()
    yu, bu, abu = sim.sample_ab(800)
    ypo, bpo, abpo = sim.sample_ab(384)
    yte, bte, abte = sim.sample_ab(500)
    # consensus feature map from unlabeled runs (oracle coordinates + calibration noise)
    ref_mz = sim.mz * (1 + rng.normal(0, 2, sim.n_feat) * 1e-6); ref_rt = sim.rt + rng.normal(0, .1, sim.n_feat)
    Vu, _, _, _ = build_matrices(sim, abu, bu, rng, ref_mz, ref_rt)
    Vpo, Ipo, shpo, acc = build_matrices(sim, abpo, bpo, rng, ref_mz, ref_rt)
    Vte, Ite, shte, _ = build_matrices(sim, abte, bte, rng, ref_mz, ref_rt)
    norm = Norm(Vu); vu, mu = norm(Vu); vpo, mpo = norm(Vpo); vte, mte = norm(Vte)
    pre = make_model(sim, ref_mz, ref_rt)
    hist = pretrain(pre, vu, mu, rng, steps=250)
    state = {k: v.clone() for k, v in pre.state_dict().items()}
    print('frac', frac, 'prep+ssl', time.time() - t0, 'ssl', hist[0], hist[-1], flush=True)
    res = dict(mode=a.mode, seed=a.seed, frac=frac, ssl_hist=hist, align_acc=acc,
               miss_idfree=float(np.isnan(Vpo).mean()), miss_idfirst=float(np.isnan(Ipo).mean()), rows=[])
    methods = ['ID-first-LR', 'ID-first+MBR-LR', 'ID-free-LR', 'FeatureFormer', 'FeatureFormer+SSL'] if a.mode == 'main' \
        else ['ID-first-LR', 'ID-first+MBR-LR', 'ID-free-LR', 'FeatureFormer+SSL']
    lib = sim.in_lib
    keep = {}

    def run(tr, te_idx, tag, n):
        for meth in methods:
            tt = time.time()
            if meth == 'ID-first-LR': p = lr_pipeline(Ipo[tr], ypo[tr], te_idx[0])
            elif meth == 'ID-first+MBR-LR': p = lr_pipeline(Vpo[tr][:, lib], ypo[tr], te_idx[1][:, lib])
            elif meth == 'ID-free-LR': p = lr_pipeline(Vpo[tr], ypo[tr], te_idx[1])
            else:
                m = make_model(sim, ref_mz, ref_rt)
                if meth.endswith('SSL'):
                    m.load_state_dict(state); m = finetune(m, vpo[tr], mpo[tr], ypo[tr], rng, steps=150, lr=1e-3, body_lr=2e-4)
                else:
                    m = finetune(m, vpo[tr], mpo[tr], ypo[tr], rng, steps=150, lr=5e-4)
                p = predict(m, te_idx[2], te_idx[3]); keep['model_' + meth] = m
            r = evaluate(p, te_idx[4]); r.update(method=meth, n=int(tr.sum()) if tr.dtype == bool else len(tr),
                                                  setting=tag, secs=time.time() - tt)
            res['rows'].append(r); keep[meth] = p; print(r, flush=True)

    ns = [24, 48, 96, 192, 384] if a.mode == 'main' else [192]
    for n in ns:
        tr = np.r_[np.where(ypo == 0)[0][:n // 2], np.where(ypo == 1)[0][:n // 2]]
        run(tr, (Ite, Vte, vte, mte, yte), 'iid', n)
    if a.mode == 'main':
        # leave-one-site-out: train sites 0-2, test unseen site 3
        tr = np.where(bpo != 3)[0]; te = bte == 3
        run(tr, (Ite[te], Vte[te], vte[te], mte[te], yte[te]), 'LOSO', len(tr))
        # interpretability: gradient x input attribution on test runs (n=384 SSL model)
        tr = np.arange(384)
        m = make_model(sim, ref_mz, ref_rt); m.load_state_dict(state)
        m = finetune(m, vpo[tr], mpo[tr], ypo[tr], rng, steps=150, lr=1e-3, body_lr=2e-4); m.eval()
        v = vte.clone().requires_grad_(True); lg = m(v, mte)
        (lg[:, 1] - lg[:, 0]).sum().backward()
        attr = (v.grad * v).abs().mean(0).detach().numpy()
        Xi, Xt_ = prep_stats(Ipo[tr], Ite); sc = StandardScaler().fit(Xi)
        from sklearn.linear_model import LogisticRegressionCV as L
        lr = L(Cs=10, cv=5, max_iter=5000, scoring='roc_auc').fit(sc.transform(Xi), ypo[tr])
        # univariate statistics on the ID-free matrix (Welch t-test)
        from scipy.stats import ttest_ind
        Xf, _ = prep_stats(Vpo, Vte) if False else (None, None)
        Vn = Vpo - np.nanmedian(Vpo, 1, keepdims=True)
        tst = ttest_ind(Vn[ypo == 1], Vn[ypo == 0], nan_policy='omit', equal_var=False)
        if a.seed == 0:
            pickle.dump(dict(attr=attr, bio=sim.bio, in_lib=sim.in_lib, effect=sim.effect, mz=sim.mz, rt=sim.rt,
                             base=sim.base, lr_coef=lr.coef_[0], ttest_p=np.asarray(tst.pvalue), ttest_t=np.asarray(tst.statistic),
                             probs={k: v for k, v in keep.items() if not k.startswith('model')}, yte=yte, bte=bte,
                             shifts=shte, batch_rt=sim.batch_rt, Vte=Vte[:200], Ite=Ite[:200]),
                        open('results/pheno2_seed0_extra.pkl', 'wb'))
        res['attr_auc'] = float(roc_auc_score(np.isin(np.arange(sim.n_feat), sim.bio), attr))
        print('attr AUROC', res['attr_auc'])
    allres.append(res)
    json.dump(allres, open(f'results/pheno2_{a.mode}_s{a.seed}.json', 'w'))
print('done')
