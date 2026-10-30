"""Fig. 6: spectrum-to-phenotype benchmark."""
import json, glob, pickle, numpy as np
from figstyle import *
from sklearn.metrics import roc_curve, roc_auc_score

rows = []
for f in sorted(glob.glob('results/pheno2_*_s*.json')):
    for r in json.load(open(f)):
        for q in r['rows']: rows.append(dict(q, frac=r['frac'], seed=r['seed'], mode=r['mode']))
M = ['ID-first-LR', 'ID-first+MBR-LR', 'ID-free-LR', 'FeatureFormer', 'FeatureFormer+SSL']
MC = {'ID-first-LR': C[3], 'ID-first+MBR-LR': C[4], 'ID-free-LR': C[2], 'FeatureFormer': C[6], 'FeatureFormer+SSL': C[0]}
def agg(mode, setting, meth, n=None, frac=0.5, key='auroc'):
    v = [r[key] for r in rows if r['mode'] == mode and r['setting'] == setting and r['method'] == meth
         and (n is None or r['n'] == n or (n == 384 and r['n'] >= 370)) and abs(r['frac'] - frac) < 1e-9]
    return (np.mean(v), np.std(v), len(v)) if v else (np.nan, np.nan, 0)

table = {}
fig = plt.figure(figsize=(7.2, 6.2)); gs = fig.add_gridspec(2, 3, hspace=0.52, wspace=0.45)
ax = fig.add_subplot(gs[0, :2]); panel(ax, 'a', -0.08)
ns = [24, 48, 96, 192, 384]
for m in M:
    mu = np.array([agg('main', 'iid', m, n)[0] for n in ns]); sd = np.array([agg('main', 'iid', m, n)[1] for n in ns])
    nn_ = [n if n != 384 else 384 for n in ns]
    ax.plot(ns, mu, 'o-', color=MC[m], ms=4, label=m); ax.fill_between(ns, mu - sd, mu + sd, color=MC[m], alpha=0.15, lw=0)
    table[m] = {n: agg('main', 'iid', m, n) for n in ns}
ax.set_xscale('log', base=2); ax.set_xticks(ns); ax.set_xticklabels(ns); ax.set_xlabel('labelled training runs')
ax.set_ylabel('test AUROC (mean ± s.d., 3 seeds)'); ax.set_title('Label efficiency (50 % of biomarkers unidentifiable)', loc='left')
ax.legend(fontsize=6.3, loc='lower right', ncol=2)
# b LOSO
ax = fig.add_subplot(gs[0, 2]); panel(ax, 'b', -0.3)
for j, m in enumerate(M):
    mu, sd, _ = agg('main', 'LOSO', m); mi, _, _ = agg('main', 'iid', m, 384)
    ax.bar(j, mu, 0.7, color=MC[m], yerr=sd, error_kw=dict(lw=0.8, capsize=2, ecolor=INK2))
    table.setdefault('LOSO', {})[m] = (mu, sd)
ax.set_xticks(range(5)); ax.set_xticklabels(['ID-1st', 'ID+MBR', 'ID-free', 'FF', 'FF+SSL'], fontsize=6.3, rotation=30)
ax.set_ylim(0.5, 1.02); ax.set_ylabel('AUROC on unseen site'); ax.set_title('Site transfer (LOSO)', loc='left')
# c sweep
ax = fig.add_subplot(gs[1, 0]); panel(ax, 'c', -0.3)
fr = [0.0, 0.25, 0.5, 0.75]
for m in ['ID-first-LR', 'ID-first+MBR-LR', 'ID-free-LR', 'FeatureFormer+SSL']:
    v = []
    for f in fr:
        mode = 'main' if f == 0.5 else 'sweep'
        v.append(agg(mode, 'iid', m, 192, f))
    ax.errorbar([f * 100 for f in fr], [x[0] for x in v], yerr=[x[1] for x in v], fmt='o-', color=MC[m], ms=4, capsize=2, lw=1.6, label=m)
    table.setdefault('sweep', {})[m] = v
ax.set_xlabel('biomarkers not identifiable (%)'); ax.set_ylabel('AUROC (n=192)'); ax.set_title('Dark-proteome sweep', loc='left')
ax.legend(fontsize=5.4, loc='lower left')
# d attribution
X = pickle.load(open('results/pheno2_seed0_extra.pkl', 'rb'))
ax = fig.add_subplot(gs[1, 1]); panel(ax, 'd', -0.3)
isb = np.isin(np.arange(len(X['attr'])), X['bio'])
a = X['attr'] / X['attr'].max()
ax.scatter(np.abs(X['effect'])[~isb] + np.random.default_rng(0).uniform(0, .03, (~isb).sum()), a[~isb], s=5, color='#b5b4ae', lw=0, label='background')
ax.scatter(np.abs(X['effect'])[isb & X['in_lib']], a[isb & X['in_lib']], s=16, color=C[0], label='biomarker (library)')
ax.scatter(np.abs(X['effect'])[isb & ~X['in_lib']], a[isb & ~X['in_lib']], s=16, marker='s', color=C[5], label='biomarker (unidentified)')
ax.set_xlabel('|true log₂ fold change|'); ax.set_ylabel('FeatureFormer attribution'); ax.set_title('Feature attribution', loc='left')
ax.legend(fontsize=5.6, loc='upper left')
table['attr_auroc'] = float(roc_auc_score(isb, X['attr']))
table['ttest_auroc'] = float(roc_auc_score(isb, -np.log10(np.nan_to_num(X['ttest_p'], nan=1.0))))
# e ROC
ax = fig.add_subplot(gs[1, 2]); panel(ax, 'e', -0.3)
for m in M:
    if m in X['probs']:
        fpr, tpr, _ = roc_curve(X['yte'][X['bte'] == 3], X['probs'][m]); ax.plot(fpr, tpr, color=MC[m], lw=1.5)
ax.plot([0, 1], [0, 1], ls='--', color=INK2, lw=0.8); ax.set_xlabel('false-positive rate'); ax.set_ylabel('true-positive rate')
ax.set_title('ROC, unseen site (seed 0)', loc='left')
fig.savefig('figures/fig6_phenotype.png'); plt.close(fig)
json.dump({k: ({str(kk): vv for kk, vv in v.items()} if isinstance(v, dict) else v) for k, v in table.items()},
          open('results/pheno_summary.json', 'w'), indent=1, default=float)
print(json.dumps({k: ({str(kk): vv for kk, vv in v.items()} if isinstance(v, dict) else v) for k, v in table.items()}, indent=0, default=float)[:3000])
