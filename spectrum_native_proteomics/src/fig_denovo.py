"""Fig. 4 (de novo benchmark) and Fig. 5 (representation & interpretability)."""
import json, pickle, numpy as np
from figstyle import *
from sklearn.manifold import TSNE
from sklearn.metrics import silhouette_score
from sklearn.preprocessing import normalize

E = json.load(open('results/denovo_eval_ev.json'))
abl = {k: json.load(open(f'results/denovo_eval_{k}.json'))['inDB']['SpecNovo'] for k in ('base', 'ma', 'ev400', 'ev800')}
logs = {k: json.load(open(f'results/trainlog_{k}.json')) for k in ('main', 'ma', 'ev')}
M = {'DB-search': C[3], 'SpectrumGraph-DP': C[2], 'SpecNovo': C[0]}
summary = {}

fig = plt.figure(figsize=(7.2, 6.4))
gs = fig.add_gridspec(2, 3, hspace=0.5, wspace=0.45)
# a learning curves
ax = fig.add_subplot(gs[0, 0]); panel(ax, 'a', -0.3)
for k, lab, c, off in [('main', 'base', C[6], 0), ('ma', '+ mass-aware', C[1], 0), ('ev', '+ ion evidence', C[0], 2000 * 128)]:
    L = logs[k]
    if k == 'ma': L = [r for r in L if r['step'] <= 2000]
    x = np.array([r['seen'] + off for r in L]) / 1e3; yv = [r['val_tok_acc'] for r in L]
    ax.plot(x, yv, color=c, label=lab, lw=1.8)
ax.set_xlabel('training spectra seen (×10³)'); ax.set_ylabel('val. token accuracy\n(teacher forcing)')
ax.legend(fontsize=6.3, loc='lower right'); ax.set_title('Training dynamics', loc='left')
# b precision-coverage
ax = fig.add_subplot(gs[0, 1]); panel(ax, 'b', -0.3)
for meth in ('SpectrumGraph-DP', 'SpecNovo'):
    cov, pr = E['inDB'][meth]['curve']; ax.plot(cov, pr, color=M[meth], label=f"{meth} (AUC {E['inDB'][meth]['pr_auc']:.2f})")
ax.set_xlabel('coverage (fraction of spectra)'); ax.set_ylabel('peptide precision'); ax.set_ylim(0, 1.02)
ax.legend(fontsize=6, loc='lower left'); ax.set_title('Precision vs coverage', loc='left')
# c database dependence
ax = fig.add_subplot(gs[0, 2]); panel(ax, 'c', -0.3)
sets = [('inDB', 'canonical'), ('SAAV', 'SAAV'), ('novel', 'novel\nproteome')]
w = 0.26
for j, meth in enumerate(M):
    vals = [E[s][meth]['peptide_recall'] for s, _ in sets]
    ax.bar(np.arange(3) + (j - 1) * w, vals, w * 0.92, color=M[meth], label=meth)
    if meth == 'DB-search':
        fi = [E[s][meth]['false_ids'] for s, _ in sets]
        ax.bar(np.arange(3) + (j - 1) * w, fi, w * 0.92, bottom=vals, color='white', edgecolor=C[7], hatch='////', lw=0.8,
               label='DB false IDs (1 % FDR)')
ax.set_xticks(range(3)); ax.set_xticklabels([l for _, l in sets]); ax.set_ylabel('peptide recall'); ax.set_ylim(0, 1.15)
ax.legend(fontsize=5.6, loc='upper right'); ax.set_title('Database dependence', loc='left')
# d,e robustness
for p, (pref, xs, xl, lab) in zip(['d', 'e'], [('pobs_', ['0.95', '0.80', '0.65', '0.50'], 'fragment observation prob.', 'Fragment completeness'),
                                         ('noise_', ['10', '40', '100', '200'], 'mean noise peaks / spectrum', 'Chemical noise')]):
    ax = fig.add_subplot(gs[1, 0 if p == 'd' else 1]); panel(ax, p, -0.3)
    for meth in M:
        v = [E[pref + x][meth]['peptide_recall'] for x in xs]
        ax.plot([float(x) for x in xs], v, 'o-', color=M[meth], ms=4, label=meth)
        summary[f'{pref}{meth}'] = v
    ax.set_xlabel(xl); ax.set_ylabel('peptide recall'); ax.set_ylim(0, 1.02); ax.set_title(lab, loc='left')
    if p == 'd': ax.invert_xaxis(); ax.legend(fontsize=6, loc='lower left')
# f ablation bars
ax = fig.add_subplot(gs[1, 2]); panel(ax, 'f', -0.3)
names = ['base', '+MA', '+MA+EV\n+51k', '+MA+EV\n+102k', '+MA+EV\n+154k']
keys = [abl['base'], abl['ma'], abl['ev400'], abl['ev800'], E['inDB']['SpecNovo']]
x = np.arange(len(keys))
ax.bar(x - 0.2, [k['aa_precision'] for k in keys], 0.38, color=C[0], label='AA precision')
ax.bar(x + 0.2, [k['peptide_recall'] for k in keys], 0.38, color=C[1], label='peptide recall')
ax.set_xticks(x); ax.set_xticklabels(names, fontsize=5.4, rotation=35, ha='right'); ax.set_ylim(0, 1.05); ax.legend(fontsize=6, loc='upper left')
ax.set_title('Ablation (in-DB)', loc='left')
fig.savefig('figures/fig4_denovo.png'); plt.close(fig)

# ---------------- Fig 5
A = pickle.load(open('results/denovo_analysis_ev.pkl', 'rb'))
emb = normalize(A['emb']); binv = normalize(A['binv']); lab = A['labels']
def r1(X):
    S = X @ X.T; np.fill_diagonal(S, -np.inf); return float(np.mean(lab[S.argmax(1)] == lab))
summary['retrieval_top1_encoder'] = r1(emb); summary['retrieval_top1_binned'] = r1(binv)
summary['silhouette_encoder'] = float(silhouette_score(emb, lab, metric='cosine'))
summary['silhouette_binned'] = float(silhouette_score(binv, lab, metric='cosine'))
summary['saliency_on_fragments_mean'] = float(np.mean(A['sal_frag'])); summary['intensity_on_fragments_mean'] = float(np.mean(A['int_frag']))
fig = plt.figure(figsize=(7.2, 6.0)); gs = fig.add_gridspec(2, 3, hspace=0.62, wspace=0.45)
for k, (X, ttl, pl) in enumerate([(emb, 'SpecNovo encoder embedding', 'a'), (binv, 'Binned spectrum (0.5 Da)', 'b')]):
    ax = fig.add_subplot(gs[0, k]); panel(ax, pl, -0.12)
    Z = TSNE(perplexity=15, random_state=0, init='pca').fit_transform(X)
    hl = lab < 6
    ax.scatter(Z[~hl, 0], Z[~hl, 1], s=5, color='#bdbcb6', lw=0)
    for q in range(6):
        ax.scatter(Z[lab == q, 0], Z[lab == q, 1], s=12, color=C[q], lw=0.4, edgecolor='white', label=A['peps'][q])
    ax.set_xticks([]); ax.set_yticks([]); ax.grid(False); ax.set_title(ttl, loc='left', fontsize=8)
    if k == 0: fig.legend(*ax.get_legend_handles_labels(), fontsize=6, loc='upper center', bbox_to_anchor=(0.36, 0.535), ncol=6, handletextpad=0.1, columnspacing=0.8)
ax = fig.add_subplot(gs[0, 2]); panel(ax, 'c', -0.3)
vals = [summary['retrieval_top1_binned'], summary['retrieval_top1_encoder'], summary['silhouette_binned'], summary['silhouette_encoder']]
ax.bar([0, 1], vals[:2], 0.7, color=[C[6], C[0]]); ax.bar([2.5, 3.5], vals[2:], 0.7, color=[C[6], C[0]])
for xx, v in zip([0, 1, 2.5, 3.5], vals): ax.text(xx, v + 0.02, f'{v:.2f}', ha='center', fontsize=6.5)
ax.set_xticks([0.5, 3]); ax.set_xticklabels(['top-1 replicate\nretrieval', 'silhouette\n(cosine)'], fontsize=6.5)
from matplotlib.patches import Patch
ax.legend(handles=[Patch(color=C[6], label='binned'), Patch(color=C[0], label='encoder')], fontsize=6, loc='upper right')
ax.set_ylim(min(0, min(vals)) - 0.05, 1.12); ax.set_title('Representation quality', loc='left')
ex = A['examples'][0]
ax = fig.add_subplot(gs[1, :2]); panel(ax, 'd', -0.07)
s = ex['sal'] / ex['sal'].max()
for m, it, f, sv in zip(ex['mz'], ex['it'], ex['isf'], s):
    ax.vlines(m, 0, it, color=C[0] if f else '#9c9b95', lw=1.0)
    ax.vlines(m, 0, -sv, color=C[1] if f else '#9c9b95', lw=1.0)
ax.axhline(0, color=INK2, lw=0.6); ax.set_ylim(-1.1, 1.1); ax.grid(False)
ax.text(0.01, 0.95, 'intensity (fragment ions blue, other grey)', transform=ax.transAxes, fontsize=6.5)
ax.text(0.01, 0.04, 'gradient×input saliency (fragment ions orange)', transform=ax.transAxes, fontsize=6.5)
ax.set_xlabel('m/z'); ax.set_title(f'Peak attribution: {ex["seq"]}', loc='left')
ax = fig.add_subplot(gs[1, 2]); panel(ax, 'e', -0.3)
ax.scatter(A['int_frag'], A['sal_frag'], s=6, color=C[0], alpha=0.6, lw=0)
ax.plot([0, 1], [0, 1], color=INK2, lw=1, ls='--'); ax.set_xlim(0, 1); ax.set_ylim(0, 1)
ax.set_xlabel('intensity share on fragments'); ax.set_ylabel('saliency share on fragments'); ax.set_title('Attribution focus', loc='left')
fig.savefig('figures/fig5_representation.png'); plt.close(fig)
json.dump(summary, open('results/denovo_summary.json', 'w'), indent=1)
print(json.dumps(summary, indent=1))
