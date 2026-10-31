"""Fig. 7: validation on real PRIDE / MassIVE-KB data."""
import json, pickle, numpy as np
from figstyle import *
from chem import fragment_ladders
from matplotlib.patches import Patch

cal = json.load(open('realdata/calibration.json'))
R = pickle.load(open('realdata/calib_rows.pkl', 'rb'))
Z = json.load(open('realdata/denovo_zeroshot.json')); K = json.load(open('realdata/denovo_calibrated.json'))
B = pickle.load(open('realdata/bsa_psms.pkl', 'rb'))
M1 = json.load(open('realdata/ms1_bottleneck.json'))
fig = plt.figure(figsize=(7.2, 7.4)); gs = fig.add_gridspec(3, 3, hspace=0.75, wspace=0.55)

# a real annotated MassIVE-KB spectrum
ax = fig.add_subplot(gs[0, :2]); panel(ax, 'a', -0.08)
rec = max(R['mk'], key=lambda r: len(r['seq']) if 9 <= len(r['seq']) <= 12 else 0)
b, y = fragment_ladders(rec['seq']); it = rec['raw_it'] / rec['raw_it'].max()
for m, s in zip(rec['mz'], it):
    c = C[0] if np.min(np.abs(b - m)) < m * 2e-5 else (C[7] if np.min(np.abs(y - m)) < m * 2e-5 else '#9c9b95')
    ax.vlines(m, 0, s, color=c, lw=1.0 if c != '#9c9b95' else 0.7)
ax.set_xlabel('m/z'); ax.set_ylabel('rel. intensity'); ax.grid(False)
ax.set_title(f'Real HCD spectrum (MassIVE-KB): {rec["seq"]}, z={rec["z"]}', loc='left')
ax.legend(handles=[Patch(color=C[0], label='b ions'), Patch(color=C[7], label='y ions'), Patch(color='#9c9b95', label='unexplained')],
          fontsize=6, loc='upper right', ncol=3)
# b score distribution BSA
ax = fig.add_subplot(gs[0, 2]); panel(ax, 'b', -0.3)
S = B['specs']; sc = np.array([s['score'] for s in S]); dc = np.array([s['decoy'] for s in S])
bins = np.linspace(0, 14, 36)
ax.hist(sc[~dc], bins=bins, color=C[0], label='target', alpha=0.9); ax.hist(sc[dc], bins=bins, color=C[1], label='decoy', alpha=0.75)
ax.set_yscale('log'); ax.set_xlabel('binomial PSM score'); ax.set_ylabel('spectra'); ax.legend(fontsize=6)
ax.set_title('PRIDE BSA search', loc='left')
# c calibration
ax = fig.add_subplot(gs[1, :2]); panel(ax, 'c', -0.08)
keys = [('p_obs_b', 'b-ion obs.'), ('p_obs_y', 'y-ion obs.'), ('frac_tic_explained', 'TIC explained'), ('y_over_b', 'y/b intensity\n(÷4)'), ('pro_effect', 'Pro effect\n(÷4)')]
srcs = [('simulator', 'simulator (default)', C[6]), ('massivekb', 'MassIVE-KB HCD (real)', C[0]), ('bsa', 'PRIDE BSA CID (real)', C[1])]
w = 0.26
for j, (k, lab, c) in enumerate(srcs):
    for i, (key, _) in enumerate(keys):
        med, lo, hi = cal[k][key]; sc_ = 0.25 if key in ('y_over_b', 'pro_effect') else 1
        ax.bar(i + (j - 1) * w, med * sc_, w * 0.92, color=c, label=lab if i == 0 else None)
        ax.vlines(i + (j - 1) * w, lo * sc_, hi * sc_, color=INK2, lw=0.8)
ax.set_xticks(range(len(keys))); ax.set_xticklabels([l for _, l in keys], fontsize=6.5); ax.set_ylabel('median (IQR)')
ax.legend(fontsize=6, loc='upper left', ncol=3); ax.set_ylim(0, 1.25); ax.set_title('Simulator calibration against real spectra', loc='left')
# d de novo on real
ax = fig.add_subplot(gs[1, 2]); panel(ax, 'd', -0.3)
cats = [('Spectrum-\ngraph DP', Z['massivekb']['SpectrumGraph-DP'], C[2]), ('SpecNovo\nzero-shot', Z['massivekb']['SpecNovo'], C[6]),
        ('SpecNovo\ncalibrated', K['massivekb']['SpecNovo'], C[0])]
x = np.arange(3)
ax.bar(x - 0.2, [c[1]['aa_precision'] for c in cats], 0.38, color=[c[2] for c in cats], alpha=0.45)
ax.bar(x + 0.2, [c[1]['peptide_recall'] for c in cats], 0.38, color=[c[2] for c in cats])
ax.set_xticks(x); ax.set_xticklabels([c[0] for c in cats], fontsize=5.8)
ax.legend(handles=[Patch(color='#bbbbbb', label='AA precision (light)'), Patch(color='#555555', label='peptide recall (dark)')], fontsize=5.6, loc='upper left')
ax.set_ylim(0, 0.85); ax.set_title(f"De novo, real HCD (n={Z['n']['massivekb']})", loc='left')
# e MS1 funnel
ax = fig.add_subplot(gs[2, :2]); panel(ax, 'e', -0.08)
t = M1['total']
vals = [t['features'], t['with_ms2'], t['identified']]; labs = ['MS1 features\n(de-isotoped)', 'sampled for MS2', 'identified\n(1 % FDR)']
ax.barh(range(3)[::-1], vals, color=[C[2], C[3], C[0]], height=0.6)
for i, v in enumerate(vals): ax.text(v * 1.15, 2 - i, f'{v:,} ({v / vals[0] * 100:.2f} %)', va='center', fontsize=7)
ax.set_xscale('log'); ax.set_yticks(range(3)[::-1]); ax.set_yticklabels(labs, fontsize=7); ax.set_xlim(10, 3e5)
ax.set_xlabel('count (6 PRIDE BSA runs, log scale)'); ax.set_title('Identification bottleneck on real MS1 data', loc='left')
# f completeness
ax = fig.add_subplot(gs[2, 2]); panel(ax, 'f', -0.3)
A = M1['alignment']
ax.bar([0, 1], [A['idfree_complete'], A['idfirst_complete']], color=[C[2], C[3]], width=0.6)
for i, v in enumerate([A['idfree_complete'], A['idfirst_complete']]): ax.text(i, v * 1.1 + 1, f'{v:,}', ha='center', fontsize=7)
ax.set_yscale('log'); ax.set_ylim(0.5, 5e3); ax.set_xticks([0, 1]); ax.set_xticklabels(['ID-free', 'ID-first'], fontsize=7)
ax.set_ylabel('features quantified in\nall 3 samples'); ax.set_title('Complete profiles', loc='left')
fig.savefig('figures/fig7_real_validation.png'); print('ok')
