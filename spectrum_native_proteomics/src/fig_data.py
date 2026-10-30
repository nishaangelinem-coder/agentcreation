"""Fig. 2: the virtual instrument - annotated MS2 spectrum, LC-MS1 maps, corpus statistics."""
import pickle, numpy as np
from figstyle import *
from chem import fragment_ladders, PROTON, H2O, peptide_mass
from simulator import CohortSimulator, simulate_ms2

D = pickle.load(open('data/ms2_corpus.pkl', 'rb'))
rng = np.random.default_rng(11)
fig = plt.figure(figsize=(7.2, 6.8))
gs = fig.add_gridspec(3, 3, height_ratios=[1.05, 1.1, 0.95], hspace=0.62, wspace=0.62)

# a) annotated spectrum
ax = fig.add_subplot(gs[0, :2]); panel(ax, 'a', -0.07)
seq = 'LVNEDTPAGLK'
mz, it, pm, z = simulate_ms2(seq, rng, p_obs=0.85, noise_peaks=40)
b, y = fragment_ladders(seq)
cols = np.array([INK2] * len(mz), dtype=object); lab = [None] * len(mz)
for k, m in enumerate(mz):
    db = np.abs(b - m); dy = np.abs(y - m)
    if db.min() < 0.02: cols[k] = C[0]; lab[k] = f'b{db.argmin() + 1}'
    elif dy.min() < 0.02: cols[k] = C[7]; lab[k] = f'y{dy.argmin() + 1}'
for m, s, c in zip(mz, it, cols):
    ax.vlines(m, 0, s, color=c, lw=1.1 if c != INK2 else 0.7, alpha=1 if c != INK2 else 0.55)
for m, s, l in zip(mz, it, lab):
    if l and s > 0.25: ax.text(m, s + 0.03, l, ha='center', fontsize=6, color=INK)
ax.set_xlim(100, pm + 20); ax.set_ylim(0, 1.18); ax.set_xlabel('m/z'); ax.set_ylabel('rel. intensity (√)')
ax.set_title(f'Simulated MS2: {seq} (z={z}, M={pm:.2f} Da)', loc='left')
ax.plot([], [], color=C[0], label='b ions'); ax.plot([], [], color=C[7], label='y ions')
ax.plot([], [], color=INK2, alpha=.55, label='noise / neutral loss'); ax.legend(loc='upper right', fontsize=6.5, ncol=3)
ax.grid(False)

# b) corpus stats: peptide length + precursor mass
ax = fig.add_subplot(gs[0, 2]); panel(ax, 'b', -0.25)
L = np.array([len(x[4]) for x in D['train'][:50000]]); Lt = np.array([len(x[4]) for x in D['tests']['inDB']])
bins = np.arange(6.5, 13.5)
ax.hist(L, bins=bins, density=True, color=C[0], rwidth=0.8, label='train (T)')
ax.hist(Lt, bins=bins, density=True, histtype='step', color=C[1], lw=2, label='test (P)')
ax.set_xlabel('peptide length (aa)'); ax.set_ylabel('density'); ax.set_ylim(0, 0.27); ax.legend(fontsize=6.3, loc='upper right', ncol=2)

# c-e) LC-MS1 maps
sim = CohortSimulator(np.random.default_rng(3), n_feat=300, n_bio=20, frac_unid=0.5, log2fc=1.0)
X, yy, bb, ab = sim.sample(120)
ext = [400, 1200, 60, 0]
ax = fig.add_subplot(gs[1, 0]); panel(ax, 'c', -0.3)
im = ax.imshow(X[0], aspect='auto', extent=ext, cmap='magma'); ax.grid(False)
ax.set_xlabel('m/z'); ax.set_ylabel('retention time (min)'); ax.set_title('LC-MS1 map (one run)', loc='left')
cb = fig.colorbar(im, ax=ax, fraction=0.05, pad=0.02); cb.set_label('log(1+I)', fontsize=6.5); cb.ax.tick_params(labelsize=6)
ax = fig.add_subplot(gs[1, 1]); panel(ax, 'd', -0.3)
diff = X[yy == 1].mean(0) - X[yy == 0].mean(0)
v = np.abs(diff).max() * 0.6
im = ax.imshow(diff, aspect='auto', extent=ext, cmap='RdBu_r', vmin=-v, vmax=v); ax.grid(False)
lib = np.isin(sim.bio, np.where(sim.in_lib)[0])
ax.scatter(sim.mz[sim.bio[lib]], sim.rt[sim.bio[lib]], s=26, facecolors='none', edgecolors=INK, lw=0.9, label='library biomarker')
ax.scatter(sim.mz[sim.bio[~lib]], sim.rt[sim.bio[~lib]], s=26, marker='s', facecolors='none', edgecolors=C[5], lw=0.9, label='unidentified biomarker')
ax.set_xlabel('m/z'); ax.set_title('Case − control', loc='left')
cb = fig.colorbar(im, ax=ax, fraction=0.05, pad=0.02); cb.ax.tick_params(labelsize=6)
ax = fig.add_subplot(gs[1, 2]); panel(ax, 'e', -0.3)
for k in range(sim.n_batches):
    ax.plot(np.linspace(0, 60, 64), X[bb == k].mean((0, 2)), color=C[k], lw=1.5, label=f'site {k + 1} (ΔRT {sim.batch_rt[k]:+.2f} min)')
ax.set_xlabel('retention time (min)'); ax.set_ylabel('TIC (a.u.)'); ax.set_title('Site effects (4 sites)', loc='left')

# f) feature map: library vs non-library
ax = fig.add_subplot(gs[2, 0:2]); panel(ax, 'f', -0.07)
ax.scatter(sim.mz[sim.in_lib], sim.rt[sim.in_lib], s=5, color=C[0], alpha=0.6, label='in identification library (70 %)', lw=0)
ax.scatter(sim.mz[~sim.in_lib], sim.rt[~sim.in_lib], s=5, color=C[3], alpha=0.8, label='not identifiable (PTM / variant / non-canonical)', lw=0)
ax.set_xlabel('m/z'); ax.set_ylabel('RT (min)'); ax.set_title('')
ax.legend(fontsize=6.3, loc='upper right', ncol=2, markerscale=2.5, frameon=True, framealpha=0.92, edgecolor='none'); ax.set_ylim(0, 74)
# g) abundance distribution
ax = fig.add_subplot(gs[2, 2]); panel(ax, 'g', -0.25)
ax.hist(ab.ravel(), bins=50, color=C[2], density=True)
ax.axvline(sim.lod, color=C[7], lw=1.5, ls='--'); ax.text(sim.lod - 0.2, ax.get_ylim()[1] * 0.85, 'LOD', color=C[7], ha='right', fontsize=7)
ax.set_xlabel('log₂ abundance'); ax.set_ylabel('density'); ax.set_title('Dynamic range', loc='left')
fig.savefig('figures/fig2_virtual_instrument.png'); print('ok')
