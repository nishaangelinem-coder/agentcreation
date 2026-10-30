"""Build the MS2 corpora: training corpus (proteome T) and benchmark test sets
(proteome P = searchable database; SAAV variants; novel proteome N)."""
import pickle, itertools
import numpy as np
from simulator import make_proteome, tryptic_digest, oxidise, mutate_saav, build_ms2_set

rng = np.random.default_rng(2026)
T = make_proteome(3000, rng)
P = make_proteome(1500, rng)
N = make_proteome(600, rng)
digest = lambda prots: sorted({p for pr in prots for p in tryptic_digest(pr)})
pepT, pepP, pepN = digest(T), digest(P), digest(N)
setP = set(pepP)
pepT = [p for p in pepT if p not in setP]
pepN = [p for p in pepN if p not in setP]
print('unique peptides  T=%d  P=%d  N=%d' % (len(pepT), len(pepP), len(pepN)))

# training corpus: 250k spectra over heterogeneous acquisition conditions
n_train = 250_000
train = []
for i in range(n_train):
    s = oxidise(pepT[rng.integers(len(pepT))], rng)
    train += build_ms2_set([s], rng, p_obs=rng.uniform(.55, .95),
                           noise_peaks=rng.uniform(10, 100), max_peaks=80)
val = build_ms2_set([oxidise(pepT[i], rng) for i in rng.integers(0, len(pepT), 2000)], rng, max_peaks=80)

def sample(peps, n):
    return [oxidise(peps[i], rng) for i in rng.choice(len(peps), n, replace=False)]

tests = {}
tests['inDB'] = build_ms2_set(sample(pepP, 3000), rng, max_peaks=80)
saav = [mutate_saav(s, rng) for s in sample(pepP, 2500)]
saav = [s for s in saav if s.replace('m', 'M') not in setP][:1500]
tests['SAAV'] = build_ms2_set(saav, rng, max_peaks=80)
tests['novel'] = build_ms2_set(sample(pepN, 1500), rng, max_peaks=80)
for po in (0.95, 0.8, 0.65, 0.5):
    tests['pobs_%.2f' % po] = build_ms2_set(sample(pepP, 1000), rng, p_obs=po, max_peaks=80)
for nz in (10, 40, 100, 200):
    tests['noise_%d' % nz] = build_ms2_set(sample(pepP, 1000), rng, noise_peaks=nz, max_peaks=80)

# search database: canonical P peptides with variable Met oxidation (<=3 Met)
db = set()
for p in pepP:
    idx = [i for i, a in enumerate(p) if a == 'M']
    if len(idx) > 3: idx = idx[:3]
    for k in range(len(idx) + 1):
        for comb in itertools.combinations(idx, k):
            q = list(p)
            for c in comb: q[c] = 'm'
            db.add(''.join(q))
db = sorted(db)
print('DB entries (targets):', len(db))
pickle.dump(dict(train=train, val=val, tests=tests, db=db,
                 n_pep=dict(T=len(pepT), P=len(pepP), N=len(pepN)),
                 n_prot=dict(T=len(T), P=len(P), N=len(N))),
            open('data/ms2_corpus.pkl', 'wb'))
