"""Calibrated MS2 corpus: simulator parameters re-fitted to real HCD statistics
(MassIVE-KB: b-ion observation 0.29, y-ion 0.67, y/b intensity 2.7, ~39
unexplained peaks, 4.4 ppm fragment error). Peptides 6-14 aa from proteome T."""
import pickle, numpy as np
from simulator import make_proteome, tryptic_digest, oxidise, simulate_ms2

rng = np.random.default_rng(2026)
T = make_proteome(3000, rng)                     # identical to data_gen.py (same seed, first call)
P = make_proteome(1500, rng)
setP = {p for pr in P for p in tryptic_digest(pr, 6, 14)}
pep = sorted({p for pr in T for p in tryptic_digest(pr, 6, 14)} - setP)
rng = np.random.default_rng(7)
out = []
for i in range(80000):
    s = oxidise(pep[rng.integers(len(pep))], rng)
    mz, it, pm, z = simulate_ms2(s, rng, noise_peaks=rng.uniform(20, 80), ppm=5.0, max_peaks=80,
                                 p_obs_b=rng.uniform(0.1, 0.6), p_obs_y=rng.uniform(0.4, 0.95),
                                 b_factor=rng.uniform(0.2, 0.6), low_mass_ions=10)
    out.append((mz, it, pm, z, s))
D = pickle.load(open('data/ms2_corpus.pkl', 'rb'))
pickle.dump(dict(train=out, val=D['val']), open('data/ms2_calibrated.pkl', 'wb'))
print(len(pep), len(out))
