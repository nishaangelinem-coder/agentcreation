"""In-silico MS2 / LC-MS1 generator (the 'virtual instrument').

Generates a synthetic proteome, tryptic peptides, realistic-ish CID/HCD
fragment spectra (b/y ladders, doubly charged fragments, neutral losses,
proline / aspartate cleavage preferences, stochastic fragment dropout, ppm
mass error, chemical-noise peaks) and LC-MS1 maps for phenotype cohorts.
"""
import numpy as np
from chem import AA_MASS, AA_FREQ, PROTON, H2O, NH3, fragment_ladders, peptide_mass

_AA19 = list(AA_FREQ.keys())
_P19 = np.array(list(AA_FREQ.values())); _P19 /= _P19.sum()


def make_proteome(n_prot, rng, len_range=(250, 650)):
    return [''.join(rng.choice(_AA19, size=rng.integers(*len_range), p=_P19))
            for _ in range(n_prot)]


def tryptic_digest(prot, min_len=7, max_len=12, missed=1):
    sites = [i + 1 for i, a in enumerate(prot[:-1]) if a in 'KR' and prot[i + 1] != 'P']
    bounds = [0] + sites + [len(prot)]
    peps = []
    for i in range(len(bounds) - 1):
        for j in range(i + 1, min(i + 2 + missed, len(bounds))):
            p = prot[bounds[i]:bounds[j]]
            if min_len <= len(p) <= max_len:
                peps.append(p)
    return peps


def oxidise(seq, rng, p=0.25):
    return ''.join('m' if (a == 'M' and rng.random() < p) else a for a in seq)


def mutate_saav(seq, rng):
    """Single amino-acid variant (not K/R so trypticity is kept)."""
    pos = rng.integers(0, len(seq) - 1)
    choices = [a for a in _AA19 if a not in ('K', 'R', seq[pos])]
    return seq[:pos] + rng.choice(choices) + seq[pos + 1:]


def simulate_ms2(seq, rng, p_obs=0.8, noise_peaks=40, ppm=10.0, max_peaks=150):
    """Return (mz, intensity, precursor_mass, charge)."""
    n = len(seq)
    z = 2 if rng.random() < (0.85 if n <= 10 else 0.6) else 3
    b, y = fragment_ladders(seq)
    pos = np.arange(1, n)  # cleavage index i: between seq[i-1] and seq[i]
    # intensity model: y > b, enhanced cleavage N-term of Pro and C-term of Asp
    base = rng.lognormal(0.0, 0.6, size=(2, n - 1))
    fac = np.ones(n - 1)
    for k, i in enumerate(pos):
        if seq[i] == 'P': fac[k] *= 4.0
        if seq[i - 1] == 'D': fac[k] *= 2.0
        if seq[i - 1] == 'P': fac[k] *= 0.5
    edge = np.minimum(pos, n - pos)
    fac *= np.where(edge <= 1, 0.35, 1.0)
    ib = base[0] * fac * 0.6
    iy = base[1] * fac * 1.0
    mz, it = [], []
    for ladder, inten in ((b, ib), (y, iy)):
        keep = rng.random(n - 1) < p_obs
        mz.append(ladder[keep]); it.append(inten[keep])
        nl = keep & (rng.random(n - 1) < 0.35)
        mz.append(ladder[nl] - H2O); it.append(inten[nl] * rng.uniform(.05, .3, nl.sum()))
        nl = keep & (rng.random(n - 1) < 0.2)
        mz.append(ladder[nl] - NH3); it.append(inten[nl] * rng.uniform(.05, .2, nl.sum()))
        if z == 3:
            dk = keep & (pos > 3) & (rng.random(n - 1) < 0.5)
            mz.append((ladder[dk] + PROTON) / 2); it.append(inten[dk] * rng.uniform(.2, .6, dk.sum()))
    mz = np.concatenate(mz); it = np.concatenate(it)
    mz = mz * (1 + rng.normal(0, ppm * 1e-6, mz.size))
    pmass = peptide_mass(seq)
    k = rng.poisson(noise_peaks)
    if k:
        mz = np.concatenate([mz, rng.uniform(100, pmass, k)])
        it = np.concatenate([it, rng.lognormal(-1.6, 0.8, k)])
    o = np.argsort(-it)[:max_peaks]
    mz, it = mz[o], it[o]
    it = np.sqrt(it / it.max())
    o = np.argsort(mz)
    pm = pmass * (1 + rng.normal(0, 3e-6))
    return mz[o].astype(np.float32), it[o].astype(np.float32), pm, z


def build_ms2_set(peptides, rng, **kw):
    out = []
    for s in peptides:
        mz, it, pm, z = simulate_ms2(s, rng, **kw)
        out.append((mz, it, pm, z, s))
    return out


# ---------------------------------------------------------------- LC-MS1 maps
class CohortSimulator:
    """LC-MS1 maps (RT x m/z images) for case/control cohorts.

    Features are peptide ions. A fraction `frac_unid` of the biomarker signal
    is carried by ions that are *not* present in the identification library
    (PTM-modified, variant or non-canonical peptides) and are therefore
    invisible to an identify-then-quantify pipeline.
    """
    RT_BINS, MZ_BINS = 64, 256
    RT_RANGE, MZ_RANGE = (0., 60.), (400., 1200.)

    def __init__(self, rng, n_feat=600, frac_lib=0.7, n_bio=24, frac_unid=0.5,
                 log2fc=0.8, n_batches=4):
        self.rng = rng
        self.n_feat = n_feat
        self.mz = rng.uniform(*self.MZ_RANGE, n_feat)
        self.rt = rng.uniform(3, 57, n_feat)
        self.z = rng.choice([2, 3], n_feat, p=[.7, .3])
        self.base = rng.normal(20, 1.6, n_feat)        # log2 abundance
        self.bio_sd = rng.uniform(0.2, 0.45, n_feat)    # biological variance
        n_lib = int(frac_lib * n_feat)
        self.in_lib = np.zeros(n_feat, bool); self.in_lib[:n_lib] = True
        n_unid = int(round(frac_unid * n_bio))
        bio_unid = rng.choice(np.arange(n_lib, n_feat), n_unid, replace=False)
        bio_id = rng.choice(np.arange(n_lib), n_bio - n_unid, replace=False)
        self.bio = np.concatenate([bio_id, bio_unid]).astype(int)
        self.effect = np.zeros(n_feat)
        self.effect[self.bio] = rng.choice([-1, 1], n_bio) * rng.uniform(.6, 1.0, n_bio) * log2fc / 0.8
        # batch effects: RT drift, global intensity, m/z calibration
        self.batch_rt = rng.normal(0, 0.6, n_batches)
        self.batch_int = rng.normal(0, 0.5, n_batches)
        self.batch_ppm = rng.normal(0, 4, n_batches)
        self.batch_feat = rng.normal(0, 0.25, (n_batches, n_feat))
        self.n_batches = n_batches
        self.lod = 17.5

    def sample_ab(self, n, prevalence=0.5):
        rng = self.rng
        y = (rng.random(n) < prevalence).astype(int)
        b = rng.integers(0, self.n_batches, n)
        ab = (self.base[None] + rng.normal(0, 1, (n, self.n_feat)) * self.bio_sd[None]
              + y[:, None] * self.effect[None] + self.batch_int[b][:, None] + self.batch_feat[b])
        return y, b, ab

    def sample(self, n, prevalence=0.5, batches=None):
        rng = self.rng
        y = (rng.random(n) < prevalence).astype(int)
        b = rng.integers(0, self.n_batches, n) if batches is None else rng.choice(batches, n)
        ab = (self.base[None] + rng.normal(0, 1, (n, self.n_feat)) * self.bio_sd[None]
              + y[:, None] * self.effect[None] + self.batch_int[b][:, None]
              + self.batch_feat[b])
        imgs = np.stack([self.render(ab[i], b[i]) for i in range(n)])
        return imgs.astype(np.float32), y, b, ab

    def render(self, ab, b):
        rng = self.rng
        R, M = self.RT_BINS, self.MZ_BINS
        img = np.zeros((R, M))
        rt = self.rt + self.batch_rt[b] + rng.normal(0, 0.15, self.n_feat)
        mz = self.mz * (1 + self.batch_ppm[b] * 1e-6)
        rt_ax = np.linspace(*self.RT_RANGE, R)
        width = 0.45
        inten = 2.0 ** (ab - 20)
        inten[ab < self.lod + rng.normal(0, .3, self.n_feat)] = 0  # missing (below LOD)
        g = np.exp(-0.5 * ((rt_ax[None] - rt[:, None]) / width) ** 2)  # F x R
        mzb = (mz - self.MZ_RANGE[0]) / (self.MZ_RANGE[1] - self.MZ_RANGE[0]) * M
        for iso, rel in enumerate((1.0, 0.55, 0.2)):
            pos = mzb + iso * (1.00336 / self.z) / (800 / M)
            lo = np.clip(np.floor(pos).astype(int), 0, M - 1)
            fr = pos - np.floor(pos)
            hi = np.clip(lo + 1, 0, M - 1)
            np.add.at(img.T, lo, (g * (inten * rel * (1 - fr))[:, None]))
            np.add.at(img.T, hi, (g * (inten * rel * fr)[:, None]))
        img += rng.gamma(1.0, 0.02, img.shape)  # chemical noise
        return np.log1p(img * 10)

    def lib_quant(self, imgs, b=None):
        """Identify-then-quantify: integrate library features at their
        *library* coordinates (fixed RT/m/z extraction windows)."""
        R, M = self.RT_BINS, self.MZ_BINS
        idx = np.where(self.in_lib)[0]
        ri = np.clip(((self.rt[idx] - self.RT_RANGE[0]) / 60 * (R - 1)).round().astype(int), 0, R - 1)
        mi = np.clip(((self.mz[idx] - 400) / 800 * M).astype(int), 0, M - 1)
        q = np.zeros((len(imgs), len(idx)))
        for d in (-1, 0, 1):
            q += imgs[:, np.clip(ri + d, 0, R - 1), mi]
        return q, idx


# ------------------------------------------------ feature-level MS1 acquisition
def acquire_features(sim, ab, b, rng, n_noise=80, p_id=0.8):
    """Simulate MS1 feature detection for one run.

    Returns observed (mz, rt, log2 intensity, true_index) for detected ions plus
    spurious features (true_index=-1), and the set of library features that
    were *identified* by data-dependent MS2 (stochastic, prob p_id)."""
    det = ab > sim.lod + rng.normal(0, .3, sim.n_feat)
    idx = np.where(det)[0]
    mz = sim.mz[idx] * (1 + (sim.batch_ppm[b] + rng.normal(0, 2, idx.size)) * 1e-6)
    rt = sim.rt[idx] + sim.batch_rt[b] + rng.normal(0, 0.15, idx.size)
    it = ab[idx] + rng.normal(0, 0.15, idx.size)
    k = rng.poisson(n_noise)
    mz = np.r_[mz, rng.uniform(*sim.MZ_RANGE, k)]; rt = np.r_[rt, rng.uniform(1, 59, k)]
    it = np.r_[it, rng.normal(18.3, .8, k)]; tid = np.r_[idx, -np.ones(k, int)]
    ident = det & sim.in_lib & (rng.random(sim.n_feat) < p_id)
    return mz, rt, it, tid, ident


def align_to_reference(mz, rt, it, ref_mz, ref_rt, ppm=15, rt_tol=0.7):
    """ID-free alignment: global RT-shift estimation + nearest-neighbour
    matching of detected features to a consensus feature map."""
    o = np.argsort(ref_mz); rm, rr = ref_mz[o], ref_rt[o]
    def match(shift, tol):
        out = np.full(len(rm), np.nan); best = np.full(len(rm), np.inf); d_all = []
        lo = np.searchsorted(rm, mz * (1 - ppm * 1e-6)); hi = np.searchsorted(rm, mz * (1 + ppm * 1e-6))
        for i in range(len(mz)):
            if hi[i] <= lo[i]: continue
            c = np.arange(lo[i], hi[i]); d = rt[i] - shift - rr[c]
            j = np.argmin(np.abs(d))
            if abs(d[j]) < tol:
                d_all.append(d[j] + shift)
                cost = abs(d[j]) + abs(mz[i] - rm[c[j]]) / rm[c[j]] * 1e5 * 0.01
                if cost < best[c[j]]: best[c[j]] = cost; out[c[j]] = it[i]
        return out, np.array(d_all)
    _, d = match(0.0, 3.0)
    shift = float(np.median(d)) if len(d) else 0.0
    v, _ = match(shift, rt_tol)
    res = np.full(len(rm), np.nan); res[o] = v
    return res, shift
