"""Real-data validation, step 1: target-decoy database search of PRIDE-derived
LC-MS/MS runs (bigbio quantms CI 'BSA' runs, LTQ Orbitrap Velos, CID ion-trap MS2)
and of MassIVE-KB annotated spectra (Casanovo sample data).

Outputs realdata/bsa_psms.pkl: one record per MS2 spectrum with the best target /
decoy hit, q-value and the peak list, plus MS1 scans for feature analysis."""
import sys, glob, pickle, re, numpy as np
from math import lgamma, log
from pyteomics import mzml, fasta, parser
from chem import AA_MASS, PROTON, H2O, peptide_mass

BSA = '/home/user/bigbio/quantms-test-datasets/testdata/lfq_ci/BSA'
FRAG_TOL = 0.4          # ion-trap CID MS2
QS = (4, 8, 12)
PREC_PPM = 10.0


def digest_db(path, min_len=7, max_len=30):
    targets = set()
    for desc, seq in fasta.read(path):
        seq = seq.replace('I', 'L').upper()
        if set(seq) - set('ACDEFGHKLMNPQRSTVWY'): seq = re.sub('[^ACDEFGHKLMNPQRSTVWY]', '', seq)
        for p in parser.cleave(seq, parser.expasy_rules['trypsin'], missed_cleavages=2, min_length=min_len, max_length=max_len):
            targets.add(p)
    peps = set()
    for p in targets:                       # variable Met oxidation (<=2)
        peps.add(p)
        idx = [i for i, a in enumerate(p) if a == 'M'][:2]
        for i in idx: peps.add(p[:i] + 'm' + p[i + 1:])
        if len(idx) == 2: q = list(p); q[idx[0]] = 'm'; q[idx[1]] = 'm'; peps.add(''.join(q))
    tset = sorted(peps)
    dec = [p[:-1][::-1] + p[-1] for p in tset]
    allp = tset + dec
    isd = np.r_[np.zeros(len(tset), bool), np.ones(len(dec), bool)]
    mass = np.array([peptide_mass(p) for p in allp])
    o = np.argsort(mass)
    return [allp[i] for i in o], mass[o], isd[o]


def theo(seq, z):
    r = np.array([AA_MASS[a] for a in seq])
    b = np.cumsum(r)[:-1] + PROTON; y = np.cumsum(r[::-1])[:-1] + H2O + PROTON
    ions = [b, y]
    if z >= 3: ions += [(b + PROTON) / 2, (y + PROTON) / 2]
    return ions


from scipy.stats import binom
TOPQ = 6


def topq_filter(mz, it, pmz, q=TOPQ):
    """Andromeda-style preprocessing: keep the q most intense peaks per 100 m/z
    window and remove the unfragmented precursor region."""
    keep = np.abs(mz - pmz) > 2.0
    mz, it = mz[keep], it[keep]
    out = []
    for w in np.unique((mz // 100).astype(int)):
        idx = np.where((mz // 100).astype(int) == w)[0]
        out.extend(idx[np.argsort(-it[idx])[:q]])
    out = np.sort(np.array(out, int))
    return mz[out], it[out]


def binom_score(seq, mz, z, tol=FRAG_TOL, q=TOPQ):
    ions = theo(seq, z)
    n = sum(len(t) for t in ions); k = 0
    for th in ions:
        j = np.clip(np.searchsorted(mz, th), 1, len(mz) - 1)
        d = np.minimum(np.abs(mz[j] - th), np.abs(mz[j - 1] - th)); k += int((d < tol).sum())
    p = q * 2 * tol / 100.0
    return float(-np.log10(binom.sf(k - 1, n, p) + 1e-300)) if k else 0.0


def hyperscore(seq, mz, it, z, tol=FRAG_TOL):
    sc, n = 0.0, []
    for th in theo(seq, z):
        j = np.clip(np.searchsorted(mz, th), 1, len(mz) - 1)
        d = np.minimum(np.abs(mz[j] - th), np.abs(mz[j - 1] - th))
        jj = np.where(np.abs(mz[j] - th) < np.abs(mz[j - 1] - th), j, j - 1)
        hit = d < tol; sc += it[jj[hit]].sum(); n.append(hit.sum())
    nb = n[0] + (n[2] if len(n) > 2 else 0); ny = n[1] + (n[3] if len(n) > 3 else 0)
    return (log(sc + 1e-9) + lgamma(nb + 1) + lgamma(ny + 1)) if sc > 0 else 0.0


def search(spec, peps, mass, isd, tol=FRAG_TOL):
    mz, it, pm, z = spec['mz'], spec['it'], spec['pmass'], spec['z']
    best = dict(t=(None, 0.0), d=(None, 0.0))
    for iso in (0, -1, 1):
        m = pm - iso * 1.003355
        lo, hi = np.searchsorted(mass, [m * (1 - PREC_PPM * 1e-6), m * (1 + PREC_PPM * 1e-6)])
        for k in range(lo, hi):
            s = max(binom_score(peps[k], spec['fmz%d' % q], z, tol, q) for q in QS)
            key = 'd' if isd[k] else 't'
            if s > best[key][1]: best[key] = (peps[k], s)
    return best


def read_runs():
    specs, ms1 = [], []
    for f in sorted(glob.glob(BSA + '/*.mzML')):
        run = f.split('/')[-1][:-5]
        for s in mzml.read(f):
            rt = float(s['scanList']['scan'][0]['scan start time'])
            if s['ms level'] == 1:
                ms1.append(dict(run=run, rt=rt, mz=np.asarray(s['m/z array']), it=np.asarray(s['intensity array'])))
                continue
            pre = s['precursorList']['precursor'][0]['selectedIonList']['selectedIon'][0]
            z = int(pre.get('charge state', 2) or 2)
            pmz = float(s['scanList']['scan'][0].get('[Thermo Trailer Extra]Monoisotopic M/Z:', 0) or 0) or float(pre['selected ion m/z'])
            mz = np.asarray(s['m/z array'], float); it = np.asarray(s['intensity array'], float)
            if len(mz) < 10: continue
            rec = dict(run=run, rt=rt, id=s['id'], z=z, pmz=pmz, pmass=(pmz - PROTON) * z,
                       mz=mz, it=np.sqrt(it / it.max()), raw_it=it)
            for q in QS: rec['fmz%d' % q] = topq_filter(mz, it, pmz, q)[0]
            specs.append(rec)
    return specs, ms1


if __name__ == '__main__':
    specs, ms1 = read_runs()
    print('MS2 spectra', len(specs), 'MS1 scans', len(ms1), flush=True)
    peps, mass, isd = digest_db(BSA + '/18Protein_SoCe_Tr_detergents_trace.fasta')
    print('DB entries (target+decoy)', len(peps), flush=True)
    for i, s in enumerate(specs):
        b = search(s, peps, mass, isd)
        # target-decoy competition per spectrum
        if b['t'][1] >= b['d'][1]: s['hit'], s['score'], s['decoy'] = b['t'][0], b['t'][1], False
        else: s['hit'], s['score'], s['decoy'] = b['d'][0], b['d'][1], True
        if i % 500 == 0: print(i, flush=True)
    sc = np.array([s['score'] for s in specs]); dc = np.array([s['decoy'] for s in specs])
    o = np.argsort(-sc); d = np.cumsum(dc[o]); t = np.cumsum(~dc[o])
    q = np.minimum.accumulate(((d + 1) / np.maximum(t, 1))[::-1])[::-1]
    qv = np.empty(len(sc)); qv[o] = q
    for s, v in zip(specs, qv): s['q'] = float(v)
    acc = [(s['hit']) for s in specs if s['q'] <= 0.01 and not s['decoy']]
    print('PSMs @1% FDR', len(acc), 'unique peptides', len(set(acc)), 'ID rate', len(acc) / len(specs))
    pickle.dump(dict(specs=specs, ms1=ms1, n_db=len(peps)), open('realdata/bsa_psms.pkl', 'wb'))
