"""Real-data validation, step 2: calibrate the virtual instrument against real
annotated spectra (MassIVE-KB HCD sample; PRIDE-derived BSA CID ion-trap PSMs)."""
import re, pickle, json, numpy as np
from pyteomics import mgf
from chem import AA_MASS, PROTON, H2O, NH3, peptide_mass

CASA = '/home/user/bigbio/casanovo/sample_data/sample_preprocessed_spectra.mgf'


def parse_seq(s):
    s = s.replace('C[Carbamidomethyl]', 'C').replace('M[Oxidation]', 'm').replace('C+57.021', 'C').replace('M+15.995', 'm')
    s = s.replace('I', 'L')
    return s if re.fullmatch('[ACDEFGHKLMNPQRSTVWYm]+', s) else None


def load_massivekb():
    out = []
    for s in mgf.read(CASA):
        seq = parse_seq(s['params']['seq'])
        if seq is None: continue
        z = int(s['params']['charge'][0]); pmz = float(s['params']['pepmass'][0])
        mz = np.asarray(s['m/z array'], float); it = np.asarray(s['intensity array'], float)
        out.append(dict(seq=seq, z=z, pmz=pmz, pmass=(pmz - PROTON) * z, mz=mz, raw_it=it, it=np.sqrt(it / it.max())))
    return out


def load_bsa(fdr=0.01):
    D = pickle.load(open('realdata/bsa_psms.pkl', 'rb')); S = D['specs']
    sc = np.array([s['score'] for s in S]); dc = np.array([s['decoy'] for s in S])
    o = np.argsort(-sc); d = np.cumsum(dc[o]); t = np.cumsum(~dc[o])
    q = np.minimum.accumulate((d / np.maximum(t, 1))[::-1])[::-1]
    qv = np.empty(len(S)); qv[o] = q
    acc = [dict(s, seq=s['hit']) for s, v in zip(S, qv) if v <= fdr and not s['decoy']]
    qp1 = np.minimum.accumulate(((d + 1) / np.maximum(t, 1))[::-1])[::-1]
    return acc, dict(n_ms2=len(S), n_psm=len(acc), n_pep=len(set(a['seq'] for a in acc)),
                     n_psm_conservative_5pct=int(((qp1 <= 0.05) & ~dc[o]).sum()), qv=qv)


def stats(rec, tol, ppm=None):
    seq, mz, it, z = rec['seq'], rec['mz'], rec['raw_it'], rec['z']
    r = np.array([AA_MASS[a] for a in seq]); n = len(seq)
    b = np.cumsum(r)[:-1] + PROTON; y = np.cumsum(r[::-1])[:-1] + H2O + PROTON
    y_by_cleavage = y[::-1]                     # y ion produced by cleavage k (same index as b)
    def match(th):
        j = np.clip(np.searchsorted(mz, th), 1, len(mz) - 1)
        jj = np.where(np.abs(mz[j] - th) < np.abs(mz[j - 1] - th), j, j - 1)
        t = tol if ppm is None else th * ppm * 1e-6
        ok = np.abs(mz[jj] - th) < t
        return ok, jj, mz[jj] - th
    ob, jb, eb = match(b); oy, jy, ey = match(y_by_cleavage)
    ib = np.where(ob, it[jb], 0); iy = np.where(oy, it[jy], 0)
    expl = set(jb[ob]) | set(jy[oy])
    for th in (b - H2O, y_by_cleavage - H2O, b - NH3, y_by_cleavage - NH3, (b + PROTON) / 2, (y_by_cleavage + PROTON) / 2):
        o, j, _ = match(th); expl |= set(j[o])
    clv = (ib + iy)
    pro = np.array([seq[k + 1] == 'P' for k in range(n - 1)])
    lowmass = mz > 150
    errs = np.r_[eb[ob] / b[ob], ey[oy] / y_by_cleavage[oy]] * 1e6
    return dict(p_obs=float((ob.sum() + oy.sum()) / (2 * (n - 1))), p_obs_b=float(ob.mean()), p_obs_y=float(oy.mean()),
                y_over_b=float(iy.sum() / max(ib.sum(), 1e-9)) if ib.sum() > 0 else np.nan,
                pro_effect=float(clv[pro].mean() / max(clv[~pro].mean(), 1e-9)) if pro.any() and clv[~pro].mean() > 0 else np.nan,
                frac_tic_explained=float(it[list(expl)].sum() / it.sum()),
                n_peaks=int(len(mz)), n_unexplained=int(len(mz) - len(expl)),
                ppm_err=errs.tolist())


def summarize(rows):
    keys = ['p_obs', 'p_obs_b', 'p_obs_y', 'y_over_b', 'pro_effect', 'frac_tic_explained', 'n_peaks', 'n_unexplained']
    out = {k: [float(np.nanmedian([r[k] for r in rows])), float(np.nanpercentile([r[k] for r in rows], 25)),
               float(np.nanpercentile([r[k] for r in rows], 75))] for k in keys}
    e = np.concatenate([r['ppm_err'] for r in rows])
    out['ppm_err_sd'] = float(np.std(e)); out['ppm_err_median'] = float(np.median(e)); out['n'] = len(rows)
    return out


if __name__ == '__main__':
    mk = load_massivekb()
    bsa, bsa_info = load_bsa()
    # simulated reference (inDB test set, same metrics on the simulator's own output)
    D = pickle.load(open('data/ms2_corpus.pkl', 'rb'))
    sim = [dict(seq=x[4], z=x[3], mz=x[0].astype(float), raw_it=x[1].astype(float) ** 2) for x in D['tests']['inDB'][:500]]
    res = dict(massivekb=summarize([stats(r, None, ppm=20) for r in mk]),
               bsa=summarize([stats(r, 0.4) for r in bsa]),
               simulator=summarize([stats(r, None, ppm=20) for r in sim]),
               bsa_info={k: v for k, v in bsa_info.items() if k != 'qv'})
    json.dump(res, open('realdata/calibration.json', 'w'), indent=1)
    pickle.dump(dict(mk=mk, bsa=bsa, mk_rows=[stats(r, None, ppm=20) for r in mk], bsa_rows=[stats(r, 0.4) for r in bsa],
                     sim_rows=[stats(r, None, ppm=20) for r in sim]), open('realdata/calib_rows.pkl', 'wb'))
    print(json.dumps(res, indent=1))
