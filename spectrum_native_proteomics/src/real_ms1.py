"""Real-data validation, step 4: the identification bottleneck on real MS1 data.
Simple centroid-linking MS1 feature detection on the BSA runs, isotope grouping,
ID-free alignment across samples, and comparison with identification-first
coverage (1 % FDR PSMs)."""
import pickle, json, numpy as np, collections
from real_calib import load_bsa

PPM = 10.0


def detect(scans, min_scans=3, topn=400):
    pts = []
    for k, s in enumerate(scans):
        mz, it = s['mz'], s['it']
        o = np.argsort(-it)[:topn]
        for m, i in zip(mz[o], it[o]): pts.append((m, k, s['rt'], i))
    pts.sort()
    feats, cur = [], [pts[0]]
    def close(c):
        c.sort(key=lambda x: x[1])
        # split on scan gaps > 2
        grp = [c[0]]
        for p in c[1:]:
            if p[1] - grp[-1][1] > 2:
                emit(grp); grp = [p]
            else: grp.append(p)
        emit(grp)
    def emit(g):
        if len(set(p[1] for p in g)) >= min_scans:
            it = np.array([p[3] for p in g]); mz = np.array([p[0] for p in g]); rt = np.array([p[2] for p in g])
            feats.append(dict(mz=float((mz * it).sum() / it.sum()), rt=float(rt[it.argmax()]), rt_lo=float(rt.min()),
                              rt_hi=float(rt.max()), apex=float(it.max()), area=float(it.sum()), n=len(g)))
    for p in pts[1:]:
        if (p[0] - cur[-1][0]) / p[0] * 1e6 < PPM: cur.append(p)
        else: close(cur); cur = [p]
    close(cur)
    return feats


def deisotope(feats):
    """Drop features that are +1/+2 isotopes (z=1..4) of a co-eluting lighter feature."""
    fs = sorted(feats, key=lambda f: f['mz']); mz = np.array([f['mz'] for f in fs]); iso = np.zeros(len(fs), bool)
    for i, f in enumerate(fs):
        for z in (1, 2, 3, 4):
            for n in (1, 2):
                t = f['mz'] + n * 1.003355 / z
                j = np.searchsorted(mz, t * (1 - PPM * 1e-6))
                while j < len(fs) and mz[j] <= t * (1 + PPM * 1e-6):
                    if abs(fs[j]['rt'] - f['rt']) < 15 and fs[j]['apex'] < f['apex'] * 3: iso[j] = True
                    j += 1
    return [f for f, x in zip(fs, iso) if not x]


def match_ms2(feats, specs):
    mz = np.array([f['mz'] for f in feats]); o = np.argsort(mz); mzs = mz[o]
    for f in feats: f['ms2'] = False; f['id'] = None
    for s in specs:
        j = np.searchsorted(mzs, s['pmz'] * (1 - PPM * 1e-6)); k = np.searchsorted(mzs, s['pmz'] * (1 + PPM * 1e-6))
        for idx in o[j:k]:
            f = feats[idx]
            if f['rt_lo'] - 30 <= s['rt'] <= f['rt_hi'] + 30:
                f['ms2'] = True
                if s.get('ident'): f['id'] = s['seq']


if __name__ == '__main__':
    D = pickle.load(open('realdata/bsa_psms.pkl', 'rb'))
    acc, info = load_bsa()
    ident = {(a['run'], a['id']): a['seq'] for a in acc}
    for s in D['specs']:
        s['ident'] = (s['run'], s['id']) in ident; s['seq'] = ident.get((s['run'], s['id']))
    runs = sorted(set(s['run'] for s in D['ms1']))
    allf = {}
    res = dict(runs={})
    for r in runs:
        sc = [s for s in D['ms1'] if s['run'] == r]
        f = deisotope(detect(sc)); match_ms2(f, [s for s in D['specs'] if s['run'] == r])
        allf[r] = f
        res['runs'][r] = dict(ms1_scans=len(sc), features=len(f), with_ms2=int(sum(x['ms2'] for x in f)),
                              identified=int(sum(x['id'] is not None for x in f)),
                              ms2=int(sum(s['run'] == r for s in D['specs'])), psms=int(sum(s['ident'] for s in D['specs'] if s['run'] == r)))
        print(r, res['runs'][r], flush=True)
    # merge fractions per sample and align across the three samples (ID-free)
    samples = sorted(set(r.split('_')[0] for r in runs))
    S = {sm: [x for r in runs if r.startswith(sm) for x in allf[r]] for sm in samples}
    ref = S[samples[0]]
    cons = sorted(ref, key=lambda x: -x['apex'])
    # consensus = union of features from all samples, clustered at 10 ppm / 60 s
    pool = sorted([dict(x, sample=sm) for sm in samples for x in S[sm]], key=lambda x: x['mz'])
    clusters = []
    for x in pool:
        placed = False
        for c in clusters[-50:]:
            if abs(c['mz'] - x['mz']) / x['mz'] * 1e6 < PPM and abs(c['rt'] - x['rt']) < 60 and x['sample'] not in c['members']:
                c['members'][x['sample']] = x; placed = True; break
        if not placed: clusters.append(dict(mz=x['mz'], rt=x['rt'], members={x['sample']: x}))
    M = np.full((len(clusters), len(samples)), np.nan); I = np.full_like(M, np.nan)
    for i, c in enumerate(clusters):
        for j, sm in enumerate(samples):
            if sm in c['members']:
                M[i, j] = np.log2(c['members'][sm]['area'])
                if c['members'][sm]['id'] is not None: I[i, j] = M[i, j]
    idrows = ~np.all(np.isnan(I), 1)
    res['alignment'] = dict(consensus_features=len(clusters), samples=samples,
                            idfree_missing=float(np.isnan(M).mean()),
                            idfree_complete=int((~np.isnan(M)).all(1).sum()),
                            idfirst_features=int(idrows.sum()),
                            idfirst_missing_within_identified=float(np.isnan(I[idrows]).mean()),
                            idfirst_complete=int((~np.isnan(I)).all(1).sum()))
    comp = (~np.isnan(M)).all(1)
    cc = np.corrcoef(M[comp].T)
    res['alignment']['idfree_pairwise_pearson'] = [float(cc[0, 1]), float(cc[0, 2]), float(cc[1, 2])]
    tot = {k: sum(v[k] for v in res['runs'].values()) for k in ('features', 'with_ms2', 'identified', 'ms2', 'psms')}
    res['total'] = tot
    print(json.dumps(res['alignment'], indent=1)); print(tot)
    json.dump(res, open('realdata/ms1_bottleneck.json', 'w'), indent=1)
    pickle.dump(dict(allf=allf, M=M, I=I, samples=samples), open('realdata/ms1_features.pkl', 'wb'))
