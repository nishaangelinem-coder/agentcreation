"""Extract the numbers quoted in the manuscript into paper/numbers.tex (LaTeX macros)
and paper/table_sweep.tex, so that the text is consistent with the data."""
import json, glob, os, numpy as np
import sys
sys.path.insert(0, '.')
from siwlab.geometry import SIWParams
from siwlab import siw

prm = SIWParams()
des = json.load(open('results/design.json'))
cal = json.load(open('results/cma_calibration.json'))
tl = json.load(open('results/throughline.json'))
M = {}


def db(x):
    return 20 * np.log10(np.maximum(np.abs(x), 1e-12))


def band(f, s11):
    """-10 dB band around 28 GHz (GHz)."""
    ok = db(s11) < -10
    i = np.argmin(abs(f - 28e9))
    if not ok[i]:
        i = np.argmin(db(s11))
    lo = i
    while lo > 0 and ok[lo - 1]:
        lo -= 1
    hi = i
    while hi < len(f) - 1 and ok[hi + 1]:
        hi += 1
    return f[lo] / 1e9, f[hi] / 1e9


f = np.array(tl['f']); i28 = np.argmin(abs(f - 28e9))
M['betaErr'] = f"{abs(tl['beta_fdtd'][i28] / tl['beta_an'][i28] - 1) * 100:.1f}"
M['aeffFD'] = f"{des['a_eff_fdtd_mm']:.2f}"
M['aeffDiff'] = f"{(1 - des['a_eff_fdtd_mm'] / prm.a_eff) * 100:.1f}"
M['lamg'] = f"{des['lamg_mm']:.2f}"
M['STLtl'] = f"{tl['S21_dB'][i28]:.2f}"
fit = des['fit']
M['Kfit'] = f"{des['K']:.1f}"
M['Lres'] = f"{des['L_res_mm']:.1f}"
Lr = np.array(fit['L']); fr_ = np.array(fit['f_res']); xr_ = np.array(fit['xo'])
m_ = xr_ == 0.2; p_, q_ = np.polyfit(Lr[m_], 1.0 / fr_[m_], 1); M['LresIso'] = f"{(1 / 28e9 - q_) / p_:.2f}"
g02 = [g for L, x, g, fr in zip(fit['L'], fit['xo'], fit['g_res'], fit['f_res']) if x == 0.2 and abs(fr - 28e9) < 1e9]
M['gslot'] = f"{np.mean(g02):.2f}"
M['gsum'] = f"{des['g_sum']:.2f}"
M['SSum'] = f"{abs((1 - des['g_sum']) / (1 + des['g_sum'])):.2f}"
M['fAir'] = f"{cal['f_air'] / 1e9:.1f}"
M['epsEff'] = f"{cal['eps_eff']:.2f}"
M['Ld'] = f"{cal['slit_L_28GHz']:.2f}"
Ld = 2 * round(cal['slit_L_28GHz'] / 2, 1)
# arrays
import os
REF_TAG = os.environ.get('REF_TAG', 'D94_base'); DEC_TAG = os.environ.get('DEC_TAG', 'D94_o40')
ref = np.load(f'results/array_{REF_TAG}_post.npz'); dec = np.load(f'results/array_{DEC_TAG}_post.npz'); decraw = np.load(f'results/array_{DEC_TAG}.npz')
d70 = np.load('results/array_nodecouple_post.npz')
fa = ref['f']; j28 = np.argmin(abs(fa - 28e9))
lo_r, hi_r = band(fa, ref['S11']); lo_d, hi_d = band(fa, dec['S11'])
M['BWref'] = f"{lo_r:.2f}--{hi_r:.2f}"; M['BWrefRel'] = f"{(hi_r - lo_r) / 28 * 100:.1f}"
M['BWdec'] = f"{lo_d:.2f}--{hi_d:.2f}"; M['fcDec'] = f"{0.5 * (lo_d + hi_d):.1f}"; M['BWdecRel'] = f"{(hi_d - lo_d) / 28 * 100:.1f}"
inb_r = (fa >= 27.5e9) & (fa <= 28.5e9); inb_d = inb_r   # fixed evaluation band 27.5-28.5 GHz (n261)
def band6(f, s11):
    ok = db(s11) < -6; i = np.argmin(db(s11)); lo = i; hi = i
    while lo > 0 and ok[lo - 1]: lo -= 1
    while hi < len(f) - 1 and ok[hi + 1]: hi += 1
    return f[lo] / 1e9, f[hi] / 1e9
l6r, h6r = band6(fa, ref['S11']); l6d, h6d = band6(fa, dec['S11'])
M['BWsixRef'] = f"{l6r:.2f}--{h6r:.2f}"; M['BWsixDec'] = f"{l6d:.2f}--{h6d:.2f}"; M['BWsixDecRel'] = f"{(h6d - l6d) / 28 * 100:.1f}"
M['SElevenDec'] = f"{db(dec['S11']).min():.1f}"; M['SElevenRef'] = f"{db(ref['S11']).min():.1f}"
M['IsoRefTwentyEight'] = f"{db(ref['S21'])[j28]:.1f}"; M['IsoDecTwentyEight'] = f"{db(dec['S21'])[j28]:.1f}"
M['IsoRefMin'] = f"{db(ref['S21'])[inb_r].max():.1f}"; M['IsoDecMin'] = f"{db(dec['S21'])[inb_d].max():.1f}"
M['IsoGain'] = f"{db(ref['S21'])[j28] - db(dec['S21'])[j28]:.1f}"
fd = ref['dft_f']; k28 = np.argmin(abs(fd - 28e9))
M['IsoSeventy'] = f"{db(d70['S21'])[j28]:.1f}"; lo7, hi7 = band(fa, d70['S11']); inb7 = inb_r
M['IsoSeventyMin'] = f"{db(d70['S21'])[inb7].max():.1f}"; M['BWseventy'] = f"{lo7:.2f}--{hi7:.2f}"; M['GainSeventy'] = f"{10 * np.log10(d70['Gr'][k28].max()):.1f}"
M['Dfin'] = f"{float(decraw['D']):.1f}"; M['Lduse'] = f"{float(decraw['Ls']):.1f}"; M['Xoff'] = f"{float(decraw['xoff']):.1f}"; M['IsoDecMax'] = f"{db(dec['S21'])[inb_d].min():.1f}"
M['GainRef'] = f"{10 * np.log10(ref['Gr'][k28].max()):.1f}"; M['GainDec'] = f"{10 * np.log10(dec['Gr'][k28].max()):.1f}"
M['EffDec'] = f"{dec['eff'][k28] * 100:.0f}"
M['GainDelta'] = f"{10 * np.log10(dec['Gr'][k28].max()) - 10 * np.log10(ref['Gr'][k28].max()):.1f}"
th = dec['theta']; ph = dec['phi']
G = dec['G'][k28]
TH, PH = np.meshgrid(th, ph, indexing='ij')
def hpbw(G, th, iph1, iph2):
    cut = np.r_[G[::-1, iph2], G[:, iph1]]; ang = np.r_[-np.degrees(th[::-1]), np.degrees(th)]
    ok = cut >= cut.max() / 2
    return ang[ok].max() - ang[ok].min()
M['HPBWe'] = f"{hpbw(G, th, np.argmin(abs(ph - np.pi/2)), np.argmin(abs(ph - 3*np.pi/2))):.0f}"
M['HPBWh'] = f"{hpbw(G, th, 0, np.argmin(abs(ph - np.pi))):.0f}"
def tilt(p):
    G = p['G'][k28]; i = np.unravel_index(G.argmax(), G.shape)
    t = np.degrees(th[i[0]]); return t * np.sign(np.cos(ph[i[1]]))
M['TiltRef'] = f"{tilt(ref):+.0f}"; M['TiltDec'] = f"{tilt(dec):+.0f}"
def ecc(p):
    Eth, Eph = p['Eth'][k28], p['Eph'][k28]
    idx = np.argmin(abs(((np.pi - ph)[None, :] - ph[:, None] + np.pi) % (2 * np.pi) - np.pi), axis=1)
    E2th = Eth[:, idx]; E2ph = -Eph[:, idx]; w = np.sin(TH)
    num = abs(np.sum((Eth * np.conj(E2th) + Eph * np.conj(E2ph)) * w)) ** 2
    den = np.sum((abs(Eth) ** 2 + abs(Eph) ** 2) * w) * np.sum((abs(E2th) ** 2 + abs(E2ph) ** 2) * w)
    return num / den
M['ECCref'] = f"{ecc(ref):.3f}"; M['ECCdec'] = f"{ecc(dec):.3f}"
ac = siw.conductor_attenuation(28e9, des['a_eff_fdtd_mm'] * 1e-3, prm.h * 1e-3, prm.er) * 0.025 * 20 * np.log10(np.e)
M['CondLoss'] = f"{ac:.2f}"
# CMA plate
try:
    r = np.load('results/cma_ref.npz'); s = np.load('results/cma_slits.npz')
    from siwlab.mom import modal_significance
    def norm_gci(data):
        cen = data['cen']; area = data['area']; gap = np.abs(cen[:, 0]) < (float(data['D']) / 2 - 2.8 - 0.2) * 1e-3
        return data['gci'] * area.sum() / area[gap].sum()
    gr_, gs_ = norm_gci(r), norm_gci(s)
    ir = np.argmin(abs(r['freqs'] - 28e9)); n = int(np.argmax(gr_[ir]))
    M['MSref'] = f"{modal_significance(r['lam'][ir, n]):.2f}"; M['GCIref'] = f"{gr_[ir, n]:.2f}"
    isd = np.argmin(abs(s['freqs'] - 28e9)); sig = modal_significance(s['lam'][isd]) > 0.5
    M['GCIslit'] = f"{gs_[isd][sig].max() if sig.any() else gs_[isd].max():.2f}"
    M['platesize'] = f"{2 * float(r['px']):.1f}~mm$\\times$~{float(r['y1']) - float(r['y0']):.0f}~mm"
except Exception as e:
    print('CMA plate data missing:', e)
    M.update(MSref='--', GCIref='--', GCIslit='--', platesize='--')
# sweep table over all array variants
def describe(p):
    D_ = float(p['D']); L = float(p['Ls']) if 'Ls' in p.files else 0.0
    cav = int(p['cav']) if 'cav' in p.files else 0; xo = float(p['xoff']) if 'xoff' in p.files else 0.0
    n = int(p['n_slits']) if 'n_slits' in p.files else 5
    if L == 0:
        return f"{D_:.1f} & none & --"
    kind = ('cont. slit' if n == 1 else ('cav. slot' if cav else 'slit'))
    return f"{D_:.1f} & {kind}{' ($x_d$=%.1f)' % xo if xo else (' (centred)' if cav else '')} & {L:.1f}"
rows = []
for fn in sorted(glob.glob('results/array_*_post.npz')):
    tag = fn[len('results/array_'):-len('_post.npz')]
    if tag in ('nodecouple_L40',):
        continue
    p = np.load(fn); raw = np.load(fn.replace('_post.npz', '.npz'))
    lo, hi = band(p['f'], p['S11']); inb = (p['f'] >= lo * 1e9) & (p['f'] <= hi * 1e9)
    inb = (p['f'] >= 27.5e9) & (p['f'] <= 28.5e9)
    rows.append((float(raw['D']), float(raw['Ls']) if 'Ls' in raw.files else 0.0, tag, describe(raw), db(p['S11'])[j28], db(p['S21'])[j28],
                 db(p['S21'])[inb].max(), 10 * np.log10(p['Gr'][k28].max())))
rows.sort(key=lambda r: (r[0], r[1], r[2]))
with open('paper/table_sweep.tex', 'w') as fh:
    fh.write("\\begin{tabular}{@{}llccccc@{}}\\toprule\n$D$ (mm) & structure & $L_d$ (mm) & $|S_{11}|$ & $|S_{21}|$ & $|S_{21}|_{\\max}$ & $G_r$ (dBi)\\\\\\midrule\n")
    for r in rows:
        fh.write(f"{r[3]} & {r[4]:.1f} & {r[5]:.1f} & {r[6]:.1f} & {r[7]:.1f}\\\\\n")
    fh.write("\\bottomrule\\end{tabular}\n")
M['SweepSentence'] = "see Table~\\ref{tab:sweep}"
with open('paper/numbers.tex', 'w') as fh:
    for k, v in M.items():
        fh.write(f"\\newcommand{{\\{k}}}{{{v}}}\n")
print(json.dumps(M, indent=1))
