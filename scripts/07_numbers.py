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
Ld = 2 * round(cal['slit_L_28GHz'] / 2, 1)   # implemented on the 0.2 mm length grid
M['Lduse'] = f"{Ld:.1f}"
# arrays
ref = np.load('results/array_nodecouple_post.npz'); dec = np.load('results/array_decoupled_post.npz')
fa = ref['f']; j28 = np.argmin(abs(fa - 28e9))
lo_r, hi_r = band(fa, ref['S11']); lo_d, hi_d = band(fa, dec['S11'])
M['BWref'] = f"{lo_r:.2f}--{hi_r:.2f}"; M['BWrefRel'] = f"{(hi_r - lo_r) / 28 * 100:.1f}"
M['BWdec'] = f"{lo_d:.2f}--{hi_d:.2f}"; M['fcDec'] = f"{0.5 * (lo_d + hi_d):.1f}"; M['BWdecRel'] = f"{(hi_d - lo_d) / 28 * 100:.1f}"
inb_r = (fa >= lo_r * 1e9) & (fa <= hi_r * 1e9); inb_d = (fa >= lo_d * 1e9) & (fa <= hi_d * 1e9)
M['IsoRefTwentyEight'] = f"{db(ref['S21'])[j28]:.1f}"; M['IsoDecTwentyEight'] = f"{db(dec['S21'])[j28]:.1f}"
M['IsoRefMin'] = f"{db(ref['S21'])[inb_r].max():.1f}"; M['IsoDecMin'] = f"{db(dec['S21'])[inb_d].max():.1f}"
M['IsoGain'] = f"{db(ref['S21'])[j28] - db(dec['S21'])[j28]:.1f}"
fd = ref['dft_f']; k28 = np.argmin(abs(fd - 28e9))
M['GainRef'] = f"{10 * np.log10(ref['Gr'][k28].max()):.1f}"; M['GainDec'] = f"{10 * np.log10(dec['Gr'][k28].max()):.1f}"
M['EffDec'] = f"{dec['eff'][k28] * 100:.0f}"
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
# slit sweep table
rows = []
for tag in ['nodecouple'] + sorted(glob.glob('results/array_decoupled*_post.npz')):
    p = np.load(tag if tag.endswith('.npz') else f'results/array_{tag}_post.npz')
    L = float(p['Ls']) if 'Ls' in p.files else 0.0
    if tag == 'nodecouple':
        L = 0.0
    lo, hi = band(p['f'], p['S11']); inb = (p['f'] >= lo * 1e9) & (p['f'] <= hi * 1e9)
    rows.append((L, db(p['S11'])[j28], db(p['S21'])[j28], db(p['S21'])[inb].max(), lo, hi, 10 * np.log10(p['Gr'][k28].max())))
rows = sorted(set(rows))
with open('paper/table_sweep.tex', 'w') as fh:
    fh.write("\\begin{tabular}{@{}lccccc@{}}\\toprule\n$L_d$ (mm) & $|S_{11}|$ @28 GHz (dB) & $|S_{21}|$ @28 GHz (dB) & max in-band $|S_{21}|$ (dB) & $-10$ dB band (GHz) & realised gain (dBi)\\\\\\midrule\n")
    for (L, s11, s21, s21m, lo, hi, g) in rows:
        fh.write(f"{'none' if L == 0 else f'{L:.1f}'} & {s11:.1f} & {s21:.1f} & {s21m:.1f} & {lo:.2f}--{hi:.2f} & {g:.1f}\\\\\n")
    fh.write("\\bottomrule\\end{tabular}\n")
sw = [r for r in rows if r[0] > 0]
best = min(sw, key=lambda r: r[2])
M['SweepSentence'] = (f"the isolation at 28~GHz is " + ", ".join(f"{r[2]:.1f}~dB for $L_d = {r[0]:.1f}$~mm" for r in sw) +
                      f", the best value being obtained for $L_d = {best[0]:.1f}$~mm" + (" (the CMA value)" if abs(best[0] - Ld) < 0.05 else f" against the CMA prediction of {Ld:.1f}~mm"))
with open('paper/numbers.tex', 'w') as fh:
    for k, v in M.items():
        fh.write(f"\\newcommand{{\\{k}}}{{{v}}}\n")
print(json.dumps(M, indent=1))
