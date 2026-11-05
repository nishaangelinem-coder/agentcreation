"""Fit the single-slot data and synthesise the 4-slot standing-wave subarray."""
import sys, json, numpy as np
sys.path.insert(0, '.')
from siwlab.slots import fit_offset_law, offset_for_conductance, design_standing_wave_array
from siwlab.geometry import SIWParams
from siwlab import siw

prm = SIWParams()
d = json.load(open('results/slot_characterization.json'))
f = np.array(d['f'])
# --- FDTD-calibrated SIW dispersion (through-line) ---
tl = json.load(open('results/throughline.json'))
ft = np.array(tl['f']); bfd = np.array(tl['beta_fdtd'])
i28 = np.argmin(abs(ft - 28e9))
k28 = 2 * np.pi * 28e9 / 3e8 * np.sqrt(prm.er)
a_eff_fd = np.pi / np.sqrt(k28 ** 2 - bfd[i28] ** 2)
lamg = 2 * np.pi / bfd[i28] * 1e3
print(f"FDTD-calibrated: a_eff = {a_eff_fd*1e3:.3f} mm (Cassivi {prm.a_eff:.3f}), lambda_g(28 GHz) = {lamg:.3f} mm")
# --- resonance vs length at the reference offset ---
cases = d['cases']
Lr, fr, gr, xr = [], [], [], []
for key, c in cases.items():
    g = np.array(c['g']); b = np.array(c['b'])
    # resonance: zero crossing of b (from + to -) nearest to the conductance peak
    ip = int(np.argmax(g))
    idx = np.where((b[:-1] > 0) & (b[1:] <= 0))[0]
    if len(idx) == 0 or g[ip] < 0.05:
        print(f"  L={c['L']} xo={c['xo']}: no resonance in band"); continue
    i = idx[np.argmin(np.abs(idx - ip))]
    f0 = f[i] - b[i] * (f[i + 1] - f[i]) / (b[i + 1] - b[i])
    g0 = np.interp(f0, f, g)
    print(f"  L={c['L']} xo={c['xo']}: f_res={f0/1e9:.2f} GHz  g_res={g0:.3f}")
    Lr.append(c['L']); fr.append(f0); gr.append(g0); xr.append(c['xo'])
Lr, fr, gr, xr = map(np.array, (Lr, fr, gr, xr))
# resonant length at 28 GHz from the xo = reference-offset family
ref_xo = 0.2
m = xr == ref_xo
p, q = np.polyfit(Lr[m], 1.0 / fr[m], 1)
L28 = (1 / 28e9 - q) / p
print(f"resonant length at 28 GHz (xo={ref_xo}): {L28:.3f} mm")
# conductance law from cases near resonance at 28 GHz (|f_res-28| < 1.5 GHz)
m2 = np.abs(fr - 28e9) < 1.5e9
K = fit_offset_law(xr[m2], gr[m2], prm.a) if m2.sum() >= 2 else fit_offset_law(xr, gr, prm.a)
print(f"g(x) = K sin^2(pi x/a) with K = {K:.2f}  (from {m2.sum()} near-resonant cases)")
N = 5
slots, g, xo = design_standing_wave_array(N, K, prm.a, L28, lamg, 'uniform', y_short=0.0)
print(f"uniform {N}-slot design: g_n = {g}  ideal offsets (mm) = {np.round(xo, 3)}")
# snap offsets to the 0.1 mm FDTD grid and report the resulting sum of conductances
xo_snap = np.round(xo / 0.1) * 0.1
g_snap = K * np.sin(np.pi * xo_snap / prm.a) ** 2
print(f"grid-snapped offsets = {xo_snap}, g_n = {np.round(g_snap, 3)}, sum = {g_snap.sum():.3f} -> "
      f"|S11| at resonance ~ {abs((1 - g_snap.sum()) / (1 + g_snap.sum())):.3f}")
slots = [(yc, np.sign(x) * xs, L, w) for (yc, x, L, w), xs in zip(slots, xo_snap)]
xo = xo_snap
des = dict(spacing_mm=5.0, first_offset_mm=2.5, g_sum=float(g_snap.sum()), a_eff_fdtd_mm=float(a_eff_fd * 1e3), lamg_mm=float(lamg), L_res_mm=float(L28), K=float(K),
           N=N, g=list(map(float, g)), xo_mm=list(map(float, xo)),
           slots=[(float(yc), float(x), float(L), 0.4) for (yc, x, L, _) in slots],
           fit=dict(L=list(map(float, Lr)), f_res=list(map(float, fr)), g_res=list(map(float, gr)), xo=list(map(float, xr))))
json.dump(des, open('results/design.json', 'w'), indent=1)
print(json.dumps(des['slots']))
