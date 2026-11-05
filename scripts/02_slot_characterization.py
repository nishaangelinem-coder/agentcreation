"""Single-slot characterisation (Elliott procedure with full-wave data):
one longitudinal slot in a matched SIW through-line; shunt admittance
y = 2(1 - S21)/S21 extracted with the slot-free through-line as reference."""
import sys, time, json, numpy as np
sys.path.insert(0, '.')
from siwlab.geometry import Layout, SIWParams, build
from siwlab.fdtd import gaussian_modulated, spectrum

cases = [(3.6, 0.9), (3.9, 0.9), (4.2, 0.9), (4.5, 0.9), (4.2, 0.5), (4.2, 0.7), (4.2, 1.1), (4.2, 1.3)]
if len(sys.argv) > 1:
    cases = [tuple(map(float, c.split(','))) for c in sys.argv[1:]]
prm = SIWParams()
ref = np.load('results/throughline.npz')
f = ref['f']
import os
res = json.load(open('results/slot_characterization.json'))['cases'] if os.path.exists('results/slot_characterization.json') else {}
for (L, xo) in cases:
    lay = Layout(xc=[0.0], y_short=None, slots=[[(0.0, xo, L, 0.4)]], y_fine=(-6.0, 6.0))
    sim, ports, dft_nodes, apg = build(lay, prm)
    src = gaussian_modulated(28e9, 6e9)
    nsteps = int(2.6e-9 / sim.dt)
    t0 = time.time()
    out = sim.run(nsteps, ports[0], src, ports, progress=5000)
    V1 = spectrum(out['rec']['P1']['V'], sim.dt, f)
    V2 = spectrum(out['rec']['P1out']['V'], sim.dt, f)
    S21 = V2 / ref['V2']
    S11 = (V1 - ref['V1']) / ref['V1']
    y = 2 * (1 - S21) / S21
    i28 = np.argmin(abs(f - 28e9))
    ires = np.argmin(abs(y.imag))
    print(f"L={L} xo={xo}: {time.time()-t0:.0f} s; at 28 GHz y = {y[i28]:.3f}; b=0 at {f[ires]/1e9:.2f} GHz, g_res={y[ires].real:.3f}")
    res[f"{L},{xo}"] = dict(L=L, xo=xo, S21=[(c.real, c.imag) for c in S21], S11=[(c.real, c.imag) for c in S11],
                            g=list(y.real), b=list(y.imag))
    json.dump(dict(f=list(f), cases=res), open('results/slot_characterization.json', 'w'))
