"""CMA step 1 (aperture formulation): characteristic modes of a single aperture
in an infinite PEC plane with air above and an effective dielectric half-space
below (magnetic-current MoM, B M = lambda G M).
 (i)  the effective permittivity of the lower half-space is calibrated so that the
      FDTD-resonant radiating slot (L_s, xo = 0.2 mm) has its aperture-mode
      resonance (CA = 180 deg) at 28 GHz;
 (ii) with that medium, the design curve 'slit length -> resonance' of the
      0.3 mm wide decoupling slit is computed and L_d(28 GHz) is read off."""
import sys, json, time, numpy as np
sys.path.insert(0, '.')
from scipy.optimize import brentq
from siwlab.mesh import plate_mesh
from siwlab.mom import RWGModel, aperture_admittance_matrix, characteristic_modes, characteristic_angle

MM = 1e-3
L_slot = float(sys.argv[1]) if len(sys.argv) > 1 else 4.0
w_slot, w_slit = 0.4, 0.3


def aperture_model(L, w, h=0.1):
    nodes, tris = plate_mesh(-w / 2 * MM, w / 2 * MM, -L / 2 * MM, L / 2 * MM, [], h * MM)
    return RWGModel(nodes, tris)


def first_mode_lambda(m, f, eps_below):
    Y = aperture_admittance_matrix(m, f, (1.0, eps_below))
    lam, M = characteristic_modes(Y, 3)      # same generalised eigenproblem with (G, B)
    return lam[0]


def resonance(m, eps_below, f_lo=14e9, f_hi=40e9):
    """First +/- sign change of the dominant eigenvalue on a coarse scan, refined by brentq."""
    g = lambda f: first_mode_lambda(m, f, eps_below)
    fs = np.arange(f_lo, f_hi + 1, 1e9)
    v = [g(f) for f in fs]
    for i in range(len(fs) - 1):
        if v[i] * v[i + 1] < 0:
            return brentq(g, fs[i], fs[i + 1], xtol=1e7)
    return np.nan


t0 = time.time()
out = {}
ms = aperture_model(L_slot, w_slot)
print(f"slot aperture mesh: {ms.nt} triangles, {ms.nE} edges")
f_air = resonance(ms, 1.0)
print(f"slot L={L_slot} mm in air/air: aperture-mode resonance {f_air/1e9:.2f} GHz ({time.time()-t0:.0f}s)")
eps_eff = brentq(lambda e: resonance(ms, e) - 28e9, 1.0, 3.0, xtol=1e-3)
print(f"==> effective permittivity of the lower half-space: eps_eff = {eps_eff:.3f}  ({time.time()-t0:.0f}s)")
out.update(slot_L=L_slot, f_air=f_air, eps_eff=eps_eff, slit_w=w_slit)
# CA(f) curves for the figure
fr = np.arange(22e9, 36.1e9, 1e9)
out['slot_ca_free'] = dict(f=list(fr), ca=[float(characteristic_angle(first_mode_lambda(ms, f, 1.0))) for f in fr])
out['slot_ca_eff'] = dict(f=list(fr), ca=[float(characteristic_angle(first_mode_lambda(ms, f, eps_eff))) for f in fr])
# slit design curve
curve = {}
for L in [3.4, 3.7, 4.0, 4.3, 4.6]:
    md = aperture_model(L, w_slit)
    frs = resonance(md, eps_eff)
    curve[f"{L:.1f}"] = dict(f_res=frs, ca=[float(characteristic_angle(first_mode_lambda(md, f, eps_eff))) for f in fr])
    print(f"slit L={L} mm: aperture-mode resonance {frs/1e9:.2f} GHz ({time.time()-t0:.0f}s)")
Lv = np.array([float(k) for k in curve]); fv = np.array([curve[k]['f_res'] for k in curve])
p, q = np.polyfit(Lv, 1.0 / fv, 1)
L28 = float((1 / 28e9 - q) / p)
# refine by root finding on the continuous model
L28 = brentq(lambda L: resonance(aperture_model(L, w_slit), eps_eff) - 28e9, Lv.min(), Lv.max(), xtol=1e-3)
print(f"==> slit length for a 28 GHz aperture-mode resonance: L_d = {L28:.3f} mm")
out['slit_curve'] = curve; out['slit_f'] = list(fr); out['slit_L_28GHz'] = L28
json.dump(out, open('results/cma_calibration.json', 'w'), indent=1)
