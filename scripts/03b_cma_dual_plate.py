"""CMA step 2: characteristic modes of the common top wall of the two SIW
subarrays (8 radiating slots) without and with the decoupling slits.
usage: python 03b_cma_dual_plate.py <tag> <D_mm> <slit_L_mm or 0> [n_slits]"""
import sys, json, time, numpy as np
sys.path.insert(0, '.')
from siwlab.mesh import plate_mesh
from siwlab.mom import RWGModel, characteristic_modes, characteristic_angle, modal_significance, track_modes

MM = 1e-3
tag = sys.argv[1]; D = float(sys.argv[2]); Ls = float(sys.argv[3])
n_slits = int(sys.argv[4]) if len(sys.argv) > 4 else 5
cal = json.load(open('results/cma_calibration.json'))
des = json.load(open('results/design.json'))
eps_eff = cal['eps_eff']; lamg = des['lamg_mm']
aps = []
for xc in (-D / 2, D / 2):
    for (yc, xo, L, w) in des['slots']:
        aps.append(((xc + xo - w / 2) * MM, (xc + xo + w / 2) * MM, (yc - L / 2) * MM, (yc + L / 2) * MM))
w_slit = cal['slit_w']
for n in range(n_slits if Ls > 0 else 0):
    yc = -des['first_offset_mm'] - n * des['spacing_mm']
    aps.append((-w_slit / 2 * MM, w_slit / 2 * MM, (yc - Ls / 2) * MM, (yc + Ls / 2) * MM))
px = D / 2 + 2.8 + 1.5  # plate half-width: outer via row + 1.5 mm
y0, y1 = -27.0, 2.0
nodes, tris = plate_mesh(-px * MM, px * MM, y0 * MM, y1 * MM, aps, 0.6 * MM)
m = RWGModel(nodes, tris, eps_bg=eps_eff)
print(f"plate {2*px:.1f} x {y1-y0:.1f} mm, {len(tris)} triangles, {m.nE} RWG edges, eps_eff={eps_eff:.3f}")
freqs = np.arange(25e9, 31.01e9, 1e9)
nm = 10
cen = m.rwg['cen']; area = m.rwg['area']
gap = np.abs(cen[:, 0]) < (D / 2 - 2.8 - 0.2) * MM   # region between the inner via rows
lam_tr = np.zeros((len(freqs), nm)); gci = np.zeros((len(freqs), nm)); J_prev = None
Jmaps = {}
t0 = time.time()
for fi, f in enumerate(freqs):
    Z = m.impedance_matrix(f)
    lam, J = characteristic_modes(Z, nm)
    R = 0.5 * (Z.real + Z.real.T)
    if J_prev is not None:
        idx = track_modes(J_prev, J, R)
        lam, J = lam[idx], J[:, idx]
    J_prev = J
    for n in range(nm):
        Jc = m.current_at_centroids(J[:, n])
        Jx = np.abs(Jc[:, 0]); Jm = np.linalg.norm(Jc, axis=1)
        gci[fi, n] = (Jx[gap] * area[gap]).sum() / (Jm * area).sum()
        lam_tr[fi, n] = lam[n]
        if abs(f - 28e9) < 1e6:
            Jmaps[n] = Jc
    print(f"f={f/1e9:.1f} GHz  lam={np.round(lam[:6], 2)}  MS={np.round(modal_significance(lam[:6]), 2)} "
          f" GCI={np.round(gci[fi, :6], 3)}  ({time.time()-t0:.0f}s)", flush=True)
np.savez(f'results/cma_{tag}.npz', freqs=freqs, lam=lam_tr, gci=gci, nodes=nodes, tris=tris, cen=cen, area=area,
         J28=np.array([Jmaps[n] for n in range(nm)]), D=D, Ls=Ls, apertures=np.array(aps), px=px, y0=y0, y1=y1)
print("saved")
