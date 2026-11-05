"""Aperture-system analysis of the dual sub-array (magnetic-current MoM):
all radiating slots (+ optional decoupling slits) are apertures in the infinite
top wall, air above, effective dielectric half-space below.
 (a) external coupling between the two sub-arrays through the aperture admittance
     matrix, with the parasitic slits eliminated by a Schur complement, versus
     slit length;
 (b) characteristic modes (B M = lambda G M) of the aperture system versus
     frequency for the chosen slit length.
usage: python 04_aperture_cma_array.py [eps_eff]"""
import sys, json, time, numpy as np
sys.path.insert(0, '.')
from siwlab.mesh import plate_mesh, rwg_edges
from siwlab.mom import RWGModel, aperture_admittance_matrix, characteristic_modes, modal_significance

MM = 1e-3
des = json.load(open('results/design.json'))
cal = json.load(open('results/cma_calibration.json'))
eps_eff = float(sys.argv[1]) if len(sys.argv) > 1 else cal['eps_eff']
D = float(sys.argv[2]) if len(sys.argv) > 2 else 7.0; w_slit = 0.3; hmax = 0.2 * MM


def aperture_system(Ld, n_slits=5):
    nodes_all, tris_all, group = [], [], []
    n0 = 0
    aps = []
    for gi, xc in enumerate((-D / 2, D / 2)):
        for (yc, xo, L, w) in des['slots']:
            aps.append((xc + xo, yc, L, w, gi))
    if Ld > 0:
        if n_slits == 1:
            yc = -des['first_offset_mm'] - (des['N'] - 1) * des['spacing_mm'] / 2
            aps.append((0.0, yc, Ld, w_slit, 2))
        else:
            for n in range(n_slits):
                aps.append((0.0, -des['first_offset_mm'] - n * des['spacing_mm'], Ld, w_slit, 2))
    for (xc, yc, L, w, g) in aps:
        nodes, tris = plate_mesh((xc - w / 2) * MM, (xc + w / 2) * MM, (yc - L / 2) * MM, (yc + L / 2) * MM, [], hmax)
        nodes_all.append(nodes); tris_all.append(tris + n0); group.append(np.full(len(tris), g)); n0 += len(nodes)
    nodes = np.vstack(nodes_all); tris = np.vstack(tris_all); group = np.concatenate(group)
    m = RWGModel(nodes, tris, chunk=256)
    egroup = group[m.rwg['tp']]          # group of each RWG edge
    return m, egroup


def coupling_metric(Y, egroup):
    """Normalised external mutual admittance between the two sub-arrays after
    eliminating the parasitic slits (Schur complement). The slot excitation is
    the designed uniform, in-phase distribution (unit magnetic current coefficient
    on every RWG edge of a slot is not meaningful, so the sub-array 'port' vector is
    the dominant aperture mode of an isolated slot, i.e. the first characteristic
    mode of a single slot, replicated on every slot)."""
    s = egroup < 2; l = egroup == 2
    Yss = Y[np.ix_(s, s)]
    if l.any():
        Yll = Y[np.ix_(l, l)]; Ysl = Y[np.ix_(s, l)]
        Yeff = Yss - Ysl @ np.linalg.solve(Yll, Ysl.T)
    else:
        Yeff = Yss
    g = egroup[s]
    # port vectors: dominant mode of each sub-array block (self admittance eigenvector)
    v = []
    for gi in (0, 1):
        k = g == gi
        Yk = Yeff[np.ix_(k, k)]
        lam, M = characteristic_modes(Yk, 1)
        vk = np.zeros(len(g), complex); vk[k] = M[:, 0]
        v.append(vk)
    y12 = v[1].conj() @ Yeff @ v[0]; y11 = v[0].conj() @ Yeff @ v[0]
    return abs(y12) / abs(y11), Yeff


t0 = time.time()
out = dict(eps_eff=eps_eff, sweep={})
f0 = 28e9
m, eg = aperture_system(0.0)
Y0 = aperture_admittance_matrix(m, f0, (1.0, eps_eff))
c0_, _ = coupling_metric(Y0, eg)
print(f"no slits: {m.nE} edges, |y12/y11| = {c0_:.4f}  ({time.time()-t0:.0f}s)")
out['sweep']['0'] = c0_
for Ld in [3.2, 3.6, 4.0, 4.4, 4.8, 5.2, 5.6]:
    m, eg = aperture_system(Ld)
    Y = aperture_admittance_matrix(m, f0, (1.0, eps_eff))
    c, _ = coupling_metric(Y, eg)
    out['sweep'][f"{Ld:.1f}"] = c
    print(f"Ld = {Ld} mm: {m.nE} edges, |y12/y11| = {c:.4f}  ({20*np.log10(c/c0_):+.1f} dB)  ({time.time()-t0:.0f}s)", flush=True)
m, eg = aperture_system(24.0, 1)
Y = aperture_admittance_matrix(m, f0, (1.0, eps_eff)); c, _ = coupling_metric(Y, eg)
out['sweep']['cont'] = c
print(f"continuous slit: |y12/y11| = {c:.4f} ({20*np.log10(c/c0_):+.1f} dB)")
json.dump(out, open(f'results/aperture_coupling_D{D:.1f}_eps{eps_eff:.2f}.json', 'w'), indent=1)
# characteristic modes of the aperture system versus frequency for the best slit length
best = min((k for k in out['sweep'] if k not in ('0', 'cont')), key=lambda k: out['sweep'][k])
Ld = float(best)
m, eg = aperture_system(Ld)
freqs = np.arange(25e9, 31.01e9, 1e9)
nm = 8
lam_tr = np.zeros((len(freqs), nm)); slitfrac = np.zeros((len(freqs), nm)); even = np.zeros((len(freqs), nm))
area = m.rwg['area']; tgroup = eg[np.arange(m.nE)]
for fi, f in enumerate(freqs):
    Y = aperture_admittance_matrix(m, f, (1.0, eps_eff))
    lam, M = characteristic_modes(Y, nm)
    for n in range(nm):
        Jc = m.current_at_centroids(M[:, n]); Jm = np.linalg.norm(Jc, axis=1)
        tg = np.zeros(m.nt, int)
        for e in range(m.nE):
            tg[m.rwg['tp'][e]] = eg[e]; tg[m.rwg['tm'][e]] = eg[e]
        slitfrac[fi, n] = (Jm[tg == 2] * area[tg == 2]).sum() / (Jm * area).sum()
        p1 = (Jm[tg == 0] * area[tg == 0]).sum(); p2 = (Jm[tg == 1] * area[tg == 1]).sum()
        even[fi, n] = min(p1, p2) / max(p1, p2, 1e-30)
        lam_tr[fi, n] = lam[n]
    print(f"f={f/1e9:.0f} GHz MS={np.round(modal_significance(lam), 2)} slitfrac={np.round(slitfrac[fi], 2)} ({time.time()-t0:.0f}s)", flush=True)
np.savez(f'results/aperture_cma_D{D:.1f}_eps{eps_eff:.2f}.npz', freqs=freqs, lam=lam_tr, slitfrac=slitfrac, even=even, Ld=Ld, eps_eff=eps_eff)
print('saved')
