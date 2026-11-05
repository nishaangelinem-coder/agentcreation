"""Dual-port 2 x (1x4) SIW slot array: FDTD S-parameters, far field, ECC.
usage: python 05_array_fdtd.py <tag> <D_mm> <slit_L_mm or 0> [n_slits] [slit_w] [T_ns]"""
import sys, time, json, numpy as np
sys.path.insert(0, '.')
from siwlab.geometry import Layout, SIWParams, build
from siwlab.fdtd import gaussian_modulated, spectrum, aperture_farfield, radiation_intensity, hemisphere_integral
from siwlab.constants import eta0

tag = sys.argv[1]; D = float(sys.argv[2]); Ls = float(sys.argv[3])
n_slits = int(sys.argv[4]) if len(sys.argv) > 4 else 5
w_slit = float(sys.argv[5]) if len(sys.argv) > 5 else 0.3
T = float(sys.argv[6]) * 1e-9 if len(sys.argv) > 6 else 3.6e-9
prm = SIWParams()
des = json.load(open('results/design.json'))
slots = [tuple(s) for s in des['slots']]
lamg = des['lamg_mm']
slits = []
if Ls > 0:
    for n in range(n_slits):
        yc = -des['first_offset_mm'] - n * des['spacing_mm']
        slits.append((0.0, yc, Ls, w_slit))
xc = [-D / 2, D / 2]
reference = (tag == 'reference')
lay = Layout(xc=xc, y_short=None if reference else 0.0, slots=[[], []] if reference else [slots, slots],
             slits=[] if reference else slits, y_fine=(-25.5, 1.0), y_src=-31.0, y_probe1=-27.0, y_min=-36.0)
sim, ports, dft_nodes, apg = build(lay, prm)
src = gaussian_modulated(28e9, 6e9)
nsteps = int(T / sim.dt)
dft_f = np.arange(26e9, 30.01e9, 0.5e9)
t0 = time.time()
out = sim.run(nsteps, ports[0], src, ports[:2], dft_freqs=dft_f, dft_nodes=dft_nodes, progress=4000)
print(f"{nsteps} steps in {time.time()-t0:.0f} s")
f = np.linspace(24e9, 32e9, 321)
V1 = spectrum(out['rec']['P1']['V'], sim.dt, f); I1 = spectrum(out['rec']['P1']['I'], sim.dt, f)
V2 = spectrum(out['rec']['P2']['V'], sim.dt, f)
np.savez(f'results/array_{tag}.npz', f=f, V1=V1, I1=I1, V2=V2, dft_f=dft_f, accEx=out['accEx'], accEy=out['accEy'],
         ex_xy=np.array(apg['ex_xy']), ex_dS=apg['ex_dS'], ey_xy=np.array(apg['ey_xy']), ey_dS=apg['ey_dS'],
         norm=ports[0]['norm'], D=D, Ls=Ls, n_slits=n_slits, w_slit=w_slit, slots=np.array(slots))
if reference:
    print("reference run saved"); sys.exit()
ref = np.load('results/array_reference.npz')
S11 = (V1 - ref['V1']) / ref['V1']; S21 = V2 / ref['V1']
# incident power at the DFT frequencies (from the reference run)
# recompute reference spectra at the dft frequencies
refV = np.interp(dft_f, f, ref['V1'].real) + 1j * np.interp(dft_f, f, ref['V1'].imag)
refI = np.interp(dft_f, f, ref['I1'].real) + 1j * np.interp(dft_f, f, ref['I1'].imag)
Pinc = np.real(refV * np.conj(refI)) / (2 * ports[0]['norm'])
th = np.linspace(0, np.pi / 2, 91); ph = np.linspace(0, 2 * np.pi, 181)
TH, PH = np.meshgrid(th, ph, indexing='ij')
res = dict(f=f, S11=S11, S21=S21, dft_f=dft_f, Pinc=Pinc, theta=th, phi=ph)
G = []; Gr = []; eff = []; Eth_all = []; Eph_all = []
for i, fd in enumerate(dft_f):
    Eth, Eph = aperture_farfield(fd, apg['ex_xy'], apg['ey_xy'], out['accEx'][i], out['accEy'][i], apg['ex_dS'], apg['ey_dS'], TH, PH)
    U = radiation_intensity(fd, Eth, Eph)
    Prad = hemisphere_integral(U, TH, PH)
    j = np.argmin(abs(f - fd))
    Pacc = Pinc[i] * (1 - abs(S11[j]) ** 2 - abs(S21[j]) ** 2)
    G.append(4 * np.pi * U / Pacc); Gr.append(4 * np.pi * U / Pinc[i]); eff.append(Prad / Pacc)
    Eth_all.append(Eth); Eph_all.append(Eph)
    print(f"f={fd/1e9:.1f} GHz |S11|={20*np.log10(abs(S11[j])):.1f} dB |S21|={20*np.log10(abs(S21[j])):.1f} dB "
          f"Gmax={10*np.log10(G[-1].max()):.2f} dBi  realized={10*np.log10(Gr[-1].max()):.2f} dBi  eff={eff[-1]*100:.1f}%  "
          f"peak at th={np.degrees(TH.ravel()[G[-1].argmax()]):.0f} ph={np.degrees(PH.ravel()[G[-1].argmax()]):.0f}")
np.savez(f'results/array_{tag}_post.npz', **res, G=np.array(G), Gr=np.array(Gr), eff=np.array(eff),
         Eth=np.array(Eth_all), Eph=np.array(Eph_all))
