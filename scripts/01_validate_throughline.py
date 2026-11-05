"""Validation 1: plain SIW through-line. Checks CPML/port behaviour, extracts the
TE10 phase constant from S21 and compares with the analytical SIW model."""
import sys, time, json, numpy as np
sys.path.insert(0, '.')
from siwlab.geometry import Layout, SIWParams, build
from siwlab.fdtd import gaussian_modulated, spectrum
from siwlab import siw

prm = SIWParams()
lay = Layout(xc=[0.0], y_short=None, slots=[[]], y_fine=(-6.0, 6.0))
sim, ports, dft_nodes, apg = build(lay, prm)
src = gaussian_modulated(28e9, 6e9)
T = 2.2e-9
nsteps = int(T / sim.dt)
t0 = time.time()
out = sim.run(nsteps, ports[0], src, ports, progress=2000)
print(f"{nsteps} steps in {time.time()-t0:.1f} s  ({(time.time()-t0)/nsteps*1e3:.1f} ms/step)")
f = np.linspace(22e9, 34e9, 241)
V1 = spectrum(out['rec']['P1']['V'], sim.dt, f); I1 = spectrum(out['rec']['P1']['I'], sim.dt, f)
V2 = spectrum(out['rec']['P1out']['V'], sim.dt, f)
np.savez('results/throughline.npz', f=f, V1=V1, I1=I1, V2=V2, t=out['t'],
         Vt1=out['rec']['P1']['V'], It1=out['rec']['P1']['I'], Vt2=out['rec']['P1out']['V'],
         norm=ports[0]['norm'], dt=sim.dt)
# S21 magnitude (incident = V1 in a matched line)
S21 = V2 / V1
Ly = (lay.y_probe2 - lay.y_probe1) * 1e-3
phase = np.unwrap(np.angle(S21))
beta_fd = -phase / Ly
# remove 2*pi ambiguity: anchor to analytical at 28 GHz
b_an = siw.beta(f, prm.a_eff * 1e-3, prm.er)
k = np.round((b_an[120] - beta_fd[120]) * Ly / (2 * np.pi))
beta_fd += 2 * np.pi * k / Ly
print("a_eff = %.3f mm, fc = %.2f GHz" % (prm.a_eff, siw.cutoff_frequency(prm.a_eff*1e-3, prm.er)/1e9))
for fi in (24e9, 26e9, 28e9, 30e9, 32e9):
    i = np.argmin(abs(f - fi))
    print(f"f={fi/1e9:.0f} GHz |S21|={20*np.log10(abs(S21[i])):.3f} dB  beta_fdtd={beta_fd[i]:.1f}  beta_an={b_an[i]:.1f} rad/m  ({100*(beta_fd[i]/b_an[i]-1):+.2f}%)")
ad = siw.dielectric_attenuation(f, prm.a_eff*1e-3, prm.er, prm.tand)
print("analytical dielectric loss over %.1f mm at 28 GHz: %.3f dB" % (Ly*1e3, 20*np.log10(np.e)*ad[120]*Ly))
json.dump(dict(f=list(f), S21_dB=list(20*np.log10(abs(S21))), beta_fdtd=list(beta_fd), beta_an=list(b_an)),
          open('results/throughline.json', 'w'))
