"""Publication figures for the manuscript (all data from results/)."""
import sys, json, os, numpy as np
sys.path.insert(0, '.')
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, Circle
from siwlab.geometry import SIWParams
from siwlab import siw

plt.rcParams.update({'font.size': 8, 'font.family': 'serif', 'axes.linewidth': 0.6, 'lines.linewidth': 1.1,
                     'legend.fontsize': 7, 'figure.dpi': 200, 'savefig.bbox': 'tight'})
OUT = 'paper/figures'
os.makedirs(OUT, exist_ok=True)
prm = SIWParams()
des = json.load(open('results/design.json'))
D = 7.0


def db(x):
    return 20 * np.log10(np.maximum(np.abs(x), 1e-12))


def fig_geometry(Ls):
    fig, (ax, ax2) = plt.subplots(2, 1, figsize=(3.6, 6.2), gridspec_kw=dict(height_ratios=[4.2, 1]))
    # top view
    ax.add_patch(Rectangle((-9, -36), 18, 39, fc='#d9d9d9', ec='none'))   # top metal
    for xc in (-D / 2, D / 2):
        for xr in (xc - prm.a / 2, xc + prm.a / 2):
            for yc in np.arange(0, -36, -prm.p):
                ax.add_patch(Circle((xr, yc), prm.d / 2, fc='k', ec='none'))
        ax.plot([xc - prm.a / 2, xc + prm.a / 2], [0, 0], 'k-', lw=2.5)
        for (yc, xo, L, w) in des['slots']:
            ax.add_patch(Rectangle((xc + xo - w / 2, yc - L / 2), w, L, fc='white', ec='k', lw=0.6))
        ax.annotate('', xy=(xc, -33.5), xytext=(xc, -36.5), arrowprops=dict(arrowstyle='->', lw=1))
        ax.text(xc, -37.3, 'P1' if xc < 0 else 'P2', ha='center', va='top', fontsize=8)
    for n in range(des['N']):
        yc = -des['first_offset_mm'] - n * des['spacing_mm']
        ax.add_patch(Rectangle((-0.15, yc - Ls / 2), 0.3, Ls, fc='white', ec='r', lw=0.8))
    ax.set_xlim(-9.5, 9.5); ax.set_ylim(-39, 4.5); ax.set_aspect('equal')
    ax.set_xlabel('x (mm)'); ax.set_ylabel('y (mm)')
    ax.text(0, 1.2, 'shorting via wall', ha='center', fontsize=7)
    ax.annotate('decoupling slits\n$L_d \\times$ 0.3 mm', xy=(0.2, -22.5), xytext=(2.4, -31.5), fontsize=7, color='r',
                arrowprops=dict(arrowstyle='->', color='r', lw=0.7))
    ax.annotate('$L_s \\times w$', xy=(-D/2 + 0.2 + 0.25, -2.5), xytext=(-9.2, -1.0), fontsize=7,
                arrowprops=dict(arrowstyle='->', lw=0.6))
    ax.annotate('$\\pm x_0$', xy=(-D/2 - 0.2, -7.5), xytext=(-9.2, -10.5), fontsize=7,
                arrowprops=dict(arrowstyle='->', lw=0.6))
    ax.annotate('', xy=(D/2 + 3.6, -2.5), xytext=(D/2 + 3.6, -7.5), arrowprops=dict(arrowstyle='<->', lw=0.6))
    ax.text(D/2 + 3.9, -5.0, '$d_s$', fontsize=7, va='center')
    ax.annotate('', xy=(-D/2, -29.5), xytext=(D/2, -29.5), arrowprops=dict(arrowstyle='<->', lw=0.6))
    ax.text(0, -28.6, '$D$', fontsize=7, ha='center')
    ax.annotate('', xy=(D/2 - prm.a/2, -26.5), xytext=(D/2 + prm.a/2, -26.5), arrowprops=dict(arrowstyle='<->', lw=0.6))
    ax.text(D/2, -25.6, '$a$', fontsize=7, ha='center')
    ax.set_title('(a) top view', fontsize=8)
    # cross-section
    ax2.add_patch(Rectangle((-9, 0), 18, prm.h, fc='#f3e6c8', ec='k', lw=0.5))
    ax2.plot([-9, 9], [0, 0], 'k-', lw=2); ax2.plot([-9, 9], [prm.h, prm.h], 'k-', lw=2)
    for xc in (-D / 2, D / 2):
        for xr in (xc - prm.a / 2, xc + prm.a / 2):
            ax2.add_patch(Rectangle((xr - prm.d / 2, 0), prm.d, prm.h, fc='k'))
        ax2.add_patch(Rectangle((xc + 0.2 - 0.2, prm.h - 0.03), 0.4, 0.08, fc='white', ec='none'))
    ax2.add_patch(Rectangle((-0.15, prm.h - 0.03), 0.3, 0.08, fc='white', ec='r', lw=0.6))
    ax2.set_xlim(-9.5, 9.5); ax2.set_ylim(-0.3, 1.3); ax2.set_yticks([0, 0.5])
    ax2.set_xlabel('x (mm)'); ax2.set_ylabel('z (mm)')
    ax2.text(0, 0.85, f'$\\varepsilon_r$ = {prm.er}, $h$ = {prm.h} mm, $d$ = {prm.d} mm, $p$ = {prm.p} mm', ha='center', fontsize=6.5)
    ax2.set_title('(b) cross-section', fontsize=8)
    fig.tight_layout()
    fig.savefig(f'{OUT}/fig_geometry.pdf'); plt.close(fig)


def fig_validation():
    tl = json.load(open('results/throughline.json'))
    f = np.array(tl['f']); bfd = np.array(tl['beta_fdtd']); ban = np.array(tl['beta_an'])
    a_fd = des['a_eff_fdtd_mm'] * 1e-3
    bfit = siw.beta(f, a_fd, prm.er)
    fig, ax = plt.subplots(1, 2, figsize=(7.0, 2.6))
    ax[0].plot(f / 1e9, ban, 'k--', label=f'analytical, $a_{{eff}}$ = {prm.a_eff:.2f} mm (Cassivi)')
    ax[0].plot(f / 1e9, bfit, 'b:', label=f'analytical, $a_{{eff}}$ = {a_fd*1e3:.2f} mm (fit)')
    ax[0].plot(f[::8] / 1e9, bfd[::8], 'ro', ms=3, mfc='none', label='FDTD (through line)')
    ax[0].set_xlabel('frequency (GHz)'); ax[0].set_ylabel(r'$\beta$ (rad/m)'); ax[0].legend(); ax[0].grid(alpha=.3)
    ax[0].set_title('(a) TE$_{10}$ dispersion', fontsize=8)
    ax[1].plot(f / 1e9, np.array(tl['S21_dB']), 'r-', label='FDTD $|S_{21}|$, 32 mm')
    ad = siw.dielectric_attenuation(f, a_fd, prm.er, prm.tand) * 0.032 * 20 * np.log10(np.e)
    ac = siw.conductor_attenuation(f, a_fd, prm.h * 1e-3, prm.er) * 0.032 * 20 * np.log10(np.e)
    ax[1].plot(f / 1e9, -ad, 'k--', label='dielectric loss (analytical)')
    ax[1].plot(f / 1e9, -ad - ac, 'k:', label='dielectric + conductor (analytical)')
    ax[1].set_ylim(-1.0, 0.2); ax[1].set_xlabel('frequency (GHz)'); ax[1].set_ylabel('dB'); ax[1].legend(); ax[1].grid(alpha=.3)
    ax[1].set_title('(b) through-line transmission', fontsize=8)
    fig.savefig(f'{OUT}/fig_validation.pdf'); plt.close(fig)


def fig_slot():
    d = json.load(open('results/slot_characterization.json')); f = np.array(d['f'])
    fig, ax = plt.subplots(1, 3, figsize=(7.4, 2.5), gridspec_kw=dict(wspace=0.45))
    for key, ls in (('3.7,0.2', '-'), ('3.9,0.2', '--'), ('4.1,0.2', '-.'), ('4.5,0.2', ':')):
        if key in d['cases']:
            c = d['cases'][key]
            ax[0].plot(f / 1e9, c['g'], 'b', ls=ls, label=f'$L_s$ = {c["L"]} mm')
            ax[0].plot(f / 1e9, c['b'], 'r', ls=ls)
    ax[0].set_xlim(24, 32); ax[0].set_xlabel('frequency (GHz)'); ax[0].set_ylabel('$g$ (blue), $b$ (red)')
    ax[0].legend(); ax[0].grid(alpha=.3); ax[0].set_title('(a) slot admittance, $x_0$ = 0.2 mm', fontsize=8)
    fit = des['fit']; L = np.array(fit['L']); fr = np.array(fit['f_res']); xo = np.array(fit['xo']); gr = np.array(fit['g_res'])
    for x, m in ((0.2, 'o'), (0.3, 's'), (0.5, '^')):
        k = xo == x
        ax[1].plot(L[k], fr[k] / 1e9, m, ms=4, label=f'$x_0$ = {x} mm')
    ax[1].axhline(28, color='k', lw=0.5, ls='--'); ax[1].set_xlabel('slot length $L_s$ (mm)'); ax[1].set_ylabel('resonant frequency (GHz)')
    ax[1].legend(); ax[1].grid(alpha=.3); ax[1].set_title('(b) resonance vs. length', fontsize=8)
    k = np.abs(fr - 28e9) < 1.5e9
    xx = np.linspace(0, 0.7, 100)
    ax[2].plot(xx, des['K'] * np.sin(np.pi * xx / prm.a) ** 2, 'k-', label=f'$K\\sin^2(\\pi x_0/a)$, $K$ = {des["K"]:.1f}')
    ax[2].plot(xo[k], gr[k], 'ro', ms=4, label='FDTD, $|f_{res}-28|<1.5$ GHz')
    ax[2].axhline(1 / des['N'], color='b', lw=0.6, ls='--'); ax[2].text(0.45, 1 / des['N'] + 0.06, f'$1/N$ = {1/des["N"]:.2f}', color='b', fontsize=7)
    ax[2].set_xlabel('slot offset $x_0$ (mm)'); ax[2].set_ylabel('resonant conductance $g_r$'); ax[2].legend(); ax[2].grid(alpha=.3)
    ax[2].set_title('(c) offset law', fontsize=8)
    fig.savefig(f'{OUT}/fig_slot.pdf'); plt.close(fig)


def fig_cma_cal():
    c = json.load(open('results/cma_calibration.json'))
    fig, ax = plt.subplots(1, 2, figsize=(7.0, 2.5))
    f = np.array(c['slot_ca_free']['f']); ax[0].plot(f / 1e9, c['slot_ca_free']['ca'], 'k--', label=f'slot {c["slot_L"]} mm, $\\varepsilon_{{bg}}$ = 1')
    f = np.array(c['slot_ca_eff']['f']); ax[0].plot(f / 1e9, c['slot_ca_eff']['ca'], 'b-', label=f'slot {c["slot_L"]} mm, $\\varepsilon_{{bg}}$ = {c["eps_eff"]:.2f}')
    ax[0].axhline(180, color='k', lw=0.5); ax[0].axvline(28, color='k', lw=0.5, ls=':')
    ax[0].set_xlabel('frequency (GHz)'); ax[0].set_ylabel('characteristic angle (deg)'); ax[0].legend(); ax[0].grid(alpha=.3)
    ax[0].set_title('(a) aperture-mode calibration', fontsize=8)
    Ls = sorted(float(k) for k in c['slit_curve']); fr = [c['slit_curve'][str(L) if str(L) in c['slit_curve'] else f'{L:.1f}']['f_res'] for L in Ls]
    ax[1].plot(Ls, np.array(fr) / 1e9, 'ro-', ms=4)
    ax[1].axhline(28, color='k', lw=0.5, ls='--'); ax[1].axvline(c['slit_L_28GHz'], color='r', lw=0.5, ls=':')
    ax[1].text(c['slit_L_28GHz'] + 0.03, 24.3, f'$L_d$ = {c["slit_L_28GHz"]:.2f} mm', color='r', fontsize=7)
    ax[1].set_xlabel('slit length $L_d$ (mm)'); ax[1].set_ylabel('CM resonance of slit mode (GHz)'); ax[1].grid(alpha=.3)
    ax[1].set_title('(b) slit design curve (CA = 180$^\\circ$)', fontsize=8)
    fig.savefig(f'{OUT}/fig_cma_cal.pdf'); plt.close(fig)


def norm_gci(data):
    """Area-normalised gap-crossing index: mean |Jx| in the channel / mean |J| on the plate."""
    cen = data['cen']; area = data['area']; D_ = float(data['D'])
    gap = np.abs(cen[:, 0]) < (D_ / 2 - 2.8 - 0.2) * 1e-3
    return data['gci'] * area.sum() / area[gap].sum()


def fig_cma_plate():
    r = np.load('results/cma_ref.npz'); s = np.load('results/cma_slits.npz')
    fig, ax = plt.subplots(1, 2, figsize=(7.0, 2.6), sharey=True)
    for data, a, ttl in ((r, ax[0], '(a) without slits'), (s, ax[1], '(b) with slits')):
        fr = data['freqs'] / 1e9; gci = norm_gci(data)
        for n in range(gci.shape[1]):
            a.plot(fr, gci[:, n], color='0.6', lw=0.7)
        a.plot(fr, gci.max(axis=1), 'r-', lw=1.6, label='maximum over the 10 modes')
        a.plot(fr, gci.mean(axis=1), 'b--', lw=1.2, label='mean over the 10 modes')
        a.set_xlabel('frequency (GHz)'); a.grid(alpha=.3); a.axvline(28, color='k', lw=0.5, ls=':')
        a.set_title(ttl, fontsize=8); a.legend(fontsize=6)
    ax[0].set_ylabel('gap-crossing index')
    fig.savefig(f'{OUT}/fig_cma_ms.pdf'); plt.close(fig)
    # current maps of the strongest gap-crossing mode at 28 GHz
    fig, ax = plt.subplots(1, 2, figsize=(5.0, 4.2))
    for data, a, ttl in ((r, ax[0], '(a) without slits'), (s, ax[1], '(b) with slits')):
        fr = data['freqs'] / 1e9; i28 = np.argmin(abs(fr - 28)); gci = norm_gci(data); n = int(np.argmax(gci[i28]))
        J = data['J28'][n]; Jm = np.linalg.norm(J, axis=1); Jm /= Jm.max()
        nodes = data['nodes'] * 1e3; tris = data['tris']
        tp = a.tripcolor(nodes[:, 0], nodes[:, 1], tris, facecolors=Jm, cmap='inferno', vmin=0, vmax=1)
        cen = data['cen'] * 1e3; k = np.arange(0, len(cen), 6)
        a.quiver(cen[k, 0], cen[k, 1], J[k, 0].real, J[k, 1].real, color='c', scale=np.abs(J).max() * 25, width=0.004)
        a.set_aspect('equal'); a.set_xlabel('x (mm)'); a.set_ylabel('y (mm)')
        a.set_title(f'{ttl}\nmode {n+1}, $\\lambda_n$ = {data["lam"][i28, n]:.2f}, GCI = {gci[i28, n]:.2f}', fontsize=7)
    fig.colorbar(tp, ax=ax, shrink=0.6, label='|J| (normalised)')
    fig.savefig(f'{OUT}/fig_cma_currents.pdf'); plt.close(fig)


def fig_aperture():
    import glob
    fig, ax = plt.subplots(1, 2, figsize=(7.0, 2.6))
    for fn, ls in zip(sorted(glob.glob('results/aperture_coupling_D7.0_eps*.json')), ('r-o', 'b--s', 'g:^')):
        c = json.load(open(fn)); sw = c['sweep']
        Ls = sorted(float(k) for k in sw if k not in ('0', 'cont')); y = [sw[f'{L:.1f}'] for L in Ls]
        ax[0].plot(Ls, 20 * np.log10(np.array(y) / sw['0']), ls, ms=4, label=f'$\\varepsilon_{{eff}}$ = {c["eps_eff"]:.2f}')
        if 'cont' in sw:
            ax[0].axhline(20 * np.log10(sw['cont'] / sw['0']), color=ls[0], lw=0.6, ls=':')
    ax[0].axhline(0, color='k', lw=0.5); ax[0].set_xlabel('slit length $L_d$ (mm)'); ax[0].set_ylabel('$|y_{12}|$ change (dB)')
    ax[0].legend(); ax[0].grid(alpha=.3); ax[0].set_title('(a) sub-array mutual admittance, 28 GHz', fontsize=8)
    fn = sorted(glob.glob('results/aperture_cma_D7.0_eps*.npz'))
    if fn:
        d = np.load(fn[-1]); fr = d['freqs'] / 1e9
        from siwlab.mom import modal_significance
        for n in range(d['lam'].shape[1]):
            sf = d['slitfrac'][np.argmin(abs(fr - 28)), n]
            ax[1].plot(fr, modal_significance(d['lam'][:, n]), color=plt.cm.coolwarm(sf), lw=1.0 + 1.5 * sf)
        sm = plt.cm.ScalarMappable(cmap='coolwarm', norm=plt.Normalize(0, 1)); sm.set_array([])
        fig.colorbar(sm, ax=ax[1], label='slit current fraction at 28 GHz')
        ax[1].set_xlabel('frequency (GHz)'); ax[1].set_ylabel('modal significance'); ax[1].set_ylim(0, 1.05); ax[1].grid(alpha=.3)
        ax[1].axvline(28, color='k', lw=0.5, ls=':'); ax[1].set_title(f'(b) aperture-system modes, $L_d$ = {float(d["Ld"]):.1f} mm', fontsize=8)
    fig.savefig(f'{OUT}/fig_aperture.pdf'); plt.close(fig)


def fig_sparams(tags):
    fig, ax = plt.subplots(1, 2, figsize=(7.0, 2.6))
    for tag, lab, ls in tags:
        p = np.load(f'results/array_{tag}_post.npz'); f = p['f'] / 1e9
        ax[0].plot(f, db(p['S11']), ls, label=lab)
        ax[1].plot(f, db(p['S21']), ls, label=lab)
    for a, t in zip(ax, ('(a) $|S_{11}|$', '(b) $|S_{21}|$ (isolation)')):
        a.set_xlabel('frequency (GHz)'); a.set_ylabel('dB'); a.grid(alpha=.3); a.legend(); a.set_title(t, fontsize=8); a.set_xlim(25, 31)
    ax[0].axhline(-10, color='k', lw=0.5, ls=':'); ax[1].axhline(-20, color='k', lw=0.5, ls=':')
    fig.savefig(f'{OUT}/fig_sparams.pdf'); plt.close(fig)


def ecc(E1th, E1ph, E2th, E2ph, th):
    w = np.sin(th)
    num = np.abs(np.sum((E1th * np.conj(E2th) + E1ph * np.conj(E2ph)) * w)) ** 2
    den = np.sum((abs(E1th) ** 2 + abs(E1ph) ** 2) * w) * np.sum((abs(E2th) ** 2 + abs(E2ph) ** 2) * w)
    return num / den


def fig_radiation(tags):
    fig, ax = plt.subplots(1, 3, figsize=(7.2, 2.5))
    for tag, lab, ls in tags:
        p = np.load(f'results/array_{tag}_post.npz'); fd = p['dft_f'] / 1e9
        i28 = np.argmin(abs(fd - 28)); G = 10 * np.log10(np.maximum(p['G'][i28], 1e-6))
        th = np.degrees(p['theta']); ph = p['phi']
        iE = np.argmin(abs(ph - np.pi / 2)); iE2 = np.argmin(abs(ph - 3 * np.pi / 2))  # y-z plane (E-plane of the array, along slots)
        iH = np.argmin(abs(ph - 0)); iH2 = np.argmin(abs(ph - np.pi))              # x-z plane
        ax[0].plot(np.r_[-th[::-1], th], np.r_[G[::-1, iE2], G[:, iE]], ls, label=lab)
        ax[1].plot(np.r_[-th[::-1], th], np.r_[G[::-1, iH2], G[:, iH]], ls, label=lab)
        # gain vs frequency and ECC (port 2 pattern = mirror image in x: phi -> pi - phi)
        Gmax = [10 * np.log10(p['Gr'][i].max()) for i in range(len(fd))]
        TH, PH = np.meshgrid(p['theta'], p['phi'], indexing='ij')
        ec = []
        for i in range(len(fd)):
            Eth, Eph = p['Eth'][i], p['Eph'][i]
            # mirror in x: phi -> pi - phi ; E_theta even, E_phi odd under this reflection for x-polarised aperture fields
            idx = np.argmin(abs(((np.pi - ph)[None, :] - ph[:, None] + np.pi) % (2 * np.pi) - np.pi), axis=1)
            E2th = Eth[:, idx]; E2ph = -Eph[:, idx]
            ec.append(ecc(Eth, Eph, E2th, E2ph, TH))
        ax[2].plot(fd, Gmax, ls, label=lab + ' (realised gain)')
        ax2 = ax[2].twinx() if tag == tags[0][0] else ax2
        ax2.semilogy(fd, ec, ls.replace('-', ':') if '-' in ls else ls, color='g', label=lab + ' (ECC)')
    for a, t in zip(ax[:2], ('(a) y-z plane, 28 GHz', '(b) x-z plane, 28 GHz')):
        a.set_xlabel(r'$\theta$ (deg)'); a.set_ylabel('gain (dBi)'); a.set_ylim(-20, 15); a.grid(alpha=.3); a.legend(); a.set_title(t, fontsize=8)
    ax[2].set_xlabel('frequency (GHz)'); ax[2].set_ylabel('realised gain (dBi)'); ax2.set_ylabel('ECC', color='g'); ax[2].grid(alpha=.3)
    ax[2].legend(loc='lower left', fontsize=6); ax[2].set_title('(c) gain and ECC', fontsize=8)
    fig.savefig(f'{OUT}/fig_radiation.pdf'); plt.close(fig)


if __name__ == '__main__':
    cal = json.load(open('results/cma_calibration.json')) if os.path.exists('results/cma_calibration.json') else None
    Ls = cal['slit_L_28GHz'] if cal else 4.0
    for fn, args in ((fig_geometry, (Ls,)), (fig_validation, ()), (fig_slot, ()), (fig_cma_cal, ()), (fig_cma_plate, ()), (fig_aperture, ()),
                     (fig_sparams, ([('nodecouple', 'without slits', 'k--'), ('decoupled', 'with CM-guided slits', 'r-')],)),
                     (fig_radiation, ([('nodecouple', 'without slits', 'k--'), ('decoupled', 'with slits', 'r-')],))):
        try:
            fn(*args); print('ok', fn.__name__)
        except Exception as e:
            print('skip', fn.__name__, '->', repr(e))
