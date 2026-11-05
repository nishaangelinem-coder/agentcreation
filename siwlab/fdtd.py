"""
Numba-accelerated 3-D FDTD solver on a nonuniform Yee grid with CPML.

Designed for planar SIW structures:
  * bottom boundary (z = 0) is PEC (ground plane), all other boundaries CPML
  * cell-based permittivity / conductivity, edge-based PEC masks
  * TE10 SIW waveguide ports (soft modal source + modal voltage/current probes)
  * running DFT of tangential E in the apertures of a PEC sheet for the
    aperture (magnetic-current) near-to-far-field transformation.

Units: SI (metres, seconds). Geometry helpers accept millimetres.
"""
import numpy as np
import numba as nb
from .constants import c0, eps0, mu0, eta0


# --------------------------------------------------------------------------
# nonuniform axis construction
# --------------------------------------------------------------------------
def graded_axis(segments):
    """Build node coordinates from consecutive segments.

    segments: list of tuples
        (x0, x1, dx)            uniform cells of size ~dx (rounded to integer count)
        (x0, x1, dx0, dx1)      geometric grading from dx0 at x0 to dx1 at x1
    Returns node coordinate array (length Ncells+1).
    """
    nodes = [segments[0][0]]
    for seg in segments:
        x0, x1 = seg[0], seg[1]
        L = x1 - x0
        if len(seg) == 3:
            n = max(1, int(round(L / seg[2])))
            d = np.full(n, L / n)
        else:
            d0, d1 = seg[2], seg[3]
            # geometric progression d0, d0 r, ..., with sum L and last ~ d1
            if abs(d1 - d0) < 1e-12:
                n = max(1, int(round(L / d0)))
                d = np.full(n, L / n)
            else:
                r = d1 / d0
                # number of cells n such that d0 (r^n - 1)/(r-1) = L with r = (d1/d0)^(1/(n-1))
                n = max(2, int(round(L / (0.5 * (d0 + d1)))))
                for _ in range(50):
                    rr = (d1 / d0) ** (1.0 / (n - 1))
                    Lt = d0 * (rr ** n - 1) / (rr - 1)
                    if Lt > L and n > 2:
                        n -= 1
                    elif Lt < L:
                        n += 1
                    else:
                        break
                rr = (d1 / d0) ** (1.0 / (n - 1))
                d = d0 * rr ** np.arange(n)
                d *= L / d.sum()
        nodes.extend(list(nodes[-1] + np.cumsum(d)))
    return np.array(nodes)


# --------------------------------------------------------------------------
# Numba kernels
# --------------------------------------------------------------------------
@nb.njit(parallel=True, fastmath=True, cache=True)
def update_E(Ex, Ey, Ez, Hx, Hy, Hz,
             cax, cbx, cay, cby, caz, cbz,
             idxd, idyd, idzd,
             pEx_y, pEx_z, pEy_x, pEy_z, pEz_x, pEz_y,
             bex, cex, kex, bey, cey, key, bez, cez, kez):
    Nx = Ex.shape[0]; Ny = Ey.shape[1]; Nz = Ez.shape[2]
    # Ex : i in [0,Nx), j in [1,Ny), k in [1,Nz)
    for i in nb.prange(Nx):
        for j in range(1, Ny):
            for k in range(1, Nz):
                dHz = (Hz[i, j, k] - Hz[i, j - 1, k]) * idyd[j]
                dHy = (Hy[i, j, k] - Hy[i, j, k - 1]) * idzd[k]
                py = bey[j] * pEx_y[i, j, k] + cey[j] * dHz
                pz = bez[k] * pEx_z[i, j, k] + cez[k] * dHy
                pEx_y[i, j, k] = py
                pEx_z[i, j, k] = pz
                Ex[i, j, k] = cax[i, j, k] * Ex[i, j, k] + cbx[i, j, k] * (
                    dHz * key[j] + py - dHy * kez[k] - pz)
    # Ey : i in [1,Nx), j in [0,Ny), k in [1,Nz)
    for i in nb.prange(1, Nx):
        for j in range(Ny):
            for k in range(1, Nz):
                dHx = (Hx[i, j, k] - Hx[i, j, k - 1]) * idzd[k]
                dHz = (Hz[i, j, k] - Hz[i - 1, j, k]) * idxd[i]
                pz = bez[k] * pEy_z[i, j, k] + cez[k] * dHx
                px = bex[i] * pEy_x[i, j, k] + cex[i] * dHz
                pEy_z[i, j, k] = pz
                pEy_x[i, j, k] = px
                Ey[i, j, k] = cay[i, j, k] * Ey[i, j, k] + cby[i, j, k] * (
                    dHx * kez[k] + pz - dHz * kex[i] - px)
    # Ez : i in [1,Nx), j in [1,Ny), k in [0,Nz)
    for i in nb.prange(1, Nx):
        for j in range(1, Ny):
            for k in range(Nz):
                dHy = (Hy[i, j, k] - Hy[i - 1, j, k]) * idxd[i]
                dHx = (Hx[i, j, k] - Hx[i, j - 1, k]) * idyd[j]
                px = bex[i] * pEz_x[i, j, k] + cex[i] * dHy
                py = bey[j] * pEz_y[i, j, k] + cey[j] * dHx
                pEz_x[i, j, k] = px
                pEz_y[i, j, k] = py
                Ez[i, j, k] = caz[i, j, k] * Ez[i, j, k] + cbz[i, j, k] * (
                    dHy * kex[i] + px - dHx * key[j] - py)


@nb.njit(parallel=True, fastmath=True, cache=True)
def update_H(Ex, Ey, Ez, Hx, Hy, Hz, dtmu,
             idx, idy, idz,
             pHx_y, pHx_z, pHy_x, pHy_z, pHz_x, pHz_y,
             bhx, chx, khx, bhy, chy, khy, bhz, chz, khz):
    Nx = Ex.shape[0]; Ny = Ey.shape[1]; Nz = Ez.shape[2]
    # Hx : i in [0,Nx], j in [0,Ny), k in [0,Nz)
    for i in nb.prange(Nx + 1):
        for j in range(Ny):
            for k in range(Nz):
                dEz = (Ez[i, j + 1, k] - Ez[i, j, k]) * idy[j]
                dEy = (Ey[i, j, k + 1] - Ey[i, j, k]) * idz[k]
                py = bhy[j] * pHx_y[i, j, k] + chy[j] * dEz
                pz = bhz[k] * pHx_z[i, j, k] + chz[k] * dEy
                pHx_y[i, j, k] = py
                pHx_z[i, j, k] = pz
                Hx[i, j, k] -= dtmu * (dEz * khy[j] + py - dEy * khz[k] - pz)
    # Hy : i in [0,Nx), j in [0,Ny], k in [0,Nz)
    for i in nb.prange(Nx):
        for j in range(Ny + 1):
            for k in range(Nz):
                dEx = (Ex[i, j, k + 1] - Ex[i, j, k]) * idz[k]
                dEz = (Ez[i + 1, j, k] - Ez[i, j, k]) * idx[i]
                pz = bhz[k] * pHy_z[i, j, k] + chz[k] * dEx
                px = bhx[i] * pHy_x[i, j, k] + chx[i] * dEz
                pHy_z[i, j, k] = pz
                pHy_x[i, j, k] = px
                Hy[i, j, k] -= dtmu * (dEx * khz[k] + pz - dEz * khx[i] - px)
    # Hz : i in [0,Nx), j in [0,Ny), k in [0,Nz]
    for i in nb.prange(Nx):
        for j in range(Ny):
            for k in range(Nz + 1):
                dEy = (Ey[i + 1, j, k] - Ey[i, j, k]) * idx[i]
                dEx = (Ex[i, j + 1, k] - Ex[i, j, k]) * idy[j]
                px = bhx[i] * pHz_x[i, j, k] + chx[i] * dEy
                py = bhy[j] * pHz_y[i, j, k] + chy[j] * dEx
                pHz_x[i, j, k] = px
                pHz_y[i, j, k] = py
                Hz[i, j, k] -= dtmu * (dEy * khx[i] + px - dEx * khy[j] - py)


@nb.njit(cache=True)
def accumulate_dft(acc, field, idx_i, idx_j, idx_k, phasors):
    # acc[nf, npts] += field[idx] * phasors[nf]
    nf = phasors.shape[0]
    n = idx_i.shape[0]
    for f in range(nf):
        ph = phasors[f]
        for p in range(n):
            acc[f, p] += field[idx_i[p], idx_j[p], idx_k[p]] * ph


# --------------------------------------------------------------------------
# CPML parameter construction for one axis
# --------------------------------------------------------------------------
def cpml_axis(nodes, npml_lo, npml_hi, dt, m=3, kappa_max=1.0, alpha_max=0.05, sigma_scale=1.0):
    """Return (b_e, c_e, k_e) for integer (E) positions and (b_h, c_h, k_h)
    for half (H) positions along an axis with node coordinates `nodes`."""
    N = len(nodes) - 1
    xe = nodes.copy()
    xh = 0.5 * (nodes[1:] + nodes[:-1])
    d = np.diff(nodes)

    def params(x):
        sig = np.zeros_like(x); kap = np.ones_like(x); alp = np.zeros_like(x)
        if npml_lo > 0:
            x0 = nodes[npml_lo]; L = x0 - nodes[0]
            dpml = d[:npml_lo].mean()
            smax = sigma_scale * 0.8 * (m + 1) / (eta0 * dpml)
            rho = np.clip((x0 - x) / L, 0, 1)
            msk = rho > 0
            sig[msk] = smax * rho[msk] ** m
            kap[msk] = 1 + (kappa_max - 1) * rho[msk] ** m
            alp[msk] = alpha_max * (1 - rho[msk])
        if npml_hi > 0:
            x1 = nodes[N - npml_hi]; L = nodes[-1] - x1
            dpml = d[N - npml_hi:].mean()
            smax = sigma_scale * 0.8 * (m + 1) / (eta0 * dpml)
            rho = np.clip((x - x1) / L, 0, 1)
            msk = rho > 0
            sig[msk] = smax * rho[msk] ** m
            kap[msk] = 1 + (kappa_max - 1) * rho[msk] ** m
            alp[msk] = alpha_max * (1 - rho[msk])
        b = np.exp(-(sig / kap + alp) * dt / eps0)
        c = np.zeros_like(x)
        msk = sig > 0
        c[msk] = sig[msk] * (b[msk] - 1) / (sig[msk] * kap[msk] + kap[msk] ** 2 * alp[msk])
        return b, c, 1.0 / kap

    be, ce, ke = params(xe)
    bh, ch, kh = params(xh)
    # inside the computational domain, b=0 (no memory), c=0, k=1
    be[np.abs(ce) == 0] = 0.0
    bh[np.abs(ch) == 0] = 0.0
    return (be.astype(np.float32), ce.astype(np.float32), ke.astype(np.float32),
            bh.astype(np.float32), ch.astype(np.float32), kh.astype(np.float32))


# --------------------------------------------------------------------------
# Simulation container
# --------------------------------------------------------------------------
class FDTD:
    def __init__(self, x, y, z, npml=(10, 10, 10, 10, 0, 10), courant=0.99, dtype=np.float32):
        """x, y, z: node coordinate arrays (m). npml: (xlo, xhi, ylo, yhi, zlo, zhi)."""
        self.x, self.y, self.z = x, y, z
        self.Nx, self.Ny, self.Nz = len(x) - 1, len(y) - 1, len(z) - 1
        self.dx, self.dy, self.dz = np.diff(x), np.diff(y), np.diff(z)
        self.xd = np.r_[self.dx[0], 0.5 * (self.dx[1:] + self.dx[:-1]), self.dx[-1]]
        self.yd = np.r_[self.dy[0], 0.5 * (self.dy[1:] + self.dy[:-1]), self.dy[-1]]
        self.zd = np.r_[self.dz[0], 0.5 * (self.dz[1:] + self.dz[:-1]), self.dz[-1]]
        dmin = min(self.dx.min(), self.dy.min(), self.dz.min())
        self.dt = courant / (c0 * np.sqrt(1 / self.dx.min() ** 2 + 1 / self.dy.min() ** 2 + 1 / self.dz.min() ** 2))
        self.npml = npml
        self.dtype = dtype
        Nx, Ny, Nz = self.Nx, self.Ny, self.Nz
        # material (cell based)
        self.epsr = np.ones((Nx, Ny, Nz))
        self.sigma = np.zeros((Nx, Ny, Nz))
        # PEC edge masks
        self.pecx = np.zeros((Nx, Ny + 1, Nz + 1), bool)
        self.pecy = np.zeros((Nx + 1, Ny, Nz + 1), bool)
        self.pecz = np.zeros((Nx + 1, Ny + 1, Nz), bool)
        self.sources = []
        self.probes = []
        self.dft_monitors = []
        # node position helpers
        self.xh = 0.5 * (x[1:] + x[:-1]); self.yh = 0.5 * (y[1:] + y[:-1]); self.zh = 0.5 * (z[1:] + z[:-1])

    # ---------------- geometry helpers (metres) ----------------
    def idx(self, axis, val):
        arr = {'x': self.x, 'y': self.y, 'z': self.z}[axis]
        return int(np.argmin(np.abs(arr - val)))

    def set_box_material(self, x0, x1, y0, y1, z0, z1, epsr, sigma=0.0):
        i0, i1 = np.searchsorted(self.xh, x0), np.searchsorted(self.xh, x1)
        j0, j1 = np.searchsorted(self.yh, y0), np.searchsorted(self.yh, y1)
        k0, k1 = np.searchsorted(self.zh, z0), np.searchsorted(self.zh, z1)
        self.epsr[i0:i1, j0:j1, k0:k1] = epsr
        self.sigma[i0:i1, j0:j1, k0:k1] = sigma

    def add_pec_sheet_z(self, k, apertures=()):
        """PEC sheet on node plane k with rectangular apertures [(x0,x1,y0,y1),...]."""
        X, Y = np.meshgrid(self.xh, self.y, indexing='ij')  # Ex nodes
        mx = np.ones(X.shape, bool)
        tol = 1e-7
        for (x0, x1, y0, y1) in apertures:
            mx &= ~((X > x0 + tol) & (X < x1 - tol) & (Y > y0 + tol) & (Y < y1 - tol))
        self.pecx[:, :, k] |= mx
        X, Y = np.meshgrid(self.x, self.yh, indexing='ij')  # Ey nodes
        my = np.ones(X.shape, bool)
        for (x0, x1, y0, y1) in apertures:
            my &= ~((X > x0 + tol) & (X < x1 - tol) & (Y > y0 + tol) & (Y < y1 - tol))
        self.pecy[:, :, k] |= my

    def add_pec_cylinder_z(self, xc, yc, r, k0, k1):
        """PEC cylinder (axis along z) between node planes k0..k1. Cell-centred
        staircase: every cell whose centre lies inside the circle becomes PEC
        (all its edges), which conserves the via cross-section area."""
        X, Y = np.meshgrid(self.xh, self.yh, indexing='ij')
        cells = (X - xc) ** 2 + (Y - yc) ** 2 <= r * r
        ci, cj = np.where(cells)
        if len(ci) == 0:
            # tiny via: at least the nearest node
            ci = np.array([np.argmin(np.abs(self.xh - xc))]); cj = np.array([np.argmin(np.abs(self.yh - yc))])
        for di, dj in ((0, 0), (1, 0), (0, 1), (1, 1)):
            self.pecz[ci + di, cj + dj, k0:k1] = True
        for dj in (0, 1):
            self.pecx[ci, cj + dj, k0:k1 + 1] = True
        for di in (0, 1):
            self.pecy[ci + di, cj, k0:k1 + 1] = True

    def add_pec_wall_y(self, j, x0, x1, k0, k1):
        """PEC wall on node plane j spanning x0..x1 and node planes k0..k1."""
        ix = (self.xh > x0) & (self.xh < x1)
        self.pecx[ix, j, k0:k1 + 1] = True
        ix = (self.x >= x0) & (self.x <= x1)
        self.pecz[ix, j, k0:k1] = True

    def add_pec_wall_x(self, i, y0, y1, k0, k1):
        iy = (self.yh > y0) & (self.yh < y1)
        self.pecy[i, iy, k0:k1 + 1] = True
        iy = (self.y >= y0) & (self.y <= y1)
        self.pecz[i, iy, k0:k1] = True

    # ---------------- SIW TE10 port ----------------
    def add_siw_port(self, name, j_src, j_probe, xc, a_eff, k0, k1, direction=+1):
        """Soft TE10 source on plane j_src and modal probes on plane j_probe.
        The mode profile sin(pi (x - (xc - a/2))/a) spans node planes k0..k1 (Ez edges k0..k1-1)."""
        xl = xc - a_eff / 2
        ix = np.where((self.x > xl) & (self.x < xl + a_eff))[0]
        prof = np.sin(np.pi * (self.x[ix] - xl) / a_eff)
        ks = np.arange(k0, k1)
        port = dict(name=name, j_src=j_src, j_probe=j_probe, ix=ix, prof=prof, ks=ks,
                    a_eff=a_eff, direction=direction)
        # probe weights: sum_{i,k} E * prof * dxd * dz
        W = prof[:, None] * self.xd[ix][:, None] * self.dz[ks][None, :]
        port['W'] = W
        port['norm'] = float(np.sum(prof[:, None] ** 2 * self.xd[ix][:, None] * self.dz[ks][None, :]))
        return port

    # ---------------- build coefficient arrays ----------------
    def build(self, f_sigma=None):
        Nx, Ny, Nz = self.Nx, self.Ny, self.Nz
        dt = self.dt
        eps = self.epsr; sig = self.sigma

        def edge_avg(arr, axis):
            # average of the 4 cells sharing an edge directed along `axis`
            a = np.zeros((Nx + 1, Ny + 1, Nz + 1))
            # pad cells
            P = np.pad(arr, 1, mode='edge')
            if axis == 0:   # Ex at (i+1/2, j, k): cells (i, j-1..j, k-1..k)
                a = 0.25 * (P[1:-1, :-1, :-1] + P[1:-1, 1:, :-1] + P[1:-1, :-1, 1:] + P[1:-1, 1:, 1:])
                return a  # shape (Nx, Ny+1, Nz+1)
            if axis == 1:
                a = 0.25 * (P[:-1, 1:-1, :-1] + P[1:, 1:-1, :-1] + P[:-1, 1:-1, 1:] + P[1:, 1:-1, 1:])
                return a  # (Nx+1, Ny, Nz+1)
            a = 0.25 * (P[:-1, :-1, 1:-1] + P[1:, :-1, 1:-1] + P[:-1, 1:, 1:-1] + P[1:, 1:, 1:-1])
            return a      # (Nx+1, Ny+1, Nz)

        self.ca = []; self.cb = []
        for ax, pec in zip((0, 1, 2), (self.pecx, self.pecy, self.pecz)):
            e = edge_avg(eps, ax) * eps0
            s = edge_avg(sig, ax)
            ca = (1 - s * dt / (2 * e)) / (1 + s * dt / (2 * e))
            cb = (dt / e) / (1 + s * dt / (2 * e))
            ca[pec] = 0.0; cb[pec] = 0.0
            self.ca.append(ca.astype(self.dtype)); self.cb.append(cb.astype(self.dtype))
        # PML
        xlo, xhi, ylo, yhi, zlo, zhi = self.npml
        self.pmlx = cpml_axis(self.x, xlo, xhi, dt)
        self.pmly = cpml_axis(self.y, ylo, yhi, dt)
        self.pmlz = cpml_axis(self.z, zlo, zhi, dt)
        # fields
        f = self.dtype
        self.Ex = np.zeros((Nx, Ny + 1, Nz + 1), f); self.Ey = np.zeros((Nx + 1, Ny, Nz + 1), f)
        self.Ez = np.zeros((Nx + 1, Ny + 1, Nz), f)
        self.Hx = np.zeros((Nx + 1, Ny, Nz), f); self.Hy = np.zeros((Nx, Ny + 1, Nz), f)
        self.Hz = np.zeros((Nx, Ny, Nz + 1), f)
        self.psiE = [np.zeros_like(self.Ex), np.zeros_like(self.Ex), np.zeros_like(self.Ey),
                     np.zeros_like(self.Ey), np.zeros_like(self.Ez), np.zeros_like(self.Ez)]
        self.psiH = [np.zeros_like(self.Hx), np.zeros_like(self.Hx), np.zeros_like(self.Hy),
                     np.zeros_like(self.Hy), np.zeros_like(self.Hz), np.zeros_like(self.Hz)]
        self.idx = (1.0 / self.dx).astype(f); self.idy = (1.0 / self.dy).astype(f); self.idz = (1.0 / self.dz).astype(f)
        self.idxd = (1.0 / self.xd).astype(f); self.idyd = (1.0 / self.yd).astype(f); self.idzd = (1.0 / self.zd).astype(f)

    # ---------------- run ----------------
    def run(self, nsteps, src_port, src_func, probe_ports, dft_freqs=(), dft_nodes=None,
            progress=None, src_amp=1.0):
        """
        src_port  : port dict for excitation (soft Ez source on j_src)
        src_func  : callable t -> amplitude
        probe_ports: list of port dicts; records V(t), I(t)
        dft_freqs : frequencies for the running DFT of aperture tangential fields
        dft_nodes : dict with 'ex': (ii,jj,kk) index arrays of Ex nodes and 'ey': (...) of Ey nodes
        Returns dict with time axis, V/I records and DFT accumulators.
        """
        f = self.dtype
        dt = self.dt
        nf = len(dft_freqs)
        rec = {p['name']: {'V': np.zeros(nsteps), 'I': np.zeros(nsteps)} for p in probe_ports}
        if nf and dft_nodes is not None:
            accEx = np.zeros((nf, len(dft_nodes['ex'][0])), np.complex128)
            accEy = np.zeros((nf, len(dft_nodes['ey'][0])), np.complex128)
            exi = tuple(np.asarray(v, np.int64) for v in dft_nodes['ex'])
            eyi = tuple(np.asarray(v, np.int64) for v in dft_nodes['ey'])
            w = 2 * np.pi * np.asarray(dft_freqs)
        bex, cex, kex, bhx, chx, khx = self.pmlx
        bey, cey, key, bhy, chy, khy = self.pmly
        bez, cez, kez, bhz, chz, khz = self.pmlz
        dtmu = f(dt / mu0)
        sp = src_port
        ix, ks, prof = sp['ix'], sp['ks'], sp['prof']
        srcprof = (prof[:, None] * np.ones(len(ks))[None, :]).astype(f)
        for n in range(nsteps):
            t = n * dt
            update_H(self.Ex, self.Ey, self.Ez, self.Hx, self.Hy, self.Hz, dtmu,
                     self.idx, self.idy, self.idz, *self.psiH,
                     bhx, chx, khx, bhy, chy, khy, bhz, chz, khz)
            update_E(self.Ex, self.Ey, self.Ez, self.Hx, self.Hy, self.Hz,
                     self.ca[0], self.cb[0], self.ca[1], self.cb[1], self.ca[2], self.cb[2],
                     self.idxd, self.idyd, self.idzd, *self.psiE,
                     bex, cex, kex, bey, cey, key, bez, cez, kez)
            # soft source on Ez
            tE = (n + 1) * dt
            self.Ez[np.ix_(ix, [sp['j_src']], ks)] += (src_amp * src_func(tE) * srcprof)[:, None, :]
            # probes
            for p in probe_ports:
                jp = p['j_probe']
                Ez = self.Ez[np.ix_(p['ix'], [jp], p['ks'])][:, 0, :]
                Hx = 0.5 * (self.Hx[np.ix_(p['ix'], [jp], p['ks'])][:, 0, :] + self.Hx[np.ix_(p['ix'], [jp - 1], p['ks'])][:, 0, :])
                rec[p['name']]['V'][n] = np.sum(Ez * p['W'])
                rec[p['name']]['I'][n] = np.sum(Hx * p['W'])
            if nf and dft_nodes is not None:
                ph = np.exp(-1j * w * tE) * dt
                accumulate_dft(accEx, self.Ex, exi[0], exi[1], exi[2], ph)
                accumulate_dft(accEy, self.Ey, eyi[0], eyi[1], eyi[2], ph)
            if progress and (n % progress == 0):
                emax = float(np.abs(self.Ez).max())
                print(f"  step {n}/{nsteps}  t={tE*1e9:.3f} ns  max|Ez|={emax:.3e}", flush=True)
                if not np.isfinite(emax):
                    raise RuntimeError("FDTD diverged")
        out = dict(t=(np.arange(nsteps) + 1) * dt, dt=dt, rec=rec)
        if nf and dft_nodes is not None:
            out['dft_freqs'] = np.asarray(dft_freqs); out['accEx'] = accEx; out['accEy'] = accEy
        return out


# --------------------------------------------------------------------------
# post-processing helpers
# --------------------------------------------------------------------------
def gaussian_modulated(f0, bw, t0=None):
    """Gaussian-modulated sinusoid; bw = half-width (Hz) of spectrum at 1/e."""
    tau = 1.0 / (np.pi * bw)
    if t0 is None:
        t0 = 4.5 * tau

    def s(t):
        return np.exp(-((t - t0) / tau) ** 2) * np.sin(2 * np.pi * f0 * t)
    return s


def spectrum(sig, dt, freqs):
    """Discrete Fourier transform of a time signal at given frequencies."""
    t = (np.arange(len(sig)) + 1) * dt
    w = 2 * np.pi * np.asarray(freqs)
    return (np.exp(-1j * np.outer(w, t)) @ sig) * dt


def aperture_farfield(freq, nodes_ex, nodes_ey, Ex_f, Ey_f, dS_ex, dS_ey, theta, phi):
    """Far-field from tangential E in apertures of an infinite PEC plane (z = const).
    Ms = -2 n x E,  L = int Ms exp(jk r.r') dS.  Returns (Etheta, Ephi) up to the
    common factor jk e^{-jkr}/(4 pi r) (so |E|^2 -> U with the k^2/(32 pi^2 eta) factor).
    nodes_*: (x, y) arrays of node positions (m)."""
    k = 2 * np.pi * freq / c0
    th = np.asarray(theta)[..., None]; ph = np.asarray(phi)[..., None]
    sx = np.sin(th) * np.cos(ph); sy = np.sin(th) * np.sin(ph)
    # Ms = 2 Ey x_hat - 2 Ex y_hat
    xE, yE = nodes_ex
    xF, yF = nodes_ey
    Lx = np.sum(2 * Ey_f * dS_ey * np.exp(1j * k * (sx * xF + sy * yF)), axis=-1)
    Ly = np.sum(-2 * Ex_f * dS_ex * np.exp(1j * k * (sx * xE + sy * yE)), axis=-1)
    th = th[..., 0]; ph = ph[..., 0]
    Lth = Lx * np.cos(th) * np.cos(ph) + Ly * np.cos(th) * np.sin(ph)
    Lph = -Lx * np.sin(ph) + Ly * np.cos(ph)
    # E_theta = -jk e/(4 pi r) L_phi ; E_phi = +jk e/(4 pi r) L_theta   (Balanis 12-10)
    Eth = -Lph; Eph = Lth
    return Eth, Eph


def radiation_intensity(freq, Eth, Eph):
    k = 2 * np.pi * freq / c0
    return k ** 2 / (32 * np.pi ** 2 * eta0) * (np.abs(Eth) ** 2 + np.abs(Eph) ** 2)


def hemisphere_integral(U, theta, phi):
    """Integrate U(theta,phi) over the upper hemisphere; theta,phi are 2-D meshgrids (ij)."""
    dth = theta[1, 0] - theta[0, 0]; dph = phi[0, 1] - phi[0, 0]
    return np.sum(U * np.sin(theta)) * dth * dph
