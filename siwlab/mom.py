"""RWG electric-field integral-equation Method of Moments for PEC surfaces and
characteristic-mode analysis (Harrington & Mautz generalized eigenproblem
X J = lambda R J).  Homogeneous background medium with relative permittivity
eps_bg (effective-medium treatment of a substrate-loaded aperture structure)."""
import numpy as np
import numba as nb
from scipy.linalg import eigh
from .constants import c0, eps0, mu0
from .mesh import rwg_edges

# 6-point Gauss rule on the unit triangle (degree 4), barycentric coords + weights
_G6 = np.array([
    [0.816847572980459, 0.091576213509771, 0.091576213509771, 0.109951743655322],
    [0.091576213509771, 0.816847572980459, 0.091576213509771, 0.109951743655322],
    [0.091576213509771, 0.091576213509771, 0.816847572980459, 0.109951743655322],
    [0.108103018168070, 0.445948490915965, 0.445948490915965, 0.223381589678011],
    [0.445948490915965, 0.108103018168070, 0.445948490915965, 0.223381589678011],
    [0.445948490915965, 0.445948490915965, 0.108103018168070, 0.223381589678011]])


def _subdivision_rule(level):
    """Barycentric centroids of the 4^level uniform sub-triangles (equal weights)."""
    pts = [((1, 0, 0), (0, 1, 0), (0, 0, 1))]
    for _ in range(level):
        new = []
        for (a, b, c) in pts:
            a, b, c = np.array(a, float), np.array(b, float), np.array(c, float)
            ab, bc, ca = (a + b) / 2, (b + c) / 2, (c + a) / 2
            new += [(a, ab, ca), (ab, b, bc), (ca, bc, c), (ab, bc, ca)]
        pts = new
    cen = np.array([(a + b + c) / 3 for (a, b, c) in pts])
    w = np.full(len(cen), 1.0 / len(cen))
    return np.c_[cen, w]


@nb.njit(parallel=True, fastmath=True, cache=True)
def _pair_integrals(p0, p1, V, rule6, ruleN, ruleS, near_d2, k, out):
    """For source triangles q in [p0,p1) (observation loop over all p), compute
    normalised double-quadrature integrals:
       out[q-p0, p, 0] = I0 = sum w w G
       out[.., 1:4]    = sum w w r G      (observation point r)
       out[.., 4:7]    = sum w w r' G     (source point r')
       out[.., 7]      = sum w w (r.r') G
    G = exp(-jkR)/(4 pi R)."""
    nt = V.shape[0]
    for qq in nb.prange(p1 - p0):
        q = p0 + qq
        vq0 = V[q, 0]; vq1 = V[q, 1]; vq2 = V[q, 2]
        cq = (vq0 + vq1 + vq2) / 3.0
        for p in range(nt):
            vp0 = V[p, 0]; vp1 = V[p, 1]; vp2 = V[p, 2]
            cp = (vp0 + vp1 + vp2) / 3.0
            d2 = (cp[0] - cq[0]) ** 2 + (cp[1] - cq[1]) ** 2 + (cp[2] - cq[2]) ** 2
            if p == q:
                rule = ruleS
            elif d2 < near_d2:
                rule = ruleN
            else:
                rule = rule6
            I0 = 0j; Irx = 0j; Iry = 0j; Irz = 0j; Ipx = 0j; Ipy = 0j; Ipz = 0j; Irr = 0j
            for a in range(rule6.shape[0]):
                wa = rule6[a, 3]
                rx = rule6[a, 0] * vp0[0] + rule6[a, 1] * vp1[0] + rule6[a, 2] * vp2[0]
                ry = rule6[a, 0] * vp0[1] + rule6[a, 1] * vp1[1] + rule6[a, 2] * vp2[1]
                rz = rule6[a, 0] * vp0[2] + rule6[a, 1] * vp1[2] + rule6[a, 2] * vp2[2]
                for b in range(rule.shape[0]):
                    wb = rule[b, 3]
                    sx = rule[b, 0] * vq0[0] + rule[b, 1] * vq1[0] + rule[b, 2] * vq2[0]
                    sy = rule[b, 0] * vq0[1] + rule[b, 1] * vq1[1] + rule[b, 2] * vq2[1]
                    sz = rule[b, 0] * vq0[2] + rule[b, 1] * vq1[2] + rule[b, 2] * vq2[2]
                    R = np.sqrt((rx - sx) ** 2 + (ry - sy) ** 2 + (rz - sz) ** 2)
                    if R < 1e-15:
                        R = 1e-15
                    g = wa * wb * np.exp(-1j * k * R) / (4.0 * np.pi * R)
                    I0 += g
                    Irx += rx * g; Iry += ry * g; Irz += rz * g
                    Ipx += sx * g; Ipy += sy * g; Ipz += sz * g
                    Irr += (rx * sx + ry * sy + rz * sz) * g
            out[qq, p, 0] = I0
            out[qq, p, 1] = Irx; out[qq, p, 2] = Iry; out[qq, p, 3] = Irz
            out[qq, p, 4] = Ipx; out[qq, p, 5] = Ipy; out[qq, p, 6] = Ipz
            out[qq, p, 7] = Irr


@nb.njit(cache=True)
def _assemble(Z, out, p0, p1, tri_edges, tri_sign, tri_free, length, omega_mu, inv_omega_eps):
    """Accumulate Z from pair integrals for source triangles q in [p0,p1)."""
    nt = out.shape[1]
    for qq in range(p1 - p0):
        q = p0 + qq
        for p in range(nt):
            I0 = out[qq, p, 0]
            Irx = out[qq, p, 1]; Iry = out[qq, p, 2]; Irz = out[qq, p, 3]
            Ipx = out[qq, p, 4]; Ipy = out[qq, p, 5]; Ipz = out[qq, p, 6]
            Irr = out[qq, p, 7]
            for a in range(3):
                m = tri_edges[p, a]
                if m < 0:
                    continue
                sm = tri_sign[p, a]
                pm = tri_free[p, a]
                for b in range(3):
                    n = tri_edges[q, b]
                    if n < 0:
                        continue
                    sn = tri_sign[q, b]
                    pn = tri_free[q, b]
                    Q1 = Irr - (pn[0] * Irx + pn[1] * Iry + pn[2] * Irz) \
                        - (pm[0] * Ipx + pm[1] * Ipy + pm[2] * Ipz) \
                        + (pm[0] * pn[0] + pm[1] * pn[1] + pm[2] * pn[2]) * I0
                    Z[m, n] += sm * sn * length[m] * length[n] * (
                        1j * omega_mu * 0.25 * Q1 - 1j * inv_omega_eps * I0)


class RWGModel:
    def __init__(self, nodes, tris, eps_bg=1.0, near_factor=2.5, chunk=192):
        self.nodes = nodes; self.tris = tris
        self.rwg = rwg_edges(nodes, tris)
        self.eps_bg = eps_bg
        self.chunk = chunk
        nt = len(tris); nE = self.rwg['nE']
        # per-triangle edge tables
        self.tri_edges = -np.ones((nt, 3), np.int64)
        self.tri_sign = np.zeros((nt, 3))
        self.tri_free = np.zeros((nt, 3, 3))
        cnt = np.zeros(nt, int)
        r = self.rwg
        for m in range(nE):
            for t, s, v in ((r['tp'][m], 1.0, r['vp'][m]), (r['tm'][m], -1.0, r['vm'][m])):
                c = cnt[t]
                self.tri_edges[t, c] = m; self.tri_sign[t, c] = s; self.tri_free[t, c] = nodes[v]
                cnt[t] += 1
        self.V = nodes[tris].astype(np.float64)
        emax = np.sqrt(4 * r['area'].max())
        self.near_d2 = (near_factor * emax) ** 2
        self.rule6 = _G6.copy()
        self.ruleN = _subdivision_rule(4)
        self.ruleS = _subdivision_rule(5)
        self.nE = nE; self.nt = nt

    def impedance_matrix(self, f):
        k = 2 * np.pi * f / c0 * np.sqrt(self.eps_bg)
        omega = 2 * np.pi * f
        Z = np.zeros((self.nE, self.nE), np.complex128)
        nt = self.nt
        for p0 in range(0, nt, self.chunk):
            p1 = min(nt, p0 + self.chunk)
            out = np.zeros((p1 - p0, nt, 8), np.complex128)
            _pair_integrals(p0, p1, self.V, self.rule6, self.ruleN, self.ruleS, self.near_d2, k, out)
            _assemble(Z, out, p0, p1, self.tri_edges, self.tri_sign, self.tri_free,
                      self.rwg['length'], omega * mu0, 1.0 / (omega * eps0 * self.eps_bg))
        return 0.5 * (Z + Z.T)

    def current_at_centroids(self, I):
        """Surface current density (A/m) at triangle centroids for coefficient vector I."""
        J = np.zeros((self.nt, 3), np.complex128)
        r = self.rwg; cen = r['cen']
        for m in range(self.nE):
            for t, s, v in ((r['tp'][m], 1.0, r['vp'][m]), (r['tm'][m], -1.0, r['vm'][m])):
                J[t] += I[m] * s * r['length'][m] / (2 * r['area'][t]) * (cen[t] - self.nodes[v])
        return J


def characteristic_modes(Z, nmodes=12, rtol=1e-9):
    """Solve X J = lambda R J with R-subspace filtering. Returns eigenvalues
    (sorted by |lambda|) and R-normalised eigencurrents (columns)."""
    R = 0.5 * (Z.real + Z.real.T); X = 0.5 * (Z.imag + Z.imag.T)
    w, U = eigh(R)
    keep = w > rtol * w.max()
    S = U[:, keep] / np.sqrt(w[keep])
    Xt = S.T @ X @ S
    lam, Vt = eigh(Xt)
    order = np.argsort(np.abs(lam))[:nmodes]
    return lam[order], S @ Vt[:, order]


def modal_significance(lam):
    return 1.0 / np.abs(1 + 1j * lam)


def characteristic_angle(lam):
    return 180.0 - np.degrees(np.arctan(lam))


def track_modes(J_prev, J_new, R_new):
    """Match columns of J_new to J_prev by R-weighted correlation. Returns index array."""
    C = np.abs(J_prev.T @ R_new @ J_new)
    n = min(C.shape)
    idx = -np.ones(J_prev.shape[1], int)
    used = set()
    for _ in range(n):
        i, j = np.unravel_index(np.argmax(C), C.shape)
        idx[i] = j; used.add(j)
        C[i, :] = -1; C[:, j] = -1
    return idx


def _assemble_generic(model, f, k, coef_vec, coef_sca):
    """Low-level: Z = sum s s l l [ j*coef_vec/4 * Q1 - j*coef_sca * I0 ] with wavenumber k."""
    Z = np.zeros((model.nE, model.nE), np.complex128)
    nt = model.nt
    for p0 in range(0, nt, model.chunk):
        p1 = min(nt, p0 + model.chunk)
        out = np.zeros((p1 - p0, nt, 8), np.complex128)
        _pair_integrals(p0, p1, model.V, model.rule6, model.ruleN, model.ruleS, model.near_d2, k, out)
        _assemble(Z, out, p0, p1, model.tri_edges, model.tri_sign, model.tri_free, model.rwg['length'], coef_vec, coef_sca)
    return 0.5 * (Z + Z.T)


def aperture_admittance_matrix(model, f, eps_sides=(1.0, 1.0)):
    """Admittance operator of apertures in an infinite PEC plane (magnetic-current
    MoM, RWG expansion of M). Each side i (relative permittivity eps_i, mu0) is a
    half-space; the image theorem doubles the free-space magnetic-current operator:
        Y_i = 2 [ j w eps_i  int int f.f G_i  -  1/(j w mu0) int int (div f)(div f) G_i ]
    H_tan = -Y M.  Returns Y = sum_i Y_i (siemens)."""
    omega = 2 * np.pi * f
    Y = np.zeros((model.nE, model.nE), np.complex128)
    for er in eps_sides:
        k = omega / c0 * np.sqrt(er)
        Y += 2 * _assemble_generic(model, f, k, omega * eps0 * er, 1.0 / (omega * mu0))
    return Y
