"""Complex non-linear least-squares (CNLS) fitting and Kramers-Kronig validation."""

from __future__ import annotations

import numpy as np
from scipy.optimize import least_squares

from .constants import EPS0

# --------------------------------------------------------------------------- dielectric


DIEL_NAMES = ("Q_ep", "n_ep", "eps_inf", "d_eps", "f_c", "alpha", "sigma")


def dielectric_model(w, K, Q, n, eps_inf, d_eps, f_c, alpha, sigma):
    """EP-CPE in series with a Cole-Cole beta-dispersion plus dc conductivity."""
    tau = 1.0 / (2 * np.pi * f_c)
    eps = eps_inf + d_eps / (1 + (1j * w * tau) ** (1 - alpha)) - 1j * sigma / (w * EPS0)
    return 1.0 / (Q * (1j * w) ** n) + K / (1j * w * EPS0 * eps)


def _pack(p):
    Q, n, eps_inf, d_eps, f_c, alpha, sigma = p
    return np.array([np.log(Q), n, np.log(eps_inf), np.log(d_eps), np.log(f_c), alpha, np.log(sigma)])


def _unpack(x):
    return np.array([np.exp(x[0]), x[1], np.exp(x[2]), np.exp(x[3]), np.exp(x[4]), x[5], np.exp(x[6])])


def fit_dielectric(f, z, K):
    """Fit the 7-parameter model with modulus weighting; returns dict of params + stats."""
    w = 2 * np.pi * f
    zabs = np.abs(z)

    def resid(x):
        zm = dielectric_model(w, K, *_unpack(x))
        r = (z - zm) / zabs
        return np.concatenate([r.real, r.imag])

    i_mid = np.argmin(np.abs(np.angle(z)))
    sigma0 = K / zabs[i_mid]
    lo = np.array([np.log(1e-9), 0.5, np.log(1.0), np.log(1.0), np.log(1e5), 0.0, np.log(1e-3)])
    hi = np.array([np.log(1e-4), 1.0, np.log(500), np.log(5e3), np.log(1e9), 0.6, np.log(10)])
    best = None
    for fc0 in (3e6, 1e7, 3e7):
        for de0 in (50.0, 200.0):
            x0 = _pack([4e-7, 0.8, 70.0, de0, fc0, 0.1, sigma0])
            x0 = np.clip(x0, lo + 1e-6, hi - 1e-6)
            sol = least_squares(resid, x0, bounds=(lo, hi), x_scale="jac", max_nfev=4000)
            if best is None or sol.cost < best.cost:
                best = sol
    p = _unpack(best.x)
    # covariance -> standard errors (in natural units via delta method)
    J = best.jac
    dof = max(1, 2 * f.size - p.size)
    s2 = 2 * best.cost / dof
    try:
        cov = np.linalg.inv(J.T @ J) * s2
        se_x = np.sqrt(np.diag(cov))
    except np.linalg.LinAlgError:
        se_x = np.full(p.size, np.nan)
    deriv = np.array([p[0], 1, p[2], p[3], p[4], 1, p[6]])
    out = dict(zip(DIEL_NAMES, p))
    out.update({f"se_{k}": v for k, v in zip(DIEL_NAMES, se_x * deriv)})
    out["chi2_red"] = s2
    out["z_fit"] = dielectric_model(w, K, *p)
    return out


# --------------------------------------------------------------------------- Randles


RANDLES_NAMES = ("Rs", "Q_dl", "n_dl", "Rct", "sigma_w")


def randles_model(w, Rs, Q, n, Rct, sw):
    zw = sw * (1 - 1j) / np.sqrt(w)
    zf = Rct + zw
    ycpe = Q * (1j * w) ** n
    return Rs + 1.0 / (ycpe + 1.0 / zf)


def fit_randles(f, z):
    w = 2 * np.pi * f
    zabs = np.abs(z)

    def resid(x):
        Rs, lQ, n, lR, ls = x
        zm = randles_model(w, Rs, np.exp(lQ), n, np.exp(lR), np.exp(ls))
        r = (z - zm) / zabs
        return np.concatenate([r.real, r.imag])

    Rs0 = z.real[np.argmax(f)]
    span = z.real.max() - Rs0
    x0 = np.array([Rs0, np.log(1e-6), 0.9, np.log(max(span * 0.3, 1.0)), np.log(100.0)])
    lo = [0, np.log(1e-9), 0.5, np.log(1e-1), np.log(1e-1)]
    hi = [1e4, np.log(1e-3), 1.0, np.log(1e7), np.log(1e6)]
    sol = least_squares(resid, np.clip(x0, np.array(lo) + 1e-6, np.array(hi) - 1e-6),
                        bounds=(lo, hi), x_scale="jac", max_nfev=4000)
    Rs, lQ, n, lR, ls = sol.x
    p = (Rs, np.exp(lQ), n, np.exp(lR), np.exp(ls))
    out = dict(zip(RANDLES_NAMES, p))
    out["z_fit"] = randles_model(w, *p)
    out["chi2_red"] = 2 * sol.cost / max(1, 2 * f.size - 5)
    return out


# --------------------------------------------------------------------------- Lin-KK


def lin_kk(f, z, c=0.85, max_M=80):
    """Linear Kramers-Kronig test (Schonleber et al., 2014) with automatic M selection.

    Model: Z = R0 + jwL + 1/(jwC) + sum_k R_k / (1 + jw tau_k), tau log-spaced.
    Returns (M, mu, residual_real, residual_imag) with residuals relative to |Z|.
    """
    w = 2 * np.pi * f
    zabs = np.abs(z)
    tmin, tmax = 1 / w.max(), 1 / w.min()
    for M in range(3, max_M + 1):
        tau = np.logspace(np.log10(tmin), np.log10(tmax), M)
        basis = 1.0 / (1 + 1j * w[:, None] * tau[None, :])
        A = np.hstack([np.ones((w.size, 1)), basis, (1j * w)[:, None], (1 / (1j * w))[:, None]])
        Aw = A / zabs[:, None]
        bw = z / zabs
        Ar = np.vstack([Aw.real, Aw.imag])
        br = np.concatenate([bw.real, bw.imag])
        scale = np.linalg.norm(Ar, axis=0)
        coef, *_ = np.linalg.lstsq(Ar / scale, br, rcond=None)
        coef = coef / scale
        Rk = coef[1:M + 1]
        pos = Rk[Rk >= 0].sum()
        neg = -Rk[Rk < 0].sum()
        mu = 1 - neg / pos if pos > 0 else 0
        if mu <= c:
            break
    zf = A @ coef
    res = (z - zf) / zabs
    return M, mu, res.real, res.imag
