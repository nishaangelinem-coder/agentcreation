"""Continuous-wave Rayleigh-Sommerfeld field of a linear phased array.

Each element is split into N_SUB point-source strips in azimuth; the
elevation (y) integral over the element height is evaluated in closed form
with the Fresnel approximation, so the field in the azimuth (x-z) plane
accounts for the finite element height without a 3-D grid.

Two propagation models are available:
  * homogeneous: c = 1540 m/s with a single attenuation coefficient;
  * layered: straight-ray path through the skin/fat/muscle/tumour stack with
    per-layer sound speed (phase), attenuation and normal-incidence
    interface transmission coefficients.
"""
import numpy as np
from scipy.special import fresnel
from scipy.signal import windows
from . import params as P

N_SUB = 4


def element_centres(n=P.N_ELEM, pitch=P.PITCH):
    return (np.arange(n) - (n - 1) / 2) * pitch


def apodization(name, n=P.N_ELEM):
    name = name.lower()
    if name in ("rect", "rectangular", "none"):
        return np.ones(n)
    if name == "hann":
        return windows.hann(n + 2)[1:-1]
    if name == "hamming":
        return windows.hamming(n)
    if name.startswith("tukey"):
        return windows.tukey(n, 0.5)
    raise ValueError(name)


def _elev_factor(r, k, H):
    """int_{-H/2}^{H/2} exp(-j k y^2 / 2r) dy  (Fresnel approximation)."""
    u = (H / 2) * np.sqrt(2 * k / (np.pi * r))
    S, C = fresnel(u)
    return np.sqrt(np.pi * r / k) * 2 * (C - 1j * S) / np.sqrt(2)


def _layer_path(z_src, z, x_src, x, layers):
    """Travel time, attenuation (Np) and transmission product along the
    straight ray from (x_src, 0) to (x, z) through horizontal layers."""
    R = np.sqrt((x - x_src)**2 + z**2)
    tt = np.zeros_like(R)
    att = np.zeros_like(R)
    T = np.ones_like(R)
    prevZ = None
    for (name, z0, z1, c, rho, a_db, G) in layers:
        dz = np.clip(np.minimum(z, z1) - z0, 0, None)
        frac = np.where(z > 0, dz / np.maximum(z, 1e-12), 0)
        seg = frac * R
        tt += seg / c
        att += P.alpha_np(a_db) * seg
        Zc = rho * c
        if prevZ is not None:
            crossed = z > z0
            T = np.where(crossed, T * 2 * Zc / (Zc + prevZ), T)
        prevZ = Zc
    return R, tt, att, T


def focal_delays(xf, zf, xe, layered=False, layers=None):
    """Transmit delays (s) that align arrivals at (xf, zf).
    If layered is False the beamformer assumes c = 1540 m/s."""
    if layered:
        _, tt, _, _ = _layer_path(0, np.full_like(xe, zf), xe,
                                  np.full_like(xe, xf), layers)
    else:
        tt = np.sqrt((xf - xe)**2 + zf**2) / P.C0_TISSUE
    return tt.max() - tt


def compute_field(x, z, **kw):
    """Complex pressure on the x-z grid defined by vectors x, z."""
    X, Z = np.meshgrid(x, z)
    return field_at(X.ravel(), Z.ravel(), **kw).reshape(X.shape)


def field_at(Xf, Zf, xf=0.0, zf=P.FOCUS_Z, apod="rect", f=P.F0,
                  pitch=P.PITCH, width=None, n_elem=P.N_ELEM,
                  model="homogeneous", alpha_db=0.54, layered_delays=False,
                  amp_err=None, phase_err=None, delay_quant=None,
                  dead=None, z_offset=None, H=P.ELEV_HEIGHT, chunk=40000):
    """Complex pressure (arbitrary units) at points (Xf, Zf).

    Hardware non-idealities (all optional, per element):
      amp_err    multiplicative amplitude error
      phase_err  additive phase error (rad)
      delay_quant  delay quantisation step (s), e.g. FPGA clock period
      dead       boolean mask of dead elements
      z_offset   element axial displacement (m) from patch curvature
    """
    if width is None:
        width = pitch - P.KERF if pitch - P.KERF > 0 else 0.9 * pitch
    xe = element_centres(n_elem, pitch)
    w_apod = apodization(apod, n_elem)
    if amp_err is not None:
        w_apod = w_apod * amp_err
    if dead is not None:
        w_apod = np.where(dead, 0.0, w_apod)
    layers = P.LAYERS + [P.DEEP]
    tau = focal_delays(xf, zf, xe, layered=layered_delays, layers=layers)
    if delay_quant:
        tau = np.round(tau / delay_quant) * delay_quant
    phase = -2 * np.pi * f * tau  # delayed firing
    if phase_err is not None:
        phase = phase + phase_err
    ze = np.zeros(n_elem) if z_offset is None else z_offset

    sub = (np.arange(N_SUB) - (N_SUB - 1) / 2) * (width / N_SUB)
    xs = (xe[:, None] + sub[None, :]).ravel()
    zs = np.repeat(ze, N_SUB)
    ws = np.repeat(w_apod * np.exp(1j * phase), N_SUB) * (width / N_SUB)

    Xf, Zf = np.asarray(Xf, float).ravel(), np.asarray(Zf, float).ravel()
    out = np.zeros(Xf.size, complex)
    k_ref = 2 * np.pi * f / P.C0_TISSUE
    for i0 in range(0, Xf.size, chunk):
        xx = Xf[i0:i0 + chunk][None, :]
        zz = Zf[i0:i0 + chunk][None, :]
        dzs = zz - zs[:, None]
        dzs = np.maximum(dzs, 1e-6)
        if model == "homogeneous":
            R = np.sqrt((xx - xs[:, None])**2 + dzs**2)
            ph = k_ref * R
            att = P.alpha_np(alpha_db, f) * R
            T = 1.0
        else:
            R, tt, att, T = _layer_path(0, dzs, xs[:, None], xx, layers)
            ph = 2 * np.pi * f * tt
        cos_t = dzs / R
        elev = _elev_factor(R, k_ref, H)
        g = (1j * k_ref / (2 * np.pi)) * np.exp(-1j * ph - att) * T \
            * elev / R * (0.5 * (1 + cos_t))   # obliquity factor
        out[i0:i0 + chunk] = (ws[:, None] * g).sum(0)
    return out


def scale_to_peak(p, target=P.P_TARGET):
    s = target / np.abs(p).max()
    return p * s, s


def fwhm(axis, profile, level=0.5):
    """Full width at `level` x maximum of |profile| (default: pressure
    FWHM, -6 dB). level=1/sqrt(2) gives the intensity FWHM (-3 dB)."""
    prof = np.abs(profile)
    i = int(np.argmax(prof))
    half = prof[i] * level
    l = i
    while l > 0 and prof[l] > half:
        l -= 1
    r = i
    while r < prof.size - 1 and prof[r] > half:
        r += 1
    xl = np.interp(half, [prof[l], prof[l + 1]], [axis[l], axis[l + 1]])
    xr = np.interp(half, [prof[r], prof[r - 1]], [axis[r], axis[r - 1]])
    return xr - xl


def mainlobe_and_psl(axis, profile):
    """-6 dB main-lobe width and peak side-lobe level (dB) of a profile."""
    prof = np.abs(profile)
    prof = prof / prof.max()
    i = int(np.argmax(prof))
    l = i
    while l > 0 and prof[l - 1] < prof[l]:
        l -= 1
    r = i
    while r < prof.size - 1 and prof[r + 1] < prof[r]:
        r += 1
    side = np.concatenate([prof[:l], prof[r + 1:]])
    psl = 20 * np.log10(side.max()) if side.size else -np.inf
    return fwhm(axis, prof), psl


def arc_pattern(theta, R=P.FOCUS_Z, **kw):
    """Field sampled on a circular arc of radius R (polar beam pattern)."""
    xs = R * np.sin(theta)
    zs = R * np.cos(theta)
    return field_at(xs, zs, **kw)
