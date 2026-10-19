"""Intensity, safety indices (MI, TIS, TIB), Pennes bioheat and acoustic
radiation force / static displacement.

The 3-D source is built by separating the azimuth-depth intensity map I(x,z)
from the unfocused elevation beam profile E(y,z) of the element height.
The Pennes equation is solved exactly in the spatial-frequency domain for a
homogeneous medium with the thermal properties of the chosen tissue, and
the radiation-force displacement uses the incompressible Kelvin solution
convolved with the body-force distribution.
"""
import numpy as np
from scipy.special import fresnel
from scipy.interpolate import RegularGridInterpolator
from . import params as P


def intensity(p, rho=P.RHO0, c=P.C0_TISSUE):
    """Pulse-average (temporal-peak) intensity in W/m^2 of a CW amplitude."""
    return np.abs(p)**2 / (2 * rho * c)


def mechanical_index(p_r_pa, f=P.F0):
    return (p_r_pa / 1e6) / np.sqrt(f / 1e6)


def derated(p_pa, z_m, f=P.F0, a_db=0.3):
    """IEC derating (0.3 dB/cm/MHz)."""
    return p_pa * 10 ** (-a_db * (f / 1e6) * (z_m * 100) / 20)


def elevation_profile(y, z, H=P.ELEV_HEIGHT, f=P.F0, c=P.C0_TISSUE):
    """|pressure| elevation profile of an unfocused strip of height H,
    normalised to 1 on the y = 0 axis at each depth (Fresnel)."""
    k = 2 * np.pi * f / c
    Y, Z = np.meshgrid(y, z, indexing="ij")
    s = np.sqrt(k / (np.pi * Z))
    S2, C2 = fresnel((H / 2 - Y) * s)
    S1, C1 = fresnel((-H / 2 - Y) * s)
    amp = np.abs((C2 - C1) - 1j * (S2 - S1))
    S0, C0 = fresnel((H / 2) * s[0:1])
    amp0 = np.abs(2 * (C0 - 1j * S0))
    return amp / amp0           # shape (ny, nz)


def build_3d(I_xz, x, z, d=0.3e-3, half_x=24e-3, half_y=24e-3, z_max=64e-3):
    """Interpolate I(x,z) onto a padded 3-D grid and apply the elevation
    profile. Returns grid vectors and I(x,y,z)."""
    gx = np.arange(-half_x, half_x, d)
    gy = np.arange(-half_y, half_y, d)
    gz = np.arange(0, z_max, d)
    interp = RegularGridInterpolator((z, x), I_xz, bounds_error=False,
                                     fill_value=0.0)
    ZZ, XX = np.meshgrid(gz, gx, indexing="ij")
    Ixz = interp(np.stack([ZZ.ravel(), XX.ravel()], -1)).reshape(gz.size,
                                                                  gx.size).T
    E = elevation_profile(gy, np.maximum(gz, 0.5e-3))  # (ny, nz)
    I3 = Ixz[:, None, :] * (E**2)[None, :, :]
    return gx, gy, gz, I3


def absorption_map(gz, f=P.F0):
    a = np.array([P.alpha_np(P.layer_at(z)[5], f) for z in gz])
    return a


def pennes_peak_vs_time(q3, d, tissue="Tumour", times=(np.inf,), at=None):
    """Temperature rise (K) from heat source q3 (W/m^3) on a 3-D grid with
    spacing d, homogeneous Pennes medium. Returns (dT_at_point[t], dT_ss)."""
    k_t, rc, wb = P.THERMAL[tissue]
    D = k_t / rc
    lam_b = wb * P.RHO_C_BLOOD / rc
    qh = np.fft.rfftn(q3)
    kx = 2 * np.pi * np.fft.fftfreq(q3.shape[0], d)
    ky = 2 * np.pi * np.fft.fftfreq(q3.shape[1], d)
    kz = 2 * np.pi * np.fft.rfftfreq(q3.shape[2], d)
    K2 = kx[:, None, None]**2 + ky[None, :, None]**2 + kz[None, None, :]**2
    lam = D * K2 + lam_b
    Tss_h = qh / rc / lam
    Tss = np.fft.irfftn(Tss_h, s=q3.shape)
    if at is None:
        at = np.unravel_index(np.argmax(Tss), Tss.shape)
    # point evaluation of transient via phase factors (exact, cheap)
    ex = np.exp(1j * kx * at[0] * d)[:, None, None]
    ey = np.exp(1j * ky * at[1] * d)[None, :, None]
    ez = np.exp(1j * kz * at[2] * d)[None, None, :]
    basis = ex * ey * ez
    wts = np.ones(kz.size)
    wts[1:] = 2.0
    if q3.shape[2] % 2 == 0:
        wts[-1] = 1.0
    N = q3.size
    out = []
    for t in times:
        fac = 1.0 if np.isinf(t) else (1 - np.exp(-lam * t))
        val = (Tss_h * fac * basis * wts[None, None, :]).sum().real / N
        out.append(val)
    return np.array(out), Tss, at


def kelvin_displacement(Fz3, d, mu=1.0):
    """Axial static displacement u_z (m) of an incompressible, homogeneous
    elastic solid (shear modulus mu) under axial body force Fz (N/m^3).
    Kelvin solution: u_z = F/(8 pi mu r) (1 + dz^2/r^2), FFT convolution."""
    nx, ny, nz = Fz3.shape
    gx = np.fft.fftfreq(2 * nx, 1 / (2 * nx)) * d
    gy = np.fft.fftfreq(2 * ny, 1 / (2 * ny)) * d
    gz = np.fft.fftfreq(2 * nz, 1 / (2 * nz)) * d
    X, Y, Z = np.meshgrid(gx, gy, gz, indexing="ij")
    r = np.sqrt(X**2 + Y**2 + Z**2)
    r_reg = np.where(r == 0, 0.62 * d, r)   # cell-averaged singularity
    G = (1 + np.where(r == 0, 1 / 3, Z**2 / r_reg**2)) / (8 * np.pi * mu * r_reg)
    Fp = np.zeros((2 * nx, 2 * ny, 2 * nz))
    Fp[:nx, :ny, :nz] = Fz3
    u = np.fft.irfftn(np.fft.rfftn(Fp) * np.fft.rfftn(G), s=Fp.shape)
    return u[:nx, :ny, :nz] * d**3


def iec_thermal_indices(W_mW, Ita_focus_mWcm2, f=P.F0):
    """Simplified IEC 62359 indices for an unscanned mode.
    TIS (large aperture) = W * f / 210 mW.MHz;
    TIB (bone at focus)  = sqrt(W * Ita) / 50 (mW, mW/cm^2)."""
    fm = f / 1e6
    TIS = W_mW * fm / 210.0
    TIB = np.sqrt(W_mW * Ita_focus_mWcm2) / 50.0
    return TIS, TIB
