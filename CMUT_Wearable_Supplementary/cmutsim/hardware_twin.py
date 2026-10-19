"""Hardware digital twin: Monte Carlo emulation of what a bench
characterisation of the fabricated patch would be expected to report.

EVERYTHING PRODUCED HERE IS SIMULATED. The outputs are predictions of
measurement outcomes (distributions with realistic error sources), not
experimental data, and must be labelled "emulated" in any publication.

Error sources modelled
  Fabrication : gap, plate-thickness and radius spread (wafer level)
  Electronics : FPGA delay quantisation, per-channel gain/phase mismatch,
                dead elements, common bias set below the weakest element
  Wearable    : conformal curvature of the patch (uncorrected / corrected)
  Measurement : needle-hydrophone spatial averaging, calibration
                uncertainty, electronic noise, scan-stage positioning error,
                thermocouple noise and viscous-heating artefact
"""
import numpy as np
from . import params as P
from . import cmut_cell as CC
from . import field as F

# ---- default non-ideality budget (1-sigma unless stated) ---------------
HW = dict(
    sigma_gap=0.03,          # relative, LPCVD sacrificial-layer thickness
    sigma_t=0.03,            # relative, plate thickness
    sigma_a=0.01,            # relative, lithographic radius
    bias_frac=0.80,          # Vdc = bias_frac * min(Vc) over the array
    sigma_ch_gain_db=0.5,    # TX channel gain mismatch
    sigma_ch_phase_deg=3.0,  # TX channel phase mismatch
    fpga_clk=200e6,          # delay resolution 5 ns
    p_dead=0.015,            # probability an element is dead (~1 per array)
    hyd_diam=0.2e-3,         # needle hydrophone active diameter
    hyd_cal_unc=0.08,        # 1-sigma calibration uncertainty (relative)
    hyd_noise_pa=150.0,      # RMS noise after averaging (Pa)
    stage_sigma=10e-6,       # positioning repeatability (m)
    tc_noise_k=0.010,        # thermocouple RMS noise (K) after averaging
    tc_viscous_frac=0.3,     # viscous-heating artefact / true rise (early time)
)


def fabricate_array(rng, n=P.N_ELEM, hw=HW):
    """Per-element geometric spread (cells in an element share a die region,
    so element-to-element variation dominates)."""
    g = P.GAP * (1 + hw["sigma_gap"] * rng.standard_normal(n))
    t = P.PLATE_T * (1 + hw["sigma_t"] * rng.standard_normal(n))
    a = P.CELL_RADIUS * (1 + hw["sigma_a"] * rng.standard_normal(n))
    Vc = np.array([CC.collapse_voltage(ai, ti, gi) for ai, ti, gi in zip(a, t, g)])
    return dict(gap=g, t=t, a=a, Vc=Vc)


def element_response(arr, hw=HW, f=P.F0):
    """Relative transmit sensitivity and phase of each element at f under a
    common bias Vdc = bias_frac * min(Vc)."""
    Vdc = hw["bias_frac"] * arr["Vc"].min()
    sens = np.array([CC.surface_pressure_per_volt(Vdc, a, g)
                     for a, g in zip(arr["a"], arr["gap"])])
    # resonance spread -> phase of a 2nd-order response with Q_LOADED
    fr = np.array([CC.resonance_immersed(a, t, V=Vdc, g0=g)
                   for a, t, g in zip(arr["a"], arr["t"], arr["gap"])])
    fr_nom = CC.resonance_immersed(P.CELL_RADIUS, P.PLATE_T,
                                   V=hw["bias_frac"] * CC.collapse_voltage())
    def ph(fr_):
        x = f / fr_
        return np.arctan2(x / P.Q_LOADED, 1 - x**2)
    dphi = ph(fr) - ph(fr_nom)
    return sens / sens.mean(), dphi, Vdc


def realise_nonidealities(rng, hw=HW, curvature_R=None, corrected=False):
    arr = fabricate_array(rng, hw=hw)
    s, dphi, Vdc = element_response(arr, hw)
    g_ch = 10 ** (hw["sigma_ch_gain_db"] * rng.standard_normal(P.N_ELEM) / 20)
    ph_ch = np.deg2rad(hw["sigma_ch_phase_deg"]) * rng.standard_normal(P.N_ELEM)
    dead = rng.random(P.N_ELEM) < hw["p_dead"]
    zoff = None
    if curvature_R:
        xe = F.element_centres()
        # convex body: element centres fall away from the flat-patch plane
        zoff = curvature_R - np.sqrt(curvature_R**2 - xe**2)
        if corrected:
            # shape-sensed delay correction (strain-gauge estimate, 50 um RMS):
            # elements displaced toward the focus have a shorter path and
            # must fire later by the estimated path difference
            zest = zoff + 50e-6 * rng.standard_normal(P.N_ELEM)
            k = 2 * np.pi * P.F0 / P.C0_TISSUE
            dR = np.hypot(xe, P.FOCUS_Z) - np.hypot(xe, P.FOCUS_Z - zest)
            ph_ch = ph_ch - k * dR
    return dict(amp_err=s * g_ch, phase_err=dphi + ph_ch, dead=dead,
                delay_quant=1 / hw["fpga_clk"], z_offset=zoff, Vdc=Vdc,
                Vc=arr["Vc"])


def hydrophone_scan(x_true, p_true_fn, rng, hw=HW, n_avg_pts=7):
    """Emulate one raster line measured with a needle hydrophone.
    p_true_fn(x) returns |p| (Pa) along the line."""
    pos = x_true + hw["stage_sigma"] * rng.standard_normal(x_true.size)
    off = np.linspace(-hw["hyd_diam"] / 2, hw["hyd_diam"] / 2, n_avg_pts)
    w = np.sqrt(np.clip(1 - (2 * off / hw["hyd_diam"])**2, 0, None))
    w /= w.sum()
    vals = sum(wi * np.abs(p_true_fn(pos + oi)) for wi, oi in zip(w, off))
    cal = 1 + hw["hyd_cal_unc"] * rng.standard_normal()
    noisy = vals * cal + hw["hyd_noise_pa"] * rng.standard_normal(vals.size)
    return np.abs(noisy)


def thermocouple_reading(dT_true, rng, n, hw=HW):
    art = hw["tc_viscous_frac"] * dT_true
    return dT_true + art + hw["tc_noise_k"] * rng.standard_normal(n), art
