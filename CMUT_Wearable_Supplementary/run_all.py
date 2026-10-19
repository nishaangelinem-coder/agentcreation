#!/usr/bin/env python3
"""Reproduce every in-silico experiment of the manuscript and run the
hardware digital twin. Usage:  python run_all.py [--quick]

Outputs
  results/data/                  raw numerical data (CSV, NPZ)
  results/tables/                tables (CSV + LaTeX)
  results/figures/               figures (PNG 600 dpi + vector PDF)
  results/emulated_measurements/ digital-twin (SIMULATED) bench outputs
  results/summary.json           every headline number
"""
import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from cmutsim import params as P
from cmutsim import cmut_cell as CC
from cmutsim import field as F
from cmutsim import bioeffects as B
from cmutsim import control as CL
from cmutsim import hardware_twin as HWT

QUICK = "--quick" in sys.argv
SEED = 20261019
ROOT = Path(__file__).resolve().parent / "results"
D_DATA, D_TAB, D_FIG, D_EMU = (ROOT / s for s in
                               ("data", "tables", "figures",
                                "emulated_measurements"))
for d in (D_DATA, D_TAB, D_FIG, D_EMU):
    d.mkdir(parents=True, exist_ok=True)

SUMMARY = {}
C = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300",
     "#6250d6", "#e34948"]
plt.rcParams.update({
    "font.family": "DejaVu Sans", "font.size": 8, "axes.titlesize": 8.5,
    "axes.labelsize": 8, "legend.fontsize": 7, "xtick.labelsize": 7,
    "ytick.labelsize": 7, "axes.spines.top": False,
    "axes.spines.right": False, "axes.grid": True, "grid.alpha": 0.25,
    "grid.linewidth": 0.5, "lines.linewidth": 1.4, "savefig.dpi": 600,
    "axes.prop_cycle": plt.cycler(color=C), "figure.dpi": 100,
})
EMU_TAG = "EMULATED (digital twin) - not experimental data"
mm = 1e3


def save_fig(fig, name, emulated=False):
    if emulated:
        fig.text(0.995, 0.005, EMU_TAG, ha="right", va="bottom", fontsize=6,
                 color="#b00020", style="italic")
    d = D_EMU if emulated else D_FIG
    fig.savefig(d / f"{name}.png", bbox_inches="tight")
    fig.savefig(d / f"{name}.pdf", bbox_inches="tight")
    plt.close(fig)


def save_table(df, name, caption, emulated=False):
    d = D_EMU if emulated else D_TAB
    df.to_csv(d / f"{name}.csv", index=False)
    tex = df.to_latex(index=False, float_format="%.4g", caption=caption,
                      label=f"tab:{name}", escape=True)
    (d / f"{name}.tex").write_text(tex)


def log(msg):
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


# grids ------------------------------------------------------------------
XM = np.arange(-18, 18.0001, 0.1 if not QUICK else 0.25) * 1e-3
ZM = np.arange(0.5, 32.0001, 0.1 if not QUICK else 0.25) * 1e-3
XL = np.arange(-6, 6.0001, 0.005) * 1e-3
ZL = np.arange(2, 34.0001, 0.02) * 1e-3


def focal_profiles(**kw):
    """Fine axial (through xf) and lateral (through the axial peak) lines."""
    xf = kw.get("xf", 0.0)
    zf = kw.get("zf", P.FOCUS_Z)
    # axial line along the steering direction
    th = np.arctan2(xf, zf)
    rr = ZL
    pa = F.field_at(rr * np.sin(th), rr * np.cos(th), **kw)
    ip = int(np.argmax(np.abs(pa)))
    zpk = rr[ip] * np.cos(th)
    xpk = rr[ip] * np.sin(th)
    pl = F.field_at(xpk + XL, np.full_like(XL, zpk), **kw)
    return dict(r=rr, pa=pa, x=xpk + XL, pl=pl, zpk=zpk, xpk=xpk)


# ======================================================================
def e01_cmut():
    log("E01 CMUT cell and element")
    Vc = CC.collapse_voltage()
    C0 = CC.clamped_capacitance()
    fvac = CC.resonance_vacuum()
    fimm0 = CC.resonance_immersed()
    fimmb = CC.resonance_immersed(V=0.8 * Vc)
    a15 = CC.radius_for_immersed_resonance(P.F0)
    f = np.linspace(0.2e6, 3.2e6, 3001)
    Z, H, bvd = CC.bvd_element(f)
    fbw, fc = CC.fractional_bandwidth(f, H, -6)
    fbw3, _ = CC.fractional_bandwidth(f, H, -3)
    ppv = CC.surface_pressure_per_volt(0.8 * Vc)
    pd.DataFrame(dict(f_Hz=f, Zabs_ohm=np.abs(Z), Zphase_deg=np.angle(Z, 1),
                      Htx_norm=H)).to_csv(D_DATA / "E01_element_response.csv",
                                          index=False)
    radii = np.linspace(20e-6, 50e-6, 61)
    sweep = pd.DataFrame(dict(
        radius_um=radii * 1e6,
        Vc_V=[CC.collapse_voltage(a) for a in radii],
        f_vac_MHz=[CC.resonance_vacuum(a) / 1e6 for a in radii],
        f_imm_0V_MHz=[CC.resonance_immersed(a) / 1e6 for a in radii],
        f_imm_08Vc_MHz=[CC.resonance_immersed(
            a, V=0.8 * CC.collapse_voltage(a)) / 1e6 for a in radii]))
    sweep.to_csv(D_DATA / "E01_radius_sweep.csv", index=False)

    tab = pd.DataFrame([
        ["Collapse voltage Vc (V)", 46.5, Vc],
        ["Clamped capacitance per cell C0 (pF)", 0.125, C0 * 1e12],
        ["Vacuum resonance (MHz)", "n.r.", fvac / 1e6],
        ["Immersed resonance, 0 V (MHz)", "n.r.", fimm0 / 1e6],
        ["Immersed resonance, 0.8 Vc (MHz)", "1.5 (implied)", fimmb / 1e6],
        ["Radius for 1.5 MHz immersed at 0.8 Vc (um)", 30.0, a15 * 1e6],
        ["Cells per element (assumed 6 mm elevation)", "n.r.", bvd["n_cells"]],
        ["Element C0 (pF)", "n.r.", bvd["C0"] * 1e12],
        ["-6 dB fractional bandwidth (%)", 51.8, fbw * 100],
        ["-3 dB fractional bandwidth (%)", "n.r.", fbw3 * 100],
        ["Surface pressure per volt at 0.8 Vc (Pa/V)", "n.r.", ppv],
    ], columns=["Quantity", "Manuscript", "Reproduced"])
    save_table(tab, "T_S1_cmut_cell", "CMUT cell and element model.")

    fig, ax = plt.subplots(1, 3, figsize=(7.16, 2.2))
    ax[0].semilogy(f / 1e6, np.abs(Z), color=C[0])
    ax[0].axvline(P.F0 / 1e6, ls="--", color="0.5", lw=0.8)
    ax[0].set(xlabel="Frequency (MHz)", ylabel="|Z| (ohm)",
              title="(a) Element input impedance")
    ax[1].plot(f / 1e6, 20 * np.log10(H), color=C[2])
    ax[1].axhline(-6, ls=":", color="0.5", lw=0.8)
    ax[1].set(xlabel="Frequency (MHz)", ylabel="Norm. output (dB)",
              ylim=(-30, 2), title=f"(b) TX response, -6 dB FBW {fbw*100:.1f}%")
    ax[2].plot(sweep.radius_um, sweep.f_imm_0V_MHz, label="immersed, 0 V")
    ax[2].plot(sweep.radius_um, sweep.f_imm_08Vc_MHz, label="immersed, 0.8 Vc")
    ax[2].axhline(1.5, ls="--", color="0.4", lw=0.8)
    ax[2].axvline(30, ls=":", color=C[7], lw=0.8)
    ax[2].annotate("design\na = 30 um", (30, 4), fontsize=6, color=C[7])
    ax[2].set(xlabel="Plate radius (um)", ylabel="Resonance (MHz)",
              title="(c) Cell resonance vs radius", ylim=(0, 6))
    ax[2].legend(frameon=False)
    fig.tight_layout()
    save_fig(fig, "Fig_S1_cmut_element")
    SUMMARY["cmut"] = dict(Vc=Vc, C0_pF=C0 * 1e12, f_imm_bias_MHz=fimmb / 1e6,
                           radius_for_1p5MHz_um=a15 * 1e6, fbw6=fbw,
                           Pa_per_V=ppv, n_cells=bvd["n_cells"],
                           elem_C0_pF=bvd["C0"] * 1e12)


# ======================================================================
def e02_focus():
    log("E02 focused field (homogeneous)")
    p = F.compute_field(XM, ZM)
    p, s = F.scale_to_peak(p)
    prof = focal_profiles()
    s_l = P.P_TARGET / np.abs(prof["pa"]).max()
    pa, pl = prof["pa"] * s_l, prof["pl"] * s_l
    lat6 = F.fwhm(prof["x"], pl)
    ax6 = F.fwhm(prof["r"], pa)
    lat3 = F.fwhm(prof["x"], pl, 2**-0.5)
    ax3 = F.fwhm(prof["r"], pa, 2**-0.5)
    Isppa = B.intensity(P.P_TARGET)
    # source surface pressure & acoustic power
    p_s = s_l * F.apodization("rect")             # Pa per element
    W = (p_s**2 / (2 * P.RHO0 * P.C0_TISSUE) * P.WIDTH * P.ELEV_HEIGHT).sum()
    Vac = p_s.max() / SUMMARY["cmut"]["Pa_per_V"]
    np.savez_compressed(D_DATA / "E02_field_homogeneous.npz", x=XM, z=ZM,
                        p_abs=np.abs(p))
    pd.DataFrame(dict(depth_mm=prof["r"] * mm, p_kPa=np.abs(pa) / 1e3)).to_csv(
        D_DATA / "E02_axial_profile.csv", index=False)
    pd.DataFrame(dict(x_mm=prof["x"] * mm, p_kPa=np.abs(pl) / 1e3)).to_csv(
        D_DATA / "E02_lateral_profile.csv", index=False)
    tab = pd.DataFrame([
        ["Focal depth of peak (mm)", 20.0, prof["zpk"] * mm],
        ["Peak focal pressure (kPa)", 12.0, P.P_TARGET / 1e3],
        ["Lateral FWHM, pressure -6 dB (mm)", 0.69, lat6 * mm],
        ["Lateral FWHM, intensity -3 dB (mm)", "-", lat3 * mm],
        ["Axial FWHM, pressure -6 dB (mm)", 3.74, ax6 * mm],
        ["Axial FWHM, intensity -3 dB (mm)", "-", ax3 * mm],
        ["Isppa (mW/cm^2)", 4.45, Isppa / 10],
        ["Ispta at 20 % duty (mW/cm^2)", 0.891, 0.2 * Isppa / 10],
        ["Required surface pressure (kPa)", "n.r.", p_s.max() / 1e3],
        ["Required AC drive at 0.8 Vc bias (V)", "n.r.", Vac],
        ["Acoustic output power, CW (mW)", "n.r.", W * 1e3],
    ], columns=["Quantity", "Manuscript", "Reproduced"])
    save_table(tab, "T_S2_focus", "Focused field at 20 mm (homogeneous).")

    fig = plt.figure(figsize=(7.16, 2.6))
    g = fig.add_gridspec(2, 2, width_ratios=[1.5, 1])
    a0 = fig.add_subplot(g[:, 0])
    im = a0.imshow(np.abs(p) / 1e3, extent=[XM[0]*mm, XM[-1]*mm, ZM[-1]*mm,
                   ZM[0]*mm], aspect="auto", cmap="inferno")
    fig.colorbar(im, ax=a0, label="|p| (kPa)")
    a0.set(xlabel="Azimuth (mm)", ylabel="Depth (mm)",
           title="(a) Focal pressure field")
    a1 = fig.add_subplot(g[0, 1])
    a1.plot(prof["x"] * mm, np.abs(pl) / 1e3)
    a1.set(xlim=(-4, 4), xlabel="Azimuth (mm)", ylabel="|p| (kPa)",
           title=f"(b) Lateral: FWHM {lat6*mm:.2f} mm (-6 dB), "
                 f"{lat3*mm:.2f} mm (-3 dB)")
    a2 = fig.add_subplot(g[1, 1])
    a2.plot(prof["r"] * mm, np.abs(pa) / 1e3, color=C[2])
    a2.set(xlabel="Depth (mm)", ylabel="|p| (kPa)",
           title=f"(c) Axial: FWHM {ax6*mm:.2f} mm (-6 dB), "
                 f"{ax3*mm:.2f} mm (-3 dB)")
    fig.tight_layout()
    save_fig(fig, "Fig_S2_focal_field")
    SUMMARY["focus"] = dict(z_peak_mm=prof["zpk"]*mm, lat6_mm=lat6*mm,
                            lat3_mm=lat3*mm, ax6_mm=ax6*mm, ax3_mm=ax3*mm,
                            Isppa_mWcm2=Isppa/10, Ispta20_mWcm2=0.2*Isppa/10,
                            p_surface_kPa=p_s.max()/1e3, Vac_V=Vac,
                            W_cw_mW=W*1e3, drive_scale=s_l)
    return p


# ======================================================================
def e03_steering():
    log("E03 beam steering")
    rows, maps = [], []
    s0 = SUMMARY["focus"]["drive_scale"]
    for ang in (-20, -10, 0, 10, 20):
        th = np.deg2rad(ang)
        kw = dict(xf=P.FOCUS_Z*np.sin(th), zf=P.FOCUS_Z*np.cos(th))
        pr = focal_profiles(**kw)
        pk = np.abs(pr["pa"]).max() * s0
        # width perpendicular to the beam axis
        u = XL
        pp = F.field_at(pr["xpk"] + u*np.cos(th), pr["zpk"] - u*np.sin(th), **kw)
        rows.append([ang, pk/1e3, 20*np.log10(pk/P.P_TARGET),
                     F.fwhm(pr["x"], pr["pl"])*mm, F.fwhm(u, pp)*mm,
                     F.fwhm(u, pp, 2**-0.5)*mm,
                     np.hypot(pr["xpk"]-kw["xf"], pr["zpk"]-kw["zf"])*mm])
        xs = np.arange(-18, 18.001, 0.2)*1e-3
        zs = np.arange(0.5, 32.001, 0.2)*1e-3
        maps.append(np.abs(F.compute_field(xs, zs, **kw))*s0)
    tab = pd.DataFrame(rows, columns=[
        "Steer (deg)", "Peak (kPa)", "Peak change (dB)",
        "Lateral FWHM along x (mm)", "FWHM perp. to beam -6 dB (mm)",
        "FWHM perp. to beam -3 dB (mm)", "Focal position error (mm)"])
    save_table(tab, "T_S3_steering", "Beam-steering performance at 20 mm range.")
    fig, axs = plt.subplots(1, 5, figsize=(7.16, 1.9), sharey=True)
    vmax = max(m.max() for m in maps)/1e3
    for a, m, ang in zip(axs, maps, (-20, -10, 0, 10, 20)):
        a.imshow(m/1e3, extent=[-18, 18, 32, 0.5], aspect="auto",
                 cmap="inferno", vmax=vmax)
        a.set(title=f"{ang:+d} deg", xlabel="Az (mm)")
        a.grid(False)
    axs[0].set_ylabel("Depth (mm)")
    fig.tight_layout()
    save_fig(fig, "Fig_S3_steering")
    SUMMARY["steering"] = tab.to_dict("list")


# ======================================================================
def e04_grating():
    log("E04 grating lobes vs pitch")
    theta = np.deg2rad(np.linspace(-89, 89, 1781))
    rows = []
    fig, ax = plt.subplots(figsize=(3.5, 2.4))
    for i, pl in enumerate((0.5, 0.75, 1.0, 1.5)):
        pitch = pl * P.LAM
        kw = dict(xf=P.FOCUS_Z*np.sin(np.deg2rad(20)),
                  zf=P.FOCUS_Z*np.cos(np.deg2rad(20)), pitch=pitch,
                  width=pitch - P.KERF)
        v = np.abs(F.arc_pattern(theta, **kw))
        vdb = 20*np.log10(v/v.max())
        sg = np.sin(np.deg2rad(20)) - 1/pl
        gl_ang = np.rad2deg(np.arcsin(sg)) if abs(sg) <= 1 else np.nan
        ml = np.abs(np.rad2deg(theta) - 20) < 6
        if np.isnan(gl_ang):
            gl = np.nan
        else:
            win = np.abs(np.rad2deg(theta) - gl_ang) < 8
            gl = vdb[win].max()
        psl = vdb[~ml].max()
        rows.append([pl, pitch*mm, gl_ang, gl, psl])
        ax.plot(np.rad2deg(theta), vdb, label=f"{pl:.2f} lambda", color=C[i])
    ax.set(xlabel="Angle on 20 mm arc (deg)", ylabel="Norm. |p| (dB)",
           ylim=(-50, 2), title="Arc beam pattern, 20 deg steer")
    ax.legend(frameon=False, title="pitch", ncol=2)
    fig.tight_layout()
    save_fig(fig, "Fig_S4_grating_lobes")
    tab = pd.DataFrame(rows, columns=[
        "Pitch (lambda)", "Pitch (mm)", "Predicted GL angle (deg)",
        "Grating-lobe level (dB)", "Peak lobe outside main lobe (dB)"])
    save_table(tab, "T_S4_grating", "Grating-lobe level versus pitch at 20 deg steer.")
    SUMMARY["grating"] = tab.to_dict("list")


# ======================================================================
def e05_apod():
    log("E05 apodization")
    rows = []
    fig, ax = plt.subplots(figsize=(3.5, 2.4))
    for i, w in enumerate(("rect", "hann", "hamming", "tukey")):
        pr = focal_profiles(apod=w)
        lw, psl = F.mainlobe_and_psl(pr["x"], pr["pl"])
        lw3 = F.fwhm(pr["x"], pr["pl"], 2**-0.5)
        gain = np.abs(pr["pa"]).max()/np.abs(focal_profiles()["pa"]).max()
        rows.append([w.capitalize(), lw*mm, lw3*mm, psl, 20*np.log10(gain)])
        v = np.abs(pr["pl"])
        ax.plot(pr["x"]*mm, 20*np.log10(v/v.max()), label=w.capitalize())
    ax.set(xlim=(-6, 6), ylim=(-60, 2), xlabel="Azimuth (mm)",
           ylabel="Norm. |p| (dB)", title="Apodization at 20 mm focus")
    ax.legend(frameon=False)
    fig.tight_layout()
    save_fig(fig, "Fig_S5_apodization")
    tab = pd.DataFrame(rows, columns=["Window", "Main lobe -6 dB (mm)",
                                      "Main lobe -3 dB (mm)",
                                      "Peak side lobe (dB)",
                                      "Focal gain vs rect (dB)"])
    save_table(tab, "T_S5_apodization", "Apodization trade-off.")
    SUMMARY["apod"] = tab.to_dict("list")


# ======================================================================
def e06_layered():
    log("E06 layered tissue")
    out = {}
    fig, axs = plt.subplots(1, 3, figsize=(7.16, 2.4))
    for j, (lab, ld) in enumerate((("homogeneous-delay", False),
                                   ("layer-corrected", True))):
        kw = dict(model="layered", layered_delays=ld)
        pr = focal_profiles(**kw)
        s = P.P_TARGET/np.abs(pr["pa"]).max()
        ins = 20*np.log10(s/SUMMARY["focus"]["drive_scale"])
        out[lab] = dict(z_peak_mm=pr["zpk"]*mm,
                        lat6_mm=F.fwhm(pr["x"], pr["pl"])*mm,
                        ax6_mm=F.fwhm(pr["r"], pr["pa"])*mm,
                        drive_increase_dB=ins, scale=s,
                        gain_vs_corrected=None)
        if ld:
            pmap = F.compute_field(XM, ZM, **kw)*s
            np.savez_compressed(D_DATA / "E06_field_layered.npz", x=XM, z=ZM,
                                p_abs=np.abs(pmap))
            out["map"] = pmap
            pd.DataFrame(dict(depth_mm=pr["r"]*mm,
                              p_kPa=np.abs(pr["pa"])*s/1e3)).to_csv(
                D_DATA / "E06_axial_layered.csv", index=False)
        axs[2].plot(pr["r"]*mm, np.abs(pr["pa"])*s/1e3, label=lab,
                    color=C[j])
    # skin-surface (near-field) peak pressure
    pm = np.abs(out.pop("map"))
    out["skin_peak_kPa"] = pm[ZM < 2e-3].max()/1e3
    out["fat_peak_kPa"] = pm[(ZM >= 2e-3) & (ZM < 10e-3)].max()/1e3
    im = axs[0].imshow(pm/1e3, extent=[XM[0]*mm, XM[-1]*mm, ZM[-1]*mm,
                       ZM[0]*mm], aspect="auto", cmap="inferno")
    for L in P.LAYERS:
        axs[0].axhline(L[2]*mm, color="w", lw=0.5, ls=":")
        axs[0].text(-17, L[1]*mm+1.5, L[0], color="w", fontsize=6)
    axs[0].set(xlabel="Azimuth (mm)", ylabel="Depth (mm)",
               title="(a) Layered model, corrected delays")
    axs[0].grid(False)
    fig.colorbar(im, ax=axs[0], label="|p| (kPa)")
    z = np.linspace(0, 28e-3, 500)
    axs[1].step(z*mm, [P.layer_at(zz)[3] for zz in z], where="post")
    axs[1].axhline(1540, ls="--", color="0.5", lw=0.8)
    axs[1].set(xlabel="Depth (mm)", ylabel="c (m/s)",
               title="(b) Sound-speed profile")
    for L in P.LAYERS:
        axs[2].axvspan(L[1]*mm, L[2]*mm, alpha=0.06, color="k")
    axs[2].axvline(20, ls="--", color=C[7], lw=0.8)
    axs[2].set(xlabel="Depth (mm)", ylabel="|p| (kPa)",
               title="(c) On-axis pressure, 12 kPa at peak")
    axs[2].legend(frameon=False)
    fig.tight_layout()
    save_fig(fig, "Fig_S6_layered_tissue")
    tab = pd.DataFrame([[k, v["z_peak_mm"], v["lat6_mm"], v["ax6_mm"],
                         v["drive_increase_dB"]] for k, v in out.items()
                        if isinstance(v, dict)],
                       columns=["Delay law", "Peak depth (mm)",
                                "Lateral FWHM (mm)", "Axial FWHM (mm)",
                                "Drive increase vs homogeneous (dB)"])
    save_table(tab, "T_S6_layered", "Effect of the tissue stack on the focus.")
    SUMMARY["layered"] = {k: (v if not isinstance(v, dict) else
                              {kk: vv for kk, vv in v.items()
                               if kk != "gain_vs_corrected"})
                          for k, v in out.items()}
    return pm


# ======================================================================
def e07_thermal(pmap):
    log("E07 Pennes bioheat (3-D spectral)")
    I_xz = B.intensity(pmap)
    d = 0.3e-3 if not QUICK else 0.5e-3
    gx, gy, gz, I3 = B.build_3d(I_xz, XM, ZM, d=d)
    alpha = B.absorption_map(gz)
    q100 = 2*alpha[None, None, :]*I3          # W/m^3 at 100 % duty
    times = np.concatenate([np.linspace(0, 240, 49)[1:], [np.inf]])
    rows, curves = [], {}
    for tissue in ("Tumour", "Fat"):
        Tt, Tss, at = B.pennes_peak_vs_time(q100, d, tissue, times)
        for duty in (0.05, 0.2, 0.5, 1.0):
            dT = Tt[-1]*duty
            rows.append([tissue, duty*100, dT*1e3, Tt[np.argmin(
                np.abs(times-240))]*duty*1e3,
                gz[at[2]]*mm])
            curves[(tissue, duty)] = Tt[:-1]*duty
        if tissue == "Tumour":
            Tmap = Tss[:, np.argmin(np.abs(gy)), :]*0.2
    tab = pd.DataFrame(rows, columns=["Thermal props", "Duty (%)",
                                      "Steady dT (mK)", "dT at 4 min (mK)",
                                      "Hot-spot depth (mm)"])
    # safety indices
    W = SUMMARY["focus"]["W_cw_mW"]
    ins = 10**(SUMMARY["layered"]["layer-corrected"]["drive_increase_dB"]/10)
    sfty = []
    for duty in (0.05, 0.2, 0.5, 1.0):
        Ita = B.intensity(P.P_TARGET)/10*duty
        TIS, TIB = B.iec_thermal_indices(W*ins*duty, Ita)
        sfty.append([duty*100, Ita, W*ins*duty, TIS, TIB])
    st = pd.DataFrame(sfty, columns=["Duty (%)", "Ispta (mW/cm^2)",
                                     "Acoustic power (mW)", "TIS", "TIB"])
    tab = tab.merge(st[["Duty (%)", "Ispta (mW/cm^2)", "TIS", "TIB"]],
                    on="Duty (%)")
    save_table(tab, "T_S7_thermal", "Pennes bioheat and IEC thermal indices.")
    tt = times[:-1]
    pd.DataFrame({f"{k[0]}_{int(k[1]*100)}pct_mK": v*1e3
                  for k, v in curves.items()} | {"t_s": tt}).to_csv(
        D_DATA / "E07_heating_transients.csv", index=False)
    fig, ax = plt.subplots(1, 2, figsize=(7.16, 2.4))
    im = ax[0].imshow(Tmap.T*1e3, extent=[gx[0]*mm, gx[-1]*mm, gz[-1]*mm,
                      gz[0]*mm], aspect="auto", cmap="magma")
    ax[0].set(xlim=(-12, 12), ylim=(30, 0), xlabel="Azimuth (mm)",
              ylabel="Depth (mm)", title="(a) Steady dT, 20 % duty (y = 0)")
    ax[0].grid(False)
    fig.colorbar(im, ax=ax[0], label="dT (mK)")
    for i, duty in enumerate((0.05, 0.2, 0.5, 1.0)):
        ax[1].plot(tt/60, curves[("Tumour", duty)]*1e3, color=C[i],
                   label=f"{int(duty*100)} % (tumour)")
        ax[1].plot(tt/60, curves[("Fat", duty)]*1e3, color=C[i], ls="--",
                   lw=0.9)
    ax[1].set(xlabel="Time (min)", ylabel="Peak dT (mK)",
              title="(b) Heating transients (dashed: fat properties)")
    ax[1].legend(frameon=False)
    fig.tight_layout()
    save_fig(fig, "Fig_S7_thermal")
    SUMMARY["thermal"] = tab.to_dict("list")
    return dict(I3=I3, gx=gx, gy=gy, gz=gz, d=d, alpha=alpha)


# ======================================================================
def e08_mi():
    log("E08 mechanical index")
    rows = []
    for pk in (10e3, 12e3, 15e3, 30e3, 50e3):
        MI = B.mechanical_index(pk)
        rows.append([pk/1e3, MI, B.intensity(pk)/10, 1.9/MI, 0.7/MI])
    tab = pd.DataFrame(rows, columns=["p (kPa)", "MI", "Isppa (mW/cm^2)",
                                      "Margin to FDA MI 1.9 (x)",
                                      "Margin to MI 0.7 nucleation (x)"])
    save_table(tab, "T_S8_mi", "Mechanical index and cavitation margin.")
    pp = np.linspace(1e3, 600e3, 300)
    ff = np.linspace(0.5e6, 3e6, 200)
    PP, FF = np.meshgrid(pp, ff)
    MI = B.mechanical_index(PP, FF)
    fig, ax = plt.subplots(figsize=(3.5, 2.5))
    cs = ax.contourf(PP/1e3, FF/1e6, MI, 20, cmap="viridis")
    fig.colorbar(cs, ax=ax, label="MI")
    for lvl, ls in ((0.5, "-"), (1.9, "--")):
        c = ax.contour(PP/1e3, FF/1e6, MI, [lvl], colors="w", linestyles=ls,
                       linewidths=0.9)
        ax.clabel(c, fmt=f"MI={lvl}", fontsize=6)
    ax.plot(12, 1.5, "o", color=C[7], ms=5)
    ax.annotate("operating point", (12, 1.5), (60, 1.7), color="w",
                fontsize=6, arrowprops=dict(arrowstyle="->", color="w"))
    ax.set(xlabel="Peak rarefactional pressure (kPa)",
           ylabel="Frequency (MHz)", title="Mechanical-index envelope")
    ax.grid(False)
    fig.tight_layout()
    save_fig(fig, "Fig_S8_mechanical_index")
    SUMMARY["mi"] = tab.to_dict("list")


# ======================================================================
def e09_radiation_force(th):
    log("E09 radiation force and displacement (Kelvin solution)")
    duty = 0.2
    c_z = np.array([P.layer_at(z)[3] for z in th["gz"]])
    Fz = 2*th["alpha"][None, None, :]*th["I3"]*duty/c_z[None, None, :]
    # crop to keep the convolution light
    ix = np.abs(th["gx"]) < 12e-3
    iy = np.abs(th["gy"]) < 12e-3
    iz = th["gz"] < 40e-3
    Fc = Fz[np.ix_(ix, iy, iz)]
    u1 = B.kelvin_displacement(Fc, th["d"], mu=1.0)   # times 1/mu
    gz = th["gz"][iz]
    cx, cy = np.argmin(np.abs(th["gx"][ix])), np.argmin(np.abs(th["gy"][iy]))
    rows = []
    axial_u = {}
    for L in P.LAYERS:
        m = (gz >= L[1]) & (gz < L[2])
        Ion = th["I3"][np.argmin(np.abs(th["gx"])),
                       np.argmin(np.abs(th["gy"])), :][iz][m]
        Fon = Fc[cx, cy, m]
        u = u1[cx, cy, m]/L[6]
        strain = np.gradient(u1[cx, cy, :]/L[6], th["d"])[m]
        rows.append([L[0], Ion.max(), Ion.max()*duty, Fon.max(),
                     u.max()*1e9, np.abs(strain).max()])
        axial_u[L[0]] = u
    tab = pd.DataFrame(rows, columns=[
        "Layer", "Peak Isppa on axis (W/m^2)", "Ispta at 20 % (W/m^2)",
        "Peak body force (N/m^3)", "Peak static displacement (nm)",
        "Peak axial strain (-)"])
    save_table(tab, "T_S9_radiation_force",
               "Acoustic radiation force and static displacement, 20 % duty.")
    fig, ax = plt.subplots(1, 2, figsize=(7.16, 2.3))
    ax[0].plot(gz*mm, Fc[cx, cy, :], color=C[6])
    ax[0].set(xlabel="Depth (mm)", ylabel="Body force (N/m^3)",
              title="(a) On-axis radiation body force, 20 % duty")
    for i, (k, u) in enumerate(axial_u.items()):
        L = [l for l in P.LAYERS if l[0] == k][0]
        m = (gz >= L[1]) & (gz < L[2])
        ax[1].plot(gz[m]*mm, u*1e9, color=C[i], label=f"{k} (G={L[6]/1e3:g} kPa)")
    ax[1].set(xlabel="Depth (mm)", ylabel="u_z (nm)",
              title="(b) Static displacement (local shear modulus)")
    ax[1].legend(frameon=False)
    fig.tight_layout()
    save_fig(fig, "Fig_S9_radiation_force")
    SUMMARY["radiation_force"] = tab.to_dict("list")


# ======================================================================
def e10_waveform():
    log("E10 stimulation waveform")
    fs = 50e6
    prf, duty = 100.0, 0.2
    t = np.arange(0, 0.05, 1/fs)
    gate = (((t*prf) % 1) < duty).astype(float)
    x = np.sin(2*np.pi*P.F0*t)*gate
    X = np.abs(np.fft.rfft(x*np.hanning(x.size)))
    fr = np.fft.rfftfreq(x.size, 1/fs)
    Xdb = 20*np.log10(X/X.max()+1e-12)
    fig, ax = plt.subplots(1, 2, figsize=(7.16, 2.2))
    ax[0].fill_between(t*1e3, -gate, gate, step="pre", color=C[0], alpha=0.35,
                       lw=0)
    ax[0].plot(t*1e3, gate, color=C[1], lw=1)
    ax[0].set(xlim=(0, 25), xlabel="Time (ms)", ylabel="Norm. drive",
              title=f"(a) Burst gating: PRF {prf:g} Hz, {duty*100:g} % duty "
                    f"({duty/prf*P.F0:.0f} cycles/burst)")
    ax[1].plot(fr/1e6, Xdb, color=C[2])
    ax[1].set(xlim=(0, 4), ylim=(-100, 3), xlabel="Frequency (MHz)",
              ylabel="Magnitude (dB)", title="(b) Spectrum")
    ins = ax[1].inset_axes([0.6, 0.55, 0.37, 0.4])
    m = np.abs(fr-P.F0) < 600
    ins.plot((fr[m]-P.F0), Xdb[m], color=C[2], lw=0.8)
    ins.set_title("PRF sidebands (Hz)", fontsize=6)
    ins.tick_params(labelsize=5)
    fig.tight_layout()
    save_fig(fig, "Fig_S10_waveform")
    SUMMARY["waveform"] = dict(cycles_per_burst=duty/prf*P.F0,
                               burst_ms=duty/prf*1e3)


# ======================================================================
def e11_control():
    log("E11 closed-loop control")
    t, y, u, r, _ = CL.simulate()
    m = CL.metrics(t, y, r)
    tn, yn, un, rn, ymn = CL.simulate(noise_kpa=0.3, quant_kpa=0.05,
                                      rng=np.random.default_rng(SEED))
    mn = CL.metrics(tn, yn, rn, band=0.05)
    pd.DataFrame(dict(t_s=t, setpoint_kPa=r, focal_kPa=y, drive=u,
                      focal_noisy_kPa=yn)).to_csv(D_DATA / "E11_control.csv",
                                                  index=False)
    # robustness: plant gain and lag +/-50 %
    rng = np.random.default_rng(SEED)
    rob = []
    nmc = 200 if not QUICK else 30
    for _ in range(nmc):
        K = rng.uniform(0.5, 1.5)
        tau = 0.06*rng.uniform(0.5, 1.5)
        th = rng.choice([0.002, 0.005, 0.010])
        tt, yy, _, rr, _ = CL.simulate(K=K, tau=tau, theta=th)
        mm_ = CL.metrics(tt, yy, rr)
        rob.append([K, tau, th, mm_["settling_s"], mm_["overshoot_pct"],
                    mm_["recovery_s"]])
    rob = pd.DataFrame(rob, columns=["K", "tau_s", "theta_s", "settling_s",
                                     "overshoot_pct", "recovery_s"])
    rob.to_csv(D_DATA / "E11_robustness_mc.csv", index=False)
    tab = pd.DataFrame([
        ["Kp / Ki / Kd", "4 / 55 / 0", "4 / 55 / 0"],
        ["Plant (K, tau, dead time)", "FOPDT (n.r.)", "1, 60 ms, 5 ms"],
        ["Settling time, 2 % (ms)", 150, round(m["settling_s"]*1e3, 1)],
        ["Overshoot (%)", 0.0, round(m["overshoot_pct"], 2)],
        ["Steady-state error (Pa)", 0.03, round(m["ss_error_kpa"]*1e3, 3)],
        ["Recovery after 40 % coupling loss (ms)", 426,
         round(m["recovery_s"]*1e3, 1)],
        ["Pressure dip during disturbance (kPa)", "n.r.",
         round(m["disturbance_dip_kpa"], 2)],
        ["Settling, 5 % band, with 0.3 kPa sensor noise (ms)", "n.r.",
         round(mn["settling_s"]*1e3, 1)],
        ["MC +/-50 % plant: settling median [95th pct] (ms)", "n.r.",
         f"{rob.settling_s.median()*1e3:.0f} [{rob.settling_s.quantile(.95)*1e3:.0f}]"],
        ["MC +/-50 % plant: overshoot max (%)", "n.r.",
         round(rob.overshoot_pct.max(), 2)],
    ], columns=["Metric", "Manuscript", "Reproduced"])
    save_table(tab, "T_S10_control", "Closed-loop focal-pressure regulation.")
    fig, ax = plt.subplots(1, 2, figsize=(7.16, 2.3))
    ax[0].plot(t, r, "--", color="0.4", lw=1, label="setpoint")
    ax[0].plot(tn, yn, color=C[0], alpha=0.35, lw=0.7, label="with sensor noise")
    ax[0].plot(t, y, color=C[0], label="focal pressure")
    ax[0].axvspan(1.0, 1.02, color=C[7], alpha=0.2)
    ax[0].annotate("40 % coupling loss", (1.0, 10), (1.2, 5), fontsize=6,
                   arrowprops=dict(arrowstyle="->"))
    ax[0].set(xlabel="Time (s)", ylabel="Focal pressure (kPa)", ylim=(0, 15),
              title="(a) Step and disturbance response")
    ax[0].legend(frameon=False, loc="lower right")
    sc = ax[1].scatter(rob.K, rob.settling_s*1e3, c=rob.tau_s*1e3, s=8,
                       cmap="viridis")
    fig.colorbar(sc, ax=ax[1], label="tau (ms)")
    ax[1].set(xlabel="Plant gain K (rel.)", ylabel="Settling time (ms)",
              title="(b) Robustness to +/-50 % plant variation")
    fig.tight_layout()
    save_fig(fig, "Fig_S11_control")
    SUMMARY["control"] = dict(nominal=m, noisy=mn,
                              mc_settling_median_ms=rob.settling_s.median()*1e3,
                              mc_overshoot_max=rob.overshoot_pct.max())


# ======================================================================
def e12_dose_and_power():
    log("E12 dose-response hypothesis and sample-size planning")
    from scipy.stats import norm, nct, t as tdist
    p = np.linspace(0, 60, 601)
    Emax, p50, n = 40.0, 25.0, 2.0       # placeholders (to be estimated)
    prolif = 100 - Emax*p**n/(p50**n + p**n)
    pd.DataFrame(dict(p_kPa=p, proliferation_pct=prolif)).to_csv(
        D_DATA / "E12_dose_response_model.csv", index=False)

    def n_per_group(delta, sd, alpha=0.05, power=0.8, k=1):
        a = alpha/k
        for nn in range(3, 200):
            df = 2*nn - 2
            ncp = delta/sd*np.sqrt(nn/2)
            tc = tdist.ppf(1-a/2, df)
            pw = 1 - nct.cdf(tc, df, ncp) + nct.cdf(-tc, df, ncp)
            if pw >= power:
                return nn
        return np.nan
    rows = []
    for eff in (10, 15, 20, 30):
        for cv in (10, 15, 20):
            rows.append([eff, cv, n_per_group(eff, cv),
                         n_per_group(eff, cv, k=4)])
    tab = pd.DataFrame(rows, columns=[
        "Effect (% reduction)", "CV of spheroid assay (%)",
        "n/group (2 groups)", "n/group (4 doses vs sham, Bonferroni)"])
    save_table(tab, "T_S11_sample_size",
               "Phase III spheroid sample size (two-sided t-test, power 0.8).")
    fig, ax = plt.subplots(figsize=(3.5, 2.4))
    ax.axvspan(10, 15, color=C[3], alpha=0.25, lw=0)
    ax.text(10.3, 101, "operating\nwindow", fontsize=6, color="#8a5d00")
    ax.plot(p, prolif, color=C[0])
    ax.axhline(80, ls="--", color=C[7], lw=0.8)
    ax.text(30, 81, "20 % suppression target", fontsize=6, color=C[7])
    ax.set(xlabel="Peak focal pressure (kPa)",
           ylabel="Proliferation (% of sham)", ylim=(55, 105),
           title="Dose-response MODEL (hypothesis, Emax=40, p50=25, n=2)")
    fig.tight_layout()
    save_fig(fig, "Fig_S12_dose_response_model")
    SUMMARY["dose"] = dict(model_at_12kPa=float(np.interp(12, p, prolif)),
                           n_20pct_cv15=int(tab.query(
                               "`Effect (% reduction)`==20 and "
                               "`CV of spheroid assay (%)`==15")
                               .iloc[0, 2]))


def e13_power_budget():
    log("E13 wearable power budget (estimate)")
    Vac = SUMMARY["focus"]["Vac_V"]
    C0 = SUMMARY["cmut"]["elem_C0_pF"]*1e-12
    duty = 0.2
    p_cv2f = P.N_ELEM*C0*Vac**2*P.F0*duty
    items = [
        ["Low-power FPGA (beamformer, sequencer, PI)", 150.0, "vendor class"],
        ["64-ch low-voltage TX drivers (linear, quiescent)", 64*3.0,
         "3 mW/ch assumed"],
        ["TX dynamic C*V^2*f*duty (computed)", p_cv2f*1e3,
         f"Vac={Vac:.2f} V, C0={C0*1e12:.0f} pF"],
        ["DC bias generator (~41 V, leakage-limited)", 40.0, "boost + LDO"],
        ["Feedback RX/ADC (8 channels, duty-cycled)", 50.0, "assumed"],
        ["MCU + BLE link", 15.0, "assumed"],
    ]
    tab = pd.DataFrame(items, columns=["Block", "Power (mW)", "Basis"])
    tot = tab["Power (mW)"].sum()
    tab.loc[len(tab)] = ["TOTAL", tot, ""]
    tab.loc[len(tab)] = ["Energy for 4 h session (Wh)", tot*4/1e3, ""]
    tab.loc[len(tab)] = ["Battery @3.7 V, 80 % usable (mAh)",
                         tot*4/3.7/0.8, ""]
    save_table(tab, "T_S12_power_budget",
               "Estimated wearable power budget (assumption-based).")
    SUMMARY["power"] = dict(total_mW=tot, battery_mAh=tot*4/3.7/0.8)


# ======================================================================
# HARDWARE DIGITAL TWIN (emulated measurements)
# ======================================================================
def hw1_fabrication():
    log("HW1 fabrication Monte Carlo")
    rng = np.random.default_rng(SEED + 1)
    n_arr = 300 if not QUICK else 40
    Vc_all, sens_all, ph_all, Vmin = [], [], [], []
    for _ in range(n_arr):
        arr = HWT.fabricate_array(rng)
        s, dphi, Vdc = HWT.element_response(arr)
        Vc_all.append(arr["Vc"])
        sens_all.append(20*np.log10(s))
        ph_all.append(np.rad2deg(dphi))
        Vmin.append(arr["Vc"].min())
    Vc_all = np.concatenate(Vc_all)
    sens_all = np.concatenate(sens_all)
    ph_all = np.concatenate(ph_all)
    Vmin = np.array(Vmin)
    # yield vs a fixed global bias (no per-array trim)
    biases = np.linspace(30, 50, 81)
    yld = [(Vmin > b).mean()*100 for b in biases]
    pd.DataFrame(dict(Vc_V=Vc_all, sens_dB=sens_all, phase_deg=ph_all)).to_csv(
        D_EMU / "HW1_element_population.csv", index=False)
    tab = pd.DataFrame([
        ["Arrays simulated (64 elements each)", n_arr],
        ["Collapse voltage mean +/- SD (V)",
         f"{Vc_all.mean():.1f} +/- {Vc_all.std():.1f}"],
        ["Collapse voltage 1st-99th percentile (V)",
         f"{np.percentile(Vc_all,1):.1f} - {np.percentile(Vc_all,99):.1f}"],
        ["Element sensitivity spread SD (dB)", f"{sens_all.std():.2f}"],
        ["Element phase spread SD (deg)", f"{ph_all.std():.2f}"],
        ["Array yield at fixed 41 V bias (no trim) (%)",
         f"{np.interp(41, biases, yld):.1f}"],
        ["Bias for 99 % array yield (V)",
         f"{np.percentile(Vmin, 1):.1f}"],
    ], columns=["Quantity", "Emulated value"])
    save_table(tab, "T_E1_fabrication_mc",
               "Fabrication Monte Carlo (EMULATED).", emulated=True)
    fig, ax = plt.subplots(1, 3, figsize=(7.16, 2.2))
    ax[0].hist(Vc_all, 50, color=C[0])
    ax[0].axvline(46.5, ls="--", color=C[7], lw=0.8)
    ax[0].set(xlabel="Collapse voltage (V)", ylabel="Elements",
              title="(a) Vc distribution")
    ax[1].hist(sens_all, 50, color=C[2])
    ax[1].set(xlabel="Relative TX sensitivity (dB)", title="(b) Sensitivity")
    ax[2].plot(biases, yld, color=C[1])
    ax[2].set(xlabel="Global DC bias (V)", ylabel="Array yield (%)",
              title="(c) No-collapse yield")
    fig.tight_layout()
    save_fig(fig, "Fig_E1_fabrication_mc", emulated=True)
    SUMMARY["hw_fab"] = dict(Vc_mean=Vc_all.mean(), Vc_sd=Vc_all.std(),
                             sens_sd_dB=sens_all.std(),
                             phase_sd_deg=ph_all.std(),
                             bias_99yield=np.percentile(Vmin, 1))


def hw2_beam_mc():
    log("HW2 beam Monte Carlo with hardware non-idealities")
    rng = np.random.default_rng(SEED + 2)
    nmc = 150 if not QUICK else 20
    s0 = SUMMARY["focus"]["drive_scale"]
    conds = [("Ideal (model)", None, False),
             ("Flat patch + HW errors", None, False),
             ("Curved R=60 mm, uncorrected", 60e-3, False),
             ("Curved R=60 mm, shape-corrected", 60e-3, True)]
    recs = []
    for name, R, corr in conds:
        nn = 1 if name.startswith("Ideal") else nmc
        for _ in range(nn):
            if name.startswith("Ideal"):
                kw = {}
            else:
                ni = HWT.realise_nonidealities(rng, curvature_R=R,
                                               corrected=corr)
                kw = {k: ni[k] for k in ("amp_err", "phase_err", "dead",
                                         "delay_quant", "z_offset")}
            pr = focal_profiles(**kw)
            pk = np.abs(pr["pa"]).max()*s0
            lw, psl = F.mainlobe_and_psl(pr["x"], pr["pl"])
            p_tgt = np.interp(P.FOCUS_Z, pr["r"], np.abs(pr["pa"]))*s0
            recs.append([name, pk/1e3, p_tgt/1e3, lw*mm, F.fwhm(pr["r"], pr["pa"])*mm,
                         (pr["zpk"]-P.FOCUS_Z)*mm, pr["xpk"]*mm, psl])
    df = pd.DataFrame(recs, columns=["Condition", "Peak (kPa)",
                                     "Pressure at 20 mm target (kPa)",
                                     "Lateral FWHM (mm)", "Axial FWHM (mm)",
                                     "Axial focus shift (mm)",
                                     "Lateral focus shift (mm)",
                                     "Peak side lobe (dB)"])
    df.to_csv(D_EMU / "HW2_beam_mc_raw.csv", index=False)
    agg = df.groupby("Condition", sort=False).agg(["mean", "std"])
    out = pd.DataFrame({"Condition": agg.index})
    for col in df.columns[1:]:
        out[col] = [f"{m:.3f} +/- {s:.3f}" if not np.isnan(s) else f"{m:.3f}"
                    for m, s in zip(agg[(col, "mean")], agg[(col, "std")])]
    save_table(out, "T_E2_beam_mc",
               f"Beam quality under hardware non-idealities, n={nmc} "
               "(EMULATED).", emulated=True)
    fig, ax = plt.subplots(1, 3, figsize=(7.16, 2.3))
    names = [c[0] for c in conds[1:]]
    for i, col in enumerate(("Pressure at 20 mm target (kPa)",
                             "Axial focus shift (mm)",
                             "Peak side lobe (dB)")):
        data = [df[df.Condition == n][col] for n in names]
        bp = ax[i].boxplot(data, widths=0.5, patch_artist=True,
                           medianprops=dict(color="k"))
        for b, cc in zip(bp["boxes"], C):
            b.set_facecolor(cc)
            b.set_alpha(0.6)
        ideal = df[df.Condition == "Ideal (model)"][col].iloc[0]
        ax[i].axhline(ideal, ls="--", color="0.3", lw=0.8)
        ax[i].set_xticks([1, 2, 3], ["flat", "curved\nuncorr.",
                                     "curved\ncorr."])
        ax[i].set_title(col, fontsize=7.5)
    fig.tight_layout()
    save_fig(fig, "Fig_E2_beam_mc", emulated=True)
    SUMMARY["hw_beam"] = out.to_dict("list")


def hw3_hydrophone():
    log("HW3 emulated hydrophone characterisation (n = 5 scans)")
    rng = np.random.default_rng(SEED + 3)
    ni = HWT.realise_nonidealities(rng)      # one 'fabricated' patch
    kw = {k: ni[k] for k in ("amp_err", "phase_err", "dead", "delay_quant")}
    s0 = SUMMARY["focus"]["drive_scale"]
    pr = focal_profiles(**kw)
    zpk = pr["zpk"]
    lat_fn = lambda x: F.field_at(x, np.full_like(x, zpk), **kw)*s0
    ax_fn = lambda z: F.field_at(np.zeros_like(z), z, **kw)*s0
    xs = np.arange(-3, 3.0001, 0.05)*1e-3
    zs = np.arange(10, 30.0001, 0.25)*1e-3
    model_lat = np.abs(F.field_at(xs, np.full_like(xs, P.FOCUS_Z)))*s0
    model_ax = np.abs(F.field_at(np.zeros_like(zs), zs))*s0
    nrep = 5
    L = np.array([HWT.hydrophone_scan(xs, lat_fn, rng) for _ in range(nrep)])
    A = np.array([HWT.hydrophone_scan(zs, ax_fn, rng) for _ in range(nrep)])
    stats = []
    for i in range(nrep):
        stats.append([i+1, L[i].max()/1e3, F.fwhm(xs, L[i])*mm,
                      F.fwhm(zs, A[i])*mm, zs[np.argmax(A[i])]*mm])
    st = pd.DataFrame(stats, columns=["Scan", "Peak (kPa)",
                                      "Lateral FWHM (mm)",
                                      "Axial FWHM (mm)", "Focus depth (mm)"])
    st.to_csv(D_EMU / "HW3_hydrophone_scans_summary.csv", index=False)
    pd.DataFrame(np.column_stack([xs*mm, L.T/1e3]),
                 columns=["x_mm"]+[f"scan{i+1}_kPa" for i in range(nrep)]
                 ).to_csv(D_EMU / "HW3_lateral_scans.csv", index=False)
    pd.DataFrame(np.column_stack([zs*mm, A.T/1e3]),
                 columns=["z_mm"]+[f"scan{i+1}_kPa" for i in range(nrep)]
                 ).to_csv(D_EMU / "HW3_axial_scans.csv", index=False)
    from scipy.stats import t as tdist
    tq = tdist.ppf(0.975, nrep-1)
    rows = []
    for col, mval in (("Peak (kPa)", 12.0),
                      ("Lateral FWHM (mm)", SUMMARY["focus"]["lat6_mm"]),
                      ("Axial FWHM (mm)", SUMMARY["focus"]["ax6_mm"]),
                      ("Focus depth (mm)", SUMMARY["focus"]["z_peak_mm"])):
        m, sd = st[col].mean(), st[col].std(ddof=1)
        rows.append([col, round(mval, 3), f"{m:.3f} +/- {sd:.3f}",
                     f"[{m-tq*sd/np.sqrt(nrep):.3f}, {m+tq*sd/np.sqrt(nrep):.3f}]",
                     f"{(m-mval)/mval*100:+.1f}"])
    tab = pd.DataFrame(rows, columns=["Quantity", "Model (ideal)",
                                      "Emulated mean +/- SD (n=5)",
                                      "95 % CI", "Deviation (%)"])
    save_table(tab, "T_E3_hydrophone",
               "Emulated needle-hydrophone characterisation (EMULATED).",
               emulated=True)
    fig, ax = plt.subplots(1, 2, figsize=(7.16, 2.3))
    for a, x, M, mdl, lab in ((ax[0], xs, L, model_lat, "Azimuth (mm)"),
                              (ax[1], zs, A, model_ax, "Depth (mm)")):
        mu, sd = M.mean(0)/1e3, M.std(0, ddof=1)/1e3
        a.fill_between(x*mm, mu-sd, mu+sd, color=C[1], alpha=0.25, lw=0,
                       label="emulated mean +/- SD (n=5)")
        a.plot(x*mm, mu, "o", ms=2, color=C[1])
        a.plot(x*mm, mdl/1e3, color=C[0], label="model")
        a.set(xlabel=lab, ylabel="|p| (kPa)")
    ax[0].set_title("(a) Lateral scan at focus")
    ax[1].set_title("(b) Axial scan")
    ax[0].legend(frameon=False, fontsize=6, loc="upper left")
    fig.tight_layout()
    save_fig(fig, "Fig_E3_hydrophone", emulated=True)
    SUMMARY["hw_hydrophone"] = tab.to_dict("list")


def hw4_thermocouple():
    log("HW4 emulated phantom thermometry")
    rng = np.random.default_rng(SEED + 4)
    dT20 = [v for k, v in zip(SUMMARY["thermal"]["Thermal props"],
                              zip(SUMMARY["thermal"]["Duty (%)"],
                                  SUMMARY["thermal"]["Steady dT (mK)"]))
            if k == "Tumour" and v[0] == 20][0][1]/1e3
    rows = []
    for p_kpa in (12, 50, 100, 200, 400):
        true = dT20*(p_kpa/12)**2*5     # 100 % duty, scales as p^2
        n = 10
        rd, art = HWT.thermocouple_reading(true, rng, n)
        lod = 3*HWT.HW["tc_noise_k"]/np.sqrt(n)
        rows.append([p_kpa, true*1e3, art*1e3, rd.mean()*1e3,
                     rd.std(ddof=1)*1e3, lod*1e3,
                     "yes" if true > lod else "no"])
    tab = pd.DataFrame(rows, columns=[
        "Focal p (kPa)", "True dT 100 % duty (mK)",
        "Viscous artefact (mK)", "Emulated reading mean (mK)",
        "SD (mK)", "LOD 3sigma (mK)", "Resolvable?"])
    save_table(tab, "T_E4_thermocouple",
               "Emulated thermocouple in hydrogel phantom (EMULATED). "
               "dT scales with p^2; thermal dose at 12 kPa is below the "
               "detection limit and must be verified by p^2 extrapolation.",
               emulated=True)
    SUMMARY["hw_thermo"] = tab.to_dict("list")


def hw5_session():
    log("HW5 emulated 4-hour closed-loop session")
    rng = np.random.default_rng(SEED + 5)
    dt = 0.01                       # one control update per burst (100 Hz)
    T = 4*3600 if not QUICK else 1800
    n = int(T/dt)
    # coupling drift: slow gel drying (-25 % over 4 h) + random walk +
    # motion events (brief 30-60 % losses)
    drift = 1 - 0.25*np.arange(n)/n
    drift *= np.exp(np.cumsum(rng.standard_normal(n))*0.0015)
    for _ in range(int(T/600)):
        i0 = rng.integers(0, n-300)
        drift[i0:i0+rng.integers(50, 300)] *= rng.uniform(0.4, 0.7)
    t, y, u, r, _ = CL.simulate(dt=dt, T=T, t_step=0.0, dist_t=T+1,
                                drift=drift, noise_kpa=0.3, quant_kpa=0.05,
                                rng=rng, theta=0.01)
    inwin = ((y >= 10) & (y <= 15))[int(5/dt):]
    within2 = (np.abs(y-12) <= 0.6)[int(5/dt):]
    ds = slice(None, None, 100)
    pd.DataFrame(dict(t_s=t[ds], focal_kPa=y[ds], drive=u[ds],
                      coupling=drift[ds])).to_csv(
        D_EMU / "HW5_4h_session_decimated.csv", index=False)
    tab = pd.DataFrame([
        ["Session length (h)", T/3600],
        ["Time within 10-15 kPa window (%)", inwin.mean()*100],
        ["Time within +/-5 % of 12 kPa (%)", within2.mean()*100],
        ["Mean focal pressure (kPa)", y[int(5/dt):].mean()],
        ["SD focal pressure (kPa)", y[int(5/dt):].std()],
        ["Drive increase end vs start (dB)",
         20*np.log10(u[-1000:].mean()/u[1000:2000].mean())],
    ], columns=["Quantity", "Emulated value"])
    save_table(tab, "T_E5_session",
               "Emulated 4 h closed-loop session with coupling drift "
               "(EMULATED).", emulated=True)
    fig, ax = plt.subplots(2, 1, figsize=(7.16, 2.8), sharex=True)
    ax[0].plot(t[ds]/3600, y[ds], lw=0.5, color=C[0])
    ax[0].axhspan(10, 15, color=C[2], alpha=0.12, lw=0)
    ax[0].set(ylabel="Focal p (kPa)", ylim=(0, 18),
              title="Emulated 4 h session: focal pressure (band = 10-15 kPa "
                    "window)")
    ax[1].plot(t[ds]/3600, drift[ds], lw=0.5, color=C[1])
    ax[1].set(xlabel="Time (h)", ylabel="Coupling gain")
    fig.tight_layout()
    save_fig(fig, "Fig_E5_4h_session", emulated=True)
    SUMMARY["hw_session"] = tab.to_dict("list")


# ======================================================================
def comparison_table():
    log("Manuscript vs reproduced vs emulated")
    f, s = SUMMARY["focus"], SUMMARY["steering"]
    c = SUMMARY["cmut"]
    th = SUMMARY["thermal"]
    dT20 = [v for k, d_, v in zip(th["Thermal props"], th["Duty (%)"],
                                  th["Steady dT (mK)"])
            if k == "Tumour" and d_ == 20][0]
    TIS20 = [v for k, d_, v in zip(th["Thermal props"], th["Duty (%)"],
                                   th["TIS"]) if k == "Tumour" and d_ == 20][0]
    hb = SUMMARY["hw_beam"]
    hy = SUMMARY["hw_hydrophone"]
    rows = [
        ["Collapse voltage (V)", "46.5", f"{c['Vc']:.1f}",
         f"{SUMMARY['hw_fab']['Vc_mean']:.1f} +/- {SUMMARY['hw_fab']['Vc_sd']:.1f}",
         "Within 11 %; depends on E and insulator (assumed)"],
        ["Immersed cell resonance (MHz)", "1.5 (implied)",
         f"{c['f_imm_bias_MHz']:.2f}", "-",
         f"MISMATCH: a ~ {c['radius_for_1p5MHz_um']:.1f} um needed for 1.5 MHz"],
        ["-6 dB fractional BW (%)", "51.8", f"{c['fbw6']*100:.1f}", "-",
         "Consistent with Q = 3.5"],
        ["Lateral FWHM (mm)", "0.69", f"{f['lat6_mm']:.2f} (-6 dB) / "
         f"{f['lat3_mm']:.2f} (-3 dB)", hy["Emulated mean +/- SD (n=5)"][1],
         "Manuscript value equals the -3 dB (intensity) width"],
        ["Axial FWHM (mm)", "3.74", f"{f['ax6_mm']:.2f} (-6 dB) / "
         f"{f['ax3_mm']:.2f} (-3 dB)", hy["Emulated mean +/- SD (n=5)"][2],
         "Manuscript value close to -3 dB width"],
        ["Peak focal pressure (kPa)", "12.0", "12.0 (set)",
         hy["Emulated mean +/- SD (n=5)"][0],
         "Hardware spread from cal. uncertainty + element errors"],
        ["Steering loss at +/-20 deg (dB)", "-0.33",
         f"{s['Peak change (dB)'][0]:.2f}", "-", ""],
        ["Isppa (mW/cm^2)", "4.45", f"{f['Isppa_mWcm2']:.2f}", "-",
         "rho*c = 1050*1540"],
        ["MI", "0.010", f"{B.mechanical_index(12e3):.4f}", "-", "Agrees"],
        ["Peak dT, 20 % duty (mK)", "0.37-0.4", f"{dT20:.3f}",
         "below thermocouple LOD", "Same order; 3-D model"],
        ["Thermal index", "0.0004 (= dT)", f"TIS {TIS20:.2e}", "-",
         "Manuscript TI equals dT numerically; use IEC 62359 TIS/TIB"],
        ["Settling time (ms)", "150",
         f"{SUMMARY['control']['nominal']['settling_s']*1e3:.0f}", "-",
         "Plant fitted (FOPDT)"],
        ["Disturbance recovery (ms)", "426",
         f"{SUMMARY['control']['nominal']['recovery_s']*1e3:.0f}", "-",
         "Depends on disturbance definition (n.r. in manuscript)"],
        ["Peak with HW errors, flat (kPa)", "-", "12.0",
         hb["Peak (kPa)"][1], "Expect ~0.5-1 dB loss on hardware"],
        ["Pressure at target, curved R=60 mm uncorrected (kPa)", "-", "12.0",
         hb["Pressure at 20 mm target (kPa)"][2],
         "Focus moves toward skin; shape sensing required"],
        ["Pressure at target, curved R=60 mm corrected (kPa)", "-", "12.0",
         hb["Pressure at 20 mm target (kPa)"][3], "Restored by correction"],
    ]
    tab = pd.DataFrame(rows, columns=["Metric", "Manuscript",
                                      "Reproduced (model)",
                                      "Emulated hardware", "Note"])
    save_table(tab, "T_S0_manuscript_vs_reproduced",
               "Manuscript values vs independent reproduction vs emulated "
               "hardware (digital twin).")


def main():
    t0 = time.time()
    e01_cmut()
    e02_focus()
    e03_steering()
    e04_grating()
    e05_apod()
    pm = e06_layered()
    th = e07_thermal(pm)
    e08_mi()
    e09_radiation_force(th)
    e10_waveform()
    e11_control()
    e12_dose_and_power()
    e13_power_budget()
    hw1_fabrication()
    hw2_beam_mc()
    hw3_hydrophone()
    hw4_thermocouple()
    hw5_session()
    comparison_table()
    SUMMARY["meta"] = dict(seed=SEED, quick=QUICK,
                           runtime_s=round(time.time()-t0, 1),
                           numpy=np.__version__, pandas=pd.__version__)
    (ROOT / "summary.json").write_text(json.dumps(SUMMARY, indent=2,
                                                  default=float))
    log(f"done in {time.time()-t0:.0f} s")


if __name__ == "__main__":
    main()
