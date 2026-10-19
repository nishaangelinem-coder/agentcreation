#!/usr/bin/env python3
"""Assemble supplementary/S1_Supplementary_Methods_and_Results.md from the
result tables and convert it to .docx with pandoc (run after run_all.py)."""
import json
import shutil
import subprocess
import tempfile
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parent
R = HERE / "results"
OUT = HERE / "supplementary" / "S1_Supplementary_Methods_and_Results.md"
S = json.loads((R / "summary.json").read_text())


def md_table(path, digits=3):
    df = pd.read_csv(path)

    def fmt(v):
        if isinstance(v, float) and v != v:
            return "-"
        if isinstance(v, str):
            return v.replace("*", "\\*").replace("|", "\\|")
        if isinstance(v, float):
            if v == 0:
                return "0"
            if abs(v) < 1e-2 or abs(v) >= 1e5:
                return f"{v:.3g}"
            return f"{v:.{digits}g}" if abs(v) < 1 else f"{v:.{digits+1}g}"
        return str(v)
    head = "| " + " | ".join(df.columns) + " |"
    sep = "|" + "---|" * len(df.columns)
    rows = ["| " + " | ".join(fmt(v) for v in r) + " |" for r in df.itertuples(index=False)]
    return "\n".join([head, sep] + rows)


def fig(rel, cap):
    return f"![{cap}](../results/{rel}){{width=100%}}\n"


f = S["focus"]
c = S["cmut"]
parts = []
A = parts.append
A(f"""---
title: "Supplementary Information"
subtitle: "A Wearable Phased-Array CMUT Ultrasound Platform for Mechanotransductive Modulation of Tumour Proliferation: An In-Silico Design and Feasibility Study"
author: "M. Nisha Angeline and S. K. Manikandan"
---

# Contents

* Note S1. Model implementation and assumptions
* Note S2. Governing equations
* Note S3. Reproduction results (Tables S0–S12, Figs S1–S12)
* Note S4. Hardware digital twin: emulated bench characterisation (Tables E1–E5, Figs E1–E5)
* Note S5. Code and data availability

> **Status of results.** Every number in this document comes from a
> simulation. Notes S1–S3 are deterministic model results. Note S4 is a
> Monte Carlo emulation of expected hardware measurements. No device has
> been fabricated and no measurement has been made. All Note S4 items are
> labelled *EMULATED*.

# Note S1. Model implementation and assumptions

The models were implemented independently in Python (NumPy/SciPy) in the
open-source package `cmutsim`. One script, `run_all.py`, regenerates every
table and figure, with a fixed random seed ({S['meta']['seed']}). A full
run takes about {S['meta']['runtime_s']:.0f} s on a single CPU core.

**Table S-A. Parameters used in addition to manuscript Tables I–II.**

| Quantity | Value | Basis |
|---|---|---|
| Plate material | Si₃N₄, E = 250 GPa, ν = 0.23, ρ = 3100 kg/m³ | Literature values for LPCVD nitride |
| Insulation layer | 0.15 µm SiO₂ (ε_r = 3.9) | Typical CMUT process |
| Element elevation height | 6 mm (unfocused) | Assumed |
| Cell pitch / cells per element | 66 µm / {c['n_cells']} | Assumed |
| Element clamped capacitance | {c['elem_C0_pF']:.1f} pF | Derived |
| Effective coupling k² (BVD) | 0.5 | Assumed (biased near 0.8 V_c) |
| Homogeneous attenuation | 0.54 dB/cm/MHz | Soft tissue |
| Field discretisation | 4 strips per element; Fresnel elevation integral | — |
| Field grid | 0.1 mm (maps); 5 µm lateral, 20 µm axial (profiles) | Convergence-checked |
| Thermal/elastic grid | 0.3 mm, 3-D, 48 × 48 × 64 mm domain | — |
| Thermal properties | IT'IS v4.1 (tumour; fat as worst case) | Ref. [18] of the manuscript |

# Note S2. Governing equations

**CMUT cell.** Clamped-plate rigidity D = Et³/12(1−ν²) and piston-equivalent
stiffness k = 64πD/a². The collapse voltage is

$$V_c=\\sqrt{{\\frac{{8k\\,g_\\mathrm{{eff}}^3}}{{27\\varepsilon_0\\pi a^2}}}},\\qquad g_\\mathrm{{eff}}=g_0+t_i/\\varepsilon_r .$$

The vacuum resonance is f = 10.21/(2πa²)·√(D/ρt). Fluid loading enters
through the added-mass factor β = 0.6689(ρ_f/ρ_p)(a/t), and the bias
softens the spring as k − ε₀AV²/(g−x)³.

**Acoustic field.** The Rayleigh–Sommerfeld integral is discretised over
element strips, and the elevation integral is evaluated in closed form
(Fresnel):

$$p(\\mathbf r)=\\sum_n w_n e^{{-j\\omega\\tau_n}}\\frac{{jk}}{{2\\pi}}\\,\\frac{{1+\\cos\\theta_n}}{{2}}\\,\\frac{{e^{{-jkR_n-\\alpha R_n}}}}{{R_n}}\\int_{{-H/2}}^{{H/2}}e^{{-jky^2/2R_n}}dy .$$

In the layered model, phase, attenuation and normal-incidence transmission
coefficients accumulate along straight rays through each layer.

**Safety indices.** I = |p|²/2ρc and MI = p_r[MPa]/√f[MHz]. TIS = W·f/210
mW·MHz and TIB = √(W·I_ta)/50, following IEC 62359.

**Bioheat.** The Pennes equation is solved exactly in k-space:

$$\\hat T(\\mathbf k,t)=\\frac{{\\hat q(\\mathbf k)}}{{\\rho c\\,\\lambda(\\mathbf k)}}\\left(1-e^{{-\\lambda t}}\\right),\\quad \\lambda=\\frac{{\\kappa}}{{\\rho c}}|\\mathbf k|^2+\\frac{{w_b\\rho_b c_b}}{{\\rho c}},\\quad q=2\\alpha I_{{ta}}.$$

**Radiation force.** The body force is F = 2αI_ta/c. The static
displacement uses the incompressible Kelvin solution
u_z = F/(8πμr)(1 + Δz²/r²), convolved by FFT.

**Control.** The plant is first-order-plus-dead-time,
G(s) = K e^(−θs)/(τs + 1), with K = 1, τ = 60 ms and θ = 5 ms (zero-order
hold of a 100 Hz estimate). The controller is a discrete PI (K_p = 4,
K_i = 55 s⁻¹) with conditional-integration anti-windup, and the plant
output is clipped at 0.

**Sample size.** A two-sided two-sample t-test, using the non-central t
distribution, at α = 0.05 (Bonferroni-corrected for 4 doses) and power 0.8.

# Note S3. Reproduction results

**Table S0. Manuscript values vs independent reproduction vs emulated hardware.**

{md_table(R/'tables/T_S0_manuscript_vs_reproduced.csv')}

**Table S1. CMUT cell and element.**

{md_table(R/'tables/T_S1_cmut_cell.csv')}

{fig('figures/Fig_S1_cmut_element.png', 'Fig. S1. (a) Element input impedance (BVD). (b) Normalised transmit response. (c) Immersed first resonance against plate radius, at 0 V and 0.8 Vc. The 1.5 MHz line crosses the biased curve at a ≈ ' + f"{c['radius_for_1p5MHz_um']:.1f}" + ' µm.')}

**Table S2. Focused field at 20 mm (homogeneous medium, rectangular apodization).**

{md_table(R/'tables/T_S2_focus.csv')}

{fig('figures/Fig_S2_focal_field.png', 'Fig. S2. (a) Pressure map at 12 kPa. (b) Lateral profile through the peak. (c) Axial profile. Widths are given at both −6 dB (pressure FWHM) and −3 dB (intensity FWHM).')}

**Table S3. Beam steering at a 20 mm focal range.**

{md_table(R/'tables/T_S3_steering.csv')}

{fig('figures/Fig_S3_steering.png', 'Fig. S3. Steered fields from −20° to +20°, on a common colour scale.')}

**Table S4. Grating lobes against pitch, 20° steer, evaluated on the 20 mm arc.**

{md_table(R/'tables/T_S4_grating.csv')}

{fig('figures/Fig_S4_grating_lobes.png', 'Fig. S4. Arc beam patterns for pitches of 0.5–1.5 λ. A grating lobe is predicted at sin θg = sin θs − λ/p.')}

**Table S5. Apodization trade-off.**

{md_table(R/'tables/T_S5_apodization.csv')}

{fig('figures/Fig_S5_apodization.png', 'Fig. S5. Lateral beam profiles at the focus for four apodization windows.')}

**Table S6. Effect of the skin–fat–muscle–tumour stack.**

{md_table(R/'tables/T_S6_layered.csv')}

{fig('figures/Fig_S6_layered_tissue.png', 'Fig. S6. (a) Field in the layered model with layer-corrected delays. (b) Sound-speed profile. (c) On-axis pressure with homogeneous and with layer-corrected delays.')}

**Table S7. Pennes bioheat (3-D) and IEC thermal indices.**

{md_table(R/'tables/T_S7_thermal.csv')}

{fig('figures/Fig_S7_thermal.png', 'Fig. S7. (a) Steady-state temperature rise at 20 % duty. (b) Peak-temperature transients over 4 min. Solid lines: tumour thermal properties. Dashed lines: fat thermal properties (worst case).')}

**Table S8. Mechanical index and cavitation margin.**

{md_table(R/'tables/T_S8_mi.csv')}

{fig('figures/Fig_S8_mechanical_index.png', 'Fig. S8. Mechanical index over the frequency–pressure plane, with the MI = 0.5 design limit and the MI = 1.9 regulatory limit marked.')}

**Table S9. Acoustic radiation force and static displacement (20 % duty).**

{md_table(R/'tables/T_S9_radiation_force.csv')}

{fig('figures/Fig_S9_radiation_force.png', 'Fig. S9. (a) On-axis radiation body force. (b) Static axial displacement in each layer, from the Kelvin solution with the local shear modulus.')}

{fig('figures/Fig_S10_waveform.png', 'Fig. S10. (a) Burst gating at a 100 Hz PRF and 20 % duty: 3000 cycles per 2 ms burst. (b) Spectrum, with the PRF sidebands shown inset.')}

**Table S10. Closed-loop regulation.**

{md_table(R/'tables/T_S10_control.csv')}

{fig('figures/Fig_S11_control.png', 'Fig. S11. (a) Step response and response to a 40 % coupling loss, with and without sensor noise. (b) Settling time across 200 plants with ±50 % uncertainty in gain and lag.')}

**Table S11. Phase III spheroid sample size.**

{md_table(R/'tables/T_S11_sample_size.csv')}

{fig('figures/Fig_S12_dose_response_model.png', 'Fig. S12. Hypothesised Hill dose–response model. The parameters are placeholders to be estimated in Phase III. This is not a prediction.')}

**Table S12. Estimated wearable power budget.**

{md_table(R/'tables/T_S12_power_budget.csv')}

# Note S4. Hardware digital twin (EMULATED results)

> The tables and figures in this Note are **emulated**. They are Monte
> Carlo predictions of what a bench characterisation would be expected to
> report, given the error budget in Supplementary file S3. They are **not
> experimental measurements**. They set the acceptance bands for the
> Phase II protocol (Supplementary file S4).

**Table E1. Fabrication Monte Carlo (EMULATED).**

{md_table(R/'emulated_measurements/T_E1_fabrication_mc.csv')}

{fig('emulated_measurements/Fig_E1_fabrication_mc.png', 'Fig. E1 (EMULATED). (a) Collapse-voltage distribution; the dashed line marks the manuscript value of 46.5 V. (b) Element sensitivity spread. (c) Array yield against a global DC bias.')}

**Table E2. Beam quality under hardware non-idealities (EMULATED).**

{md_table(R/'emulated_measurements/T_E2_beam_mc.csv')}

{fig('emulated_measurements/Fig_E2_beam_mc.png', 'Fig. E2 (EMULATED). Monte Carlo results for a flat patch, an uncorrected curved patch (R = 60 mm) and a shape-corrected curved patch. The dashed lines show the ideal model.')}

**Table E3. Emulated needle-hydrophone characterisation, n = 5 scans (EMULATED).**

{md_table(R/'emulated_measurements/T_E3_hydrophone.csv')}

{fig('emulated_measurements/Fig_E3_hydrophone.png', 'Fig. E3 (EMULATED). Emulated hydrophone scans (mean ± SD, n = 5) against the model.')}

**Table E4. Emulated phantom thermometry (EMULATED).**

{md_table(R/'emulated_measurements/T_E4_thermocouple.csv')}

**Table E5. Emulated 4 h closed-loop session (EMULATED).**

{md_table(R/'emulated_measurements/T_E5_session.csv')}

{fig('emulated_measurements/Fig_E5_4h_session.png', 'Fig. E5 (EMULATED). Focal pressure over a 4 h session with gel drying, a random-walk drift in coupling, and motion events.')}

# Note S5. Code and data availability

The code (`cmutsim`, `run_all.py`), all raw numerical data (CSV/NPZ) and
the figures are released under the MIT licence. They will be deposited in
a public repository with a DOI before publication (see file S5).
""")
OUT.write_text("\n".join(parts))
print("wrote", OUT)
subprocess.run(["pandoc", OUT.name, "-o", OUT.with_suffix(".docx").name,
                "--resource-path=."], cwd=OUT.parent, check=True)
print("wrote", OUT.with_suffix(".docx"))

# PDF preview: pandoc -> standalone HTML (MathML) -> headless Chromium
CSS = """body{font-family:"DejaVu Serif",serif;font-size:10.5pt;line-height:1.35}
h1{font-size:14pt;color:#1f3b63;margin-top:18pt} h1.title{font-size:16pt;text-align:center}
p.subtitle,p.author{text-align:center}
table{border-collapse:collapse;font-size:8.5pt;margin:6pt 0;width:100%}
th,td{border-bottom:0.5pt solid #999;padding:2pt 4pt;text-align:left;vertical-align:top}
th{border-bottom:1pt solid #000} img{max-width:100%} figcaption{font-size:9pt}
blockquote{border-left:3pt solid #b00020;margin-left:0;padding-left:8pt}
@page{size:A4;margin:18mm}"""
chromium = shutil.which("chromium") or "/opt/pw-browsers/chromium"
if Path(chromium).exists():
    with tempfile.TemporaryDirectory() as td:
        css = Path(td) / "s.css"
        css.write_text(CSS)
        html = Path(td) / "s1.html"
        subprocess.run(["pandoc", OUT.name, "-s", "--mathml",
                        "--embed-resources", "--resource-path=.",
                        "-c", str(css), "-o", str(html)], cwd=OUT.parent,
                       check=True)
        subprocess.run([chromium, "--headless", "--no-sandbox", "--disable-gpu",
                        "--no-pdf-header-footer",
                        f"--print-to-pdf={OUT.with_suffix('.pdf')}", str(html)],
                       check=True, capture_output=True)
    print("wrote", OUT.with_suffix(".pdf"))
