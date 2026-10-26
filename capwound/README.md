# capwound — Wideband bioimpedance & electrochemistry of CAP-treated infected wound phantoms

An open, reproducible **in silico** framework, with a Springer Nature–format manuscript, for

> *Wideband Bioimpedance and Electrochemical Characterization of Cold Atmospheric
> Plasma-Treated Infected Wound Phantoms: A Multiphysics In Silico Study*

> **Important:** every result in `results/`, `figures/` and the manuscript comes from a
> simulation (virtual phantoms). Nothing here is a wet-lab measurement. Before a Q1
> submission the model needs experimental calibration and validation; the pipeline
> already accepts measured spectra (see *Using real data*).

## What is modelled

| Module | Physics |
|---|---|
| `plasma_chemistry.py` | 1-D implicit diffusion–reaction of H2O2, NO2⁻, NO3⁻, H⁺ in biofilm + agarose wound bed; catalase sink; phosphate-buffer pH |
| `microbiology.py` | depth-resolved Weibull (Mafart) inactivation; membrane permeabilisation |
| `dielectric.py` | single-shell bacterium, Maxwell–Garnett biofilm mixture, Debye gel |
| `impedance.py` | interdigitated electrode, 10 Hz–100 MHz, electrode polarisation CPE, stray C, noise |
| `electrochem.py` | SPCE: Randles EIS, finite-difference cyclic voltammetry (Butler–Volmer + iR), Prussian-blue H2O2 amperometry |
| `fitting.py` | CNLS Cole–Cole+CPE and Randles fits, Lin-KK test (Schönleber 2014) |
| `stats.py` | ANOVA/Tukey, Spearman, nested-CV PLS + VIP |
| `pipeline.py` | 9 groups × n = 6 virtual phantoms with phantom-to-phantom variability |
| `io.py` | load and fit **measured** spectra from CSV |

## Quick start

```bash
pip install -r requirements.txt
python -m pytest -q tests          # 8 verification tests
python scripts/run_all.py          # ~1 min: results/ + figures/ (PDF + 600-dpi PNG)
cd manuscript && ./build.sh        # compiles main.pdf (needs a TeX Live with latexmk)
```

## Manuscript (Springer Nature `sn-jnl`, formatted for *Medical & Biological Engineering & Computing*)

See `manuscript/MBEC_submission_checklist.md` for the journal-specific items.

`manuscript/main.tex` is written for the official Springer Nature LaTeX template
(`\documentclass[pdflatex,sn-mathphys-num]{sn-jnl}`). The class file is not
redistributed here. Download the template from Springer Nature's LaTeX author support page and copy
`sn-jnl.cls` and `sn-mathphys-num.bst` into `manuscript/`, and the official class is used
automatically. Without it the draft compiles with `sn-jnl-compat.cls`, a local preview shim
that is **not** the official class.

## Using real data

Put one CSV per spectrum (`freq_Hz,Zreal_ohm,Zimag_ohm`) in a folder together with a
`samples.csv` (`file,group,replicate,log10_cfu`), then:

```python
from capwound.io import analyse_folder
df = analyse_folder("data/measured", K=<your cell constant>)
```

## Key simulated findings (seed 2026, n = 6)

* 30–300 s CAP → 0.31–2.11 log10 reduction; gas flow alone has no effect.
* β-dispersion Δε 182 → 78, f_c 16.6 → 27.1 MHz; f_c is the best single marker (η² = 0.973).
* Sterile phantoms show a 5 % ionic/evaporative conductivity artefact. Use dispersion parameters, not |Z| alone.
* Biofilm H2O2 is 16× lower than in sterile phantoms (catalase scavenging).
* R_ct falls 40 % as CAP unblocks the electrode; nested PLS RMSECV = 0.52 log.
