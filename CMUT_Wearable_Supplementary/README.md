# Supplementary Material — Wearable Phased-Array CMUT Mechanostimulation Platform


This package does three things:

1. **Reproduction.** It re-implements, independently, every in-silico experiment
   in the manuscript: CMUT cell, focused field, steering, grating lobes,
   apodization, layered tissue, Pennes bioheat, MI/TI, radiation force,
   waveform, closed-loop control and the dose-response model.
2. **Hardware digital twin.** It runs a Monte Carlo emulation of the bench
   results a fabricated patch would be expected to give, covering fabrication
   spread, electronics errors, body curvature and measurement uncertainty.
   **These are emulated, not measured** (see `supplementary/S3_Emulated_Data_Disclosure.md`).


## Folder layout

```
CMUT_Wearable_Supplementary/
├── README.md                     this file
├── run_all.py                    one command reproduces every result
├── make_supplementary_doc.py     builds S1 (.md/.docx/.pdf) from the results
├── requirements.txt  LICENSE  CITATION.cff
├── cmutsim/                      model library
│   ├── params.py                 Tables I–II + stated assumptions
│   ├── cmut_cell.py              collapse voltage, resonance, BVD model
│   ├── field.py                  Rayleigh–Sommerfeld field, layered paths, apodization
│   ├── bioeffects.py             intensity, MI, TIS/TIB, 3-D Pennes, Kelvin displacement
│   ├── control.py                FOPDT plant + PI with anti-windup
│   └── hardware_twin.py          Monte Carlo hardware / measurement emulation
├── results/
│   ├── figures/                  Fig_S1–S12  (PNG 600 dpi + vector PDF)
│   ├── tables/                   T_S0–T_S12  (CSV + LaTeX)
│   ├── data/                     raw numerical data (CSV, NPZ field maps)
│   ├── emulated_measurements/    Fig_E*, T_E*, raw twin outputs (EMULATED)
│   └── summary.json              every headline number
└── supplementary/
    ├── S1_Supplementary_Methods_and_Results.md / .docx / .pdf   (upload this)
    ├── S2_Manuscript_Consistency_Report.md
    ├── S3_Emulated_Data_Disclosure.md
    ├── S4_Phase_II_Bench_Protocol.md
    └── S5_Data_and_Code_Availability.md
```

## Reproduce

```bash
pip install -r requirements.txt
python run_all.py            # full resolution (~3–5 min)
python run_all.py --quick    # coarse grids, ~1.5 min (for checking only)
python make_supplementary_doc.py   # rebuild the Supplementary Information file
```

All random draws use the fixed seed `20261019`, so reruns give identical
numbers.

## Figure and table map

| Supplementary item | Manuscript counterpart | Content |
|---|---|---|
| Table S0 | Table XII | Manuscript vs reproduced vs emulated hardware |
| Fig S1 / Table S1 | Fig 3–4, Table I | CMUT element impedance, TX response, resonance vs radius |
| Fig S2 / Table S2 | Fig 7, Table III | Focal field, -6 dB and -3 dB widths, drive and power |
| Fig S3 / Table S3 | Fig 8, Table IV | Steering ±20° |
| Fig S4 / Table S4 | Fig 9, Table V | Grating lobes vs pitch (arc pattern) |
| Fig S5 / Table S5 | Fig 10, Table VI | Apodization |
| Fig S6 / Table S6 | Fig 11–12, Table II | Layered tissue, aberration, insertion loss |
| Fig S7 / Table S7 | Fig 13–14, Table VII | 3-D Pennes, TIS/TIB |
| Fig S8 / Table S8 | Fig 15, Table VIII | MI envelope |
| Fig S9 / Table S9 | Fig 16, Table IX | Radiation force, Kelvin displacement |
| Fig S10 | Fig 17 | Burst waveform and spectrum |
| Fig S11 / Table S10 | Fig 18, Table X | Closed loop + robustness Monte Carlo |
| Fig S12 / Table S11 | Fig 21 | Dose-response hypothesis + Phase III sample size |
| Table S12 | Section II (">4 h") | Wearable power budget |
| Fig/Table E1–E5 | — (new) | Digital twin: fabrication, beam MC, hydrophone, thermometry, 4 h session |
