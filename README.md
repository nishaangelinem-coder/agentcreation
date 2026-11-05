# siwlab — open-source design environment for a low-profile SIW slot array with characteristic-mode-guided decoupling (28 GHz)

This repository contains a self-contained, open-source (Python/NumPy/SciPy/Numba) electromagnetic design
environment and the manuscript *"Low-Profile SIW Slot Array With Characteristic-Mode-Guided Decoupling for
28 GHz 5G mmWave Front Ends"* prepared in the Elsevier `elsarticle` format for
*AEU – International Journal of Electronics and Communications*.

## What is inside

| Path | Content |
|---|---|
| `siwlab/siw.py` | SIW design equations: effective width (Cassivi, Xu–Wu), cutoff, TE10 dispersion, dielectric/conductor attenuation, design rules |
| `siwlab/slots.py` | Standing-wave longitudinal-slot array synthesis (Elliott first-order procedure with full-wave slot data) |
| `siwlab/mesh.py` | Structured triangular mesher for metal plates with rectangular apertures |
| `siwlab/mom.py` | RWG electric-field integral-equation MoM (Numba) and characteristic-mode solver (X J = λ R J) |
| `siwlab/fdtd.py` | Numba-accelerated 3-D FDTD: nonuniform Yee grid, CPML, TE10 SIW ports, aperture near-to-far-field |
| `siwlab/geometry.py` | Parametric builder of the dual-port 2×(1×4) SIW slot array (vias, slots, slits, ports) |
| `scripts/01_validate_throughline.py` | FDTD validation on an SIW through line (dispersion vs. analytical model) |
| `scripts/02_slot_characterization.py`, `02b_synthesis.py` | Single-slot admittance extraction and array synthesis |
| `scripts/03a_cma_calibration.py` | Aperture-mode CMA: effective-medium calibration and slit design curve |
| `scripts/03b_cma_dual_plate.py` | Electric-current CMA of the finite shared top wall (gap-crossing index) |
| `scripts/04_aperture_cma_array.py` | Aperture-system MoM/CMA of the two sub-arrays with parasitic slits (Schur-complement coupling metric) |
| `scripts/05_array_fdtd.py` | Full dual-port array: S-parameters, far field, gain, efficiency, ECC |
| `scripts/06_figures.py` | Publication figures |
| `results/` | Raw and post-processed data (JSON/NPZ) produced by the scripts |
| `paper/` | LaTeX source (`elsarticle`), bibliography, figures and compiled PDF |

## Reproducing

```bash
pip install -r requirements.txt
python scripts/01_validate_throughline.py
python scripts/02_slot_characterization.py          # ~2.5 min per case, 4 cores
python scripts/02b_synthesis.py
python scripts/03a_cma_calibration.py 4.5
python scripts/03b_cma_dual_plate.py ref 7.0 0
python scripts/03b_cma_dual_plate.py slits 7.0 <L_slit>
python scripts/04_aperture_cma_array.py
python scripts/05_array_fdtd.py reference 7.0 0
python scripts/05_array_fdtd.py nodecouple 7.0 0
python scripts/05_array_fdtd.py decoupled 7.0 <L_slit>
python scripts/06_figures.py
python scripts/07_numbers.py
python -m pytest tests -q
cd paper && latexmk -pdf main.tex
```

All solvers are written from scratch; no commercial EM software is required. The FDTD was verified against the
analytical SIW dispersion relation and the MoM/CMA against the resonance of a thin strip dipole.
