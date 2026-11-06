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

## Main results (FDTD, port 1 excited, 28 GHz)

| Configuration | $|S_{11}|$ | $|S_{21}|$ | worst $|S_{21}|$ 27.5–28.5 GHz | realised gain |
|---|---|---|---|---|
| D = 7.0 mm, no decoupling structure | −9.0 dB | −15.3 dB | −14.9 dB | 10.9 dBi |
| D = 7.0 mm, top-wall slits 3.6…4.8 mm or continuous | −9…−10 dB | −14.8…−13.9 dB | — | — |
| D = 9.4 mm, parasitic cavities, no slots | −9.1 dB | −17.5 dB | −17.2 dB | 10.1 dBi |
| D = 9.4 mm, CM-guided offset cavity-backed slots (4.0 mm, 0.8 mm offset) | −8.2 dB | **−23.5 dB** | **−20.5 dB** | **11.9 dBi** |

Key physical finding: slits etched in the top wall over the sub-cut-off channel between adjacent SIWs cannot
resonate (the channel acts as a short at 0.05 λ0 depth) and therefore do not decouple the sub-arrays, whatever
their length; a parasitic aperture must be backed by a cavity/open region, which the aperture-mode CMA makes
explicit. See `paper/main.pdf`.

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
python scripts/05_array_fdtd.py reference D=7.0
python scripts/05_array_fdtd.py nodecouple D=7.0
python scripts/05_array_fdtd.py decoupled D=7.0 Ls=4.0          # top-wall slits (negative result)
python scripts/05_array_fdtd.py reference D=9.4
python scripts/05_array_fdtd.py D94_base D=9.4 cav=1
python scripts/05_array_fdtd.py D94_o40 D=9.4 cav=1 Ls=4.0 xoff=0.8   # final design
python scripts/06_figures.py
python scripts/07_numbers.py
python -m pytest tests -q
cd paper && latexmk -pdf main.tex
```

All solvers are written from scratch; no commercial EM software is required. The FDTD was verified against the
analytical SIW dispersion relation and the MoM/CMA against the resonance of a thin strip dipole.
