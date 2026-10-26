# DT-GCN: Digital-Twin-Assisted Graph Learning for Stealthy FDIA Detection & Localization

Code, results and the paper for an ICISPES 2027 submission
(track: *AI/ML Applications, Digital Twins, Data Analytics, IoT, and Cybersecurity in Power & Energy Systems*).

**Paper:** `paper/main.pdf` (IEEE conference format; LaTeX source in `paper/`).

## What is in here

| Path | Content |
|---|---|
| `fdia_dt/grid.py` | AC measurement model h(x), analytic Jacobian, Newton–Raphson power flow, WLS state estimator |
| `fdia_dt/data.py` | 15-min load/PV scenarios, physics-based digital twin, AC-stealthy and naive FDIA generator, twin-residual features |
| `fdia_dt/models.py` | Dual-head GCN (detection + localization), MLP, training loop |
| `fdia_dt/run_experiments.py` | Benchmark of 8 detectors on IEEE 14/118-bus systems (5 seeds) |
| `fdia_dt/sensitivity.py` | Sensitivity to digital-twin forecast error σ_f |
| `fdia_dt/make_tables.py`, `fdia_dt/figures.py` | LaTeX tables and figures for the paper |
| `results/` | CSV results, logs, `digest.json` with every number used in the paper |
| `figures/` | PDF/PNG figures |

## Reproduce

```bash
python3 -m venv .venv && . .venv/bin/activate
pip install -r requirements.txt
./run_all.sh            # ~1.5 h on a 4-core CPU
cd paper && pdflatex main && bibtex main && pdflatex main && pdflatex main
```

Datasets (`data/*.pkl`) are regenerated deterministically by `run_all.sh` (seeds fixed).
