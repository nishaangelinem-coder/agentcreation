# Microplastics & Neuroinflammation — In-Silico Pipeline

Computational companion to the systematic review *"Microplastics and
Neuroinflammation: Does Polymer Type Matter?"* It tests whether a polymer's
predicted neuroinflammatory potential tracks its **leachate chemistry** acting on
neuroinflammation receptors (AhR, TLR4–MD2, NLRP3) rather than its inert backbone
or particle size.

> **Framing.** Every ranking here is a **falsifiable computational prediction**
> for wet-lab testing, not established fact. Binding affinity is strength, **not**
> activation direction. Virtual-cell kinetic parameters are **illustrative**, not
> fitted. Nothing is fabricated: steps that could not run in the execution
> environment are reported as such (see `results/summary.md` and
> `results/docking_status.md`).

## Pipeline

```
polymer → leachate SMILES (S1) → Tox21 ML triage (S2) → receptor docking (S3)
        → virtual microglia ODE cell (S4) → ranking + figures (S5)
```

| Stage | Script | Status in this run |
|---|---|---|
| S1 leachate inventory | `src/s1_leachates.py` | ✅ 56 RDKit-validated compounds, ≥10/polymer |
| S2 Tox21 ML triage | `src/s2_ml_triage.py` | ✅ scaffold-split RandomForest; real ROC-AUC/AUPRC |
| S3 receptor docking | `src/s3_docking.py` | ⚠️ not run — PDB/AlphaFold/UniProt unreachable (no fabricated affinities) |
| S4 virtual microglia cell | `src/s4_virtual_cell.py` + `models/microglia.py` | ✅ ODE index (illustrative kinetics) |
| S5 ranking & figures | `src/s5_rank_and_plot.py` | ✅ 3 figures + `summary.md` |

## Reproduce

```bash
pip install -r requirements.txt     # or: conda env create -f environment.yml
python src/s1_leachates.py
python src/s2_ml_triage.py          # downloads+caches Tox21 to data/
python src/s3_docking.py            # honest no-op without receptors/Vina
python src/s4_virtual_cell.py
python src/s5_rank_and_plot.py
```

## Headline result (this run)

Predicted neuroinflammatory potential: **PU > PVC > PE > PS > PP**
(PU driven by aromatic-amine hydrolysis products TDA/MDA; PVC by BPA/phthalates).
Scaffold-split triage: NR-AhR ROC-AUC 0.84 / AUPRC 0.55. See `results/summary.md`.

## Key outputs

- `data/leachates.csv`, `data/experimental_ranking.csv` (author fills)
- `results/ml_scores.csv`, `results/ml_metrics.csv`, `results/ml_calibration.md`
- `results/docking.csv` (+ `docking_status.md`)
- `results/neuroinflammation_index.csv`, `results/summary.md`, `results/METHODS_PARAGRAPH.md`
- `results/figures/fig1..3.{png,svg}`
- `paper/` — IEEE-format manuscript (`.tex` + readable `.md`)

## To enable Stage-3 docking

Place `data/receptors/{AhR,TLR4_MD2,NLRP3}.pdbqt`, `pip install vina meeko`, and
rerun `src/s3_docking.py`; S4 then blends the docking term automatically.
