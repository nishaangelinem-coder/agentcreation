# Spectrum-native proteomics benchmark

Code, results, figures and manuscript for:

> **From identification-first to spectrum-native proteomics: a controlled benchmark of transformer models that learn directly from MS1 and MS2 signal.**
> M. Nisha Angeline (Velalar College of Engineering and Technology)

Everything runs on a CPU (PyTorch, NumPy, SciPy, scikit-learn, matplotlib).

**Note:** the main benchmark is simulated by the included *virtual instrument*. It is validated on real data: PRIDE-derived BSA LC-MS/MS runs (bigbio/quantms CI data) and MassIVE-KB annotated spectra (Casanovo sample data). The PRIDE archive itself was not reachable from the compute environment.

## Layout

| Path | Contents |
|---|---|
| `src/chem.py` | monoisotopic masses, token vocabulary, b/y ladders |
| `src/simulator.py` | virtual instrument: proteomes, tryptic digestion, MS2 fragmentation model, LC-MS1 cohorts, ID-free alignment |
| `src/models.py` | `SpecNovo` (de novo transformer: mass-aware + ion-evidence), `FeatureFormer` (ProFormer-style MS1 transformer + SSL), `MS1Former` (image ViT, exploratory only) |
| `src/baselines.py` | target–decoy database search (hyperscore-like), spectrum-graph DP de novo |
| `src/metrics.py` | mass-aligned AA / peptide metrics, precision–coverage |
| `src/data_gen.py` | builds the 250k-spectrum training corpus and the 11 test sets (`data/ms2_corpus.pkl`, ~160 MB, not committed) |
| `src/train_denovo.py`, `src/eval_denovo.py` | SpecNovo training / full evaluation |
| `src/pheno2.py`, `src/run_pheno2.py` | spectrum-to-phenotype benchmark (main + dark-proteome sweep) |
| `src/fig_*.py`, `src/make_tables.py`, `src/build_paper.py`, `src/style_docx.py` | figures, tables, manuscript build |
| `src/real_*.py`, `src/data_gen_calibrated.py`, `src/fig_real.py` | real-data validation: search, calibration, de novo transfer, MS1 bottleneck, Fig. 7 |
| `realdata/` | real-data results (JSON, PSMs, features) |
| `results/` | JSON metrics, training logs, model checkpoints (`*.pt`), analysis pickles |
| `figures/` | Fig. 1–7 (300 dpi PNG) |
| `paper/` | `manuscript.pdf`, `manuscript.docx`, `manuscript.tex`, `manuscript.md`, `tables.md` |

## Reproduce

Run each command from this directory.

1. Build the corpus:
   ```bash
   PYTHONPATH=src python3 src/data_gen.py
   ```
2. Train the base model:
   ```bash
   PYTHONPATH=src python3 src/train_denovo.py --steps 3500 --tag main
   ```
3. Train the mass-aware model:
   ```bash
   PYTHONPATH=src python3 src/train_denovo.py --steps 3000 --save_at 1000,2000 --mass_aware 1 --tag ma
   ```
4. Train the ion-evidence model, starting from the mass-aware checkpoint:
   ```bash
   PYTHONPATH=src python3 src/train_denovo.py --steps 1200 --save_at 400,800 --mass_aware 1 --evidence 1 --init results/specnovo_ma_s2000.pt --lr 3e-4 --tag ev
   ```
5. Run the full evaluation of the final model:
   ```bash
   PYTHONPATH=src python3 src/eval_denovo.py --ckpt results/specnovo_ev.pt --tag ev
   ```
6. Run the ablation evaluations:
   ```bash
   PYTHONPATH=src python3 src/eval_denovo.py --ckpt results/specnovo_main.pt --tag base --full 0 --mass_aware 0 --evidence 0
   ```
   ```bash
   PYTHONPATH=src python3 src/eval_denovo.py --ckpt results/specnovo_ma_s2000.pt --tag ma --full 0 --evidence 0
   ```
   ```bash
   PYTHONPATH=src python3 src/eval_denovo.py --ckpt results/specnovo_ev_s400.pt --tag ev400 --full 0
   ```
   ```bash
   PYTHONPATH=src python3 src/eval_denovo.py --ckpt results/specnovo_ev_s800.pt --tag ev800 --full 0
   ```
7. Run the phenotype benchmark. Repeat each command for seeds 0, 1 and 2:
   ```bash
   PYTHONPATH=src python3 src/run_pheno2.py --mode main --seed 0
   ```
   ```bash
   PYTHONPATH=src python3 src/run_pheno2.py --mode sweep --seed 0 --fracs 0.0,0.25,0.75
   ```
8. Build the figures:
   ```bash
   for f in fig_schematic fig_data fig_denovo fig_pheno; do PYTHONPATH=src python3 src/$f.py; done
   ```
9. Build the tables and the paper:
   ```bash
   PYTHONPATH=src python3 src/make_tables.py && PYTHONPATH=src python3 src/build_paper.py && python3 src/style_docx.py
   ```

## Real-data validation

All commands run from this directory.

1. Get the real data:
   ```bash
   git clone --depth 1 https://github.com/bigbio/quantms-test-datasets /home/user/bigbio/quantms-test-datasets
   ```
   ```bash
   git clone --depth 1 https://github.com/Noble-Lab/casanovo /home/user/bigbio/casanovo
   ```
2. Install the mzML parsing dependencies:
   ```bash
   pip install psims lxml
   ```
3. Run the target–decoy search of the PRIDE BSA runs:
   ```bash
   PYTHONPATH=src python3 src/real_search.py
   ```
4. Calibrate the simulator against the real spectra:
   ```bash
   PYTHONPATH=src python3 src/real_calib.py
   ```
5. Run the zero-shot de novo test on real spectra:
   ```bash
   PYTHONPATH=src python3 src/real_denovo.py --tag zeroshot
   ```
6. Generate the calibrated simulation corpus:
   ```bash
   PYTHONPATH=src python3 src/data_gen_calibrated.py
   ```
7. Fine-tune SpecNovo on the calibrated corpus:
   ```bash
   PYTHONPATH=src python3 src/train_denovo.py --steps 800 --save_at 400 --mass_aware 1 --evidence 1 --init results/specnovo_ev.pt --lr 1e-4 --warmup 50 --corpus data/ms2_calibrated.pkl --n_train 80000 --tag cal
   ```
8. Re-test the fine-tuned model on real spectra:
   ```bash
   PYTHONPATH=src python3 src/real_denovo.py --ckpt results/specnovo_cal.pt --tag calibrated
   ```
9. Compute bootstrap confidence intervals, the MS1 bottleneck and Fig. 7:
   ```bash
   PYTHONPATH=src python3 src/real_bootstrap.py && PYTHONPATH=src python3 src/real_ms1.py && PYTHONPATH=src python3 src/fig_real.py
   ```

| Real-data measurement | Value |
|---|---|
| PRIDE BSA: PSMs at 1 % FDR / MS2 spectra | 64 / 3,135 |
| MS1 features identified | 48 / 25,613 (0.19 %) |
| Complete profiles, ID-free vs ID-first | 1,342 vs 1 |
| Real HCD de novo peptide recall (n = 119): spectrum graph | 0.218 |
| SpecNovo zero-shot | 0.076 |
| SpecNovo calibrated | 0.151 |

## Headline measurements

| Measurement | Value |
|---|---|
| DB search recall, canonical / variant / novel peptides | 0.999 / 0.000 / 0.000 |
| SpecNovo peptide recall, canonical / variant / novel | 0.416 / 0.389 / 0.387 |
| SpecNovo vs spectrum-graph precision–coverage AUC | 0.663 vs 0.630 |
| Ion-evidence ablation: peptide recall | 0.005 → 0.416 |
| Ion-evidence ablation: token accuracy | 0.582 → 0.864 |
| Phenotype AUROC, n = 192, 50 % dark biomarkers: ID-first / ID+MBR / ID-free / FeatureFormer | 0.648 / 0.932 / 0.975 / 0.767 |
| Unseen-site AUROC, FeatureFormer from scratch vs with SSL | 0.651 vs 0.728 |
