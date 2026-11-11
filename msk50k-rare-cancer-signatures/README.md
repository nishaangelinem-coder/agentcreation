# Mutational signatures in three rare cancers from MSK-IMPACT 50K

Pipeline, results and manuscript for a mutational-signature analysis of angiosarcoma (ASI), thymic epithelial
tumours (THYM) and adrenocortical carcinoma (ACC) in the public MSK-IMPACT 50K cohort (54,331 tumours;
cBioPortal study `msk_impact_50k_2026`, Bandlamudi et al., Cancer Cell 2026).

## Layout

| Path | Content |
|---|---|
| `env/environment.yml` | Conda environment (Python 3.11, SigProfiler suite, lifelines, bedtools) |
| `env/environment.lock.yml` | Exact package versions used |
| `code/01_acquire.py` … `10_figures.py` | Analysis phases, run in order (`run_all.sh`) |
| `results/` | Aggregate outputs (JSON/TSV): provenance, cohort flow, gates, statistics |
| `figures/` | Manuscript figures (300 dpi PNG) |
| `manuscript/` | Springer-style manuscript (`.docx`) and its build script |

## Reproduce

```bash
micromamba create -f env/environment.yml
```

```bash
micromamba run -n msk50k ./run_all.sh
```

Seeds are fixed at 42. Raw data are downloaded by `01_acquire.py`; SHA-256 checksums are written to
`results/provenance.json`.

## Data licence

The MSK-IMPACT 50K data are released under **CC BY-NC-ND 4.0**. Raw data and tumour-level derived tables are
therefore not stored in this repository; they are regenerated locally by the scripts. Cite the original study
and cBioPortal when using the data.

## Quality gates (pre-specified)

| Gate | Outcome |
|---|---|
| ≥50 analysable tumours per group | Met for ASI (82); not met for ACC (48) or THYM (23) |
| Simulation cosine ≥0.90 | Met only at ≥50 SNVs per tumour; cohort median is 2 |
| Mean pairwise Cohen's κ ≥0.60 | 0.585, not met |
| Proportional hazards | Violated in the pooled model; satisfied within groups |
