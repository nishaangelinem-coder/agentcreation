# Methods paragraph (for the manuscript)

**In-silico neuroinflammation screening pipeline.** For each polymer (PE, PP, PS,
PVC, PU) a curated inventory of 56 documented additives, residual monomers and
commonly adsorbed contaminants (11–12 per polymer) was assembled with canonical
SMILES validated and standardised in RDKit (2025.x); every structure parsed
cleanly and no entry was discarded. Toxicological triage used the public Tox21
benchmark (7,823 structures after RDKit featurisation) with three
inflammation/stress-relevant endpoints — NR-AhR (aryl hydrocarbon receptor),
SR-ARE (Nrf2/antioxidant response) and SR-p53 (genotoxic stress). Each endpoint
was modelled with a RandomForest classifier (300 trees, balanced class weights)
on 2,048-bit Morgan fingerprints (radius 2), evaluated under a deterministic
Bemis–Murcko **scaffold split** (80/20). Because the positive classes are rare,
performance is reported as both ROC-AUC and the imbalance-aware AUPRC: NR-AhR
ROC-AUC 0.84 / AUPRC 0.55 (test base rate 0.14), SR-ARE 0.72 / 0.47 (0.24), and
SR-p53 0.76 / 0.35 (0.11). Models refit on all labelled Tox21 rows then scored
every leachate. Receptor docking (AutoDock Vina) into AhR (AlphaFold P35869),
TLR4–MD2 (PDB 3FXI) and NLRP3 (PDB 6NPY), with documented 22.5–25 Å boxes, was
implemented but **not executed** in this run because RCSB PDB, AlphaFold DB and
UniProt were unreachable from the execution environment; no affinities were
fabricated, and the receptor-engagement signal therefore derives from the Tox21
triage (NR-AhR directly for the AhR axis; SR-ARE/SR-p53 for oxidative and
genotoxic stress). Per polymer, the top-3 mean probability per endpoint was
combined (weights NR-AhR 0.50, SR-ARE 0.30, SR-p53 0.20) into a receptor drive
*u* ∈ [0,1] feeding a minimal microglial ODE model (receptor → NF-κB → {TNF-α,
IL-6}; **illustrative, non-fitted kinetics**) integrated to steady state
(SciPy LSODA); the normalised steady-state cytokine sum is the neuroinflammation
index. Each result CSV is reproducible from a single `src/` script.

**Reproducibility parameters.** Fingerprint: Morgan r=2, 2048 bits. Split:
deterministic scaffold, train fraction 0.80, seed 42. Classifier:
RandomForest(n_estimators=300, class_weight="balanced"). Endpoint weights:
{NR-AhR:0.50, SR-ARE:0.30, SR-p53:0.20}; aggregation: top-3 mean. ODE drive
u∈[0,1]; integrator LSODA, rtol 1e-8, atol 1e-10, t_end 100. Docking boxes:
AhR 22.5³ Å, TLR4–MD2 25³ Å, NLRP3 24³ Å (pending structure access).
