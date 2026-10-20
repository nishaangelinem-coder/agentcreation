# Stage 3 — docking status

**Docking was NOT executed. No affinities were fabricated.**

## Why

- Receptor structures missing: AhR, TLR4-MD2, NLRP3. RCSB PDB (files.rcsb.org), AlphaFold DB (alphafold.ebi.ac.uk) and UniProt (rest.uniprot.org) were all unreachable from the execution sandbox (egress proxy returned 403/no-tunnel), so they could not be fetched automatically. Stage files locally into data/receptors/ to run.
- AutoDock Vina Python bindings (`pip install vina`) are not installed in the sandbox.

## What a full run would do

With the receptor structures and AutoDock Vina present, this stage docks the Stage-2-prioritised leachates into each receptor and records the best affinity per (ligand, receptor), **paired with the Stage-2 Tox21 label** so that binding strength is never confused with activation direction.

## Receptor targets and boxes (documented for reproducibility)

| receptor | structure source | box size (A) |
|---|---|---|
| AhR | AlphaFold:P35869 (PAS-B LBD) | 22.5x22.5x22.5 |
| TLR4-MD2 | PDB:3FXI (TLR4-MD2-LPS) | 25.0x25.0x25.0 |
| NLRP3 | PDB:6NPY (NACHT/LRR) | 24.0x24.0x24.0 |

## Consequence for the pipeline

Stage 4 falls back to the Stage-2 Tox21 ML signal as the receptor-engagement proxy (NR-AhR is a *direct* Tox21 endpoint for the AhR axis; SR-ARE and SR-p53 cover oxidative/genotoxic stress). TLR4-MD2 and NLRP3 docking remain **pending structure access** and are reported as a limitation, not as a result.
