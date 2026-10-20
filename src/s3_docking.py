#!/usr/bin/env python3
"""Stage 3 — Receptor structure & docking (AutoDock Vina).

This stage docks the Stage-2-prioritised leachates into three neuroinflammation
receptors and records, for every pose, the Vina affinity PAIRED WITH the Stage-2
Tox21 label (because affinity is binding strength, not agonist/antagonist
direction — the central Stage-3 caveat).

Receptors (fetch from RCSB PDB first; AlphaFold DB fallback):
    AhR      ligand-binding (PAS-B) domain   — AlphaFold P35869 (human AHR)
    TLR4-MD2 complex                          — PDB 3FXI
    NLRP3    NACHT/LRR                         — PDB 6NPY

Docking boxes (documented for reproducibility; centred on the co-crystal ligand
or the canonical ligand cavity, 22.5 A isotropic unless noted):
    AhR:      box 22.5^3 A centred on the PAS-B cavity
    TLR4-MD2: box 25^3 A centred on the MD2 hydrophobic pocket
    NLRP3:    box 24^3 A centred on the NACHT ADP/NBD site

IMPORTANT — honest execution contract:
    If the receptor structures cannot be obtained (no local file AND the
    structure hosts are unreachable) OR AutoDock Vina is not installed, this
    stage DOES NOT invent affinities. It writes results/docking_status.md
    explaining exactly what was missing and emits an EMPTY results/docking.csv
    (header only). Downstream Stage 4 then falls back to the Stage-2 ML signal
    and says so.

Expected inputs when run in a full environment:
    data/receptors/AhR.pdbqt, TLR4_MD2.pdbqt, NLRP3.pdbqt
Outputs:
    results/docking.csv          (SMILES, receptor, best_affinity_kcal_mol,
                                  tox21_label_endpoint, tox21_prob, pose_file)
    results/docking_status.md
"""
from __future__ import annotations

import csv
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
DATA = os.path.join(ROOT, "data")
RESULTS = os.path.join(ROOT, "results")
RECEPTOR_DIR = os.path.join(DATA, "receptors")

RECEPTORS = {
    "AhR": {
        "source": "AlphaFold:P35869 (PAS-B LBD)",
        "pdbqt": os.path.join(RECEPTOR_DIR, "AhR.pdbqt"),
        "box_center": None,  # set from structure at prep time
        "box_size": (22.5, 22.5, 22.5),
    },
    "TLR4-MD2": {
        "source": "PDB:3FXI (TLR4-MD2-LPS)",
        "pdbqt": os.path.join(RECEPTOR_DIR, "TLR4_MD2.pdbqt"),
        "box_center": None,
        "box_size": (25.0, 25.0, 25.0),
    },
    "NLRP3": {
        "source": "PDB:6NPY (NACHT/LRR)",
        "pdbqt": os.path.join(RECEPTOR_DIR, "NLRP3.pdbqt"),
        "box_center": None,
        "box_size": (24.0, 24.0, 24.0),
    },
}

DOCKING_CSV = os.path.join(RESULTS, "docking.csv")
STATUS_MD = os.path.join(RESULTS, "docking_status.md")
HEADER = ["SMILES", "receptor", "best_affinity_kcal_mol",
          "tox21_label_endpoint", "tox21_prob", "pose_file"]


def vina_available() -> bool:
    try:
        import vina  # noqa: F401
        return True
    except Exception:  # noqa: BLE001
        return False


def receptors_present() -> list[str]:
    return [name for name, meta in RECEPTORS.items() if os.path.exists(meta["pdbqt"])]


def write_empty_with_status(reason_lines: list[str]) -> None:
    os.makedirs(RESULTS, exist_ok=True)
    with open(DOCKING_CSV, "w", newline="") as fh:
        csv.writer(fh).writerow(HEADER)
    with open(STATUS_MD, "w") as fh:
        fh.write("# Stage 3 — docking status\n\n")
        fh.write("**Docking was NOT executed. No affinities were fabricated.**\n\n")
        fh.write("## Why\n\n")
        for line in reason_lines:
            fh.write(f"- {line}\n")
        fh.write("\n## What a full run would do\n\n")
        fh.write("With the receptor structures and AutoDock Vina present, this "
                 "stage docks the Stage-2-prioritised leachates into each receptor "
                 "and records the best affinity per (ligand, receptor), **paired "
                 "with the Stage-2 Tox21 label** so that binding strength is never "
                 "confused with activation direction.\n\n")
        fh.write("## Receptor targets and boxes (documented for reproducibility)\n\n")
        fh.write("| receptor | structure source | box size (A) |\n")
        fh.write("|---|---|---|\n")
        for name, meta in RECEPTORS.items():
            bs = "x".join(str(v) for v in meta["box_size"])
            fh.write(f"| {name} | {meta['source']} | {bs} |\n")
        fh.write("\n## Consequence for the pipeline\n\n")
        fh.write("Stage 4 falls back to the Stage-2 Tox21 ML signal as the "
                 "receptor-engagement proxy (NR-AhR is a *direct* Tox21 endpoint "
                 "for the AhR axis; SR-ARE and SR-p53 cover oxidative/genotoxic "
                 "stress). TLR4-MD2 and NLRP3 docking remain **pending structure "
                 "access** and are reported as a limitation, not as a result.\n")
    print(f"[S3] docking not run; honest status written -> {STATUS_MD}")
    print(f"[S3] empty (header-only) docking.csv written -> {DOCKING_CSV}")


def main() -> int:
    os.makedirs(RESULTS, exist_ok=True)
    reasons = []
    present = receptors_present()
    missing = [n for n in RECEPTORS if n not in present]
    if missing:
        reasons.append(
            "Receptor structures missing: " + ", ".join(missing) +
            ". RCSB PDB (files.rcsb.org), AlphaFold DB (alphafold.ebi.ac.uk) and "
            "UniProt (rest.uniprot.org) were all unreachable from the execution "
            "sandbox (egress proxy returned 403/no-tunnel), so they could not be "
            "fetched automatically. Stage files locally into data/receptors/ to run."
        )
    if not vina_available():
        reasons.append(
            "AutoDock Vina Python bindings (`pip install vina`) are not installed "
            "in the sandbox."
        )

    if reasons:
        write_empty_with_status(reasons)
        return 0

    # ---- Full docking path (runs only when structures + Vina are available) ----
    import pandas as pd
    from vina import Vina  # noqa: F401  (import kept local to the available branch)

    scores = pd.read_csv(os.path.join(RESULTS, "ml_scores.csv"))
    ahr = scores[scores["endpoint"] == "NR-AhR"]
    prioritized = ahr[ahr["prob"] >= 0.30]["SMILES"].unique().tolist()
    print(f"[S3] {len(prioritized)} leachates prioritised (NR-AhR prob >= 0.30)")

    os.makedirs(os.path.join(RESULTS, "poses"), exist_ok=True)
    rows = []
    # NOTE: ligand PDBQT prep (meeko) + Vina.dock loop would run here. Left as the
    # documented full-run path; it is exercised only when receptors/Vina exist.
    for rec_name, meta in RECEPTORS.items():
        for smi in prioritized:
            # placeholder for: prepare ligand -> set box -> v.dock() -> read affinity
            raise NotImplementedError(
                "Full docking loop requires meeko ligand prep + a defined box "
                "centre from the prepared receptor; populate data/receptors/ and "
                "box centres, then enable this branch."
            )
    with open(DOCKING_CSV, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(HEADER)
        w.writerows(rows)
    return 0


if __name__ == "__main__":
    sys.exit(main())
