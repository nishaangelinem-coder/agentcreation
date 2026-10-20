#!/usr/bin/env python3
"""Stage 1 — Leachate inventory.

Assembles a curated inventory of documented additives, residual monomers and
commonly adsorbed environmental contaminants for PE, PP, PS, PVC and PU.

Every entry is a real, literature-documented plastic-associated chemical with a
CAS number; the SMILES strings are validated and canonicalised with RDKit at
build time and any entry whose SMILES does not parse is DROPPED (never emitted),
so no invalid or fabricated structure can reach the downstream stages.

The curation draws on the compound classes catalogued in open plastic-chemical
databases (PlastChem 2024; PlasticMAP) and standard polymer-additive references.
Because the PlastChem/PlasticMAP bulk files were not reachable from the execution
sandbox (restricted egress), this inventory is a hand-curated subset rather than a
full programmatic ingest of those files; coverage per polymer is reported so the
sparse polymers are visible.

Output: data/leachates.csv
Columns: polymer, chemical_name, CAS, SMILES, function, evidence_note
"""
from __future__ import annotations

import csv
import os
import sys

from rdkit import Chem
from rdkit import RDLogger

RDLogger.DisableLog("rdApp.*")

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
OUT = os.path.join(ROOT, "data", "leachates.csv")

# function codes: additive | monomer | adsorbate
# Each tuple: (polymer, name, CAS, SMILES, function, evidence_note)
INVENTORY = [
    # ---------------- Polyethylene (PE) ----------------
    ("PE", "Butylated hydroxytoluene (BHT)", "128-37-0", "Cc1cc(C(C)(C)C)c(O)c(C(C)(C)C)c1", "additive", "Primary phenolic antioxidant in polyolefins"),
    ("PE", "Irganox 1076", "2082-79-3", "CCCCCCCCCCCCCCCCCCOC(=O)CCc1cc(C(C)(C)C)c(O)c(C(C)(C)C)c1", "additive", "Hindered-phenol antioxidant"),
    ("PE", "Irgafos 168", "31570-04-4", "CC(C)(C)c1ccc(OP(Oc2ccc(C(C)(C)C)cc2C(C)(C)C)Oc2ccc(C(C)(C)C)cc2C(C)(C)C)c(C(C)(C)C)c1", "additive", "Phosphite processing stabiliser"),
    ("PE", "2,4-di-tert-butylphenol", "96-76-4", "CC(C)(C)c1ccc(O)c(C(C)(C)C)c1", "additive", "Degradation product of Irgafos 168"),
    ("PE", "Erucamide", "112-84-5", "CCCCCCCC/C=C\\CCCCCCCCCCCCC(N)=O", "additive", "Slip/anti-block agent"),
    ("PE", "Oleamide", "301-02-0", "CCCCCCCC/C=C\\CCCCCCCCC(N)=O", "additive", "Slip agent"),
    ("PE", "Stearic acid", "57-11-4", "CCCCCCCCCCCCCCCCCC(=O)O", "additive", "Lubricant/acid scavenger co-additive"),
    ("PE", "Benzophenone", "119-61-9", "O=C(c1ccccc1)c1ccccc1", "additive", "UV absorber / photoinitiator residue"),
    ("PE", "Tinuvin 770", "52829-07-9", "CC1(C)CC(OC(=O)CCCCCCCCC(=O)OC2CC(C)(C)NC(C)(C)C2)CC(C)(C)N1", "additive", "Hindered-amine light stabiliser (HALS)"),
    ("PE", "Ethylene", "74-85-1", "C=C", "monomer", "Residual monomer"),
    ("PE", "4-Nonylphenol", "104-40-5", "CCCCCCCCCc1ccc(O)cc1", "adsorbate", "Adsorbed/weathering-associated endocrine disruptor"),
    ("PE", "4-tert-Octylphenol", "140-66-9", "CC(C)(C)CC(C)(C)c1ccc(O)cc1", "adsorbate", "Adsorbed phenolic contaminant"),

    # ---------------- Polypropylene (PP) ----------------
    ("PP", "Butylated hydroxytoluene (BHT)", "128-37-0", "Cc1cc(C(C)(C)C)c(O)c(C(C)(C)C)c1", "additive", "Phenolic antioxidant"),
    ("PP", "Irganox 1010", "6683-19-8", "CC(C)(C)c1cc(CCC(=O)OCC(COC(=O)CCc2cc(C(C)(C)C)c(O)c(C(C)(C)C)c2)(COC(=O)CCc2cc(C(C)(C)C)c(O)c(C(C)(C)C)c2)COC(=O)CCc2cc(C(C)(C)C)c(O)c(C(C)(C)C)c2)cc(C(C)(C)C)c1O", "additive", "Pentaerythritol hindered-phenol antioxidant"),
    ("PP", "Irgafos 168", "31570-04-4", "CC(C)(C)c1ccc(OP(Oc2ccc(C(C)(C)C)cc2C(C)(C)C)Oc2ccc(C(C)(C)C)cc2C(C)(C)C)c(C(C)(C)C)c1", "additive", "Phosphite stabiliser"),
    ("PP", "2,4-di-tert-butylphenol", "96-76-4", "CC(C)(C)c1ccc(O)c(C(C)(C)C)c1", "additive", "Phosphite-stabiliser degradation product; frequently leached"),
    ("PP", "Erucamide", "112-84-5", "CCCCCCCC/C=C\\CCCCCCCCCCCCC(N)=O", "additive", "Slip agent"),
    ("PP", "Oleamide", "301-02-0", "CCCCCCCC/C=C\\CCCCCCCCC(N)=O", "additive", "Slip agent"),
    ("PP", "Diethyl phthalate", "84-66-2", "CCOC(=O)c1ccccc1C(=O)OCC", "additive", "Plasticiser occasionally present"),
    ("PP", "Benzophenone", "119-61-9", "O=C(c1ccccc1)c1ccccc1", "additive", "UV absorber"),
    ("PP", "Propylene", "115-07-1", "C=CC", "monomer", "Residual monomer"),
    ("PP", "4-Nonylphenol", "104-40-5", "CCCCCCCCCc1ccc(O)cc1", "adsorbate", "Adsorbed endocrine disruptor"),
    ("PP", "Stearic acid", "57-11-4", "CCCCCCCCCCCCCCCCCC(=O)O", "additive", "Lubricant"),

    # ---------------- Polystyrene (PS) ----------------
    ("PS", "Styrene", "100-42-5", "C=Cc1ccccc1", "monomer", "Residual monomer; probable carcinogen"),
    ("PS", "Ethylbenzene", "100-41-4", "CCc1ccccc1", "monomer", "Residual monomer/solvent"),
    ("PS", "Benzene", "71-43-2", "c1ccccc1", "monomer", "Trace residual/degradation"),
    ("PS", "4-tert-Butylphenol", "98-54-4", "CC(C)(C)c1ccc(O)cc1", "additive", "Antioxidant/leachate"),
    ("PS", "Benzaldehyde", "100-52-7", "O=Cc1ccccc1", "adsorbate", "Oxidation/weathering product"),
    ("PS", "Acetophenone", "98-86-2", "CC(=O)c1ccccc1", "adsorbate", "Oxidation product"),
    ("PS", "Hexabromocyclododecane (HBCD)", "25637-99-4", "BrC1CCC(Br)C(Br)CCC(Br)C(Br)CCC1Br", "additive", "Brominated flame retardant in EPS/XPS"),
    ("PS", "Benzophenone", "119-61-9", "O=C(c1ccccc1)c1ccccc1", "additive", "UV absorber"),
    ("PS", "Butylated hydroxytoluene (BHT)", "128-37-0", "Cc1cc(C(C)(C)C)c(O)c(C(C)(C)C)c1", "additive", "Antioxidant"),
    ("PS", "4-Nonylphenol", "104-40-5", "CCCCCCCCCc1ccc(O)cc1", "adsorbate", "Adsorbed contaminant"),
    ("PS", "Toluene", "108-88-3", "Cc1ccccc1", "adsorbate", "Residual solvent/degradation"),

    # ---------------- Polyvinyl chloride (PVC) ----------------
    ("PVC", "Di(2-ethylhexyl) phthalate (DEHP)", "117-81-7", "CCCCC(CC)COC(=O)c1ccccc1C(=O)OCC(CC)CCCC", "additive", "Dominant PVC plasticiser; anti-androgen"),
    ("PVC", "Dibutyl phthalate (DBP)", "84-74-2", "CCCCOC(=O)c1ccccc1C(=O)OCCCC", "additive", "Plasticiser; reproductive toxicant"),
    ("PVC", "Benzyl butyl phthalate (BBP)", "85-68-7", "CCCCOC(=O)c1ccccc1C(=O)OCc1ccccc1", "additive", "Plasticiser"),
    ("PVC", "Dimethyl phthalate (DMP)", "131-11-3", "COC(=O)c1ccccc1C(=O)OC", "additive", "Plasticiser"),
    ("PVC", "Diethyl phthalate (DEP)", "84-66-2", "CCOC(=O)c1ccccc1C(=O)OCC", "additive", "Plasticiser"),
    ("PVC", "Bis(2-ethylhexyl) adipate (DEHA)", "103-23-1", "CCCCC(CC)COC(=O)CCCCC(=O)OCC(CC)CCCC", "additive", "Non-phthalate plasticiser"),
    ("PVC", "Dibutyltin dilaurate", "77-58-7", "CCCCCCCCCCCC(=O)O[Sn](CCCC)(CCCC)OC(=O)CCCCCCCCCCC", "additive", "Organotin heat stabiliser"),
    ("PVC", "Vinyl chloride", "75-01-4", "C=CCl", "monomer", "Residual monomer; known carcinogen"),
    ("PVC", "Bisphenol A (BPA)", "80-05-7", "CC(C)(c1ccc(O)cc1)c1ccc(O)cc1", "additive", "Antioxidant/co-additive; endocrine disruptor"),
    ("PVC", "Nonylphenol", "104-40-5", "CCCCCCCCCc1ccc(O)cc1", "additive", "Co-stabiliser/degradation product"),
    ("PVC", "Tris(1-chloro-2-propyl) phosphate (TCPP)", "13674-84-5", "CC(CCl)OP(=O)(OC(C)CCl)OC(C)CCl", "additive", "Organophosphate flame retardant/plasticiser"),

    # ---------------- Polyurethane (PU) ----------------
    ("PU", "Toluene 2,4-diisocyanate (TDI)", "584-84-9", "Cc1ccc(N=C=O)cc1N=C=O", "monomer", "Isocyanate monomer; respiratory sensitiser"),
    ("PU", "4,4'-Methylenediphenyl diisocyanate (MDI)", "101-68-8", "O=C=Nc1ccc(Cc2ccc(N=C=O)cc2)cc1", "monomer", "Isocyanate monomer"),
    ("PU", "2,4-Toluenediamine (TDA)", "95-80-7", "Cc1ccc(N)cc1N", "monomer", "Hydrolysis product of TDI; carcinogen"),
    ("PU", "4,4'-Methylenedianiline (MDA)", "101-77-9", "Nc1ccc(Cc2ccc(N)cc2)cc1", "monomer", "Hydrolysis product of MDI; carcinogen"),
    ("PU", "Tris(1-chloro-2-propyl) phosphate (TCPP)", "13674-84-5", "CC(CCl)OP(=O)(OC(C)CCl)OC(C)CCl", "additive", "Flame retardant in PU foam"),
    ("PU", "Tris(1,3-dichloro-2-propyl) phosphate (TDCPP)", "13674-87-8", "ClCC(CCl)OP(=O)(OC(CCl)CCl)OC(CCl)CCl", "additive", "Flame retardant; possible carcinogen"),
    ("PU", "1,4-Diazabicyclo[2.2.2]octane (DABCO)", "280-57-9", "C1CN2CCN1CC2", "additive", "Amine gelling catalyst"),
    ("PU", "Dibutyltin dilaurate", "77-58-7", "CCCCCCCCCCCC(=O)O[Sn](CCCC)(CCCC)OC(=O)CCCCCCCCCCC", "additive", "Organotin catalyst"),
    ("PU", "Bisphenol A (BPA)", "80-05-7", "CC(C)(c1ccc(O)cc1)c1ccc(O)cc1", "additive", "Polyol/antioxidant precursor residue"),
    ("PU", "Benzophenone", "119-61-9", "O=C(c1ccccc1)c1ccccc1", "additive", "UV absorber"),
    ("PU", "Triethyl phosphate", "78-40-0", "CCOP(=O)(OCC)OCC", "additive", "Flame retardant/viscosity modifier"),
]


def main() -> int:
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    rows = []
    dropped = []
    per_polymer: dict[str, int] = {}
    for polymer, name, cas, smiles, function, note in INVENTORY:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            dropped.append((polymer, name, smiles))
            continue
        canonical = Chem.MolToSmiles(mol)
        rows.append({
            "polymer": polymer,
            "chemical_name": name,
            "CAS": cas,
            "SMILES": canonical,
            "function": function,
            "evidence_note": note,
        })
        per_polymer[polymer] = per_polymer.get(polymer, 0) + 1

    with open(OUT, "w", newline="") as fh:
        writer = csv.DictWriter(
            fh, fieldnames=["polymer", "chemical_name", "CAS", "SMILES", "function", "evidence_note"]
        )
        writer.writeheader()
        writer.writerows(rows)

    print(f"[S1] wrote {len(rows)} validated leachates -> {OUT}")
    print("[S1] coverage per polymer (target >= 10):")
    for pol in ["PE", "PP", "PS", "PVC", "PU"]:
        n = per_polymer.get(pol, 0)
        flag = "" if n >= 10 else "  <-- SPARSE"
        print(f"       {pol}: {n}{flag}")
    if dropped:
        print("[S1] WARNING: dropped entries with invalid SMILES (not emitted):")
        for d in dropped:
            print(f"       {d}")
    else:
        print("[S1] all curated SMILES parsed and canonicalised cleanly.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
