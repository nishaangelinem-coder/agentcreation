#!/usr/bin/env python
"""Phase 3b: trinucleotide-opportunity profile of the MSK-IMPACT capture territory (per panel version) vs genome/exome.
Capture territory is approximated by the coding exons (+/-10 bp) of each panel's genes (RefSeq Select, hg19),
because the vendor capture BED is not distributed with the public study."""
import gzip, json, numpy as np, pandas as pd, py2bit
from common import *
tb = py2bit.open(f"{REF}/hg19.2bit")
LUT = np.full(256, 4, np.int64)
for i, b in enumerate("ACGT"): LUT[ord(b)] = i; LUT[ord(b.lower())] = i
TRIS = [a + b + c for a in "ACGT" for b in "ACGT" for c in "ACGT"]
PYR = sorted({t if t[1] in "CT" else revcomp(t) for t in TRIS})
def count64(seq):
    x = LUT[np.frombuffer(seq.encode(), np.uint8)]
    if len(x) < 3: return np.zeros(64, np.int64)
    a, b, c = x[:-2], x[1:-1], x[2:]; ok = (a < 4) & (b < 4) & (c < 4)
    return np.bincount((a * 16 + b * 4 + c)[ok], minlength=64)
def collapse(c64):
    s = pd.Series(0.0, index=PYR)
    for t, n in zip(TRIS, c64): s[t if t[1] in "CT" else revcomp(t)] += n
    return s
main = [f"chr{i}" for i in list(range(1, 23)) + ["X", "Y"]]
g = np.zeros(64, np.int64)
for c in main: g += count64(tb.sequence(c))
genome = collapse(g)
cols = ["bin", "name", "chrom", "strand", "txStart", "txEnd", "cdsStart", "cdsEnd", "exonCount", "exonStarts", "exonEnds", "score", "name2", "s1", "s2", "fr"]
def cds_table(fn):
    t = pd.read_csv(f"{REF}/{fn}", sep="\t", names=cols); t = t[t.chrom.isin(main) & (t.cdsEnd > t.cdsStart)]
    return t
sel, cur = cds_table("ncbiRefSeqSelect.txt.gz"), cds_table("ncbiRefSeqCurated.txt.gz")
def intervals(row, flank=10):
    out = []
    for s, e in zip(map(int, row.exonStarts.strip(",").split(",")), map(int, row.exonEnds.strip(",").split(","))):
        s, e = max(s, int(row.cdsStart)), min(e, int(row.cdsEnd))
        if e > s: out.append((str(row.chrom), int(s - flank), int(e + flank)))
    return out
def merge(iv):
    iv = sorted(iv); out = []
    for c, s, e in iv:
        if out and out[-1][0] == c and s <= out[-1][2]: out[-1][2] = max(out[-1][2], e)
        else: out.append([c, s, e])
    return out
ALIAS = {"ABRAXAS1": "FAM175A", "MRE11": "MRE11A", "NSD2": "WHSC1", "NSD3": "WHSC1L1", "PRKN": "PARK2", "KMT5A": "SETD8", "ELOC": "TCEB1", "CARS1": "CARS",
         "TENT5C": "FAM46C", "EMSY": "C11orf30", "MYCL": "MYCL1", "H3C1": "HIST1H3A", "KNSTRN": "C15orf23", "GSDME": "DFNA5", "H1-2": "HIST1H1C", "MRTFA": "MKL1"}
exome = np.zeros(64, np.int64)
for _, r in sel.iterrows():
    pass
ex_iv = merge([x for _, r in sel.iterrows() for x in intervals(r, 0)])
for c, s, e in ex_iv: exome += count64(tb.sequence(c, s - 1, e + 1))
exome = collapse(exome)
opp = {"genome": genome, "exome_refseq_cds": exome}; info = {"exome_cds_bp": int(sum(e - s for _, s, e in ex_iv))}
terr = {}
for p in ["IMPACT341", "IMPACT410", "IMPACT468", "IMPACT505"]:
    genes = [l.strip() for l in open(f"{REF}/{p}.genes.txt") if l.strip()]; iv = []; miss = []
    for gname in genes:
        hit = None
        for tab in (sel, cur):
            for nm in (gname, ALIAS.get(gname)):
                h = tab[tab.name2 == nm] if nm else tab.iloc[0:0]
                if len(h): hit = h.assign(l=h.cdsEnd - h.cdsStart).sort_values("l").iloc[-1]; break
            if hit is not None: break
        if hit is None: miss.append(gname)
        else: iv += intervals(hit)
    iv = merge(iv); terr[p] = iv
    c = np.zeros(64, np.int64)
    for ch, s, e in iv: c += count64(tb.sequence(ch, s - 1, e + 1))
    opp[p] = collapse(c); info[p] = {"genes": len(genes), "unmatched_genes": miss, "territory_bp": int(sum(e - s for _, s, e in iv))}
    pd.DataFrame(iv).to_csv(f"{REF}/{p}.cds_flank10.hg19.bed", sep="\t", header=False, index=False)
O = pd.DataFrame(opp); O.to_csv(f"{RES}/trinucleotide_opportunity_counts.tsv", sep="\t")
P = O / O.sum(); P.to_csv(f"{RES}/trinucleotide_opportunity_proportions.tsv", sep="\t")
info["cosine_panel_vs_genome"] = {p: cosine(P[p], P.genome) for p in P.columns[1:]}
info["cosine_505_vs_other_panels"] = {p: cosine(P[p], P.IMPACT505) for p in ["IMPACT341", "IMPACT410", "IMPACT468", "exome_refseq_cds"]}
info["max_fold_enrichment_505_vs_genome"] = (P.IMPACT505 / P.genome).sort_values().iloc[[0, -1]].round(3).to_dict()
# empirical check: share of observed cohort SNVs that fall inside the modelled territory (IMPACT505)
snv = pd.read_parquet(f"{RES}/snv_annotated.parquet", columns=["Chromosome", "Start_Position", "Variant_Classification"])
import bisect
idx = {}
for c, s, e in terr["IMPACT505"]: idx.setdefault(c, []).append((s, e))
def inside(c, p):
    L = idx.get("chr" + str(c), []); i = bisect.bisect_right(L, (p, 10**12)) - 1
    return i >= 0 and L[i][0] < p <= L[i][1]
ins = np.array([inside(c, p) for c, p in zip(snv.Chromosome, snv.Start_Position)])
info["frac_all_SNVs_inside_modelled_territory"] = float(ins.mean())
info["outside_by_class"] = snv.Variant_Classification[~ins].value_counts().head(6).to_dict()
json.dump(info, open(f"{RES}/panel_opportunity_info.json", "w"), indent=2); print(json.dumps(info, indent=2))
