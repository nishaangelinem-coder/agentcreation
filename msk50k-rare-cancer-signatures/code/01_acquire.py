#!/usr/bin/env python
"""Phase 1 - data acquisition with provenance (cBioPortal REST API + public datahub)."""
import hashlib, json, os, subprocess, sys, datetime, requests, tarfile
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RAW = os.path.join(ROOT, "data", "raw"); REF = os.path.join(ROOT, "data", "ref")
API = "https://www.cbioportal.org/api"; HUB = "https://datahub.assets.cbioportal.org"
STUDY = "msk_impact_50k_2026"
EXTERNAL = ["acc_tcga_pan_can_atlas_2018", "thym_tcga_pan_can_atlas_2018", "angs_painter_2025", "angs_project_painter_2018"]

def sha256(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""): h.update(b)
    return h.hexdigest()

def fetch(url, dest):
    if not os.path.exists(dest):
        subprocess.run(["curl", "-sS", "-L", "--fail", "--retry", "4", "-o", dest, url], check=True)
    return dest

prov = {"access_date_utc": datetime.datetime.utcnow().isoformat() + "Z", "studies": {}}
for sid in [STUDY] + EXTERNAL:
    tgz = os.path.join(RAW, sid + ".tar.gz")
    try: fetch(f"{HUB}/{sid}.tar.gz", tgz)
    except subprocess.CalledProcessError:
        print("unavailable:", sid); prov["studies"][sid] = "unavailable"; continue
    out = os.path.join(RAW, "study")
    if not os.path.isdir(os.path.join(out, sid)):
        with tarfile.open(tgz) as t: t.extractall(out, filter="data")
    meta = requests.get(f"{API}/studies/{sid}", timeout=60).json()
    prov["studies"][sid] = {"tarball_sha256": sha256(tgz), "name": meta.get("name"), "pmid": meta.get("pmid"),
                            "citation": meta.get("citation"), "api_allSampleCount": meta.get("allSampleCount")}
# main-study integrity checks
S = os.path.join(RAW, "study", STUDY)
prov["maf_sha256"] = sha256(os.path.join(S, "data_mutations.txt"))
ids = [l.split("\t")[0] for l in open(os.path.join(S, "data_clinical_sample.txt")) if not l.startswith("#")][1:]
prov["n_samples_file"] = len(ids); prov["n_missing_sampleId"] = sum(1 for i in ids if not i.strip())
api_n = requests.get(f"{API}/studies/{STUDY}/samples", params={"projection": "META"}, timeout=60).headers.get("total-count")
prov["n_samples_api"] = int(api_n) if api_n else None
prov["license"] = "CC BY-NC-ND 4.0 (per study meta)"
# datahub commit hash of the study meta file (GitHub API may be unreachable -> recorded as None)
try:
    r = requests.get("https://api.github.com/repos/cBioPortal/datahub/commits", params={"path": f"public/{STUDY}", "per_page": 1}, timeout=30)
    prov["datahub_commit"] = r.json()[0]["sha"]
except Exception as e: prov["datahub_commit"] = None
# gene panels via API
for p in ["IMPACT341", "IMPACT410", "IMPACT468", "IMPACT505"]:
    g = requests.get(f"{API}/gene-panels/{p}", timeout=60).json()
    with open(os.path.join(REF, f"{p}.genes.txt"), "w") as f: f.write("\n".join(sorted(x["hugoGeneSymbol"] for x in g["genes"])) + "\n")
    prov.setdefault("panels", {})[p] = len(g["genes"])
json.dump(prov, open(os.path.join(ROOT, "results", "provenance.json"), "w"), indent=2)
print(json.dumps(prov, indent=2))
assert prov["n_samples_file"] == 54331 and prov["n_missing_sampleId"] == 0
