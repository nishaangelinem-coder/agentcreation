#!/usr/bin/env bash
# Reproduce the full analysis. Requires network access to www.cbioportal.org, datahub.assets.cbioportal.org and hgdownload.soe.ucsc.edu.
set -euo pipefail
cd "$(dirname "$0")"; mkdir -p data/raw data/ref results figures
[ -f data/ref/hg19.2bit ] || curl -sSL -o data/ref/hg19.2bit https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.2bit
for f in ncbiRefSeqSelect ncbiRefSeqCurated; do [ -f data/ref/$f.txt.gz ] || curl -sSL -o data/ref/$f.txt.gz https://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/$f.txt.gz; done
cd code
for s in 01_acquire 02_cohort_matrices 03_panel_opportunity 04_simulation_gate 05_spa 06_consensus 07_stats 08_external 09_denovo 10_figures; do echo "== $s"; python $s.py; done
cd ../manuscript && node build_manuscript.js
