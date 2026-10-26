#!/usr/bin/env bash
# Copy the latest figures/tables from the analysis into the manuscript bundle and compile.
set -euo pipefail
cd "$(dirname "$0")"
mkdir -p figures tables
cp ../figures/fig*.pdf figures/
cp ../results/table_groups.tex ../results/table_anova.tex tables/
latexmk -pdf -interaction=nonstopmode -halt-on-error main.tex
