#!/bin/bash
# Full pipeline: data generation -> benchmarks -> sensitivity -> tables & figures.
set -e
cd "$(dirname "$0")"
source .venv/bin/activate
[ -f data/ieee14.pkl ] || python -m fdia_dt.gen_all ieee14
[ -f data/ieee118.pkl ] || python -m fdia_dt.gen_all ieee118
for c in ieee14 ieee118; do python -m fdia_dt.run_experiments $c 2>&1 | grep -v -i numba | tee results/run_$c.log; done
for c in ieee14 ieee118; do python -m fdia_dt.sensitivity $c 2>&1 | grep -v -i numba | tee results/sens_$c.log; done
python -m fdia_dt.make_tables > /dev/null
python -m fdia_dt.figures
echo ALL_DONE
