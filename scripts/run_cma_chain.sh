#!/bin/sh
cd /home/user/agentcreation
while ! grep -q slit_L_28GHz results/cma_calibration.json 2>/dev/null; do sleep 20; done
LD=$(python3 -c "import json;print(round(json.load(open('results/cma_calibration.json'))['slit_L_28GHz'],1))")
python3 scripts/03b_cma_dual_plate.py ref 7.0 0 > results/cma_ref.log 2>&1
python3 scripts/03b_cma_dual_plate.py slits 7.0 $LD > results/cma_slits.log 2>&1
echo done > results/cma_chain.log
