#!/bin/sh
# FDTD chain: undecoupled array, then decoupled arrays with slit lengths around the CMA-predicted value
cd /home/user/agentcreation
python3 scripts/05_array_fdtd.py nodecouple 7.0 0 > results/array_nodecouple.log 2>&1
while ! grep -q slit_L_28GHz results/cma_calibration.json 2>/dev/null; do sleep 20; done
LD=$(python3 -c "import json;print(round(json.load(open('results/cma_calibration.json'))['slit_L_28GHz'],1))")
echo "L_d = $LD" > results/fdtd_chain.log
python3 scripts/05_array_fdtd.py decoupled 7.0 $LD > results/array_decoupled.log 2>&1
for dl in -0.4 0.4; do
  L=$(python3 -c "print(round($LD+($dl),1))")
  python3 scripts/05_array_fdtd.py decoupled_L$L 7.0 $L > results/array_decoupled_L$L.log 2>&1
done
echo done >> results/fdtd_chain.log
