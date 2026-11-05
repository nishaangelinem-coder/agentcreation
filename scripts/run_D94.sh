#!/bin/sh
cd /home/user/agentcreation
cp results/array_reference.npz results/array_reference_D70.npz
python3 scripts/05_array_fdtd.py reference D=9.4 > results/array_reference_D94.log 2>&1
cp results/array_reference.npz results/array_reference_D94.npz
python3 scripts/05_array_fdtd.py D94_base D=9.4 cav=1 > results/array_D94_base.log 2>&1 &
python3 scripts/05_array_fdtd.py D94_c40 D=9.4 cav=1 Ls=4.0 > results/array_D94_c40.log 2>&1
wait
python3 scripts/05_array_fdtd.py D94_o40 D=9.4 cav=1 Ls=4.0 xoff=0.8 > results/array_D94_o40.log 2>&1 &
python3 scripts/05_array_fdtd.py D94_c36 D=9.4 cav=1 Ls=3.6 > results/array_D94_c36.log 2>&1
wait
echo doneD94 >> results/fdtd_chain.log
