#!/bin/sh
cd /home/user/agentcreation
while pgrep -f run_D94.sh > /dev/null; do sleep 20; done
python3 scripts/05_array_fdtd.py D94_long_o40 D=9.4 cav=2 Ls=4.0 xoff=0.8 > results/array_D94_long_o40.log 2>&1 &
python3 scripts/05_array_fdtd.py D94_long_base D=9.4 cav=2 > results/array_D94_long_base.log 2>&1
wait
echo doneD94b >> results/fdtd_chain.log
