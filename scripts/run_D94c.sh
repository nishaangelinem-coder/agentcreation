#!/bin/sh
cd /home/user/agentcreation
while pgrep -f "run_D94b.sh|run_D94.sh" > /dev/null; do sleep 20; done
python3 scripts/05_array_fdtd.py D94_o44 D=9.4 cav=1 Ls=4.4 xoff=0.8 > results/array_D94_o44.log 2>&1 &
python3 scripts/05_array_fdtd.py D94_o36 D=9.4 cav=1 Ls=3.6 xoff=0.8 > results/array_D94_o36.log 2>&1
wait
python3 scripts/05_array_fdtd.py D94_o40x12 D=9.4 cav=1 Ls=4.0 xoff=1.2 > results/array_D94_o40x12.log 2>&1
echo doneD94c >> results/fdtd_chain.log
