#!/bin/sh
cd /home/user/agentcreation
[ -f results/array_D94_o44_post.npz ] || python3 scripts/05_array_fdtd.py D94_o44 D=9.4 cav=1 Ls=4.4 xoff=0.8 > results/array_D94_o44.log 2>&1 &
[ -f results/array_D94_o36_post.npz ] || python3 scripts/05_array_fdtd.py D94_o36 D=9.4 cav=1 Ls=3.6 xoff=0.8 > results/array_D94_o36.log 2>&1
wait
[ -f results/array_D94_o40x12_post.npz ] || python3 scripts/05_array_fdtd.py D94_o40x12 D=9.4 cav=1 Ls=4.0 xoff=1.2 > results/array_D94_o40x12.log 2>&1 &
[ -f results/array_D94_long_o40_post.npz ] || python3 scripts/05_array_fdtd.py D94_long_o40 D=9.4 cav=2 Ls=4.0 xoff=0.8 > results/array_D94_long_o40.log 2>&1
wait
echo doneopt >> results/fdtd_chain.log
