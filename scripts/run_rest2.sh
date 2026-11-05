#!/bin/sh
cd /home/user/agentcreation
while pgrep -f run_rest.sh > /dev/null; do sleep 20; done
[ -f results/array_decoupled_L4.8_post.npz ] || python3 scripts/05_array_fdtd.py decoupled_L4.8 7.0 4.8 > results/array_decoupled_L4.8.log 2>&1
[ -f results/array_decoupled_cont_post.npz ] || python3 scripts/05_array_fdtd.py decoupled_cont 7.0 24.0 1 > results/array_decoupled_cont.log 2>&1
echo done2 >> results/fdtd_chain.log
