#!/bin/sh
cd /home/user/agentcreation
[ -f results/cma_ref.npz ] || python3 scripts/03b_cma_dual_plate.py ref 7.0 0 > results/cma_ref.log 2>&1
[ -f results/cma_slits.npz ] || python3 scripts/03b_cma_dual_plate.py slits 7.0 4.0 > results/cma_slits.log 2>&1
echo done > results/cma_chain.log
