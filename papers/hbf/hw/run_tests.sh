#!/bin/sh
# Runs every functional test; prints one RESULT line per DUT.
set -e
cd "$(dirname "$0")"
mkdir -p build
run() { iverilog -g2012 -o build/sim_$1 -DDUT=$1 -DN=$2 -DSIGNED=$3 -DLAT=$4 $5 tb/tb_mult.v rtl/$1.v && vvp -n build/sim_$1 | grep -E "RESULT|MISMATCH" | sed "s/^/$1 /"; }
run bw_array8    8 1 0
run bw_dadda8    8 1 0
run u_array4     4 0 0
run bw_array8_p3 8 1 3 -DPIPE
run bw_dadda8_p2 8 1 2 -DPIPE
run bw_wallace8  8 1 0
run booth8       8 1 0
run booth8_p2    8 1 2 -DPIPE
run behav_mult8  8 1 0
python3 gen_fp8_vectors.py >/dev/null
iverilog -g2012 -o build/sim_fp8 tb/tb_fp8_mac.v rtl/fp8_mac.v rtl/u_array4.v && vvp -n build/sim_fp8 | grep -E "RESULT|MISMATCH" | sed "s/^/fp8_mac /"
iverilog -g2012 -o build/sim_sys tb/tb_systolic.v rtl/systolic_array.v rtl/int8_pe.v rtl/bw_array8.v && vvp -n build/sim_sys | grep -E "RESULT|MISMATCH" | sed "s/^/systolic_array /"
for N in 4 16; do for d in bw_array bw_wallace bw_dadda booth behav_mult; do run $d$N $N 1 0; done; done
run bw_array16_p4 16 1 4 -DPIPE
for N in 8 16; do run bw_dadda${N}_ks $N 1 0; run booth${N}_ks $N 1 0; done
