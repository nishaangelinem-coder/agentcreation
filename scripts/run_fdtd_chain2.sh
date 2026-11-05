#!/bin/sh
cd /home/user/agentcreation
[ -f results/array_nodecouple_post.npz ] || python3 scripts/05_array_fdtd.py nodecouple 7.0 0 > results/array_nodecouple.log 2>&1
for L in 4.0 3.6 4.4; do
  tag=decoupled; [ "$L" != "4.0" ] && tag=decoupled_L$L
  [ -f results/array_${tag}_post.npz ] || python3 scripts/05_array_fdtd.py $tag 7.0 $L > results/array_$tag.log 2>&1
done
echo done >> results/fdtd_chain.log
