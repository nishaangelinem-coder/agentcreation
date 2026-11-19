#!/usr/bin/env bash
# Full measurement campaign. Run as root (SCHED_FIFO + mlockall) from the project root.
#   REPS (default 3), DUR probe-run seconds (default 100), CPU core for the monitor threads (default 1)
set -u
REPS=${REPS:-3}; DUR=${DUR:-100}; CPU=${CPU:-1}; OUT=${OUT:-results/campaign}
LOAD_ARGS=${LOAD_ARGS:-"8 2 1"}          # cpu-bound, memory-churn, file-I/O threads
mkdir -p "$OUT"
{ date -u; uname -a; nproc; grep -m1 "model name" /proc/cpuinfo 2>/dev/null; cat /sys/kernel/realtime 2>/dev/null || echo "realtime: no";
  cat /proc/sys/kernel/sched_rt_runtime_us 2>/dev/null; ${CC:-gcc} --version | head -1; echo "load: $LOAD_ARGS"; } > "$OUT/platform.txt" 2>&1

run() { # name load extra-args...
  local name=$1 load=$2; shift 2; local lp=""
  if [ "$load" = stress ]; then ./build/loadgen $LOAD_ARGS 100000 & lp=$!; sleep 2; fi
  ./build/monitor --cpu "$CPU" --out "$OUT/$name" "$@" | tee -a "$OUT/log.txt"
  if [ -n "$lp" ]; then kill "$lp" 2>/dev/null; wait "$lp" 2>/dev/null; sleep 1; fi
}

# A. latency probe: 2 (scheduling) x 2 (alarm path) x 2 (load), REPS repetitions, blocked by repetition
for rep in $(seq 1 "$REPS"); do for load in idle stress; do for pol in fifo other; do for al in event poll; do
  run "probe_${pol}_${al}_${load}_r${rep}" $load --scenario probe --duration "$DUR" --policy $pol --alarm $al --seed "$rep"
done; done; done; done

# B. watchdog fault injection (runaway processing thread), varied phase
i=0; for t in 8.003 8.007 8.011 8.017 8.019; do i=$((i+1)); for pol in fifo other; do
  run "stall_${pol}_event_stress_r${i}" stress --scenario normal --duration 12 --policy $pol --alarm event --inject-stall $t --seed "$i"
done; done

# C. clinical scenario in real time
run clinical_fifo_event_idle_r1   idle   --scenario clinical --policy fifo  --alarm event --seed 11
run clinical_fifo_event_stress_r1 stress --scenario clinical --policy fifo  --alarm event --seed 12
run clinical_fifo_event_stress_r2 stress --scenario clinical --policy fifo  --alarm event --seed 13
run clinical_other_event_stress_r1 stress --scenario clinical --policy other --alarm event --seed 12
date -u >> "$OUT/platform.txt"; echo CAMPAIGN_DONE | tee -a "$OUT/log.txt"
