# RT Patient-Monitor Simulator

A real-time patient-monitor simulator with **measured alarm-latency deadlines** and a
**hazard log styled after ISO 14971 / IEC 62304**. Synthetic ECG (500 Hz) and SpO2/PPG
(100 Hz) are acquired, processed (Pan-Tompkins-style QRS detection, ratio-of-ratios SpO2)
and turned into prioritised, latching alarms by a fixed-priority thread architecture.

The code is plain POSIX C11 (pthreads, `SCHED_FIFO`, priority-inheritance mutexes,
`clock_nanosleep`, `mlockall`) so that one source tree targets **QNX Neutrino** and **Linux**.

> **Status - read this first**
> * This is a research/teaching simulator, **not a medical device**. The hazard log is
>   *styled after* the standards; it is not evidence of conformity.
> * All results in `results/` and in the paper were measured on **stock Linux in a virtual
>   machine**. The QNX build target (`make qnx`) has **not been compiled or run** yet.
> * Signals are synthetic. Ventricular fibrillation is not classified.

## Quick start (Linux)

```bash
make                 # builds build/monitor, build/loadgen, build/test_detect
make verify          # off-target verification suite -> results/verification.csv
```

Live demo with a terminal display (run as root for real-time scheduling):

```bash
sudo ./build/monitor --scenario clinical --ui
```

The 210 s clinical scenario plays: normal rhythm, asystole, tachycardia (180 bpm),
bradycardia (32 bpm), desaturation (80 %), leads-off. The banner turns red / yellow / cyan
for high / medium / low priority alarms.

## Reproducing the paper

```bash
sudo ./scripts/run_experiments.sh       # ~1 h: 24 probe runs, 10 fault injections, 4 clinical runs
python3 scripts/analyze.py              # statistics, figures, LaTeX tables, run-time verdicts
python3 scripts/gen_hazard_log.py       # docs/HAZARD_LOG.md + traceability from executed tests
cd paper && latexmk -pdf main.tex
```

Python needs `numpy scipy pandas matplotlib`.

## Running on QNX (not yet done - contributions welcome)

1. Install QNX SDP 8.0 (a free non-commercial licence is available from QNX) and source
   `qnxsdp-env.sh`.
2. `make qnx QNX_ARCH=aarch64le` (Raspberry Pi 4 image) or `QNX_ARCH=x86_64` (QEMU/VM image).
3. Copy `build/` and `scripts/run_experiments.sh` to the target and run the campaign as root.
4. Copy `results/campaign/` back and run the two Python scripts on a host.

Things to check on the first port: thread priorities above 63 need root or the
`PROCMGR_AID_PRIORITY` ability; the CPU-binding call in `src/rt_util.c` (`ThreadCtl` run
mask) is untested; `loadgen` writes to `/tmp`.

## Layout

| Path | Content |
|---|---|
| `src/siggen.c` | ECG + PPG source with scripted scenarios and ground truth |
| `src/ecg_proc.c` | QRS detector, heart rate, asystole / tachy / brady / leads-off conditions |
| `src/spo2_proc.c` | SpO2 estimator, desaturation and low-perfusion conditions |
| `src/alarm.c` | Alarm state machine: priorities, latching, acknowledge (no OS calls) |
| `src/monitor.c` | Real-time application: 5 threads, logging, fault injection, terminal UI |
| `src/chan.c`, `src/rt_util.c` | Bounded queue with PI mutex; scheduling / time helpers |
| `src/loadgen.c` | CPU, memory and I/O interference load |
| `tests/test_detect.c` | Off-target verification suite (27 checks) |
| `docs/hazard_log.csv`, `docs/requirements.csv` | Source tables of the hazard analysis |
| `docs/HAZARD_LOG.md`, `docs/traceability.csv` | Generated report (do not edit by hand) |
| `results/` | Raw campaign data, summaries, verification verdicts |
| `paper/` | Manuscript (LaTeX), generated tables and figures |

## `monitor` options

```
--policy fifo|other     fixed-priority real-time (default) or time-sharing baseline
--alarm  event|poll     event-driven alarm path (default) or polled (--poll-ms N, default 5)
--scenario clinical|probe|normal
--duration S  --seed N  --cpu N|-1 (core binding)  --mlock 0|1
--inject-stall T        runaway loop in the processing thread at time T (watchdog test)
--ui                    terminal display
--out PREFIX            writes PREFIX_{alarms,samples,truth}.csv and PREFIX_summary.json
```

Alarm latency is `t_fire - t_release(k*)`: from the nominal release instant of the sample
that makes the condition true to the instant the annunciator has been driven.
