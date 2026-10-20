# Near-flash multiplier study (Section 7 of the HBF paper)

Open-source flow: Python generator → structural Verilog → Icarus Verilog verification →
Yosys/ABC synthesis (SkyWater 130 nm, ASAP7) → unit-delay gate-level power → system model.

| Step | Command |
|---|---|
| Fetch and prepare cell libraries (~100 MB, not committed) | `python3 prepare_libs.py` |
| Regenerate a multiplier | `python3 gen_mult.py NAME N {array\|wallace\|dadda\|booth}[+ks] {s\|u} [cuts]` |
| Functional regression (26 suites) | `./run_tests.sh` |
| Synthesis + power sweep (27 designs x 2 libs x 2 modes x 2 data profiles) | `python3 syn_flow.py [design ...]` |
| Merge parallel-worker outputs | `python3 merge_results.py` |
| Figures 7–16 and table data | `python3 make_hw_figures.py` |

- `rtl/` — generated multipliers (array, Wallace, Dadda, Booth-4, ±Kogge–Stone, pipelined; 4/8/16 bit),
  `fp8_mac.v` (E4M3, exact 48-bit accumulation), `int8_pe.v`, `systolic_array.v`
- `tb/` — exhaustive / random / cycle-accurate testbenches; `gen_fp8_vectors.py` builds exact golden vectors
- `nmc_model.py` — S1 (HBM, no HBF) / S2 (HBF, far compute) / S3 (HBF, near-flash MAC) energy, capacity, throughput
- `results/` — `functional_tests.txt`, `syn_results.json` (108 synthesis records), `structure_stats.json`, `hw_tables.json`

Requires `iverilog`, `yosys` (tested with 0.33), Python 3 with matplotlib/numpy.
