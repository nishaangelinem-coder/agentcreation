# Supplementary Material S1: in-silico case study

`S1_closed_loop_cap_simulation.py` reproduces Section 7, Table 4 and Fig. 2 of the manuscript.

```
pip install numpy
python S1_closed_loop_cap_simulation.py        # writes S1_simulation_results.json
```

- Deterministic: master seed 2026. Plants are drawn once and shared by all three arms. Each closed-loop arm has its own child random stream (`numpy.random.SeedSequence.spawn`).
- Model: τ dT/dt = −(T − T_amb) + K·P, forward Euler with dt = 0.1 s, CEM43 per Sapareto & Dewey (1984).
- Dose error is computed from the **true** temperature. The controller stops the treatment using CEM43 from the **measured** temperature.
- The parameters are illustrative and not fitted to any device. Do not cite these results as device performance.

Tested with Python 3.13 and NumPy 2.5.
