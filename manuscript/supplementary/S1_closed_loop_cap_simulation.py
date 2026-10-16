"""Supplementary Material S1: illustrative Monte Carlo case study (Section 7).

Closed-loop versus open-loop delivery of a CEM43 thermal dose under a cold
atmospheric plasma (CAP) jet, using a generic first-order lumped thermal model.

    tau * dT/dt = -(T - T_amb) + K * P(t),   0 <= P <= P_max

The parameters are illustrative and NOT fitted to any device. Results show
mechanisms only and must not be cited as device performance.

Usage:  python S1_closed_loop_cap_simulation.py [--out S1_simulation_results.json]
Deterministic: the master seed is 2026 and each arm uses its own child stream.
"""

import argparse
import json

import numpy as np

SEED = 2026
N_TREATMENTS = 1000
DT = 0.1  # s

# Nominal plant
T_AMB_NOM = 32.0  # degC
TAU_NOM = 8.0  # s
K_NOM = 3.0  # degC / W
P_MAX = 5.0  # W

# Treatment specification
CEM_TARGET = 1.0  # min
T_SAFETY = 45.0  # degC, hard limit (true temperature)
T_SET = 44.0  # degC, closed-loop set-point
T_INTERLOCK = 44.7  # degC, supervisory interlock on the measured temperature
T_MAX_TIME = 300.0  # s, hard stop for any arm

# Open-loop calibration
P_OPEN = 4.0  # W -> 44 degC steady state on the nominal plant

# PI controller
KP = 1.2  # W / degC
KI = 0.15  # W / (degC s)

SENSORS = {
    "closed_nominal": {"lag": 1.0, "noise": 0.2, "bias_sd": 0.3},
    "closed_degraded": {"lag": 3.0, "noise": 0.5, "bias_sd": 1.0},
}

GAIN_BINS = [(0.75, 0.85), (0.85, 0.95), (0.95, 1.05), (1.05, 1.15), (1.15, 1.25)]


def cem43_rate(temp):
    """CEM43 accumulation rate in min per s (Sapareto and Dewey, 1984)."""
    r = np.where(temp >= 43.0, 0.5, 0.25)
    return r ** (43.0 - temp) / 60.0


def plant_step(temp, power, k, tau, t_amb):
    return temp + DT * (-(temp - t_amb) + k * power) / tau


def open_loop_duration():
    """Time for the nominal plant at P_OPEN to accumulate CEM_TARGET."""
    temp, dose, t = T_AMB_NOM, 0.0, 0.0
    while dose < CEM_TARGET:
        dose += float(cem43_rate(temp)) * DT
        temp = plant_step(temp, P_OPEN, K_NOM, TAU_NOM, T_AMB_NOM)
        t += DT
    return round(t, 1)


def draw_plants(rng, n):
    k_mult = rng.uniform(0.75, 1.25, n)
    tau = TAU_NOM * rng.uniform(0.8, 1.2, n)
    t_amb = rng.normal(T_AMB_NOM, 1.0, n)
    return k_mult, tau, t_amb


def run_open_loop(k_mult, tau, t_amb, t_open):
    n = len(k_mult)
    temp = t_amb.copy()
    dose = np.zeros(n)
    t_peak = temp.copy()
    steps = int(round(t_open / DT))
    for _ in range(steps):
        dose += cem43_rate(temp) * DT
        temp = plant_step(temp, P_OPEN, K_NOM * k_mult, tau, t_amb)
        t_peak = np.maximum(t_peak, temp)
    return dose, t_peak, np.full(n, steps * DT)


def run_closed_loop(k_mult, tau, t_amb, sensor, rng):
    n = len(k_mult)
    bias = rng.normal(0.0, sensor["bias_sd"], n)
    temp = t_amb.copy()
    t_lag = temp.copy()  # lagged (filtered) true temperature
    integ = np.zeros(n)
    dose_true = np.zeros(n)
    dose_meas = np.zeros(n)
    t_peak = temp.copy()
    t_end = np.full(n, T_MAX_TIME)
    active = np.ones(n, dtype=bool)
    alpha = DT / (sensor["lag"] + DT)
    k = K_NOM * k_mult
    for step in range(int(round(T_MAX_TIME / DT))):
        noise = rng.normal(0.0, sensor["noise"], n)
        meas = t_lag + bias + noise
        err = T_SET - meas
        u = KP * err + KI * integ
        power = np.clip(u, 0.0, P_MAX)
        # Anti-windup: integrate only when the actuator is not saturated
        # in the direction of the error.
        unsat = (u > 0.0) & (u < P_MAX)
        windup_ok = unsat | ((u >= P_MAX) & (err < 0)) | ((u <= 0.0) & (err > 0))
        integ = np.where(active & windup_ok, integ + err * DT, integ)
        power = np.where(meas >= T_INTERLOCK, 0.0, power)
        power = np.where(active, power, 0.0)

        dose_true += np.where(active, cem43_rate(temp) * DT, 0.0)
        dose_meas += np.where(active, cem43_rate(meas) * DT, 0.0)

        temp = plant_step(temp, power, k, tau, t_amb)
        t_lag = t_lag + alpha * (temp - t_lag)
        t_peak = np.where(active, np.maximum(t_peak, temp), t_peak)

        done = active & (dose_meas >= CEM_TARGET)
        t_end = np.where(done, (step + 1) * DT, t_end)
        active &= ~done
        if not active.any():
            break
    return dose_true, t_peak, t_end


def summarise(dose, t_peak, t_end, k_mult):
    err = (dose - CEM_TARGET) / CEM_TARGET * 100.0
    out = {
        "median_abs_dose_error_pct": float(np.median(np.abs(err))),
        "within_10pct_pct": float(np.mean(np.abs(err) <= 10.0) * 100.0),
        "exceeding_45C_pct": float(np.mean(t_peak > T_SAFETY) * 100.0),
        "median_time_s": float(np.median(t_end)),
        "dose_error_p5_p95_pct": [float(np.percentile(err, 5)), float(np.percentile(err, 95))],
        "by_gain_bin": [],
    }
    for lo, hi in GAIN_BINS:
        sel = (k_mult >= lo) & (k_mult < hi if hi < 1.25 else k_mult <= hi)
        e = err[sel]
        out["by_gain_bin"].append({
            "bin": f"{lo:.2f}-{hi:.2f}",
            "n": int(sel.sum()),
            "median_pct": float(np.median(e)),
            "p10_pct": float(np.percentile(e, 10)),
            "p90_pct": float(np.percentile(e, 90)),
        })
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="S1_simulation_results.json")
    args = ap.parse_args()

    master = np.random.SeedSequence(SEED)
    plant_ss, nominal_ss, degraded_ss = master.spawn(3)
    k_mult, tau, t_amb = draw_plants(np.random.default_rng(plant_ss), N_TREATMENTS)

    t_open = open_loop_duration()
    results = {
        "seed": SEED,
        "n_per_arm": N_TREATMENTS,
        "open_loop_power_W": P_OPEN,
        "open_loop_duration_s": t_open,
        "note": "Illustrative generic lumped model; parameters not fitted to any device.",
        "arms": {},
    }
    results["arms"]["open_loop"] = summarise(*run_open_loop(k_mult, tau, t_amb, t_open), k_mult)
    for name, ss in (("closed_nominal", nominal_ss), ("closed_degraded", degraded_ss)):
        res = run_closed_loop(k_mult, tau, t_amb, SENSORS[name], np.random.default_rng(ss))
        results["arms"][name] = summarise(*res, k_mult)

    with open(args.out, "w") as fh:
        json.dump(results, fh, indent=2)
    for name, arm in results["arms"].items():
        print(f"{name:16s} med|err| {arm['median_abs_dose_error_pct']:6.1f}%  "
              f"within10 {arm['within_10pct_pct']:5.1f}%  >45C {arm['exceeding_45C_pct']:5.1f}%  "
              f"t {arm['median_time_s']:5.1f}s")


if __name__ == "__main__":
    main()
