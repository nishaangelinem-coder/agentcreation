#!/usr/bin/env python3
"""System-level model: where should the multiply-accumulate engines live?

S1  Conventional, no HBF : weights in HBM, MAC array on the accelerator die.
S2  HBF, far compute     : weights in HBF, streamed over the interposer PHY to the accelerator.
S3  HBF, near-flash MAC  : weights consumed by MAC engines in the HBF base die; only
                           activations and partial sums cross the interposer.

Per-bit data-movement energies are explicit, swept assumptions (DEFAULTS). MAC energy,
area and clock come from the synthesis results (results/syn_results.json).
"""
import json, os, sys, math, itertools
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, ".."))
import hbf_model as hm

SYN = json.load(open(os.path.join(HERE, "results/syn_results.json")))

DEFAULTS = dict(
    e_hbm=4.0,      # pJ/bit  HBM DRAM access + I/O to accelerator (HBM2-class ~4 pJ/bit [25])
    e_nand=2.0,     # pJ/bit  NAND array sense + page-buffer read (assumption, swept 0.5-10)
    e_tsv=0.1,      # pJ/bit  core-die -> base-die TSV transfer
    e_phy=0.5,      # pJ/bit  base die -> accelerator over interposer (advanced-package link)
    e_sram=0.2,     # pJ/bit  staging into accelerator SRAM / operand delivery
    e_buf=0.05,     # pJ/bit  page buffer -> near-flash MAC operand delivery on the base die
    bw_stack=0.8e12,  # B/s external HBF read bandwidth per stack
    p_budget=5.0,   # W  power budget for compute on one HBF base die
    a_budget=20.0,  # mm^2 base-die area available for compute
    gamma=1.0,      # internal-array bandwidth / external interface bandwidth
)

MULTS8 = ["bw_array8", "bw_wallace8", "bw_dadda8", "bw_dadda8_ks", "booth8", "booth8_ks", "behav_mult8"]

def rec(design, lib, mode="area"):
    return SYN[f"{design}.{lib}.{mode}"]

def ff_overhead(lib):
    return max(r["seq_overhead_ps"] for r in SYN.values() if r["lib"] == lib)

def t_clk_ps(r):
    """Registered clock period: combinational delay + flop clk->Q + setup."""
    return r["comb_delay_ps"] + ff_overhead(r["lib"]) if r["seq_overhead_ps"] == 0 else r["t_clk_ps"]

def e_op(r, profile="gaussian"):
    """Energy per operation (J), dynamic at the activity profile plus leakage over one cycle."""
    t = t_clk_ps(r) * 1e-12
    return (r[f"e_cycle_J_{profile}"] + r["p_leak_W"] * t) / r["ops"]

def pe_variant(mult, lib, mode="area", profile="gaussian"):
    """INT8 MAC PE built around a given multiplier: measured PE minus its array multiplier
    plus the chosen multiplier (accumulator, registers and pass-through are shared)."""
    pe, arr, m = rec("int8_pe", lib, mode), rec("bw_array8", lib, mode), rec(mult, lib, mode)
    area = pe["area_um2"] - arr["area_um2"] + m["area_um2"]
    energy = e_op(pe, profile) - e_op(arr, profile) + e_op(m, profile)
    # cycle: max(multiplier stage, accumulate stage) ; accumulate stage = PE cycle when mult is array
    acc_stage = max(0.0, t_clk_ps(pe) - (arr["comb_delay_ps"] + ff_overhead(lib)))
    tclk = max(m["comb_delay_ps"] + ff_overhead(lib), acc_stage, 1.0)
    return dict(area_um2=area, e_mac_J=energy, t_clk_ps=tclk)

def fp8_variant(lib, mode="area", profile="gaussian"):
    r = rec("fp8_mac", lib, mode)
    return dict(area_um2=r["area_um2"], e_mac_J=e_op(r, profile), t_clk_ps=t_clk_ps(r))

# ------------------------------------------------------------------ energy per token
def move_energy_per_bit(scn, p):
    if scn == "S1": return p["e_hbm"] + p["e_sram"]
    if scn == "S2": return p["e_nand"] + p["e_tsv"] + p["e_phy"] + p["e_sram"]
    if scn == "S3": return p["e_nand"] + p["e_tsv"] + p["e_buf"]
    raise ValueError(scn)

def token_energy(model, B, scn, e_mac_J, p=DEFAULTS, d_model=7168):
    """Joules per generated token for the weight-GEMV part of one decode step (batch B)."""
    m = hm.MODELS[model]
    W = hm.weight_bytes_touched(m, B)                  # bytes read per step
    e_move = W * 8 * move_energy_per_bit(scn, p) * 1e-12
    e_comp = m["P_act"] * B * e_mac_J                   # one MAC per active weight per token
    e_act = 0.0
    if scn == "S3":   # per layer ~12*d_model elements cross the PHY: FP8 inputs in, FP32 partial sums out
        e_act = m["layers"] * B * 12 * d_model * (1 + 4) * 8 * p["e_phy"] * 1e-12
    return (e_move + e_comp + e_act) / B, dict(move=e_move / B, compute=e_comp / B, act=e_act / B)

# ------------------------------------------------------------------ near-flash capacity
def nmc_batch_limit(unit, p=DEFAULTS):
    """Largest decode batch whose MAC demand fits the base-die power and area budgets.
    MAC rate per stack = internal bandwidth (B/s) x B; units = rate / f_clk."""
    bw = p["bw_stack"] * p["gamma"]
    f = 1e12 / unit["t_clk_ps"]
    per_B_power = bw * unit["e_mac_J"]                       # W per unit of batch
    per_B_area = bw / f * unit["area_um2"] * 1e-6            # mm^2 per unit of batch
    return min(p["p_budget"] / per_B_power, p["a_budget"] / per_B_area), per_B_power, per_B_area

def decode_with_nmc(model, cfg, S, N, Bcap, gamma, slo=0.050):
    """Best tokens/s per package with weight traffic served at gamma x HBF bandwidth,
    restricted to batches the near-flash engines can sustain (B <= Bcap)."""
    old = hm.HBF_STACK["bw"]
    hm.HBF_STACK["bw"] = old * gamma
    try:
        best, B = None, 1
        while B <= max(1, int(Bcap)):
            r = hm.decode_step(model, cfg, B, S, N)
            if r is None: break
            if r["t"] <= slo and (best is None or r["thr"] > best["thr"]):
                best = dict(r, B=B)
            B += 1 if B < 16 else max(1, B // 8)
        return best
    finally:
        hm.HBF_STACK["bw"] = old

if __name__ == "__main__":
    out = {}
    for lib in ("sky130", "asap7"):
        for mlt in MULTS8:
            u = pe_variant(mlt, lib)
            bmax, pw, ar = nmc_batch_limit(u)
            out[f"pe|{lib}|{mlt}"] = dict(u, Bmax=bmax, W_per_B=pw, mm2_per_B=ar)
        u = fp8_variant(lib); bmax, pw, ar = nmc_batch_limit(u)
        out[f"pe|{lib}|fp8_mac"] = dict(u, Bmax=bmax, W_per_B=pw, mm2_per_B=ar)
    print(json.dumps(out, indent=1))
