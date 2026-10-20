"""Analytical model for HBM / HBF memory configurations in LLM inference.
All device parameters are stated assumptions (see Table 2 of the paper)."""
import json, math

GB, TB = 1e9, 1e12
PEAK = 1.5e15            # effective FP8 FLOP/s per package (2.5 PFLOP/s peak x 0.6 MFU)
HBM_STACK = dict(cap=36*GB, bw=1.0*TB)
HBF_STACK = dict(cap=512*GB, bw=0.8*TB)
USABLE = 0.9             # fraction of capacity usable for weights + KV

CONFIGS = {              # (HBM stacks, HBF stacks) on an 8-site package
    "HBM-8":        (8, 0),
    "HBM-6/HBF-2":  (6, 2),
    "HBM-4/HBF-4":  (4, 4),
    "HBM-2/HBF-6":  (2, 6),
}

MODELS = {
    "Dense-70B": dict(P=70e9, P_act=70e9, shared=70e9, E=1, k=1, layers=80,
                      kvb=2*80*8*128*1, attn_dim=64*128),
    "MoE-1T":    dict(P=1000e9, P_act=32e9, shared=12e9, E=384, k=8, layers=61,
                      kvb=576*61*2, attn_dim=128*128),
}

def tiers(cfg):
    h, f = CONFIGS[cfg]
    return (h*HBM_STACK["cap"], h*HBM_STACK["bw"], f*HBF_STACK["cap"], f*HBF_STACK["bw"])

def weight_bytes_touched(m, B):
    if m["E"] == 1:
        return m["P"]
    frac = 1 - (1 - m["k"]/m["E"])**B
    return m["shared"] + (m["P"] - m["shared"]) * frac

ALPHA = 10e-6            # per-layer collective latency (s) when sharded over N>1 packages

def decode_step(mname, cfg, B, S, N=1):
    """One decode iteration for batch B at context S, sharded over N packages."""
    m = MODELS[mname]
    hc, hb, fc, fb = tiers(cfg)
    W = weight_bytes_touched(m, B) / N
    KV = B*S*m["kvb"] / N
    flops = (2*m["P_act"]*B + 4*S*m["attn_dim"]*m["layers"]*B) / N
    t_comm = 2*m["layers"]*ALPHA if N > 1 else 0.0
    if fc == 0:
        if m["P"]/N + KV > USABLE*hc: return None
        t_m, t_c = (W + KV)/hb, flops/PEAK
        t = max(t_m, t_c) + t_comm
        return dict(t=t, thr=B/t/N, bound="HBM" if t_m >= t_c else "compute")
    if m["P"]/N > USABLE*fc or KV > USABLE*hc: return None
    parts = dict(HBF=W/fb, HBM=KV/hb, compute=flops/PEAK)
    b = max(parts, key=parts.get)
    t = parts[b] + t_comm
    return dict(t=t, thr=B/t/N, bound=b)

def best_under_slo(mname, cfg, S, N, slo):
    best = None
    B = 1
    while B <= 8192:
        r = decode_step(mname, cfg, B, S, N)
        if r is None: break
        if r["t"] <= slo and (best is None or r["thr"] > best["thr"]):
            best = dict(r, B=B)
        B = int(B*1.25)+1 if B > 4 else B+1
    return best

def max_batch(mname, cfg, S):
    m = MODELS[mname]; hc, hb, fc, fb = tiers(cfg)
    if fc == 0: return None
    return int(USABLE*hc // (S*m["kvb"]))

def prefill(mname, cfg, L):
    m = MODELS[mname]; hc, hb, fc, fb = tiers(cfg)
    flops = 2*m["P_act"]*L + 2*L*L*m["attn_dim"]*m["layers"]
    W = weight_bytes_touched(m, L)
    bw = fb if fc else hb
    return dict(t_comp=flops/PEAK, t_mem=W/bw, Lstar=PEAK*W/(bw*2*m["P_act"]))

def lifetime_days(cap, pe, rate, waf=1.0):
    return pe*cap/(rate*waf*86400)

def parallel_units(bw, tR, page):
    return bw*tR/page

if __name__ == "__main__":
    for S in [8192, 32768]:
        for mn in MODELS:
            for N in [1, 2, 4, 8, 16]:
                for cfg in CONFIGS:
                    r = best_under_slo(mn, cfg, S, N, 0.050)
                    print(S, mn, N, cfg, None if r is None else (r["B"], round(r["t"]*1e3,1), round(r["thr"]), r["bound"]))
    for mn in MODELS:
        for cfg in ["HBM-8", "HBM-4/HBF-4", "HBM-2/HBF-6"]:
            print("prefill", mn, cfg, {k: round(v,4) for k,v in prefill(mn, cfg, 8192).items()})
