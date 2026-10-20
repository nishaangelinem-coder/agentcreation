#!/usr/bin/env python3
"""Figures and table data for the near-flash multiplier study (Section 7 of the paper)."""
import json, os, sys
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch

HERE = os.path.dirname(os.path.abspath(__file__)); os.chdir(HERE)
sys.path.insert(0, HERE); sys.path.insert(0, "..")
import nmc_model as nm
import hbf_model as hm

FIG = "../figures"
plt.rcParams.update({"font.family": "DejaVu Serif", "font.size": 9, "axes.grid": True, "grid.alpha": .3,
                     "axes.spines.top": False, "axes.spines.right": False, "savefig.dpi": 300, "savefig.bbox": "tight"})
SYN = nm.SYN
STRUCT = json.load(open("results/structure_stats.json"))
LIBN = {"sky130": "SkyWater 130 nm", "asap7": "ASAP7 7 nm"}
ARCH8 = [("bw_array8", "Array (BW)"), ("bw_wallace8", "Wallace"), ("bw_dadda8", "Dadda"),
         ("bw_dadda8_ks", "Dadda+KS"), ("booth8", "Booth-4"), ("booth8_ks", "Booth-4+KS"), ("behav_mult8", "Tool a*b")]
COL = ["#1f4e79", "#2e86c1", "#16a085", "#7dcea0", "#e67e22", "#f5b041", "#7f8c8d"]

def metrics(r, profile="gaussian"):
    t = nm.t_clk_ps(r)
    e = nm.e_op(r, profile)
    eu = nm.e_op(r, "uniform")
    f = 1e12 / t
    return dict(area=r["area_um2"], cells=r["cells"], delay=r["comb_delay_ps"], tclk=t, fmax_ghz=f / 1e9,
                e_pj=e * 1e12, eu_pj=eu * 1e12, leak_uw=r["p_leak_W"] * 1e6,
                p_mw=(r[f"e_cycle_J_{profile}"] * f + r["p_leak_W"]) * 1e3,
                adp=r["area_um2"] * t * 1e-3,             # um^2 * ns
                edp=e * 1e12 * t * 1e-3,                  # pJ * ns
                gops_mm2=r["ops"] * f / 1e9 / (r["area_um2"] * 1e-6),
                tops_w=1.0 / (e * 1e12),                  # TOPS/W = 1 / (pJ per op); one MAC = one op
                lat=r["lat"])

tables = {}

# ---------------- Fig 7: base-die engine block diagram ----------------
fig, ax = plt.subplots(figsize=(7.6, 3.9)); ax.axis("off"); ax.set_xlim(0, 116); ax.set_ylim(0, 56)
def bx(x, y, w, h, t, fc, fs=7.2, ec="#333", ls="-"):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.3", fc=fc, ec=ec, lw=1, ls=ls))
    ax.text(x + w / 2, y + h / 2, t, ha="center", va="center", fontsize=fs)
def arr(x1, y1, x2, y2, c="#333", lw=1.1):
    ax.annotate("", xy=(x2, y2), xytext=(x1, y1), arrowprops=dict(arrowstyle="->", lw=lw, color=c))
bx(2, 44, 96, 9, "3D NAND core dies (16-high): ~1,000 independently sensed sub-arrays per stack\n(read-only weight regions, sealed after load)", "#f5cba7")
for i, x in enumerate([6, 30, 54, 78]):
    arr(x + 8, 44, x + 8, 39, "#c0392b")
ax.text(50, 41, "TSV array (≈0.1 pJ/bit)", ha="center", fontsize=7, color="#c0392b")
ax.add_patch(FancyBboxPatch((2, 2), 96, 35, boxstyle="round,pad=0.4", fc="#f4f6f7", ec="#555", lw=1.2, ls="--"))
ax.text(50, 35.2, "HBF logic base die", ha="center", fontsize=8.5, style="italic")
for i, x in enumerate([5, 29, 53, 77]):
    bx(x, 27, 18, 5.5, f"Page buffers +\nLDPC decode (ch {4*i}-{4*i+3})", "#fdebd0", fs=6.5)
    arr(x + 9, 27, x + 9, 22.5, "#c0392b")
    bx(x, 12.5, 18, 9.5, "Systolic MAC tile\nR×C PEs\n(array / Booth\nmultipliers)", "#d6eaf8", fs=6.5)
    arr(x + 9, 12.5, x + 9, 9.5)
bx(5, 4, 42, 5, "Activation SRAM (x, broadcast to rows)", "#e8f8f5", fs=6.8)
bx(53, 4, 42, 5, "Partial-sum buffer + FP32 reduce", "#e8f8f5", fs=6.8)
for x in [14, 38]:
    arr(26, 9, x, 12.5, "#16a085", .8)
ax.annotate("", xy=(104, 6.5), xytext=(96, 6.5), arrowprops=dict(arrowstyle="<->", lw=1.2, color="#1f4e79"))
ax.text(105, 6.5, "to accelerator\n(interposer):\nactivations in,\npartial sums out", fontsize=6.5, ha="left", va="center", color="#1f4e79")
fig.savefig(f"{FIG}/fig7_nmc_engine.png"); plt.close(fig)

# ---------------- Fig 8: 8-bit architecture comparison (both nodes) ----------------
fig, axs = plt.subplots(2, 3, figsize=(7.4, 4.6))
keys = [("area", "Area (µm²)"), ("tclk", "Registered clock period (ps)"), ("e_pj", "Energy / multiply (pJ)"),
        ("adp", "Area-delay product (µm²·ns)"), ("edp", "Energy-delay product (pJ·ns)"), ("tops_w", "Efficiency (TOPS/W)")]
for ax, (k, lab) in zip(axs.flat, keys):
    for j, lib in enumerate(("sky130", "asap7")):
        base = metrics(nm.rec("bw_array8", lib))[k]
        vals = [metrics(nm.rec(d, lib))[k] / base for d, _ in ARCH8]
        ax.bar(np.arange(len(ARCH8)) + (j - .5) * .38, vals, .38, color=["#34495e", "#e67e22"][j], label=LIBN[lib])
    ax.axhline(1, c="k", lw=.7, ls=":")
    ax.set_xticks(range(len(ARCH8))); ax.set_xticklabels([n for _, n in ARCH8], rotation=55, ha="right", fontsize=6.8)
    ax.set_title(lab + "\n(normalised to array)", fontsize=7.6)
axs[0, 0].legend(fontsize=6.5)
fig.tight_layout(); fig.savefig(f"{FIG}/fig8_mult_compare.png"); plt.close(fig)

# ---------------- Fig 9: width scaling ----------------
fams = [("bw_array", "Array (BW)"), ("bw_wallace", "Wallace"), ("bw_dadda", "Dadda"), ("booth", "Booth-4"), ("behav_mult", "Tool a*b")]
fig, axs = plt.subplots(1, 3, figsize=(7.4, 2.7))
for ax, (k, lab) in zip(axs, [("area", "Area (µm²)"), ("delay", "Critical path (ps)"), ("e_pj", "Energy / multiply (pJ)")]):
    for (f, n), c in zip(fams, COL[:3] + COL[4:5] + COL[6:]):
        xs = [4, 8, 16]
        ax.plot(xs, [metrics(nm.rec(f"{f}{w}", "asap7"))[k] for w in xs], "o-", c=c, lw=1.4, ms=3.5, label=n)
    ax.set_xscale("log", base=2); ax.set_xticks([4, 8, 16]); ax.set_xticklabels(["4", "8", "16"])
    ax.set_yscale("log"); ax.set_xlabel("Operand width n (bits)"); ax.set_title(lab + " — ASAP7", fontsize=8)
axs[0].legend(fontsize=6.3)
fig.tight_layout(); fig.savefig(f"{FIG}/fig9_width_scaling.png"); plt.close(fig)

# ---------------- Fig 10: energy-delay Pareto ----------------
fig, axs = plt.subplots(1, 2, figsize=(7.4, 3.0))
group = {"bw_array": ("Array", "o"), "bw_wallace": ("Wallace", "s"), "bw_dadda": ("Dadda", "^"), "booth": ("Booth-4", "D"), "behav_mult": ("Tool a*b", "x")}
for ax, lib in zip(axs, ("sky130", "asap7")):
    for (pre, (lab, mk)), c in zip(group.items(), [COL[0], COL[1], COL[2], COL[4], COL[6]]):
        pts = []
        for key, r in SYN.items():
            if r["lib"] != lib or r["n"] != 8 or not r["design"].startswith(pre):
                continue
            if r["design"] in ("int8_pe", "fp8_mac", "systolic_4x8", "u_array4"): continue
            m = metrics(r); pts.append((m["tclk"], m["e_pj"], r["design"], r["mode"]))
        if pts:
            ax.scatter([p[0] for p in pts], [p[1] for p in pts], marker=mk, c=c, s=22, label=lab)
            for p in pts:
                if "_p" in p[2] and p[3] == "area":
                    ax.annotate("pipelined", (p[0], p[1]), fontsize=5.8, xytext=(3, 3), textcoords="offset points")
    ax.set_xlabel("Registered clock period (ps)"); ax.set_ylabel("Energy / multiply (pJ)"); ax.set_yscale("log")
    ax.set_title(f"8-bit multipliers, {LIBN[lib]}\n(area- and delay-driven synthesis)", fontsize=8)
axs[0].legend(fontsize=6.3)
fig.tight_layout(); fig.savefig(f"{FIG}/fig10_pareto.png"); plt.close(fig)

# ---------------- Fig 11: operand statistics ----------------
fig, axs = plt.subplots(1, 2, figsize=(7.4, 2.7))
for ax, lib in zip(axs, ("sky130", "asap7")):
    u = [metrics(nm.rec(d, lib))["eu_pj"] for d, _ in ARCH8]
    g = [metrics(nm.rec(d, lib))["e_pj"] for d, _ in ARCH8]
    x = np.arange(len(ARCH8))
    ax.bar(x - .2, u, .4, color="#95a5a6", label="Uniform random operands")
    ax.bar(x + .2, g, .4, color="#c0392b", label="Gaussian (LLM-weight-like)")
    ax.set_xticks(x); ax.set_xticklabels([n for _, n in ARCH8], rotation=55, ha="right", fontsize=6.8)
    ax.set_ylabel("Energy / multiply (pJ)"); ax.set_title(LIBN[lib], fontsize=8)
axs[0].legend(fontsize=6.3)
fig.tight_layout(); fig.savefig(f"{FIG}/fig11_operand_stats.png"); plt.close(fig)

# ---------------- table data: multipliers ----------------
order = [d for d, _ in ARCH8] + ["bw_array8_p3", "bw_dadda8_p2", "booth8_p2"]
for lib in ("sky130", "asap7"):
    rows = []
    for d in order:
        for mode in ("area", "delay"):
            r = nm.rec(d, lib, mode); m = metrics(r)
            rows.append(dict(design=d, arch=r["arch"], mode=mode, **m))
    tables[f"mult8_{lib}"] = rows
w_rows = []
for f, n in fams:
    for w in (4, 8, 16):
        for lib in ("sky130", "asap7"):
            m = metrics(nm.rec(f"{f}{w}", lib)); w_rows.append(dict(fam=n, n=w, lib=lib, **m))
for d in ("bw_dadda16_ks", "booth16_ks", "bw_array16_p4"):
    for lib in ("sky130", "asap7"):
        r = nm.rec(d, lib); m = metrics(r); w_rows.append(dict(fam=r["arch"], n=16, lib=lib, **m))
tables["width"] = w_rows
mac_rows = []
for d in ("u_array4", "int8_pe", "fp8_mac", "systolic_4x8"):
    for lib in ("sky130", "asap7"):
        if f"{d}.{lib}.area" not in SYN: continue
        r = nm.rec(d, lib); m = metrics(r); mac_rows.append(dict(design=d, arch=r["arch"], lib=lib, ops=r["ops"], **m))
tables["macs"] = mac_rows
tables["struct"] = STRUCT

# ---------------- Fig 12: pipelining ----------------
pairs = [("bw_array8", "bw_array8_p3", "Array 8b\n3-stage"), ("bw_dadda8", "bw_dadda8_p2", "Dadda 8b\n2-stage"),
         ("booth8", "booth8_p2", "Booth 8b\n2-stage"), ("bw_array16", "bw_array16_p4", "Array 16b\n4-stage")]
fig, axs = plt.subplots(1, 3, figsize=(7.4, 2.6))
for j, lib in enumerate(("sky130", "asap7")):
    sp = [metrics(nm.rec(b, lib))["fmax_ghz"] / metrics(nm.rec(a, lib))["fmax_ghz"] for a, b, _ in pairs]
    ar = [metrics(nm.rec(b, lib))["area"] / metrics(nm.rec(a, lib))["area"] for a, b, _ in pairs]
    en = [metrics(nm.rec(b, lib))["e_pj"] / metrics(nm.rec(a, lib))["e_pj"] for a, b, _ in pairs]
    for ax, v in zip(axs, (sp, ar, en)):
        ax.bar(np.arange(len(pairs)) + (j - .5) * .38, v, .38, color=["#34495e", "#e67e22"][j], label=LIBN[lib])
for ax, t in zip(axs, ("Fmax gain (×)", "Area overhead (×)", "Energy / multiply (×)")):
    ax.axhline(1, c="k", lw=.7, ls=":"); ax.set_xticks(range(len(pairs))); ax.set_xticklabels([p[2] for p in pairs], fontsize=6.5)
    ax.set_title(t, fontsize=8)
axs[0].legend(fontsize=6.3)
fig.tight_layout(); fig.savefig(f"{FIG}/fig12_pipelining.png"); plt.close(fig)

# ---------------- System level: S1 / S2 / S3 ----------------
P = dict(nm.DEFAULTS)
acc = nm.pe_variant("bw_array8", "asap7")          # accelerator-die MAC (advanced node)
sys_rows = []
for model in ("Dense-70B", "MoE-1T"):
    for B in (1, 8, 64):
        for scn, unit, lab in [("S1", acc, "S1: HBM, MAC on accelerator (no HBF)"),
                               ("S2", acc, "S2: HBF, MAC on accelerator"),
                               ("S3", nm.pe_variant("bw_array8", "sky130"), "S3: HBF + near-flash MAC (130 nm base die)"),
                               ("S3", nm.pe_variant("bw_array8", "asap7"), "S3: HBF + near-flash MAC (7 nm base die)")]:
            e, parts = nm.token_energy(model, B, scn, unit["e_mac_J"], P)
            sys_rows.append(dict(model=model, B=B, scn=lab, e_mJ=e * 1e3, **{k: v * 1e3 for k, v in parts.items()}))
tables["sys_energy"] = sys_rows

fig, axs = plt.subplots(1, 2, figsize=(7.4, 3.1))
for ax, model in zip(axs, ("Dense-70B", "MoE-1T")):
    rows = [r for r in sys_rows if r["model"] == model]
    labs = sorted(set(r["scn"] for r in rows), key=lambda s: [r["scn"] for r in rows].index(s))
    x = np.arange(3); w = .2
    for j, lab in enumerate(labs):
        rr = [r for r in rows if r["scn"] == lab]
        mv = [r["move"] for r in rr]; cp = [r["compute"] for r in rr]; ac = [r["act"] for r in rr]
        ax.bar(x + (j - 1.5) * w, mv, w, color=COL[j], label=lab.split(":")[0] + (": " + lab.split(":")[1] if True else ""))
        ax.bar(x + (j - 1.5) * w, cp, w, bottom=mv, color=COL[j], hatch="////", ec="w", lw=0)
    ax.set_yscale("log"); ax.set_xticks(x); ax.set_xticklabels(["B = 1", "B = 8", "B = 64"])
    ax.set_ylabel("Weight-GEMV energy per token (mJ)"); ax.set_title(model, fontsize=8.5)
axs[1].legend(fontsize=5.8, loc="upper right")
fig.tight_layout(); fig.savefig(f"{FIG}/fig13_system_energy.png"); plt.close(fig)

# near-flash capacity: Bmax per multiplier/node
cap_rows = []
for lib in ("sky130", "asap7"):
    for d, n in ARCH8 + [("fp8_mac", "FP8 MAC")]:
        u = nm.fp8_variant(lib) if d == "fp8_mac" else nm.pe_variant(d, lib)
        bmax, pw, ar = nm.nmc_batch_limit(u, P)
        cap_rows.append(dict(lib=lib, design=d, name=n, Bmax=bmax, W_per_B=pw, mm2_per_B=ar,
                             e_mac_fJ=u["e_mac_J"] * 1e15, area=u["area_um2"], tclk=u["t_clk_ps"],
                             Bmax_power=P["p_budget"] / pw, Bmax_area=P["a_budget"] / ar))
tables["nmc_capacity"] = cap_rows

fig, axs = plt.subplots(1, 2, figsize=(7.4, 3.0))
ax = axs[0]
for j, lib in enumerate(("sky130", "asap7")):
    rr = [r for r in cap_rows if r["lib"] == lib]
    ax.bar(np.arange(len(rr)) + (j - .5) * .38, [r["Bmax"] for r in rr], .38, color=["#34495e", "#e67e22"][j], label=LIBN[lib])
ax.set_yscale("log"); ax.set_xticks(range(len(ARCH8) + 1)); ax.set_xticklabels([n for _, n in ARCH8] + ["FP8 MAC"], rotation=55, ha="right", fontsize=6.6)
ax.set_ylabel("Max. decode batch served in base die"); ax.set_title("(a) 5 W, 20 mm² per stack, 0.8 TB/s", fontsize=8); ax.legend(fontsize=6.3)
ax = axs[1]
pb = np.linspace(1, 15, 30)
for lib, ls in (("sky130", "--"), ("asap7", "-")):
    for d, c in (("bw_array8", COL[0]), ("booth8", COL[4]), ("fp8_mac", "#8e44ad")):
        u = nm.fp8_variant(lib) if d == "fp8_mac" else nm.pe_variant(d, lib)
        ys = [nm.nmc_batch_limit(u, dict(P, p_budget=x))[0] for x in pb]
        ax.plot(pb, ys, ls, c=c, lw=1.3, label=f"{dict(ARCH8 + [('fp8_mac', 'FP8 MAC')])[d]}, {lib}")
ax.set_yscale("log"); ax.set_xlabel("Base-die compute power budget (W)"); ax.set_ylabel("Max. batch")
ax.set_title("(b) Sensitivity to thermal budget", fontsize=8); ax.legend(fontsize=5.6, ncol=2)
fig.tight_layout(); fig.savefig(f"{FIG}/fig14_nmc_capacity.png"); plt.close(fig)

# sensitivity heatmaps: S3 / S2 and S3 / S1 energy at B = 1 (MoE-1T)
en = np.array([0.5, 1, 2, 4, 6, 10]); ep = np.array([0.25, 0.5, 1, 1.5, 2, 3])
u3 = nm.pe_variant("bw_array8", "asap7")
fig, axs = plt.subplots(1, 2, figsize=(7.4, 3.0))
heat = {}
for ax, ref in zip(axs, ("S2", "S1")):
    Z = np.zeros((len(en), len(ep)))
    for i, a in enumerate(en):
        for j, b in enumerate(ep):
            p = dict(P, e_nand=a, e_phy=b)
            e3 = nm.token_energy("MoE-1T", 1, "S3", u3["e_mac_J"], p)[0]
            er = nm.token_energy("MoE-1T", 1, ref, acc["e_mac_J"], p)[0]
            Z[i, j] = e3 / er
    heat[ref] = Z.tolist()
    im = ax.imshow(Z, origin="lower", cmap="RdYlGn_r", vmin=0, vmax=max(1.5, Z.max()), aspect="auto")
    for i in range(len(en)):
        for j in range(len(ep)):
            ax.text(j, i, f"{Z[i, j]:.2f}", ha="center", va="center", fontsize=6.3)
    ax.set_xticks(range(len(ep))); ax.set_xticklabels(ep); ax.set_yticks(range(len(en))); ax.set_yticklabels(en)
    ax.set_xlabel("Interposer PHY energy (pJ/bit)"); ax.set_ylabel("NAND array read energy (pJ/bit)")
    ax.set_title(f"S3 / {ref} energy per token (MoE-1T, B = 1)", fontsize=8); ax.grid(False)
fig.tight_layout(); fig.savefig(f"{FIG}/fig15_sensitivity.png"); plt.close(fig)
tables["heat"] = dict(e_nand=en.tolist(), e_phy=ep.tolist(), **heat)

# decode throughput with near-flash compute vs internal bandwidth gamma
thr_rows = []
gammas = [1, 1.5, 2, 3, 4]
fig, ax = plt.subplots(figsize=(6.0, 3.1))
base = {}
for N in (1, 4):
    r8 = hm.best_under_slo("MoE-1T", "HBM-8", 32768, N, .050)
    base[N] = r8["thr"] if r8 else 0
series = [("HBM-4/HBF-4", 1, "asap7", COL[0], "-"), ("HBM-4/HBF-4", 4, "asap7", COL[4], "-"),
          ("HBM-2/HBF-6", 4, "asap7", "#c0392b", "-"), ("HBM-4/HBF-4", 4, "sky130", COL[4], "--")]
for cfg, N, lib, c, ls in series:
    u = nm.pe_variant("bw_array8", lib)
    ys = []
    for g in gammas:
        bcap = nm.nmc_batch_limit(u, dict(P, gamma=g))[0]
        r = nm.decode_with_nmc("MoE-1T", cfg, 32768, N, bcap, g)
        ys.append(r["thr"] if r else 0)
        thr_rows.append(dict(cfg=cfg, N=N, lib=lib, gamma=g, Bcap=bcap, thr=r["thr"] if r else 0, B=r["B"] if r else 0,
                             t_ms=r["t"] * 1e3 if r else None))
    ax.plot(gammas, ys, "o" + ls, c=c, lw=1.4, ms=3.5, label=f"{cfg}, N={N}, MACs in {LIBN[lib]} base die")
ax.axhline(base[4], c="k", ls=":", lw=.9); ax.text(1.02, base[4] * 1.08, f"HBM-8, N=4 (no HBF): {base[4]:.0f} tok/s/pkg", fontsize=6.5)
ax.set_xlabel("Internal-array / external bandwidth ratio γ"); ax.set_ylabel("Tokens/s per package @ 50 ms")
ax.set_yscale("log"); ax.set_title("MoE-1T decode with near-flash MAC engines (32k context)", fontsize=8.5); ax.legend(fontsize=6)
fig.savefig(f"{FIG}/fig16_nmc_throughput.png"); plt.close(fig)
tables["nmc_thr"] = thr_rows
tables["nmc_thr_base"] = base
tables["defaults"] = P
tables["ff_overhead"] = {lib: nm.ff_overhead(lib) for lib in ("sky130", "asap7")}
json.dump(tables, open("results/hw_tables.json", "w"), indent=1)
print("figures + tables written")
