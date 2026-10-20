import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
import numpy as np, copy, json
import hbf_model as hm

plt.rcParams.update({"font.family": "DejaVu Serif", "font.size": 9, "axes.grid": True,
                     "grid.alpha": .3, "axes.spines.top": False, "axes.spines.right": False,
                     "savefig.dpi": 300, "savefig.bbox": "tight"})
C = {"HBM-8": "#1f4e79", "HBM-6/HBF-2": "#2e86c1", "HBM-4/HBF-4": "#e67e22", "HBM-2/HBF-6": "#c0392b"}
results = {}

# ---------- Fig 1: architecture ----------
fig, ax = plt.subplots(figsize=(7.2, 3.6)); ax.axis("off"); ax.set_xlim(0, 100); ax.set_ylim(-4, 50)
def box(x, y, w, h, t, fc, ec="#333", fs=8, bold=False):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.3", fc=fc, ec=ec, lw=1))
    ax.text(x+w/2, y+h/2, t, ha="center", va="center", fontsize=fs, weight="bold" if bold else "normal")
ax.add_patch(FancyBboxPatch((2, 6), 62, 40, boxstyle="round,pad=0.4", fc="#f4f6f7", ec="#555", lw=1.2, ls="--"))
ax.text(33, 47.5, "2.5D package (silicon interposer / bridge)", ha="center", fontsize=8.5, style="italic")
box(21.5, 17, 23, 18, "Accelerator die\nTensor cores + SRAM\nHBM PHY | HBF PHY\nflash-aware DMA and\nprefetch engine", "#d6eaf8", fs=7.5)
for i, y in enumerate([36, 26]):
    box(5, y-2, 13, 8, "HBM stack\n(KV cache,\nactivations)", "#aed6f1", fs=7)
for i, y in enumerate([36, 26]):
    box(47.5, y-2, 14.5, 8, "HBF stack\n(weights,\nprefix cache)", "#f5cba7", fs=7)
box(5, 9, 13, 8, "HBM stack", "#aed6f1", fs=7); box(47.5, 9, 14.5, 8, "HBF stack", "#f5cba7", fs=7)
for y in [38, 28, 13]:
    ax.annotate("", xy=(22, min(max(y, 20), 32)), xytext=(18, y), arrowprops=dict(arrowstyle="<->", lw=1, color="#1f4e79"))
    ax.annotate("", xy=(44, min(max(y, 20), 32)), xytext=(47.5, y), arrowprops=dict(arrowstyle="<-", lw=1, color="#c0392b"))
ax.text(3, 1.5, "HBM: ~1 TB/s per stack, unlimited R/W", fontsize=7, color="#1f4e79", ha="left")
ax.text(3, -2.0, "HBF: ~0.8 TB/s read per stack, read-mostly", fontsize=7, color="#c0392b", ha="left")
box(70, 36, 27, 8, "Host CPU + LPDDR/DDR\n(scheduler, tokenizer)", "#e8f8f5", fs=7.5)
box(70, 24, 27, 8, "CXL memory pool\n(KV spill, warm prefixes)", "#e8f8f5", fs=7.5)
box(70, 12, 27, 8, "NVMe SSD / object store\n(checkpoints, cold KV)", "#e8f8f5", fs=7.5)
for y in [40, 28, 16]:
    ax.annotate("", xy=(64.5, 26), xytext=(70, y), arrowprops=dict(arrowstyle="<->", lw=.8, color="#555"))
ax.text(83.5, 47, "Off-package tiers", ha="center", fontsize=8.5, style="italic")
fig.savefig("figures/fig1_architecture.png"); plt.close(fig)

# ---------- Fig 2: roofline ----------
fig, ax = plt.subplots(figsize=(6.0, 3.6))
I = np.logspace(-0.5, 4.5, 400)
for bw, lab, c in [(8e12, "HBM-8 (8 TB/s)", C["HBM-8"]), (4.8e12, "HBF-6 tier (4.8 TB/s)", C["HBM-2/HBF-6"]),
                   (3.2e12, "HBF-4 tier (3.2 TB/s)", C["HBM-4/HBF-4"])]:
    ax.loglog(I, np.minimum(hm.PEAK, bw*I)/1e12, c=c, lw=1.6, label=lab)
    ax.axvline(hm.PEAK/bw, c=c, ls=":", lw=.8)
pts = {"Decode B=1": 2, "Decode B=32": 64, "Decode B=256": 512, "Prefill L=512": 1024, "Prefill L=8192": 16384}
offs = {"Decode B=1": (4, -11), "Decode B=32": (4, -11), "Decode B=256": (-20, -30), "Prefill L=512": (4, 8), "Prefill L=8192": (-30, -14)}
for k, v in pts.items():
    ax.scatter([v], [min(hm.PEAK, 3.2e12*v)/1e12], c="k", s=14, zorder=5)
    ax.annotate(k, (v, min(hm.PEAK, 3.2e12*v)/1e12), textcoords="offset points", xytext=offs[k], fontsize=7)
ax.set_xlabel("Arithmetic intensity (FLOP / weight byte, FP8)"); ax.set_ylabel("Attainable TFLOP/s")
ax.set_title("Roofline of a 1.5 PFLOP/s package with HBM and HBF ceilings", fontsize=9)
ax.legend(fontsize=7, loc="lower right"); fig.savefig("figures/fig2_roofline.png"); plt.close(fig)

# ---------- Fig 3: MoE decode TPOT vs batch + expert touch ----------
fig, (a1, a2) = plt.subplots(1, 2, figsize=(7.2, 3.0))
Bs = np.unique(np.logspace(0, 3, 40).astype(int))
for cfg in hm.CONFIGS:
    ts = [(b, hm.decode_step("MoE-1T", cfg, b, 32768, 4)) for b in Bs]
    ts = [(b, r["t"]*1e3) for b, r in ts if r]
    if ts: a1.plot(*zip(*ts), c=C[cfg], lw=1.6, label=cfg)
a1.axhline(50, c="k", ls="--", lw=.8); a1.text(1.2, 54, "TPOT SLO = 50 ms", fontsize=7)
a1.set_xscale("log"); a1.set_yscale("log"); a1.set_xlabel("Decode batch size B"); a1.set_ylabel("Time per output token (ms)")
a1.set_title("(a) MoE-1T, N = 4 packages, S = 32k", fontsize=8.5); a1.legend(fontsize=6.5)
m = hm.MODELS["MoE-1T"]
for E, k, lab in [(384, 8, "E=384, k=8 (MoE-1T)"), (64, 4, "E=64, k=4"), (128, 2, "E=128, k=2")]:
    a2.plot(Bs, 100*(1-(1-k/E)**Bs), lw=1.6, label=lab)
a2.set_xscale("log"); a2.set_xlabel("Decode batch size B"); a2.set_ylabel("Experts touched per layer (%)")
a2.set_title("(b) Expert-weight traffic grows with batch", fontsize=8.5); a2.legend(fontsize=6.5)
fig.tight_layout(); fig.savefig("figures/fig3_moe_decode.png"); plt.close(fig)

# ---------- Fig 4: best throughput under SLO vs N ----------
fig, axs = plt.subplots(1, 2, figsize=(7.2, 3.0), sharey=False)
Ns = [1, 2, 4, 8, 16]
for ax, S in zip(axs, [8192, 32768]):
    w = .2
    for j, cfg in enumerate(hm.CONFIGS):
        vals = []
        for N in Ns:
            r = hm.best_under_slo("MoE-1T", cfg, S, N, .050)
            vals.append(r["thr"] if r else 0)
            results[f"slo|{S}|{N}|{cfg}"] = None if r is None else dict(B=r["B"], t=round(r["t"]*1e3, 1), thr=round(r["thr"]), bound=r["bound"])
        ax.bar(np.arange(len(Ns))+(j-1.5)*w, vals, w, color=C[cfg], label=cfg)
    ax.set_xticks(range(len(Ns))); ax.set_xticklabels([str(n) for n in Ns]); ax.set_xlabel("Packages N")
    ax.set_ylabel("Tokens/s per package @ 50 ms"); ax.set_yscale("log"); ax.set_ylim(30, 2e4)
    ax.set_title(f"MoE-1T, context S = {S//1024}k", fontsize=8.5)
axs[0].legend(fontsize=6.5, loc="upper left"); fig.tight_layout(); fig.savefig("figures/fig4_slo_throughput.png"); plt.close(fig)

# sensitivity: HBF bandwidth 1.6 TB/s/stack and shared-weights-in-HBM
for bw in [0.8, 1.6]:
    hm.HBF_STACK["bw"] = bw*hm.TB
    for cfg in ["HBM-4/HBF-4", "HBM-2/HBF-6"]:
        for N in [1, 4]:
            r = hm.best_under_slo("MoE-1T", cfg, 32768, N, .050)
            results[f"sens|{bw}|{cfg}|{N}"] = None if r is None else dict(B=r["B"], t=round(r["t"]*1e3, 1), thr=round(r["thr"]), bound=r["bound"])
hm.HBF_STACK["bw"] = 0.8*hm.TB

# ---------- Fig 5: endurance ----------
fig, ax = plt.subplots(figsize=(6.0, 3.4))
rates = np.logspace(-3, 2, 200)  # GB/s per stack
for pe, lab, c in [(3000, "TLC, 3k P/E", "#c0392b"), (30000, "pSLC, 30k P/E (1/3 capacity)", "#e67e22"), (100000, "SLC-class, 100k P/E (1/3 capacity)", "#2e86c1")]:
    cap = 512*hm.GB if pe == 3000 else 512*hm.GB/3
    ax.loglog(rates, hm.lifetime_days(cap, pe, rates*hm.GB, waf=2.0)/365, c=c, lw=1.6, label=lab)
ax.axhline(5, c="k", ls="--", lw=.8); ax.text(1.2e-3, 6, "5-year service life", fontsize=7)
for x, t in [(0.003, "weight refresh\n1x/day"), (1.0, "KV write\n(sustained)"), (20, "KV swap-in\n(burst)")]:
    ax.axvline(x, c="#777", ls=":", lw=.8); ax.text(x*1.1, 2e-3, t, fontsize=6.5)
ax.set_xlabel("Sustained write rate per HBF stack (GB/s), WAF = 2"); ax.set_ylabel("Lifetime (years)")
ax.set_title("Endurance-limited lifetime of a 512 GB HBF stack", fontsize=9); ax.legend(fontsize=7)
fig.savefig("figures/fig5_endurance.png"); plt.close(fig)
for pe, cap in [(3000, 512), (30000, 512/3)]:
    for r in [0.003, 1.0, 20.0]:
        results[f"life|{pe}|{r}"] = round(hm.lifetime_days(cap*hm.GB, pe, r*hm.GB, 2.0), 1)

# ---------- Fig 6: parallel read units ----------
fig, ax = plt.subplots(figsize=(6.0, 3.4))
tR = np.logspace(0, 2, 100)*1e-6
for page, c in [(4096, "#2e86c1"), (16384, "#e67e22"), (65536, "#c0392b")]:
    ax.loglog(tR*1e6, hm.parallel_units(0.8e12, tR, page), c=c, lw=1.6, label=f"{page//1024} KiB read unit")
ax.set_xlabel("Array read latency tR (µs)"); ax.set_ylabel("Concurrent reads in flight per stack")
ax.set_title("Little's law: concurrency needed to sustain 0.8 TB/s per stack", fontsize=9)
ax.legend(fontsize=7); fig.savefig("figures/fig6_littles_law.png"); plt.close(fig)
for page in [4096, 16384]:
    for t in [3, 20, 50]:
        results[f"little|{page}|{t}"] = round(hm.parallel_units(0.8e12, t*1e-6, page))

json.dump(results, open("results.json", "w"), indent=1)
print(json.dumps(results, indent=1))
