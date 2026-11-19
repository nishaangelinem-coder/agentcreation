#!/usr/bin/env python3
"""Analyse a measurement campaign: statistics, run-time verification verdicts, LaTeX tables, figures.
usage: analyze.py [campaign_dir] [paper_dir]"""
import glob, json, math, os, re, sys
import numpy as np, pandas as pd
from scipy import stats
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

CAMP = sys.argv[1] if len(sys.argv) > 1 else "results/campaign"
PAPER = sys.argv[2] if len(sys.argv) > 2 else "paper"
FS, DEADLINE_US = 500, 10000.0
os.makedirs(f"{PAPER}/tables", exist_ok=True); os.makedirs(f"{PAPER}/figs", exist_ok=True)
plt.rcParams.update({"font.size": 9, "font.family": "serif", "axes.spines.top": False, "axes.spines.right": False,
                     "axes.grid": True, "grid.color": "#e4e3df", "grid.linewidth": 0.6, "axes.axisbelow": True,
                     "legend.frameon": False, "pdf.fonttype": 42, "lines.linewidth": 1.6})
CFG = [("fifo", "event"), ("fifo", "poll"), ("other", "event"), ("other", "poll")]
LBL = {("fifo", "event"): "FP + event (proposed)", ("fifo", "poll"): "FP + polling",
       ("other", "event"): "TS + event", ("other", "poll"): "TS + polling"}
COL = {("fifo", "event"): "#2a78d6", ("fifo", "poll"): "#1baf7a", ("other", "event"): "#eb6834", ("other", "poll"): "#eda100"}
LS = {("fifo", "event"): "-", ("fifo", "poll"): "--", ("other", "event"): "-.", ("other", "poll"): ":"}

def runs(kind):
    out = []
    for f in sorted(glob.glob(f"{CAMP}/{kind}_*_summary.json")):
        m = re.match(rf".*/{kind}_(fifo|other)_(event|poll)_(idle|stress)_r(\d+)_summary.json", f)
        if m: out.append((m.group(1), m.group(2), m.group(3), int(m.group(4)), f[:-len("_summary.json")]))
    return out

def pct(a, q): return float(np.percentile(a, q)) if len(a) else float("nan")
def cliffs_delta(x, y):
    u = stats.mannwhitneyu(x, y, alternative="two-sided").statistic
    return 2.0 * u / (len(x) * len(y)) - 1.0
def cp_upper(k, n, conf=0.95):          # one-sided Clopper-Pearson upper bound on a proportion
    return 1.0 if k >= n else float(stats.beta.ppf(conf, k + 1, n - k))
def fmt(v, d=0): return f"{v:,.{d}f}".replace(",", "\\,")

num = {}            # numbers exported as LaTeX macros
verdicts = []       # run-time verification rows
def verdict(tid, desc, metric, value, op, lim):
    ok = value <= lim if op == "<=" else value >= lim if op == ">=" else value == lim
    verdicts.append((tid, desc, metric, f"{value:.4f}", f"{op} {lim:g}", "PASS" if ok else "FAIL"))

# ------------------------------------------------------------------ A. latency probe
lat, jit, pipe, cacq, cproc, meta, resid = {}, {}, {}, {}, {}, {}, {}
for pol, al, load, rep, base in runs("probe"):
    key = (pol, al, load)
    a = pd.read_csv(base + "_alarms.csv"); a = a[(a.event == "fire") & (a.alarm == "LEADS_OFF")]
    sf = base + "_samples.csv"; s = pd.read_csv(sf if os.path.exists(sf) else sf + ".gz")
    resid.setdefault(key, []).append(a.latency_us.values - s.jitter_ns.values[a.k.values] / 1e3)
    s = s[s.k >= FS]            # skip first second (start-up)
    j = json.load(open(base + "_summary.json"))
    lat.setdefault(key, []).append(a.latency_us.values)
    jit.setdefault(key, []).append(s.jitter_ns.values / 1e3); pipe.setdefault(key, []).append(s.pipeline_ns.values / 1e3)
    cacq.setdefault(key, []).append(s.c_acq_ns.values / 1e3); cproc.setdefault(key, []).append(s.c_proc_ns.values / 1e3)
    m = meta.setdefault(key, {"drops": 0, "wc_alarm": 0.0, "reps": 0, "dur": 0.0})
    m["drops"] += j["sample_drops"] + j["alarm_ev_drops"]; m["wc_alarm"] = max(m["wc_alarm"], j["wcet_alarm_ns"] / 1e3)
    m["reps"] += 1; m["dur"] += j["duration_s"]

rows, tab = [], []
for load in ("idle", "stress"):
    for c in CFG:
        key = c + (load,)
        if key not in lat: continue
        L = np.concatenate(lat[key]); J = np.concatenate(jit[key]); P = np.concatenate(pipe[key])
        miss = int((L > DEADLINE_US).sum()); repmax = [float(x.max()) for x in lat[key]]
        r = dict(load=load, policy=c[0], alarm=c[1], n=len(L), reps=meta[key]["reps"], min=L.min(), p50=pct(L, 50), mean=L.mean(),
                 p99=pct(L, 99), p999=pct(L, 99.9), max=L.max(), miss=miss, miss_pct=100.0 * miss / len(L),
                 miss_ub95_pct=100.0 * cp_upper(miss, len(L)), rep_max_min=min(repmax), rep_max_max=max(repmax),
                 jit_p50=pct(J, 50), jit_p999=pct(J, 99.9), jit_max=J.max(), jit_gt1ms=int((J > 1000).sum()), n_samples=len(J),
                 pipe_p50=pct(P, 50), pipe_p999=pct(P, 99.9), pipe_max=P.max(), pipe_miss=int((P > DEADLINE_US).sum()),
                 drops=meta[key]["drops"], resid_p50=pct(np.concatenate(resid[key]), 50), resid_p999=pct(np.concatenate(resid[key]), 99.9),
                 resid_max=float(np.concatenate(resid[key]).max()),
                 late_share=float(np.mean((np.concatenate(lat[key]) - np.concatenate(resid[key]))[L > 1000] / L[L > 1000])) if (L > 1000).any() else float("nan"),
                 n_gt1ms=int((L > 1000).sum()))
        rows.append(r)
res = pd.DataFrame(rows); res.to_csv("results/latency_summary.csv", index=False)

with open(f"{PAPER}/tables/tab_latency.tex", "w") as f:
    for load in ("idle", "stress"):
        f.write(f"\\multicolumn{{9}}{{l}}{{\\textit{{{'No interference load' if load == 'idle' else 'Overload (8 CPU-bound, 2 memory-churn, 1 file-I/O threads)'}}}}}\\\\\n")
        for _, r in res[res.load == load].iterrows():
            f.write(f"{LBL[(r.policy, r.alarm)]} & {fmt(r.n)} & {fmt(r.p50)} & {fmt(r['mean'])} & {fmt(r.p99)} & {fmt(r.p999)} & {fmt(r['max'])} & "
                    f"{int(r.miss)} ({r.miss_pct:.2f}\\,\\%) & {r.miss_ub95_pct:.3f}\\,\\% \\\\\n")
        if load == "idle": f.write("\\midrule\n")
with open(f"{PAPER}/tables/tab_jitter.tex", "w") as f:
    for load in ("idle", "stress"):
        f.write(f"\\multicolumn{{9}}{{l}}{{\\textit{{{'No interference load' if load == 'idle' else 'Overload'}}}}}\\\\\n")
        for _, r in res[res.load == load].iterrows():
            f.write(f"{LBL[(r.policy, r.alarm)]} & {fmt(r.n_samples)} & {fmt(r.jit_p50, 1)} & {fmt(r.jit_p999)} & {fmt(r.jit_max)} & "
                    f"{fmt(r.pipe_p50, 1)} & {fmt(r.pipe_p999)} & {fmt(r.pipe_max)} & {int(r.pipe_miss)} \\\\\n")
        if load == "idle": f.write("\\midrule\n")

def get(pol, al, load): 
    q = res[(res.policy == pol) & (res.alarm == al) & (res.load == load)]
    return q.iloc[0] if len(q) else None
for pol in ("fifo", "other"):
    for al in ("event", "poll"):
        for load in ("idle", "stress"):
            r = get(pol, al, load)
            if r is None: continue
            tag = {"fifo": "Fp", "other": "Ts"}[pol] + {"event": "Ev", "poll": "Po"}[al] + {"idle": "Idle", "stress": "Str"}[load]
            num[tag + "N"] = fmt(r.n); num[tag + "Med"] = fmt(r.p50); num[tag + "Max"] = fmt(r["max"]); num[tag + "MaxMs"] = f"{r['max'] / 1e3:.2f}"
            num[tag + "PNNN"] = fmt(r.p999); num[tag + "Miss"] = str(int(r.miss)); num[tag + "MissPct"] = f"{r.miss_pct:.2f}"
            num[tag + "UB"] = f"{r.miss_ub95_pct:.3f}"; num[tag + "JitMax"] = fmt(r.jit_max); num[tag + "JitPNNN"] = fmt(r.jit_p999)
            num[tag + "ResMax"] = fmt(r.resid_max); num[tag + "ResPNNN"] = fmt(r.resid_p999); num[tag + "NGtOne"] = str(int(r.n_gt1ms))
            num[tag + "LateShare"] = "n/a" if r.late_share != r.late_share else f"{100 * r.late_share:.0f}"
            num[tag + "PipeMiss"] = str(int(r.pipe_miss)); num[tag + "PipeMissPpm"] = f"{1e6 * r.pipe_miss / r.n_samples:.0f}"
            num[tag + "PipeMaxMs"] = f"{r.pipe_max / 1e3:.1f}"; num[tag + "JitMaxMs"] = f"{r.jit_max / 1e3:.1f}"; num[tag + "ResP"] = fmt(r.resid_p50)
            num[tag + "PipeMax"] = fmt(r.pipe_max); num[tag + "JitOver"] = str(int(r.jit_gt1ms)); num[tag + "NSamp"] = fmt(r.n_samples)

# hypothesis tests under overload
tests = {}
def cmp(a, b, name):
    if a in lat and b in lat:
        x, y = np.concatenate(lat[a]), np.concatenate(lat[b])
        p = stats.mannwhitneyu(x, y, alternative="two-sided").pvalue
        tests[name] = dict(p=p, cliffs_delta=cliffs_delta(x, y), n1=len(x), n2=len(y))
        num[name + "Delta"] = f"{tests[name]['cliffs_delta']:.3f}"; num[name + "P"] = "<10^{-300}" if p < 1e-300 else f"{p:.1e}"
cmp(("other", "event", "stress"), ("fifo", "event", "stress"), "TestPolicy")
cmp(("fifo", "poll", "stress"), ("fifo", "event", "stress"), "TestPath")
cmp(("fifo", "event", "stress"), ("fifo", "event", "idle"), "TestLoad")
cmp(("other", "event", "stress"), ("other", "event", "idle"), "TestLoadTs")

r = get("fifo", "event", "stress")
if r is not None:
    verdict("TC-RT-01", f"Alarm latency, proposed configuration under overload ({int(r.n)} alarms, {int(r.reps)} runs): maximum", "L_max_ms", r["max"] / 1e3, "<=", 10.0)
    verdict("TC-RT-01b", f"Pipeline latency of every sample (release to end of processing), proposed configuration under overload ({int(r.n_samples)} samples): samples above 10 ms", "late_samples", float(r.pipe_miss), "==", 0)
    verdict("TC-RT-02", f"Acquisition release jitter, proposed configuration under overload ({int(r.n_samples)} samples): maximum", "J_max_ms", r.jit_max / 1e3, "<=", 1.0)
    verdict("TC-RT-03", "Samples or alarm events lost, proposed configuration, all probe runs", "lost",
            float(sum(meta[k]["drops"] for k in meta if k[0] == "fifo" and k[1] == "event")), "==", 0)

# WCET observations and response-time bound (proposed configuration, overload)
key = ("fifo", "event", "stress")
if key in cacq:
    Ca, Cp, Cal = np.concatenate(cacq[key]).max(), np.concatenate(cproc[key]).max(), meta[key]["wc_alarm"]
    supp = {}
    for f in glob.glob("results/supp/*_summary.json"):
        j = json.load(open(f))
        for k in ("wcet_wdg_ns", "wcet_ui_cs_ns"): supp[k] = max(supp.get(k, 0), j.get(k, 0))
    Cw, B = supp.get("wcet_wdg_ns", 0) / 1e3, supp.get("wcet_ui_cs_ns", 0) / 1e3
    Rb = B + Ca + Cp + Cal; R = Rb                      # fixed-point iteration with watchdog interference (T = 20 ms)
    for _ in range(50): R = Rb + math.ceil(R / 20000.0) * Cw
    Lmax = r["max"]
    wc = dict(C_acq=Ca, C_proc=Cp, C_alarm=Cal, C_wdg=Cw, B_ui=B, R_app=R, L_max=Lmax, platform=Lmax - R,
              C_acq_mean=float(np.concatenate(cacq[key]).mean()), C_proc_mean=float(np.concatenate(cproc[key]).mean()),
              C_acq_pq=pct(np.concatenate(cacq[key]), 99.99), C_proc_pq=pct(np.concatenate(cproc[key]), 99.99),
              resid_max=float(r.resid_max), jit_max=float(r.jit_max))
    json.dump(wc, open("results/wcet_rta.json", "w"), indent=1)
    for k, v in wc.items(): num["Wc" + k.replace("_", "")] = fmt(v, 1)
    num["WcUtil"] = f"{100.0 * (wc['C_acq_mean'] + wc['C_proc_mean']) / 2000.0:.1f}"
    num["WcShare"] = f"{100.0 * R / Lmax:.0f}"

# figure: complementary CDF of alarm latency
fig, axs = plt.subplots(1, 2, figsize=(7.0, 2.9), sharey=True)
for ax, load, ttl in zip(axs, ("idle", "stress"), ("(a) no interference load", "(b) overload")):
    for c in CFG:
        key = c + (load,)
        if key not in lat: continue
        x = np.sort(np.concatenate(lat[key])) / 1e3; y = 1.0 - np.arange(len(x)) / len(x)
        ax.step(x, y, where="post", color=COL[c], ls=LS[c], label=LBL[c])
    ax.axvline(10, color="#444", lw=0.9); ax.text(10.6, 0.5, "10 ms deadline", rotation=90, va="center", fontsize=7.5, color="#444")
    ax.set_xscale("log"); ax.set_yscale("log"); ax.set_xlabel("alarm latency $L$ (ms)"); ax.set_title(ttl, fontsize=9, loc="left")
    ax.set_xlim(0.01, 300); ax.set_ylim(5e-5, 1.3)
axs[0].set_ylabel("fraction of alarms with latency > $L$"); axs[1].legend(loc="lower left", fontsize=7.5)
fig.tight_layout(); fig.savefig(f"{PAPER}/figs/fig_latency_ccdf.pdf"); fig.savefig(f"{PAPER}/figs/fig_latency_ccdf.png", dpi=200); plt.close(fig)

# figure: release jitter CCDF under overload + per-run maxima
fig, axs = plt.subplots(1, 2, figsize=(7.0, 2.8))
for c in (("fifo", "event"), ("other", "event")):
    key = c + ("stress",)
    if key in jit:
        x = np.sort(np.concatenate(jit[key])) / 1e3; x = np.maximum(x, 1e-3); y = 1.0 - np.arange(len(x)) / len(x)
        axs[0].step(x, y, where="post", color=COL[c], ls=LS[c], label={"fifo": "FP (SCHED_FIFO)", "other": "TS (time-sharing)"}[c[0]])
axs[0].axvline(1, color="#444", lw=0.9); axs[0].text(1.07, 0.4, "1 ms limit", rotation=90, va="center", fontsize=7.5, color="#444")
axs[0].set_xscale("log"); axs[0].set_yscale("log"); axs[0].set_xlabel("acquisition release jitter (ms)"); axs[0].set_ylabel("fraction of samples above")
axs[0].set_title("(a) release jitter under overload", fontsize=9, loc="left"); axs[0].legend(fontsize=7.5, loc="lower left"); axs[0].set_ylim(5e-6, 1.3)
xs, labs = [], []
for i, c in enumerate(CFG):
    for jx, load in enumerate(("idle", "stress")):
        key = c + (load,)
        if key not in lat: continue
        mx = [x.max() / 1e3 for x in lat[key]]
        axs[1].scatter([i + (jx - 0.5) * 0.3] * len(mx), mx, s=26, color=COL[c], marker="o" if load == "idle" else "s",
                       facecolors="none" if load == "idle" else COL[c], zorder=3)
axs[1].scatter([], [], s=26, color="#555", facecolors="none", marker="o", label="no load"); axs[1].scatter([], [], s=26, color="#555", marker="s", label="overload")
axs[1].axhline(10, color="#444", lw=0.9); axs[1].set_yscale("log"); axs[1].set_xticks(range(4))
axs[1].set_xticklabels(["FP\nevent", "FP\npolling", "TS\nevent", "TS\npolling"]); axs[1].set_ylabel("per-run maximum latency (ms)")
axs[1].set_title("(b) worst alarm latency in each run", fontsize=9, loc="left"); axs[1].legend(fontsize=7.5, loc="upper left")
fig.tight_layout(); fig.savefig(f"{PAPER}/figs/fig_jitter.pdf"); fig.savefig(f"{PAPER}/figs/fig_jitter.png", dpi=200); plt.close(fig)

# ------------------------------------------------------------------ B. watchdog fault injection
st = {}
for pol, al, load, rep, base in runs("stall"):
    a = pd.read_csv(base + "_alarms.csv"); a = a[(a.event == "fire") & (a.alarm == "SYS_STALL")]
    st.setdefault(pol, []).append(a.latency_us.values[0] / 1e3 if len(a) else float("nan"))
with open(f"{PAPER}/tables/tab_stall.tex", "w") as f:
    for pol in ("fifo", "other"):
        if pol in st:
            v = np.array(st[pol]); nd = int(np.isnan(v).sum())
            f.write(f"{'FP + event (proposed)' if pol == 'fifo' else 'TS + event'} & {len(v)} & {nd} & {np.nanmin(v):.1f} & {np.nanmean(v):.1f} & {np.nanmax(v):.1f} \\\\\n")
            num[("Fp" if pol == "fifo" else "Ts") + "StallMax"] = f"{np.nanmax(v):.1f}"; num[("Fp" if pol == "fifo" else "Ts") + "StallMin"] = f"{np.nanmin(v):.1f}"
if "fifo" in st:
    v = np.array(st["fifo"])
    verdict("TC-WDG-01", f"Runaway processing thread injected under overload ({len(v)} runs): undetected stalls", "undetected", float(np.isnan(v).sum()), "==", 0)
    verdict("TC-WDG-01b", "Runaway processing thread: maximum time from stall start to SYS_STALL annunciation", "t_max_ms", float(np.nanmax(v)), "<=", 200.0)

# ------------------------------------------------------------------ C. clinical scenario in real time
EXPECT = [("ASYSTOLE", 1, "last"), ("LOW_PERF", 1, "start"), ("TACHY", 3, "start"), ("BRADY", 5, "start"), ("DESAT", 7, "start"),
          ("SPO2_LOW", 7, "start"), ("LEADS_OFF", 9, "start")]
LIMIT = {"ASYSTOLE": 4.5, "LOW_PERF": 10, "TACHY": 10, "BRADY": 10, "DESAT": 10, "SPO2_LOW": 15, "LEADS_OFF": 0.504}
clin, clin_err, clin_lat = {}, 0, {}
for pol, al, load, rep, base in runs("clinical"):
    a = pd.read_csv(base + "_alarms.csv"); t = pd.read_csv(base + "_truth.csv"); fires = a[a.event == "fire"]
    tag = (pol, load)
    if pol == "fifo" and (len(fires) != len(EXPECT) or set(fires.alarm) != {e[0] for e in EXPECT}): clin_err += 1
    for name, seg, ref in EXPECT:
        q = fires[fires.alarm == name]
        if not len(q): continue
        k0 = t.k_last_r[seg] if ref == "last" else t.k_start[seg]
        clin.setdefault(name, {}).setdefault(tag, []).append((q.k.values[0] - k0) / FS)
        clin_lat.setdefault(tag, []).append(q.latency_us.values[0])
with open(f"{PAPER}/tables/tab_clinical.tex", "w") as f:
    worst_margin = 1e9
    for name, seg, ref in EXPECT:
        if name not in clin: continue
        allfp = [v for tg, vs in clin[name].items() if tg[0] == "fifo" for v in vs]
        ts = clin[name].get(("other", "stress"), [float("nan")])
        pr = {"ASYSTOLE": "High", "TACHY": "High", "BRADY": "High", "DESAT": "High", "SPO2_LOW": "Medium", "LEADS_OFF": "Medium", "LOW_PERF": "Low"}[name]
        refs = {"last": "last R wave", "start": "onset"}[ref]
        f.write(f"\\texttt{{{name.replace('_', chr(92) + '_')}}} & {pr} & {refs} & {min(allfp):.3f}--{max(allfp):.3f} & {ts[0]:.3f} & $\\le${LIMIT[name]:g} \\\\\n")
        worst_margin = min(worst_margin, LIMIT[name] - max(allfp))
    fp_l = [v for tg, vs in clin_lat.items() if tg[0] == "fifo" for v in vs]; ts_l = [v for tg, vs in clin_lat.items() if tg[0] == "other" for v in vs]
    if fp_l: num["ClinFpLatMax"] = fmt(max(fp_l)); num["ClinFpN"] = str(len(fp_l))
    if ts_l: num["ClinTsLatMax"] = fmt(max(ts_l)); num["ClinTsLatMaxMs"] = f"{max(ts_l) / 1e3:.2f}"
if clin:
    nfp = len([1 for p, a_, l, r_, b in runs("clinical") if p == "fifo"])
    verdict("TC-RT-04", f"Clinical scenario in real time, proposed configuration ({nfp} runs): runs with missing, extra or repeated alarms", "bad_runs", float(clin_err), "==", 0)
    verdict("TC-RT-04b", "Clinical scenario in real time: smallest margin to the time-to-alarm limits", "margin_s", worst_margin, ">=", 0.0)
    if fp_l: verdict("TC-RT-04c", "Clinical scenario in real time: maximum alarm-system latency", "L_max_ms", max(fp_l) / 1e3, "<=", 10.0)

# figure: clinical scenario overview from the off-target trace (same code, same seed as run r1 idle)
if os.path.exists("results/clinical_trace.csv"):
    d = pd.read_csv("results/clinical_trace.csv"); tt = d.k / FS
    fig, axs = plt.subplots(3, 1, figsize=(7.0, 4.6), gridspec_kw={"height_ratios": [1.15, 1, 1]})
    w = (tt > 17.5) & (tt < 26.0); axs[0].plot(tt[w], d.ecg[w], color="#2a78d6", lw=0.8)
    ka = d.k[d.asys == 1].min() / FS; kl = 20.244
    q = d[(d.qrs == 1) & w]; 
    axs[0].axvline(ka, color="#d03b3b", lw=1.2); axs[0].text(ka + 0.08, 1.05, "ASYSTOLE alarm", color="#333", fontsize=7.5)
    lastq = ka - 4.0
    axs[0].annotate("", xy=(ka, -0.32), xytext=(lastq, -0.32), arrowprops=dict(arrowstyle="<->", color="#555", lw=0.8))
    axs[0].text((ka + lastq) / 2, -0.27, "4.0 s since last detected QRS", ha="center", fontsize=7.5, color="#333")
    axs[0].set_ylabel("ECG (mV)"); axs[0].set_xlabel("time (s)"); axs[0].set_title("(a) asystole onset", fontsize=9, loc="left"); axs[0].set_ylim(-0.45, 1.3)
    hr = d.hr.where(d.hr_ok == 1); axs[1].plot(tt[::25], hr[::25], color="#2a78d6", lw=1.2)
    axs[1].axhline(150, color="#777", lw=0.7, ls="--"); axs[1].axhline(40, color="#777", lw=0.7, ls="--"); axs[1].set_ylabel("heart rate (bpm)")
    axs[1].set_title("(b) heart rate, limits 40/150 bpm; shaded: alarm condition true", fontsize=9, loc="left")
    sp = d.spo2.where((d.spo2_ok == 1) & (d.ir >= 0)); m = d.ir >= 0
    axs[2].plot(tt[m], sp[m], color="#2a78d6", lw=1.2); axs[2].axhline(90, color="#777", lw=0.7, ls="--"); axs[2].axhline(85, color="#777", lw=0.7, ls="--")
    axs[2].set_ylabel("SpO$_2$ (%)"); axs[2].set_xlabel("time (s)"); axs[2].set_ylim(76, 100)
    axs[2].set_title("(c) SpO$_2$, limits 90/85 %", fontsize=9, loc="left")
    def shade(ax, col, name, y):
        c = d[col].values; on = np.flatnonzero(np.diff(np.r_[0, c, 0]) == 1); off = np.flatnonzero(np.diff(np.r_[0, c, 0]) == -1)
        for a0, b0 in zip(on, off):
            ax.axvspan(a0 / FS, b0 / FS, color="#e34948", alpha=0.13, lw=0); ax.text(a0 / FS + 0.5, y, name, fontsize=7, color="#333", va="top")
    shade(axs[1], "asys", "ASYSTOLE", 195); shade(axs[1], "tachy", "TACHY", 120); shade(axs[1], "brady", "BRADY", 195); shade(axs[1], "lo", "LEADS_OFF", 195)
    shade(axs[2], "low", "SPO2_LOW", 93.5); shade(axs[2], "perf", "LOW_PERF", 93.5)
    c = d["desat"].values; on = np.flatnonzero(np.diff(np.r_[0, c]) == 1)
    for a0 in on: axs[2].axvline(a0 / FS, color="#d03b3b", lw=1.0); axs[2].text(a0 / FS - 0.8, 93.5, "DESAT", fontsize=7, color="#333", ha="right", va="top")
    axs[1].set_ylim(0, 200); axs[1].set_xlim(0, 210); axs[2].set_xlim(0, 210)
    fig.tight_layout(); fig.savefig(f"{PAPER}/figs/fig_clinical.pdf"); fig.savefig(f"{PAPER}/figs/fig_clinical.png", dpi=200); plt.close(fig)

# ------------------------------------------------------------------ outputs
pd.DataFrame(verdicts, columns=["test_id", "description", "metric", "value", "limit", "verdict"]).to_csv("results/verification_rt.csv", index=False)
json.dump({"tests": tests, "stall_ms": st, "clinical_s": {k: {f"{a}/{b}": v for (a, b), v in d_.items()} for k, d_ in clin.items()}}, open("results/stats.json", "w"), indent=1, default=float)
with open(f"{PAPER}/tables/numbers.tex", "w") as f:
    for k, v in sorted(num.items()): f.write(f"\\newcommand{{\\N{k}}}{{{v}}}\n")
print(res[["load", "policy", "alarm", "n", "p50", "p999", "max", "miss", "jit_p999", "jit_max", "jit_gt1ms", "pipe_miss", "resid_p50", "resid_p999", "resid_max", "n_gt1ms", "late_share", "drops"]].round(1).to_string(index=False))
for v in verdicts: print(v)
print(json.dumps(tests, indent=1, default=float))
