#!/usr/bin/env python3
"""Build the hazard log and traceability matrix from the source tables and the *executed* verification results.
Inputs : docs/hazard_log.csv, docs/requirements.csv, results/verification.csv, results/verification_rt.csv, src/monitor.c
Outputs: docs/HAZARD_LOG.md, docs/traceability.csv, paper/tables/tab_hazard.tex, paper/tables/tab_trace.tex"""
import csv, os, re, datetime
import pandas as pd

def risk_class(s, p):
    r = s * p
    return "Unacceptable" if r >= 12 else "Tolerable" if r >= 6 else "Acceptable"

# ---- TC-INSP-01: automated code inspection of the real-time thread bodies
src = open("src/monitor.c").read()
def body(fn):
    m = re.search(r"\n(?:static )?[\w \*]+\b" + fn + r"\([^)]*\)\s*\{", src); i = m.end(); depth = 1
    while depth: depth += {"{": 1, "}": -1}.get(src[i], 0); i += 1
    return src[m.end():i]
banned = re.compile(r"\b(malloc|calloc|realloc|free|fopen|fclose|fprintf|printf|fputs|fwrite|usleep|sleep)\s*\(")
viol = sum(len(banned.findall(body(f))) for f in ("th_acq", "th_proc", "th_alarm", "th_wdg", "post_alarm", "alarm_handle", "log_alarm"))
pi_ok = all(re.search(p, src) for p in (r"chan_init\(&ch_samp,[^;]*cfg\.fifo\)", r"chan_init\(&ch_alarm,[^;]*cfg\.fifo\)",
                                         r"rt_mutex_init\(&pt\.mu, cfg\.fifo\)", r"rt_mutex_init\(&snap\.mu, cfg\.fifo\)"))
mlock_ok = "rt_lock_memory()" in src and "PTHREAD_PRIO_INHERIT" in open("src/rt_util.c").read()
insp = [("TC-INSP-01", "Automated inspection: heap, stdio or sleep calls inside real-time thread bodies; PI on all shared mutexes; memory locking",
         "violations", f"{viol + (0 if pi_ok else 1) + (0 if mlock_ok else 1):.4f}", "== 0", "PASS" if viol == 0 and pi_ok and mlock_ok else "FAIL")]

ver = pd.concat([pd.read_csv("results/verification.csv"), pd.read_csv("results/verification_rt.csv"),
                 pd.DataFrame(insp, columns=["test_id", "description", "metric", "value", "limit", "verdict"])], ignore_index=True)
ver.to_csv("results/verification_all.csv", index=False)
V = {r.test_id: r for r in ver.itertuples()}
def tests_for(prefixes):                      # 'TC-ASY-01' also covers its sub-results TC-ASY-01b, TC-ASY-01c
    out = []
    for p in prefixes:
        out += [t for t in V if t == p or (t.startswith(p) and len(t) == len(p) + 1 and t[-1].isalpha())]
    return sorted(set(out))
def status(tests):
    if not tests: return "NOT RUN"
    if any(t not in V for t in tests): return "NOT RUN"
    return "PASS" if all(V[t].verdict == "PASS" for t in tests) else "FAIL"

req = pd.read_csv("docs/requirements.csv"); haz = pd.read_csv("docs/hazard_log.csv").fillna("")
R = {}
for r in req.itertuples():
    want = r.verification.split("|"); ts = tests_for(want)
    missing = [w for w in want if w not in V]            # a planned test without an executed result
    R[r.id] = dict(text=r.requirement, tests=ts, status="NOT RUN" if missing else status(ts))

rows, trace = [], []
for h in haz.itertuples():
    reqs = h.requirements.split("|"); ts = sorted({t for q in reqs for t in R[q]["tests"]})
    st = "NOT RUN" if any(R[q]["status"] == "NOT RUN" for q in reqs) else status(ts); ri, rr = h.sev * h.p_init, h.sev * h.p_res
    ver_state = "Verified" if st == "PASS" else "OPEN: control not verified" if st == "FAIL" else "Not run"
    rows.append(dict(h=h, reqs=reqs, tests=ts, st=st, ri=ri, rr=rr, ci=risk_class(h.sev, h.p_init), cr=risk_class(h.sev, h.p_res), vs=ver_state))
    for q in reqs:
        for t in R[q]["tests"]: trace.append((h.id, h.sw_item, q, t, V[t].verdict if t in V else "NOT RUN"))
pd.DataFrame(trace, columns=["hazard", "software_item", "requirement", "test", "verdict"]).to_csv("docs/traceability.csv", index=False)

n_pass, n_fail = int((ver.verdict == "PASS").sum()), int((ver.verdict == "FAIL").sum())
with open("docs/HAZARD_LOG.md", "w") as f:
    f.write(f"""# Hazard log and traceability report

*Generated {datetime.date.today().isoformat()} by `scripts/gen_hazard_log.py` from executed verification results. Do not edit by hand.*

> **Scope and status.** This document is *styled after* the risk-management file required by ISO 14971 and the
> traceability chain of IEC 62304 (hazardous situation -> software item -> software cause -> risk control measure ->
> verification). It belongs to a research/teaching **simulator**, not a medical device. It has not been reviewed by a
> clinical or regulatory expert, severity/probability ratings are illustrative engineering judgements with no field
> data behind them, and nothing here constitutes evidence of conformity to any standard.

## 1. Risk acceptability scheme

Severity S: 1 negligible, 2 minor, 3 serious, 4 critical, 5 catastrophic. Probability of harm P: 1 improbable ... 5 frequent.
Following IEC 62304, the probability that a software fault *occurs* is taken as 1; P rates the likelihood that the
resulting hazardous situation leads to harm, before and after the risk control measures.
Risk index R = S x P: **>= 12 unacceptable**, **6-10 tolerable** (reduce as far as possible), **<= 5 acceptable**.

Software safety classification (IEC 62304, illustrative): the alarm path (acquisition, processing, alarm manager,
watchdog) is treated as **Class C**; the display/operator thread as Class B, segregated by scheduling class and priority.

## 2. Hazard log

Verification summary: **{n_pass} PASS, {n_fail} FAIL** out of {len(ver)} executed checks.

""")
    for r in rows:
        h = r["h"]
        f.write(f"### {h.id} - {h.hazard}\n\n| | |\n|---|---|\n| Software cause | {h.cause} |\n| Hazardous situation | {h.hazardous_situation} |\n"
                f"| Harm | {h.harm} |\n| Software item | `{h.sw_item}` |\n| Initial risk | S{h.sev} x P{h.p_init} = {r['ri']} ({r['ci']}) |\n"
                f"| Risk control measures | {h.controls} |\n| Control type | {h.control_type} |\n| Requirements | {', '.join(r['reqs'])} |\n"
                f"| Verification | {', '.join(f'{t} ({V[t].verdict})' for t in r['tests'])} |\n"
                f"| Residual risk | S{h.sev} x P{h.p_res} = {r['rr']} ({r['cr']}) |\n| Status | **{r['vs']}** |\n"
                + (f"| Residual-risk note | {h.residual_note} |\n" if h.residual_note else "") + "\n")
    f.write("## 3. Software requirements and verification status\n\n| ID | Requirement | Tests | Status |\n|---|---|---|---|\n")
    for q, d in R.items(): f.write(f"| {q} | {d['text']} | {', '.join(d['tests'])} | **{d['status']}** |\n")
    f.write("\n## 4. Verification results (as executed)\n\n| Test | Description | Metric | Value | Limit | Verdict |\n|---|---|---|---|---|---|\n")
    for v in ver.itertuples(): f.write(f"| {v.test_id} | {v.description} | {v.metric} | {float(v.value):g} | {v.limit} | **{v.verdict}** |\n")
    f.write("\n## 5. Traceability matrix\n\nSee `docs/traceability.csv` (hazard, software item, requirement, test, verdict).\n")

os.makedirs("paper/tables", exist_ok=True)
tex = lambda s: str(s).replace("_", "\\_").replace("%", "\\%").replace("&", "\\&").replace(">=", "$\\ge$").replace("<=", "$\\le$").replace("+/-", "$\\pm$")
def short(ts):                               # collapse TC-ASY-01, -01b, -01c -> TC-ASY-01
    return sorted({re.sub(r"(?<=\d)[a-z]$", "", t) for t in ts})
with open("paper/tables/tab_hazard.tex", "w") as f:
    for r in rows:
        h = r["h"]
        f.write(f"{h.id} & {tex(h.hazard)}: {tex(h.cause[0].lower() + h.cause[1:])} & {r['ri']} & {tex(h.controls)} & "
                f"{tex(', '.join(q.replace('SRS-', '') for q in r['reqs']))} & {tex(', '.join(t.replace('TC-', '') for t in short(r['tests'])))} & {r['rr']} & "
                f"{'\\checkmark' if r['st'] == 'PASS' else '\\textbf{open}' if r['st'] == 'FAIL' else 'not run'} \\\\\n")
with open("paper/tables/tab_trace.tex", "w") as f:
    f.write(f"\\newcommand{{\\NVerTotal}}{{{len(ver)}}}\\newcommand{{\\NVerPass}}{{{n_pass}}}\\newcommand{{\\NVerFail}}{{{n_fail}}}\n"
            f"\\newcommand{{\\NHaz}}{{{len(rows)}}}\\newcommand{{\\NReq}}{{{len(R)}}}\\newcommand{{\\NTraceLinks}}{{{len(trace)}}}\n"
            f"\\newcommand{{\\NHazOpen}}{{{sum(1 for r in rows if r['st'] == 'FAIL')}}}\n"
            f"\\newcommand{{\\NHazUnaccInit}}{{{sum(1 for r in rows if r['ci'] == 'Unacceptable')}}}\n"
            f"\\newcommand{{\\NHazUnaccRes}}{{{sum(1 for r in rows if r['cr'] == 'Unacceptable')}}}\n"
            f"\\newcommand{{\\NFailList}}{{{tex(', '.join(ver[ver.verdict == 'FAIL'].test_id)) or 'none'}}}\n")
print(f"hazards={len(rows)} requirements={len(R)} checks={len(ver)} pass={n_pass} fail={n_fail}")
for r in rows: print(r["h"].id, r["ri"], r["ci"], "->", r["rr"], r["cr"], r["vs"])
print("unlinked tests:", sorted(set(V) - {t for d in R.values() for t in d["tests"]}))
