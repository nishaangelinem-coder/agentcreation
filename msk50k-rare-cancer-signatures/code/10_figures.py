#!/usr/bin/env python
"""Manuscript figures (300 dpi PNG)."""
import json, numpy as np, pandas as pd, matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch
from lifelines import KaplanMeierFitter
from common import *
PAL = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#6250d6", "#e34948"]
FC = dict(zip(FAMILIES, PAL)); FC.update({"Other": "#898781", "AMBIGUOUS": "#c3c2b7", "Other called": "#898781"})
INK, MUT, GRID = "#0b0b0b", "#52514e", "#e1e0d9"
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 8, "axes.edgecolor": "#c3c2b7", "axes.labelcolor": INK, "xtick.color": MUT, "ytick.color": MUT,
                     "axes.spines.top": False, "axes.spines.right": False, "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.5, "axes.axisbelow": True,
                     "legend.frameon": False, "figure.dpi": 100, "savefig.dpi": 300, "savefig.bbox": "tight", "axes.titlesize": 9, "axes.titleweight": "bold", "axes.titlelocation": "left"})
J = lambda f: json.load(open(f"{RES}/{f}"))
flow, cons, st, ext = J("cohort_flow.json")["flow"], J("consensus_summary.json"), J("stats_summary.json"), J("external_validation.json")
b = pd.read_csv(f"{RES}/basket_sample_calls.tsv", sep="\t")
def tag(ax, s): ax.text(-0.12, 1.08, s, transform=ax.transAxes, fontsize=11, fontweight="bold", va="bottom", color=INK)
# ---------------- Fig 1 flow
fig, ax = plt.subplots(figsize=(9.5, 7.0)); ax.axis("off"); ax.set_xlim(0, 10); ax.set_ylim(0, 10)
def box(x, y, w, h, t, fc="#ffffff"):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.02,rounding_size=0.12", fc=fc, ec="#52514e", lw=0.8)); ax.text(x + w / 2, y + h / 2, t, ha="center", va="center", fontsize=7, color=INK)
def arrow(x1, y1, x2, y2): ax.annotate("", (x2, y2), (x1, y1), arrowprops=dict(arrowstyle="-|>", color="#52514e", lw=0.8))
box(2.5, 9.0, 5, 0.8, f"MSK-IMPACT 50K study (msk_impact_50k_2026)\n{flow['all_samples']:,} tumours, {flow['all_patients']:,} patients; all tumour–normal matched")
arrow(5, 9.0, 5, 8.45); box(2.0, 7.55, 6, 0.9, "Basket selection by OncoTree code + detailed cancer type\n" + "\n".join(["; ".join(f"{g} {flow['basket_samples'][g]} samples / {flow['basket_patients'][g]} patients" for g in ["ASI", "THYM"]), f"ACC {flow['basket_samples']['ACC']} samples / {flow['basket_patients']['ACC']} patients"]))
box(8.3, 7.7, 1.65, 0.6, f"Excluded: {flow['all_samples'] - sum(flow['basket_samples'].values()):,}\nother cancer types", "#f4f3ef"); arrow(8.0, 8.0, 8.3, 8.0)
xs = {"ASI": 0.2, "THYM": 3.5, "ACC": 6.8}; names = {"ASI": "Angiosarcoma (ASI)", "THYM": "Thymic epithelial tumour (THYM)", "ACC": "Adrenocortical carcinoma (ACC)"}
for g, x in xs.items():
    fc_ = flow["fit_class"][g]; arrow(5, 7.55, x + 1.5, 6.95)
    box(x, 6.15, 3, 0.8, f"{names[g]}\nOne tumour per patient: n = {flow['after_dedup'][g]}\n(duplicates removed: {flow['basket_samples'][g] - flow['after_dedup'][g]})")
    arrow(x + 1.5, 6.15, x + 1.5, 5.6); box(x, 4.4, 3, 1.2, f"SNV-count fit class\nFULL (≥4 SNVs): {fc_['FULL']}\nSPARSE (2–3 SNVs): {fc_['SPARSE']}\nDROP (<2 SNVs): {fc_['DROP']}")
    ok = flow["gate_ge50_per_group"][g]; arrow(x + 1.5, 4.4, x + 1.5, 3.85)
    box(x, 2.85, 3, 1.0, f"Analysable per-sample: n = {flow['analysable'][g]}\nGate (≥50): {'PASS' if ok else 'FAIL'}\nHypermutators (>50 SNVs): {flow['hypermutators'].get(g, 0)}", "#e6f4ec" if ok else "#fbe9e4")
    arrow(x + 1.5, 2.85, x + 1.5, 2.3); n_os = int(b[(b.GROUP == g)].OS_MONTHS.notna().sum())
    box(x, 1.3, 3, 1.0, f"Pooled-spectrum analysis: n = {flow['after_dedup'][g]}\n({cons['pooled'][g]['n_snv']} SNVs)\nSurvival analysis: n = {n_os}", "#eaf2fc")
fig.savefig(f"{FIG}/Fig1_flow.png"); plt.close()
# ---------------- Fig 2 panel bias + simulation gate
P = pd.read_csv(f"{RES}/trinucleotide_opportunity_proportions.tsv", sep="\t", index_col=0); R = pd.read_csv(f"{RES}/simulation_gate.tsv", sep="\t")
fig, axs = plt.subplots(1, 3, figsize=(11, 3.2), gridspec_kw={"width_ratios": [1.5, 1, 1], "wspace": 0.32})
ax = axs[0]; r = (P.IMPACT505 / P.genome).sort_values(); ax.bar(range(32), np.log2(r.values), color=np.where(r.values > 1, PAL[1], PAL[0]), width=0.72)
ax.set_xticks(range(32)); ax.set_xticklabels(r.index, rotation=90, fontsize=6, family="monospace"); ax.set_ylabel("log2 opportunity ratio\n(panel territory / genome)"); ax.axhline(0, color="#52514e", lw=0.8)
ax.set_title("Trinucleotide opportunity, IMPACT505 vs genome"); ax.grid(axis="x", visible=False); tag(ax, "a")
ax = axs[1]; A = R[R.scenario == "ALL"]
for (name, lab), c in zip([("uncorrected_genome", "No correction"), ("exome_SPA", "Exome frame"), ("panel_opportunity", "Panel frame")], [PAL[7], PAL[0], PAL[2]]):
    d = A[A.correction == name]; ax.plot(d.N, d.cos_family, "-o", color=c, lw=2, ms=4, label=lab)
ax.axhline(0.90, color="#52514e", lw=0.8, ls="--"); ax.text(2.6, 0.905, "gate 0.90", fontsize=7, color=MUT, va="bottom")
ax.axvspan(2, 4, color="#eda100", alpha=0.18, lw=0); ax.text(3, 0.42, "cohort\nmedian\n1–2 SNVs", fontsize=6.5, ha="center", color=MUT)
ax.set_xscale("log"); ax.set_xlabel("SNVs per synthetic tumour"); ax.set_ylabel("Mean cosine, true vs inferred\nfamily exposures"); ax.set_ylim(0.4, 1.01); ax.legend(loc="lower right", fontsize=7); ax.set_title("Simulation gate (200 tumours per N)"); tag(ax, "b")
ax = axs[2]; d = R[(R.correction == "panel_opportunity") & (R.scenario != "ALL")]
for sc, c in [("UV", FC["UV"]), ("APOBEC", FC["APOBEC"]), ("Clock", FC["Clock"]), ("HRD", FC["HRD"])]:
    x = d[d.scenario == sc]; ax.plot(x.N, x.nnls_dominant_acc, "-o", color=c, lw=2, ms=4, label=sc)
ax.set_xscale("log"); ax.set_xlabel("SNVs per synthetic tumour"); ax.set_ylabel("Dominant-process accuracy"); ax.set_ylim(0, 1.03); ax.legend(loc="lower right", fontsize=7, title="Simulated process", title_fontsize=7); ax.set_title("Dominant-process recovery"); tag(ax, "c")
fig.savefig(f"{FIG}/Fig2_panel_bias_simulation.png"); plt.close()
# ---------------- Fig 3 pooled spectra
S = pd.read_csv(f"{RES}/pooled_spectra.tsv", sep="\t", index_col=0); SC = dict(zip(SUBS, ["#2a78d6", "#0b0b0b", "#e34948", "#898781", "#1baf7a", "#e87ba4"]))
keys = ["ASI: Cutaneous/H&N", "ASI: Breast", "ASI: Visceral/soft tissue/other", "ACC", "THYM"]
fig, axs = plt.subplots(len(keys), 1, figsize=(9, 8.2), sharex=True, gridspec_kw={"hspace": 0.42})
for ax, k in zip(axs, keys):
    v = S[k] / S[k].sum(); ax.bar(range(96), v.values * 100, color=[SC[c[2:5]] for c in CHANNELS], width=0.8); ax.grid(axis="x", visible=False); ax.set_ylabel("% of SNVs")
    p = cons["pooled"][k]; ax.set_title(f"{k.replace('ASI: ', 'Angiosarcoma – ').replace('ACC', 'Adrenocortical carcinoma').replace('THYM', 'Thymic epithelial tumours')}   (n = {p['n_patients']} patients, {p['n_snv']} SNVs)"); ax.set_xlim(-1, 96)
    for i in range(1, 6): ax.axvline(i * 16 - 0.5, color=GRID, lw=0.8)
for i, s in enumerate(SUBS): axs[0].text(i * 16 + 7.5, axs[0].get_ylim()[1] * 1.0, s, ha="center", va="top", fontsize=8, fontweight="bold", color=INK); axs[0].plot([i * 16, i * 16 + 15], [axs[0].get_ylim()[1] * 0.86] * 2, color=SC[s], lw=3, solid_capstyle="butt")
axs[-1].set_xticks(range(96)); axs[-1].set_xticklabels([c[0] + c[2] + c[6] for c in CHANNELS], rotation=90, fontsize=4, family="monospace"); axs[-1].set_xlabel("Trinucleotide context (pyrimidine frame)")
fig.savefig(f"{FIG}/Fig3_pooled_spectra.png"); plt.close()
# ---------------- Fig 4 composition + consensus
fig, axs = plt.subplots(1, 2, figsize=(11.5, 4.2), gridspec_kw={"width_ratios": [1.25, 1], "wspace": 0.62})
ax = axs[0]; keys4 = ["ASI: Cutaneous/H&N", "ASI: Breast", "ASI: Visceral/soft tissue/other", "ACC", "THYM: thymoma", "THYM: thymic carcinoma"]; fams = ["UV", "Clock", "MMR", "HRD", "APOBEC"]
y0 = np.arange(len(keys4))[::-1]; h = 0.15
for j, f in enumerate(fams):
    v = np.array([[cons["pooled"][k]["families"][f][q] for q in ("frac", "lo", "hi")] for k in keys4]); yy = y0 + (2 - j) * h
    ax.errorbar(v[:, 0], yy, xerr=[v[:, 0] - v[:, 1], v[:, 2] - v[:, 0]], fmt="o", color=FC[f], ms=5, lw=1.5, capsize=0, label=f)
ax.set_yticks(y0); ax.set_yticklabels([f"{k.replace('ASI: ', 'ASI – ').replace('THYM: ', 'THYM – ')}\n({cons['pooled'][k]['n_snv']} SNVs)" for k in keys4], fontsize=7.5); ax.set_xlim(-0.02, 1); ax.set_xlabel("Fraction of pooled SNVs attributed (NNLS, panel frame; 95% bootstrap CI)")
ax.legend(ncol=5, loc="upper center", bbox_to_anchor=(0.45, -0.16), fontsize=7.5, handletextpad=0.2, columnspacing=1); ax.set_title("Pooled process composition"); ax.grid(axis="y", visible=False); tag(ax, "a")
for y in y0[:-1] - 0.5: ax.axhline(y, color=GRID, lw=0.8)
ax = axs[1]; an = b[b.FIT_CLASS != "DROP"].copy(); an["STR"] = np.where(an.GROUP == "ASI", "ASI – " + an.ASI_SITE.str.replace("Visceral/soft tissue/other", "Visceral/other"), an.GROUP)
order = ["ASI – Cutaneous/H&N", "ASI – Breast", "ASI – Visceral/other", "ACC", "THYM"]; cats = ["UV", "Clock", "MMR", "APOBEC", "POLE", "Other", "AMBIGUOUS"]
T = pd.crosstab(an.STR, an.CONSENSUS).reindex(index=order, columns=cats, fill_value=0); left = np.zeros(len(order)); yy = np.arange(len(order))[::-1]
for c in cats:
    ax.barh(yy, T[c].values, left=left, color=FC[c], height=0.6, edgecolor="white", lw=1.5, label="Ambiguous (no ≥2/3 majority)" if c == "AMBIGUOUS" else c)
    for yv, l, w in zip(yy, left, T[c].values):
        if w >= 3: ax.text(l + w / 2, yv, str(w), ha="center", va="center", fontsize=7, color="#0b0b0b" if c in ("AMBIGUOUS", "Clock", "POLE") and c != "Clock" else "white")
    left += T[c].values
ax.set_yticks(yy); ax.set_yticklabels([f"{o} (n = {int(T.loc[o].sum())})" for o in order], fontsize=7.5); ax.set_xlabel("Analysable tumours (≥2 SNVs)"); ax.grid(axis="y", visible=False)
ax.legend(ncol=4, loc="upper center", bbox_to_anchor=(0.4, -0.16), fontsize=7, handlelength=1, handletextpad=0.3, columnspacing=0.9); ax.set_title("Per-tumour consensus dominant process"); tag(ax, "b")
fig.savefig(f"{FIG}/Fig4_composition_consensus.png"); plt.close()
# ---------------- Fig 5 concordance + positive controls
fig, axs = plt.subplots(1, 3, figsize=(11, 3.3), gridspec_kw={"wspace": 0.45})
ax = axs[0]; sets = [("basket_ge2", "Basket\n≥2 SNVs"), ("basket_plus_controls_SPA_set", "Basket +\ncontrols"), ("pan_cohort_ge20_(M2~M3 only)", "≥20 SNVs")]; prs = [("M1_SPA~M3_cosine", "SPA vs cosine", PAL[0]), ("M1_SPA~M2_MSKv2", "SPA vs MSK v2", PAL[1]), ("M2_MSKv2~M3_cosine", "MSK v2 vs cosine", PAL[2])]
for j, (p, lab, c) in enumerate(prs):
    v = [cons["kappa"][s][p]["called_only"]["kappa"] for s, _ in sets]; xx = np.arange(3) + (j - 1) * 0.25; ax.bar(xx, v, width=0.22, color=c, label=lab)
    for x_, v_ in zip(xx, v): ax.text(x_, v_ + 0.01, f"{v_:.2f}", ha="center", fontsize=6.5, color=INK)
ax.axhline(0.6, color="#52514e", ls="--", lw=0.8); ax.text(2.45, 0.61, "gate 0.60", fontsize=7, color=MUT, ha="right"); ax.set_xticks(range(3)); ax.set_xticklabels([l for _, l in sets]); ax.set_ylabel("Cohen's κ (called samples)"); ax.set_ylim(0, 0.95)
ax.legend(fontsize=7, loc="upper left", ncol=1); ax.set_title("Inter-method concordance"); ax.grid(axis="x", visible=False); tag(ax, "a")
ax = axs[1]; ctr = [("UV_in_cutaneous_melanoma_vs_non_skin", "UV\nmelanoma"), ("MMR_in_MSI_instable_vs_stable", "MMR\nMSI-H"), ("APOBEC_in_bladder_vs_colorectal", "APOBEC\nbladder"), ("Tobacco_in_NSCLC_vs_colorectal", "Tobacco\nNSCLC")]
for j, (m, lab, c) in enumerate([("M2_MSKv2", "MSK v2", PAL[1]), ("M3_cosine", "Cosine", PAL[2]), ("CONSENSUS", "Consensus", PAL[6])]):
    xx = np.arange(4) + (j - 1) * 0.25; ax.bar(xx, [cons["controls"][k][m]["sens"] for k, _ in ctr], width=0.22, color=c, label=lab)
    ax.plot(xx, [cons["controls"][k][m]["spec"] for k, _ in ctr], "D", color=INK, ms=3.5, label="Specificity" if j == 2 else None)
ax.set_xticks(range(4)); ax.set_xticklabels([l for _, l in ctr], fontsize=7); ax.set_ylabel("Sensitivity (bars) / specificity (◆)"); ax.set_ylim(0, 1.25); ax.set_yticks(np.arange(0, 1.01, 0.2)); ax.legend(fontsize=7, ncol=4, loc="upper center", handlelength=1, columnspacing=0.8, handletextpad=0.3)
ax.set_title("Pan-cohort positive controls (≥4 SNVs)"); ax.grid(axis="x", visible=False); tag(ax, "b")
ax = axs[2]; ks = [">=4", ">=10", ">=20"]; ax.plot(range(3), [cons["AUC_UVfraction_SKCM_vs_rest"][k] for k in ks], "-o", color=FC["UV"], lw=2, label="UV fraction → cutaneous melanoma")
ax.plot(range(3), [cons["AUC_MMRfraction_vs_MSI"][k]["auc"] for k in ks], "-o", color=FC["MMR"], lw=2, label="MMR fraction → MSI-unstable")
for i, k in enumerate(ks): ax.text(i, cons["AUC_UVfraction_SKCM_vs_rest"][k] + 0.008, f"{cons['AUC_UVfraction_SKCM_vs_rest'][k]:.2f}", ha="center", fontsize=7); ax.text(i, cons["AUC_MMRfraction_vs_MSI"][k]["auc"] - 0.022, f"{cons['AUC_MMRfraction_vs_MSI'][k]['auc']:.2f}", ha="center", fontsize=7)
ax.set_xticks(range(3)); ax.set_xticklabels(["≥4", "≥10", "≥20"]); ax.set_xlabel("Minimum SNVs per tumour"); ax.set_ylabel("ROC AUC"); ax.set_ylim(0.78, 1.0); ax.set_xlim(-0.3, 2.3); ax.legend(fontsize=7, loc="lower right"); ax.set_title("Discrimination vs SNV count"); tag(ax, "c")
fig.savefig(f"{FIG}/Fig5_concordance_controls.png"); plt.close()
# ---------------- Fig 6 survival
s = pd.read_csv(f"{RES}/survival_table.tsv", sep="\t"); fig, axs = plt.subplots(1, 3, figsize=(12.5, 3.5), gridspec_kw={"wspace": 0.3, "width_ratios": [1, 1, 0.9]})
def km(ax, d, lab, c):
    k = KaplanMeierFitter().fit(d["T"], d.E); sf = k.survival_function_.iloc[:, 0]; ci = k.confidence_interval_
    ax.step(sf.index, sf.values, where="post", color=c, lw=2, label=f"{lab} (n = {len(d)}, {int(d.E.sum())} deaths)"); ax.fill_between(ci.index, ci.iloc[:, 0], ci.iloc[:, 1], step="post", color=c, alpha=0.12, lw=0)
ax = axs[0]
for g, c in [("THYM", PAL[2]), ("ACC", PAL[1]), ("ASI", PAL[0])]: km(ax, s[s.GROUP == g], g, c)
ax.set_xlim(0, 72); ax.set_ylim(0, 1.02); ax.set_xlabel("Months from sequencing"); ax.set_ylabel("Overall survival"); ax.legend(fontsize=7, loc="lower left"); ax.set_title(f"Tumour group (log-rank P = {st['survival']['logrank_groups_p']:.4f})"); tag(ax, "a")
ax = axs[1]; d = s[(s.GROUP == "ASI") & (s.FIT_CLASS != "DROP")]; km(ax, d[d.UV == 1], "UV-dominant", FC["UV"]); km(ax, d[d.UV == 0], "Other / ambiguous", "#898781")
ax.set_xlim(0, 72); ax.set_ylim(0, 1.02); ax.set_xlabel("Months from sequencing"); ax.legend(fontsize=7, loc="lower left"); ax.set_title(f"ASI by consensus UV (log-rank P = {st['survival']['ASI_UV_logrank']['p']:.2f})"); tag(ax, "b")
ax = axs[2]; C = {c["label"][0:2].strip(":"): c for c in st["survival"]["cox"]}; items = [("B", "UV", "ASI: UV-dominant vs not"), ("C", "CUT", "ASI: cutaneous/H&N site"), ("C", "MET", "ASI: metastatic specimen"), ("D", "AGE", "ACC: age (per year)"), ("D", "log10TMB", "ACC: log10(TMB+1)"), ("E", "log10TMB", "THYM: log10(TMB+1)"), ("A", "G_THYM", "THYM vs ASI (pooled model)")]
for i, (m, v, lab) in enumerate(items[::-1]):
    h_ = C[m]["hr"][v]; ax.plot([h_["exp(coef) lower 95%"], h_["exp(coef) upper 95%"]], [i, i], color=INK, lw=1.5); ax.plot(h_["exp(coef)"], i, "s", color=PAL[0], ms=6)
    ax.text(40, i, f"{h_['exp(coef)']:.2f} ({h_['exp(coef) lower 95%']:.2f}–{h_['exp(coef) upper 95%']:.2f})", va="center", fontsize=6.5, color=INK)
ax.axvline(1, color="#52514e", lw=0.8, ls="--"); ax.set_xscale("log"); ax.set_xlim(0.15, 600); ax.set_xticks([0.25, 1, 4, 16]); ax.set_xticklabels(["0.25", "1", "4", "16"]); ax.set_yticks(range(len(items))); ax.set_yticklabels([l for _, _, l in items[::-1]], fontsize=7)
ax.set_xlabel("Hazard ratio (95% CI), multivariable Cox"); ax.grid(axis="y", visible=False); ax.set_title("Adjusted hazard ratios"); ax.text(-0.75, 1.08, "c", transform=ax.transAxes, fontsize=11, fontweight="bold"); p_=ax.get_position(); ax.set_position([p_.x0 + 0.075, p_.y0, p_.width - 0.075, p_.height])
fig.savefig(f"{FIG}/Fig6_survival.png"); plt.close()
# ---------------- Fig 7 external validation
fig, axs = plt.subplots(1, 3, figsize=(11, 3.3), gridspec_kw={"wspace": 0.42}); rng = np.random.default_rng(1)
for ax, (study, ttl, t) in zip(axs[:2], [("angs_painter_2025", "Angiosarcoma Project 2025 (WES)", "a"), ("angs_project_painter_2018", "Angiosarcoma Project 2020 (WES)", "b")]):
    e = pd.read_csv(f"{RES}/external_{study}.tsv", sep="\t", index_col=0); e = e[e.n_wes >= 20]; e["pt"] = e.index.str.extract(r"^(.*?)(?:[-_]T?\d+)?$")[0].values
    for i, (hn, lab) in enumerate([(True, "Head/neck/face/scalp"), (False, "Other sites")]):
        v = e[e.head_neck == hn].UV_frac; ax.scatter(i + rng.uniform(-0.16, 0.16, len(v)), v, s=22, color=FC["UV"] if hn else "#898781", edgecolor="white", lw=0.6, alpha=0.9); ax.plot([i - 0.25, i + 0.25], [v.median()] * 2, color=INK, lw=2)
    ax.set_xticks([0, 1]); ax.set_xticklabels(["Head/neck/\nface/scalp", "Other sites"]); ax.set_xlim(-0.6, 1.6); ax.set_ylim(-0.03, 1.03); ax.set_ylabel("UV fraction of exome SNVs (NNLS)"); ax.grid(axis="x", visible=False)
    ax.set_title(f"{ttl}\nFisher P = {ext[study]['fisher_p']:.1e} (UV-dominant × site)"); tag(ax, t)
ax = axs[2]; labs = [("angs_painter_2025", "ANGS\n2025"), ("angs_project_painter_2018", "ANGS\n2020"), ("acc_tcga_pan_can_atlas_2018", "TCGA\nACC"), ("thym_tcga_pan_can_atlas_2018", "TCGA\nTHYM")]
cr = [ext[k]["panel_downsampling"]["cosine_call_rate"] for k, _ in labs]; ac = [ext[k]["panel_downsampling"]["cosine_acc_when_called"] for k, _ in labs]; xx = np.arange(4)
ax.bar(xx - 0.2, cr, width=0.36, color=PAL[0], label="Call rate"); ax.bar(xx + 0.2, ac, width=0.36, color=PAL[2], label="Accuracy when called")
for x_, a_, c_, (k, _) in zip(xx, ac, cr, labs): ax.text(x_ + 0.2, a_ + 0.015, f"{a_:.2f}", ha="center", fontsize=7); ax.text(x_ - 0.2, c_ + 0.015, f"{c_:.2f}", ha="center", fontsize=7)
ax.set_xticks(xx); ax.set_xticklabels([f"{l}\n(n = {ext[k]['panel_downsampling']['n_panel_ge2']})" for k, l in labs], fontsize=6); ax.set_ylim(0, 1.3); ax.set_yticks(np.arange(0, 1.01, 0.2)); ax.legend(fontsize=7, ncol=2, loc="upper center"); ax.grid(axis="x", visible=False)
ax.set_title("In-silico panel down-sampling of exomes\n(cosine classifier vs exome truth)"); tag(ax, "c")
fig.savefig(f"{FIG}/Fig7_external_validation.png"); plt.close(); print("figures done")
