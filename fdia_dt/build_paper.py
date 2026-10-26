"""Fill %%KEY%% placeholders in the paper sources from results/digest.json."""
import json
import re

d = json.load(open("results/digest.json"))
a, b = d["ieee14"], d["ieee118"]
P = "DT-GCN (proposed)"
loc118 = {m: v["loc_f1"] for m, v in b.items() if "loc_f1" in v and m != P}
vals = {
    "BDDNAIVE14": d["ieee14_naive_recall"]["BDD (chi-square)"],
    "BDDNAIVE118": d["ieee118_naive_recall"]["BDD (chi-square)"],
    "BDDST14": d["ieee14_bdd_stealthy_recall"],
    "BDDST118": d["ieee118_bdd_stealthy_recall"],
    "F14": a[P]["f1"], "AUC14": a[P]["auc"],
    "F118": b[P]["f1"], "FPR118": b[P]["fpr"], "AUC118": b[P]["auc"],
    "GCNRAW118": b["GCN (raw)"]["f1"],
    "DTTH118": b["DT-Threshold"]["f1"], "DTMLP118": b["DT-MLP"]["f1"],
    "K118": "/".join(f"{v:.1f}" for v in d["ieee118_by_k"][P]),
    "K118RF": "/".join(f"{v:.1f}" for v in d["ieee118_by_k"]["RandomForest (raw)"]),
    "LOCF118": b[P]["loc_f1"], "LOCEX118": b[P]["loc_exact"],
    "LOCF118NEXT": max(loc118.values()),
    "LOCF14": a[P]["loc_f1"], "LOCF14MLP": a["DT-MLP"]["loc_f1"], "LOCEX14": a[P]["loc_exact"],
    "INF118": f'{d["ieee118_infer_ms"][P]:.3f}',
    "SENS": open("paper/sens_text.tex").read().strip(),
}


def fill(txt):
    def rep(m):
        v = vals[m.group(1)]
        return f"{v:.2f}" if isinstance(v, float) else str(v)
    return re.sub(r"%%([A-Z0-9]+)%%", rep, txt)


for src, dst in [("paper/sec_results.tex", "paper/sec_results_filled.tex"),
                 ("paper/abstract.tex", "paper/abstract_filled.tex")]:
    out = fill(open(src).read())
    assert "%%" not in out, src
    open(dst, "w").write(out)
print("filled", len(vals), "values")
