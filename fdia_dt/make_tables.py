"""Generate LaTeX tables and a JSON digest of key numbers from results/."""
import json
import pickle
import numpy as np
import pandas as pd

CASES = [("ieee14", "IEEE 14-bus"), ("ieee118", "IEEE 118-bus")]
digest = {}


def fmt(m, s, pct=True, best=False):
    k = 100 if pct else 1
    txt = f"{m * k:.2f}" + (f"$\\pm${s * k:.2f}" if s > 0 else "")
    return f"\\textbf{{{txt}}}" if best else txt


def table_detection():
    lines = []
    means = {c: pd.read_csv(f"results/mean_{c}.csv", index_col=0) for c, _ in CASES}
    stds = {c: pd.read_csv(f"results/std_{c}.csv", index_col=0) for c, _ in CASES}
    methods = list(means["ieee14"].index)
    cols = ["acc", "f1", "fpr", "auc"]
    for mt in methods:
        cells = [mt.replace("(proposed)", "(prop.)")]
        for c, _ in CASES:
            for k in cols:
                col = means[c][k]
                best = (col.idxmin() if k == "fpr" else col.idxmax()) == mt
                cells.append(fmt(col[mt], stds[c].loc[mt, k], best=best))
        lines.append(" & ".join(cells) + r" \\")
    open("paper/tab_detection.tex", "w").write("\n".join(lines))
    for c, _ in CASES:
        digest[c] = {mt: {k: round(float(means[c].loc[mt, k]) * 100, 2) for k in means[c].columns
                          if not np.isnan(means[c].loc[mt, k])} for mt in methods}
        digest[c + "_std"] = {mt: {k: round(float(stds[c].loc[mt, k]) * 100, 2) for k in stds[c].columns
                                   if not np.isnan(stds[c].loc[mt, k])} for mt in methods}


def table_localisation():
    lines = []
    means = {c: pd.read_csv(f"results/mean_{c}.csv", index_col=0) for c, _ in CASES}
    stds = {c: pd.read_csv(f"results/std_{c}.csv", index_col=0) for c, _ in CASES}
    methods = [m for m in means["ieee14"].index if not np.isnan(means["ieee14"].loc[m, "loc_f1"])]
    cols = ["loc_rec", "loc_f1", "loc_exact"]
    for mt in methods:
        cells = [mt.replace("(proposed)", "(prop.)")]
        for c, _ in CASES:
            for k in cols:
                col = means[c].loc[methods, k]
                cells.append(fmt(col[mt], stds[c].loc[mt, k], best=col.idxmax() == mt))
        lines.append(" & ".join(cells) + r" \\")
    open("paper/tab_localisation.tex", "w").write("\n".join(lines))


def table_k_and_time():
    lines = []
    for c, title in CASES:
        s = pickle.load(open(f"results/store_{c}.pkl", "rb"))
        te = s["te"]; k = s["kind"][te]; N = s["N"][te]
        nk = N.sum(1)
        row = {}
        for m in ["RandomForest (raw)", "GCN (raw)", "DT-MLP", "DT-GCN (proposed)"]:
            pred = s["scores"][m] > 0.5
            row[m] = [pred[(k == 1) & (nk == j)].mean() for j in (1, 2, 3)]
        row["naive"] = {m: float((s["scores"][m] > 0.5)[k == 2].mean()) for m in row}
        row["naive"]["BDD (chi-square)"] = float((s["scores"]["BDD (chi-square)"] > s["thr"])[k == 2].mean())
        row["stealthy_bdd"] = float((s["scores"]["BDD (chi-square)"] > s["thr"])[k == 1].mean())
        digest[c + "_by_k"] = {m: [round(v * 100, 2) for v in r] for m, r in row.items() if isinstance(r, list)}
        digest[c + "_naive_recall"] = {m: round(v * 100, 2) for m, v in row["naive"].items()}
        digest[c + "_bdd_stealthy_recall"] = round(row["stealthy_bdd"] * 100, 2)
        digest[c + "_infer_ms"] = {m: round(v * 1000, 4) for m, v in s["timing"].items()}
        digest[c + "_test_counts"] = {"normal": int((k == 0).sum()), "stealthy": int((k == 1).sum()),
                                      "naive": int((k == 2).sum())}
        for j in (1, 2, 3):
            cells = [title if j == 1 else "", str(j)] + [f"{row[m][j - 1] * 100:.2f}" for m in
                                                         ["RandomForest (raw)", "GCN (raw)", "DT-MLP",
                                                          "DT-GCN (proposed)"]]
            lines.append(" & ".join(cells) + r" \\")
        lines.append(r"\midrule" if c == "ieee14" else "")
    open("paper/tab_k.tex", "w").write("\n".join(lines))


def sens_digest():
    for c, _ in CASES:
        try:
            df = pd.read_csv(f"results/sens_{c}.csv")
        except FileNotFoundError:
            continue
        g = df.groupby(["method", "level"])[["f1", "loc_f1"]].mean().round(4) * 100
        digest[c + "_sens"] = {f"{m}@{lv}": v for (m, lv), v in g.to_dict("index").items()}


if __name__ == "__main__":
    table_detection(); table_localisation(); table_k_and_time(); sens_digest()
    json.dump(digest, open("results/digest.json", "w"), indent=1)
    print(json.dumps(digest, indent=1))
