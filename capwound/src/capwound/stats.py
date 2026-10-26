"""Statistical analysis helpers."""

from __future__ import annotations

import numpy as np
import pandas as pd
from scipy import stats
from sklearn.cross_decomposition import PLSRegression
from sklearn.model_selection import GridSearchCV, KFold, LeaveOneOut
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler


def anova_tukey(df: pd.DataFrame, value: str, group: str, order=None):
    order = order or list(dict.fromkeys(df[group]))
    samples = [df.loc[df[group] == g, value].to_numpy() for g in order]
    F, p = stats.f_oneway(*samples)
    k, N = len(samples), sum(len(s) for s in samples)
    ss_b = sum(len(s) * (s.mean() - np.concatenate(samples).mean()) ** 2 for s in samples)
    ss_t = ((np.concatenate(samples) - np.concatenate(samples).mean()) ** 2).sum()
    eta2 = ss_b / ss_t
    tk = stats.tukey_hsd(*samples)
    rows = []
    for i in range(k):
        for j in range(i + 1, k):
            rows.append({"A": order[i], "B": order[j],
                         "diff": samples[j].mean() - samples[i].mean(),
                         "p_adj": tk.pvalue[i, j]})
    return {"F": F, "p": p, "eta2": eta2, "df_between": k - 1, "df_within": N - k,
            "pairs": pd.DataFrame(rows)}


def spearman_table(df: pd.DataFrame, features, target):
    rows = []
    for f in features:
        rho, p = stats.spearmanr(df[f], df[target])
        rows.append({"feature": f, "rho": rho, "p": p})
    return pd.DataFrame(rows)


def nested_pls(X, y, max_comp=8, seed=0):
    """Leave-one-out outer loop, 5-fold inner loop for the number of PLS components."""
    y = np.asarray(y)
    outer = LeaveOneOut()
    pred = np.empty_like(y, dtype=float)
    chosen = []
    for tr, te in outer.split(X):
        pipe = make_pipeline(StandardScaler(), PLSRegression(scale=False))
        gs = GridSearchCV(pipe, {"plsregression__n_components": list(range(1, max_comp + 1))},
                          cv=KFold(5, shuffle=True, random_state=seed),
                          scoring="neg_root_mean_squared_error")
        gs.fit(X[tr], y[tr])
        pred[te] = gs.predict(X[te]).ravel()
        chosen.append(gs.best_params_["plsregression__n_components"])
    rmse = np.sqrt(np.mean((pred - y) ** 2))
    q2 = 1 - np.sum((pred - y) ** 2) / np.sum((y - y.mean()) ** 2)
    final = make_pipeline(StandardScaler(), PLSRegression(n_components=int(np.median(chosen)), scale=False))
    final.fit(X, y)
    return {"pred": pred, "rmsecv": rmse, "q2": q2, "n_comp": int(np.median(chosen)),
            "model": final}


def vip_scores(pls: PLSRegression):
    t = pls.x_scores_
    w = pls.x_weights_
    q = pls.y_loadings_
    p, h = w.shape
    s = np.diag(t.T @ t @ q.T @ q).reshape(h, -1)
    total = s.sum()
    wn = (w / np.linalg.norm(w, axis=0)) ** 2
    return np.sqrt(p * (wn @ s).ravel() / total)


def linear_lod(x, y):
    """Calibration slope/intercept, R^2 and 3.3*s/slope detection limit."""
    res = stats.linregress(x, y)
    resid = y - (res.intercept + res.slope * x)
    s = np.std(resid, ddof=2)
    return {"slope": res.slope, "intercept": res.intercept, "r2": res.rvalue ** 2,
            "p": res.pvalue, "lod": 3.3 * s / abs(res.slope), "s_res": s}
