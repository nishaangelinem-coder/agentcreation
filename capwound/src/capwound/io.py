"""Load measured spectra so the fitting/statistics pipeline runs on real lab data.

Expected CSV columns (header row required):
    freq_Hz, Zreal_ohm, Zimag_ohm          (one file per sample)
Optional sample sheet ``samples.csv`` with columns: file, group, replicate, log10_cfu
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from .fitting import fit_dielectric, fit_randles, lin_kk


def read_spectrum(path):
    d = pd.read_csv(path)
    d = d.sort_values("freq_Hz")
    return d["freq_Hz"].to_numpy(float), (d["Zreal_ohm"] + 1j * d["Zimag_ohm"]).to_numpy()


def analyse_folder(folder, K=600.0, kind="bioimpedance"):
    """Fit every spectrum listed in ``folder/samples.csv``; returns a DataFrame."""
    folder = Path(folder)
    sheet = pd.read_csv(folder / "samples.csv")
    rows = []
    for _, s in sheet.iterrows():
        f, z = read_spectrum(folder / s["file"])
        M, mu, rr, ri = lin_kk(f, z)
        res = fit_dielectric(f, z, K) if kind == "bioimpedance" else fit_randles(f, z)
        res.pop("z_fit")
        rows.append({**s.to_dict(), **res, "kk_rms": float(np.sqrt(np.mean(rr ** 2 + ri ** 2)))})
    return pd.DataFrame(rows)
