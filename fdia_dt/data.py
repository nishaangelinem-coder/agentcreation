"""Scenario, digital-twin and stealthy AC-FDIA dataset generation."""
import numpy as np
from .grid import Grid

DT_STEP_H = 0.25  # 15-min resolution


def load_profile(T, rng):
    """Normalised system load shape: daily double peak, weekly cycle, AR(1) noise."""
    h = (np.arange(T) * DT_STEP_H) % 24
    day = (np.arange(T) * DT_STEP_H // 24).astype(int)
    shape = (0.62 + 0.18 * np.exp(-((h - 9.5) ** 2) / 8) + 0.28 * np.exp(-((h - 19.5) ** 2) / 6)
             - 0.05 * np.exp(-((h - 3.5) ** 2) / 6))
    weekly = np.where(day % 7 >= 5, 0.9, 1.0)
    ar = np.zeros(T)
    for k in range(1, T):
        ar[k] = 0.95 * ar[k - 1] + rng.normal(0, 0.008)
    return shape * weekly * (1 + ar), h


def pv_profile(h, T, rng):
    clear = np.clip(np.sin(np.pi * (h - 6) / 13), 0, None)
    day = (np.arange(T) * DT_STEP_H // 24).astype(int)
    cloud_day = rng.uniform(0.5, 1.0, day.max() + 1)[day]
    cloud = np.clip(cloud_day + rng.normal(0, 0.08, T), 0.1, 1.0)
    return clear * cloud


def operating_points(g: Grid, T, rng, pv_share=0.25, bus_sigma=0.03):
    """True bus injections for T steps. Returns Pd, Qd (T,n) and Pg (T,ngen)."""
    sysload, h = load_profile(T, rng)
    n = g.n
    loadbus = np.where(g.Pd0 > 0)[0]
    # PV plants at ~30% of load buses
    pvbus = rng.choice(loadbus, max(1, int(0.3 * len(loadbus))), replace=False)
    pvcap = np.zeros(n)
    pvcap[pvbus] = pv_share * g.Pd0.sum() / len(pvbus)
    pv = pv_profile(h, T, rng)
    # bus-level idiosyncratic AR(1) variation
    eps = np.zeros((T, n))
    for k in range(1, T):
        eps[k] = 0.9 * eps[k - 1] + rng.normal(0, bus_sigma * np.sqrt(1 - 0.81), n)
    mult = sysload[:, None] * (1 + eps) / 0.85
    Pd = mult * g.Pd0[None, :] - pv[:, None] * pvcap[None, :]
    Qd = mult * g.Qd0[None, :]
    # generators follow net load proportionally (slack balances losses)
    net = Pd.sum(1) / g.Pd0.sum()
    Pg = net[:, None] * g.Pg0[None, :]
    return Pd, Qd, Pg, h


def stealthy_attack(g: Grid, Va, Vm, rng, k_range=(1, 3), mag=(0.001, 0.02), vfrac=0.5):
    """AC stealthy FDIA: a = h(x+c) - h(x) on a connected region of k buses."""
    k = rng.integers(k_range[0], k_range[1] + 1)
    cand = g.nonref
    start = rng.choice(cand)
    region = [start]
    while len(region) < k:
        nb = np.unique(np.where(g.adj[region].sum(0) > 0)[0])
        nb = [b for b in nb if b not in region and b != g.ref]
        if not nb:
            break
        region.append(rng.choice(nb))
    region = np.array(region)
    m = rng.uniform(*mag)
    c_a = np.zeros(g.n)
    c_v = np.zeros(g.n)
    c_a[region] = rng.choice([-1, 1]) * m * rng.uniform(0.8, 1.0, len(region))
    if rng.random() < vfrac:
        c_v[region] = rng.choice([-1, 1]) * rng.uniform(0.2, 0.5) * m
    a = g.h(Va + c_a, Vm + c_v) - g.h(Va, Vm)
    label = np.zeros(g.n)
    label[region] = 1
    return a, label, m


def naive_attack(g: Grid, rng, sig):
    """Unstructured random bias on a few measurements (detectable by BDD)."""
    a = np.zeros(g.nz)
    idx = rng.choice(g.nz, rng.integers(3, 8), replace=False)
    a[idx] = rng.choice([-1, 1], len(idx)) * rng.uniform(5, 15, len(idx)) * sig[idx]
    return a


def digital_twin(g: Grid, Pd, Qd, Pg, rng, fsig):
    """Twin state from forecast injections (very-short-term forecast error fsig)."""
    ef = 1 + rng.normal(0, fsig, Pd.shape)
    Va, Vm, ok = g.powerflow(Pd * ef, Qd * ef, Pg * (1 + rng.normal(0, fsig / 4, Pg.shape)))
    return Va, Vm


def generate(name, T=6000, attack_rate=0.5, naive_rate=0.1, fsig=0.02, seed=0, verbose=True):
    rng = np.random.default_rng(seed)
    g = Grid(name)
    sig = g.sigma()
    Pd, Qd, Pg, hours = operating_points(g, T, rng)
    out = {k: [] for k in ["z", "Va_hat", "Vm_hat", "Va_dt", "Vm_dt", "J", "y", "ynode",
                           "mag", "kind", "Va_true", "Vm_true", "fc_Pd", "fc_Qd", "fc_Pg"]}
    Va_prev = Vm_prev = None
    for t in range(T):
        Va, Vm, ok = g.powerflow(Pd[t], Qd[t], Pg[t])
        z = g.h(Va, Vm) + rng.normal(0, sig)
        # attacker's knowledge: SE on clean telemetry
        Va_c, Vm_c, J_c, _ = g.wls(z, sig, Va_prev, Vm_prev)
        u = rng.random()
        ynode = np.zeros(g.n)
        mag, kind = 0.0, 0
        if u < attack_rate:
            a, ynode, mag = stealthy_attack(g, Va_c, Vm_c, rng)
            kind = 1
        elif u < attack_rate + naive_rate:
            a = naive_attack(g, rng, sig)
            kind = 2
        else:
            a = None
        if a is not None:
            z = z + a
            Va_h, Vm_h, J, _ = g.wls(z, sig, Va_c, Vm_c)
        else:
            Va_h, Vm_h, J = Va_c, Vm_c, J_c
        Va_d, Vm_d = digital_twin(g, Pd[t], Qd[t], Pg[t], rng, fsig)
        for k, v in zip(["z", "Va_hat", "Vm_hat", "Va_dt", "Vm_dt", "J", "y", "ynode", "mag", "kind",
                         "Va_true", "Vm_true"],
                        [z, Va_h, Vm_h, Va_d, Vm_d, J, int(kind > 0), ynode, mag, kind, Va, Vm]):
            out[k].append(v)
        Va_prev, Vm_prev = Va_c, Vm_c
        if verbose and t % 1000 == 0:
            print(f"  [{name}] step {t}/{T}", flush=True)
    data = {k: np.array(v) for k, v in out.items() if v}
    data["Pd"], data["Qd"], data["Pg"], data["hours"] = Pd, Qd, Pg, hours
    data["sig"] = sig
    return g, data


def redo_twin(g, data, fsig, seed):
    """Recompute the twin for a different forecast-error level (sensitivity study)."""
    rng = np.random.default_rng(seed)
    T = len(data["y"])
    Va_d = np.zeros((T, g.n))
    Vm_d = np.zeros((T, g.n))
    for t in range(T):
        Va_d[t], Vm_d[t] = digital_twin(g, data["Pd"][t], data["Qd"][t], data["Pg"][t], rng, fsig)
    d = dict(data)
    d["Va_dt"], d["Vm_dt"] = Va_d, Vm_d
    return d


def node_features(g, data):
    """Per-bus DT-residual features (T, n, 6)."""
    Va_h, Vm_h = data["Va_hat"], data["Vm_hat"]
    Va_d, Vm_d = data["Va_dt"], data["Vm_dt"]
    z = data["z"]
    n, m = g.n, g.m
    T = len(z)
    zdt = np.stack([g.h(Va_d[t], Vm_d[t]) for t in range(T)])
    r = z - zdt
    rP = r[:, n:2 * n]
    rQ = r[:, 2 * n:3 * n]
    rV = r[:, :n]
    rPf = np.abs(r[:, 3 * n:3 * n + m])
    # incident-line flow residual aggregated to buses
    inc = np.zeros((m, n))
    inc[np.arange(m), g.f] = 1
    inc[np.arange(m), g.t] = 1
    rflow = rPf @ inc
    # angle deviation relative to neighbours removes the common drift from the slack
    dA = Va_h - Va_d
    deg = g.adj.sum(1)
    dA_rel = dA - (dA @ g.adj) / np.maximum(deg, 1)
    dV = Vm_h - Vm_d
    return np.stack([dA, dA_rel, dV, rP, rQ, rflow], axis=-1).astype(np.float32), rV


def raw_node_features(g, data):
    """Per-bus raw telemetry (no twin): V, P_inj, Q_inj, incident flows (T, n, 4)."""
    z = data["z"]
    n, m = g.n, g.m
    inc = np.zeros((m, n))
    inc[np.arange(m), g.f] = 1
    inc[np.arange(m), g.t] = 1
    Pf = np.abs(z[:, 3 * n:3 * n + m]) @ inc
    return np.stack([z[:, :n], z[:, n:2 * n], z[:, 2 * n:3 * n], Pf], axis=-1).astype(np.float32)
