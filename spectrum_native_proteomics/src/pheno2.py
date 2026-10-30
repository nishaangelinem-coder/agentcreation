"""Spectrum-to-phenotype benchmark (feature-token formulation).

Pipelines compared
  ID-first      : DDA identifications of library peptides -> quant matrix ->
                  median normalisation -> min imputation -> L2 logistic regr.
  ID-free-LR    : all MS1 features aligned to a consensus map (no MS2 IDs)
                  -> same statistics.
  FeatureFormer : ProFormer-style transformer on ID-free feature tokens,
                  trained from scratch.
  FeatureFormer+SSL : same, pretrained with masked-intensity modelling on
                  unlabeled runs (MassNet-style corpus) then fine-tuned.
"""
import numpy as np, torch, torch.nn.functional as F
from sklearn.linear_model import LogisticRegressionCV
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import roc_auc_score, balanced_accuracy_score
from simulator import CohortSimulator, acquire_features, align_to_reference
from models import FeatureFormer


def build_matrices(sim, ab, b, rng, ref_mz, ref_rt):
    n = len(ab)
    Vfree = np.full((n, sim.n_feat), np.nan); Vid = np.full((n, sim.n_feat), np.nan)
    shifts = np.zeros(n); align_acc = []
    for i in range(n):
        mz, rt, it, tid, ident = acquire_features(sim, ab[i], b[i], rng)
        Vfree[i], shifts[i] = align_to_reference(mz, rt, it, ref_mz, ref_rt)
        # ID-first: only identified library peptides are quantified
        obs = {t: v for t, v in zip(tid, it) if t >= 0}
        for j in np.where(ident)[0]:
            Vid[i, j] = obs.get(j, np.nan)
        truth = np.full(sim.n_feat, np.nan); truth[tid[tid >= 0]] = it[tid >= 0]
        both = ~np.isnan(truth) & ~np.isnan(Vfree[i])
        align_acc.append(np.mean(np.isclose(truth[both], Vfree[i][both])) if both.any() else 0)
    return Vfree, Vid[:, sim.in_lib], shifts, float(np.mean(align_acc))


def prep_stats(Vtr, Vte):
    keep = np.mean(~np.isnan(Vtr), 0) >= 0.3
    Vtr, Vte = Vtr[:, keep], Vte[:, keep]
    med = lambda V: V - np.nanmedian(V, 1, keepdims=True)
    Vtr, Vte = med(Vtr), med(Vte)
    mn = np.nanmin(Vtr, 0) - 0.5
    Vtr = np.where(np.isnan(Vtr), mn, Vtr); Vte = np.where(np.isnan(Vte), mn, Vte)
    return Vtr, Vte


def lr_pipeline(Vtr, ytr, Vte):
    Xtr, Xte = prep_stats(Vtr, Vte)
    sc = StandardScaler().fit(Xtr)
    cv = int(min(5, np.bincount(ytr).min()))
    lr = LogisticRegressionCV(Cs=10, cv=cv, max_iter=5000, scoring='roc_auc')
    lr.fit(sc.transform(Xtr), ytr)
    return lr.predict_proba(sc.transform(Xte))[:, 1]


class Norm:
    """Per-run median centring + global standardisation for token values."""
    def __init__(self, V):
        c = V - np.nanmedian(V, 1, keepdims=True)
        c = np.where(np.isnan(c).all(0, keepdims=True), 0.0, c)
        self.mu = np.nanmean(c, 0); self.sd = np.nanstd(c, 0) + 1e-3
        self.mu = np.nan_to_num(self.mu); self.sd = np.nan_to_num(self.sd, nan=1.0)
    def __call__(self, V):
        c = (V - np.nanmedian(V, 1, keepdims=True) - self.mu) / self.sd
        miss = np.isnan(c)
        return torch.tensor(np.nan_to_num(c), dtype=torch.float32), torch.from_numpy(miss)


def make_model(sim, ref_mz, ref_rt):
    return FeatureFormer(sim.n_feat, ref_mz, ref_rt)


def pretrain(model, v, miss, rng, steps=400, bs=64, lr=1e-3):
    opt = torch.optim.AdamW(model.parameters(), lr=lr, weight_decay=0.05)
    sch = torch.optim.lr_scheduler.OneCycleLR(opt, max_lr=lr, total_steps=steps, pct_start=0.1)
    hist = []; model.train()
    for s in range(steps):
        i = torch.from_numpy(rng.integers(0, len(v), bs))
        loss = model.ssl_loss(v[i], miss[i])
        opt.zero_grad(); loss.backward(); opt.step(); sch.step()
        if s % 25 == 0 or s == steps - 1: hist.append((s, loss.item()))
    return hist


def finetune(model, v, miss, y, rng, steps=200, bs=32, lr=5e-4, body_lr=None):
    hp = list(model.head.parameters()) + list(model.tok_head.parameters())
    head = {id(p) for p in hp}
    body = [p for p in model.parameters() if id(p) not in head]
    groups = [dict(params=body, lr=body_lr or lr), dict(params=hp, lr=lr)]
    opt = torch.optim.AdamW(groups, weight_decay=0.05)
    sch = torch.optim.lr_scheduler.OneCycleLR(opt, max_lr=[g['lr'] for g in groups], total_steps=steps, pct_start=0.1)
    yt = torch.from_numpy(y).long(); model.train()
    for s in range(steps):
        i = torch.from_numpy(rng.choice(len(v), min(bs, len(v)), replace=False))
        vv = v[i] + 0.1 * torch.randn_like(v[i])                 # intensity jitter
        drop = (torch.rand_like(vv) < 0.1) | miss[i]               # feature dropout
        loss = F.cross_entropy(model(vv, drop), yt[i], label_smoothing=0.1)
        opt.zero_grad(); loss.backward(); opt.step(); sch.step()
    return model


@torch.no_grad()
def predict(model, v, miss):
    model.eval()
    return torch.cat([F.softmax(model(v[i:i + 128], miss[i:i + 128]), -1)[:, 1] for i in range(0, len(v), 128)]).numpy()


def evaluate(p, y):
    return dict(auroc=float(roc_auc_score(y, p)), bacc=float(balanced_accuracy_score(y, p > 0.5)))
