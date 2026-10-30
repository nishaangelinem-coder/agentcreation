"""Spectrum-to-phenotype benchmark on simulated LC-MS1 cohorts.

Compares an identify-then-quantify pipeline with spectrum-native models
(pixel logistic regression, CNN, MS1Former from scratch, MS1Former with
masked-autoencoder pretraining on unlabeled runs)."""
import time, json, sys, argparse
import numpy as np, torch, torch.nn as nn, torch.nn.functional as F
from sklearn.linear_model import LogisticRegressionCV
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import roc_auc_score, balanced_accuracy_score
from simulator import CohortSimulator
from models import MS1Former


class SmallCNN(nn.Module):
    def __init__(self):
        super().__init__()
        c = lambda i, o: nn.Sequential(nn.Conv2d(i, o, 3, padding=1), nn.BatchNorm2d(o), nn.GELU(), nn.MaxPool2d(2))
        self.f = nn.Sequential(c(1, 16), c(16, 32), c(32, 64), c(64, 64))
        self.h = nn.Sequential(nn.Dropout(0.3), nn.Linear(64, 2))

    def forward(self, x):
        return self.h(self.f(x.unsqueeze(1)).mean((2, 3)))


def augment(x, rng):
    sh = int(rng.integers(-2, 3))
    x = torch.roll(x, sh, 1)
    return x * float(rng.uniform(0.9, 1.1))


def train_clf(model, X, y, rng, steps=200, lr=3e-4, bs=32, wd=0.05, layer_lr=None):
    params = model.parameters() if layer_lr is None else layer_lr
    opt = torch.optim.AdamW(params, lr=lr, weight_decay=wd)
    mx = [g['lr'] for g in opt.param_groups]
    sched = torch.optim.lr_scheduler.OneCycleLR(opt, max_lr=mx, total_steps=steps, pct_start=0.1)
    Xt, yt = torch.from_numpy(X), torch.from_numpy(y).long()
    model.train()
    for s in range(steps):
        idx = rng.choice(len(X), min(bs, len(X)), replace=False)
        out = model(augment(Xt[idx], rng))
        loss = F.cross_entropy(out, yt[idx], label_smoothing=0.1)
        opt.zero_grad(); loss.backward(); opt.step(); sched.step()
    return model


@torch.no_grad()
def predict(model, X):
    model.eval()
    return torch.cat([F.softmax(model(torch.from_numpy(X[i:i + 128])), -1)[:, 1]
                      for i in range(0, len(X), 128)]).numpy()


def pretrain_mae(X, rng, steps=500, bs=48):
    m = MS1Former()
    opt = torch.optim.AdamW(m.parameters(), lr=1e-3, weight_decay=0.05)
    sched = torch.optim.lr_scheduler.OneCycleLR(opt, max_lr=1e-3, total_steps=steps, pct_start=0.1)
    Xt = torch.from_numpy(X); hist = []
    m.train()
    for s in range(steps):
        idx = rng.integers(0, len(X), bs)
        loss = m.mae_loss(augment(Xt[idx], rng))
        opt.zero_grad(); loss.backward(); opt.step(); sched.step()
        if s % 50 == 0 or s == steps - 1: hist.append((s, loss.item()))
    return m, hist


def evaluate(p, y):
    return dict(auroc=float(roc_auc_score(y, p)), bacc=float(balanced_accuracy_score(y, p > 0.5)))


def id_first(sim, Xtr, ytr, Xte, raw_tr, raw_te):
    qtr, _ = sim.lib_quant(raw_tr); qte, _ = sim.lib_quant(raw_te)
    # per-sample median normalisation (standard label-free quant practice)
    qtr = qtr - np.median(qtr, 1, keepdims=True); qte = qte - np.median(qte, 1, keepdims=True)
    sc = StandardScaler().fit(qtr)
    lr = LogisticRegressionCV(Cs=8, cv=min(5, int(min(np.bincount(ytr)))), max_iter=3000, penalty='l2', scoring='roc_auc')
    lr.fit(sc.transform(qtr), ytr)
    return lr.predict_proba(sc.transform(qte))[:, 1]


def pixel_lr(Xtr, ytr, Xte):
    f = lambda X: X.reshape(len(X), 32, 2, 128, 2).mean((2, 4)).reshape(len(X), -1)
    sc = StandardScaler().fit(f(Xtr))
    lr = LogisticRegressionCV(Cs=8, cv=min(5, int(min(np.bincount(ytr)))), max_iter=3000, scoring='roc_auc')
    lr.fit(sc.transform(f(Xtr)), ytr)
    return lr.predict_proba(sc.transform(f(Xte)))[:, 1]


def finetune(pre, Xtr, ytr, rng):
    m = MS1Former(); m.load_state_dict(pre.state_dict())
    head = list(m.head.parameters())
    hid = {id(p) for p in head}
    body = [p for p in m.parameters() if id(p) not in hid]
    groups = [dict(params=body, lr=1e-4), dict(params=head, lr=1e-3)]
    return train_clf(m, Xtr, ytr, rng, lr=1e-3, layer_lr=groups)
