"""Detectors: DT-GCN (proposed), ablations and baselines."""
import numpy as np
import torch
import torch.nn as nn
import torch.nn.functional as F


def norm_adj(A):
    A = A + np.eye(len(A))
    d = 1 / np.sqrt(A.sum(1))
    return torch.tensor(d[:, None] * A * d[None, :], dtype=torch.float32)


class GCNLayer(nn.Module):
    def __init__(self, i, o):
        super().__init__()
        self.lin = nn.Linear(i, o)
        self.skip = nn.Linear(i, o, bias=False)

    def forward(self, x, A):
        return F.gelu(A @ self.lin(x) + self.skip(x))


class GCNDetector(nn.Module):
    """Graph-level detection + node-level localisation."""

    def __init__(self, A, fin, hid=64, layers=3, drop=0.1):
        super().__init__()
        self.register_buffer("A", norm_adj(A))
        dims = [fin] + [hid] * layers
        self.layers = nn.ModuleList([GCNLayer(a, b) for a, b in zip(dims[:-1], dims[1:])])
        self.drop = nn.Dropout(drop)
        self.node_head = nn.Linear(hid, 1)
        self.graph_head = nn.Sequential(nn.Linear(2 * hid, hid), nn.GELU(), nn.Linear(hid, 1))

    def forward(self, x):
        for l in self.layers:
            x = self.drop(l(x, self.A))
        node = self.node_head(x).squeeze(-1)
        g = torch.cat([x.mean(1), x.max(1).values], -1)
        return self.graph_head(g).squeeze(-1), node


class MLPDetector(nn.Module):
    def __init__(self, n, fin, hid=256, drop=0.1):
        super().__init__()
        self.net = nn.Sequential(nn.Flatten(), nn.Linear(n * fin, hid), nn.GELU(), nn.Dropout(drop),
                                 nn.Linear(hid, hid), nn.GELU(), nn.Dropout(drop), nn.Linear(hid, n + 1))

    def forward(self, x):
        o = self.net(x)
        return o[:, 0], o[:, 1:]


def train_torch(model, Xtr, ytr, Ntr, Xva, yva, Nva, seed, epochs=60, lr=2e-3, bs=128, lam=1.0):
    torch.manual_seed(seed)
    np.random.seed(seed)
    opt = torch.optim.AdamW(model.parameters(), lr=lr, weight_decay=1e-4)
    sched = torch.optim.lr_scheduler.CosineAnnealingLR(opt, epochs)
    Xtr, ytr, Ntr = map(lambda a: torch.tensor(a, dtype=torch.float32), (Xtr, ytr, Ntr))
    Xva_t = torch.tensor(Xva, dtype=torch.float32)
    pos = Ntr.mean().clamp(min=1e-3)
    pw = ((1 - pos) / pos).clamp(max=20)
    best, best_state, hist = 1e9, None, []
    for ep in range(epochs):
        model.train()
        perm = torch.randperm(len(Xtr))
        for i in range(0, len(Xtr), bs):
            b = perm[i:i + bs]
            g, nd = model(Xtr[b])
            loss = F.binary_cross_entropy_with_logits(g, ytr[b]) + lam * F.binary_cross_entropy_with_logits(
                nd, Ntr[b], pos_weight=pw)
            opt.zero_grad()
            loss.backward()
            opt.step()
        sched.step()
        model.eval()
        with torch.no_grad():
            g, nd = model(Xva_t)
            vl = (F.binary_cross_entropy_with_logits(g, torch.tensor(yva, dtype=torch.float32))
                  + lam * F.binary_cross_entropy_with_logits(nd, torch.tensor(Nva, dtype=torch.float32),
                                                             pos_weight=pw)).item()
        hist.append(vl)
        if vl < best:
            best = vl
            best_state = {k: v.clone() for k, v in model.state_dict().items()}
    model.load_state_dict(best_state)
    return model, hist


@torch.no_grad()
def predict_torch(model, X):
    model.eval()
    g, nd = model(torch.tensor(X, dtype=torch.float32))
    return torch.sigmoid(g).numpy(), torch.sigmoid(nd).numpy()
