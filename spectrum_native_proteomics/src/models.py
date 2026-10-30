"""Spectrum-native transformer models.

SpecNovo   : encoder-decoder transformer that translates MS2 peaks -> peptide.
MS1Former  : ViT-style transformer on LC-MS1 maps with masked-patch
             (MAE-style) self-supervised pretraining and a phenotype head.
"""
import math
import torch
import torch.nn as nn
import torch.nn.functional as F
import numpy as np
from chem import VOCAB, PAD, EOS, TOK_MASS, H2O, PROTON


class MassEncoder(nn.Module):
    """Sinusoidal m/z encoding with log-spaced wavelengths (1e-3 .. 1e4 Da)."""
    def __init__(self, d, lmin=1e-3, lmax=1e4):
        super().__init__()
        n = d // 2
        base = lmin / (2 * math.pi)
        scale = lmax / lmin
        self.register_buffer('freq', 1.0 / (base * scale ** (torch.arange(n) / (n - 1))))

    def forward(self, m):
        x = m.unsqueeze(-1) * self.freq
        return torch.cat([torch.sin(x), torch.cos(x)], -1)


class SpecNovo(nn.Module):
    """mass_aware=True adds (i) complementary-ion encoding of every peak and
    (ii) prefix-mass conditioning of decoder tokens."""
    def __init__(self, d=128, heads=4, enc_layers=3, dec_layers=3, ff=512, drop=0.1, max_len=16,
                 mass_aware=False, evidence=False):
        super().__init__()
        self.evidence = evidence
        if evidence:
            self.ev_mlp = nn.Sequential(nn.Linear(8, 32), nn.GELU(), nn.Linear(32, 1))
        self.d = d
        self.mass_aware = mass_aware
        self.menc = MassEncoder(d)
        self.int_proj = nn.Linear(1, d)
        self.peak_mlp = nn.Sequential(nn.Linear(2 * d if mass_aware else d, d), nn.GELU(), nn.Linear(d, d))
        if mass_aware:
            self.prefix_proj = nn.Linear(2 * d, d)
        el = nn.TransformerEncoderLayer(d, heads, ff, drop, batch_first=True, norm_first=True)
        self.encoder = nn.TransformerEncoder(el, enc_layers)
        self.charge = nn.Embedding(5, d)
        self.aa = nn.Embedding(VOCAB, d, padding_idx=PAD)
        self.pos = nn.Embedding(max_len + 2, d)
        dl = nn.TransformerDecoderLayer(d, heads, ff, drop, batch_first=True, norm_first=True)
        self.decoder = nn.TransformerDecoder(dl, dec_layers)
        self.head = nn.Linear(d, VOCAB)
        self.register_buffer('tok_mass', torch.tensor(TOK_MASS, dtype=torch.float32))

    def encode(self, mz, it, mask, pmass=None):
        e = self.menc(mz)
        if self.mass_aware:
            comp = pmass.unsqueeze(1) + 2 * PROTON - mz          # complementary ion m/z
            e = torch.cat([e, self.menc(comp)], -1)
        x = self.peak_mlp(e) + self.int_proj(it.unsqueeze(-1))
        return self.encoder(x, src_key_padding_mask=mask)

    def decode(self, mem, mask, pmass, charge, tokens, peaks=None):
        prec = self.menc(pmass) + self.charge(charge)
        x = torch.cat([prec.unsqueeze(1), self.aa(tokens)], 1)
        L = x.size(1)
        x = x + self.pos(torch.arange(L, device=x.device))
        if self.mass_aware:
            cum = torch.cat([torch.zeros_like(pmass).unsqueeze(1), self.tok_mass[tokens].cumsum(1)], 1)
            rem = pmass.unsqueeze(1) - H2O - cum
            x = x + self.prefix_proj(torch.cat([self.menc(cum + PROTON), self.menc(rem)], -1))
        causal = torch.triu(torch.full((L, L), float('-inf'), device=x.device), 1)
        h = self.decoder(x, mem, tgt_mask=causal, memory_key_padding_mask=mask)
        logits = self.head(h)
        if self.evidence and peaks is not None:
            logits = logits + self.ion_evidence(peaks, mask, pmass, tokens)
        return logits

    def ion_evidence(self, peaks, mask, pmass, tokens, tol=0.02):
        """DeepNovo-style fragment evidence: for every candidate residue a at
        every step, the observed intensity at its theoretical b/y ion m/z."""
        mz, it = peaks
        cum = torch.cat([torch.zeros_like(pmass).unsqueeze(1), self.tok_mass[tokens].cumsum(1)], 1)  # B,L
        nxt = cum.unsqueeze(-1) + self.tok_mass[3:].view(1, 1, -1)                             # B,L,A
        pm = pmass.view(-1, 1, 1)
        b = nxt + PROTON; y = pm - nxt + PROTON
        th = torch.stack([b, y, b - H2O, y - H2O, (b + PROTON) / 2, (y + PROTON) / 2,
                          cum.unsqueeze(-1).expand_as(b) + PROTON, (pm - cum.unsqueeze(-1) + PROTON).expand_as(b)], -1)  # B,L,A,8
        d = (th.unsqueeze(-1) - mz.view(mz.size(0), 1, 1, 1, -1)).abs()                     # B,L,A,8,P
        w = (d < tol).float() * it.view(it.size(0), 1, 1, 1, -1) * (~mask).float().view(mask.size(0), 1, 1, 1, -1)
        f = w.amax(-1)                                                                       # B,L,A,8
        ev = self.ev_mlp(f).squeeze(-1)                                                      # B,L,A
        out = torch.zeros(ev.size(0), ev.size(1), VOCAB, device=ev.device)
        out[..., 3:] = ev
        return out

    def forward(self, mz, it, mask, pmass, charge, tokens):
        mem = self.encode(mz, it, mask, pmass)
        return self.decode(mem, mask, pmass, charge, tokens, peaks=(mz, it))

    @torch.no_grad()
    def beam_search(self, mz, it, mask, pmass, charge, beam=5, max_len=14, tol=0.05):
        B = mz.size(0)
        mem = self.encode(mz, it, mask, pmass)
        K = beam
        mem = mem.repeat_interleave(K, 0); mmask = mask.repeat_interleave(K, 0)
        pk = (mz.repeat_interleave(K, 0), it.repeat_interleave(K, 0))
        pm = pmass.repeat_interleave(K, 0); ch = charge.repeat_interleave(K, 0)
        tokens = torch.zeros(B * K, 0, dtype=torch.long)
        score = torch.full((B, K), float('-inf')); score[:, 0] = 0
        score = score.view(-1)
        cum = torch.zeros(B * K)
        done = torch.zeros(B * K, dtype=torch.bool)
        tm = self.tok_mass
        aa_ids = torch.arange(3, VOCAB)
        for t in range(max_len + 1):
            logits = self.decode(mem, mmask, pm, ch, tokens, peaks=pk)[:, -1]
            lp = F.log_softmax(logits, -1)
            target = pm - H2O
            nxt = cum.unsqueeze(1) + tm.unsqueeze(0)
            allowed = torch.zeros_like(lp, dtype=torch.bool)
            allowed[:, 3:] = nxt[:, 3:] <= target.unsqueeze(1) + tol
            if t == max_len:
                allowed[:, 3:] = False
            allowed[:, EOS] = ((cum - target).abs() < tol) | (t == max_len)
            lp = lp.masked_fill(~allowed, float('-inf'))
            # finished beams: only PAD with zero cost
            lp[done] = float('-inf'); lp[done, PAD] = 0.0
            cand = (score.unsqueeze(1) + lp).view(B, -1)
            top, idx = cand.topk(K, -1)
            src = (idx // VOCAB) + torch.arange(B).unsqueeze(1) * K
            tok = idx % VOCAB
            src = src.view(-1); tok = tok.view(-1)
            tokens = torch.cat([tokens[src], tok.unsqueeze(1)], 1)
            score = top.view(-1)
            cum = cum[src] + tm[tok]
            done = done[src] | (tok == EOS) | (tok == PAD)
            if done.all():
                break
        # pick best beam per spectrum
        score = score.view(B, K)
        best = score.argmax(-1)
        toks = tokens.view(B, K, -1)[torch.arange(B), best]
        lens = ((toks >= 3)).sum(-1).clamp(min=1)
        conf = score[torch.arange(B), best] / (lens + 1)
        cumb = cum.view(B, K)[torch.arange(B), best]
        ok = (cumb - (pmass - H2O)).abs() < tol
        conf = torch.where(ok, conf, conf - 10.0)   # mass-inconsistent beams ranked last
        return toks, conf


# ------------------------------------------------------------------ MS1Former
class MS1Former(nn.Module):
    def __init__(self, img=(64, 256), patch=(8, 16), d=128, heads=4, layers=4, ff=384,
                 drop=0.1, dec_d=64, dec_layers=2):
        super().__init__()
        self.patch = patch
        self.nr, self.nc = img[0] // patch[0], img[1] // patch[1]
        self.N = self.nr * self.nc
        pdim = patch[0] * patch[1]
        self.embed = nn.Linear(pdim, d)
        self.pos = nn.Parameter(torch.randn(1, self.N, d) * .02)
        self.cls = nn.Parameter(torch.zeros(1, 1, d))
        el = nn.TransformerEncoderLayer(d, heads, ff, drop, batch_first=True, norm_first=True)
        self.enc = nn.TransformerEncoder(el, layers)
        self.norm = nn.LayerNorm(d)
        # MAE decoder
        self.dec_embed = nn.Linear(d, dec_d)
        self.mask_tok = nn.Parameter(torch.zeros(1, 1, dec_d))
        self.dec_pos = nn.Parameter(torch.randn(1, self.N, dec_d) * .02)
        dl = nn.TransformerEncoderLayer(dec_d, 4, dec_d * 4, 0.0, batch_first=True, norm_first=True)
        self.dec = nn.TransformerEncoder(dl, dec_layers)
        self.dec_out = nn.Linear(dec_d, pdim)
        self.head = nn.Sequential(nn.LayerNorm(d), nn.Dropout(0.2), nn.Linear(d, 2))

    def patchify(self, x):
        B = x.size(0); pr, pc = self.patch
        x = x.view(B, self.nr, pr, self.nc, pc).permute(0, 1, 3, 2, 4)
        return x.reshape(B, self.N, pr * pc)

    def features(self, x):
        p = self.embed(self.patchify(x)) + self.pos
        p = torch.cat([self.cls.expand(x.size(0), -1, -1), p], 1)
        h = self.norm(self.enc(p))
        return h

    def forward(self, x):
        h = self.features(x)
        return self.head(0.5 * (h[:, 0] + h[:, 1:].mean(1)))

    def mae_loss(self, x, ratio=0.75):
        B = x.size(0)
        tgt = self.patchify(x)
        tgt = (tgt - tgt.mean(-1, keepdim=True)) / (tgt.std(-1, keepdim=True) + 1e-4)
        p = self.embed(self.patchify(x)) + self.pos
        nkeep = int(self.N * (1 - ratio))
        noise = torch.rand(B, self.N)
        ids = noise.argsort(1); restore = ids.argsort(1)
        keep = ids[:, :nkeep]
        pk = torch.gather(p, 1, keep.unsqueeze(-1).expand(-1, -1, p.size(-1)))
        pk = torch.cat([self.cls.expand(B, -1, -1), pk], 1)
        h = self.norm(self.enc(pk))[:, 1:]
        h = self.dec_embed(h)
        full = torch.cat([h, self.mask_tok.expand(B, self.N - nkeep, -1)], 1)
        full = torch.gather(full, 1, restore.unsqueeze(-1).expand(-1, -1, full.size(-1)))
        out = self.dec_out(self.dec(full + self.dec_pos))
        m = torch.ones(B, self.N); m[:, :nkeep] = 0
        m = torch.gather(m, 1, restore)
        return (((out - tgt) ** 2).mean(-1) * m).sum() / m.sum()


class FeatureFormer(nn.Module):
    """ProFormer-style transformer over MS1 feature tokens (ID-free).

    token_j = E_id[j] + MLP(sin(m/z_j), sin(RT_j)) + value_j * w + missing_j * u
    Self-supervised objective: masked intensity modelling on unlabeled runs.
    """
    def __init__(self, n_feat, mz, rt, d=64, heads=4, layers=2, ff=128, drop=0.1):
        super().__init__()
        self.idemb = nn.Embedding(n_feat, d)
        nn.init.normal_(self.idemb.weight, std=0.02)
        menc = MassEncoder(d)
        with torch.no_grad():
            coord = torch.cat([menc(torch.tensor(mz, dtype=torch.float32)),
                               menc(torch.tensor(rt, dtype=torch.float32) * 10)], -1)
        self.register_buffer('coord', coord)
        self.cproj = nn.Linear(2 * d, d)
        self.W = nn.Parameter(torch.randn(n_feat, d))   # FT-Transformer numeric embedding
        self.miss = nn.Parameter(torch.zeros(d))
        self.mask_tok = nn.Parameter(torch.zeros(d))
        self.cls = nn.Parameter(torch.zeros(1, 1, d))
        el = nn.TransformerEncoderLayer(d, heads, ff, drop, batch_first=True, norm_first=True)
        self.enc = nn.TransformerEncoder(el, layers)
        self.norm = nn.LayerNorm(d)
        self.recon = nn.Linear(d, 1)
        self.head = nn.Sequential(nn.Dropout(0.1), nn.Linear(d, 2))
        self.tok_head = nn.Linear(d, 2)   # additive token readout (zero-init)
        nn.init.zeros_(self.tok_head.weight); nn.init.zeros_(self.tok_head.bias)

    def tokens(self, v, missing, masked=None):
        t = self.idemb.weight + 0.1 * self.cproj(self.coord)
        x = t.unsqueeze(0) + v.unsqueeze(-1) * self.W.unsqueeze(0) * (~missing).unsqueeze(-1).float() \
            + missing.unsqueeze(-1).float() * self.miss
        if masked is not None:
            x = torch.where(masked.unsqueeze(-1), t.unsqueeze(0) + self.mask_tok, x)
        return x

    def features(self, v, missing, masked=None):
        x = self.tokens(v, missing, masked)
        x = torch.cat([self.cls.expand(x.size(0), -1, -1), x], 1)
        return self.norm(self.enc(x))

    def forward(self, v, missing):
        h = self.features(v, missing)
        # CLS logit + additive per-feature evidence (sum over feature tokens)
        return self.head(h[:, 0]) + self.tok_head(h[:, 1:]).sum(1)

    def ssl_loss(self, v, missing, ratio=0.3):
        masked = (torch.rand_like(v) < ratio) & ~missing
        h = self.features(v, missing, masked)[:, 1:]
        pred = self.recon(h).squeeze(-1)
        return ((pred - v) ** 2 * masked).sum() / masked.sum().clamp(min=1)
