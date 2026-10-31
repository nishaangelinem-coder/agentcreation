import sys, time, json, pickle, argparse
import numpy as np, torch, torch.nn.functional as F
from chem import TOK, PAD, EOS
from models import SpecNovo

ap = argparse.ArgumentParser()
ap.add_argument('--n_train', type=int, default=250000)
ap.add_argument('--steps', type=int, default=3500)
ap.add_argument('--tag', default='main')
ap.add_argument('--threads', type=int, default=4)
ap.add_argument('--save_at', default='1500')
ap.add_argument('--mass_aware', type=int, default=0)
ap.add_argument('--evidence', type=int, default=0)
ap.add_argument('--init', default='')
ap.add_argument('--lr', type=float, default=5e-4)
ap.add_argument('--corpus', default='data/ms2_corpus.pkl')
ap.add_argument('--warmup', type=int, default=200)
a = ap.parse_args()
torch.set_num_threads(a.threads); torch.manual_seed(0)
D = pickle.load(open(a.corpus, 'rb'))
train = D['train'][:a.n_train]; val = D['val'][:1000]
P = 80

def collate(items):
    B = len(items)
    mz = np.zeros((B, P), np.float32); it = np.zeros((B, P), np.float32)
    mask = np.ones((B, P), bool); tok = np.zeros((B, 16), np.int64)
    pm = np.zeros(B, np.float32); ch = np.zeros(B, np.int64)
    for i, (m, s, p, z, seq) in enumerate(items):
        n = len(m); mz[i, :n] = m; it[i, :n] = s; mask[i, :n] = False
        t = [TOK[c] for c in seq] + [EOS]; tok[i, :len(t)] = t
        pm[i] = p; ch[i] = z
    return [torch.from_numpy(x) for x in (mz, it, mask, pm, ch, tok)]

model = SpecNovo(mass_aware=bool(a.mass_aware), evidence=bool(a.evidence))
if a.init:
    print(model.load_state_dict(torch.load(a.init), strict=False))
opt = torch.optim.AdamW(model.parameters(), lr=a.lr, weight_decay=1e-4)
warm = a.warmup
sched = torch.optim.lr_scheduler.LambdaLR(opt, lambda s: min(1, (s + 1) / warm) * 0.5 * (1 + np.cos(np.pi * min(s, a.steps) / a.steps)))
rng = np.random.default_rng(0)
vb = collate(val)
log = []
save_at = {int(x) for x in a.save_at.split(',') if x}
t0 = time.time(); step = 0
while step < a.steps:
    for idx in np.array_split(rng.permutation(len(train)), len(train) // 128):
        model.train()
        mz, it, mask, pm, ch, tok = collate([train[i] for i in idx])
        inp = tok[:, :-1].clone(); inp[inp == EOS] = PAD
        logits = model(mz, it, mask, pm, ch, inp)
        loss = F.cross_entropy(logits.reshape(-1, logits.size(-1)), tok.reshape(-1),
                               ignore_index=PAD, label_smoothing=0.05)
        opt.zero_grad(); loss.backward()
        torch.nn.utils.clip_grad_norm_(model.parameters(), 1.0)
        opt.step(); sched.step(); step += 1
        if step % 100 == 0 or step == 1:
            model.eval()
            with torch.no_grad():
                mz, it, mask, pm, ch, tok = vb
                inp = tok[:, :-1].clone(); inp[inp == EOS] = PAD
                lg = model(mz, it, mask, pm, ch, inp)
                vl = F.cross_entropy(lg.reshape(-1, lg.size(-1)), tok.reshape(-1), ignore_index=PAD).item()
                acc = ((lg.argmax(-1) == tok) & (tok != PAD)).sum().item() / (tok != PAD).sum().item()
            log.append(dict(step=step, train_loss=loss.item(), val_loss=vl, val_tok_acc=acc,
                            seen=step * 128, time=time.time() - t0))
            print(log[-1], flush=True)
            json.dump(log, open(f'results/trainlog_{a.tag}.json', 'w'))
        if step in save_at:
            torch.save(model.state_dict(), f'results/specnovo_{a.tag}_s{step}.pt')
        if step >= a.steps: break
torch.save(model.state_dict(), f'results/specnovo_{a.tag}.pt')
print('done', time.time() - t0)
