"""Brute force: random assignments of signals to e0..e7, choice of output component, sign and deformation g."""
import numpy as np, time, sys
from multiprocessing import Pool
from sim import *

TRAIN = ['constant', 'steps', 'ramps', 'clouds1']
TEST = ['clouds2', 'slow_ramp_hot', 'temp_drift']
N = int(sys.argv[1]) if len(sys.argv) > 1 else 20000

rng = np.random.default_rng(7)
mapping = np.array([rng.permutation(len(SIGNALS))[:8] for _ in range(N)])
k = rng.integers(0, 8, N)
s = rng.choice([-1, 1], N)
gc = np.where(rng.random((N, 1)) < 0.5, 0.0, rng.normal(0, 0.3, (N, 8)))

def run(args):
    prof, idx = args
    c = OctonionAssociator(mapping[idx], k[idx], s[idx], gc[idx])
    return simulate(c, len(idx), profile(prof)), c.agreement()

def evaluate(profs, idx):
    chunks = np.array_split(idx, 4)
    with Pool(16) as p:
        res = p.map(run, [(pr, ch) for pr in profs for ch in chunks])
    eta = np.array([np.concatenate([r[0] for r in res[i*4:(i+1)*4]]) for i in range(len(profs))])
    agr = np.array([np.concatenate([r[1] for r in res[i*4:(i+1)*4]]) for i in range(len(profs))])
    return eta, agr

if __name__ == '__main__':
    t0 = time.time()
    ref = {p: simulate(PerturbObserve(1), 1, profile(p))[0] for p in TRAIN + TEST}
    fix = {p: simulate(FixedDuty(1), 1, profile(p))[0] for p in TRAIN + TEST}
    eta, agr = evaluate(TRAIN, np.arange(N))
    score = eta.mean(0)
    print(f'{N} candidates, training {time.time()-t0:.0f}s')
    print(f'P&O training mean {np.mean([ref[p] for p in TRAIN])*100:.2f} %, fixed duty {np.mean([fix[p] for p in TRAIN])*100:.2f} %')
    for q in (50, 90, 99, 99.9, 100):
        print(f'  percentile {q:5}: {np.percentile(score, q)*100:.2f} %')
    print(f'  candidates better than fixed duty: {(score > np.mean([fix[p] for p in TRAIN])).sum()}')
    print(f'  candidates better than P&O: {(score > np.mean([ref[p] for p in TRAIN])).sum()}')

    top = np.argsort(-score)[:20]
    teta, tagr = evaluate(TEST, top)
    print('\nTOP 20 (training -> test on unseen profiles), agreement of decisions with P&O')
    print(f'{"P&O":>60s} {np.mean([ref[p] for p in TEST])*100:6.2f} %')
    print(f'{"fixed duty":>60s} {np.mean([fix[p] for p in TEST])*100:6.2f} %')
    for j, i in enumerate(top):
        m = ' '.join(f'e{c}={SIGNALS[mapping[i][c]]}' for c in range(8))
        print(f'{m:48s} out=e{k[i]} s={s[i]:+d} g={"def" if gc[i].any() else "1  "} '
              f'train {score[i]*100:6.2f}%  test {teta[:, j].mean()*100:6.2f}% '
              f'(min {teta[:, j].min()*100:6.2f})  agree w. P&O {np.nanmean(tagr[:, j])*100:5.1f}%')
    np.savez('search.npz', mapping=mapping, k=k, s=s, gc=gc, score=score, agr=agr, top=top, teta=teta, tagr=tagr)
    print(f'done in {time.time()-t0:.0f}s')
