"""Unrestricted step: both the octonion controller and P&O may jump anywhere."""
import numpy as np
from multiprocessing import Pool
from sim import *
from shading import pv_shaded, p_mpp_shaded, shading_profile
import adaptive as A

GAINS = np.logspace(-5, 5, 21)
TRAIN = ['constant', 'steps', 'ramps', 'clouds1']
TEST = ['clouds2', 'temp_drift']
nc = len(A.C_map); ng = len(GAINS)
M = np.repeat(np.array(A.C_map), ng, 0); K = np.repeat(A.C_k, ng); S = np.repeat(A.C_s, ng); GA = np.tile(GAINS, nc)
NS_ = np.logspace(-5, 0, 21)

def job(a):
    what, kind, idx, D0, seed = a
    pr, kw = (shading_profile(), dict(pvf=pv_shaded, pmf=p_mpp_shaded)) if what == 'shade' else (profile(what), {})
    if kind == 'okt_free':
        c = OctonionAssociator(M[idx], K[idx], S[idx], mode='free', gain=GA[idx]); n = len(idx)
    elif kind == 'vspo_free':
        c = VariableStepPO(len(idx), NS_[idx], cap=None); n = len(idx)
    elif kind == 'po':
        c = PerturbObserve(1); n = 1
    elif kind == 'scan':
        c = ScanPO(1); n = 1
    eta, V = simulate(c, n, pr, D0=D0, log=True, seed=seed, **kw)
    return eta, V

def batch(profs, kind, idx, D0s=(0.75, 0.3), seed=12345):
    jobs = [(p, kind, ch, D0, seed) for p in profs for D0 in D0s for ch in np.array_split(idx, 4)]
    with Pool(16) as pool:
        r = pool.map(job, jobs)
    eta = np.array([np.concatenate([x[0] for x in r[i*4:(i+1)*4]]) for i in range(len(profs)*len(D0s))])
    return eta, [np.concatenate([x[1] for x in r[i*4:(i+1)*4]], 1) for i in range(len(profs)*len(D0s))]

if __name__ == '__main__':
    eo, _ = batch(TRAIN, 'okt_free', np.arange(len(M)))
    ev, _ = batch(TRAIN, 'vspo_free', np.arange(len(NS_)))
    so = eo.mean(0).reshape(nc, ng); sv = ev.mean(0)
    best = np.argsort(-so.max(1))[:5]; bv = sv.argmax()
    print('tuning (training, mean of 4 profiles x 2 starting points):')
    print(f'  P&O unrestricted step: N={NS_[bv]:.1e} -> {sv[bv]*100:.2f} %')
    for c in best:
        print(f'  {A.C_name[c]:18s} gain {GAINS[so[c].argmax()]:9.3g} -> {so[c].max()*100:.2f} %')

    ts, Gs, Ts = shading_profile(); pms, vms = p_mpp_shaded(Gs, Ts)
    def hits(V):
        h = []
        for seg in range(6):
            i = slice(seg*500+400, seg*500+500)
            h.append(np.median(np.abs(V[i] - vms[i])) < 3)
        return np.mean(h)

    print(f'\n{"":30s} {"clouds2":>8s} {"temp.drift":>11s} {"SHADE":>8s} {"global peak":>13s}')
    rows = [('P&O fixed step', 'po', np.arange(1)), ('P&O unrestricted step', 'vspo_free', np.array([bv])),
            ('P&O + full-curve sweep', 'scan', np.arange(1))]
    rows += [(A.C_name[c] + ' (free)', 'okt_free', np.array([c*ng + so[c].argmax()])) for c in best]
    for name, kind, idx in rows:
        e, _ = batch(TEST, kind, idx)
        es, Vs = zip(*[batch(['shade'], kind, idx, seed=sd) for sd in range(4)])
        estin = np.mean([x[:, 0].mean() for x in es]) * 100
        hh = np.mean([hits(V[0][:, 0]) for _, V in [(None, V) for V in Vs]] + [hits(V[1][:, 0]) for _, V in [(None, V) for V in Vs]])
        m = lambda i: (e[2*i, 0] + e[2*i+1, 0]) / 2 * 100
        print(f'{name:30s} {m(0):7.2f}% {m(1):10.2f}% {estin:7.2f}% {hh*100:12.0f} %')
