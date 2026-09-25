"""Variable step: octonion controller (step size from |associator|) vs variable-step P&O."""
import numpy as np, time
from multiprocessing import Pool
from sim import *
from shading import pv_shaded, p_mpp_shaded, shading_profile

d = np.load('search.npz')
faith = [SIGNALS.index(x) for x in ['V', 'I', 'P', 'I/V', 'T', 'dV', 'dP', 'Vbus']]
C_map = [faith] * 16 + list(d['mapping'][d['top']])
C_k = list(range(8)) * 2 + list(d['k'][d['top']])
C_s = [1] * 8 + [-1] * 8 + list(d['s'][d['top']])
C_name = [f'faithful e{k} s={s:+d}' for s in (1, -1) for k in range(8)] + [f'search #{j+1}' for j in range(20)]
GAINS = np.logspace(-3, 4, 15)
NS_ = np.logspace(-5, -1, 15)
nc, ng = len(C_map), len(GAINS)
M = np.repeat(np.array(C_map), ng, 0); K = np.repeat(C_k, ng); S = np.repeat(C_s, ng); GA = np.tile(GAINS, nc)

TRAIN = ['constant', 'steps', 'ramps', 'clouds1']
TEST = ['clouds2', 'temp_drift', 'shade']

def job(a):
    prof, D0, which, idx = a
    if prof == 'shade':
        pr, kw = shading_profile(), dict(pvf=pv_shaded, pmf=p_mpp_shaded)
    else:
        pr, kw = profile(prof), {}
    if which == 'okt':
        c = OctonionAssociator(M[idx], K[idx], S[idx], mode='var', gain=GA[idx][:, None][:, 0])
        n = len(idx)
    elif which == 'vspo':
        c = VariableStepPO(len(idx), NS_[idx]); n = len(idx)
    else:
        c = PerturbObserve(1); n = 1
    eta, V = simulate(c, n, pr, D0=D0, log=True, **kw)
    return eta, V

def run(profs, which, idx, D0s=(0.75, 0.3)):
    jobs = [(p, D0, which, ch) for p in profs for D0 in D0s for ch in np.array_split(idx, 4)]
    with Pool(16) as pool:
        res = pool.map(job, jobs)
    eta = np.array([np.concatenate([r[0] for r in res[i*4:(i+1)*4]]) for i in range(len(profs) * len(D0s))])
    V = [np.concatenate([r[1] for r in res[i*4:(i+1)*4]], 1) for i in range(len(profs) * len(D0s))]
    return eta, V

def settle(V, prof='constant'):
    t, G, T = profile(prof); pm, _ = p_mpp(G, T)
    P = V * pv_current(V, G[:, None], T[:, None])
    bad = P / pm[:, None] < 0.98
    last = np.array([np.where(b)[0][-1] + 1 if b.any() else 0 for b in bad.T])
    return np.where(last < len(t), last * TS, np.inf)

if __name__ == '__main__':
    t0 = time.time()
    # 1) tuning on the training profiles
    eo, _ = run(TRAIN, 'okt', np.arange(len(M)))
    ev, _ = run(TRAIN, 'vspo', np.arange(len(NS_)))
    so = eo.mean(0).reshape(nc, ng); sv = ev.mean(0)
    best_g = so.argmax(1); bv = sv.argmax()
    order = np.argsort(-so.max(1))[:6]
    sel = np.array([c * ng + best_g[c] for c in order])
    print(f'tuning done in {time.time()-t0:.0f}s')
    print(f'variable-step P&O: best N={NS_[bv]:.1e}, training {sv[bv]*100:.2f} %')
    for c in order:
        print(f'  {C_name[c]:18s} best gain {GAINS[best_g[c]]:8.3g}  training {so[c].max()*100:.2f} %')

    # 2) test: unseen profiles incl. shading + speed
    rows = [('P&O fixed step', run(TEST + ['constant'], 'po', np.arange(1))),
            ('P&O variable step', run(TEST + ['constant'], 'vspo', np.array([bv])))]
    eo2, Vo2 = run(TEST + ['constant'], 'okt', sel)
    for j, c in enumerate(order):
        rows.append((C_name[c] + ' (var)', (eo2[:, [j]], [v[:, [j]] for v in Vo2])))
    print(f'\n{"":24s} {"clouds2":>8s} {"temp.drift":>10s} {"SHADE":>8s}   time to 98 % (start 12 V / 33.6 V)')
    for name, (e, V) in rows:
        e = e[:, 0] * 100
        m = lambda i: (e[2*i] + e[2*i+1]) / 2
        st = [settle(V[6][:, 0][:, None])[0], settle(V[7][:, 0][:, None])[0]]
        print(f'{name:24s} {m(0):7.2f}% {m(1):9.2f}% {m(2):7.2f}%   ' + ' / '.join(f'{x:.2f}s' if np.isfinite(x) else 'never' for x in st))
    print(f'done in {time.time()-t0:.0f}s')
