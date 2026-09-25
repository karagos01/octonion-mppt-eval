"""Rise time + steady-state efficiency; shading over 8 noise seeds and which peak the controller ends on."""
import numpy as np
from multiprocessing import Pool
from sim import *
from shading import pv_shaded, p_mpp_shaded, shading_profile, panel_v
import adaptive as A

# selection: the best from tuning (taken from the output of adaptive.py)
CANDS = {'P&O fixed step': ('po', None), 'P&O variable step': ('vspo', 1.9e-3)}
for name, g in [('search #7', 0.0316), ('search #13', 0.01), ('search #19', 0.0316), ('search #1', 0.00316)]:
    c = A.C_name.index(name)
    CANDS[name + ' (var)'] = ('okt', (A.C_map[c], A.C_k[c], A.C_s[c], g))
CANDS['search #1 (fixed step)'] = ('okt_sign', (A.C_map[A.C_name.index('search #1')], A.C_k[A.C_name.index('search #1')], A.C_s[A.C_name.index('search #1')], None))

def make(kind, p, n):
    if kind == 'po': return PerturbObserve(n)
    if kind == 'vspo': return VariableStepPO(n, np.full(n, p))
    m, k, s, g = p
    if kind == 'okt': return OctonionAssociator([m] * n, [k] * n, [s] * n, mode='var', gain=g)
    return OctonionAssociator([m] * n, [k] * n, [s] * n)

def job(a):
    name, what, D0, seed = a
    kind, p = CANDS[name]
    if what == 'shade':
        pr = shading_profile(); eta, V = simulate(make(kind, p, 1), 1, pr, D0=D0, log=True, pvf=pv_shaded, pmf=p_mpp_shaded, seed=seed)
    else:
        pr = profile('constant'); eta, V = simulate(make(kind, p, 1), 1, pr, D0=D0, log=True, seed=seed)
    return name, what, D0, seed, eta[0], V[:, 0]

if __name__ == '__main__':
    seeds = range(8)
    jobs = [(n, w, D0, sd) for n in CANDS for w in ('konst', 'shade') for D0 in (0.75, 0.3) for sd in seeds]
    with Pool(16) as pool:
        res = pool.map(job, jobs)
    t, G, T = profile('constant'); pm, _ = p_mpp(G, T)
    ts, Gs, Ts = shading_profile(); pms, vms = p_mpp_shaded(Gs, Ts)
    print(f'{"":26s} {"rise 12V":>9s} {"rise 33.6V":>11s} {"steady":>9s} {"SHADE":>8s} {"+-(seeds)":>9s} {"global peak":>15s}')
    for name in CANDS:
        r = [x for x in res if x[0] == name]
        rise = {}
        for D0 in (0.75, 0.3):
            v = [x[5] for x in r if x[1] == 'konst' and x[2] == D0]
            rt = []
            for V in v:
                P = V * pv_current(V, G, T); ok = np.where(P / pm >= 0.98)[0]
                rt.append(ok[0] * TS if len(ok) else np.inf)
            rise[D0] = np.median(rt)
        steady = np.mean([np.mean((x[5] * pv_current(x[5], G, T))[-1500:] / pm[-1500:]) for x in r if x[1] == 'konst']) * 100
        es = np.array([x[4] for x in r if x[1] == 'shade']) * 100
        # on the global peak = during the last 2 s of each 10 s segment the voltage is within 3 V of the global MPP
        hits = []
        for x in r:
            if x[1] != 'shade': continue
            for seg in range(6):
                i = slice(seg * 500 + 400, seg * 500 + 500)
                hits.append(np.median(np.abs(x[5][i] - vms[i])) < 3)
        print(f'{name:26s} {rise[0.75]:8.2f}s {rise[0.3]:10.2f}s {steady:8.2f}% {es.mean():7.2f}% {es.std():8.2f}% {np.mean(hits)*100:13.0f} %')
