import numpy as np
from multiprocessing import Pool
from sim import *
from analyze import symbolic

d = np.load('search.npz')
mapping, k, s, gc, score = d['mapping'], d['k'], d['s'], d['gc'], d['score']
PROFS = ['constant', 'steps', 'ramps', 'clouds1', 'clouds2', 'temp_drift']

# 1) which signals do the successful candidates use (better than fixed duty in training)?
FIXED = np.mean([simulate(FixedDuty(1), 1, profile(p))[0] for p in PROFS[:4]])
good = score > FIXED
print(f'fixed-duty baseline on the training profiles: {FIXED*100:.2f} %')
print(f'successful candidates: {good.sum()} of {len(score)}')
print('signal | share among the successful ones  (67 % would be chance)')
for i, sg in enumerate(SIGNALS):
    has = (mapping == i).any(1)
    print(f'  {sg:5s} {has[good].mean()*100:5.1f} %')
has_dp = (mapping == SIGNALS.index('dP')).any(1)
print(f'successful without dP: {(good & ~has_dp).sum()}')

# 2) strict test of the top 20: start at 12 V and 33.6 V, with and without noise
top = d['top']
def run(a):
    prof, D0, noise = a
    c = OctonionAssociator(mapping[top], k[top], s[top], gc[top])
    return simulate(c, len(top), profile(prof), noise=noise, D0=D0)
jobs = [(p, D0, nz) for nz in (True, False) for p in PROFS for D0 in (0.75, 0.3)]
with Pool(16) as pool:
    res = np.array(pool.map(run, jobs))
half = len(PROFS) * 2
ns, nn = res[:half].mean(0) * 100, res[half:].mean(0) * 100
po = [np.mean([simulate(PerturbObserve(1), 1, profile(p), noise=nz, D0=D0)[0] for p in PROFS for D0 in (0.75, 0.3)]) * 100 for nz in (True, False)]
print(f'\nstrict test (6 profiles x starts at 12 V and 33.6 V):  P&O {po[0]:.2f} % with noise, {po[1]:.2f} % noise-free')
for j, i in enumerate(top):
    print(f'  #{j+1:2d} out=e{k[i]}  with noise {ns[j]:6.2f} %   noise-free {nn[j]:6.2f} %')

# 3) what does the best candidate compute?
i = top[0]
print('\n#1:', ' '.join(f'e{c}={SIGNALS[mapping[i][c]]}' for c in range(8)), f'out=e{k[i]} s={s[i]:+d}')
print('   expression:', symbolic(mapping[i], k[i]))
