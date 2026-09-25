import numpy as np, time
from sim import *

PROFS = ['constant', 'steps', 'ramps', 'clouds1', 'clouds2', 'slow_ramp_hot']
t0 = time.time()
profs = {p: profile(p) for p in PROFS}
print(f'{"":22s}' + ''.join(f'{p[:11]:>12s}' for p in PROFS))
for name, mk in [('P&O', PerturbObserve), ('fixed duty', FixedDuty), ('random walk', RandomWalk)]:
    print(f'{name:22s}' + ''.join(f'{simulate(mk(1), 1, profs[p])[0]*100:11.2f}%' for p in PROFS))

# faithful reading of the comments in causal_associator.c:
# e0 potential=V, e1-e3 "vector fields"=I,P,I/V, e4 entropy=T, e5-e7 gradients & DC bus = dV,dP,Vbus
faith = [SIGNALS.index(s) for s in ['V', 'I', 'P', 'I/V', 'T', 'dV', 'dP', 'Vbus']]
for k in range(8):
    for s in (1, -1):
        for gname, gc in [('g=1', None), ('g thermal', np.array([[0.3, -0.2, 0.1, 0.25, -0.15, 0.05, -0.3, 0.2]]))]:
            c = OctonionAssociator([faith], [k], [s], g_coef=gc)
            eta = [simulate(c, 1, profs[p])[0] * 100 for p in PROFS]
            print(f'faithful e{k} s={s:+d} {gname:10s}' + ''.join(f'{e:11.2f}%' for e in eta))
print(f'{time.time()-t0:.1f}s')
