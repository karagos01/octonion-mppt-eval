"""Why is the octonion controller slower? Share of steps in the correct direction + strength of the useful term."""
import numpy as np
from sim import *
from octonion import assoc

faith = [SIGNALS.index(x) for x in ['V', 'I', 'P', 'I/V', 'T', 'dV', 'dP', 'Vbus']]
d = np.load('search.npz'); i1 = d['top'][0]
cands = {'P&O': lambda: PerturbObserve(1),
         'faithful e5': lambda: OctonionAssociator([faith], [5], [-1]),
         'search #1': lambda: OctonionAssociator([d['mapping'][i1]], [d['k'][i1]], [d['s'][i1]])}
t, G, T = pr = profile('constant'); pm, vm = p_mpp(G, T)

print('share of steps in the correct direction (towards the MPP), start at 12 V, median of 8 noise seeds:')
for name, mk in cands.items():
    far, near = [], []
    for sd in range(8):
        _, V = simulate(mk(), 1, pr, log=True, D0=0.75, seed=sd)
        V = V[:, 0]; step = np.sign(np.diff(V)); want = np.sign(vm[:-1] - V[:-1])
        ok = step == want
        m = np.abs(V[:-1] - vm[:-1]) > 2          # far from the MPP
        far.append(ok[m].mean()); near.append(ok[~m].mean())
    print(f'  {name:12s} far from MPP {np.median(far)*100:5.1f} %   near MPP {np.median(near)*100:5.1f} %')

# how large is the useful term inside the associator compared with the rest?
c = OctonionAssociator([faith], [5], [-1])
A_tot, A_use = [], []
prev = []
Vp = Ip = Pp = None
_, V = simulate(OctonionAssociator([faith], [5], [-1]), 1, pr, log=True, D0=0.75)
V = V[:, 0]; I = pv_current(V, G, T); P = V * I
rng = np.random.default_rng(12345)
Vm_ = V + rng.normal(0, 0.05, V.size); Im_ = I + rng.normal(0, 0.02, V.size); Pm_ = Vm_ * Im_
sig = np.stack([np.ones_like(V), Vm_/VOC, Im_/ISC, Pm_/PSTC,
                np.r_[0, np.diff(Vm_)]/(DSTEP*VBUS), np.r_[0, np.diff(Im_)]/0.05, np.r_[0, np.diff(Pm_)]/1.0,
                (T-25)/50, np.full_like(V, VBUS/VOC), 1-V/VBUS, (Im_/ISC)/np.maximum(Vm_/VOC, .05),
                np.r_[0, np.diff(1-V/VBUS)]/DSTEP], 1)[:, faith]
X, Y, Z = sig[:-2], sig[1:-1], sig[2:]
Atot = assoc(X, Y, Z)[:, 5]
a, b = Y - X, Z - Y
use = 2 * (a[:, 3]*b[:, 6] - a[:, 6]*b[:, 3]) * (VBUS/VOC)     # the term Vbus*(dG x d(dP))
print(f'\nin associator e5: |useful term| / |whole| = {np.median(np.abs(use)/np.maximum(np.abs(Atot),1e-12)):.2f}')
print(f'the sign of the whole matches the sign of the useful term in {np.mean(np.sign(use)==np.sign(Atot))*100:.0f} % of steps')
# how much noise the second difference adds
print(f'\nnoise in dP (1st difference):    {np.std(np.diff(Pm_)-np.diff(P)):.3f} W')
print(f'noise in d(dP) (2nd difference): {np.std(np.diff(Pm_,2)-np.diff(P,2)):.3f} W')
