"""Does the controller know which way the MPP is? Start far from it, constant irradiance, with and without noise."""
import numpy as np
from sim import *

faith = [SIGNALS.index(x) for x in ['V', 'I', 'P', 'I/V', 'T', 'dV', 'dP', 'Vbus']]
cands = {
    'P&O': lambda: PerturbObserve(1),
    'random walk': lambda: RandomWalk(1),
    'faithful e1 s=-1': lambda: OctonionAssociator([faith], [1], [-1]),
    'faithful e5 s=-1': lambda: OctonionAssociator([faith], [5], [-1]),
    'faithful e0, thermal g': lambda: OctonionAssociator([faith], [0], [1], [[0.3, -0.2, 0.1, 0.25, -0.15, 0.05, -0.3, 0.2]]),
    'search #2, e3': lambda: OctonionAssociator([[SIGNALS.index(x) for x in ['dD', 'dV', 'I', 'dI', '1', 'T', 'D', 'P']]], [3], [-1]),
}
t, G, T = pr = profile('constant')
pm, vm = p_mpp(G, T)
print(f'Vmpp = {vm[-1]:.2f} V, Pmpp = {pm[-1]:.1f} W')
for noise in (True, False):
    print(f'\n--- {"with measurement noise" if noise else "NO noise"} ---')
    for D0, lab in [(0.75, 'start at 12 V'), (0.30, 'start at 33.6 V')]:
        print(f'  {lab}:')
        for name, mk in cands.items():
            eta, V = simulate(mk(), 1, pr, noise=noise, log=True, D0=D0)
            vs = V[[int(1/TS), int(5/TS), int(15/TS), -1], 0]
            print(f'    {name:22s} η={eta[0]*100:6.2f}%   V(1s,5s,15s,60s) = ' + ' '.join(f'{v:6.2f}' for v in vs))
