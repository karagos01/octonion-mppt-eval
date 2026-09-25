"""Behaviour during a cloud passage: recovery time after an edge and voltage error on a fast ramp."""
import numpy as np
from sim import *

faith = [SIGNALS.index(x) for x in ['V', 'I', 'P', 'I/V', 'T', 'dV', 'dP', 'Vbus']]
d = np.load('search.npz'); i1 = d['top'][0]
cands = {'P&O fixed step': lambda: PerturbObserve(1),
         'P&O variable step': lambda: VariableStepPO(1, 1.9e-3),
         'faithful e5': lambda: OctonionAssociator([faith], [5], [-1]),
         'best from the search': lambda: OctonionAssociator([d['mapping'][i1]], [d['k'][i1]], [d['s'][i1]]),
         'best from search (step from memory)': lambda: OctonionAssociator([d['mapping'][i1]], [d['k'][i1]], [d['s'][i1]], mode='var', gain=0.00316)}

for prof in ('cloud_edges', 'fast_ramp'):
    t, G, T = pr = profile(prof)
    pm, vm = p_mpp(G, T)
    print(f'\n=== profile {prof} ===')
    hdr = 'recovery after edge [s]' if prof == 'cloud_edges' else 'voltage error from Vmpp [V]'
    print(f'{"":33s} {"η":>7s}  {hdr}')
    for name, mk in cands.items():
        etas, mets = [], []
        for sd in range(8):
            eta, V = simulate(mk(), 1, pr, log=True, D0=0.5, seed=sd)
            V = V[:, 0]; P = V * pv_current(V, G, T); etas.append(eta[0])
            if prof == 'cloud_edges':
                rec = []
                for e in np.arange(5, 60, 5):
                    i0 = int(e / TS); w = slice(i0, i0 + int(4.9 / TS))
                    ok = np.where(P[w] / pm[w] > 0.98)[0]
                    rec.append(ok[0] * TS if len(ok) else 4.9)
                mets.append(np.mean(rec))
            else:
                mets.append(np.mean(np.abs(V - vm)))
        print(f'{name:33s} {np.mean(etas)*100:6.2f}%  {np.median(mets):8.2f}')
