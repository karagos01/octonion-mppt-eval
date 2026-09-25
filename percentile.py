"""Where would the published configuration rank among the 20 000 randomly assembled ones?
Requires search.npz (produced by search.py)."""
import numpy as np
from sim import *

d = np.load('search.npz'); score = d['score']
TRAIN = ['constant', 'steps', 'ramps', 'clouds1']      # same conditions as the search (start at D=0.5)
profs = {p: profile(p) for p in TRAIN}
faith = [SIGNALS.index(x) for x in ['V', 'I', 'P', 'I/V', 'T', 'dV', 'dP', 'Vbus']]
gdef = [[0.3, -0.2, 0.1, 0.25, -0.15, 0.05, -0.3, 0.2]]   # an arbitrary thermal deformation


def ev(k, s, gc, D0=0.5):
    c = OctonionAssociator([faith], [k], [s], gc)
    return np.mean([simulate(c, 1, profs[p], D0=D0)[0] for p in TRAIN]) * 100


print(f'{"configuration":34s} {"eta":>7s}  percentile  better random ones')
for name, k, s, gc in [('as published: e0, g = 1', 0, 1, None),
                       ('e0 with an arbitrary deformation g', 0, 1, gdef),
                       ('best faithful: e5, s = -1', 5, -1, None)]:
    v = ev(k, s, gc)
    print(f'{name:34s} {v:6.2f} %  {(score*100 < v).mean()*100:7.1f} %  {(score*100 > v).sum():6d} of {len(score)}')
print(f'{"P&O (reference)":34s}  99.81 %')
print(f'\nas published, e0 with deformation, start at 12 V: {ev(0, 1, gdef, D0=0.75):.1f} %')
