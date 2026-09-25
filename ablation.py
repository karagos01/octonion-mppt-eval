"""Where does the faithful e5 get its sense of direction from? Ablation of signals and terms."""
import numpy as np
from sim import *

base = ['V', 'I', 'P', 'I/V', 'T', 'dV', 'dP', 'Vbus']
PROFS = ['constant', 'steps', 'ramps', 'clouds1', 'clouds2', 'temp_drift']
profs = {p: profile(p) for p in PROFS}

def score(mk, D0s=(0.75, 0.3), noise=True):
    return np.mean([simulate(mk(), 1, profs[p], noise=noise, D0=d)[0] for p in PROFS for d in D0s]) * 100

print('mean over 6 profiles x 2 starting points (12 V and 33.6 V); with noise / noise-free')
print(f'{"P&O":45s} {score(lambda: PerturbObserve(1)):6.2f}%  {score(lambda: PerturbObserve(1), noise=False):6.2f}%')
print(f'{"fixed duty":45s} {score(lambda: FixedDuty(1)):6.2f}%')
print(f'{"random walk":45s} {score(lambda: RandomWalk(1)):6.2f}%')
m = [SIGNALS.index(x) for x in base]
f = lambda mm: (lambda: OctonionAssociator([mm], [5], [-1]))
print(f'{"faithful e5 (all signals)":45s} {score(f(m)):6.2f}%  {score(f(m), noise=False):6.2f}%')
for i, name in enumerate(base):
    if name in ('T', 'Vbus'): continue
    mm = list(m); mm[i] = SIGNALS.index('1') if name != 'V' else SIGNALS.index('1')
    # ablation = the signal is replaced by a constant (its increments are 0)
    print(f'{"  without " + name + " (replaced by a constant)":45s} {score(f(mm)):6.2f}%  {score(f(mm), noise=False):6.2f}%')

# e5 = 2*(T*(a_I b_G - a_G b_I) + T*(a_P b_dP - a_dP b_P) - Vbus*(a_I b_P - a_P b_I) + Vbus*(a_G b_dP - a_dP b_G))
# (G = I/V; a = increment t-2 -> t-1, b = increment t-1 -> t). Each term is tried on its own as a controller.
class Term(OctonionAssociator):
    def __init__(self, fn):
        super().__init__([m], [5], [-1]); self.fn = fn
    def step(self, Vm, Im, Tm):
        S = self._state(Vm, Im, Tm); self.hist = (self.hist + [S])[-3:]
        if len(self.hist) == 3:
            X, Y, Z = self.hist; a, b = Y - X, Z - Y
            A = self.fn(a, b)
            self.dir = np.where(np.abs(A) > 1e-12, -np.sign(A), self.dir)
        self.D = np.clip(self.D + self.dir * DSTEP, D_MIN, D_MAX); return self.D
I_, P_, G_, dP_ = 1, 2, 3, 6
terms = {
    'T·(ΔI×ΔG)': lambda a, b: a[:, I_]*b[:, G_] - a[:, G_]*b[:, I_],
    'T·(ΔP×ΔdP)': lambda a, b: a[:, P_]*b[:, dP_] - a[:, dP_]*b[:, P_],
    '−Vbus·(ΔI×ΔP)': lambda a, b: -(a[:, I_]*b[:, P_] - a[:, P_]*b[:, I_]),
    'Vbus·(ΔG×ΔdP)': lambda a, b: a[:, G_]*b[:, dP_] - a[:, dP_]*b[:, G_],
}
print('\nindividual terms of the e5 expression, each used alone as the controller:')
for n_, fn in terms.items():
    print(f'{"  " + n_:45s} {score(lambda: Term(fn)):6.2f}%  {score(lambda: Term(fn), noise=False):6.2f}%')
