"""Analysis of a single candidate: symbolic form of the output + behaviour of the trajectory."""
import sys, numpy as np, sympy as sp
from sim import *
from octonion import C

def symbolic(mapping, k):
    const = {'1', 'Vbus', 'T'}          # constants / slow quantities: increments = 0
    X = [sp.Symbol(SIGNALS[m]) for m in mapping]
    a = [0 if SIGNALS[m] in const else sp.Symbol('a_' + SIGNALS[m]) for m in mapping]
    b = [0 if SIGNALS[m] in const else sp.Symbol('b_' + SIGNALS[m]) for m in mapping]
    mul = lambda u, v: [sum(u[i] * v[j] * int(C[i, j, kk]) for i in range(8) for j in range(8) if C[i, j, kk]) for kk in range(8)]
    A = [sp.expand(p - q) for p, q in zip(mul(mul(X, a), b), mul(X, mul(a, b)))]
    subs = {sp.Symbol('1'): 1}
    return sp.factor(A[k].subs(subs))

def behaviour(mapping, k, s, gc=None, prof='steps'):
    c = OctonionAssociator([mapping], [k], [s], None if gc is None else [gc])
    t, G, T = pr = profile(prof)
    eta, V = simulate(c, 1, pr, log=True)
    _, vm = p_mpp(G, T)
    D = 1 - V[:, 0] / VBUS
    dirs = np.sign(np.diff(D))
    print(f'  η={eta[0]*100:.2f} %, agree w. P&O {c.agreement()[0]*100:.1f} %')
    for ti in (1, 5, 14, 20, 29, 35, 44, 50, 59):
        i = int(ti / TS)
        print(f'  t={ti:2d}s  V={V[i,0]:6.2f}  Vmpp={vm[i]:6.2f}  last 16 steps: '
              + ''.join('+' if d > 0 else '-' if d < 0 else '0' for d in -dirs[i-16:i]))

if __name__ == '__main__':
    cands = {
        'faithful e1': ([SIGNALS.index(x) for x in ['V', 'I', 'P', 'I/V', 'T', 'dV', 'dP', 'Vbus']], 1, -1),
        'faithful e5': ([SIGNALS.index(x) for x in ['V', 'I', 'P', 'I/V', 'T', 'dV', 'dP', 'Vbus']], 5, -1),
        'search #2': ([SIGNALS.index(x) for x in ['dD', 'dV', 'I', 'dI', '1', 'T', 'D', 'P']], 3, -1),
    }
    for name, (m, k, s) in cands.items():
        print(f'\n== {name}: out e{k}, s={s:+d}')
        print('  expression:', symbolic(m, k))
        behaviour(m, k, s)
    print('\n== P&O for comparison')
    c = PerturbObserve(1); t, G, T = pr = profile('steps')
    eta, V = simulate(c, 1, pr, log=True); _, vm = p_mpp(G, T)
    dirs = np.sign(np.diff(1 - V[:, 0] / VBUS))
    for ti in (1, 5, 14, 20, 29):
        i = int(ti / TS)
        print(f'  t={ti:2d}s  V={V[i,0]:6.2f}  Vmpp={vm[i]:6.2f}  last 16 steps: '
              + ''.join('+' if d > 0 else '-' if d < 0 else '0' for d in -dirs[i-16:i]))
