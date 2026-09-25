"""Noise-free decision table: at voltage V and after the last three moves (d0,d1,d2), where does the
controller send the voltage? Compared with P&O (reverse when the power dropped)."""
import numpy as np, itertools
from sim import *
from octonion import assoc

G_, T_ = 1000.0, profile('constant')[2][-1]
_, vm = p_mpp(np.array([G_]), np.array([T_])); vm = vm[0]
faith = [SIGNALS.index(x) for x in ['V', 'I', 'P', 'I/V', 'T', 'dV', 'dP', 'Vbus']]
dv = DSTEP * VBUS

def sig(V, Vp, Ip, Pp, D, dD):
    I = pv_current(np.array([V]), G_, T_)[0]; P = V * I
    return np.array([1, V/VOC, I/ISC, P/PSTC, (V-Vp)/dv, (I-Ip)/0.05, (P-Pp)/1.0, (T_-25)/50, VBUS/VOC, D,
                     (I/ISC)/max(V/VOC, .05), dD/DSTEP]), I, P

def decide(k, s, Vend, moves):
    """moves = voltage directions (+1 = voltage rising) of the last 3 steps; returns the new direction"""
    Vs = [Vend - dv * sum(moves[i:]) for i in range(4)] + [Vend]   # V_{t-3}..V_t
    Vs = Vs[1:] if len(Vs) > 4 else Vs
    Vs = [Vend - dv*(moves[0]+moves[1]+moves[2]), Vend - dv*(moves[1]+moves[2]), Vend - dv*moves[2], Vend]
    Vp, Ip, Pp, Dp = Vs[0], pv_current(np.array([Vs[0]]), G_, T_)[0], None, 1 - Vs[0]/VBUS
    Pp = Vp * Ip; dDp = 0.0; states = []
    for V in Vs[1:]:
        D = 1 - V / VBUS
        st, I, P = sig(V, Vp, Ip, Pp, D, D - Dp)
        states.append(st[faith]); Vp, Ip, Pp, Dp = V, I, P, D
    A = assoc(*[np.array(x) for x in states])[k]
    dD = s * np.sign(A)
    return -int(dD)            # duty up = voltage down

print(f'Vmpp = {vm:.2f} V.  Table: new voltage direction (+ up, - down); correct is up below the MPP, down above it')
Vgrid = [12, 18, 22, 25, 26, 27, 28, 30, 33, 35]
pats = [(1,1,1), (-1,1,1), (1,-1,1), (-1,-1,1), (1,1,-1), (-1,1,-1), (1,-1,-1), (-1,-1,-1)]
for k, s, name in [(5, -1, 'faithful e5'), (1, -1, 'faithful e1')]:
    print(f'\n{name}   moves(d0,d1,d2) \\ V: ' + ' '.join(f'{v:4d}' for v in Vgrid))
    good = tot = 0
    for p in pats:
        row = []
        for V in Vgrid:
            d = decide(k, s, V, p)
            ok = (d > 0) == (V < vm)
            good += ok; tot += 1
            row.append(('+' if d > 0 else '-' if d < 0 else '0') + (' ' if ok else '✗'))
        print(f'  {str(p):30s}      ' + ' '.join(f'{r:>4s}' for r in row))
    print(f'  correct direction: {good}/{tot}')
