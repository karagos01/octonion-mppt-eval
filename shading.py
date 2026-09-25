"""Module of 3 substrings of 20 cells, each with a bypass diode -> multi-peak P(V) under shading."""
import numpy as np
from sim import *

def _vsub(I, G, T):
    Vt = 1.380649e-23 * (T + 273.15) / 1.602176634e-19
    a = N_ID * (NS / 3) * Vt
    isc_t = ISC * (1 + ALPHA_I * (T - 25))
    Iph = isc_t * G / 1000.0
    I0 = isc_t / (np.exp(VOC / 3 * (1 + BETA_V * (T - 25)) / a) - 1)
    arg = (Iph - I) / I0 + 1
    V = np.where(arg > 0, a * np.log(np.maximum(arg, 1e-300)) - I * RS / 3, -np.inf)
    return np.maximum(V, -0.5)            # bypass diode

def panel_v(I, G3, T):
    return sum(_vsub(I, G3[..., j], T) for j in range(3))

def pv_shaded(V, G3, T):
    """current at voltage V (bisection; V(I) is decreasing)"""
    V = np.asarray(V, float)
    lo, hi = np.zeros_like(V), np.full_like(V, ISC * 1.1 * G3.max() / 1000)
    for _ in range(40):
        mid = (lo + hi) / 2
        above = panel_v(mid, G3, T) > V
        lo, hi = np.where(above, mid, lo), np.where(above, hi, mid)
    return np.where(panel_v(np.zeros_like(V), G3, T) > V, (lo + hi) / 2, 0.0)

def p_mpp_shaded(G, T):
    """global MPP by a dense grid"""
    Vg = np.linspace(0.5, VBUS, 800)
    pm, vm = np.zeros(len(G)), np.zeros(len(G))
    cache = {}
    for i in range(len(G)):
        key = (tuple(np.round(G[i])), round(T[i], 1))
        if key not in cache:
            P = Vg * pv_shaded(Vg, G[i][None, :].repeat(len(Vg), 0), T[i])
            cache[key] = (P.max(), Vg[P.argmax()])
        pm[i], vm[i] = cache[key]
    return pm, vm

def shading_profile(name='shade'):
    t = np.arange(0, 60.0, TS)
    pats = [(1000, 1000, 1000), (1000, 1000, 300), (1000, 250, 250), (900, 900, 500), (1000, 1000, 1000), (400, 1000, 1000)]
    G = np.array([pats[min(int(ti // 10), 5)] for ti in t], float)
    T = np.full(t.size, 45.0)
    return t, G, T

if __name__ == '__main__':
    t, G, T = shading_profile()
    Vg = np.linspace(0.5, 40, 400)
    for p in sorted(set(map(tuple, G)), key=lambda x: -sum(x)):
        P = Vg * pv_shaded(Vg, np.array(p)[None, :].repeat(len(Vg), 0), 45.0)
        # local peaks
        pk = [i for i in range(1, len(P) - 1) if P[i] > P[i-1] and P[i] >= P[i+1] and P[i] > 5]
        print(f'{str(p):20s} peaks: ' + ', '.join(f'{P[i]:.0f} W @ {Vg[i]:.1f} V' for i in pk))
