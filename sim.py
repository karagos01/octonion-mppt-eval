"""PV panel + boost converter + MPPT controllers, vectorised over N controllers at once."""
import numpy as np
from octonion import C

# --- PV panel: single-diode model, typical 60-cell ~250 W module ---
ISC, VOC, NS, N_ID, RS, RSH = 8.9, 37.6, 60, 1.2, 0.35, 250.0
ALPHA_I, BETA_V = 0.0005, -0.0032      # temperature coefficients of Isc (+), Voc (-) [1/K]
PSTC = 250.0
VBUS = 48.0                            # boost into a 48 V battery: Vpv = Vbus*(1-D)
TS = 0.02                              # MPPT period, 50 Hz
DSTEP = 0.004                          # duty step (~0.19 V), same for every controller
D_MIN, D_MAX = 0.0, 0.9


def pv_current(V, G, T):
    Vt = 1.380649e-23 * (T + 273.15) / 1.602176634e-19
    a = N_ID * NS * Vt
    isc_t = ISC * (1 + ALPHA_I * (T - 25))
    Iph = isc_t * G / 1000.0
    I0 = isc_t / (np.exp(VOC * (1 + BETA_V * (T - 25)) / a) - 1)
    I = np.array(Iph, dtype=float) * np.ones_like(V)
    for _ in range(20):                # Newton
        ex = np.exp(np.minimum((V + I * RS) / a, 80))
        f = Iph - I0 * (ex - 1) - (V + I * RS) / RSH - I
        df = -I0 * RS / a * ex - RS / RSH - 1
        I = I - f / df
    return np.maximum(I, 0.0)          # blocking diode


def p_mpp(G, T):
    """True MPP by golden-section search (P(V) is unimodal without shading)."""
    lo, hi = np.zeros_like(G), np.full_like(G, VOC * 1.1)
    r = (np.sqrt(5) - 1) / 2
    for _ in range(60):
        m1, m2 = hi - r * (hi - lo), lo + r * (hi - lo)
        p1, p2 = m1 * pv_current(m1, G, T), m2 * pv_current(m2, G, T)
        hi = np.where(p1 > p2, m2, hi)
        lo = np.where(p1 > p2, lo, m1)
    v = (lo + hi) / 2
    return v * pv_current(v, G, T), v


# --- irradiance [W/m2] and cell temperature profiles ---
# Explicit seeds: the random cloud realisations must not depend on the profile names,
# so that renaming a profile cannot silently change any published number.
PROFILE_SEEDS = {'clouds1': 4001462442, 'clouds2': 2005412112}
def profile(name, dur=60.0):
    t = np.arange(0, dur, TS)
    rng = np.random.default_rng(PROFILE_SEEDS.get(name, 0))
    if name == 'constant':
        G = np.full_like(t, 1000.0)
    elif name == 'steps':
        G = np.select([t < 15, t < 30, t < 45], [1000.0, 300.0, 800.0], 150.0)
    elif name == 'ramps':                       # EN 50530 style
        G = np.interp(t, [0, 5, 20, 30, 35, 45, 60], [100, 100, 1000, 1000, 300, 300, 900])
    elif name.startswith('clouds'):              # random clouds, smoothed noise
        x = np.cumsum(rng.normal(size=t.size)) * 0.06
        x = np.convolve(x - x.mean(), np.ones(25) / 25, 'same')
        G = np.clip(650 + 330 * np.tanh(x), 80, 1100)
    elif name == 'slow_ramp_hot':
        G = np.interp(t, [0, 60], [200, 1000])
    elif name == 'cloud_edges':                       # sharp cloud edges: 1000 <-> 200 W/m2 in 100 ms
        G = np.where((t // 5).astype(int) % 2 == 0, 1000.0, 200.0)
        G = np.interp(t, t, G)
        for e in np.arange(5, 60, 5):           # soften the edge to 100 ms
            m = (t >= e) & (t < e + 0.1)
            G[m] = np.interp(t[m], [e, e + 0.1], [G[np.searchsorted(t, e) - 1], G[np.searchsorted(t, e + 0.1)]])
    elif name == 'fast_ramp':                # +-1000 W/m2 within 2 s, repeatedly
        G = 550 + 450 * np.sign(np.sin(2 * np.pi * t / 8)) * np.minimum(np.abs(np.sin(2 * np.pi * t / 8)) * 4, 1)
    elif name == 'temp_drift':              # constant irradiance, ambient 0 -> 50 C (moves Vmpp)
        G = np.full_like(t, 900.0)
    else:
        raise ValueError(name)
    # cell temperature: NOCT model with thermal time constant tau = 20 s
    tamb = 35.0 if 'hot' in name else 20.0
    if name == 'temp_drift':
        tamb = np.interp(t, [0, 60], [0.0, 50.0])
    Tt = tamb + G / 800 * 25
    T = np.empty_like(Tt); T[0] = Tt[0]
    for k in range(1, t.size):
        T[k] = T[k - 1] + (Tt[k] - T[k - 1]) * TS / 20.0
    return t, G, T


# --- controllers: step(Vm, Im, Tm) -> new duty cycle D (vector of N) ---
class PerturbObserve:
    def __init__(self, n):
        self.D = np.full(n, 0.5); self.dir = np.ones(n); self.P0 = np.zeros(n)

    def step(self, Vm, Im, Tm):
        P = Vm * Im
        self.dir = np.where(P < self.P0, -self.dir, self.dir)
        self.P0 = P
        self.D = np.clip(self.D + self.dir * DSTEP, D_MIN, D_MAX)
        return self.D


class VariableStepPO:
    """classical variable-step P&O: step = N*|dP/dV|, clipped to 1x...10x DSTEP"""
    def __init__(self, n, N, cap=10.0):
        self.D = np.full(n, 0.5); self.dir = np.ones(n); self.P0 = np.zeros(n); self.V0 = np.zeros(n)
        self.N = np.asarray(N, float); self.cap = cap

    def step(self, Vm, Im, Tm):
        P = Vm * Im
        self.dir = np.where(P < self.P0, -self.dir, self.dir)
        slope = np.abs(P - self.P0) / np.maximum(np.abs(Vm - self.V0), 0.05)
        self.P0, self.V0 = P, Vm
        step = self.N * slope
        step = np.clip(step, DSTEP, self.cap * DSTEP) if self.cap else np.maximum(step, DSTEP)
        self.D = np.clip(self.D + self.dir * step, D_MIN, D_MAX)
        return self.D


class ScanPO:
    """the proper answer to shading: P&O + a full-curve sweep every `period` s, jump to the highest peak"""
    def __init__(self, n, period=10.0, pts=40):
        self.D = np.full(n, 0.5); self.dir = np.ones(n); self.P0 = np.zeros(n)
        self.grid = np.linspace(0.05, 0.85, pts); self.i = 0; self.period = period; self.pts = pts
        self.best = np.zeros(n); self.bestD = np.full(n, 0.5)

    def step(self, Vm, Im, Tm):
        P = Vm * Im
        k = int(round(self.period / TS))
        ph = self.i % k
        if ph < self.pts:                      # sweep phase
            if ph > 0:
                better = P > self.best
                self.best = np.where(better, P, self.best)
                self.bestD = np.where(better, self.grid[ph - 1], self.bestD)
            if ph == 0:
                self.best = np.zeros_like(P)
            self.D = np.full_like(self.D, self.grid[ph])
        elif ph == self.pts:                   # jump to the best point found
            self.D = self.bestD.copy()
        else:                                  # ordinary P&O
            self.dir = np.where(P < self.P0, -self.dir, self.dir)
            self.D = np.clip(self.D + self.dir * DSTEP, D_MIN, D_MAX)
        self.P0 = P; self.i += 1
        return self.D


class FixedDuty:
    """dumb reference: constant-voltage method, duty fixed at 0.76*Voc regardless of the start"""
    own_duty = True          # simulate() must not overwrite this controller's duty with D0

    def __init__(self, n):
        self.D = np.full(n, 1 - 0.76 * VOC / VBUS)

    def step(self, Vm, Im, Tm):
        return self.D


class RandomWalk:
    """random walk with the same step: what a controller that knows nothing achieves"""
    def __init__(self, n, seed=1):
        self.D = np.full(n, 0.5); self.rng = np.random.default_rng(seed)

    def step(self, Vm, Im, Tm):
        self.D = np.clip(self.D + DSTEP * self.rng.choice([-1, 1], self.D.size), D_MIN, D_MAX)
        return self.D


# signals that can be assigned to octonion components (dimensionless, order 1)
SIGNALS = ['1', 'V', 'I', 'P', 'dV', 'dI', 'dP', 'T', 'Vbus', 'D', 'I/V', 'dD']
CFLAT = C.reshape(64, 8)


def omul(a, b, g):
    return (a[:, :, None] * b[:, None, :]).reshape(len(a), 64) @ CFLAT * g


class OctonionAssociator:
    """X, Y, Z = three consecutive states ("causal memory"), output = component k of the associator.

    mapping: (N,8) indices into SIGNALS for e0..e7
    k:       (N,) which associator component drives the duty cycle, s: (N,) sign
    g_coef:  (N,8) metric deformation g_i = 1 + g_coef_i * (T-25)/50 ("degradation history")
    mode:    'sign' -> dD = s*DSTEP*sign(A_k) (holds the previous direction when A_k ~ 0)
             'var'  -> magnitude from |A_k|, 'free' -> no limit at all
    """
    def __init__(self, mapping, k, s, g_coef=None, mode='sign', gain=None, hold=True):
        n = len(mapping)
        self.map, self.k, self.s = np.asarray(mapping), np.asarray(k), np.asarray(s, float)
        self.gc = np.zeros((n, 8)) if g_coef is None else np.asarray(g_coef)
        self.mode, self.gain, self.hold = mode, gain, hold
        self.D = np.full(n, 0.5); self.dir = np.ones(n)
        self.prev = None; self.hist = []
        self.A_log = []
        self.agree = []; self.P0 = None; self.last_dir = np.ones(n)

    def _state(self, Vm, Im, Tm):
        P = Vm * Im
        if self.prev is None:
            self.prev = (Vm, Im, P, self.D.copy(), np.zeros_like(Vm))
        V0, I0, P0, D0, _ = self.prev
        dD = self.D - D0
        sig = np.stack([np.ones_like(Vm), Vm / VOC, Im / ISC, P / PSTC,
                        (Vm - V0) / (DSTEP * VBUS), (Im - I0) / 0.05, (P - P0) / 1.0,
                        (Tm - 25) / 50, np.full_like(Vm, VBUS / VOC), self.D,
                        (Im / ISC) / np.maximum(Vm / VOC, 0.05), dD / DSTEP], axis=1)
        self.prev = (Vm, Im, P, self.D.copy(), dD)
        return np.take_along_axis(sig, self.map, axis=1)

    def step(self, Vm, Im, Tm):
        S = self._state(Vm, Im, Tm)
        self.hist = (self.hist + [S])[-3:]
        if len(self.hist) < 3:
            A = np.zeros(len(S))
        else:
            X, Y, Z = self.hist
            g = 1 + self.gc * ((Tm - 25) / 50)[:, None]
            assoc = omul(omul(X, Y, g), Z, g) - omul(X, omul(Y, Z, g), g)
            A = assoc[np.arange(len(S)), self.k]
        self.A_log.append(A)
        P = Vm * Im
        if self.P0 is not None:   # what would P&O do in this situation?
            po = np.where(P < self.P0, -self.last_dir, self.last_dir)
            self.po_dir = po
        self.P0 = P
        if self.mode == 'sign':
            nz = np.abs(A) > 1e-12
            new = self.s * np.sign(A)
            self.dir = np.where(nz, new, self.dir if self.hold else 0.0)
            dD = self.dir * DSTEP
        elif self.mode == 'free':   # no limit on the step size whatsoever
            dD = self.s * self.gain * A
        else:   # 'var': direction from the sign, step size from |A| (the memory)
            nz = np.abs(A) > 1e-12
            self.dir = np.where(nz, self.s * np.sign(A), self.dir)
            dD = self.dir * np.clip(self.gain * np.abs(A), DSTEP, 10 * DSTEP)
        if hasattr(self, 'po_dir'):
            moved = dD != 0
            self.agree.append(np.where(moved, np.sign(dD) == self.po_dir, np.nan))
        self.last_dir = np.where(dD != 0, np.sign(dD), self.last_dir)
        self.D = np.clip(self.D + dD, D_MIN, D_MAX)
        return self.D

    def agreement(self):
        return np.nanmean(np.array(self.agree), axis=0)


def simulate(ctrl, n, prof, noise=True, log=False, D0=0.5, pvf=None, pmf=None, seed=12345):
    """Returns the tracking efficiency eta = E / E_mpp for each of the N controllers."""
    t, G, T = prof
    pvf = pvf or pv_current
    pm, _ = (pmf or p_mpp)(G, T)
    rng = np.random.default_rng(seed)
    D = np.full(n, float(D0))
    if getattr(ctrl, 'own_duty', False):
        D = ctrl.D.copy()    # controller with its own operating point (e.g. constant voltage)
    else:
        ctrl.D = D.copy()
    E = np.zeros(n); Vlog = []
    for i in range(t.size):
        V = VBUS * (1 - D)
        I = pvf(V, G[i], T[i])
        E += V * I
        if log: Vlog.append(V.copy())
        Vm = V + (rng.normal(0, 0.05, n) if noise else 0)
        Im = I + (rng.normal(0, 0.02, n) if noise else 0)
        D = ctrl.step(Vm, Im, np.full(n, T[i]))
    eta = E / pm.sum()
    return (eta, np.array(Vlog)) if log else eta
