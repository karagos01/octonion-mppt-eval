"""Complete octonion product (Fano plane) + metric deformation as in causal_associator.c.

All functions are vectorised: inputs have shape (..., 8).
"""
import numpy as np

# Fano plane (Baez, "The Octonions"): e_i e_j = e_k for these oriented triples
FANO = [(1, 2, 4), (2, 3, 5), (3, 4, 6), (4, 5, 7), (5, 6, 1), (6, 7, 2), (7, 1, 3)]


def _structure_tensor():
    C = np.zeros((8, 8, 8))
    C[0, 0, 0] = 1.0
    for i in range(1, 8):
        C[0, i, i] = C[i, 0, i] = 1.0   # e0 is the unit
        C[i, i, 0] = -1.0               # e_i^2 = -1
    for i, j, k in FANO:
        for a, b, c in ((i, j, k), (j, k, i), (k, i, j)):
            C[a, b, c] = 1.0            # e_a e_b = e_c
            C[b, a, c] = -1.0           # anticommutativity of the imaginary units
    return C


C = _structure_tensor()


def mul(a, b, g=None):
    """(ab)_k = g_k * sum_ij a_i b_j C_ijk. g=None -> standard octonions.

    The component-wise deformation g_k is the only sensible generalisation of line 24 in
    causal_associator.c, where e0 is multiplied by history.g[0] (the author left the rest unwritten).
    """
    r = np.einsum('...i,...j,ijk->...k', a, b, C)
    return r if g is None else r * g


def assoc(x, y, z, g=None):
    """[X, Y, Z] = (XY)Z - X(YZ)"""
    return mul(mul(x, y, g), z, g) - mul(x, mul(y, z, g), g)


def conj(a):
    return a * np.array([1, -1, -1, -1, -1, -1, -1, -1.0])


if __name__ == '__main__':
    rng = np.random.default_rng(0)
    x, y, z = rng.normal(size=(3, 100000, 8))
    n = lambda v: np.linalg.norm(v, axis=-1)
    E = np.eye(8)

    checks = {
        'e_i*e_i = -1 (i=1..7)': max(abs(mul(E[i], E[i])[0] + 1) for i in range(1, 8)),
        'normed algebra |xy| = |x||y|': np.max(abs(n(mul(x, y)) - n(x) * n(y)) / (n(x) * n(y))),
        'alternativity [x,x,y] = 0': np.max(abs(assoc(x, x, y))),
        'alternativity [x,y,y] = 0': np.max(abs(assoc(x, y, y))),
        'flexibility [x,y,x] = 0': np.max(abs(assoc(x, y, x))),
        'scalar part of the associator = 0': np.max(abs(assoc(x, y, z)[:, 0])),
        'power associativity (xx)x = x(xx)': np.max(abs(mul(mul(x, x), x) - mul(x, mul(x, x)))),
    }
    for k, v in checks.items():
        print(f'{k:40s} max error {v:.2e}  {"OK" if v < 1e-9 else "!!"}')
    print(f'{"associator nonzero in general":40s} mean |[x,y,z]| = {n(assoc(x, y, z)).mean():.3f}')

    # his code: g rescales the result. For constant g the associator just scales by g^2, properties hold.
    g = np.full(8, 0.9)
    print(f'{"g = 0.9 const: associator scalar":40s} max {np.max(abs(assoc(x, y, z, g)[:, 0])):.2e}')
    g = 1 + 0.05 * rng.normal(size=8)
    print(f'{"g non-const (+-5%): assoc. scalar":40s} mean {np.mean(abs(assoc(x, y, z, g)[:, 0])):.3f}'
          f'  vs imag. parts {np.mean(abs(assoc(x, y, z, g)[:, 1:])):.3f}')
    print(f'{"g non-const: [x,x,y]":40s} mean {np.mean(n(assoc(x, x, y, g))):.3f}')
