#!/usr/bin/env python3
"""Do the component labels in MPPT Sec. 2.1 carry any information?

Sec. 2.1 of the preprint (v1.1) reads:

    "Let O be a state octonion mapping dimensions e_0 ... e_7 to physical
     potentials (e.g., scalar voltage, vector current, thermodynamic entropy,
     and electromotive gradients)."

That is the only place the physics enters the octonion layer, so the whole
proposal rests on it.  Three things are measured here.

1. How many relabellings of e_1..e_7 leave the multiplication table bitwise
   identical.  If a relabelling is an automorphism, then two different physical
   assignments define the same algebra, and the difference between them is not
   observable in anything the algebra computes -- the labels are decoration.
   The basis automorphism group is measured directly rather than quoted, by
   checking which of the 5040 permutations preserve the Fano incidence structure
   and then which sign vectors complete each of them.

2. Whether the stated assignment is dimensionally possible at all.  An octonion
   is an algebra over R: all eight components are added to each other, so they
   must share one physical dimension.  Reading Sec. 2.1 literally gives each
   component a different unit, and every product then equates incompatible ones.

3. What e_0 is.  It is not a slot for a physical quantity; it is the
   multiplicative identity.

The counts are invariant across valid orientations of the Fano plane, which the
script confirms on the orientation used throughout this repository and on the
one the same author later published as SOTP Eq. (3).
"""
import itertools

import numpy as np

from octonion import C as C_REPO

# SOTP Eq. (3), verbatim, for the cross-check.  A different orientation of the
# same seven lines from the one in octonion.py, and also a valid octonion algebra.
SOTP = [(1, 2, 3), (1, 4, 5), (2, 4, 6), (3, 4, 7), (1, 7, 6), (2, 5, 7), (3, 6, 5)]


def from_triads(triads):
    """Structure tensor from seven oriented triads, e_0 as the unit."""
    C = np.zeros((8, 8, 8))
    C[0, 0, 0] = 1.0
    for i in range(1, 8):
        C[0, i, i] = C[i, 0, i] = 1.0
        C[i, i, 0] = -1.0
    for i, j, k in triads:
        for a, b, c in ((i, j, k), (j, k, i), (k, i, j)):
            C[a, b, c] = 1.0
            C[b, a, c] = -1.0
    return C


def sign_index(C):
    """(sign, index) form: e_i e_j = sign[i,j] * e_index[i,j]."""
    idx = np.argmax(np.abs(C), axis=2)
    sgn = np.array([[C[i, j, idx[i, j]] for j in range(8)] for i in range(8)])
    return sgn, idx


def automorphisms(C):
    """Relabellings e_i -> s_i * e_{p(i)} (i >= 1, e_0 fixed) that preserve C.

    Returns (pure permutations, index-preserving permutations, full group order).
    """
    sgn, idx = sign_index(C)
    pure, struct, order = [], [], 0
    for perm in itertools.permutations(range(1, 8)):
        p = [0] + list(perm)
        # the index condition does not involve the signs: p must carry the
        # product e_i e_j onto the product e_p(i) e_p(j)
        if any(idx[p[i], p[j]] != p[idx[i, j]] for i in range(8) for j in range(8)):
            continue
        struct.append(perm)
        for bits in range(128):
            s = [1] + [1 if (bits >> (m - 1)) & 1 == 0 else -1 for m in range(1, 8)]
            if all(sgn[p[i], p[j]] == sgn[i, j] * s[i] * s[j] * s[idx[i, j]]
                   for i in range(8) for j in range(8)):
                order += 1
                if bits == 0:
                    pure.append(perm)
    return pure, struct, order


print("1. Relabellings of e_1..e_7 that leave the multiplication table identical")
print(f"   {'table':<34s} {'pure':>6s} {'up to sign':>11s} {'group order':>12s}"
      f" {'distinct':>9s}")
for name, C in (("this repository (octonion.py)", C_REPO),
                ("SOTP Eq. (3), same author", from_triads(SOTP))):
    pure, struct, order = automorphisms(C)
    print(f"   {name:<34s} {len(pure):>6d} {len(struct):>11d} {order:>12d}"
          f" {5040 // len(struct):>9d}")

pure, struct, order = automorphisms(C_REPO)
print(f"\n   Of the 7! = {5040} ways to assign seven physical quantities to")
print(f"   e_1..e_7, each one gives an algebra identical to {len(struct)-1} of the others,")
print(f"   so the 5040 assignments collapse to {5040 // len(struct)} distinct algebras."
      f"  {len(pure)} of the")
print("   relabellings need no sign change at all.")
print(f"   The full basis automorphism group has order {order} = 8 x {len(struct)};")
print("   over the reals the automorphism group of O is the 14-dimensional")
print("   exceptional Lie group G_2, which acts transitively on basic triples, so")
print("   no invariant of the algebra distinguishes e_1 from e_4 in the first")
print("   place.  Naming one of them 'thermodynamic entropy' adds no constraint")
print("   that any computation can see.")

# ---- 2. dimensional consistency of the stated assignment ------------------
print("\n2. Units, taking Sec. 2.1 literally")
# SI base dimensions (kg, m, s, A, K)
DIM = {"V":   np.array([1, 2, -3, -1, 0]),
       "A":   np.array([0, 0, 0, 1, 0]),
       "J/K": np.array([1, 2, -2, 0, -1]),
       "V/m": np.array([1, 1, -3, -1, 0]),
       "1":   np.array([0, 0, 0, 0, 0])}
unit = {0: "V", 1: "A", 2: "A", 3: "A", 4: "J/K", 5: "V/m", 6: "V/m", 7: "V/m"}
print("   e_0 = V (scalar voltage), e_1..e_3 = A (vector current),")
print("   e_4 = J/K (thermodynamic entropy), e_5..e_7 = V/m (electromotive gradient)")

sgn, idx = sign_index(C_REPO)
products = [(i, j, int(idx[i, j])) for i in range(8) for j in range(8)]
mixed = [(i, j, k) for i, j, k in products if i and j and i != j]
bad, shown = [], 0
for i, j, k in mixed:
    if not np.array_equal(DIM[unit[i]] + DIM[unit[j]], DIM[unit[k]]):
        bad.append((i, j, k))
        if shown < 5:
            shown += 1
            s = "-" if sgn[i, j] < 0 else "+"
            print(f"     e_{i} e_{j} = {s}e_{k}"
                  f"   requires  [{unit[i]}][{unit[j]}] = [{unit[k]}]")
print(f"   products of two distinct imaginary units            {len(mixed)} of 56")
print(f"   of those, dimensionally inconsistent as labelled    {len(bad)}")
sums = {tuple(DIM[unit[i]]) for i in range(8)}
print(f"   distinct dimensions among the eight components      {len(sums)}")

# is ANY assignment of dimensions to e_0..e_7 consistent?  One linear system per
# base dimension: d_i + d_j - d_k = 0 for every product e_i e_j = +-e_k.
A = np.array([[(1 if t == i else 0) + (1 if t == j else 0) - (1 if t == k else 0)
               for t in range(8)] for i, j, k in products], dtype=float)
rank = np.linalg.matrix_rank(A)
print(f"\n   constraints d_i + d_j = d_k from the full table      {A.shape[0]}")
print(f"   rank of that system in the 8 unknowns d_0..d_7       {rank} of 8")
print(f"   dimension of the solution space                      {8 - rank}")
print("   -> the only assignment of physical dimensions to the eight components")
print("      that the multiplication table admits is d_0 = ... = d_7 = 0: every")
print("      component dimensionless.  It follows from two entries of the table")
print("      alone -- e_0 e_i = e_i forces d_0 = 0, and e_i e_i = -e_0 then")
print("      forces d_i = 0.  An octonion is an algebra over R, its components")
print("      are added to each other, and there is no octonion over R whose")
print("      components carry different units.  Sec. 2.1 does not describe a")
print("      richer object; it describes an expression that cannot be evaluated.")

# ---- 3. what e_0 is -------------------------------------------------------
e0 = np.zeros(8); e0[0] = 1.0
prod = np.einsum("i,j,ijk->k", e0, e0, C_REPO)
print(f"\n3. e_0 is not a slot: e_0 e_0 = {prod.astype(int).tolist()}, i.e. e_0 itself.")
print("   e_0 is the multiplicative identity, so calling it 'scalar voltage'")
print("   asserts V^2 = V.  The mapping has eight names and seven places, and")
print("   the seven places are interchangeable.")
