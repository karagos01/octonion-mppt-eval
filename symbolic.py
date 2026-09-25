"""What does associator component k actually compute for the 'faithful' signal assignment?
By trilinearity and alternativity, [X,Y,Z] with Y = X+a, Z = Y+b equals [X, a, b]
(X = state, a = first increment, b = next increment)."""
import sympy as sp
from octonion import C

names = ['V', 'I', 'P', 'G', 'T', 'dV', 'dP', 'Vb']   # G = I/V (conductance)
X = sp.symbols([n for n in names])
a = sp.symbols(['a_' + n for n in names])      # increment between t-2 and t-1
b = sp.symbols(['b_' + n for n in names])      # increment between t-1 and t
zero = {a[4]: 0, b[4]: 0, a[7]: 0, b[7]: 0}    # T moves slowly, Vbus is constant

def mul(u, v):
    return [sp.expand(sum(u[i] * v[j] * C[i, j, k] for i in range(8) for j in range(8) if C[i, j, k]))
            for k in range(8)]

Y = [X[i] + a[i] for i in range(8)]
Z = [Y[i] + b[i] for i in range(8)]
A = [sp.expand(p - q) for p, q in zip(mul(mul(X, Y), Z), mul(X, mul(Y, Z)))]
check = [sp.expand(p - q) for p, q in zip(mul(mul(X, a), b), mul(X, mul(a, b)))]
print('[X,Y,Z] == [X,a,b]:', all(sp.simplify(A[k] - check[k]) == 0 for k in range(8)))
for k in range(8):
    e = sp.factor(A[k].subs(zero))
    print(f'e{k}: {e}')
