# Octonion MPPT — independent numerical evaluation

Reproducible evaluation of the claim that the octonion associator `[X,Y,Z] = (XY)Z − X(YZ)` can
control a photovoltaic MPPT converter, as proposed in M. Mazgal, *Intensional Maxwellian Formalism
in Power Electronics* (Zenodo, DOI 10.5281/zenodo.22877912 / 10.5281/zenodo.22914238) and the
repository `Causal-Octonion-MPPT.`

The write-up is in [`PAPER.md`](PAPER.md). Short version: the complete octonion product is
implemented and verified, the missing definitions are reconstructed as faithfully as possible,
20 000 further configurations are searched, and none outperforms Perturb & Observe. The component
the original code computes (`e0`) is identically zero for any input.

> Czech version: [`README.cs.md`](README.cs.md) and [`PAPER.cs.md`](PAPER.cs.md).

## Requirements

```
python3 -m pip install -r requirements.txt   # numpy, sympy
```

No GPU needed. `search.py` uses 16 processes (`multiprocessing.Pool(16)`); lower it if you have
fewer cores.

## Reproducing the results

```
./run_all.sh          # everything, ~20 min on 16 cores
```

or individually, in this order (later scripts consume `search.npz`):

| Script | What it produces | Runtime |
|---|---|---|
| `octonion.py` | verification of the algebra (norm, alternativity, scalar part of the associator) | 3 s |
| `baseline.py` | P&O, fixed duty, random walk and all 16 faithful variants on 6 profiles | 4 min |
| `symbolic.py` | proof that `[X,Y,Z] = [X,a,b]`, and each component written out symbolically | 20 s |
| `convergence.py` | convergence from 12 V and 33.6 V, with and without measurement noise | 1 min |
| `policy.py` | noise-free decision table: does the controller know where the MPP is? | 30 s |
| `search.py` | search over 20 000 configurations, writes `search.npz` | 14 min |
| `strict.py` | strict re-test of the best twenty, signal enrichment among the successful ones | 2 min |
| `percentile.py` | where the published configuration ranks among the 20 000 random ones | 30 s |
| `ablation.py` | which signal and which single term carries the information | 3 min |
| `adaptive.py` | step size derived from the associator magnitude, gain tuning | 1 min |
| `adaptive2.py` | rise time, steady-state efficiency, partial shading over 8 noise seeds | 2 min |
| `free.py` | unrestricted step size; reference P&O with periodic full-curve scan | 3 min |
| `clouds.py` | 100 ms cloud edges and ±1000 W/m²-in-2 s ramps | 2 min |
| `why_slow.py` | share of steps in the correct direction; noise in first vs second difference | 1 min |

## Modules

- `octonion.py` — Fano-plane structure tensor, product, associator, component-wise deformation `g`,
  self-tests.
- `sim.py` — single-diode PV model, boost stage, irradiance/temperature profiles, controllers
  (P&O, variable-step P&O, P&O with full-curve scan, fixed duty, random walk, octonion associator in
  sign / variable / free modes), simulation loop vectorised over N controllers.
- `shading.py` — three substrings with bypass diodes, multi-peak P(V) curve, global MPP.

## Key numbers

| Controller | with noise | noise-free |
|---|---|---|
| P&O, fixed step | **99.35 %** | **99.45 %** |
| random walk | 50.9 % | – |
| constant voltage (duty at 0.76·Voc) | 97.6 % | – |
| as published, output `e0` | 49–94 % | 20–78 % |
| best faithful reconstruction (`e5`) | 96.6 % | 97.2 % |
| one single term of it, without octonions | 97.1 % | 96.2 % |
| best of 20 000 configurations | 98.4 % | 90.7 % |

Mean over 6 irradiance/temperature profiles × 2 starting points far from the MPP. Of 20 000
configurations, 404 beat a fixed duty cycle and none beat P&O.

The configuration actually published (output `e0`, undeformed metric, as the released code implies)
reaches 20.3 % and ranks below 16 470 of the 20 000 randomly assembled configurations — assigning the
eight quantities to components by throwing dice beats it four times out of five (`percentile.py`).

## Companion evaluations

- [`xternary-eval`](https://github.com/karagos01/xternary-eval) — the 2-bit LLM inference engine
- [`causal-trilogy-eval`](https://github.com/karagos01/causal-trilogy-eval) — the CQFT / PCTP / SOTP trilogy of September 2026
- [`openql-eval`](https://github.com/karagos01/openql-eval) — the openQL / openOL tensor matrix architecture of October 2026

## Licence

Code (all `*.py` and `run_all.sh`): MIT, see `LICENSE`.
Text of `PAPER.md` and `PAPER.cs.md`: CC BY 4.0.
