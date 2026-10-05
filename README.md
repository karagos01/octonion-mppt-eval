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
| `labels.py` | whether the Sec. 2.1 component labels carry information: basis automorphisms, dimensional consistency | 1 s |
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
- `census.sh` — **the only script that needs the network.** Downloads all ten of the author's
  Zenodo deposits, converts them with `pdftotext` and counts the keywords behind the Maxwell
  attribution claim (PAPER.md §5).

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

Dice do as well because the slots are interchangeable. Of the 7! = 5040 ways to map seven physical
quantities onto `e1…e7`, **168** leave the multiplication table bitwise identical, so each assignment
is algebraically indistinguishable from 167 others and the 5040 collapse to **30** distinct algebras;
the full basis automorphism group has order **1344** (`labels.py`, same four numbers on this
repository's orientation and on the one later published as SOTP Eq. 3). Taking the units of §2.1
literally, **42 of 56** products of distinct imaginary units equate incompatible dimensions, and the
constraint system `d_i + d_j = d_k` over the whole table has rank **8 of 8** — so the only
dimensionally consistent octonion over ℝ is the one whose eight components are all dimensionless.
`e0` is not a slot at all: it is the multiplicative identity, so calling it "scalar voltage" asserts
V² = V.

The *Maxwellian* in the title is uncited, and `census.sh` shows the attribution tracks the algebra
rather than the history: Maxwell is named in the four deposits that pair him with octonions
(2, 1, 4 and 5 times) and **zero** times in CQFT, PCTP and SOTP, where the algebra is quaternionic.
"Treatise", "1873", "1865" and "Heaviside" occur **zero** times in all ten deposits. Quaternions are
in the 1873 *Treatise* at §§618–619 — two pages of Hamiltonian operator notation in which Maxwell
never multiplies two quaternions — and octonions appear nowhere in Maxwell at all.

## Companion evaluations

- [`xternary-eval`](https://github.com/karagos01/xternary-eval) — the 2-bit LLM inference engine
- [`causal-trilogy-eval`](https://github.com/karagos01/causal-trilogy-eval) — the CQFT / PCTP / SOTP trilogy of September 2026
- [`openql-eval`](https://github.com/karagos01/openql-eval) — the openQL / openOL tensor matrix architecture of October 2026
- [`cymatic-eval`](https://github.com/karagos01/cymatic-eval) — the cymatic stomatal stimulation and pulsed light proposal of October 2026

## Licence

Code (all `*.py` and `run_all.sh`): MIT, see `LICENSE`.
Text of `PAPER.md` and `PAPER.cs.md`: CC BY 4.0.
