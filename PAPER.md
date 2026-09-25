# Does the Octonion Associator Help MPPT? An Independent Numerical Evaluation

**Author:** karagos01 · **Status:** Working paper / technical comment · **Date:** 2026-09-25 ·
**Czech version:** [PAPER.cs.md](PAPER.cs.md)

**Comment on:** M. Mazgal, *Intensional Maxwellian Formalism in Power Electronics: Causal Octonion
Control of Distributed Multi-Phase Converters*, v1.0 (DOI 10.5281/zenodo.22877912) and v1.1
(DOI 10.5281/zenodo.22914238), and the accompanying repository `Causal-Octonion-MPPT.`

---

## Abstract

A recent preprint proposes controlling a photovoltaic MPPT converter through the associator of a
non-associative octonion algebra, `[X, Y, Z] = (XY)Z − X(YZ)`, used as a projection of a "hidden
causal memory" onto the PWM duty cycle. The proposal specifies neither the multiplication rules of
the claimed deformed algebra, nor the physical quantities entering the eight octonion components,
nor the mapping from the associator to the duty cycle, and it reports no measurements. This comment
supplies the missing definitions in the most faithful way we could reconstruct them, implements the
complete octonion product, and evaluates the resulting controllers in a photovoltaic simulation
against Perturb & Observe (P&O). We searched 20 000 assignments of measured signals to octonion
components, output components, signs and metric deformations. No configuration outperformed P&O.
The best ones reduce, after symbolic expansion, to the P&O decision rule "change of power × change
of direction". The single component that the published source code actually computes, `e0`, is
identically zero for any input; taken with the undeformed metric that the released code implies, the
published configuration ranks below 16 470 of the 20 000 randomly assembled ones. We also show that the associator of three consecutive states equals
`[X, a, b]`, a trilinear function of the state and the last two increments, which explains both the
measured 3× slower convergence and the failure under fast irradiance ramps: the useful signal is
obtained one difference later, with 74 % more measurement noise. All code is provided.

## 1. What is claimed and what is published

The preprint proposes a two-layer architecture: a "hidden causal memory" layer in octonion algebra
(O) accumulating thermal drift, magnetic saturation and material fatigue, and an observable tensor
layer producing the duty cycle, bridged by the octonion associator. Claimed consequences include
pre-emptive hardware protection, resolution of magnetic saturation and thermal degradation, and
wire-free phase-locking of a converter swarm over the shared DC bus.

The published artifacts contain:

| Needed to evaluate the proposal | Present? |
|---|---|
| multiplication rules of the deformed algebra (all 8 components) | no |
| update rule for the deformation ("metric") `g` | no |
| physical quantity, unit and scaling of each component `e0…e7` | labels only |
| definition of the associator arguments X, Y, Z | no |
| mapping from associator to duty cycle | v1.0: one dimensionally inconsistent equation; v1.1: none |
| MPP search procedure | no |
| measurements, figures, comparison with any standard method | no |

In the source file `causal_associator.c`, `OctonionMultiply` computes only component `e[0]`;
components `e[1]…e[7]` are replaced by the comment `// ...` and are left uninitialised, so the
returned associator is undefined behaviour rather than a number.

## 2. Method

**Algebra.** We implement the full octonion product from the Fano plane
(triples (1,2,4), (2,3,5), (3,4,6), (4,5,7), (5,6,1), (6,7,2), (7,1,3)) as a structure tensor,
and verify on 10⁵ random triples that the implementation satisfies, to machine precision
(max error ≤ 2.8·10⁻¹⁴): `e_i² = −1`, multiplicativity of the norm `|xy| = |x||y|`, alternativity
`[x,x,y] = [x,y,y] = 0`, flexibility `[x,y,x] = 0`, power associativity, and that the scalar
component of the associator vanishes identically. Following the published code, the deformation is
applied component-wise to each product, `(ab)_k → g_k·(ab)_k`.

**Plant.** Single-diode model of a 250 W, 60-cell module (Isc 8.9 A, Voc 37.6 V, Rs 0.35 Ω,
Rsh 250 Ω), boost stage into a 48 V battery (`Vpv = Vbus(1−D)`), quasi-static per control step,
MPPT update rate 50 Hz, duty step 0.004 (≈0.19 V), measurement noise σ_V = 50 mV, σ_I = 20 mA.
True MPP per sample is obtained by golden-section search, so efficiencies are referenced to the
actual maximum, not to a fixed set point. Six irradiance/temperature profiles are used (constant,
steps, EN 50530-style ramps, two random cloud profiles, and a temperature drift from 0 °C to 50 °C
ambient), plus a partial-shading model of three substrings with bypass diodes, and a hard test with
100 ms cloud edges and ±1000 W/m² ramps within 2 s.

**Controllers.** The octonion controller takes X, Y, Z as three consecutive state octonions
("causal memory"), and drives the duty cycle from one component `A_k` of the associator, in three
variants: sign only (fixed step, the fairest comparison to P&O), magnitude from `|A_k|` (1× to 10×
the base step, gain tuned), and unrestricted step. Baselines: P&O with fixed step, P&O with
variable step ∝ |ΔP/ΔV| (gain tuned identically), P&O with periodic full-curve scan, fixed duty,
and a random walk with the same step size.

**Search.** 20 000 random candidates over: assignment of 8 distinct signals out of 12 (constant,
V, I, P, ΔV, ΔI, ΔP, T, Vbus, D, I/V, ΔD) to components `e0…e7`; output component; sign; and
either `g = 1` or a random thermal deformation. Candidates were scored on four training profiles
and the best twenty re-tested on unseen profiles.

## 3. Results

### 3.1 Tracking efficiency

Strict test: mean over 6 profiles × 2 starting points far from the MPP (12 V and 33.6 V).

| Controller | with noise | noise-free |
|---|---|---|
| P&O, fixed step | **99.35 %** | **99.45 %** |
| random walk, same step | 50.9 % | – |
| constant voltage (duty fixed at 0.76·Voc) | 97.6 % | – |
| as published, output `e0` | 49–94 % | 20–78 % |
| most faithful reconstruction, output `e5` | 96.6 % | 97.2 % |
| one single term of it, no octonions | 97.1 % | 96.2 % |
| best of 20 000 searched configurations | 98.4 % | 90.7 % |

Of the 20 000 configurations, 404 beat a fixed duty cycle and **none** beat P&O
(median 48.8 %, 99th percentile 98.9 %, maximum 99.78 % against 99.81 % for P&O, evaluated at a
starting point near the MPP). Thirteen of the best twenty lose more than 10 percentage points (up to 61) when
the measurement noise is removed; several then run to open-circuit voltage, where the output power is
zero. Their tracking is therefore partly noise-driven.

Note the constant-voltage baseline: a controller that measures nothing and simply holds the duty
cycle at 0.76·Voc reaches 97.6 %, i.e. more than the most faithful reconstruction of the proposal
(96.6 %) and more than all but 404 of the 20 000 searched configurations. On a single module whose MPP
voltage moves by about 4 V over the whole temperature range, that baseline is hard to beat, which is
why we report P&O as the reference rather than the naive one.

**Sensitivity to the random profiles.** At the easy operating point the margin is small enough to
depend on the particular cloud realisation. With a different random seed for the two cloud profiles,
one configuration out of 20 000 edges past P&O by 0.02 percentage points when started near the MPP
(99.80 % against 99.79 %), while in the strict test from 12 V and 33.6 V the same configuration stays
0.7 points behind P&O with noise and 8.5 points behind without it. The comparison that is stable is
therefore the strict test, which is what we report; the random seeds of the profiles are pinned in the
code so that no published number depends on them.

The configuration actually published deserves a separate line. With the output taken from component
`e0` and the undeformed metric that the released code implies — no update rule for the deformation is
given anywhere — it reaches 20.3 % tracking efficiency under the conditions of the search, ranking
below 16 470 of the 20 000 randomly assembled configurations (11.6th percentile). Assigning the eight
physical quantities to components by throwing dice outperforms it four times out of five. Supplying an
arbitrary thermal deformation lifts the same component to 97.0 % when the controller starts near the
MPP, but the same controller reaches only 68.1 % from a 12 V start on the same profiles, and 48.8 %
under constant irradiance, where it never leaves the neighbourhood of its starting point. The
improvement is therefore a property of the starting point, not of the method.

### 3.2 Where the working configurations get their information

For the faithful reconstruction, symbolic expansion of component `e5` yields eight bilinear terms.
Used alone as the controller, the term `Vbus·(ΔG × Δ(ΔP))` (G = I/V) reaches 97.1 %, i.e. slightly
more than the complete octonion expression; the remaining terms only add noise. Replacing ΔP by a
constant drops the controller to 33.9 %, replacing I/V drops it to 25.0 %. For the best searched
configuration, the expansion is `2·(Vbus·(ΔI × ΔP) − (ΔP × Δ(ΔD)))`, whose second term is the P&O
rule: power change multiplied by the change of travel direction.

### 3.3 Why it converges more slowly

Because the associator is trilinear and alternating, three consecutive states give
`[X, Y, Z] = [X, a, b]` with `a = Y − X`, `b = Z − Y` (verified symbolically), i.e. a function of
the state and the last two increments. The decision therefore rests on second differences, whose
noise is 74 % larger than that of the first difference (1.66 W against 0.95 W in our setup), while
the informative part of the determinant largely cancels when the controller moves at a constant
step. Measured share of steps taken in the correct direction while more than 2 V away from the MPP:

| Controller | correct direction | net progress per cycle |
|---|---|---|
| P&O | **97.1 %** | 0.94 step |
| faithful reconstruction (e5) | 64.6 % | 0.29 step |
| best searched configuration | 68.7 % | 0.37 step |

Alternativity adds a structural drawback: `[X, a, a] = 0`, so the output vanishes exactly while the
controller is climbing steadily in one direction, and only reappears when the trajectory bends.

### 3.4 Step size from "memory", and unrestricted steps

Deriving the step size from `|A_k|` (gain tuned over ten orders of magnitude) gives a rise time from
12 V of 0.43 s for the best configuration, against 0.38 s for the classical variable-step P&O; the
remaining configurations need 3.8 to 4.7 s. Removing the step limit altogether degrades the octonion
controllers from 98–99 % to 55–63 %, because the magnitude of the associator has no units and no
calibrated relation to the distance from the MPP. Tuning gave the unrestricted P&O the smallest
possible step, i.e. the tuner itself rejected large jumps.

### 3.5 Fast transients and partial shading

On 100 ms cloud edges, most controllers need no recovery at all: the MPP voltage moves by about 1 V
over the whole irradiance range, while 30 K of temperature move it by 4 V. The faithful
reconstruction is the exception, dropping to 72.6 % on a ±1000 W/m²-in-2 s ramp with a mean voltage
error of 8.55 V, because the irradiance change injects a large spurious signal into the second
difference.

Under partial shading (three substrings, bypass diodes; e.g. 146 W at 17.6 V versus 80 W at 30.5 V
for one shaded substring) no configuration located the global peak, which is expected: a function of
the last two increments carries no information about a peak 10 V away. For reference, P&O with a
periodic full-curve scan found the global peak in 100 % of segments and reached 95.7 % under shading
against 87 % for plain P&O.

### 3.6 Cost

One associator requires four octonion products, i.e. 256 multiply-accumulate operations plus eight
subtractions per control step, against one comparison for P&O.

## 4. Limitations

We evaluate our reconstruction, not the author's method, because the latter is not defined; a
different assignment of quantities to components may behave differently, which is precisely the
problem — the spread across assignments is 20 % to 99 %. The plant model is quasi-static and covers
a single module without converter dynamics; multi-phase current sharing, the swarm synchronisation
claim and the hardware-protection claim are not tested, as no testable description of them exists.
Our results say nothing about whether an octonion layer could be useful in some other formulation.

## 5. Conclusion

For the MPP search itself, the octonion associator adds nothing. Where a configuration works, it
works because the expansion happens to contain the change of power, which is the quantity P&O uses
directly, one difference earlier and with half the noise. The component the published code computes
is identically zero. No configuration among 20 000 outperformed a method whose principle dates from the late 1960s.

For the proposal to become testable, three things are needed: (i) the multiplication table of the
claimed algebra and the update rule of its deformation, with parameters that can be identified from
measurements; (ii) the physical quantity, unit and scaling of each component, the definition of
X, Y, Z, and the mapping to the duty cycle; (iii) one measurable claim compared against a standard
method. A useful self-test for any proposed formalism is whether it reproduces the ordinary case:
write the plain converter in it and show that the standard equations follow.

A closing note of method. Statements produced with the help of a large language model read as
confident whether or not they are true; the model will elaborate an unfounded premise as fluently as
a sound one. That makes independent verification, of the kind attempted here, part of the work
rather than an optional afterthought — and it is cheap: the whole evaluation is a few hundred lines
of Python. A formalism drafted with such a model can be eight-dimensional and still dimensionally
inconsistent.

## Data and code availability

All code, models and experiment scripts are available at https://github.com/karagos01/octonion-mppt-eval and reproduce every
number in this comment (`run_all.sh`; the 20 000-configuration search takes about 14 minutes on
16 cores).

## References

1. M. Mazgal, *Intensional Maxwellian Formalism in Power Electronics: Causal Octonion Control of
   Distributed Multi-Phase Converters*, Zenodo, 2026. DOI 10.5281/zenodo.22877912 (v1.0),
   DOI 10.5281/zenodo.22914238 (v1.1).
2. J. C. Baez, *The Octonions*, Bulletin of the American Mathematical Society 39(2), 2002, 145–205.
3. T. Esram, P. L. Chapman, *Comparison of Photovoltaic Array Maximum Power Point Tracking
   Techniques*, IEEE Transactions on Energy Conversion 22(2), 2007, 439–449.
4. N. Femia, G. Petrone, G. Spagnuolo, M. Vitelli, *Optimization of Perturb and Observe Maximum
   Power Point Tracking Method*, IEEE Transactions on Power Electronics 20(4), 2005, 963–973.
5. H. Patel, V. Agarwal, *Maximum Power Point Tracking Scheme for PV Systems Operating Under
   Partially Shaded Conditions*, IEEE Transactions on Industrial Electronics 55(4), 2008, 1689–1698.
6. D. P. Hohm, M. E. Ropp, *Comparative Study of Maximum Power Point Tracking Algorithms*,
   Progress in Photovoltaics: Research and Applications, 2003.
7. EN 50530, *Overall efficiency of grid connected photovoltaic inverters* (dynamic MPPT efficiency
   test profiles).
