# Does the Octonion Associator Help MPPT? An Independent Numerical Evaluation

**Author:** karagos01 · **Status:** Working paper / technical comment · **Date:** 2026-09-25 ·
**Czech version:** [PAPER.cs.md](PAPER.cs.md)
**Licence:** this text CC BY 4.0 · the accompanying code MIT

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
obtained one difference later, with 74 % more measurement noise.

Two results do not depend on that reconstruction at all. First, the component labels carry no
information: of the 5040 ways to assign seven physical quantities to the imaginary units, 168 leave
the multiplication table bitwise identical, the basis automorphism group has order 1344, and the 5040
assignments collapse to 30 distinct algebras — the same four counts on two different published
orientations. Taking the stated units literally, 42 of 56 products equate incompatible dimensions,
and the constraint system over the whole table has rank 8 of 8, so the only dimensionally consistent
octonion over ℝ is the one whose components are all dimensionless; `e0` is in any case the
multiplicative identity rather than a slot. Second, the *Maxwellian* of the title is uncited and the
attribution tracks the algebra rather than the history: across all ten of the author's Zenodo
deposits Maxwell is named only in the four that pair him with octonions, zero times in the three
quaternionic ones, and "Treatise", "1873", "1865" and "Heaviside" occur zero times anywhere. All code
is provided.

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

## 4. The component labels carry no information

Section 2.1 of v1.1 is the only place where physics enters the octonion layer, so everything in the
proposal rests on one sentence:

> "Let O be a state octonion mapping dimensions e_0 … e_7 to physical potentials (e.g., scalar
> voltage, vector current, thermodynamic entropy, and electromotive gradients)."

`labels.py` measures three things about that sentence. None of them depends on our reconstruction of
the controller, so this section stands independently of Section 3.

**The seven imaginary units are interchangeable.** A physical assignment is a bijection from seven
quantities onto e_1 … e_7, so there are 7! = 5040 of them. Two assignments describe the same algebra
whenever the relabelling between them is an automorphism of the multiplication table. Counting those
directly, on the orientation used throughout this repository and on the one the same author later
published as SOTP Eq. (3):

| | measured |
|---|---|
| relabellings leaving the table bitwise identical | **21 of 5040** |
| …allowing a sign flip on each unit | **168 of 5040** |
| order of the full basis automorphism group | **1344** = 8 × 168 |
| distinct algebras reachable by relabelling | **30** |

Both tables give the same four numbers, so the result is a property of the octonions and not of a
particular orientation. Each of the 5040 assignments is algebraically identical to 167 of the others.
The deeper reason is that Aut(𝕆) over the reals is the 14-dimensional exceptional group G₂, which acts
transitively on basic triples: no invariant of the algebra distinguishes e_1 from e_4. Calling one
component "thermodynamic entropy" and another "electromotive gradient" therefore adds no constraint
that any computation can see. It is consistent with our search result in Section 3.2 — where a
configuration worked, it worked because of *which signal* was fed in, never because of *which slot*
it was fed into.

**No assignment of units is possible at all.** Taking the sentence literally, with e_0 in volts,
e_1 … e_3 in amperes, e_4 in J/K and e_5 … e_7 in V/m, **42 of the 56** products of two distinct
imaginary units equate incompatible dimensions: e_1 e_2 = e_4 requires [A][A] = [J/K], e_1 e_3 = e_7
requires [A][A] = [V/m], and so on. That is every one of them.

The stronger statement is that no relabelling or rescaling rescues it. Write one unknown dimension
exponent d_i per component and collect the constraint d_i + d_j = d_k from all 64 entries of the
table. The resulting system has **rank 8 in 8 unknowns**, so its solution space has dimension **0**:
the only dimensionally consistent octonion over ℝ is the one in which all eight components are
dimensionless. Two entries suffice to see why — e_0 e_i = e_i forces d_0 = 0, and e_i e_i = −e_0 then
forces d_i = 0. An octonion is an algebra over ℝ; its eight components are added to one another, so
they must share one dimension. Section 2.1 does not describe a richer object than a vector of eight
real numbers. It describes an expression that cannot be evaluated.

**e_0 is not a slot.** e_0 e_0 = e_0: it is the multiplicative identity. Assigning it "scalar
voltage" asserts V² = V. The mapping offers eight names for seven interchangeable places and one
place that is not a place.

This is the dimensional inconsistency referred to in Section 1 as "one dimensionally inconsistent
equation" in v1.0, located at its source rather than in a downstream formula.

## 5. Which algebra is attributed to Maxwell

The title of both versions is *Intensional Maxwellian Formalism*, and the abstract of v1.0 describes
the work as "restoring James Clerk Maxwell's original hypercomplex architecture". The claim is
load-bearing: it supplies the historical warrant for choosing a non-associative algebra in the first
place. It is also uncited — the preprint has no reference to any work of Maxwell's, and no
bibliography at all.

`census.sh` downloads all ten of the author's Zenodo deposits, converts each to text and counts the
keywords. It is the only script here that needs the network; the counts below are from a run on
3 October 2026.

| deposit | DOI suffix | Maxwell | octonion | quaternion | Treatise | 1873 | 1865 | Heaviside |
|---|---|---|---|---|---|---|---|---|
| acoustic projector | 22876757 | 2 | 19 | 0 | 0 | 0 | 0 | 0 |
| X-Ternary | 22877345 | 1 | 2 | 0 | 0 | 0 | 0 | 0 |
| MPPT v1.0 | 22877912 | 4 | 13 | 0 | 0 | 0 | 0 | 0 |
| MPPT v1.1 | 22914238 | 5 | 12 | 0 | 0 | 0 | 0 | 0 |
| COBAR | 22923154 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| momentum drive | 22934058 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| CQFT | 22959886 | **0** | 0 | 10 | 0 | 0 | 0 | 0 |
| PCTP | 22962188 | **0** | 0 | 15 | 0 | 0 | 0 | 0 |
| SOTP | 22962557 | **0** | 20 | 8 | 0 | 0 | 0 | 0 |
| openQL / openOL | 23112560 | 2 | 3 | 7 | 0 | 0 | 0 | 0 |

Three patterns fall out of the table.

**Maxwell is named only where octonions are.** The four deposits that pair Maxwell with an
engineering application all attribute to him a two-layer architecture whose first layer is octonionic:
"Following Maxwell's formalism, the system's underlying physical models are structured into two
distinct layers. The first layer is the Hidden Causal Memory, utilizing Octonion Algebra" (acoustic
projector); "a two-layer formalism inspired by Maxwell's original equations … Layer 1 (Hidden Causal
Memory): The system's causal history and interaction sequence are stored using non-associative
octonion algebra" (X-Ternary); "restoring James Clerk Maxwell's original hypercomplex architecture"
(MPPT v1.0).

**Where the algebra is quaternionic, Maxwell disappears.** CQFT and PCTP are built on quaternions and
name him zero times. SOTP, the octonionic member of the same trilogy, also names him zero times. In
place of a historical warrant, CQFT substitutes a different one — S. Adler's quaternionic quantum
mechanics — and argues that "this historical failure was not due to the nature of quaternions, but
rather the flawed axiom of scalar fungibility". An authority that attaches exactly where octonions
appear, and falls away when the algebra changes, is being selected to fit the conclusion rather than
consulted.

**Nothing is cited.** "Treatise", "1873", "1865" and "Heaviside" occur zero times in all ten
deposits. Whatever Maxwell is said to have originally written, no deposit says where.

What the historical record contains is narrower than the claim and points at the other algebra. The
1865 *Dynamical Theory of the Electromagnetic Field* is twenty scalar component equations in twenty
variables, with no hypercomplex algebra of any kind. Quaternions appear in the 1873 *Treatise* in
§§618–619, the last two articles of Chapter IX, under the head "Quaternion Expressions for the
Electromagnetic Equations" — scarcely two pages in a two-volume work of roughly a thousand. §618 is a
named list of the vectors and scalars already in implicit use; §619 states that *if* "vector" and
"scalar" are read as the vector and scalar parts of quaternions and *if* ∇ is taken as
quaternion-valued, then the equations already derived can be notated as (A)–(L). Only the operators
`S.` and `V.` occur there. Maxwell never multiplies two quaternions anywhere in the book: the full
product PQ does not appear. That is Hamilton's operator notation, not quaternion algebra, and the
real part of the quaternion carries no physics.

Maxwell's own position is on the record. He "endeavoured to avoid any process demanding from the
reader a knowledge of the Calculus of Quaternions", recommended "the introduction of the ideas, as
distinguished from the operations and methods of Quaternions", and wrote to Tait that he wanted to
"leaven my book with Hamiltonian ideas without casting the operations into a Hamiltonian form". He
also raised a physical objection: a vector quaternion squares to minus the square of its length, so a
kinetic energy written quaternionically comes out negative. In Chapter X the Gothic symbols appear
stripped of their quaternionic reading, and after that they are gone.

So the strong form of the claim is unavailable in either algebra. Octonions — Graves 1843, Cayley
1845 — appear nowhere in Maxwell, not even as a notational remark. Nor was anything suppressed:
those two pages are arguably the most consequential in the book, because Heaviside and Gibbs built
modern vector analysis out of them, and Maxwell himself coined "gradient", "convergence" and "curl"
in his 1871 essay on the classification of physical quantities. What was dropped from the quaternion
form was the real part, which Maxwell never filled with physics.

For completeness, openQL / openOL (3 October 2026) is the first deposit to attach Maxwell to
quaternions rather than octonions — "4D quaternion mechanics (for computing phase shifts, rotations,
and Maxwell's equations)". That is the algebra for which the historical claim has at least a weak
basis. It remains uncited and carries no mechanism, but the attribution now points at the right
algebra. That deposit is evaluated separately in `openql-eval`.

## 6. Limitations

We evaluate our reconstruction, not the author's method, because the latter is not defined; a
different assignment of quantities to components may behave differently, which is precisely the
problem — the spread across assignments is 20 % to 99 %. The plant model is quasi-static and covers
a single module without converter dynamics; multi-phase current sharing, the swarm synchronisation
claim and the hardware-protection claim are not tested, as no testable description of them exists.
Our results say nothing about whether an octonion layer could be useful in some other formulation.

## 7. Conclusion

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
16 cores).  Sections 4 and 5 are the cheap ones: `labels.py` reproduces every count in Section 4 in
about a second and needs nothing but numpy, and `census.sh` reproduces the table in Section 5 by
downloading the ten deposits from Zenodo — the only script here that uses the network.

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
8. J. C. Maxwell, *A Treatise on Electricity and Magnetism*, Clarendon Press, 1873, §§618–619
   ("Quaternion Expressions for the Electromagnetic Equations").
9. J. C. Maxwell, *A Dynamical Theory of the Electromagnetic Field*, Philosophical Transactions of
   the Royal Society 155, 1865, 459–512.
10. J. C. Maxwell, *On the Mathematical Classification of Physical Quantities*, Proceedings of the
   London Mathematical Society 3, 1871, 224–233 (where "gradient", "convergence" and "curl" are
   introduced).
11. N. Wheeler, *Theories of Maxwellian Design*, Reed College, 1998 (on §§618–619, the
   Maxwell–Tait correspondence, and the negative-kinetic-energy objection).
12. J. M. Chappell, A. Iqbal, J. G. Hartnett, D. Abbott, *The Vector Algebra War: A Historical
   Perspective*, IEEE Access 4, 2016, 1997–2004. arXiv:1509.00501.
13. A. Hurwitz, *Über die Composition der quadratischen Formen von beliebig vielen Variabeln*,
   Nachrichten der Gesellschaft der Wissenschaften zu Göttingen, 1898, 309–316.
