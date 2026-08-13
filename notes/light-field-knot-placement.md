# The light field's knot placement: a grid of constants

Two decisions arrived together and only one was argued for. The interpolant changed from a
value-only cubic with solved slopes to a Hermite carrying a value and a slope at every knot;
and the placement changed from adaptive refinement on value error to `u_k · h_max` with 65
uniform fractions. This states what the placement should be, and prices every alternative
against the three constraints that are actually in tension.

**The interpolant is settled and stays.** Report 03 §3.1 is the argument: a spline that solves
a tridiagonal system for its slopes has no local support, so every knot value influences every
query and its transpose is a banded solve of run-dependent width. Hermite gives four non-zeros
per query. Nothing below revisits that.

A **second** reason for it turned up here, and it is what makes the rest of this note possible:
**a Hermite interpolant's spans are local, so extending the grid upward changes no span a query
reads.** On a solved-slope spline an appended knot changes every slope and therefore every
query. So the interpolant change is not merely compatible with a fixed grid — it is the thing
that admits one. The two decisions are separable in one direction only.

Everything below is measured on the model. Where a number disagrees with a report, the
disagreement is the finding.

---

## 1. What the placement has to satisfy

Three constraints, and each candidate wins a different one.

**Gradient safety.** Knot positions are passive: report 03 §3.3 is right that dropping their
derivative is the correct treatment, because moving a knot changes the interpolant and not the
interpolated function. But a *dropped* channel is only harmless if it is *small*, and where the
positions are built from `h_max` it is not. The chain `tallest cohort's height → every knot
position → every crown integral` is re-formed and dropped at every stage.

**Accuracy.** What consumes the field is a crown integral: cohort `k` reads
`M(h_k) = ∫₀^{h_k} L(z) q(z,h_k) dz`. Substituting `t = u^η` gives `M = ∫₀¹ L(h t^{1/η}) 2(1−t) dt`,
so the weight is trivial and the sampling is concentrated at the crown top. **That, and not a
uniform norm over the profile, is the metric a placement should be judged on** — it is also the
metric `resource_spline.h`'s own comment already used.

**Forward cost.** An equally spaced grid indexes by arithmetic instead of a search, which is
what the current placement was chosen for.

The field's structure follows from one fact. With `η = 12`, `Q̃(ν) = (1−ν^η)²` is 0.9995 at
`ν = 0.5` and 0.867 at `ν = 0.8`, so **each cohort's shading transition occupies the top ~40% of
its own height**, and `A(z)` has a curvature break at every distinct cohort top. The field is
flat well below the shortest cohort and structured throughout the range where cohorts live.

---

## 2. The dropped position channel, measured

Report 03 §3.3's falsifier, run per cohort: hold the profile fixed, move one cohort's height,
and read what the interpolant does. The difference is the term the passive positions drop.
Normalised as the response to a *relative* change in that height, as a fraction of the field's
range, so it is dimensionless and step-independent. Shaded stand, eight cohorts, `h_max = 18`,
`L(0) = 0.055`.

| placement | knots | 18.0 | 15.5 | 13.0 | 10.5 | 8.1 | 5.6 | 3.1 | 0.6 | quad-sum |
|---|---|---|---|---|---|---|---|---|---|---|
| uniform × `h_max` | 65 | 2.52e-01 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 2.52e-01 |
| breaks + fill | 65 | 1.75e-01 | 3.09e-01 | 2.77e-01 | 4.09e-02 | 1.51e-02 | 6.85e-03 | 5.31e-02 | 7.75e-02 | **4.62e-01** |
| fixed absolute, Δ=0.25 | 75 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | **0** |
| fixed absolute, Δ=0.10 | 183 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | **0** |

Three readings.

**A grid of constants drops nothing, exactly.** Not small — `0.000e+00`, on every cohort.

**Tying knots to cohort tops is the worst of the three, and the mechanism is visible.** It
drops a term for *every* cohort rather than one, so its quadrature sum is 1.8× the
canopy-following grid's. This confirms the note that breaks+fill costs position sensitivity, and
locates the cost: it is eight channels where uniform × `h_max` has one.

**One row is understated and should not be read as adaptive's number.** The adaptive placement
was driven from an unperturbed patch, so only the domain's endpoint moved. Real adaptive
positions depend on the whole profile, so its true channel has a term per cohort like
breaks+fill.

### 2.1 End to end

Forward tangent of the right-hand side against a plain-double difference of the rates, which
re-runs the field build and so carries the chain the tangent drops. Shaded stand at the darkest
shade the tangent survives (min read light 0.203).

| placement | tallest height column | median other height column | ratio |
|---|---|---|---|
| uniform × `h_max`, 65 | **2.600e-03** | 1.531e-06 | 1698× |
| fixed absolute, Δ=0.10 | 3.032e-05 | 1.533e-06 | 19.8× |
| fixed absolute, Δ=0.25 | 8.022e-05 | 1.534e-06 | 52.3× |

Only the tallest cohort's column moves; every other height column sits at 1.53e-06 in all
three. One construct changed, one number moved.

**The fixed grid's remaining 3e-05 is a property of the fixture, not of the design.** It appears
only when `h_max` sits *exactly* on a knot. Off a knot:

| Δ | `h_max` | `h_max/Δ` | straddles a knot | tallest | median other |
|---|---|---|---|---|---|
| 0.10 | 18.000 | 180.00 | yes | 3.032e-05 | 1.533e-06 |
| 0.10 | 18.037 | 180.37 | no | **1.801e-06** | 1.507e-06 |
| 0.10 | 18.500 | 185.00 | yes | 3.155e-05 | 2.211e-06 |
| 0.10 | 18.617 | 186.17 | no | **1.714e-06** | 1.445e-06 |
| 0.25 | 18.000 | 72.00 | yes | 8.022e-05 | 1.534e-06 |
| 0.25 | 18.037 | 72.15 | no | **1.801e-06** | 8.032e-07 |

**With the canopy off a knot the tallest cohort's column is indistinguishable from every other
one.** The severed channel is not reduced, it is gone. Stated as the ratio the first table used,
confirmed from a clean process on a fixture whose `h_max` is not a multiple of Δ:

| placement | tallest | median other | ratio |
|---|---|---|---|
| uniform × `h_max`, 65 | 2.578e-03 | 1.212e-06 | **2127×** |
| fixed absolute, Δ=0.10 | 1.801e-06 | 1.507e-06 | **1.2×** |

### 2.2 What is left is a kink, and it is smaller than what it replaces

Two things coincide when `h_max` sits on a knot: the count changes, and the tallest cohort's
curvature break crosses a knot. Holding the count fixed separates them — and the ceil rule
(count 76→77) and a held count of 100 give **identical numbers to every digit**. So the count
contributes exactly nothing and the break crossing is the whole of it.

Walking the tallest cohort's height across a knot, Δ=0.25:

| placement | value continuity (2nd difference) | `d(error)/dh` range | jump |
|---|---|---|---|
| fixed absolute | 3.573e-07 | **[−9.6e-07, 1.413e-04]** | 1.411e-04 |
| uniform × `h_max` | 8.505e-08 | **[1.801e-04, 2.558e-04]** | 9.196e-06 |

The value stays continuous, so this is a kink and not a jump. And the comparison is the design
argument: **the fixed grid's error in `d(crown mean)/dh` is ≤ 9.6e-07 away from a crossing and
spikes to 1.4e-04 at one; the canopy-following grid's is 1.8e-04 to 2.6e-04 everywhere and
never approaches zero.** A grid of constants replaces a *persistent* first-order bias with an
*isolated* kink of comparable size — 200× more accurate almost everywhere, and never worse. A
crossing is measure-zero in state space, and a tangent at one returns the one-sided derivative
of the branch in force, which is what report 08 §4 says a transpose ought to be the transpose
of.

---

## 3. Accuracy, and the honest result

Consumer metric: max over cohorts of the relative error in crown-mean light. Floor cohort at
0.4 m throughout, so the smallest plant's requirement is in every row. Knot counts in brackets.

| `h_max` | canopy 65 | canopy 257 | fixed 0.25 | fixed 0.10 | fixed 0.05 |
|---|---|---|---|---|---|
| 1.5 | 7.17e-08 | 7.14e-10 | 5.06e-06 (9) | 2.97e-06 (18) | 2.55e-07 (33) |
| 6.0 | 1.09e-06 | 7.18e-08 | 2.90e-05 (27) | 2.93e-06 (63) | 2.58e-07 (123) |
| 18.0 | 3.97e-05 | 7.25e-07 | 2.38e-05 (75) | 1.42e-06 (183) | 1.18e-07 (363) |
| 35.0 | **6.25e-04** | 4.88e-06 | 5.69e-05 (143) | **1.09e-06** (353) | 1.18e-07 (703) |

**Neither family dominates per knot, and the fixed grid is the worse of the two on a young
stand.** At 1.5 m the incumbent is 7.17e-08 against Δ=0.10's 2.97e-06 — **40× better**, and it
stays ahead to about 12 m. The canopy-following grid holds constant *relative* spacing, which
suits a stand whose cohorts are all of similar height; the fixed grid holds constant *absolute*
spacing, which is what a 0.4 m seedling under a 35 m canopy needs. The crossover is where the
stand's height range opens up.

**What does separate them is the shape of the error over a run.** The incumbent's grows four
orders as the stand matures — 7.17e-08 to 6.25e-04 — and it degrades exactly as the stand
closes into the shade where the field carries the competition. Fixed Δ=0.10 is **flat at about
1e-06 across a 23× range of canopy height.** That is what report 03 §2 means by an error budget:
a quantity whose accuracy varies by four orders over the run does not have one.

**Where the canopy grid is genuinely more efficient, stated carefully because it is easy to
overclaim in either direction.** On FF16's offspring production — a whole-run functional, and
the sharpest one available — both families converge to the same limit (the two finest agree to
7.57e-06), and the canopy grid reaches it on far fewer knots:

| policy | knots | vs limit | | policy | knots | vs limit |
|---|---|---|---|---|---|---|
| canopy 65 | 65 | −2.178e-03 | | fixed 0.25 | 54 | −1.663e-02 |
| canopy 129 | 129 | +1.098e-04 | | fixed 0.10 | 131 | −1.072e-03 |
| canopy 257 | 257 | −1.885e-06 | | fixed 0.05 | 261 | −1.682e-04 |
| canopy 513 | 513 | +8.316e-06 | | fixed 0.025 | 519 | −5.239e-05 |

**So on FF16 the fixed grid needs roughly 3–4× the knots for equal accuracy.** That is the real
price, and §4 is why it is worth paying.

**But per knot is not the comparison that decides anything, and read against the incumbent the
same table goes the other way.** Canopy 65 is **−2.18e-03** from the limit; fixed Δ=0.10 is
**−1.07e-03**. The fixed grid is twice as close on the metric that most favours the canopy
family. Collecting every comparison against the shipped placement rather than against a matched
count:

| metric | canopy 65 (shipped) | fixed Δ=0.10 | |
|---|---|---|---|
| crown mean, `h_max`=35 | 6.25e-04 | 1.09e-06 | 570× better |
| crown mean, `h_max`=18 | 3.97e-05 | 1.42e-06 | 28× better |
| crown mean, `h_max`=6 | 1.09e-06 | 2.93e-06 | 2.7× worse |
| crown mean, `h_max`=1.5 | 7.17e-08 | 2.97e-06 | 40× worse |
| FF16 offspring, vs limit | −2.18e-03 | −1.07e-03 | 2× better |
| dropped position channel | 2.52e-01 | **0** | exact |

**The only place the proposal loses is a young canopy, and there both errors are ≤ 3e-06.** The
shipped placement is also not converged on FF16, so this is a re-blessing rather than a
bit-identical swap either way.

---

## 4. Cost, and why the knots are affordable

Build and run cost, `h_max = 18`, eight cohorts:

| what | canopy 65 | fixed 0.10 | ceiling 25 m | ceiling 40 m | ceiling 60 m |
|---|---|---|---|---|---|
| knots | 65 | 183 | 251 | 401 | 601 |
| one field build | 2.75e-05 s | 1.78× | 1.85× | 2.04× | 2.24× |
| 20-year SCM run (106 / 401 knots) | 1.631 s | **0.98×** | — | 1.06× | — |
| block inputs, `2K + layers` | 185 | 421 | 557 | 857 | 1257 |

Sweep cost against input width, one cohort's block:

| | knots | inputs | reverse sweep | forward tangent |
|---|---|---|---|---|
| canopy 65 | 65 | 185 | 1.533e-02 s (1.00×) | 2.067e-01 s (1.00×) |
| canopy 129 | 129 | 313 | 1.492e-02 s (**0.97×**) | 3.537e-01 s (**1.71×**) |

Two measurements matter here.

**The knot count is nearly free forward.** A 20-year run at 106 fixed knots ran at 0.98× the
65-knot canopy grid, and going from 65 to 129 canopy knots cost 3%. The leaf solve dominates,
exactly as report 00 says, so the field build is not where a run's time goes.

**The reverse sweep is flat in K and the forward tangent is not.** At 313 inputs against 185,
one cohort's reverse sweep ran at 0.97× and its forward Jacobian at 1.71×. That is reverse
mode's own property — sweep cost scales with output rows, not inputs (report 01 §3) — and its
consequence here is worth stating plainly: **a larger K costs the reference instrument, not the
product.** A forward tangent is a test-time reference on small fixtures, so paying 3–7× there
to remove a 1698× defect in the thing being refereed is the right trade.

**One consequence is not free and belongs in the plan.** `n_cohort_reads() = 2K + layers`, so
the recorded block's input vector widens with K while the number a cohort actually *reads* stays
bounded by the quadrature rule at ~4 non-zeros per query (report 03 §3.1, report 07 §1). A
larger K widens the gap between what is materialised and what is read, which moves report 07
§7's item 2 — use the sparsity of the light coupling — from an economy to something the design
leans on.

---

## 5. The design

**A grid of constants, uniform in absolute height, spanning `[0, ceiling]`.**

```
x_k = k · Δ,   k = 0 … ⌈ceiling / Δ⌉
```

- **Positions are constants of the run.** `∂x_k/∂(state) ≡ 0` is a fact, not a treatment, so
  there is no term to drop and report 03 §3.3's channel does not exist.
- **The count is a constant of the run**, which is what report 05 §5 relies on when it says K
  does not vary within a run, and what `n_cohort_reads()` needs in order to be answerable before
  the first build (§7).
- **The spacing stays equal**, so `hermite_interpolator`'s arithmetic index survives intact —
  the property the current placement was chosen for.
- **Knots above the canopy are exactly inert.** `Q̃` and `Q̃'` both vanish at `ν = 1`, so above
  `h_max` the field is exactly `(1, 0)` and its derivative with respect to every cohort height
  is exactly zero. Measured: a 401-knot ceiling grid and a 183-knot append-only grid give
  **bit-identical rates** (`max |diff| = 0.000e+00`) and a leaf-area census agreeing to eight
  digits over a 20-year run. They can therefore be filled without calling the reduction at all,
  which makes the ceiling's extra knots free in the build as well as inert in the answer.

**Δ is a control parameter, defaulting to 0.10 m**, which buys ~1e-06 consumer error flat across
the run. Δ=0.05 buys ~1e-07 for twice the knots if a metric ever needs it.

**Δ must not be derived from the seed height.** That quantity now carries a derivative (report
04 §6.1), and a grid position built from it would drop a term in precisely the direction report
04 just recovered. The old `tol`/`nbase`/`max_depth` slots are the natural place for Δ, and
taking them for it retires three arguments that are currently accepted and ignored.

**The ceiling is preferred over append-only**, on §7's evidence rather than on taste.

---

## 6. What is refuted

Each of these was tried and measured, and each is worth recording so it is not re-derived.

**A geometric grid, spacing proportional to height.** The argument for it is that each cohort's
transition band is a fixed *fraction* of its own height, so the requirement is scale-invariant.
Measured, it is 50× worse than uniform-absolute at matched count — `r=1.10` at 76 knots gives a
consumer error of 5.94e-03 against 2.66e-05 for Δ=0.25 at 75. A fixed relative spacing spends
most of its knots in the flat sub-canopy and is left too coarse across the tallest cohort's
band, which is both the widest feature and the one the crown integrals weight most.

**Snapping one knot to `h_max`.** If the canopy grid's advantage were its knot on the tallest
break, a fixed grid with one snapped knot would recover it for one channel instead of 65.
Measured, it recovers **nothing** — consumer error 1.197e-06 either way — and reintroduces a
dropped channel of 3.16e-02. Strictly worse, and the hypothesis that motivated it is false.

**Interpolating the optical depth instead of the transmittance.** Beer's law is currently applied
before interpolation, so the knots hold `L = exp(−A)`. Holding `A` and exponentiating at the
query is identical at Δ=0.10 (1.197e-06 both) and better by at most 30%, only at `L(0) ≈ 5e-05` —
shade the gradient already refuses. It would move an `exp` from K knots per build to every
query and change what the reduction's transpose is a transpose of, for that.

**The interpolant's undershoot.** Report 03 §4 has the light floor and the interpolant's
monotonicity guard sitting on one lever, `k_I · LAI`. The floor half is right; the undershoot
half is not reachable. Measured at five shade levels down to `L(0) = 6e-29`, the minimum
interpolated value tracks the minimum exact value and is never negative. The reason is
structural: `L` is monotone increasing in `z` and the knots carry its exact slopes, so the
first span cannot dip below the ground knot. **The `std::max(S(0.0), …)` clamp in
`get_value_at_height` never fires**, and the two clamps are not in fact co-triggered.

**Reverting to adaptive placement.** Beyond report 01 §7 item 8's prohibition on a
state-dependent node count, it is also *worse on the consumer metric* than the incumbent on a
shaded stand — 1.749e-04 against 4.845e-05 — because it chases value error, which lives at the
canopy top, and under-resolves the small cohorts whose crown means are relative.

---

## 6.5 What the change costs the suite

Measured by forcing the placement at run time and running the files that read this field most
directly. Baseline under the current placement: 692 assertions, no failures — so the switch
itself is what moves them, not the scaffolding.

| file | pass | fail | what fails |
|---|---|---|---|
| `test-canopy-methods.R` | 92 | 94 | see below |
| `test-strategy-ff16.R` | 52 | 3 | `offspring arrival`, the bit-identity guard |
| `test-environment.R` | 33 | 0 | |
| `test-strategy-ff16-reference-comparison.R` | 17 | 0 | |
| `test-scm.R` | 140 | 0 | |
| `test-patch.R` | 167 | 0 | |
| `test-environment-TF24.R` | 94 | 0 | |

**The 94 is one test and about five numbers.** 92 of them are a single assertion —
`x == u * h$height_max` — evaluated once per recorded step in
`"the light interpolant's knot positions are run-constant"`. That test asserts the *incumbent's*
placement rule, so it is rewritten rather than re-blessed, and **its replacement is stronger**:
positions become `k · Δ`, identical at every step and independent of `height_max` altogether,
which is what the test's own name asks for. The rest is one deep-crown SCM baseline, one
Beer's-law knot assertion, and FF16's three offspring-arrival values.

So the migration is: one structural test rewritten, and a re-blessing of ~5 numbers at the
1e-06 (TF24 census) to 1e-03 (FF16 offspring) level. §3's convergence table is the argument
that the new numbers are the better ones — the incumbent sits 2.2e-03 from the limit on FF16.

---

## 7. Two things found by breaking them

**Append-only cannot answer `knot_count()` before the first build, and callers need it to.**
Making `knot_count()` report the laid-out count rather than the planned one produced heap
corruption across the suite: `n_cohort_reads() = 2K + layers` sizes buffers that are filled
later. Under a ceiling the count is known at construction. Under append-only it is not, and it
changes mid-run. **That is the argument for the ceiling** — not the accuracy, which is
bit-identical, but that a run-constant K is what the surrounding code is written against.

**A partial and a total, again.** The fixed grid's tallest column reads 3e-05 or 1.8e-06
depending only on whether the fixture's `h_max` happens to be a multiple of Δ. A rung that
measured one of those and not the other would draw opposite conclusions. Report 08 §6's
fixture assertions should include **`h_max` not a multiple of Δ**, for the same reason they
already include non-commensurate heights.

---

## 8. What this does not fix, and what it is blocked behind

**The shaded-stand non-finite tangent bounds what can be refereed.** Driven into shade the
forward tangent of the right-hand side returns NaN, and the threshold is sharp: clean at a
minimum read light of 0.203, non-finite at 0.143 and below. The structure is sharper than the
existing record has it — **38 of 73 rate rows are wholly NaN and the count is identical
(2774/5329) at three different shade levels**, so it is a switch rather than a degradation. The
affected rows are height, fecundity, storage, survival-weighted offspring and five soil rates.
This is not caused by the placement: it is present under both, unchanged. It means the
end-to-end validation above stops at `L(0) ≈ 0.2`, and a realistically closed canopy
(`L(0) ≈ 0.05`, which is report 06 §5's reference stand) cannot be refereed by a tangent at all.

**Two placement-independent column disagreements** showed up beside it and are almost certainly
report 08 §4.7's trap two rather than defects: a `storage` column at 2.9e-01 and two
`log_density` columns at ~2.5e-03, identical under all three placements. The difference's step
is `max(|state|, 1) · rel`, which is absolute for a small reserve, so the row it returns is not
a derivative of anything.

---

## 9. What would falsify this

- **The consumer metric is the wrong one.** Everything above ranks placements on crown-mean
  light because that is what a cohort reads. If a metric turns out to be sensitive to the
  field's *slope* independently — the reverse pass contracts `λ_Λ'` against `∂_z E^comp` — then
  slope error needs its own budget, and breaks+fill's 26× advantage there (5.87e-04 against
  1.53e-02 at 65 knots) would have to be weighed against its eight dropped channels.
- **Δ=0.10 is not enough for a metric someone wants.** The flat ~1e-06 is measured on crown-mean
  light and on a leaf-area census. FF16's offspring production sits at 1.07e-03 there, and
  whether that matters is a modelling judgement, not a numerical one.
- **The knot count is not free at production width.** Every cost figure here is at eight to
  forty cohorts. The build is `O(K × cohorts)`, so a stand at production width with a 60 m
  ceiling is the case to re-time before committing to Δ=0.05.
- **A cohort top crossing a knot is not measure-zero in practice.** The kink of §2.2 is harmless
  because the canopy crosses a knot at isolated times. If a stand can *sit* at a crossing — a
  stalled canopy, a fixed-point iteration — the kink is permanent and the argument changes.
- **`n_cohort_reads()` widening is not affordable.** If the block's input materialisation turns
  out to matter at production width, the sparsity work of report 07 §7 item 2 becomes a
  precondition of this change rather than a companion to it.
