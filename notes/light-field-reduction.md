# The light field's cost, and why the knot count stopped being the question

The knot placement is settled: knots at fixed heights, positions constants of
the run, dropped position channel exactly zero. What was not settled is what
that costs, and the answer given so far — 1.52 times develop's run time — was
treated as the price of correctness.

It is not. It is the price of a reduction that costs knots times cohorts, and
that reduction does not have to cost that.

---

## Triage: 3

Expensive to reverse. It changes every number the model produces, so it lands
with a `scientific_version` bump and a re-blessing, and the requirement arrived
as a solution-verb — "resolve the dyadic lattice's discontinuities" — which is a
mechanism, not an outcome.

---

## Requirements ledger

Every quantity below is measured on this branch or on `develop`, on real FF16
stands taken from runs rather than constructed. The scripts are named against
each line.

**R1 — the adjoint differentiates the function the forward model computes.**
The channel a passive knot position drops, as a fraction of the field's range,
against the field's own error (`d8-anatomy.R`):

| knots | dropped channel | field error | ratio |
|---|---|---|---|
| 33 | 2.489e-01 | 2.198e-03 | 113x |
| 129 | 5.811e-02 | 1.331e-04 | 437x |
| 1025 | 7.626e-03 | 1.623e-06 | 4698x |

The channel falls as 1/K and the error as K^-2.2, so **refinement makes the
channel relatively worse** and no knot count rescues a canopy-tied grid. On a
lattice of constants the same measurement is `0.000e+00` at every spacing.

**R2 — crown-mean light, per cohort, at every stage: 1e-6 relative.**
develop's adaptive placement delivers 2.25e-06 to 8.66e-05 depending on the
stand, and at one stand leaves a cohort with no knot at all inside its own crown
(`d9-pareto.R`). The landed lattice at 0.05 delivers 1.01e-05 to 1.58e-05, and
its worst cohort is always the shortest.

**R3 — run time no worse than develop's.** The landed lattice is 1.52x.

**R4 — the field is continuous in the state, and the recorded computation's
shape does not depend on an active value.**

**R5 — the recorded block's input width is 2K + layers, at 480 slots per knot.**
Linear in the knot count, so knots are not free to the reverse pass even when
they are free to the forward one.

**Scarce resource.** Kernel evaluations in the field reduction: the build
evaluates every cohort's crown at every knot, so it costs `N x K` per stage, and
that is **62% to 88% of an FF16 right-hand side, rising with stand size**
(`d15-buildshare.R`). Everything else in this ledger is downstream of whether
that product stays a product.

---

## The floor

Keep the landed uniform lattice and accept the cost. It meets R1 exactly and R4
by construction.

It fails **R2** at the shortest cohort — 1.55e-05 against 1e-6 at lifetime 10,
because a 0.05 lattice puts three knots inside a recruit's crown where it puts
fifty inside the canopy's — and it fails **R3** at 1.52x against 1.0.

Those two failures are the same failure. The lattice cannot be refined to fix
the recruit because refinement is what costs 1.52x.

---

## Candidates

**A [first thought] — grade the lattice, spacing proportional to height.**
Move: re-shape the scarce resource. Commitment: constant *relative* resolution,
positions still constants. Pays for R2 at the recruit by moving knots down from
the canopy, where they were over-spent.

*Refuted by measurement.* At matched count it fails at the other end: at
lifetime 10, 176 knots graded gives the tallest cohort 4.16e-05 where 179 knots
uniform gives it 1.26e-07. The canopy top needs fine *absolute* resolution
because half the cohort tops are piled within a metre of it, so the field's
structure there is set by the spacing between cohort tops and not by any one
cohort's height. Eliminated on R2.

**B — dyadic lattice, level change smoothed by a blend.**
Move: make a discrete change continuous. Commitment: positions drawn from one
fixed lattice, resolution following the canopy by level.

*Refuted before the blend was designed.* A dyadic lattice **is** a uniform
lattice at whatever spacing its level picked, and its accuracy is uniform's to
three digits (lifetime 40: dyadic 175 knots 2.24e-04, uniform 0.100 176 knots
2.24e-04). It buys a bounded count and nothing else, and pays four field
discontinuities for it. There is nothing to blend *for*. Eliminated on R2.

**C — leave the placement alone and make a knot cheap.**
Move: attack the unit cost rather than the count. Commitment: the field's build
does not cost knots times cohorts.

The crown kernel expands:

```
(1 - (z/h)^eta)^2  =  1 - 2 (z/h)^eta + (z/h)^(2 eta)
```

so every z-dependence is `z^eta` or `z^(2 eta)` and every cohort-dependence is
`c_j`, `c_j h_j^-eta`, `c_j h_j^-2eta`. Three running sums over cohorts ordered
by height give every knot at once. Pays for R3 directly and for R2 by making
refinement affordable.

**Winner: C.** A and B are eliminated on R2 by measurement, not by taste. C is
the only candidate that touches the ledger line the other two are downstream of.

---

## The commitment

**The field's build cost is the number of knots plus the number of cohorts, not
their product.**

*Kept true by:* the build takes the whole knot vector and the whole cohort list
and returns the whole field, in one pass. There is no per-knot entry point on
the build path for a caller to put in a loop, so the quadratic form is not
expressible without reintroducing a function that no longer exists.
`compute_competition_and_slope(z)` stays for consumers that genuinely want one
height — the tests referee against it — but the environment no longer calls it.

---

## Kill question

*The assumption whose falsity makes this unnecessary:* that the reduction is a
large share of the work. If the leaf solve dominated, a faster field build would
buy a few percent and the 1.52x would have to come out of the knot count
instead.

*Verdict: survives, and only one of the two models needs it to.*

| | build, share of one right-hand side | swept, right-hand side becomes |
|---|---|---|
| FF16, lifetime 4 / 10 / 40 | 62% / 77% / 88% | 0.39 / 0.24 / 0.13 |
| TF24, lifetime 4 / 10 | 13% / 16% | 0.88 / 0.84 |

On FF16 the share **rises** with stand size, because the leaf work is linear in
cohorts and the build is quadratic; swept, the right-hand side is 2.6x to 7.7x
faster. On TF24 the leaf solve is 1.8 ms against the build's 0.26 ms, so the
same change buys 1.19x.

**That asymmetry is the answer to R3 rather than a weakness in it.** The 1.52x
this design exists to remove was measured on FF16, which is exactly where the
build dominates. On TF24 the knot count was never the cost: going from 166 to
315 knots moves the right-hand side by about 5%, so the lattice was close to
free there already. The two models fail R3 for different reasons and the same
change answers both.

---

## What survives deletion

- **The three running sums** — R3. Remove any one and the kernel no longer
  expands.
- **The downward sweep, carrying sums already scaled by the current height** —
  R3 and the arithmetic. Written with unscaled suffix sums the design is
  numerically dead: `h^-eta` at eta = 12 spans 1e19 over a real stand's heights,
  so the sums annihilate their own small terms. Scaled, every term is at most
  `c_j` because the sum runs only over cohorts with `h_j >= z`.
- **The general per-height reduction** — it is the referee, and the fallback for
  the shading models that do not expand.

Nothing else is added. No new type, no new control, no new configuration.

---

## What this settles

- Knot count stops being a design variable in the forward model. The build for a
  108-cohort, 348-knot stand goes from 37,584 kernel evaluations to 456.
- The accuracy requirement can be met by refining, because refining is nearly
  free. `d11-required.R` gives the requirement directly: about 0.022 near the
  seed height, which is fixed, and 0.006 to 0.015 of the canopy at the top.
- The dyadic lattice is not built, and the blend it needed is not designed.
- The knots-on-breaks placement is not built: at matched count it is *worse*
  than uniform for the consumer metric (0.4x to 0.9x), because the crown
  integral averages the localised break error away.
- Holding optical depth instead of light is not built: measured at 1.00x on
  every stand, because the stand's optical depth is only about 1.5 and Beer's
  law is not a strong enough nonlinearity to matter at that depth.

---

## What this makes hard

**The fast path is kernel-specific.** It expands the smooth Yokozawa profile and
nothing else. `FlatTopSoftBox` carries its own smoothstep and would need either
its own expansion or the general reduction; `FlatTopBox` already refuses to
build a light field at all. If a fourth profile arrives, it arrives with a
decision to make rather than a free ride.

**One digit of arithmetic.** The expansion computes `1 - 2u + u^2` where the
direct form computes `(1 - u)^2`, so it cancels where `u` approaches 1. Measured
against the direct sum's own rounding floor, the absolute error in A is 3.6e-16
to 7.5e-14 where `eps * A(0)` is 1.4e-16 to 6.3e-15 — within about one digit,
across stands from 8 to 300 cohorts. Absolute error in A is what matters,
because it is relative error in `exp(-A)`; the relative error in A itself is
large only at the canopy top, where A goes to zero and no consumer reads it.

**R5 is untouched.** The recorded block still declares `2K + layers` inputs at
480 slots per knot, so knots stay linear in the reverse pass even once they are
nearly free in the forward one. A crown integral under a fixed rule touches at
most `n+1` values and `n+1` slopes whatever the canopy does, so the width is
available; taking it is separate work and this design does not do it.

---

## Kill condition

**A production model whose competition profile does not expand**, or a
configuration where the leaf solve is so dominant that the build is a few
percent. The first hands off to the general reduction, which is kept for exactly
that reason. The second hands off to candidate A — if a knot must be paid for
after all, then where it goes matters again, and `d11-required.R` is the curve
to place against.

---

## The design

One reduction, evaluated as a sweep from the canopy down.

Let `T1(z) = sum_j c_j (z/h_j)^eta` and `T2(z)` the same at `2 eta`, both over
cohorts with `h_j >= z`, and `S0(z)` the plain sum of `c_j` over the same set.
Then

```
A(z)  = S0 - 2 T1 + T2
A'(z) = 2 eta (T2 - T1) / z
```

and moving from a knot to the next one below, at `z' < z`:

```
T1 <- T1 (z'/z)^eta                      rescale what is already in
T2 <- T2 (z'/z)^(2 eta)
for each cohort with z' <= h_j < z:      admit the ones now in range
    u  <- (z'/h_j)^eta
    S0 <- S0 + c_j ;  T1 <- T1 + c_j u ;  T2 <- T2 + c_j u^2
```

Every term is at most `c_j`, because the sums run over `h_j >= z`. Each cohort
is admitted once and each knot is visited once, so the pass is `O(N + K)` with
two powers per knot and one per cohort.

Three things the implementation has to hold, each of which is a way to lose it:

**The cohort order is a sorted view, not the storage order.** Heights invert on
this coordinate — a younger cohort can overtake an older one — so a sweep that
walks cohorts in storage order admits them at the wrong knots. The census
already needs a sorted view for the same reason.

**One set of sums per species.** `eta` is a strategy parameter, so cohorts of one
species share it and cohorts of different species do not. The pass is
`O(N + K x species)`.

**The slope comes from the same sums as the value.** `A'` is formed from the
`T1` and `T2` that produced `A`, in the same pass, so value and slope cannot
come from two merges that differ in their last bits.

The transpose is the same shape. `T1` and `T2` are linear recurrences in the
knot index, so their adjoint is one reverse pass; the cohort rows factor as
`(z_i/h_j)^eta = z_i^eta h_j^-eta`, so each cohort's height adjoint reads a
prefix sum over knots that is formed once. `O(N + K)` both ways, with the same
rescaling for the same reason.

---

## What would falsify this

- **The build is not the share it is measured to be on a production
  configuration.** The number here is FF16 at three stand sizes; a configuration
  whose leaf work dominates moves the whole argument.
- **The sweep and the direct sum disagree by more than the direct sum's own
  rounding**, on a stand the model reaches. The check needs no reference: build
  the field both ways and compare, and the scale to beat is `eps * A(0)`.
- **The transpose's adjoint does not match a forward tangent of the sweep.**
  That is the same identity the reduction transposes already answer to, and it
  needs no gradient reference.
- **A stand exists where the requirement curve is not two-humped** — where the
  spacing a cohort needs is not set by the seed height at the bottom and the
  canopy at the top. Then the placement question reopens on a different shape.
