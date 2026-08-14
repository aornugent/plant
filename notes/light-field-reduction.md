# Building the light field in one pass

The field is built by evaluating, at every knot, a sum over every cohort. That
product is 62 to 88 percent of an FF16 right-hand side and it does not have to
be a product: the crown profile is a polynomial in `u^eta`, so three running
sums over cohorts ordered by height give every knot at once.

This is the implementation. The alternatives it was chosen over, and the
measurements that eliminated them, are in the closing section.

---

## What is being claimed, and what it rests on

The reduction is a trapezium integral over the size distribution of one
quantity per cohort:

```
A(z) = sum_k w_k d_k k_I a_k Q(z / h_k)          Q(u) = (1 - u^eta)^2, u <= 1
```

`w_k` are the quadrature weights, `d_k` the node's density, `a_k` its
competition effect and `h_k` its height. **The only `z`-dependence is `Q`**, so
writing `c_k = w_k d_k k_I a_k` the field is `sum_k c_k Q(z/h_k)` with `c_k`
independent of the height it is read at. That is the whole basis of the sweep,
and it holds in both coordinates the forward model supports, because the
abscissa only ever enters `w_k`.

Expanding,

```
Q(u)     = 1 - 2 u^eta + u^(2 eta)
dQ/dz    = (-2 eta u^eta + 2 eta u^(2 eta)) / z
```

so with `S = sum c_k`, `T1 = sum c_k u_k^eta` and `T2 = sum c_k u_k^(2 eta)`,
each over the cohorts whose crowns reach `z`,

```
A(z)  = S - 2 T1 + T2
A'(z) = 2 eta (T2 - T1) / z
```

**The sums must be carried already scaled by the height they were last read
at.** Written as `z^eta sum c_k h_k^-eta`, the factor `h^-eta` at `eta = 12`
spans 1e19 over one stand's heights and the sum loses its own small terms.
Carried scaled, every term is at most `c_k`, because the sums run only over
cohorts with `h_k >= z`. Measured against the per-cohort sum's own rounding on
a real 93-cohort stand: value agrees to 4.4e-16 where `eps * A(0)` is 7.5e-16,
slope to 3.6e-15 on a range of 17.5.

Two facts about the model make the truncation free rather than a special case.
`Q` and `q` both vanish at `u = 1`, so a cohort admitted at the first knot below
its own crown top contributes nothing at the knot above it, and the join is C1 —
which is why the branch carries no derivative. And a cohort admitted once stays
admitted, because every knot below is also below its crown.

---

## What the model actually does today, which the sweep has to reproduce

Read from the forward model rather than assumed; two of these were wrong on a
first reading and are the reason the sweep is specified this precisely.

**The abscissa is not the birth date by default.** `abscissa_of` returns the
introduction time in the birth-date coordinate and **minus the height**
otherwise, and the shipped SCM default is the height coordinate. The reverse
sweep refuses any other, so both are live: the weights differ between them and
nothing else does.

**The trapezium is accumulated as a running sum**, `(x0 - x1)(f1 + f0)` over
consecutive nodes, halved once at the end, with the inflow boundary node
appended as a closing interval. Re-associating that as `sum_k w_k f_k` is the
same number in exact arithmetic and a different one in floating point, so the
sweep moves the model's numbers and needs a re-blessing.

**The early exit is not a truncation of the sum.** The loop stops at the first
node below the query height; every node past it contributes `f = 0` at both ends
of its interval, so including them adds exact zeros. The sweep may therefore
admit every cohort without changing the value.

**The field is stored in double and the reduction is not, and conflating those
deletes a channel.** The knots are a `ResourceSpline<double>` and the
environments are not templated, so what the field *stores* is plain double.
But `Individual` holds `Internals<value_type>`, and the strategy evaluates the
profile at those internals — so the reduction reads the cohort's **active**
height through the cached `height_inverse`, and a tangent seeded on a height
propagates into the field through it.

Taking the height passively in the sweep is therefore not a free simplification:
it severs every cohort's height channel out of the field. It is also invisible
in the value, which agrees to rounding either way, and invisible to every check
that does not differentiate — measured, the whole non-gradient suite passes with
the channel gone, and the ladder's tangent-against-difference margin moves from
1.0e-08 to 2.9e-04.

So the height is active where the profile consumes it, and passive only where
the model already reads it passively: the quadrature abscissa, the descent's
knot-over-knot ratio, and the comparison deciding which knot a cohort enters at.
On the height coordinate the abscissa and the height are the same quantity read
two ways, which looks wrong and is what the per-height reduction does.

**The amplitude is already on the node.** `k_I * aux(COMPETITION_EFFECT)` is the
scale and `aux(HEIGHT_INVERSE)` the reciprocal crown top, both cached, so the
sweep needs no allometry.

---

## The design

**Commitment: the field's build cost is knots plus cohorts, not their product.**

*Kept true by:* the build takes the whole knot vector and returns the whole
field. There is no per-knot entry point on the build path to put in a loop, so
the quadratic form is not expressible without reintroducing a function that the
environment no longer calls.

### Where each piece goes

**`CanopyShape` owns the expansion**, beside `Q_and_q`, because it owns `Q`. A
sum built from a second, separately written copy of the profile would agree with
the first nowhere except by accident, and the value would not show it. It gains
the three running sums and the two operations on them, and a predicate saying
whether its profile is a polynomial in `u^eta` — the box profiles are not, and
one of them cannot build a light field at all.

**`Species` owns the sweep**, because it owns the nodes, the quadrature weights
and the strategy. It gains one entry point that fills value and slope at every
knot, and keeps the per-height reduction for the transpose, for R, and as the
fallback when the profile does not expand.

**`Patch` hands the whole grid down** instead of a per-height callback, and
**`ResourceSpline` asks once** instead of per knot. The environments apply
Beer's law over the filled vectors rather than inside a per-height lambda.

### The sweep

Knots descend; cohorts are admitted in decreasing height order; the weights come
from the abscissa in whichever coordinate is in force.

```
for each knot z, from the top down:
    descend the sums to z            T1 *= r ; T2 *= r * r,  r = (z/z_prev)^eta
    admit every cohort with h_k >= z not yet admitted:
        u = (z / h_k)^eta
        S += c_k ; T1 += c_k u ; T2 += c_k u * u
    value = S - 2 T1 + T2
    slope = 2 eta (T2 - T1) / z
```

Each cohort is admitted once and each knot visited once: two powers per knot,
one per cohort.

### Four things that lose it

**The ground knot has no ratio.** `A'` divides by `z`, which is zero at the
lowest knot, and the limit there is not reachable from the sums at `eta = 1`.
That knot is evaluated by the per-height path — one knot out of hundreds — which
is also what the existing code does for its own limit.

**The cohort order is a sorted view, not the storage order.** Heights invert on
the birth-date coordinate, so a sweep walking nodes as stored admits them at the
wrong knots. The height coordinate keeps a sorted view for the same reason.

**One set of sums per species.** `eta` is a strategy parameter, so the pass is
knots plus cohorts per species, not across them.

**Value and slope come from the same sums in the same pass**, so they cannot be
two merges differing in their last bits.

---

## What this costs

**A re-blessing.** The re-association moves every number the model produces.
That is a `scientific_version` bump and a snapshot pass, not a silent change.

**The fast path is profile-specific.** It expands the smooth Yokozawa profile.
The soft box carries its own smoothstep and takes the fallback; the hard box
already refuses to build a light field.

**One digit of arithmetic**, from `1 - 2u + u^2` cancelling where `(1 - u)^2`
does not, and from the descent compounding over the knots. Measured within a
digit of the direct sum's own rounding on stands from 8 to 300 cohorts.

**The transpose is untouched and stays a product.** It is written per height,
so it remains correct — it transposes the same function — but the reverse pass
keeps its knots-times-cohorts cost. The rows factor the same way and a reverse
pass over the same recurrence would collapse it; that is separate work.

**The recorded block still declares `2K + layers` inputs**, at 480 slots per
knot. Knots are nearly free forward and still linear backward.

---

## What it delivers, measured

FF16, the same probe run against the same stands before and after:

| lifetime | cohorts | knots | build share | right-hand side | speed-up |
|---|---|---|---|---|---|
| 4 | 86 | 87 | 62% -> 16% | 2.40e-4 -> 9.5e-5 s | 2.5x |
| 10 | 93 | 178 | 77% -> 19% | 3.55e-4 -> 1.05e-4 s | 3.4x |
| 40 | 108 | 348 | 88% -> 24% | 7.95e-4 -> 1.25e-4 s | 6.4x |

A 40-year run goes from 0.60 s to 0.16 s. TF24's leaf solve is 1.8 ms against the
build's 0.26 ms, so the same change buys about 1.19x there — which is the point
rather than a weakness in it: the knot count was only ever expensive where the
build dominates, and on TF24 going from 166 to 315 knots already moves the
right-hand side by about 5 percent.

Agreement with the per-height reduction, on a field rebuilt at the state it is
read at:

| | value | slope | eps * A(0) |
|---|---|---|---|
| FF16, three stands | 2.8e-16 to 7.8e-16 | 6.7e-16 to 1.0e-15 | 3.2e-16 to 3.8e-16 |
| TF24 | 9.4e-16 | 1.8e-15 | 4.5e-16 |

And the ladder, against the branch this was cut from: tangent against a
difference 1.02e-08 on both, the transpose off the soil rates 5.97e-16 on both,
the transpose on the soil rates 3.35e-16 before and 4.23e-14 after, which is the
re-association and uses none of its budget.

## Three ways this went wrong, and what each needed to be seen

Recorded because the arithmetic was right from the first compile and every defect
was in the seam around it. Each was invisible to every check cheaper than the one
that caught it.

**The inflow boundary node was admitted unconditionally.** The per-height
reduction declines the closing interval when the node it would close from
contributes nothing, which happens when the newest cohorts carry no density yet.
Seeing it needs a stand in that state; on any other it is a fixed point of the
mistake.

**The reduction was not divided by patch area.** The per-height entry points
divide at the patch level, and filling the grid from the species directly walks
around that. Seeing it needs a patch whose area is not one — every other test
runs at the default, where dividing by one hides it.

**The cohort height was taken passively.** That deletes the field's height
channel. Seeing it needs a differentiated path: the value is identical, so the
non-gradient suite passes in full.

---

## The alternatives, and what eliminated them

Each was measured, not argued.

**Grade the lattice, spacing proportional to height.** At matched count it
starves the canopy top exactly as a uniform lattice starves the recruit — at
lifetime 10, 176 graded knots give the tallest cohort 4.16e-05 where 179 uniform
knots give it 1.26e-07. Half the cohort tops sit within a metre of the canopy,
so the structure there is set by the spacing between cohort tops and not by any
one cohort's height.

**A dyadic lattice, its level change smoothed by a blend.** A dyadic lattice is
a uniform lattice at whatever spacing its level picked, and matches one to three
digits (lifetime 40: dyadic 175 knots 2.24e-04, uniform 0.100 176 knots
2.24e-04). It buys a bounded count and pays field discontinuities for it. There
is nothing to blend for.

**Knots on the curvature breaks.** The breaks dominate the pointwise error — at
spacing 0.05 the spans holding a cohort top carry 2.2e-05 against 5.1e-07
elsewhere — but the crown integral averages that away, and at matched count the
placement is 0.4x to 0.9x, which is worse. It would also have required knot
positions to be a function of the state.

**Hold optical depth instead of light at the knots.** Measured at 1.00x on every
stand. The stand's optical depth is about 1.5, so Beer's law is not a strong
enough nonlinearity for the choice to matter.

**Store the grid from an earlier pass and reuse it.** This makes the
discretisation independent of the state by construction rather than by
`to_passive`, which is the right idea. The landed lattice already achieves it
and more cheaply: its positions are `k * spacing`, constants of every run rather
than of one, so they do not move under a trait perturbation either, and nothing
has to be stored per step.

---

## What would falsify this

- **The sweep and the per-height reduction disagree by more than the per-height
  reduction's own rounding**, on a stand the model reaches. The scale to beat is
  `eps * A(0)`; the check needs no reference.
- **The build is not the share it is measured to be** on a production
  configuration. The numbers here are FF16 and TF24 at three stand sizes.
- **A profile arrives that is neither a polynomial in `u^eta` nor a box.** Then
  the fallback is the production path and the sweep is dead weight.
- **The re-blessed numbers move by more than the re-association can explain.**
  A re-association moves the last bits; anything larger is a different change
  wearing its clothes.
