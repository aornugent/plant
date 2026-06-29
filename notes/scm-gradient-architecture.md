# SCM / stand gradient architecture (#472 scope B, calibration-facing)

Design agreed 2026-06-29 (Dan + AD spike). This is the target shape for exposing
reverse-mode trait gradients of emergent stand outputs as a *calibration-ready*
capability, superseding the single-purpose `offspring_production_gradient()`.

## Principle: gradient is a capability *on the existing object*, parameterized by metric

Not a parallel `*_gradient()` family. The thing you already call gains an optional
gradient output. Three entry points, same pattern:

| Object / fn | Output differentiated | Use |
|---|---|---|
| `run_scm` / SCM | `d(metric)/d(traits)` per species, `metric ∈ {offspring_production, biomass, LAI, size-distribution moment, …}` | stand-level calibration / selection gradients |
| `Leaf` | `d(profit)/d(traits)` | standalone leaf calibration to data |
| `grow_individual_to_size` | `d(size-at-age / growth-rate)/d(traits)` | growth-rate optimisation / calibration |

## The engine lives OUTSIDE the SCM object

Decisive fact: the gradient is two-pass, and **pass 2 (the AD replay) never runs the
SCM solver.** It consumes only *harvested data* from a resident run:

- the frozen step schedule,
- the frozen per-RK-stage resident environment,
- each cohort's birth step / weight / survival.

So the engine is coupled to the *output* of a resident run, not to the solver. That
output is common across SCM variants (single-patch, multi-patch, stochastic,
continuous). Putting the engine inside `SCM<T,E>` would (a) couple it to a thing it
doesn't use, (b) force every new SCM variant to re-own/inherit it, and (c) drag the
strategy-specific AD kernels + XAD into the generic SCM TU.

Three **orthogonal axes** — keep them orthogonal by keeping the engine outside:

| Axis | Supplies | Examples |
|---|---|---|
| SCM variant | the harvest | single-patch, multi-patch, stochastic, continuous |
| Strategy | the rate kernel | FF16, TF24 |
| Metric | the scalar accumulator | offspring_production, biomass, LAI, size-dist moment |

A new SCM variant becomes differentiable by emitting the same harvest — nothing else.
It then gets every strategy × every metric for free.

### The seam: `ResidentHarvest`

Each SCM variant exposes a cheap `harvest()` returning a common struct (it already has
the pieces: `step_history`, `environment_history`, per-species cohort node times /
patch densities / survival-at-birth / `pr_survival`). The engine consumes the harvest;
ergonomics (`run_scm(collect_gradient = …)` or a free `stand_gradient(scm, …)`) are a
thin forwarding seam — UX does not require the engine to live in SCM.

## Metrics are generic and symmetric — no lock-in, none privileged

The capability must stay as generic/flexible as possible; no metric (not even
`offspring_production`) is privileged. Two layers give that:

1. **Generic reduction core.** Most stand metrics are a weighted reduction over the
   replayed cohorts, `metric = Σ_i w_i · f(state_i)` (optionally over censuses). The
   engine differentiates a *supplied* `(w, f)` — `offspring_production`, `LAI`,
   `biomass`, a size-distribution moment `Σ density_i·h_iᵏ` are all symmetric
   instantiations. `offspring_production` is just one registered entry, not a default.
   Adding a metric = adding a `(w, f)`; additive, no core change.
2. **Escape hatch for non-reduction metrics.** For anything that is not a simple
   reduction (quantiles, ratios, bespoke statistics a downstream package invents),
   expose the per-cohort **state × trait Jacobian** (heights, masses, density/survival,
   area_leaf, accumulator states, per saved census). Any smooth metric's gradient then
   composes downstream by chain rule — `plant` never needs to know the metric. Same
   boundary as "likelihoods live downstream."

## Many metrics from one baseline

Eventual use connects gradients to data and to *likelihoods*, with several likelihood
terms per stand. So the engine runs **multiple metrics against one baseline**:

- one pass-1 resident run → one harvest,
- one AD forward replay → one recorded tape,
- one cheap reverse (adjoint) sweep **per metric** over that tape.

Cost ≈ replay + M·(cheap sweep), not M replays. Primary output = a **metrics × traits
Jacobian** + the metric values. Reverse-per-metric is optimal while metrics ≤ traits
(the usual case); forward mode wins only when metrics ≫ traits (note, don't build yet).

## Boundary: plant returns the Jacobian; likelihoods are downstream

plant computes `d(metric_k)/d(trait)`. A downstream calibration package composes the
data/likelihood side `dL_k/d(metric_k)` by chain rule and stacks likelihood terms. Data
never enters plant — that is what lets many likelihood terms share one stand baseline.

```
SCM variant ──harvest──▶ [ gradient engine: metrics × traits Jacobian ] ──▶ downstream likelihood layer (other pkg)
   (any)                  (strategy kernel × metric, outside SCM)              dL/dmetric · dmetric/dtrait
```

## Build order

1. Generic reduction engine (reuse the existing replay): differentiates a supplied set
   of `(w, f)` reductions over one forward replay, metrics × traits Jacobian out, one
   reverse sweep per reduction. Seed the registry with `offspring_production`, `LAI`,
   `biomass`, size moments — all symmetric, none privileged. Plus the state × trait
   Jacobian escape hatch for downstream-defined metrics. FF16 first, TF24 mirror. No SCM
   templating, no `scm.h` change.
2. Formalize the `ResidentHarvest` seam so multi-patch/stochastic variants plug in.
3. `Leaf` and `grow_individual_to_size` gradient methods on the same pattern.
4. Phase-E tests + guide demonstrations: each entry point as a calibration/optimisation
   *setup* (no full calibration / likelihood / adaptive dynamics in plant itself).
