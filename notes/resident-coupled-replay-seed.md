# Seed: coupled-replay resident gradient (FF16 first) — #472 scope B

Handoff for a fresh chat. Branch `spike-ff16-scm-emergent`. Read this + the COURSE
CORRECTION in `notes/resident-gradient-scope.md` first, then `scm-gradient-architecture.md`.

## The task in one line

Build the **coupled** resident total-gradient of emergent stand metrics: all alive cohorts
stepped together through the frozen schedule, the canopy light an active function of the
**trait-dependent stand state** (canopy cohort heights `h_i` AND densities `density_i`
respond to θ, not just `area_leaf`). Output: `d(metric)/d(trait)` for all 28 FF16 traits,
metrics ∈ {LAI, biomass, size_moment} (offspring stays the frozen invasion gradient).

## Why (what the first cut got wrong)

R0–R1 shipped (`6a6b1ef0`) but over-froze: it held the background-canopy `h_i` and
`density_i` at θ₀ and let only `area_leaf(θ,h_i)` move — so only the leaf-area channel of
the feedback (C-27 sign flip), and only `a_l1`/`a_l2` got any feedback. The agreed design
(Dan, 2026-06-29) is that the canopy state itself responds to θ, so **every trait feeds
back** (a trait changes growth/mortality → heights/densities → the canopy everyone reads).
That requires re-evolving the stand under θ = the coupled replay.

## The design (freeze geometry = freeze the GRID, not the stand)

Freeze: the **ODE step schedule** (`step_history`), the **cohort-introduction schedule**
(birth steps, birth weights / patch density / pr_survival), and the **env-spline knot
x-positions** (the heights at which light is sampled — odelia #32's differentiable spline:
frozen knots, active values). Recompute: the **stand state** (every cohort's height,
log-density, demographic state) by replaying all cohorts together, and from it the **active
knot y-values** (light levels) at each RK stage.

Pass 2 is a **manual fixed-schedule `Patch<...,ad>` stepper** (a free function over a
`vector<Node<...,ad>>` + the frozen schedule — NOT `Patch` itself; its ctor eagerly runs
`compute_environment`/`reset`). Per Cash-Karp RK stage:
1. From the current ACTIVE stand, build the active competition spline at the FROZEN knot
   x-positions: `competition(z_k) = (1/area)ˑtrapezium_i( density_iˑk_Iˑarea_leaf(θ,h_i)ˑ
   Q(z_k/h_i; eta) )`, heights/densities/area_leaf all active. Fill an odelia
   `basic_interpolator<ad_t>` → `light(z) = exp(-competition(z))`, read by interpolation
   (matching `FF16_Environment::get_environment_at_height`).
2. Each cohort reads its crown light from that spline; compute its rates (deep-crown net →
   `ff16_compute_rates_from_net`), incl. the census `log_density` rate.
3. Introduce new cohorts at their frozen birth steps with active establishment.
4. Step all together (`ff16_cashkarp_replay` generalised to the whole-stand state).

At θ₀ this must reproduce the SCM's per-stage env + emergent metrics bit-closely (target
the mutant replay's ~3e-14; the first-cut value gate). No value-anchoring needed — the
spline is genuinely active.

## Cost / AD strategy (the key correction)

Build the canopy spline **ONCE per RK stage** and share it across cohorts (each reads via
O(log N) interpolation). Do NOT re-loop the whole stand inside each cohort's crown integral
at every GK point — that was the first cut's O(N²·steps·stages) trap that blew the reverse
tape to tens of GB and crashed. Built once per stage, the coupled tape is
~O(N·steps·stages) ≈ the frozen monolithic tape (~1.3s measured), so **ONE reverse sweep
per metric gives all 28 traits** (the scope's "coupled tape is affordable" point). Reverse
is the right tool here (many traits, few metrics) — unlike the first cut where the frozen
canopy decoupled cohorts and made forward-per-trait cheap.

Sanity cross-check: the coupled result's leaf-area part should reproduce the committed
first-cut sign flip (dLAI/da_l1: −0.53 frozen → +126 leaf-area-only); the coupled number
will differ by the height/density channel.

## What already exists (reuse)

- **Per-RK-stage harvest** (`inst/include/plant/patch.h`): `stand_height_stage_history`,
  `stand_competition_stage_history` ([step][stage 0..5][cohort], species 0), captured in
  `cache_RK45_step` aligned 1:1 with `environment_history`'s 6 per-stage envs; RcppR6-
  exposed. Useful for the θ₀ baseline check + extracting frozen knot x-positions. NOTE:
  species-0 only — generalise to all species (cross-species canopy) for multi-species.
- **`ff16_cashkarp_replay`** (`inst/include/plant/models/ff16_production_kernel.h:481`):
  generic Cash-Karp driver over a frozen step schedule, templated on State + deriv + axpy.
  Generalise State to the whole-stand vector.
- **`ff16_resident_light_at`** (same header, ~599): reconstructs resident light from a
  frozen stand active in traits (direct sum — use the TRAPEZIUM form to match the SCM).
- **deep-crown net kernels**: `ff16_net_from_components`, `ff16_compute_rates_from_net`,
  `ff16_area_leaf`, `ff16_canopy_q`, `ff16_assimilation_leaf`, `QK::integrate_ad<S>`.
- **Engine + API** (`src/ff16_emergent.cpp`, `R/emergent_gradient.R`): `stand_gradient(...,
  feedback=c("frozen","resident"))`, `ff16_stand_gradient_impl`, `ff16_harvest`. The
  resident branch currently does the leaf-area-only graft — replace its census path with
  the coupled replay (keep `feedback="frozen"` untouched).
- **odelia differentiable spline (#32)**: frozen knot x-positions + active knot values —
  the mechanism for step 1. Check the odelia version in use exposes `basic_interpolator`
  templated on the scalar; the env light spline is `ResourceSpline` over
  `AdaptiveInterpolator` (`inst/include/plant/resource_spline.h`,
  `inst/include/plant/models/ff16_environment.h:150` `compute_environment`).

## Validation plan

- **R0 gate**: coupled DOUBLE replay reproduces SCM per-stage env + emergent metrics to
  ~3e-14 at θ₀ (extract frozen knot x-positions from `environment_history`; compare the
  active-spline values at θ₀ to `get_environment_at_height`).
- **R1**: coupled reverse sweep → resident dLAI/dθ etc. for ALL 28 traits; validate vs a
  two-pass FD that **re-runs the coupled reconstruction** (perturb θ, re-evolve the coupled
  stand, re-reduce) — NOT a full SCM re-run. Expect all 28 nonzero (vs 2 in the first cut).
- Reuse `scripts/ad_resident_gradient.R` as the harness skeleton.

## Pitfalls learned (first cut)

- The reverse-mode `adj` tape is odelia's single compiled global; SCOPE it (`{ }`) so it's
  destroyed before any forward (`xad::fwd`) pass, and don't nest it with `compute_dh0`'s
  tape concerns. Mixing a live adj tape with fwd ops caused nondeterministic heap crashes.
- Re-looping the whole stand per GK point per cohort = O(N²) tape blow-up → crash. Build
  the canopy once per stage.
- `nohup … &` returns the echo's exit code, not Rscript's; a segfault shows as truncated
  output. Append a sentinel (`echo FINISHED`) and wait for it.
- RcppR6 objects in `environment_history` hold external pointers — can't `saveRDS` the
  harvest across processes; re-run the SCM in-process.

## Build / test workflow

- Interface change (RcppR6_classes.yml or new `[[Rcpp::export]]` arg): `make rebuild`
  (or `Rscript -e 'Rcpp::compileAttributes()'` then `make full_compile`).
- C++-only: `make compile`. Then `devtools::load_all(".", compile=FALSE)`.
- Memory note: `make compile`/`full_compile` build optimized; `load_all` alone is -O0.
- Full test suite before PR (generic tests loop over helper-plant.R strategy lists).
