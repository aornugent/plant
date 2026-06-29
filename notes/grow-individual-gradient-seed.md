# Seed: FF16 trait gradient of `grow_individual_to_size` / `_to_height` — #472 scope B

Handoff for a fresh chat. Branch `spike-ff16-scm-emergent` (commit `c408babf` and
later). This is the **last FF16 surface without a trait gradient**; once it lands,
every FF16 entry point is gradient-enabled and the branch is a good PR boundary (Dan's
plan, 2026-06-30).

Read first: the guide `overstorey-staging/guides/autodiff-trait-gradients.qmd` (esp.
steps 3-5 "demographic rates / integrate through time / one plant to a stand" and "The
calibration interface" — the same `gradient as an optional output of the call you
already make` pattern). The coupled-replay sections are NOT needed here — this surface
has **no resident feedback** (the environment is given/fixed), so it is much simpler
than the SCM emergent work just finished.

## The task in one line

Given `grow_individual_to_size(individual, sizes, size_name, env, ...)` (FF16, a single
plant grown in a FIXED `env` up to target size(s)), return the derivative of its outputs
w.r.t. FF16 traits — primarily `d(state at target size)/d(theta)` and
`d(time-to-reach-size)/d(theta)` — as an optional output of the same call, exact via
reverse-mode AD, validated against finite differences.

## What this is (and why it's the easy one)

`grow_individual_to_size` (`R/individual.R:23`) integrates ONE plant's demographic ODE
in a fixed environment with the adaptive solver (`IndividualRunner` + `OdeRunner`), then
**bracket + bisect** (`grow_individual_bracket` / `grow_individual_bisect`) to find the
time `t*` at which `size_name` (e.g. height) hits each target, returning the ODE state
at `t*`. No canopy feedback, no cross-cohort coupling — it is the per-plant,
fixed-environment problem.

The gradient is a two-pass replay + an implicit-function step for the stopping time:
- **Pass 1 (double, the call as-is):** discover the time-to-size `t*` and the
  integrator's step schedule (the existing adaptive run).
- **Pass 2 (taped replay):** replay the demographic ODE over the frozen schedule with
  the trait active, to `t*`. Then
  `d(state at t*)/d(theta) = ∂state/∂theta|_{t*} + state_dt(t*) · d(t*)/d(theta)`,
  where `d(t*)/d(theta)` comes from the **implicit function theorem** on the stopping
  condition `height(t*, theta) = target`:
  `d(t*)/d(theta) = − (∂height/∂theta at t*) / height_dt(t*)`.
  (If a caller wants state-as-a-function-of-SIZE — e.g. growth rate, age, mass AT a
  fixed height — that is exactly the IFT term; for height itself the two cancel.)

## Reuse (most of the machinery exists)

- **`ff16_grow_demography<S>`** (`inst/include/plant/models/ff16_production_kernel.h:353`)
  already integrates the full FF16 demographic vector over fixed RK4 steps in a fixed
  **crown-top** light and reverse-differentiates `d(state at age t_end)/d(trait)`. This
  is the time-integrated gradient — the bridge to build on. It is exercised by
  `tests/testthat/test-ff16-grow-trajectory-ad.R` (a `Rcpp::sourceCpp` smoke test of
  `d(final height)/d(trait)`; skipped in `load_all` sessions, runs against an installed
  `plant.so`).
- **Deep-crown (the FF16 default).** `ff16_grow_demography` uses crown-top
  (`ff16_compute_rates_crown_top`). The default assimilation is deep-crown (moving-node
  GK integral). `deep_net<S>` / `deep_height_dt<S>` reading a **fixed** `FF16_Environment`
  already exist in `src/ff16_emergent.cpp` (the offspring path uses them); the
  grow-to-size kernel should read the given `env`'s light field the same way (a fixed
  light spline, no reconstruction). Decide up front whether to support crown-top only or
  deep-crown too (deep-crown is the realistic case; crown-top is the cheap first gate).
- **IFT pattern** for `d(t*)/d(theta)`: same shape as the `height_0` seedling-size IFT
  (`compute_dh0` in `src/ff16_emergent.cpp`) and the #539 leaf root-find — differentiate
  the converged stopping condition, never the bisection loop.
- **Adaptive vs fixed steps.** The SCM emergent path replays the adaptive Cash-Karp
  schedule faithfully (`ff16_cashkarp_replay`); `ff16_grow_demography` uses fixed RK4.
  For grow-to-size, the schedule is a single plant's adaptive run — harvest its step
  history (or, since Dan is happy to fix node/step schedules for gradients, run the
  discovery pass with a fixed step and replay that). Match whatever the double pass used
  so the replay reproduces `t*`/state to ~1e-12 before differentiating (the R0 gate).

## Suggested build order

1. **R0 gate:** a double-precision replay (fixed env, frozen schedule) reproducing
   `grow_individual_to_size`'s `state` and `time` at the target to ~1e-12. Crown-top
   first, then deep-crown over the fixed env.
2. **R1 crown-top:** reverse sweep for `d(state at t*)/d(theta)` and `d(t*)/d(theta)`
   (with the IFT stopping-time term), all 28 traits, vs a two-pass FD. Reuse
   `ff16_grow_demography` + add the IFT term + the `height_0` IFT for the seedling start.
3. **R1 deep-crown:** swap crown-top for `deep_net`/`deep_height_dt` over the fixed env.
4. **Wire a first-class R API + compile into `plant.so`** (like
   `offspring_production_gradient` / `stand_gradient`): e.g.
   `grow_individual_to_size_gradient(individual, sizes, size_name, env, traits=NULL)` →
   per-size × component × trait derivatives. Plain-R CI test (no sourceCpp), matching the
   established convention. Add a short guide subsection ("the same pattern on
   `grow_individual_to_size`").
5. **TF24** is a follow-up (the leaf opt re-solves per stage — use the harvested-profit
   tangent-linear replay from the TF24 emergent section of the guide); FF16 first.

## Pitfalls / notes carried over

- Build: C++-only → `cd <pkg root>` then `make compile`; new `[[Rcpp::export]]` or
  RcppR6 field → `make attributes` / `make RcppR6` then `make compile`; then
  `load_all(".", compile=FALSE)`. (Don't `cd src/` — the persisted shell cwd then breaks
  `make`.) Build optimized THEN `load_all(compile=FALSE)`; `load_all` alone is -O0.
- FD reference: a too-small step is corrupted by any hidden solver's noise floor;
  validate at the noise-optimal step (the recurring lesson — see the guide's "Lessons,
  distilled"). The bisection `t*` has its own tolerance; tighten it for the R0 gate.
- Reverse `adj` tape is odelia's single compiled global — SCOPE it (`{ }`); don't mix a
  live adj tape with a fwd pass.
- Run the WHOLE plant suite before the PR (generic tests loop over helper-plant.R
  strategy lists). Current state: 2474 pass / 0 fail / 9 skip (the skips are sourceCpp
  AD tests that need an installed `plant.so`, not `load_all`).

## Context: what's already done on this branch (the round-out just shipped)

FF16 resident gradient is complete and wired: `stand_gradient(feedback="resident")`
gives the coupled resident total for census metrics (single-species) AND the
cross-species total `d(total-stand metric)/d(theta_s)` for multi-species (fixed-schedule;
guarded error on stiff refined schedules — see the guide). Active birth-env is on by
default; the frozen-grid gap is decomposed (~half grid response, ~half birth-env).
`offspring_production_gradient`, `stand_gradient`, `stand_state_jacobian` all live for
FF16 (+ TF24 offspring). Status memory: `coupled-resident-gradient-status`. After
grow-individual lands → PR `spike-ff16-scm-emergent` (the FF16 AD surface is complete).
