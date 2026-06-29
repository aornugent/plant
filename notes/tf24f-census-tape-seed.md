# Seed: TF24f reverse-mode AD census tape (the refine step) — #472 scope B

Handoff for a fresh chat. Branch `spike-ff16-scm-emergent` (commits `deef6823` →
`b331b7d9`). The TF24f gradient **prototype phase is complete** — every surface has a
working, validated FD gradient (or exact AD for `g'`). This seed is the **refine step**:
replace the per-trait FD census gradient with ONE reverse-mode AD sweep over the replay.

Read first: `notes/tf24-stand-gradient-scope.md` (the whole thing; §8-§9 are this
session's findings) and the TF24 section of
`overstorey-staging/guides/autodiff-trait-gradients.qmd` (the harvest + linearised-profit
tangent-linear replay — TF24's net production has no adjoint tape, so you harvest the
converged leaf operating point and replay a tapeable expression). The shipped TF24
offspring tape `src/tf24_emergent.cpp` (`tf24_stand_gradient_impl`, `tf24_replay_full`,
`harvest_at`, `inj`) is the structure to extend.

## The task in one line

Return a metrics×traits Jacobian `d(census metric)/d(theta)` (LAI / biomass /
size_moment) for a **frozen-resident** (rare-mutant) TF24f stand in one reverse sweep,
matching `tf24f_census_gradient_fd` (the FD gate already built) to ~1% (the recon noise
floor). Then the resident/coupled and individual tapes follow.

## What exists to build on (this session's deliverables)

- **R0 census recon (the faithful replay):** `tf24f_census_recon_impl`
  (`src/tf24f_emergent.cpp`) — a double-precision 7-state replay {5 demog, tracked collar
  `opt_root_psi_state`, `log_density`} over the SCM's frozen Cash-Karp schedule, driving
  the REAL TF24f leaf at the tracked collar. Reproduces the SCM stand bit-faithfully at
  patch lifetime ≤4. The tape replaces the real-leaf calls with a harvested + linearised
  expression (as TF24 offspring does), but must reproduce THIS replay.
- **FD gates (the validation targets):** `tf24f_census_gradient_fd` (frozen) and
  `tf24f_resident_census_gradient_fd` (resident total, via run_scm re-runs) in
  `R/tf24f_emergent_gradient.R`; tests in `test-tf24f-census-gradient.R`.
- **Exact-AD `g'`:** `TF24f_Strategy::growth_rate_gradient_height_ad` + `Leaf::dprofit_dPPFD`
  (`src/tf24f_strategy.cpp`, `src/leaf_model.cpp`); `node_gradient_exact_ad` now exposed in
  the R `Control`. Validated native-correct.

## The one hard ingredient (why this is the heavy refine, not a copy of TF24 offspring)

**TF24f's tracked collar is strongly θ-dependent, and at `k_acclim=1` it lags the optimum
— so the envelope theorem does NOT zero its contribution** (it does for TF24, which is why
TF24 offspring is easy). Profit depends on the collar; `∂profit/∂collar = dprofit_dpsi ≠ 0`.
So you must carry the collar as a **tapeable state** with rate `k_acclim·dprofit_dpsi`, and
that rate's θ-response needs **second derivatives of profit** (a "curvature harvest"):

- collar rate linearised: `dprofit_dpsi ≈ dprofit_dpsi0 + (d²profit/dpsi²)(collar−collar0)
  + (d²profit/dpsi dh)(h−h0) + Σ_k (d²profit/dpsi dθ_k)(θ_k−θ_k0)`.
- profit (for the demographic rates) linearised in (h, collar, θ): the collar term
  `dprofit_dpsi0·(collar−collar0)` is what the offspring tape omits — add it (collar is taped).
- `log_density` rate `= −g' − mort_dt`; `g' = d(height_dt)/dh` — its θ-derivative needs
  `d²profit/dh²` and `d²profit/dh dθ` (the same curvature obstacle the guide flags for the
  TF24 census). Either harvest these too, OR reuse the exact-AD `g'` machinery
  (`growth_rate_gradient_height_ad`) and harvest its θ-sensitivity.

Harvest the second derivatives by **FD in the double pass** (perturb collar / h / θ_k,
recompute `leaf.dprofit_droot_collar_psi` and the profit, difference) — no leaf templating
needed. Cost: the `d²profit/dpsi dθ_k` term is ~one leaf-gradient eval per trait per stage;
for the prototype validate a few traits, not all 27. (Alternative, cleaner-but-bigger:
template the leaf-at-collar eval on the XAD type so the collar propagates analytically with
no curvature harvest — the leaf is currently a plain-double class with hand-coded
sensitivities, so this is a real leaf-model project.)

## Suggested build order

1. **Curvature harvest, validated standalone.** In the double pass, harvest per (cohort,
   stage) at the tracked collar: `{profit0, dprofit_dh, dprofit_dθ_k, dprofit_dpsi0,
   d²profit/dpsi², d²profit/dpsi dh, d²profit/dpsi dθ_k, d²profit/dh²}`. Validate each FD
   term against an independent FD before taping (the recurring lesson — pin the harvest
   first).
2. **Tape the 7-state replay** (mirror `tf24_replay_full` + add collar + log_density;
   census reduction from `ff16_emergent.cpp`'s `census_reduce` / trapezium + new_node
   tail, already mirrored in `tf24f_census_recon_impl`). One reverse sweep per metric.
   Validate `d(census)/dθ` vs `tf24f_census_gradient_fd` to ~1%.
3. **Wire `stand_gradient(..., strat="TF24f")`** census path (today `stand_gradient` does
   FF16 + TF24 offspring; add TF24f census). First-class R API, plain-R CI test.
4. **Individual grow-to-size AD tape** — same collar curvature harvest, no canopy/density;
   validate vs `tf24f_grow_individual_to_size_gradient_fd`.
5. **Resident/coupled AD** — swap the FF16 coupled per-cohort rate call for the TF24f
   harvested-leaf-at-tracked-collar eval inside the existing coupled engine; validate vs
   `tf24f_resident_census_gradient_fd`. Fixed node schedule (TF24f is stiff).

## Pitfalls / notes carried over

- **`Rcpp::as<>` env round-trip is NOT faithful** for the crown-sampled light above cohort
  heights. It bit this session twice (a "validation trap": an exact-AD `g'` that was
  correct natively read as 18% wrong through an `as<>` env; and it caps the census recon's
  long-horizon fidelity). **Validate env-dependent quantities INSIDE run_scm (native env),
  never through an R→C++ env round-trip.** Longer-term fix = move the replay/recon into C++
  driving native envs (memory: `move-gradient-machinery-to-cpp`).
- **Fidelity floor / horizon:** the R0 gate is bit-faithful at lifetime ≤4; lifetime 5
  recon LAI drifts ~0.4% (the `as<>` round-trip, NOT `g'` or the collar — both are sound).
  Keep gates at lifetime ≤4.
- **Drought shutdown:** `TF24f::solve_leaf`'s AD acclimation gradient now gates on
  `prepare_collar_solve` (returns 0 at shutdown) — fixed the lifetime≥8 abort. Any new code
  that calls `Leaf::dprofit_droot_collar_psi` must respect the same shutdown gate (the
  transport-spline derivative extrapolates at `psi_crit`).
- **Build:** `cd <pkg root>`; C++-only → `make compile`; new `[[Rcpp::export]]` →
  `make full_compile` (runs compileAttributes); new RcppR6 field → `make rebuild`; then
  `load_all(".", compile=FALSE)`. Don't `cd src/`. Build optimized THEN load_all(compile=FALSE).
- **Reverse `adj` tape** is odelia's single compiled global — SCOPE it; tape one cohort at
  a time where independent (offspring is per-cohort; census couples cohorts via the
  trapezium → one tape over all cohorts, then one reverse sweep per metric).
- Run the WHOLE plant suite before any PR (generic tests loop over helper-plant.R; TF24f is
  intentionally excluded from those generic lists — it is a 6-state variant).

## Context: prototype status (what shipped this session)

TF24f now has FD prototypes for frozen census, resident (total) census, and individual
grow-to-size, plus exact-AD `g'`; the resident gradient demonstrably flips the sign of
`d(LAI)/d(lma)` (+13.6 frozen → −2.3 resident) — the canopy feedback. Status memory:
`tf24f-census-r0-gate`. The refine phase is the AD tapes (this seed) + the C++ native-env
migration.
