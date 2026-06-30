# AD gradient machinery: refactor, consolidate, optimize

## Context

Reverse-mode automatic differentiation of emergent trait gradients is **built and
validated end-to-end** for FF16, TF24, and TF24f on branch `spike-ff16-scm-emergent`
(101 commits, ~21,700 net lines over `develop`). The prototype phase is done; this plan
is the **refactor + optimize phase**. It does **not** add modelling capability — it
makes the validated machinery faster, smaller, more correct, and CI-covered.

Three problems motivate it:

1. **The R-side `Rcpp::as<>` env round-trip is the correctness ceiling.** The gradient
   is two-pass: pass 1 runs a resident `run_scm(save_RK45_cache=TRUE)` and harvests a
   frozen schedule + per-RK-stage environment; pass 2 is an AD replay that never re-runs
   the solver. Pass-1 harvest happens **in R** (`ff16_harvest`, `tf24f_harvest`, …),
   pulling `patch$environment_history` across the boundary and rebuilding each env via
   `Rcpp::as<>`. That reconstruction is **not faithful** for the crown-sampled light
   above cohort heights — it is the root cause of the TF24f census-recon fidelity floor
   at patch lifetime > 4 and of a validation trap (a native-correct AD `g'` read as "18%
   wrong" only through a round-tripped env). It also rebuilds the whole RcppR6 patch on
   every access (O(stand size), ~1600× slower per `R/emergent_gradient.R:70-73`).

2. **Duplication.** Three near-identical replay engines (`src/ff16_emergent.cpp` 1851,
   `src/tf24_emergent.cpp` 479, `src/tf24f_emergent.cpp` 2266) plus three R harvest
   seams. The Cash-Karp stepper and state containers are *already* shared; what is
   triplicated is the **orchestration glue** (the `Frozen` struct, `build_frozen`/`as<>`
   marshalling, per-strategy deriv/axpy/seed, `census_reduce`, coupled-canopy
   reconstruction).

3. **No benchmarks or CI coverage for the AD path.** There is no timing harness for the
   gradient calls, and 7 FF16 AD test files `skip()` with "AD tape symbols unavailable
   in this load_all session" — so the FF16 AD surface is not exercised in CI. ~31
   intermediate `ad_*`/`coupled_*`/`recon_*` validation scripts accumulated during the
   prototype.

**Intended outcome:** the AD path drives native environments (round-trip gone, the
long-horizon fidelity floor dissolved), one strategy-parameterized engine instead of
three, a benchmark + regression net that proves every step preserves the validated
Jacobians, and the spike merged to `develop` as a clean reviewable PR stack.

### Decision on PR #541 / merge sequencing (resolved with Dan)

**Do NOT gate this work on landing #541 first; close #541 without merging.** #541
(`spike-ff16-hierarchy`, the `<T,E,S>` templating) is already the ancestor base of this
spike — nothing depends on it being in `develop`, and its commits remain in the spike
history whether or not the PR is open. The only value of keeping it open was reviewer
ergonomics, which is better served by re-cutting fresh PRs from the *final* tree.
**Action: close PR #541 (do not merge)** — it isn't serving a purpose as a live review
object. We do all the work in place on the spike (no rebase, no merge conflict), then
decompose the final tree into a stacked PR sequence for `develop` (Phase 5), where a
foundation slice (the templating) can be cut fresh as PR #1. *(Closing the PR is a
mutating GitHub action — done on plan approval, not during planning.)*

---

## Phase 1 — Safety net: benchmarks + correctness-regression fixture (do first)

Establish the baseline *before* touching anything, so every later step is provably
non-regressing in both timing and value. Templates the existing
`scripts/bench_ab.R` / `scripts/bench_tf24.R` conventions.

- **`scripts/bench_gradient.R`** — same-machine A/B harness (mirrors `bench_ab.R`:
  `path`+`label` args, `pkgload::load_all(compile=FALSE)`, warm-up, one `RESULT|…` line
  per build, interleave this-branch vs a worktree). The novel requirement: **time the
  three sub-costs separately** by calling internals directly —
  - `run_ms` = resident `run_scm(save_RK45_cache=TRUE)` (context),
  - `harvest_ms` = the R-side `*_harvest()` call alone (**the refactor's target**),
  - `impl_ms` = the C++ `*_impl` replay+sweep alone,
  - headline `harvest_frac = harvest_ms/(harvest_ms+impl_ms)` and a `val=<digest>`
    bit-identity check (à la `bench_tf24.R:31`).
  Canonical cases lifted from test fixtures (fixed schedules): `ff16_frozen` (28 traits,
  4 metrics), `ff16_resident_coupled`, `ff16_resident_ms`, `ff16_offspring`,
  `tf24f_census`, `tf24f_resident`, `grow_individual`. `BENCH_REPS`/`BENCH_CASES` env
  vars; optional `BENCH_SCALING=1` to factor `impl_ms ≈ replay(traits) + M·sweep`.

- **Correctness-regression fixture (AD-vs-AD baseline).** Snapshot the current validated
  Jacobians (values + `d(metric)/d(trait)` matrices) to
  `tests/testthat/fixtures/gradient-baseline.rds` for every distinct engine
  (frozen / coupled-single / coupled-ms / state-jacobian / offspring / grow, across
  FF16/TF24/TF24f). A new `tests/testthat/test-gradient-regression.R` reloads it and
  asserts two-tier tolerance: **bit-identical** (`1e-12` rel) for pure-relocation steps,
  **noise-floor** (`~5e-6` rel) for the coupled/ms paths that legitimately reorder FP
  sums. This is complementary to the existing AD-vs-FD tests (those pin AD to physics at
  ~1%; this pins AD to its own validated self at machine precision). Use fixed node
  schedules so snapshot/check are reproducible.

- **Extend `scripts/profile-benchmarks.R`** with gradient cases (`grad_ff16_frozen`,
  `grad_ff16_resident`, `grad_tf24f_census`) reusing its Rprof + native `/usr/bin/sample`
  machinery. Per the profile-plant skill: sample **ON** for hotspot localisation, **OFF**
  for A/B ratios; only same-session ratios are trustworthy.

- **Timing history:** `notes/profile-gradient-2026-06-30.md`, columns
  `date | step+sha | case | run_ms | harvest_ms | impl_ms | public_ms | harvest_frac |
  cum_speedup | incr_speedup | sample | fixture(PASS/FAIL + worst rel_dev) | notes`.
  No speedup row is recorded without its correctness verdict attached.

**Validates:** the fixture itself is the validation infrastructure; confirm it passes on
the current HEAD and that `harvest_ms + impl_ms ≈ public_ms`.

---

## Phase 2 — Native-env harvest in C++ (kills the `Rcpp::as<>` round-trip)

The highest-value change: it has a **correctness payoff** (removes the lifetime-floor and
the validation trap) and shrinks the per-engine surface Phase 3 then unifies. The
round-trip is structurally avoidable — confirmed: the replay reads the env only as a
*double* value + analytic derivative (`src/ff16_emergent.cpp:125-128`), lifting just those
scalars into the AD type; the AD type never touches the env. A faithful native env
**pointer** is sufficient and exact.

**Approach (honours `notes/scm-gradient-architecture.md`: engine stays OUTSIDE the SCM
object).** Do **not** add a `collect_gradient` run mode to `SCM::run()` — the harvest is
already captured during the normal `save_RK45_cache=TRUE` run. Instead:

- Add `inst/include/plant/gradient/resident_harvest.h` with a `ResidentHarvest` struct
  and a free `harvest_resident(const Patch<T,E,S>&, int species)` builder that gathers,
  **by const-ref/pointer into the Patch's own storage**, the pieces R rebuilds today:
  `step_history`, `&environment_history[n][s]` (native pointers, not copies),
  `stand_*_stage_history`, `stand_newnode_*`, and per-species cohort metadata. Move the
  two remaining R computations — trapezoid node weights (`emergent_gradient.R:94-97`) and
  the per-RK-stage `pr_survival` matrix (`:99-103`) — into this C++ builder.
- Change each engine's `Frozen.eh` from owning `std::vector<…Environment>` to holding
  `const Environment*` borrowed from the Patch. **Delete the `Rcpp::as<>` loop**
  (`ff16_emergent.cpp:450`, and the TF24/TF24f equivalents). This line is the round-trip.
- The C++ entry takes the **live `Patch&`** from the RcppR6 external pointer (env pointers
  are only valid while the SCM C++ object lives; RcppR6 keeps it alive in R). The R API
  (`stand_gradient()` etc.) becomes a thin forwarder passing the live SCM xptr down.

**Sub-steps, each bit-checkable against Phase 1's fixture:**
1. **2a** — `ResidentHarvest` + C++ gather + the two moved arithmetic computations, still
   building `Frozen` from copied envs. **Bit-identical** to current `ff16_harvest`
   outputs (assert `tw`, `ppsurv`, sampled env values equal). Isolates the gather move.
2. **2b** — switch `Frozen.eh` to borrowed native pointers; delete `as<>`. Bit-identical
   at lifetime ≤ 4; at **lifetime 5–8** the recon LAI must now match
   `compute_competition(0)` to ~1e-6 and TF24f must no longer abort at the spline edge —
   add new long-horizon gates to `test-tf24f-census-gradient.R` that were previously
   impossible. **This is the headline correctness win.**
3. **2c** — repeat for TF24 / TF24f env entries.

**Files:** new `inst/include/plant/gradient/resident_harvest.h`; `inst/include/plant/patch.h`
(public accessors already exist, `:153`, `:176-218`); `src/ff16_emergent.cpp`,
`src/tf24_emergent.cpp`, `src/tf24f_emergent.cpp` (`Frozen`, `build_frozen`, `attach_*`);
`R/emergent_gradient.R`, `R/tf24_emergent_gradient.R`, `R/tf24f_emergent_gradient.R`
(harvests shrink to forwarders). `make rebuild` if an export signature changes.

**Validates:** all `test-{ff16,tf24,tf24f}-*gradient*.R` stay green; the new lifetime-5/8
TF24f gates pass; the Phase-1 fixture confirms ≤4-horizon Jacobians unchanged; record the
`harvest_frac` drop in the timing-history file.

---

## Phase 3 — Unify the three replay engines (maximal dedup, per Dan)

Collapse the triplicated orchestration into one strategy-parameterized engine. The
stepper (`ff16_cashkarp_replay`) and state containers are already shared across all three
(`tf24_emergent.cpp` and `tf24f_emergent.cpp` already `#include ff16_production_kernel.h`).
The engines differ in `replay_cohort_final` at exactly three points: the env→net call, the
state-vector width (5/6/7), and the rate fill.

**Minimal trait-class interface** each strategy supplies (mostly *binding* things already
templated — `FF16ProdPars<S>`/`TF24ProdPars<S>` and `*_compute_rates_from_net` already
exist with matching shapes):
```
struct StrategyTraits<Strat,S> {
  using ProdPars  = …;            // existing
  using LifeState = …;            // strategy-supplied width (TF24f adds collar + log_density)
  static S        net_at(ProdPars, Env*, LifeState, /*harvest@n,stage*/);
  static Rates    rates_from_net(ProdPars, height, area_leaf, net);
  static LifeState seed(ProdPars, Env* birth_env, h0, …);   // establishment / IFT collar birth
};
```
TF24f's collar-as-6th-state + curvature harvest fits as a **hook**, not a fork: the collar
is an extra `LifeState` component with rate `k_acclim·dprofit_dpsi`; the curvature harvest
is a *harvest producer* orthogonal to the replay template. **Do not** merge TF24's
envelope-zeroed collar with TF24f's taped collar — they are different math; keep them as
two `seed`/`net_at` implementations behind the hook.

**Ordered dedup (each gated by Phase-1 fixture); Dan chose to push through all ranks:**
1. **Shared `census_reduce` + `(w,f)` metric dispatch** — pure arithmetic, currently
   hand-copied identically (`ff16:587`, `tf24f:226`). **Bit-identical.** (~300 lines)
2. **Generic `replay_cohort_final<Strat,S>`** via the trait class (frozen/invasion path).
   FF16/TF24 offspring are bit-exact references. **Bit-identical.** (~400-600 lines)
3. **Generic coupled-canopy reconstruction** — the Yokozawa trapezium
   `Q=(1-(z/h)^eta)^2` is identical across strategies (`coupled_comp_at`/`build_interp`
   vs `tf24f_comp_at`/`build_canopy`); only the per-cohort rate differs. Validated to the
   ~1e-4 FD-noise floor + sign-flip invariants (not bit-identical — the coupled tape is
   numerically delicate). (~500-800 lines)
4. **Multi-species engine** (`assemble_metrics_coupled_ms`/`FrozenMS`). **Highest risk** —
   the MS stiffness/conditioning guards (joint-env-drift > 1e-3 gate, finiteness guard,
   tape scoping) must be carried through **verbatim**; a refactor that drops them silently
   re-opens divergence. Validate against the MS-R0/MS-cross-species-R1 gates on fast fixed
   schedules.

**Files:** `inst/include/plant/models/ff16_production_kernel.h` (shared `census_reduce`,
generic `replay_cohort_final`, trait base), new
`inst/include/plant/gradient/replay_engine.h` (strategy-agnostic orchestration); the three
`src/*_emergent.cpp` shrink to: include the generic engine + supply the strategy trait
class + the `[[Rcpp::export]]` entry points.

**Validates:** ranks 1-2 bit-for-bit; rank 3 at ~1e-4 + sign-flips; rank 4 against the MS
gates. Run the **whole** test suite each step (generic tests loop over
`helper-plant.R` strategy lists; remember TF24f is the non-5-state variant excluded
there). Record line-count reduction + timing per step.

---

## Phase 4 — Cleanup: fold validation scripts into tests, cover the AD surface in CI

Per Dan: **fold the still-valuable validation scripts into proper tests**, then delete the
scripts (rather than deleting outright or leaving them as loose scripts).

- **Migrate** the meaningful checks from `scripts/coupled_r0_check.R`,
  `coupled_r1_check.R`, `recon_static_check.R`, `coupled_vs_full_scm_fd.R`,
  `coupled_multispp_r0.R`, `coupled_multispp_r1.R`, `coupled_gap_decompose.R`,
  `coupled_birthenv_channel.R`, etc. into `testthat` tests (or fold into the existing
  `test-{ff16,tf24f}-*gradient*.R`). Delete the migrated scripts. Keep a small curated set
  of `ad_*` demos as runnable examples; the `overstorey-staging` guide `.qmd` remains the
  user-facing narrative.
- **Make the AD tests actually run in CI (enabler — otherwise folded tests just skip).**
  The 7 FF16 AD test files skip with "AD tape symbols unavailable in this load_all
  session"; TF24f tests already run live. Diagnose the asymmetry and make the FF16 AD
  symbols available under the test load path — tie to the **odelia DLL load-ordering**
  fix (plant needs `importFrom(odelia,…)` so odelia's `.onLoad` runs before plant's
  `useDynLib`) and the build-optimized-then-`load_all(compile=FALSE)` path. The aim:
  folded + existing AD tests execute in CI rather than skip. (Aligns with roadmap Phase C:
  "make the headline result CI-testable in plain R".)

**Validates:** the migrated tests pass under `make test` / `devtools::test()` *without*
skipping; the AD surface is genuinely covered. Confirm no orphaned `source(scripts/…)`
references remain.

---

## Phase 5 — Decompose for merge to `develop`

PR #541 is **closed (not merged)** at the start of this work — its commits stay as spike
ancestors. Re-cut the final spike tree into a stacked, reviewable PR sequence (this is
*merge* decomposition, independent of the in-place writing order above):

1. **Foundation** (`<T,E,S>` templating + exact AD growth-rate gradient) — additive,
   bit-identical to current package; cut fresh from the final tree, the clean foundation
   reviewers vet first (replaces the role the old #541 would have played).
2. **Native-env harvest** (Phase 2) — the correctness fix; land even if later phases slip.
3. **Engine unification** (Phase 3 ranks 1-2 bit-identical; rank 3 noise-gated; rank 4 MS).
4. **Cleanup + CI coverage** (Phase 4).

Each PR carries its slice of the timing-history table and its fixture verdict.

---

## Verification (end-to-end)

- **Build discipline:** C++-only change → `make compile`; new export / RcppR6 field →
  `make rebuild`; then `pkgload::load_all(".", compile=FALSE)` (optimized -O2, **never**
  `devtools::load_all()` which is -O0). Run the **whole** plant suite before each PR
  (`make test`), not just `test-strategy-*` — generic tests loop over `helper-plant.R`.
- **Correctness gate (every step):** `Rscript scripts/gradient_fixture.R check` (or the
  `test-gradient-regression.R` testthat run) — bit-identical for relocation/rank-1-2
  steps, noise-floor for coupled/ms. Plus the existing AD-vs-FD tests stay green.
- **The Phase-2 headline proof:** new lifetime-5 and lifetime-8 TF24f census gates that
  were impossible under the round-trip now pass to ~1e-6 vs `compute_competition(0)`.
- **Timing (same-session A/B):** build this branch + a worktree of the pre-step ref,
  interleave `Rscript scripts/bench_gradient.R <path> <label>` runs, compare `RESULT`
  lines; expect `harvest_frac` → ~0 and `public_ms` to drop by the old harvest fraction on
  the frozen/offspring/grow cases. Native `/usr/bin/sample` (OFF for ratios) confirms time
  left the R harvest. Record cumulative + incremental speedup in
  `notes/profile-gradient-2026-06-30.md`.

## Critical files

- `R/emergent_gradient.R` (`ff16_harvest` :61, `ff16_harvest_ms` :129, the O(stand) comment
  :70-73, `stand_gradient` :229) — the R-side round-trip + arithmetic being moved to C++.
- `R/tf24_emergent_gradient.R`, `R/tf24f_emergent_gradient.R` — TF24/TF24f harvests + FD gates.
- `src/ff16_emergent.cpp` (reference engine: `Frozen` :66-117, `build_frozen` :442-462,
  `replay_cohort_final` :225-256, coupled :643-921), `src/tf24_emergent.cpp`,
  `src/tf24f_emergent.cpp` (7-state collar replay :71, curvature harvest, coupled canopy).
- `inst/include/plant/patch.h` (native `environment_history` :153, `stand_*_stage_history`
  :176-218 — the harvest source).
- `inst/include/plant/models/ff16_production_kernel.h` (shared stepper/state — where the
  generic engine + trait class live); `tf24_production_kernel.h`.
- New: `inst/include/plant/gradient/resident_harvest.h`, `…/gradient/replay_engine.h`.
- New: `scripts/bench_gradient.R`, `tests/testthat/fixtures/gradient-baseline.rds`,
  `tests/testthat/test-gradient-regression.R`, `notes/profile-gradient-2026-06-30.md`;
  extend `scripts/profile-benchmarks.R`.
