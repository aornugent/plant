# Gradient (AD) benchmark + timing history — #472 scope B refactor+optimize

Timing history for the reverse-mode AD trait-gradient path across the refactor
(harvest → C++, engine unification). Convention follows the profile-plant skill
(`.claude/skills/profile-plant/SKILL.md`): **only same-session interleaved A/B
ratios are trustworthy** — machine drift across sessions is large, so never compare
absolute ms across dates. Every timing row carries its correctness verdict
(`fixture` column) — no speedup is recorded without its regression check.

## Tooling

- **`scripts/bench_gradient.R`** — same-machine A/B harness. Times three sub-costs
  separately so the harvest→C++ win is attributable:
  `run_ms` (resident `run_scm` solve, context) · `harvest_ms` (R-side `*_harvest()`,
  the refactor target) · `impl_ms` (C++ `*_impl` replay+sweep) ·
  `harvest_frac = harvest_ms/(harvest_ms+impl_ms)` · `val=<digest>` (Jacobian digest
  for same-session bit-identity). Env: `BENCH_REPS` (default 7), `BENCH_CASES`,
  `BENCH_SCALING=1`.
- **`scripts/gradient_fixture.R` + `tests/testthat/test-gradient-regression.R`** —
  the AD-vs-AD baseline. `gradient-baseline.rds` pins every engine's Jacobian;
  bit-identical (1e-12) for frozen/offspring/state-jac/grow/census/tf24-offspring,
  noise-floor (5e-6) for coupled/ms. Regenerate after an intentional change:
  `Rscript --no-init-file scripts/gradient_fixture.R snapshot`.
- **`scripts/profile-benchmarks.R`** — `PLANT_PROFILE_GRADIENT=1` adds
  `grad_ff16_frozen`, `grad_ff16_resident`, `grad_tf24f_census` cases (SCM built once
  outside the timed repeat). Native `/usr/bin/sample` ON for hotspot localisation,
  OFF for A/B ratios.

## A/B protocol (for each refactor step)

```sh
make compile                                                   # this branch
git worktree add -f /private/tmp/plant-base <pre-step-sha>
( cd /private/tmp/plant-base && make compile )
for r in 1 2; do
  Rscript --no-init-file scripts/bench_gradient.R "$(pwd)"            this-r$r
  Rscript --no-init-file scripts/bench_gradient.R /private/tmp/plant-base base-r$r
done
Rscript --no-init-file scripts/gradient_fixture.R check          # correctness gate
```

## Baseline (2026-06-30, `b3b4e188`, single-build indicative)

> Indicative only — `BENCH_REPS=3`, **not** an interleaved A/B (machine drift not
> controlled). Establishes order-of-magnitude and the headline `harvest_frac`
> finding; refactor-step rows below must use the interleaved protocol above.
> Apple Darwin 25.5.0, optimized `make compile` build.

| case | run_ms | harvest_ms | impl_ms | public_ms | harvest_frac | fixture |
|---|---:|---:|---:|---:|---:|---|
| ff16_frozen (28tr × 4met) | 850 | 85 | 1861 | 1911 | 0.044 | PASS |
| ff16_offspring (28tr) | 827 | 83 | 654 | 738 | 0.113 | PASS |
| ff16_resident_coupled (3tr × 3met) | 21 | 31 | 75 | 109 | 0.292 | PASS |
| ff16_resident_ms (2tr × 2met) | 29 | 29 | 175 | 284 | 0.142 | PASS |
| tf24f_census (4tr × 3met) | 39 | 47 | 108 | 154 | 0.303 | PASS |
| tf24f_resident (4tr × 2met) | 39 | 47 | 170 | 216 | 0.217 | PASS |
| grow_individual (28tr) | — | n/a | n/a | 8 | n/a | PASS |

### Reading of the baseline (sets the optimize priorities)

- **`harvest_frac` is 4–30%, not dominant.** The C++ `*_impl` (forward replay + one
  reverse sweep per metric) dominates wall-time in every case. So moving the R harvest
  into C++ (Phase 2) is primarily the **correctness** fix (kills the `Rcpp::as<>` env
  round-trip / the long-horizon fidelity floor); the wall-time it reclaims is modest
  (largest where the impl is cheapest: the coupled / TF24f cases at ~20–30%).
- **The impl is the wall-time lever.** `ff16_frozen` spends 1861 ms in the impl for
  28 traits × 4 metrics. Use `BENCH_SCALING=1` to split `impl ≈ replay(traits) +
  M·sweep`: if the per-metric sweep is cheap and flat, the forward-replay/tape build is
  the hotspot (the dedup in Phase 3, and any tape-scoping / quadrature work, target
  this); if it grows with M, the reverse sweep is. Localise with the native-sample
  section of `profile-benchmarks.R`.
- **Expectation after Phase 2:** `harvest_ms` → ~0; `public_ms` drops by roughly the
  old `harvest_frac` on each case. If it does not drop, the native sample confirms the
  time was never in the R harvest — a falsifiable check.

## Finding: the env round-trip is faithful (Phase 2 is a refactor, not a fidelity fix)

A lossy-vs-native A/B of the TF24f census recon (`tf24f_census_recon_impl` via the R
env list vs `tf24f_census_recon_native` reading `patch.environment_history` directly):

| H | SCM LAI | lossy recon (rel err) | native recon (rel err) |
|---|---|---|---|
| 4 | 1.61802617 | 6.26e-7 | 6.26e-7 (identical) |
| 5 | 1.59417777 | 3.98e-3 | 3.98e-3 (identical) |
| 8 | 1.65127320 | 1.19e-2 | 1.19e-2 (identical), no abort |

`native == lossy` bit-for-bit at every horizon ⇒ the `Rcpp::as<>` env round-trip is
**faithful** for the census path and is **not** the lifetime>4 floor. The floor is the
g' backward-FD near-cancellation + frozen-schedule replay. Phase 2a/2b are therefore
**bit-identical refactors** (remove a real round-trip + the env-validation hazard,
unlock perf headroom) — not a fidelity fix. Lifting the census floor = exact-AD g' /
refined-schedule replay (separate item). See memory `ad-env-roundtrip-faithful-for-census`.

Note: Phase 2a/2b nativized only the ENV extraction; the O(stand) `ppsurv` loop + `tw`
still run in R, so the `harvest_frac` perf win is **not yet realized** — that needs the
full native harvest (move `ppsurv`/`tw`/birth into C++), a follow-up.

## Refactor steps

| date | step (sha) | case | harvest_ms | impl_ms | public_ms | harvest_frac | cum × | incr × | sample | fixture |
|---|---|---|---:|---:|---:|---:|---:|---:|---|---|
| 2026-06-30 | baseline (b3b4e188) | — | see table above | | | | 1.00 | — | off | PASS |
| 2026-06-30 | 2a FF16 native-env (d69014b4) | ff16_frozen | (env only) | — | — | — | — | — | off | bit-identical |
| 2026-06-30 | 2b TF24f native-env (735d02c8) | tf24f_census | (env only) | — | — | — | — | — | off | bit-identical |
