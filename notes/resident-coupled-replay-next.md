# Seed: rounding out the FF16 coupled resident gradient — #472 scope B (R1+)

Handoff for a fresh chat. Branch `spike-ff16-scm-emergent`. The coupled resident
replay (R0 + R1) is **built, validated, committed, and wired into the public API**.
Read `notes/resident-coupled-replay-seed.md` (the original design) and the guide
section "The resident total gradient: the coupled replay" in
`overstorey-staging/guides/autodiff-trait-gradients.qmd` first.

## STATUS UPDATE (2026-06-29 round-out session): steps 1–6 addressed

Steps 1, 2, 3, 5, 6 are **done**; step 4 (multi-species) is **substantially done**
(harvest + engine + cross-species gradient validated on fixed schedules; the stiff
*refined* production schedule has an R0 robustness gap, so the public API still
guards multi-species to `feedback="frozen"`). Detail per step at the bottom of this
file under "ROUND-OUT OUTCOMES". The original step list and design notes below are
kept for context.

## What landed (commits b8bad961 → 0a5d073e)

The course-corrected **coupled** replay: all alive cohorts re-evolved together over the
frozen schedule; canopy light reconstructed each RK stage from the ACTIVE stand
(heights AND log-densities respond to θ) via an odelia `basic_interpolator<S>` at the
frozen knot x-positions. Built once per stage, shared across cohorts.

- **Harvest** (`patch.h`): boundary `new_node` per RK stage
  (`stand_newnode_{height,competition}_stage_history`), snapshotted INSIDE
  `compute_environment` before `compute_rates` mutates it. RcppR6-exposed.
- **Engine** (`src/ff16_emergent.cpp`): `coupled_comp_at` / `deep_net_coupled` /
  `deep_height_dt_coupled` / `rkck_one_step` / `assemble_metrics_coupled` /
  `attach_coupled`. Exports `ff16_coupled_metrics_impl` (R0 double) +
  `ff16_coupled_gradient_impl` (R1 reverse).
- **API** (`R/emergent_gradient.R`): `stand_gradient(feedback="resident")` routes the
  census metrics (LAI/biomass/size_moment) to the coupled engine; offspring stays
  frozen (invasion). `ff16_harvest` returns `nn_h`/`nn_c`.
- **Tests**: `test-ff16-stand-gradient.R` has a coupled-resident case (AD vs
  coupled-recon FD; offspring stays frozen). Scripts: `recon_static_check.R`,
  `coupled_r0_check.R`, `coupled_r1_check.R`.

**Validated**: static env recon 4.4e-16 (machine-exact); R0 metrics 1e-8…1e-10, env
drift ~7e-6 at z=0; R1 AD-vs-coupled-recon-FD ~1e-4 (noise-optimal step). 27/28 traits
feed back (vs 2 frozen-geometry first cut); sign flips. 28 traits × 3 metrics in 3.9s.

## Q: does it work with multiple species? — NO (single-species only)

The per-RK-stage harvest (`stand_*_stage_history`, `stand_newnode_*`) is **species[0]
only** (see `patch.h:163,175`). `attach_coupled` reconstructs competition from the
re-evolved cohorts of the ONE differentiated species + the frozen boundary node — it
does NOT include other species' canopy contribution. So `stand_gradient(feedback=
"resident", species=k)` on a multi-species stand is **wrong** (the reconstructed light
omits the other species). `feedback="frozen"` is fine multi-species (it reads the joint
harvested env directly). Multi-species coupled = the cross-species resident Jacobian
`d(metric)/d(θ of species s')` — genuinely new (the mutant cross term was zero).

## Steps to round out FF16 (suggested order)

1. **Decompose the frozen-geometry gap (measure first).** Coupled `d(LAI)/dθ` vs a full
   `run_scm` FD differs ~20–33% for geometry traits (`lma` 33%, `a_l1` 20%), ~2% for
   `a_p1` (see the guide's "Honest scope" + `scripts/birthenv_test.R` logic). That gap
   bundles THREE frozen channels — quantify each before fixing:
   - **frozen adaptive step sizes** (a full re-run re-steps; the replay reuses θ₀ steps),
   - **frozen knot x-positions** (the spline grid),
   - **frozen birth-env establishment** (below).
   Build an FD that freezes the schedule (`refine_schedule=FALSE`) AND pins the step
   sizes / knots to isolate the birth-env piece. If step-size freezing dominates, the
   coupled gradient is doing the right thing and the "gap" is the grid response (a
   different, arguably-not-wanted quantity); if birth-env dominates, do #2.

2. **Active birth-env establishment.** In `assemble_metrics_coupled`, cohort initial
   conditions (`net0`/`pr_estab`/`mort0`/`g0`/`logd0`) use the FROZEN harvested birth
   env (`deep_net` over `F.eh[rn-1][5]`). Make it active: reconstruct the birth env from
   the already-alive cohorts (born `< rn`) via `build_interp` BEFORE adding cohort i
   (exact at θ₀, non-circular — the alive set at start of step rn == harvested stand at
   (rn-1,5), already proven by the static check), and compute the initial conditions
   with `deep_net_coupled`. Re-run `coupled_r0_check.R` (R0 must stay 1e-8) and the
   birth-env FD test.

3. **Active boundary node.** The boundary `new_node` is frozen (harvested ce_b/h_b) in
   the reconstruction — a tiny ground-level channel. Make it active: it is the pending
   seed, density = `birth_rate·pr_estab/g0` in the current reconstructed env (the
   `new_node_census` logic). Likely negligible; measure before bothering.

4. **Multi-species coupled (the cross-species Jacobian).** Generalise the per-RK-stage
   harvest to ALL species (`patch.h`: loop species, not `species[0]`) — additive, RcppR6
   nested one level deeper. Then `attach_coupled`/`build_interp` reconstruct the JOINT
   canopy from all species' re-evolved cohorts; re-evolve every species' cohorts
   together; differentiate w.r.t. one species' traits (or all). Validate R0 multi-species
   first (env recon over the joint stand), then the cross-species gradient vs full-SCM FD.

5. **Wire `feedback` default per metric + guard multi-species.** Decide census→resident,
   offspring→frozen as the documented default (scope decision #1). Until #4, `stand_
   gradient(feedback="resident")` should STOP with a clear error if the stand has >1
   species (currently it silently uses species[0]'s harvest — add the guard).

6. **CI + guide.** Promote a small coupled-resident FD check to a fast CI test (the
   short-schedule one already added). Update the guide's "Honest scope" once the gap is
   decomposed / birth-env is active.

## Pitfalls / notes carried over

- Reverse `adj` tape is odelia's single compiled global — SCOPE it (`{ }`); don't mix a
  live adj tape with a fwd pass. (The coupled gradient already scopes its tape.)
- The coupled reconstruction carries ~1e-8 VALUE noise (cohort-height-crossing sort
  discontinuities + re-evolution drift), so FD-validate at step ~1e-4 (noise-optimal),
  not 1e-6. AD is exact; FD is the noise-limited reference.
- Build: C++-only → `make compile`; new export/RcppR6 field → `make attributes`/`make
  RcppR6` then `make compile`; then `load_all(".", compile=FALSE)`. 137GB RAM — the
  coupled tape (~1e8 nodes) is not a memory concern here.
- Run the WHOLE plant suite before any PR (generic tests loop over helper-plant.R
  strategy lists).

## ROUND-OUT OUTCOMES (2026-06-29)

**Step 1 — decompose the gap. DONE.** `scripts/coupled_gap_decompose.R` (ODE-tolerance
sweep) + `scripts/coupled_birthenv_channel.R`. The coupled-AD vs full-SCM-FD gap splits
~half **grid response** (frozen adaptive step sizes + knot positions; vanishes as
ode_tol→0: lma 33%→14%, a_p1 1.8%→0.5%) and ~half **birth-env + boundary**. Grid
response is arguably *not* the wanted derivative (it tracks where the solver lays its
grid, not the biology). Conclusion: birth-env is material → did step 2; boundary is not
→ skipped step 3.

**Step 2 — active birth-env. DONE & SHIPPED (default on).** `assemble_metrics_coupled`
now reconstructs the birth canopy from the already-alive cohorts at each birth step
(gated by `Frozen::coupled_active_birthenv`, plumbed through `attach_coupled` +
`ff16_coupled_{metrics,gradient}_impl(..., active_birthenv=true)`). R0 stays exact at
θ0 (LAI rel 2.6e-9 — the alive set reproduces the harvested stand, non-circular). The
birth-env channel closes part of the gap (lma 33→29%, a_l1 20→16%); it slightly
overshoots the tiny a_p1 gap (which is ~all grid response). R1 AD-vs-FD still agrees
(the one >5e-3 entry, dLAI/dhmat, is a ~1e-8-magnitude derivative = noise floor).

**Step 3 — active boundary node. MEASURED, SCOPED OUT.** Boundary ce ≈ 0.4% of the
canopy sum and only in the lowest trapezium segment; its active-vs-frozen *derivative*
channel is far smaller — negligible vs the 10–30% grid/birth-env channels. Left frozen
(its circular density↔canopy dependency isn't worth <1e-3 of the gradient).

**Step 4 — multi-species coupled. SUBSTANTIALLY DONE.**
- *Harvest (patch.h):* added all-species per-RK-stage fields
  `stand_{height,competition}_stage_history_all` [step][stage][species][cohort] +
  `stand_newnode_{height,competition}_stage_history_all` [step][stage][species] +
  per-species `newnode_*_env_snapshot_all`, RcppR6-exposed (4D codegen works). Additive;
  species-0 single-species fields untouched. **Gate A** (`scripts/coupled_multispp_r0_static.R`):
  joint static recon = Σ_species per-species-trapezium/area reproduces the SCM env to
  3.3e-7 (single-species is 4e-16; the ~1e-6 is float summation order in the ~180-node
  joint canopy — the harvested per-node effect k_I·area_leaf is deterministic from
  heights, so the harvest is provably complete).
- *Engine (`src/ff16_emergent.cpp`):* `FrozenMS` + `assemble_metrics_coupled_ms` +
  `build_frozen_ms`; exports `ff16_coupled_metrics_ms_impl` (R0) +
  `ff16_coupled_gradient_ms_impl(..., target)` (R1, cross-species). Joint canopy =
  Σ_species per-species trapezium (NOT a merged sort — matches
  `Patch::compute_competition = Σ_s Species_s::compute_competition/area`). Metrics are
  TOTAL-stand (Σ_species). Differentiate one species' traits; ALL species lifted to the
  tape (canopy couples them), only the target's traits + its h0-IFT registered.
- *Validation:* **Gate B** (`scripts/coupled_multispp_r0.R`, fixed schedules): joint R0
  env drift 1e-6…1e-4, total LAI ≈ compute_competition(0) to ~1e-3. **Gate C**
  (`scripts/coupled_multispp_r1.R`): cross-species AD == coupled-recon FD (LAI clean to
  ~1e-2; mass/size metrics noise-limited — step sweep shows their FD bounces 10x around
  the stable AD), and the cross-species feedback is large/nonzero (dLAI/dlma total
  −0.47 vs target-frozen −1.03). CI: 3 new tests in `test-ff16-stand-gradient.R`
  (guard, MS-R0, MS-cross-species-R1) on fast fixed schedules.
- *WIRED INTO THE PUBLIC API (2026-06-30):* `stand_gradient(feedback="resident")` on a
  >1-species stand now routes to the cross-species engine and returns
  d(TOTAL-stand metric)/d(θ of `species`). A cheap double R0 pass gates it: if the joint
  re-evolution diverges (env drift > 1e-2 or non-finite) it raises a clear error asking
  for a fixed/uniform node schedule, instead of returning NaN. `ff16_harvest_ms` (R) does
  the all-species gather. Single-species + `feedback="frozen"` paths unchanged.
- *REMAINING (the one open item):* on the stiff **refined** production schedule (tiny
  ~1e-5 clustered early steps) the joint re-evolution goes stiff and **diverges** at a
  specific step (instrumented: env drift jumps 0→0.83 at one step when a single cohort's
  log-density runs away to ±inf via the backward-FD g' term; nalive≈179). This is a
  REPLAY artifact, NOT a model behaviour: a normal `run_scm` does not blow up (it steps
  adaptively against its own exact canopy — those 1e-5 steps ARE it controlling this
  stiffness); the replay reuses the frozen step sizes against a ~1e-3-drifted re-evolved
  canopy, which tips the cohort past the stiffness threshold. Absent from single-species
  refined (7e-6) and from fixed multi-species (≤1e-4). So the guard + fixed-schedule
  recommendation is the shipped resolution (Dan: "happy to take node schedules as fixed
  for the gradient calculations"). Possible future hardening if wanted: a log-density-rate
  cap mirroring the SCM's `check_initial_density_rates`, or sub-stepping the clustered
  early steps; not needed for fixed-schedule calibration use.

**Step 5 — feedback default + guard. DONE.** Public default stays `"frozen"` (safe,
backward-compatible no-feedback derivative); docs now state the scope decision
(census→`"resident"` is the total gradient, offspring always invasion). Multi-species
guard added in `stand_gradient` (clear error, points to `"frozen"`).

**Step 6 — CI + guide. DONE.** Guide "Honest scope" rewritten (grid-vs-birth-env
decomposition, birth-env now active, boundary negligible, single-species-only +
cross-species follow-up). New scripts listed above; 3 fast CI tests added.
