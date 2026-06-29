# Seed: rounding out the FF16 coupled resident gradient — #472 scope B (R1+)

Handoff for a fresh chat. Branch `spike-ff16-scm-emergent`. The coupled resident
replay (R0 + R1) is **built, validated, committed, and wired into the public API**.
Read `notes/resident-coupled-replay-seed.md` (the original design) and the guide
section "The resident total gradient: the coupled replay" in
`overstorey-staging/guides/autodiff-trait-gradients.qmd` first.

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
