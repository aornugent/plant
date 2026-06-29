# Scope: TF24 / TF24F stand-metric trait gradients — feasibility (#472 scope B)

Scoping note (no code). Branch `spike-ff16-scm-emergent`. Question from Dan: now that
the FF16 stand-metric gradients are done (frozen + resident census, single- and
multi-species), how well can the same work for **TF24**, given the leaf hydraulics once
suggested reverse-mode AD "might not work" — and would restricting to **TF24F** make it
tractable? Read the guide's TF24 section
(`overstorey-staging/guides/autodiff-trait-gradients.qmd`) for the leaf-gradient story.

## TL;DR

1. **Reverse mode is NOT lost on TF24 — it already ships.** `tf24_offspring_production_
   gradient()` is a committed reverse-mode emergent gradient through the live TF24 SCM.
   The leaf optimiser/root-finds are *never taped*; they run in a double discovery pass
   and only their converged sensitivities (envelope theorem + IFT) are injected into a
   linearised, fully-tapeable replay. So the worry is already answered for the offspring
   metric; the open question is only the **census** metrics.
2. **The coupled-canopy machinery ports from FF16 for free.** In TF24/TF24F `area_leaf`
   is pure allometry (`(height/a_l1)^(1/a_l2)`), and competition is
   `density·k_I·area_leaf·Q` — the leaf optimiser affects only the demographic *rates*,
   **not the light-field geometry**. So every piece of the FF16 coupled resident /
   multi-species canopy reconstruction (`coupled_comp_at`, `build_interp`, the
   per-species trapezium, the joint multi-species canopy, the active-knot spline) reuses
   unchanged. That is the bulk of `stand_gradient`, and it is identical work.
3. **The one genuinely new ingredient is the census number density's growth-rate
   gradient** `g' = d(height_dt)/d(height)`. It is the reason the guide lists TF24
   census as a follow-up. It is solvable (see §3), and **TF24F makes it clean** because
   TF24F removes the optimiser and already has exact analytic leaf derivatives.
4. **Recommendation: do TF24F, not TF24, for the census/resident stand gradients** — it
   is the cleaner, higher-fidelity, less-hacky path, and it aligns with the acclimation
   research direction. TF24 (instantaneous-optimum) census is possible but needs a
   second-order "curvature harvest" and carries TF24's ~5e-7 optimiser-tolerance floor.

## 1. The leaf solvers, and why reverse mode works anyway

Per TF24 net-production evaluation (from the leaf model):
- **Golden-section optimiser** over the root-collar potential `opt_root_psi`
  (`Leaf::find_root_collar_psi`, `util::golden_section_max`, tol `GSS_tol_abs`).
- nested **soil→collar root-find** `find_root_psi` (TOMS748, water-balance residual).
- nested **stomatal `psi_stem→ci` root-find** `psi_stem_to_ci` (TOMS748).
- (`E_from_Soil_to_Root_Collar` is closed-form, not iterative.)

None of these is taped. The shipped principle (the guide's "never differentiate a
solver — differentiate its converged point"): a double pass discovers the operating
point; then `Leaf::dprofit_d{vcmax25,g1,beta2,b,c,K_s,...}` (forward-mode AD over the
analytic algebra + **IFT** at `psi_stem→ci` + analytic spline derivatives for transport,
+ **envelope theorem** at the optimised collar) give the converged trait sensitivities,
which are injected into the tape. `tf24_emergent.cpp` harvests, per RK stage per cohort,
`{profit, dprofit_dh, k_max, E_up, and the 10 dprofit_dθ}` and replays a **linearised
profit** `profit(h,θ) = profit₀ + dprofit_dh·(h−h₀) + Σ_k inj_k·(θ_k−θ_k0)` — leaf-opt-
free and fully reverse-tapeable. Fidelity ~5e-7 (the optimiser tolerance), vs FF16's
~1e-14 closed-form. So reverse mode is demonstrably fine; it is the *operating-point
harvest* that makes it fine.

## 2. What ports from FF16 unchanged (the canopy)

`TF24_Strategy::competition_effect()` returns `area_leaf(height)` (allometry), and the
SCM canopy is `Σ_species Σ_cohorts density·k_I·area_leaf·Q / area` — structurally
identical to FF16. Therefore, with a TF24 prod-pars harvest in place:
- the generic `(w,f)` reduction engine, the census height-trapezium + pending-seed tail,
- the coupled resident replay's per-stage canopy reconstruction (active heights +
  log-densities → light spline at frozen knots),
- the multi-species joint canopy (sum of per-species trapeziums) and the cross-species
  Jacobian, and even the refined-schedule stiffness guard,

all reuse the FF16 code paths with only the per-cohort *rate* call swapped (FF16's
closed-form `deep_net` → TF24's harvested-profit leaf evaluation). The geometry, the
trapezium arithmetic, the knot handling, the tape structure are unchanged.

## 3. The one new obstacle: the census density's g′

Census metrics (LAI / biomass / size_moment) are size-distribution integrals over the
replayed **number density** `density_i = exp(log_density_i)`, where
`d(log_density)/dt = −g' − mortality_dt` and `g' = ∂(height_dt)/∂height`. The SCM forms
`g'` by a backward finite difference of `height_dt` over a 1e-6 height step — which
**re-solves the leaf optimiser at `h−1e-6`**. `height_dt ∝ net = profit·area_leaf`, so
`g'` carries `∂(profit)/∂height` *and its trait sensitivity*.

The current TF24 harvest linearises profit **linear in height** (`dprofit_dh` constant),
so its second derivative is zero → it reproduces `g'` only to first order, missing the
curvature `∂²profit/∂h²`. That is exactly why the guide says the census density "needs a
leaf-optimisation cross-sensitivity the linearised harvest does not differentiate
faithfully." Two fixes, both additive:

- **(TF24) curvature harvest.** Harvest the leaf operating point + its trait
  sensitivities at *both* `h` and `h−GEPS` (the two points the SCM's `g'` uses), per
  stage. Then the linearised `g'` and `d(g')/dθ` reproduce the SCM's own scheme exactly.
  Cost: one extra double-precision leaf solve per stage (the forward TF24 emergent path
  already does a central FD of the leaf opt for within-trajectory `d(net)/d(height)`, so
  the evaluate-at-perturbed-height machinery exists). Fidelity floor stays ~5e-7.
- **(TF24F) no harvest hack needed — see §4.**

## 4. Why TF24F is the right target

TF24F (acclimation, #525/#527) replaces TF24's per-step collar optimiser with a **6th
ODE state** `opt_root_psi_state` evolving by gradient ascent
`dψ/dt = k_acclim·dprofit_dψ`, evaluating the leaf at the *tracked* collar
(`evaluate_root_collar_psi`, no golden-section). Consequences for AD:

- **The golden-section optimiser is gone** from the trajectory (the least AD-friendly
  piece; only re-run once at birth to seed the state).
- **The soil→collar root-find `find_root_psi` is gone** from the per-eval path (closed-
  form `E_from_Soil_to_Root_Collar`). Only the `psi_stem→ci` root-find remains — and it
  is **IFT-able** (already handled inside `dprofit_droot_collar_psi`).
- **Exact analytic leaf derivatives already exist** for TF24F's own acclimation tracking:
  `Leaf::dprofit_droot_collar_psi` (fwd-AD + IFT + spline deriv) and the full
  `dprofit_dθ` family. Forward-mode acclimation gradients of the TF24F system already
  exist (#527).

So TF24F's trajectory is a clean 6-state ODE whose rates are **analytically
differentiable** (the one root-find handled by IFT). The census `g'` problem dissolves:
because the collar is a *replayed state* with an analytic rate and the leaf-at-given-
collar derivatives are analytic, the within-trajectory height-sensitivity propagates on
the tape directly (carry `opt_root_psi_state` as the 6th replayed state) — **no
constant-slope linearisation, no curvature harvest**. TF24F stand gradients can be as
clean as FF16's, modulo injecting the single `psi_stem→ci` IFT sensitivity and the
~optimiser-free fidelity (the leaf eval is now analytic, so the floor should be far
better than TF24's 5e-7 — closer to the spline/IFT precision).

Caveats to flag, not blockers:
- TF24F is **not bit-identical to TF24** — it tracks the optimum lazily (lags by a
  gain-dependent amount). Restricting gradients to TF24F means differentiating the
  *acclimation* model, not the instantaneous-optimum model. Decision for Dan: if TF24F
  is the model to calibrate/analyse (acclimation is a stated research priority), this is
  a feature, not a compromise.
- The `evaluate_root_collar_psi` **clamp** to the feasible interval `[bound_a,bound_b]`
  is a `max/min` — a kink. At the clamp boundary the derivative is one-sided; harvest the
  active branch (frozen boolean, the established pattern for the net-production clamp /
  mortality guard) so the taped expression is straight-line.
- TF24F adds a state, so the coupled replay carries 6 states/cohort not 5; trivial.

## 4b. The single-plant `grow_individual_to_size` surface (FF16 done; TF24 analogue)

FF16 now has `grow_individual_to_size_gradient()` — a trait gradient of one plant grown
in a **fixed** environment to a target size (the guide's "The same pattern on a single
plant"). It is the lightest gradient: no resident feedback, so it is just the two-pass
frozen-schedule replay + the **IFT stopping-time step** on `size(t*,θ)=target`
(`dt*/dθ = −(∂size/∂θ)/size_dt`), returning `d(t*)/dθ` and the total `d(state at t*)/dθ`.
The TF24 analogue would be wanted for the same use case (e.g.
`optimise_individual_rate_at_size_by_trait`, the leaf/whole-plant growth optimisations).

How it ports to TF24/TF24F — it is **strictly easier than the TF24 census** above,
because there is *no canopy and no density*:
- **No `g′` / number-density problem at all.** The single-plant grow tracks only the
  demographic state along one trajectory; there is no `log_density` rate, so the
  census-density curvature obstacle of §3 simply does not arise.
- **The replay machinery is the existing TF24 emergent path.** TF24's
  `tf24_offspring_production_gradient` already replays one cohort's demographic ODE over
  a frozen Cash-Karp schedule with the leaf opt harvested per stage and a linearised,
  tapeable profit injected. The grow-to-size gradient is that same per-cohort replay,
  driven by a *fixed env* (no per-stage resident light to harvest — the env is given) and
  stopped by the IFT step instead of integrated to the census. So it reuses the harvested
  -profit replay wholesale; only the FF16 `deep_net` rate call is swapped for the TF24
  leaf evaluation, exactly as for the offspring gradient.
- **The stopping-time IFT needs `size_dt(t*)` only** — `height_dt` at the operating
  point, which the harvest already produces (it is the `height_dt` whose within-trajectory
  `d/dh` the emergent path already central-differences). No new sensitivity.
- **TF24 fidelity floor stays ~5e-7** (the optimiser tolerance), as for offspring; **TF24F
  is clean** (analytic leaf rates, collar carried as the 6th replayed state — the same
  reason it is the right census target).

So the per-plant grow-to-size gradient for TF24 is low-risk follow-on work that should be
done *with* (or before) the TF24F census, sharing the harvested-profit replay; it does
**not** need the curvature harvest §3 requires. Build order: a TF24 grow-to-size R0
(double replay reproduces `grow_individual_to_size` to the optimiser floor) then R1 (AD
vs frozen-schedule FD), mirroring `test-ff16-grow-individual-gradient.R` and
`scripts/ad_grow_individual_gradient.R`.

## 5. Feasibility verdict by metric × feedback

| target | TF24 | TF24F |
|---|---|---|
| `offspring_production`, frozen (invasion) | **DONE** (`tf24_offspring_production_gradient`) | inherits TF24's; minor (tracked-state seed) |
| census (LAI/biomass/size_moment), frozen | needs **curvature harvest** for g′; ~5e-7 floor | **clean** (analytic rates; carry collar state on tape) |
| census, **resident** (coupled, self-shading) | curvature harvest + per-stage leaf-opt in the re-evolution (stiff, many solves) | **clean + cheaper** (no optimiser; 1 leaf eval/step) |
| census, **multi-species** (cross-species) | as above; canopy ports free | as above; canopy ports free |
| `stand_state_jacobian` (escape hatch) | per-cohort, harvested profit | per-cohort, analytic |
| `grow_individual_to_size` (single plant, fixed env) | reuses harvested-profit replay + IFT stop; no g′/density; ~5e-7 floor | **clean** (analytic rates + IFT stop); FF16 done |

The canopy/geometry half is free for both. The only real cost is the leaf-opt sensitivity
in the rates+density — manageable for TF24 (curvature harvest), **clean for TF24F**.

## 6. Suggested build order (if we proceed — recommend TF24F)

1. **TF24F frozen census R0 gate.** Double replay (6 states incl. collar) reproducing the
   SCM's census density/heights to the optimiser-free floor; confirm LAI ==
   `compute_competition(0)`. This proves the collar-state replay is faithful before any AD.
2. **TF24F frozen census R1.** Reverse sweep for d(census metric)/dθ (all 27 traits),
   carrying `opt_root_psi_state` on the tape, IFT-injecting the `psi_stem→ci`
   sensitivity; validate vs a two-pass FD over the same replay.
3. **TF24F resident census (coupled).** Swap the FF16 per-cohort rate call for the TF24F
   leaf-at-tracked-collar evaluation inside the existing coupled engine; the canopy code
   is unchanged. Validate R0 (joint env) then R1 (AD vs recon-FD), on a **fixed** node
   schedule (the FF16 refined-schedule stiffness guard applies — TF24F is stiffer still,
   so fixed schedules are the expected mode anyway).
4. **TF24F multi-species** falls out of the FF16 multi-species path (canopy is generic).
5. **(Optional) TF24 instantaneous-optimum census** via the curvature harvest, only if
   the non-acclimation model's census gradient is specifically wanted.

## 7. Risks / unknowns to probe first (cheap experiments)

- **Confirm the collar-state replay is faithful** (step 1) before investing in AD — TF24F
  tracks lazily, so the replayed `opt_root_psi_state` must match the SCM's stored state;
  if the gain `k_acclim` makes the trajectory schedule-sensitive, pin it.
- **Measure the leaf-eval fidelity floor** for TF24F (analytic path) vs TF24 (5e-7): a
  quick `dprofit_droot_collar_psi` AD-vs-FD check at a census operating point tells us the
  achievable census-gradient precision.
- **Clamp frequency.** How often does `evaluate_root_collar_psi` hit the `[bound_a,bound_b]`
  clamp along a real trajectory? Frequent clamping means many one-sided-derivative stages
  to harvest as frozen branches.
- **Cost.** TF24/TF24F coupled re-evolution does a leaf eval per cohort per stage; on a
  realistic stand that is many solves in the (double) discovery pass. TF24F (no optimiser)
  is ~3-4× cheaper per eval than TF24. Benchmark before committing to the resident path.

## 8. Progress — step 1 (TF24F frozen census R0 gate) DONE

Built on `spike-ff16-scm-emergent`. The R0 gate from §6 step 1 is implemented and green:

- **`src/tf24f_emergent.cpp` — `tf24f_census_recon_impl`** (double only, no XAD): a
  per-cohort replay of the **7-state** system {5 demog, tracked collar
  `opt_root_psi_state`, `log_density`} over the SCM's frozen Cash-Karp schedule. It
  drives the **real** TF24F leaf at the tracked collar (set `tracked_root_psi_`, call the
  inherited `net_mass_production_dt` → `solve_leaf` evaluates, not optimises — crown-centre
  = one analytic leaf solve/stage). Seeds each cohort as `Node::compute_initial_conditions`
  does (collar₀ from the one birth-time optimiser run via `initializing_`; establishment
  pr_estab at the tracked collar₀; logd₀ = log(birth_rate·pr_estab/g₀)). `g'` replicates
  the SCM's own scheme exactly: backward FD of `height_dt`, abs step
  `node_gradient_eps = 1e-6`, collar held, height → h−eps (defaults: dir −1, exact_ad off).
  Census reduction = the height-trapezium of `density·kI·area_leaf` (+ pending-seed tail),
  identical geometry to FF16 / `Species::compute_competition`.
- **`R/tf24f_emergent_gradient.R`** — `tf24f_harvest` (TF24f-guarded; returns the frozen
  env+schedule, birth steps, params, `k_acclim`/`use_ad_gradient`) and the internal
  `tf24f_census_recon` validation entry.
- **`tests/testthat/test-tf24f-census-gradient.R`** — the gate (6 checks, plain R/CI).

**Result (the §7 risk-1 check: is the collar-state replay faithful? — YES).** On the gate
stand (lifetime 4, 9 nodes, crown-centre, `k_acclim=1`), the recon reproduces the SCM's
stored stand essentially bit-exactly: heights/collar to ~1e-6 (8 of 9 cohorts to ~1e-14),
`log_density` to ~1e-6, **LAI vs `compute_competition(0)` to 6e-7**. The collar trajectory
is tracked faithfully, confirming the 6-state replay before any AD.

**[CORRECTED — see §9] Fidelity-floor measurements (§7 risk-2).** Scaling the patch
(tallest cohort in parens):

| lifetime | max │dH/H│ | max │d log_density│ | LAI reldiff |
|---|---|---|---|
| 3 (5.8 m) | 3.3e-7 | 1.9e-6 | 2.0e-7 |
| 4 (7.2 m) | 8.2e-7 | 3.9e-6 | 6.3e-7 |
| 5 (8.4 m) | 1.5e-6 | **6.1e-3** | **−4.0e-3** |

At lifetime 5 the oldest cohort's `log_density` jumps to 6e-3 (LAI 0.4%) while heights
stay ~1e-6. At **lifetime 8 the run aborts** ("Extrapolation disabled…"): the tracked
collar drifts outside the leaf transport spline domain. So R0/R1 work should stay at
lifetime ≤ 4–5; longer horizons need the collar drift bounded first (clamp/spline-domain
guard), independent of the gradient work.

## 9. The `g'` finite difference, and what using AD for it actually showed

Dan asked: shouldn't the census `g' = d(height_dt)/d(height)` use AD, not the SCM's
1e-6 backward FD — and that FD is in **every** `run_scm`, not just the gradient. Both are
true; the FD lives in `Node::growth_rate_gradient` (all strategies) and is the default
because `control.node_gradient_exact_ad` defaults off. The exact-AD hook
(`Strategy::growth_rate_gradient_height_ad`, #537 A1) is implemented only by **FF16**;
**TF24/TF24f/K93 inherit the base `NA`** and silently fall back to FD even with the flag
on. (Two separate things: the R1 *trait* gradient already uses AD — it differentiates
*through* whatever `g'` expression the SCM used, no θ-FD. The FD here is purely the
*height* derivative defining `log_density`.)

**Built this session.** `TF24f_Strategy::growth_rate_gradient_height_ad` — analytic
`d(height_dt)/d(height)` at the fixed tracked collar: forward-AD of the demographic growth
kernel with `profit` carrying its height derivative summed over the leaf channels height
moves. A channel decomposition showed `d(profit)/d(height)` is dominated by the crown-
centre **light** channel (~110% of the total), with **kmax** (`∝1/h`) a ~−10% correction
and the feared **E_up / root-network** term negligible (exactly 0 above the 1.5 m rooting
clamp, ≤1.6% for seedlings — so no root-network derivative is needed). New leaf sensitivity
`Leaf::dprofit_dPPFD` (the light channel; a `which==4` case of `dprofit_dphoto` through
`electron_transport`). **Exposed `node_gradient_exact_ad` in the R `Control` interface**
(it was C++-only since #537 A1) — so `run_scm` can now switch on the exact gradient for
FF16 **and** TF24f for the first time.

**What it showed, and how the two follow-ups (Dan: "do 1 and 2") resolved.**
1. **Correct in smooth / fixed light** — matches a fine FD of the live growth rate to
   ~1e-6 (`test-tf24f-growth-gradient-ad.R`).
2. **(2) Real-env `g'` — RESOLVED (it was correct all along).** An *in-`run_scm`*
   (native-env) check found the analytic `dprofit_dh` matches a native profit-FD to
   ~1e-4 and the AD `g'` matches a native height_dt-FD to ~1e-5 **even for the tallest
   cohort** (h=8.4: AD 0.61761 vs FD 0.61759). The earlier "~18% off" was an **artifact
   of validating through an `Rcpp::as<>`-round-tripped environment** (my removed
   `decomp`/`gpcheck` probes and the R-side `Individual` check all passed the env through
   `as<>`, which is NOT faithful for the crown-sampled light above cohort heights — it
   gave a different growth rate than the native env, 0.617 vs 0.442). Lesson: validate
   env-dependent gradients **inside the run**, never through an R→C++ env round-trip.
3. **(1) Lifetime-≥8 abort — FIXED (a real #527 bug).** Plain TF24 runs to lifetime 12;
   only TF24f aborted ("Extrapolation disabled"), and only with the **AD** acclimation
   gradient (`use_ad_gradient = TRUE`; the FD gradient ran fine to lifetime 12). Cause:
   under drought, `Leaf::prepare_collar_solve` returns false and sets the **shutdown**
   operating point (collar = −`psi_crit`); the FD branch of `TF24f::solve_leaf` gates on
   that and returns 0, but the **AD branch did not** — it called
   `Leaf::dprofit_droot_collar_psi` at the shutdown collar, where the non-extrapolating
   transport spline's *derivative* is out of domain → hard abort. Fix:
   `solve_leaf`'s AD branch now gates on `prepare_collar_solve` exactly like the FD
   branch (shutdown → `dprofit_dpsi_ = 0`). TF24f now runs with the default exact
   gradient to **lifetime 16** (`test-tf24f-growth-gradient-ad.R`, drought-robustness
   case). NB this was **not** collar "runaway": TF24's optimiser also reaches collar
   magnitudes ~2.6 under stress — the deep collar is physical; only the AD gradient's
   spline-edge handling was missing.

**Remaining census-floor note (separate from both).** The lifetime-5 census *recon*
floor (LAI ~0.4%) is **not** the `g'` method and **not** the collar trajectory — both
are now sound. It is the `tf24f_census_recon` engine reading the resident light through
the same unfaithful `Rcpp::as<>` round-trip (point 2): faithful enough for the value-
based height trajectory (heights match ~1e-6) but not for the slope-sensitive `g'`/
`log_density` of the tallest cohort at long horizons. The R0 gate stays at lifetime ≤ 4
(bit-faithful). For R1, prefer driving the replay from native envs (or accept the gate's
≤4 horizon); the gradient itself is AD-exact for the density the SCM computes.

## 10. Progress — the frozen-census R1 AD tape (build-order step 2) DONE

The refine step the seed (`notes/tf24f-census-tape-seed.md`) calls for is implemented and
green: the per-trait FD census gradient is replaced by ONE reverse-mode AD sweep per
metric over the 7-state replay.

- **`src/tf24f_emergent.cpp` — `tf24f_census_gradient_ad_impl`** + R entry
  `tf24f_census_gradient_ad` (`R/tf24f_emergent_gradient.R`). A double DISCOVERY pass
  harvests, per cohort per RK stage, the trait-independent leaf operating point at the
  tracked collar — at BOTH the forward height `h` and the backward `h−GEPS` the SCM's `g'`
  reads — then a second pass replays a leaf-opt-free, fully tapeable 7-state expression
  {5 demog, collar, log_density} and takes one adjoint sweep per metric. The census
  trapezium couples cohorts, so all cohorts share one tape; one reverse sweep per metric.
- **The collar curvature harvest (the one hard ingredient).** The tracked collar is a
  taped STATE with rate `k_acclim·dprofit_dpsi`, linearised with the FD-harvested second
  derivatives `d²profit/dpsi²`, `d²profit/dpsi dh`, `d²profit/dpsi dθ_k`; profit itself
  is linearised in (h, collar, θ) — the `dprofit_dpsi0·(collar−collar0)` term the offspring
  tape omits is added (collar is taped). The collar's birth seed (= optimum) is injected
  by IFT `d(collar0)/dθ = −(d²p/dpsi dθ)/(d²p/dpsi²)`. The per-trait sensitivities (`inj`,
  `d²p/dpsi dθ`) are FD over independently-built ±θ-perturbed strategies, matching the FD
  gate (which perturbs the same `pars` vector). Prototype trait subset, not all 27.
- **`g'` reproduces the SCM's backward-FD scheme** by differencing the two tapeable
  `height_dt` expressions at `h` and `h−GEPS` (as `ff16_emergent.cpp`'s census tape does) —
  no `d²profit/dh²` harvest needed. (The `node_gradient_exact_ad` `g'` path is a follow-up;
  the gate stand uses FD `g'`, which this matches.)

**Validation (`test-tf24f-census-gradient.R`, R1 tape test, lifetime-4 gate stand).** AD
values match the recon/SCM bit-for-bit (linearisation exact at θ0). AD Jacobian matches
`tf24f_census_gradient_fd` to **≤0.8%** (within the recon noise floor) over
{LAI, biomass, size_moment} × {vcmax_25, lma, a_l1, K_s}. The AD sits at the convergent
limit of the well-behaved FD steps (1e-4/1e-5) and removes the step sensitivity — FD at
1e-6 degrades (LAI/lma 13.20 vs AD 13.55) from the leaf-solve noise floor amplified by the
small step. Full plant suite green (FAIL 0, PASS 2582).

**Step 3 (first-class API) DONE.** `stand_gradient(scm, …)` now dispatches the TF24f
strategy: the census metrics under `feedback = "frozen"` route to `tf24f_census_gradient_ad`;
`feedback = "resident"` and `offspring_production` stop with a clear "follow-up" message.
Test in `test-tf24f-census-gradient.R`; commit `fcbfc60d`.

**Remaining (the seed's later steps).** Step 4: individual grow-to-size AD tape (same collar
curvature harvest, no canopy/density; validate vs `tf24f_grow_individual_to_size_gradient_fd`).
Step 5: resident/coupled AD (swap the FF16 coupled per-cohort rate for the TF24f
harvested-leaf-at-tracked-collar eval; fixed schedule; validate vs
`tf24f_resident_census_gradient_fd`). The C++-native-env migration
(memory `move-gradient-machinery-to-cpp`) still gates lifetime >4 fidelity.

## Bottom line for Dan

Reverse-mode AD is **not** a loss for TF24 — it already ships for the offspring metric,
and the principle (harvest the converged operating point, inject IFT/envelope
sensitivities, never tape the solver) generalises. The stand **census** metrics are the
only new work, and the hard part (`g'`) is **clean on TF24F** because TF24F removes the
optimiser and already carries exact analytic leaf derivatives for its acclimation
tracking. Recommendation: target **TF24F** for the census/resident/multi-species stand
gradients; it is more faithful, cheaper, and aligned with the acclimation direction.
Restrict TF24 (instantaneous optimum) to the already-shipped offspring gradient unless
its census gradient is specifically required (then a one-extra-leaf-solve curvature
harvest delivers it at the ~5e-7 floor).
