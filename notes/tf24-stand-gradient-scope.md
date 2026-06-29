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

## 5. Feasibility verdict by metric × feedback

| target | TF24 | TF24F |
|---|---|---|
| `offspring_production`, frozen (invasion) | **DONE** (`tf24_offspring_production_gradient`) | inherits TF24's; minor (tracked-state seed) |
| census (LAI/biomass/size_moment), frozen | needs **curvature harvest** for g′; ~5e-7 floor | **clean** (analytic rates; carry collar state on tape) |
| census, **resident** (coupled, self-shading) | curvature harvest + per-stage leaf-opt in the re-evolution (stiff, many solves) | **clean + cheaper** (no optimiser; 1 leaf eval/step) |
| census, **multi-species** (cross-species) | as above; canopy ports free | as above; canopy ports free |
| `stand_state_jacobian` (escape hatch) | per-cohort, harvested profit | per-cohort, analytic |

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
