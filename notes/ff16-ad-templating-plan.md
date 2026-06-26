# FF16 scalar-templating plan (#472 scope B / #537)

Design record for making FF16 outputs differentiable w.r.t. traits by reverse-mode
AD. The strategy is to add a **third template axis** — the scalar type `S`
(default `double`) — to the `<T,E>` hierarchy. Throughout: **additive**, every
templated class/kernel keeps a `using X = X<...,double>` alias (or a default
template argument) so the existing package + R boundary stay bit-identical; AD
lives only in C++ (`S = xad::adj<double>::active_type`).

Proven before this branch (reused, not redone): `Internals` templating (K93
spike), the differentiable spline/`Interpolator` (odelia #32, merged), the
two-pass orchestration + frozen-quadrature replay, IFT at root-finds (A2, #539).

## Milestone A — single-plant trait gradient (no ODE, no demography)
Differentiate FF16 net mass production / height growth at a fixed size in a fixed
(double) light environment w.r.t. traits. Smallest end-to-end AD result; exercises
the real physiology kernel.
1. `Internals` → `basic_internals<S>` (+ `using Internals = basic_internals<double>`).
2. `FF16_Pars` → `basic_FF16_Pars<S>` (+ `FF16_Pars` alias; R/RcppR6 stays double).
3. FF16 physiology kernel templated on `S` (mass cascade, respiration/turnover,
   assimilation_leaf, net production). Double methods delegate to it.
4. **Quadrature — frozen-replay only.** Do NOT template the adaptive `QK::integrate`.
   A small templated fixed-weight accumulator replays recorded nodes/weights.
5. Validate `∂(net_mass_production)/∂trait` vs central FD (~1e-8).

## Milestone B — resident stand output via the two-pass (no SCM)
Differentiate an emergent stand quantity w.r.t. a trait, with the resident light
environment as a frozen-knot active-value spline.
6. `Individual<T,E,S>` — state/aux accessors in `S`; the ODE (de)serialisation
   stays the `double` boundary.
7. `compute_competition` accumulates `S` then feeds the spline knot values
   (frozen `x`, active `y`) — the resident self-shading coupling.
8. Wire the two-pass driver onto the real FF16 light env + assimilation.
9. Validate vs FD of the whole double two-pass.

## Milestone C — end-to-end through the SCM (the big one)
10. `Node<T,E,S>` (carries `Individual<…,S>`; demographic density state stays
    `double` — mixed serialisation at the ODE boundary).
11. `Species`/`Patch`/`SCM` templated on `S`; the **ODE state boundary** is the
    main coupling — `state_type = std::vector<double>` (odelia) stays `double`;
    replay a fixed step schedule (two-pass) so the taped pass is branch-free.
12. Differentiate an emergent SCM output (fitness / equilibrium density) w.r.t. a
    trait; validate. Mutant path: background env stays `double`, only the query is
    active → the active-query spline overload.

## Hardest couplings (ranked)
1. **ODE state boundary** (`Individual::ode_state`/`set_ode_state`): odelia's
   `iterator = vector<double>::iterator` stays double; the AD path steps a
   replay loop manually instead of going through the odelia solver.
2. **FF16_Pars / R split**: template internally, expose only `<double>` to RcppR6.
3. **assimilation integrand**: must return `S`; solved by the frozen-replay
   accumulator (step 4), not a QK rewrite.
4. **compute_competition → spline**: accumulate `S`, feed frozen-`x`/active-`y`
   spline (odelia #32).

## What does NOT need templating
- The adaptive `QK`/`QAG` controller, `AdaptiveInterpolator::refine`, the RKCK
  stepper — all pass-1 (double) schedule discovery, frozen and replayed.
- The RcppR6 / R boundary — binds the `<double>` specialisations only.
- odelia's ODE solver internals — external; the `double` state boundary is the
  contract. (The AD replay sidesteps it with a fixed-step manual integrator.)

---

## Status on branch `spike-ff16-hierarchy` (PR #541)

Milestones A and B-core landed earlier (`Internals`/`FF16_Pars`/kernel templating,
deep-crown frozen replay, A1 exact growth gradient; PRs #539/#540). Milestone C
was built up additively as a series of bit-identical increments — the full test
suite stays green (and the double path bit-identical) at every step:

| Incr | Piece | Validation |
|---|---|---|
| 12–14 | `Node` / `Species` / `Patch` templated on `S=double` | suite bit-identical |
| 15 | FF16 `update_dependent_aux<S>` + `area_leaf<S>` → a live `Node<…,ad>` constructs | d(area_leaf)/d(height) ~2e-10 |
| 16 | `FF16Rates<S>` + `ff16_compute_rates_crown_top<S>` (full demographic rate fill) | 5 rates bit-exact vs live; d(fec)/d(a_p1) ~1e-12 |
| 17 | `FF16State<S>` + `ff16_grow_demography<S>` (full 5-state trajectory) | d(lifetime fecundity@T)/d(lma) ~1e-8 |
| 18 | `ff16_replay_cohort<S>` (frozen-schedule two-pass replay primitive) | multi-trait d(stand LAI)/d{lma,a_p1} ~2–5e-10 |
| 19 | `ff16_replay_cohort_active_light<S>` (within-cohort crown-light feedback) | d(height@T)/d(lma) w/ feedback ~5e-11 |
| 20 | `ff16_resident_light_at<S>` (full self-shading: profile responds to trait) | d(focal net)/d(a_l1) through active-knot spline ~2e-11 |

Every kernel lives in `inst/include/plant/models/ff16_production_kernel.h`. A
runnable demonstration of the gradients (each checked against finite differences)
is `scripts/ad_gradient_examples.R`.

**Remaining (production integration, not feasibility):** drive pass-1 from the
**live SCM** — harvest the real node-introduction schedule and the per-RK-sub-step
environment cache (`Patch::environment_history`, the `save_RK45_cache` machinery),
build the frozen per-cohort crown-light schedule, and call the replay kernels
under AD. The public `run_scm(collect=TRUE)` output samples at schedule events, not
integration sub-steps, so this is a C++-side entry rather than R plumbing.
