# Resident-gradient build — scope (#472 scope B, the next unit)

Design agreed 2026-06-29 (Dan + AD spike), after the emergent-gradient speed work.
This scopes the **resident total gradient** of emergent stand metrics, the piece
deliberately deferred when the committed `stand_gradient` shipped only the
mutant/frozen-resident gradient. Builds on `scm-gradient-architecture.md` and the
active-knot prototypes C-27/C-28.

## The two gradients (recap)

For an emergent metric `M` and a trait `θ`:

```
dM/dθ  =  ∂M/∂θ |_env        +        (δM/δenv) · (denv/dθ)
          └ mutant / invasion ┘        └─ resident feedback (NEW) ─┘
```

- **Mutant (frozen-resident):** env held fixed; cohorts independent on the tape.
  SHIPPED (`offspring_production_gradient`, `stand_gradient`, escape hatch). Correct
  & complete for `offspring_production` (= the invasion/selection gradient).
- **Resident total:** the canopy + soil water co-vary with θ. The right quantity for
  **LAI / biomass / size-moment** as resident ecosystem outcomes (and for resident
  community calibration). NOT built. C-27 showed the feedback term can *dominate and
  flip the sign* for allometric traits — so the mutant reading of LAI/biomass is not
  just incomplete, it can be qualitatively wrong.

## The design: freeze geometry, recompute knot values (Dan's, confirmed)

Avoid re-running anything adaptive. From the resident (pass-1) run, **freeze**:

- the ODE step times + the cohort-introduction schedule,
- the resource-spline **knot x-positions** (the heights at which light is sampled;
  the soil depth/time grid for water).

Then in pass 2 (the AD replay), recompute the **knot y-values** — light levels, water
levels — as *active* functions of the (trait-dependent) stand state. This is exactly
odelia #32's differentiable spline (frozen knot positions, active knot values).

What this buys: no adaptive stepping, no schedule refinement, no spline-knot
re-placement — the geometry is fixed, only the values flow the trait.

## What changes vs the mutant replay

| | mutant (shipped) | resident (this build) |
|---|---|---|
| environment | frozen `double` constant | **active** (knot y-values respond to θ) |
| cohort coupling | independent | **coupled** through the shared per-stage env |
| replay shape | each cohort alone over frozen env | **all cohorts stepped together** through the frozen schedule (a manual fixed-schedule Patch stepper) |
| tape | per-cohort (tiny) | one **coupled** tape (the whole stand) |
| per-cohort speed trick | applies | does NOT apply |

The coupled tape is the thing that *looked* expensive — but today's measurement shows
the monolithic tape is ~1.3s, so the coupled replay is affordable. The replay becomes
a templated mini-SCM: introduce cohorts at frozen birth steps, step all alive cohorts
together; at each RK stage reconstruct the active env from the current (active) stand,
each cohort reads its crown light/water from it. (This is the manual `Patch<...,ad>`
fixed-schedule stepper the earlier notes anticipated, now actually needed.)

## Two coupling channels — different structure

### Light (FF16 + TF24) — prototyped
- Canopy = `ResourceSpline` over height; `light(z) = exp(-competition(z))`,
  `competition(z) = trapezium_i ce_i · Q(z/h_i)`.
- C-27/C-28 reconstruct it active: `ce_i = C_i · area_leaf_i`, `C_i` frozen,
  `area_leaf_i` active in the trait; fill an odelia `basic_interpolator` at the frozen
  knot heights with these active values. Validated: C-27 (static census, sign flip vs
  mutant), C-28 (time-integrated over a focal cohort, ~1e-4 faithful).
- Harvest already present: `Patch::stand_height_history`,
  `stand_competition_history` (per ODE step, **species[0] only**), R-exposed.

### Water (TF24 only) — NOT built
- Soil moisture is a multi-depth state vector (`soil_number_of_depths`, default 5)
  evolved as ODE aux variables, driven by stand transpiration + infiltration; `psi_soil`
  derived from it. NOT a spatial spline — a *temporal* soil-water balance.
- Needs: a per-stage harvest of the soil-water state + the stand-total transpiration
  (active in θ via leaf area × leaf-level transpiration), and the water-balance ODE
  differentiated under frozen geometry (freeze depths/times, recompute moisture active).
- This is the genuinely new, hardest channel; nothing reusable from C-27/28.

## Cross-cutting requirements

1. **Per-RK-stage stand harvest** (not just per-step). C-28's per-step reconstruction
   was ~1e-4 faithful; matching the SCM bit-closely (as the mutant replay does, 3e-14)
   needs the stand state (every cohort's height + leaf area, + water state) at each of
   the 6 Cash-Karp stages, like `environment_history`. New harvest fields.
2. **Multi-species:** the active env is the joint canopy of ALL species; the harvest is
   currently `species[0]`. Generalise to all species (a cross-species Jacobian becomes
   nonzero here — the distinct quantity from the mutant cross term that was zero).
3. **Baseline self-consistency:** at θ₀ the recomputed active env must reproduce the
   resident's frozen env (else the value drifts). Decide the tolerance / whether to
   anchor the reconstruction to the frozen env at θ₀.

## Build order (FIRST CUT = R0–R2, decided 2026-06-29)

- **R0 — coupled stepper + baseline check (FF16, light).** Build the all-cohorts-
  together fixed-schedule `Patch<...,ad>` replay over the frozen schedule, reconstructing
  the active light spline per RK stage from the **per-RK-stage** harvested stand state.
  First goal: at θ₀ the coupled double replay reproduces the SCM's emergent metrics + the
  per-stage light bit-closely (the consistency gate, target the mutant replay's ~3e-14,
  not C-28's ~1e-4). No gradient yet.
- **R1 — resident gradient, FF16 census.** One coupled reverse sweep per metric →
  resident dLAI/dθ, dbiomass/dθ, d(size-moment)/dθ. Validate vs an active-resident
  two-pass FD (perturb θ, re-run the resident-reconstruction, re-reduce). This is where
  the C-27 sign flip becomes a committed number.
- **R2 — TF24 light (END OF FIRST CUT).** Same coupled stepper with the harvested
  leaf-opt operating points; the active light enters net via the leaf's assimilation, so
  harvest the extra envelope sensitivity d(profit)/d(light) at each operating point.
  Validate + iterate before taking on water.
- **R3 — TF24 soil water (deferred past first cut).** The new water channel: per-stage
  water-state harvest + active soil-water balance under frozen geometry.
- **R4 — tests + guide (with R2/R3).** Phase-E calibration-setup demo in the guide.

## Decisions (2026-06-29, Dan)

1. **API surface — DECIDED:** a `feedback` argument on `stand_gradient`
   (`feedback = c("frozen", "resident")`, frozen = the shipped mutant path). One unified
   entry point. `offspring_production` stays `frozen` by default (invasion fitness); the
   census metrics are the resident-flavour consumers. (Exact default per-metric: settle
   at R1 — likely census→resident, offspring→frozen.)
2. **First cut — DECIDED:** R0–R2, i.e. **FF16 + TF24 light**. Stop and validate before
   the TF24 water channel (R3).
3. **Faithfulness — DECIDED:** **per-RK-stage** stand harvest (bit-close to the SCM, the
   mutant replay's ~3e-14 standard), NOT the cheaper per-step C-28 reconstruction. New
   per-RK-stage harvest fields are part of R0.

## R0 implementation decision — VALUE-ANCHORED reconstruction (2026-06-29, AD spike)

Discovered while building R0: the baseline-consistency gate (#3, ~3e-14) and the
coupled-tape worry both dissolve under a **value-anchored** reconstruction, the same
trick `deep_net` already uses for the z-linearisation:

```
light_active(z) = light_frozen(z)            // exact harvested env VALUE (3e-14 by construction)
                + [recon(z; θ) − value(recon(z; θ))]   // adds the θ-derivative, ZERO value
```

- `recon(z; θ)` = the C-27/C-28 trapezium `exp(−Σ_i C_i·area_leaf(θ,h_i)·Q(z/h_i))` over
  the **frozen** per-RK-stage stand (heights h_i, weights C_i = ce_i/area_leaf(θ₀,h_i)
  frozen), area_leaf active in the allometric trait.
- The added term has VALUE 0, so every cohort reads **exactly** the frozen env — the
  resident metric VALUES are bit-identical to the shipped frozen engine (the R0 gate is
  automatic). Only the DERIVATIVE channel θ→canopy→light is new.
- Consequence: cohorts stay **decoupled** given θ (each reads a light that depends on θ
  through all residents' area_leaf, but not on the other cohorts' replayed states). No
  coupled monolithic tape needed — the existing per-cohort / one-recording machinery is
  reused unchanged. The "coupled tape" the scope feared is unnecessary for the
  freeze-geometry design.
- First-order correctness: the focal cohort's own height→light channel is already carried
  by `deep_net`'s frozen z-linearisation `ld·(z−z₀)`; the recon term adds θ→light at frozen
  z. The dropped cross term d²light/(dz dθ) is second order, so dM/dθ is first-order exact.
- Which traits get a resident feedback: the canopy-light formula at FROZEN geometry is
  `competition(z) = Σ_i density_i·k_I·area_leaf(a_l1,a_l2,h_i)·Q(z/h_i; eta)`. Of the **28
  differentiable traits** (lma…a_dG2 — see `field_names()`), only **a_l1, a_l2** appear in
  it (via area_leaf). `k_I` and raw `eta` DO appear but are strategy-level CONSTANTS, not
  in the trait vector — so there is nothing to graft for them. ⇒ under freeze-geometry the
  graft of a_l1/a_l2 is **complete** for the differentiable set; the other 26 traits do
  not enter the light at fixed geometry, so resident == frozen for them *exactly* (not an
  approximation). (Correction to the R0-R1 commit message, which wrongly listed k_I/eta as
  a next step — they are not traits.)
- The genuinely all-28-with-feedback resident — where lma, a_p1, … also shift the canopy
  by reshaping heights/densities over time — requires the GEOMETRY to respond, i.e. the
  coupled re-evolving replay (deferred, R3+). Freeze-geometry deliberately drops it.
- R1 validates AD vs **FD over the same reconstruction** (perturb θ, recompute area_leaf,
  re-reduce) — so the reconstruction form (trapezium) IS the definition of the committed
  resident number; per scope decision #1's "census→resident".

API: `feedback = c("frozen","resident")` on `stand_gradient` (frozen = shipped path).

## R0 starting points (next session)

- New harvest: per-RK-stage stand state. `environment_history` already stores the env at
  each of the 6 stages; mirror it with the **stand state** (each cohort's height +
  area_leaf, all species) at each stage. Extend the `cache_ode_step` machinery in
  `patch.h` (alongside `stand_height_history`, which is per-step + species[0] only).
- The coupled stepper: a free function over `vector<Node<...,ad>>` + the frozen schedule
  (NOT `Patch<...,ad>` itself — its ctor eagerly runs `compute_environment`/`reset`); at
  each stage build the active light interpolator from the current active stand, each
  cohort reads crown light via `get_environment_at_height` analogue on the active spline.
- Reuse: C-27/C-28 reconstruction (`scripts/ad_self_shading_live.R`,
  `ad_self_shading_timeint.R`), `ff16_cashkarp_replay`, the odelia differentiable spline
  (#32, frozen knots + active values).
