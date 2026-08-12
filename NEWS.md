## Plant (development version)

### Breaking changes

These change the R-facing interface and require updating downstream code. Each
entry gives the `old -> new` migration; the `plant-update-interface` skill
(`.claude/skills/plant-update-interface/`) reads this section to migrate
products using plant.

* **`Leaf$set_physiology()` takes `root_network`, not `root_carbon_per_leaf_area`,
  and `Leaf()` no longer takes `beta_R_H` or `beta_R_V`.** Requires
  phylloptim >= 0.2.0 (phylloptim #33). Migration:
  `l$set_physiology(root_carbon_per_leaf_area = x, ..., soil_depth = d, ...)` ->
  `l$set_physiology(root_network = phylloptim::root_network_from_carbon(x, soil_depth = d, beta_R_H = ..., beta_R_V = ...), ..., soil_depth = d, ...)`;
  drop `beta_R_H`/`beta_R_V` from `Leaf()` calls. `RootNetwork()` builds a network
  from resistances directly, for a caller who has those rather than a carbon
  profile.

  The leaf's supply solve reads two vectors -- `r_R_H_min` and `r_R_V_sum` -- and
  nothing in it touches root carbon, the 1/3 : 2/3 root split, the layer thickness
  or either `beta_R_*` constant. Those are a root-ARCHITECTURE model, which is now
  `TF24_Strategy`'s: `beta_R_H` and `beta_R_V` are still TF24 parameters at the
  same values, and `net_mass_production_dt` calls
  `phylloptim::root_network_from_carbon` into a strategy member before each
  `set_physiology`. Same move `leaf_specific_conductance_max` already made -- plant
  computes `kmax` from height and passes a scalar, because which
  conductance-versus-height model is in force is not the leaf's business.

  **No TF24 results change.** Verified bit-identical over 18 operating points
  (6 heights x 3 soil-moisture profiles) x 26 states, rates and auxiliaries,
  against `develop`; phylloptim's own 288-point golden file is bit-identical too.
  `scientific_version` is therefore unchanged.

* **A census gradient's trait columns are named per species: `"lma"` -> `"1.lma"`.**
  Migration:
  * `colnames(stand_gradient(scm)$gradient)`  -> same, now `"<species>.<parameter>"`
  * `stand_gradient(scm, traits = "lma")`     -> `stand_gradient(scm, traits = "1.lma")`
  * `census_trait_names_tf24(scm)`            -> same, now prefixed
  * `g$gradient[, "lma"]`                     -> `g$gradient[, "1.lma"]`, or
    `g$gradient[, sub("^[0-9]+\\.", "", colnames(g$gradient)) == "lma"]` for every
    species' column

  Concatenating each species' parameter names with no prefix gave `S * P` columns
  with every name repeated `S` times. Character indexing resolves a name to its
  *first* match, so a multi-species gradient silently returned **species one's
  column** for every named parameter, and the unknown-parameter check could not
  see it because the name was present. A bare name now refuses, naming the
  convention. The prefix is applied at every species count, including one, so a
  single-species caller is not written against a shape that changes when a second
  species arrives. Found by asserting the columns are unique; it was live in three
  of this suite's own checks, one of which was comparing a quantity summed over
  both species against species one's column.

* **`stand_gradient_unanswered()` is now empty, so no trait is refused by name.**
  Migration: none mechanical, but the *behaviour* changed and a caller relying on
  the refusal will now get an answer. It listed thirteen traits the leaf supplied
  no derivative for; the leaf now supplies rows for all of them, and eleven come
  back live. The remaining two, `psi_crit` and `root_psi_crit`, are exactly zero
  at an interior operating point by complementary slackness -- they set the dry
  bound of a feasible interval the point is inside -- so they are declared zeros
  rather than refusals, and they carry the whole row at a pin. The refusal
  mechanism is kept, and matches on the parameter rather than the column, so a
  trait that loses its row is refused rather than reported as zero.

* **`run_stochastic_collect()`'s environment field is `env`, not `light_env`.**
  Migration: `out$light_env -> out$env`. The old name was never produced by
  anything — `StochasticPatch::r_get_state()` had its environment leg commented
  out and the collector asked for a name that did not exist — so every element of
  that field was `NULL` and nothing can have depended on its contents. Code that
  tests for the name, or indexes the collected list positionally, needs updating.
  `env` is what the patch reports and what `run_scm()`'s collected output already
  calls it. See the corresponding entry under Minor changes.

* **`TF24_Environment`'s `atm_kpa` driver now defaults to 101.3 kPa, not 100.5.**
  `scientific_version` for TF24 goes 7 -> 8, TF24f 7.1 -> 8.1. Migration: none
  required, but read this if you have TF24 results on disk.
  * The leaf model's ppm -> Pa conversion was the hard-coded constant `0.1013`, which
    is 101.3 kPa in disguise (`1e-6 * 101300 Pa`). The driver said 100.5, so the
    model's conductance side responded to 100.5 while Gamma*, Kc, Ko, Km and the ci
    root-find bounds silently assumed 101.3. Deriving the conversion from `atm_kpa`
    (phylloptim #15) made the model self-consistent and turned that disagreement into a
    **+2.4%** shift in output.
  * 100.5 came from `34d46ac2` (#446), an interface refactor that does not mention
    atmospheric pressure, with no recorded rationale, while every leaf-level test used
    101.3. Pinning the driver to the value the rest of the model already assumed
    reduces the **net** movement of this whole branch to **-0.20%**, and restores
    every pinned baseline to develop's own values -- including the seeded stochastic
    TF24 counts, which match exactly. (Those counts were subsequently re-derived when
    the stochastic arrival schedule started scaling with patch area — see Minor
    changes — so the file no longer shows the matching values. The equivalence that
    was established here still holds; it is just no longer the thing the test pins.)
  * `atm_kpa` remains a driver: set it per site if you are modelling altitude.
    `env$extrinsic_drivers_set_constant("atm_kpa", <kPa>)`.

* **The collar bracket is now clamped to `root_psi_crit`** (phylloptim #24, #584).
  `scientific_version` for TF24 goes 6 -> 7, TF24f 6.1 -> 7.1. The clamp compared a
  magnitude against a signed potential, so it could never bind and the solver
  optimised over a root-collar potential the root system cannot supply. The window is
  **1.2 MPa wide at the defaults** (`psi_crit` 7.0855 vs `root_psi_crit` 5.8703), so
  water-limited runs change; no mesic run does, and the standard SCM scenario is
  bit-identical. No migration: no name or signature changes.

* **Water potential now has one representation everywhere: positive magnitudes in
  MPa** (phylloptim #25). `scientific_version` for TF24 goes 5 -> 6, TF24f 5.1 -> 6.1.
  Migration:
  * `leaf$root_collar_psi_` -> `leaf$opt_root_psi_`, **and its sign flips**: it is
    now the positive magnitude. Renamed rather than reused deliberately, so an old
    script gets an error instead of quietly reading the wrong sign.
  * the **`opt_root_psi` aux changes sign** for the same reason, and now agrees with
    TF24f's `opt_root_psi_state`, which always held the magnitude. Before, the aux
    reported the signed potential while the state was negated back from it — an
    inconsistency in plant's own outputs. Any stored output or plot that negated the
    aux must stop doing so.
  * `l$E_from_Soil_to_Root_Collar(collar, psi_soil)`, `l$find_root_psi(...)`,
    `l$find_psi_stem_from_psi_root(...)`, `l$dE_from_soil_dpsi_collar(...)` and
    `l$transpiration_to_psi_stem(E, psi_upstream)` keep their signatures but now
    take **positive magnitudes**, and the bracket ends of `find_root_psi` swap
    (wettest layer first, `psi_crit` second). Passing the old signed values raises
    an error rather than returning a wrong number — the leaf package validates it.
  * outputs move by about **-3.2e-4 relative** (one-species SCM offspring production
    83.9026 -> 83.8761). Not an equation change: the rewrite is exactly
    sign-symmetric, but boost's TOMS748 iterates depend on the bracket's
    orientation, and the SCM amplifies the resulting 1-3 ULP. Every pinned test
    value and the exact stochastic TF24 counts pass unchanged.

* The TF24 leaf gas-exchange and hydraulics model now comes from the standalone
  header-only [`phylloptim`](https://github.com/traitecoevo/phylloptim) package instead of
  a copy in this repo (phylloptim #9). `plant::Leaf` is an alias for `phylloptim::Leaf`, so
  C++ consumers that `#include <plant/leaf_model.h>` and read public `Leaf` members
  keep compiling unchanged. **`scientific_version` for TF24 goes 4 -> 5 and TF24f
  4.1 -> 5.1: TF24 output moves by about +2.4%** — see `tf24_strategy.h` for the
  measurement and its attribution. Migration:
  * `l$set_physiology(area_leaf =, mass_root_prop =, rho =, a_bio =, PPFD =,
    psi_soil =, soil_depth =, leaf_specific_conductance_max =, atm_vpd =, ca =,
    sapwood_volume_per_leaf_area =, leaf_temp =, atm_o2_kpa =, atm_kpa =)`
    -> `l$set_physiology(root_carbon_per_leaf_area =, PPFD =, psi_soil =,
    soil_depth =, leaf_specific_conductance_max =, atm_vpd =, ca =, leaf_temp =,
    atm_o2_kpa =, atm_kpa =)`. **This is a semantic change, not just a rename:**
    `area_leaf`, `rho`, `a_bio` and `sapwood_volume_per_leaf_area` were dead
    stores and are gone, and `root_carbon_per_leaf_area` is the old
    `mass_root_prop` **divided by `area_leaf`**. The leaf is purely intensive now;
    uptake is exactly homogeneous in that ratio, so passing the old absolute
    carbon silently gives a root system too weak by a factor of `1/area_leaf`.
  * `l$area_leaf_`, `l$rho_`, `l$a_bio_`, `l$sapwood_volume_per_leaf_area_`
    (removed) -> no equivalent; all four were assigned and never read.
  * `plant::umol_per_mol_to_Pa` (C++, removed) -> `Leaf::umol_per_mol_to_Pa_`, now
    derived per call from `atm_kpa` rather than hard-coded at 0.1013.
  * `l$initialize_integrator(rule, tol)` keeps its name but now sets the tolerance
    of the leaf package's header-only adaptive Simpson quadrature rather than
    configuring plant's compiled QAG; only `transpiration_full_integration`, a
    spline-fidelity diagnostic, was ever affected.
  * Requires **odelia at master** (commit `d8235d1` or later), not just the
    `>= 0.2.0` in DESCRIPTION: `Patch::ode_rates` is non-const since #585 and
    odelia 0.2.0's `r_ode_rates` takes the system by `const&`. plant's `develop`
    does not compile against 0.2.0 either — this is not new here, but it will bite
    anyone building from a released odelia.

* Penman-Monteith leaf energy balance added to TF24/TF24f behind an
  opt-in gate, **default off** (#523). With the gate off, runtime behaviour and
  outputs are bit-identical to before (verified against leaf-level and SCM
  baselines) and `scientific_version` is unchanged — but the change adds fields
  to the R-facing `TF24_Pars` and a new environment driver, so any snapshot /
  round-trip test that encodes the full TF24 parameter set or driver list will
  see new entries. Migration:
  * `TF24_Pars` gains two fields (default off / inert): `pars$use_energy_balance`
    (0 = off = today's `Tleaf = Tair`; non-zero = on) and `pars$d` (characteristic
    leaf dimension, m, for the aerodynamic resistance). No action needed unless
    you assert the exact `pars` field set.
  * `TF24_Environment` gains a `wind_speed` extrinsic driver (default `2.0`
    m s⁻¹), set like any other driver:
    `env$extrinsic_drivers_set_constant("wind_speed", U0)` /
    `env$extrinsic_drivers_set_variable("wind_speed", x =, y =)`.
  * The `Leaf` submodel exposes `use_energy_balance_`, `d_`, `wind_speed_`
    (settable) and `Tair_`/`Rn_`/`ra_` (readable) for leaf-level experimentation.
  * To enable: set `pars$use_energy_balance <- 1` (and optionally `pars$d`,
    the `wind_speed` driver) before running; leaf temperature is then solved from
    the energy balance and fed to the Farquhar temperature scaling.

* Strategy biological parameters are now stored in a nested `pars` sub-object
  rather than as flat fields on the strategy (#410). This applies to every
  strategy type (`FF16_Strategy`, `K93_Strategy`, `TF24_Strategy`). Migration —
  read/write any biological parameter through `$pars`:
  * `s$lma`    -> `s$pars$lma`
  * `s$rho`    -> `s$pars$rho`
  * `s$hmat`   -> `s$pars$hmat`
  * `s$eta`    -> `s$pars$eta`
  * `s$k_I`    -> `s$pars$k_I`
  * `s$S_D`    -> `s$pars$S_D`
  * …and likewise for every other strategy parameter (all FF16/TF24 allometry,
    production, mortality and hydraulic parameters; all K93 `b_*`/`c_*`/`d_*`
    and `height_0`). Construction by flat field also moves under `pars`:
    `FF16_Strategy(lma = 0.1)` -> `s <- FF16_Strategy(); s$pars$lma <- 0.1`
    (or build from traits via `add_strategies()`).
  * Top-level strategy fields are unchanged: `control`, `birth_rate_x`,
    `birth_rate_y`, `is_variable_birth_rate`, `collect_all_auxiliary`.
  * A `<Strategy>_Pars` constructor is exported per strategy (`FF16_Pars()`,
    `K93_Pars()`, `TF24_Pars()`).
* The trait -> strategy interface was renamed to be clearer and pipe-friendly,
  with the `Parameters` object now the first argument (#410). The old names
  remain as deprecated shims for one release. Migration:
  * `expand_parameters(traits, p, hyperpar, birth_rate_list =, keep_existing_strategies =)`
    -> `add_strategies(p, traits, hyperpar =, birth_rate =, keep_existing =)`
  * `mutant_parameters(traits, p, …)` -> `add_mutant(p, traits, …)`
  * `strategy_list(traits, p, hyperpar, birth_rate_list)`
    -> `generate_strategy(p, traits, hyperpar =, birth_rate =)`
  * argument `birth_rate_list` -> `birth_rate`; `keep_existing_strategies` ->
    `keep_existing`.

* SCM cohort-refinement now happens entirely in C++, and the R interface is
  consolidated onto a single `run_scm()`. `run_scm_collect()`, `run_scm_error()`
  and `build_schedule()` have been removed (#408, #459, #462). Migration:
  * `build_schedule(p, …)` -> `scm <- run_scm(p, …, refine_schedule = TRUE)`,
    then read the refined parameters as `scm$parameters`. The former
    `attr(p_new, "offspring_production")` side-channel is now
    `scm$offspring_production`.
  * `run_scm_collect(p, …)` -> `run_scm(p, …, collect = TRUE)`
  * `run_scm_error(p, …)` -> set `scm$collect_errors <- TRUE` then read
    `scm$combined_node_errors`
* The argument order of `run_scm()` changed; `use_ode_times` is now the last
  argument (after the new `refine_schedule` and `collect`).
* Numeric `Control` defaults are now set in the C++ `Control()` constructor
  (the pragmatic "fast" settings), and the two R preset helpers were removed
  (#463). Raw `Control()` now means **fast**, not accurate. Migration:
  * `fast_control()` -> `Control()` (or its lowercase alias `control()`)
  * `scm_base_control()` -> `Control()` / `control()`
  * for the old tight-tolerance behaviour -> `control_accurate()`
  * `scm_base_parameters()` is unaffected (it builds a `Parameters`, not a
    `Control`).
* Removed the TF24/TF24f parameter `hk_s` from the public interface and from
  the `Leaf` constructor; hydraulic cost now depends on `g1_TF24`, `beta2`, and
  the vulnerability curve terms only. Migration:
  * `TF24_Strategy()$pars$hk_s` -> removed (no equivalent)
  * `TF24f_Strategy()$pars$hk_s` -> removed (no equivalent)
  * `make_TF24_hyperpar(..., B_hks1 =, B_hks2 =)` -> adjust `g1_TF24` ~ `rho` scaling (new; no `hk_s` equivalent)
  * `make_TF24f_hyperpar(..., B_hks1 =, B_hks2 =)` -> adjust `g1_TF24` ~ `rho` scaling (new; no `hk_s` equivalent)
  * `Leaf(..., hk_s =)` -> remove `hk_s` argument
* The ODE solver and interpolator core were spun out into the standalone
  [odelia](https://github.com/traitecoevo/odelia) package, which plant now
  depends on (#456, #464). Plant's internal solver/interpolator headers
  (`inst/include/plant/ode_solver/*`, `interpolator.*`, `tk_spline.*`) were
  removed in favour of `odelia::ode::Solver` / `odelia/interpolator.hpp`.
  Migration: code that linked against plant's C++ numerics headers, or used the
  old `OdeRunner`/interpolator R6 types directly, must switch to odelia. Pure-R
  callers of `run_scm()` / `run_stochastic()` are unaffected.

The following earlier changes (since the 2.0.0 release) are also breaking and
were not previously recorded here:

* All fitness/equilibrium functionality was removed from plant; it now lives in
  the separate `regnans` package (#388). Removed `fitness_landscape()`,
  `solve_max_fitness()`, `viable_fitness()`, `fundamental_fitness()`,
  `assembly_parameters()`, `equilibrium_birth_rate()` and the `equilibrium_*`
  parameters/controls.
* The `Environment` object was removed from `Parameters`; pass it directly as
  the `env` argument to run functions, e.g.
  `run_scm(p, env = Environment("FF16"))` (#315). Likewise `Control` was removed
  from `Parameters` and is passed as the `ctrl` argument (#314).
* SCM & environment interface simplified (#446): environment construction is now
  the single `Environment("FF16")` constructor — removed `make_environment()`,
  `FF16_make_environment()`, `FF16_fixed_environment()`, `test_environment()`,
  and the helpers `scm_patch()`, `make_scm_integrate()`, `scm_state()`,
  `patch_to_internals()`, `scm_to_internals()`, `species_to_internals()`,
  `first()`, `last()`, `modify_list()`. State collection moved to C++ and
  `tidy_patch` runs by default.
* Node/Species/SCM competition methods renamed (#448):
  * `Species$competition_effects` -> `Species$compute_competition_effect_by_nodes`
  * `Species$competition_effects_error` -> `Species$compute_competition_effect_by_nodes_error`
  * `SCM$competition_effect_error` -> `SCM$compute_competition_effect_error_by_node_for_species_i`
  * the `Node$competition_effect` getter was removed.
* The exported-function surface was reduced (#420): the per-model R6 generators
  (`FF16_*`, `K93_*`, `TF24_*` `Node`/`Patch`/`SCM`/`Species`/`Stochastic*`) and
  many utilities (`Internals`, `OdeControl`, `bounds`, `check_bounds`,
  `clamp_domain`, `nlsolve`, `splinefun_log`/`splinefun_loglog`, `strategy`,
  `strategy_default`, `validate`, …) are no longer exported. Rename:
  `make_transparent()` -> `util_colour_set_opacity()`.
* The `Cohort` class/concept was renamed `Node` throughout — classes, methods,
  filenames, and R-side `coh` variables (#335).
* "Seed rain" terminology renamed for clarity (#297):
  * `seed_rain[_in]` -> `birth_rate`
  * `seed_rain[_out]` -> `net_reproduction_ratio` / `offspring_production`
  * `add_seeds()` -> `introduce_new_cohort()`
* Canopy/light renamed and generalised (#384): the `canopy` class became
  `Resource_spline` (`canopy.h` -> `resource_spline.h`); `compute_canopy` ->
  `light_availability$compute_environment`; `canopy_light_*` controls ->
  `light_availability_spline_*`; `get_canopy_at_height()` / `canopy_openness()`
  -> `get_value_at_height()`. The standalone `assimilation` class and its
  adaptive-integration options were removed.
* Extrinsic-driver API reworked (#340): per-driver `Environment` methods replaced
  by an `ExtrinsicDrivers` object at `env$extrinsic_drivers`, with
  `evaluate("rainfall", t)` / `evaluate_range("rainfall", c(...))`.
* ODE-stepper methods renamed (#413): `advance()` -> `advance_adaptive()`;
  `node_schedule_ode_times()` -> `ode_times()`. The corresponding `Parameters`
  field was also renamed: `p$node_schedule_ode_times` -> `p$ode_times`.
* Growth-optimisation routine generalised (#382):
  `FF16_solve_max_size_growth_rate_at_height()` ->
  `optimise_individual_rate_at_size_by_trait()` (new `size` / `size_name` args),
  with a height wrapper `optimise_individual_rate_at_height_by_trait()`.
* Removed the `FF16w` strategy, superseded by TF24 (`FF16w` -> `TF24`) (#438),
  and the `FF16r` strategy (#439); removed the experimental soil component from
  `FF16` (#441).
* Replaced the Leaf cost parameter `g1_TF24` with the Medlyn stomatal-conductance
  parameters `g0` (default 0.022) and `g1` (default 2.57) on the TF24 `Leaf`
  (#450, #451).
* `Disturbance` refactored into a `Disturbance_Regime` base with
  `No_Disturbance_Regime` / `Weibull_Disturbance_Regime` subclasses (#301);
  disturbance configuration moved from `Environment` to `Patch` (#290); the
  light-extinction coefficient `k_I` moved from `Parameters` onto the strategy
  (#293, #302).
* The stochastic runner's node-schedule interface was renamed to mirror `SCM`
  (#217, #506). Migration:
  * `StochasticPatchRunner$schedule`              -> `StochasticPatchRunner$node_schedule` (get and `<-` assignment)
  * `StochasticPatchRunner$set_schedule_times(x)` -> `StochasticPatchRunner$set_node_schedule_times(x)`
* `run_stochastic_collect()` no longer returns an `offspring_production` element
  in its result list — it is not defined for the finite-population model and was
  never populated (#498, #506). Migration:
  * `run_stochastic_collect(...)$offspring_production` -> removed (no equivalent)
* `Patch` computes its rates when they are read rather than when its state is
  set, so reading `$ode_rates` now evaluates and the separate compute entry point
  is gone (#585). Migration:
  * `patch$compute_rates()` (removed) -> `patch$ode_rates`, which computes at the
    state currently loaded and returns the result
  * Semantic change: `patch$ode_rates` was a cheap read of stored values and is
    now a full right-hand-side evaluation. `StochasticPatch$compute_rates()` is
    unchanged.

### New features

* **Every reduction over the size distribution now integrates over the
  coordinate that distribution is a density in.** Scientific version: FF16 and
  K93 1 -> 2, TF24 8 -> 9.

  A reduction over the size distribution is a quadrature of a density, so its
  trapezium widths are gaps in the coordinate the density is carried in. Three
  of them took gaps in height whatever the coordinate was: the fused
  value-and-slope reduction that builds the light field, the census, and the
  transposes of both resource reductions. `Species::compute_competition()` and
  `Species::consumption_rate()` already followed the coordinate, so the fused
  reduction and its own value accessor disagreed by 100% relative on a
  birth-date run while agreeing bit for bit on a height run.

  **This moves results on the birth-date coordinate and on no other.** Measured
  on the full-lifetime deep-crown FF16 anchor, offspring production is
  16.884586 over heights -- unchanged to the last bit -- and 17.172004 over
  birth dates, i.e. **+1.70%**. The height-coordinate arithmetic is bit-identical
  because the abscissa is minus the height and negation is exact.

  Two consequences worth knowing. The fused reduction and the value reduction
  now agree bit for bit on both coordinates. And on the birth-date coordinate
  the closing interval runs from the youngest cohort's birth date to the current
  time, so the field is a function of the time as well as of the state: a field
  left behind by an earlier stage and read at a later clock differs from a fresh
  build. `Species::consumption_rate()` has always had this property.

  `integrate_over_size_distribution()` genuinely does integrate over height, so
  it now orders its grid before integrating rather than assuming the rows arrive
  in height order.

* **`Control$node_density_in_birth_date`** (default `FALSE`) carries the SCM's
  size distribution as a density in birth date instead of in height.

  The transport equation's compression term is the total derivative of the growth
  rate along a cohort's own trajectory, which equals `∂g/∂h` only when growth is
  a function of size. TF24's reserve gate breaks that: the finite-difference
  probe in `Node::growth_rate_gradient` moves height while holding *absolute*
  carbon fixed, so it shifts the reserve fraction `r = S/S_max`, whereas a cohort
  actually grows with `r` roughly constant. The probe is accurate about a
  quantity the plant never experiences, so this is a different derivative rather
  than a worse approximation of the right one.

  In birth-date coordinates the density rate is mortality alone (nothing moves an
  individual along the birth-date axis), the birth density is
  `birth_rate·pr_estab` with no division by the growth rate, and both resource
  integrals run over introduction times.

  With the flag off the solver is bit-identical to before. With it on:

  - FF16 and K93 have size-only growth, so both coordinates must converge to the
    same answer, and do. Relative difference in offspring production over
    successive halvings of the node spacing: FF16 `1.0e-2`, `4.2e-3`, `1.2e-3`,
    `2.9e-4`; K93 `2.6e-3`, `5.9e-4`, `1.4e-4`, `3.5e-5`. That is the ~2nd order
    of the trapezium rule, and confirms the coordinate change is a change of
    coordinates.
  - The gap is almost entirely the *height* coordinate's error. Across those same
    four schedules FF16's birth-date answer moves 18.6956 → 18.6996 while its
    height answer climbs 18.5071 → 18.6941 toward it. At the default schedule the
    birth-date answer is already within `2e-4` of converged and the height answer
    is `1.0e-2` away — roughly 50x more accurate for the same number of nodes.
  - TF24, whose growth is not a function of size alone, is the case the two
    coordinates genuinely disagree on, and the diagnostic is that refining the
    schedule does *not* close the gap the way it does for FF16 and K93. At
    `lma = 0.0825`, `hmat = 5`, `birth_rate = 20`,
    `max_patch_lifetime = 30`, over the same three schedules: height
    474 → 587 → 696, birth date 3495 → 3729 → 3840, a ratio of 7.4 → 6.3 → 5.5.
    Both are still moving, so neither figure is a converged value; the point is
    that they are not converging *to each other*, which is what a wrong
    compression term looks like as against a coarse quadrature.

  Multi-species runs behave the same way. The established two-species FF16 and
  three-species K93 cases converge per species at the same ~2nd order (FF16
  `4.4e-3`/`7.1e-3` → `2.8e-4`/`6.0e-4`; K93 all three species `3.1e-3`–`6.2e-3`
  → `1.9e-4`–`3.7e-4`), and K93's competitively marginal first species — ~90x
  below the dominant one — converges no worse in relative terms than the
  dominant ones. So the shared competition profile does not disadvantage a
  marginal species under the coordinate change.

  **TF24 two-species changes the ecological outcome, not just the numbers.** On
  the two-species case of `test-strategy-tf24.R` (`lma` 0.0825 / 0.10,
  `max_patch_lifetime = 30`), the height coordinate excludes the slower species
  (offspring production `1.4e-4` against `503` for the faster) while the
  birth-date coordinate has them coexisting at comparable abundance (`1004`
  against `2707`), and refining the schedule does not move either toward the
  other. The conclusion recorded in that test — that reserve-gated growth
  (#517) largely excludes the slower species — is therefore coordinate
  dependent, and needs re-deriving before it is relied on. *(Re-derived and
  retracted; see the entry under Minor changes & bug fixes below.)*

  `Node::growth_rate_gradient` is not called, which also removes one leaf solve
  per cohort per Runge-Kutta stage.

  **What R receives is unchanged in meaning.** The coordinate is an internal
  choice, so `log_density` is converted back to a density in *height* before it
  leaves C++ — `Species$log_densities`, `Patch$state`, and therefore
  `run_scm(collect = TRUE)`, `tidy_outputs.R`'s `density = exp(log_density)`,
  `interpolate_to_heights()` and the plots all keep their existing meaning. The
  conversion is `N = ν / |dh/dτ|`, with the Jacobian formed by central
  differences of adjacent node heights against their birth dates (the boundary
  node supplies the extra point at the young end). Measured against an otherwise
  identical height-coordinate run, the median node agrees to `1.7e-3` in log
  space — 0.17% in density — improving to `4.4e-4` and `1.1e-4` over two
  schedule refinements. See Known issues for where this is weakest.

  The boundary node needs no differencing at all: `dh/dτ = −g(H₀)` at birth, so
  `compute_initial_conditions()` now records the birth growth rate on every node
  (`Node$growth_rate_at_birth`) and the Jacobian uses it there directly. It is
  exact *and* current for that node, because the boundary node is re-evaluated
  every step; for an introduced node the recorded rate is frozen at its own
  birth and so cannot serve as its present-day Jacobian. This removes a 31%
  error on the boundary node — it sits a whole introduction interval from its
  neighbour, which is the worst case for a one-sided difference. Recording it
  also lets the two coordinates' boundary conditions be checked against each
  other: `exp(log_density) · g(H₀)` on the height path must equal
  `birth_rate · pr_estab`, which the birth-date path carries directly, and that
  is now a test.

  The quantity actually integrated is reported alongside rather than lost:
  `Species$log_densities_state` (and a `log_density_state` row in `Patch$state`,
  present only on the birth-date path), plus `Species$height_jacobian` for the
  conversion itself. `Node$log_density` stays unconverted, since forming
  `|dh/dτ|` needs the neighbouring nodes and a `Node` does not have them.
  `export_patch_state()` resumes from `patch$ode_state`, which is the raw state,
  so resume is unaffected by any of this — pinned by a test.

  Two invariants the birth-date axis needs, which the height axis did not:

  - The boundary node's birth date is the current time, so
    `Patch::compute_environment()` refreshes it before building the profile.
    `compute_rates()` (which stamps it) runs *after* the `set_ode_state()` that
    rebuilds the environment, so reading the stamp there would use the previous
    derivs call's time. The measured effect on FF16 offspring production is below
    `1e-6` — the boundary node carries almost no leaf area — but it made the
    spatial quadrature a function of the ODE step size.
  - Introduction times are the quadrature grid, so repeated ones span zero width
    and drop out of the integrals. A scheduled run cannot produce them; a patch
    seeded or resumed *without* per-node times gives every node the boundary
    node's birth date, which would silently zero the competition profile.
    `Patch::check_birth_dates_distinct()` now rejects that with an actionable
    message. A faithful `export_patch_state()` / `set_initial_state()`
    round-trip carries the times and is unaffected.

  `Species::compute_competition()` keeps the sorted-grid fallback of #574 for the
  height coordinate only: that path integrates in height, so sending the
  birth-date coordinate down it would swap coordinates mid-run. The birth-date
  abscissa cannot invert (introduction times are fixed at birth), but the
  *height* early exit is skipped when the height ordering has broken, since a
  node below the query height can then be followed by a taller one.
* **A dry TF24f patch no longer aborts the whole run on the ci root-find.**
  `Leaf::dprofit_droot_collar_psi` — TF24f's exact AD/IFT gradient — called
  `psi_stem_to_ci()` before testing for hydraulic shut-down. In shut-down,
  `psi_upstream >= psi_stem` (both positive magnitudes), so
  `gc = const · transpiration` goes negative, the residual stops crossing zero
  over (Γ*, ca], and TOMS748 throws *"a and b do not bracket the root"* rather
  than returning non-finite — which is what the existing `isfinite` guard on the
  next line was written to catch, and cannot see. One dry patch therefore killed
  the run. The condition is now tested *before* the solve, matching what
  `set_leaf_states_rates_from_psi_stem()` has always done, and returns a zero
  gradient. Reproduced at 5 soil layers, θ = 0.005–0.03 with 1 m yr⁻¹ rainfall
  (ψ_stem = 1.23 against ψ_upstream = 5.92 MPa). Behaviour changes only in states
  that previously threw, so no scientific version bump.
* **Adaptive interpolation now names a non-finite target instead of blaming
  resolution.** `AdaptiveInterpolator::check_err()` compares against NaN, and
  every NaN comparison is false, so a single NaN or Inf from the target made its
  interval permanently unacceptable: refinement halved the spacing to
  `max_depth` and then reported *"Interpolated function as refined as currently
  possible"*. That message cost real debugging time on a TF24 patch. A
  non-finite value now fails immediately, naming the point and the value; the
  resolution-limit message additionally reports the spacing reached, the limit,
  and the tolerances missed. A new `test_adaptive_interpolator()` test hook
  drives the refiner from R (the production caller passes a C++ lambda).
* **The scenario scorecard now reports `persists`.** Whether a strategy replaces
  itself, R0 >= 1, which is a different question from whether the run completed.
  It matters because `status`/`outcome` test `total > 0`: at the current baseline
  five of the eight hydraulic scenarios return R0 between 2e-15 and 6e-14 —
  numerically extinct — and were all recorded as `"persisted"`. Only **1 of 8**
  scenarios persists. That is the main reason the gateway discriminates so little
  (3/8, with no expected failure failing). Reported alongside, not folded into,
  the existing classification: the scenario CSV's "Model failure" means the model
  *breaks numerically* (#549/#550), not that the strategy dies out, and the
  blessed baseline diff is defined on the existing columns.
  `scenario_summary()` gains `n_persists` and tolerates scorecards recorded
  before the column existed.
* **A hydraulically shut-down plant no longer draws water from the soil
  (`TF24@v4`, `TF24f@v4.1`).** Both shut-down exits in
  `Leaf::find_root_collar_psi` set `profit_` directly and bypass
  `profit_psi_stem_TF`, and the first returns before any
  `E_from_Soil_to_Root_Collar` call in that solve. Because `Leaf` is a value
  member reused across every `compute_rates` call for an individual, every leaf
  output they did not assign kept the **previous step's** value. That was not
  just a reporting problem: `soil_consumption_` feeds
  `TF24_Strategy::evapotranspiration_dt` and hence the patch water balance, so a
  plant whose stomata had closed carried on extracting its last wet-step uptake.
  Measured: an individual moved from θ = 0.30 to θ = 0.02 (past ψ_crit) reported
  `E_up_` and `transpiration` *identical* to its wet step. Both exits now zero
  the transport chain (`transpiration_`, `stom_cond_CO2_`, `E_up_`,
  `soil_consumption_`), and recovery on rewetting is tested. Note the water
  budget still *closed* in the buggy state — what was recorded as depleted was
  what was removed — so the conservation tests could not catch this; only the
  physics was wrong. Because water-limited runs change (scenario gateway:
  offspring production moves by up to 5e-3 relative on 5 of 8 scenarios, with
  every success/failure classification unchanged), the TF24 scientific version
  is bumped **3 → 4** and `TF24f` tracks to **4.1**. This invalidates `logpile`
  caches for water-limited TF24 runs, which is the intended, safe direction.
* **Root hydraulic parameters are now settable from R.** `root_b`, `root_c`,
  `root_psi_crit` and `rooting_depth_max` move into `TF24_Pars` (they were fixed
  members of `TF24_Strategy` and a file-static constant in
  `src/tf24_strategy.cpp`, so unreachable from R). Root shutoff was pinned at
  ψ ≈ 5.87 MPa, too conservative for taxa that operate below it (e.g. *Acacia
  aneura*), and rooting depth at 1.5 m. **Defaults are unchanged**, so this is
  an interface change only and the TF24 scientific version does *not* move.
  Two cautions: `root_psi_crit` is derived from `root_b` and `root_c` exactly as
  `psi_crit` is from `b` and `c`, so set it whenever you set those; and
  `rooting_depth_max` beyond the soil column depth (`TF24_Environment$depth`,
  default 1.5 m) gains nothing, because the layers do not exist.
* **`check_driver_interpolation()`** reports how a driver's control points
  survive interpolation — evaluated range, fraction of the series that goes
  negative, undershoot area, and the interpolated vs supplied integral — and
  warns when a non-negative series interpolates negative. Extrinsic drivers use
  a cubic spline, which is a poor fit for intermittent forcing: a realistic
  daily rainfall series with a ~10% wet-day fraction evaluates negative at ~45%
  of points, reaching −5.7 m yr⁻¹. Worth running on any site-forcing series
  before committing to a long run.
* **`assimilation` is now reported for TF24/TF24f.** The auxiliary variable was
  listed in `aux_names()` but never written — no index was resolved and no
  `set_aux` call referenced it — so the slot reported whatever `Internals`
  happened to hold, and carbon uptake was unavailable as a model output. It now
  carries net CO₂ assimilation at the optimal operating point, per unit leaf
  area (µmol CO₂ m⁻² s⁻¹). **Net, not gross:** `Leaf::assim_colimited()`
  subtracts dark respiration, so gross = `assimilation` + R_d with
  R_d = 0.015·vcmax at the acclimated vcmax. It is integrated over the crown
  under `deep-crown` shading alongside the other leaf outputs, and the two
  hydraulic shut-down exits now set it explicitly (to −R_d) instead of leaving a
  stale probe value, keeping `profit == assimilation − hydraulic cost` true in
  every branch.

* **NSC storage pool for TF24 (`TF24@v3`, `TF24f@v3.1`).** TF24 now carries a
  non-structural-carbohydrate storage state so growth and mortality respond to
  *buffered* carbon rather than instantaneous net production (#517, #554).
  Growth/reproduction are reserve-gated (a smooth logistic on relative reserves
  `r = S/S_max`), and mortality is now `a_dG1·exp(-a_dG2·r)` — bounded in
  `[a_dG1·e⁻ᵃ_dG2, a_dG1]` — which is the root-cause fix for the SCM
  cohort-density blow-up (#550). New parameters `a_st1`/`a_st2`/`a_st3` (storage
  capacity per unit sapwood, growth half-on reserve fraction, birth fill).
  Because this changes the simulation output for identical inputs, the TF24
  scientific version is bumped **2 → 3**; `TF24f` inherits the state and its
  compound version auto-tracks to **3.1**. This invalidates the `logpile` cache
  for both models. See the staging guide in `overstorey-staging/guides/`.
* **Per-model scientific versioning.** Each model now carries a scientific
  version — an integer that is independent of the package `Version` and is
  bumped only when the model's equations or default parameters change the
  simulation output for identical inputs (not for refactors, performance, or
  interface changes). Read it with `model_version("FF16")` (an integer) or
  `model_id("FF16")` (`"FF16@v1"`). The number is authored as the
  `scientific_version` constant on each strategy class in
  `inst/include/plant/models/*_strategy.h`, next to the equations it versions.
  Downstream tools (e.g. `logpile`) use it to decide when archived simulations
  must be re-run: reruns follow scientific changes, not every software release.
  A drift-guard test (`tests/testthat/test-model-version.R`) fails when a
  model's default parameters change without a version bump. Starting versions:
  `FF16@v1`, `K93@v1`, `TF24@v2` (a published result used the pre-versioning
  "v1" science). `TF24f`, being a fast *approximation* of TF24, carries a
  compound version `"<TF24 version>.<approximation revision>"` (`TF24f@v2.1`):
  the major component auto-tracks TF24 so a TF24 change also invalidates TF24f,
  and the minor tracks changes to the approximation itself.
* `run_scm()` can start a patch from **pre-existing nodes** rather than always
  growing from empty — to resume an exported run or to seed an arbitrary initial
  size distribution at patch age 0 (#499, revives #304). New
  `export_patch_state(scm, step)` captures a patch's full state; the initial
  condition rides on the `Parameters` object (`initial_state`,
  `n_initial_cohorts`, …) so `run_scm(p)` just works and stays serialisable.
* Selectable **canopy shading models** for FF16/TF24 via `control$shading_model`
  (#490, #417): `""` (each strategy's own default), `deep-crown`, `mean-light`,
  `crown-centre`, `flat-top-soft-box`, `flat-top-box`, `ppa` (with
  `ppa_layer_optical_depth` / `ppa_layer_smoothing`). Dispatch is bound once in
  `prepare_strategy()`; deep-crown is bit-for-bit unchanged.
* Alternative **fixed-step forward-Euler** ODE integration alongside the
  adaptive Cash-Karp solver, selected with `control$fixed_time_step` (years;
  `0` = adaptive RKCK, the default) (#489). Combining `fixed_time_step > 0`
  with a mutant run / `save_RK45_cache` / `use_ode_times` errors clearly.
* Multi-layer **root water uptake** for the TF24 strategy (soil→root→stem→leaf
  hydraulic pathway following Potkay et al. 2021), replacing the previous
  single-value soil-water treatment (#488).
* Added the **TF24 strategy** — a leaf-level water-use/hydraulics model with a
  `TF24_Environment` whose soil-layer count and parameters can be set after
  construction (#445), including a Medlyn stomatal-conductance model on the leaf
  (`solve_medlyn_ci_*`, `medlyn_model_gs`) (#450).
* New **mutant-fitness method**: caches resident environments at each ODE step so
  mutant fitness reuses the residents' resource shadow. Adds a `save_history` /
  `save_RK45_cache` control, `environment_history` / `patch_step_history` on the
  SCM, and a `mutant_parameters()` method (#362, #379).
* Per-species **time-varying birth rates** via the extrinsic-driver mechanism:
  `set_constant_birth_rate(p, i, k)`, `set_interpolated_birth_rate(p, i, x, y)`,
  and `Parameters` fields `birth_rate_x` / `birth_rate_y` /
  `is_constant_birth_rate` (#334).
* `run_scm()` collected output now exports the extrinsic environment drivers per
  timestep (`out$env[[t]]$canopy`, `…$rainfall`); new drivers export
  automatically (#347).
* `expand_state()` generalised to apply to any strategy (#443).
* New strategy parameter `recruitment_decay` — exponential decline of
  establishment probability with patch age (#330).
* Per-individual `consumption_rates` (water extraction across cohorts/soil
  layers) feeding the water-enabled environment's soil balance (#329); FF16
  rainfall getters/setters + spline exposed to R (#324); environment auxiliary
  variables (#323) and environment state variables (#305) exposed/added.
* Added an HTML report (plots + analyses) for the FF16 strategy (#350).

* **The census gradient has a reference over a whole trajectory for the first
  time: a forward tangent of the same run, stepped at the sizes it recorded.**
  `SCM::census_trait_tangent<Metrics>(direction, value)`, exported for the ladder
  as `ladder_trajectory_tangent_tf24(scm, direction)`. One weight per gradient
  column, so a unit vector returns one exact column of the census-by-trait
  Jacobian and a mixed one returns a contraction.

  It replays the recorded step **sizes** rather than the times, and rather than a
  controller of its own. A tangent run left to choose its own steps differentiates
  the controller, which the model does not contain; and a size differenced back
  out of two recorded times is not the size that was taken, since
  `fl(fl(t + h) - t) != h`. The replay lands on the run's own census to
  **3.2e-12** relative, which is what says it is the same function, and it is the
  comparison's floor.

  The reference traverses both reductions, the stage recursion and the
  introduction boundary, and none of the transposes under test are on its path.
  With the two defects below closed, the whole 3 x 88 gradient agrees with it to
  **3.7e-06** relative to the matrix, and per column to 1.2e-04 or better for
  every column the inflow boundary node does not carry.

* **The inflow boundary node now carries its own row in both reduction
  transposes, and the soil rows are right for the first time.** Both forward
  reductions integrate a trapezium from the boundary node up -- it is the size
  distribution's lower grid point -- while both transposes scattered onto the
  introduced nodes only, because the boundary node holds no ODE row to scatter to.
  Its draw and its shading do depend on the state, so the rows the reductions feed
  read low by that node's share.

  A quadrature node is not a term that can be dropped: its weight is shared with
  its neighbour, so the transposes were getting the neighbour's weight right and
  losing the node's own integrand. Both now scatter over `size() + 1` slots, the
  boundary node last within its species, and the block loop sweeps it too --
  `Patch::reduction_node_count()` is the index space that widening created.

  The two rows it has no ODE slot for are read out rather than dropped in silence:
  `Patch::boundary_node_adjoint` carries its height, which birth size imposes
  passive, and its density, which is the inflow condition. Naming them is what
  keeps the remaining channel measurable.

  Measured against both references, which agree with each other, on one
  right-hand-side evaluation:
  * the soil-rate rows against the tangent: **1.06e-05 -> 5.38e-08**, and
    1.14e-05 -> 4.67e-08 on the one-cohort patch;
  * the whole-state contraction against a plain-double difference:
    **3.01e-06 -> 8.91e-07**;
  * cohort rows and trait rows off the soil block unchanged at round-off
    (1.99e-16 and 2.39e-16).

  What is left at 5.4e-08 is the leaf's own supplied-row gap -- a difference
  re-solves the collar where the transpose grafts -- and not a missing channel.

  **One index hazard is worth recording, because it is the failure report 02 §4.8
  predicts for this shape.** Two index spaces now exist, the reduction's grid
  points and the ODE rows, and they differ by one node per species.
  `Patch::offspring_adjoint` seeded through the first while indexing the second,
  which put every fecundity seed in the wrong slot: the cohort-row residual went
  from 2e-16 to **3.5e+04**. Not a wrong number but a wrong address, exactly as
  predicted, and caught by a check that was already there.

* **Both faults the trajectory rung has to inject are injected, and both are
  caught only by the reference.** A rung whose faults have not been injected has
  not been climbed, so these are the evidence rather than a by-product. Each is a
  deliberately broken build, measured against the four-node stand and reverted.
  * **An active value held across a tape clear** -- the census recording's patch
    and tape hoisted out of the per-metric loop, which is the tempting economy,
    since rebinding is not free. Caught at **169x** the contraction's tolerance
    and up to **5264x** on a coordinate column.
  * **The lost stage term** -- the reverse pass's sum over every later stage
    replaced by the immediate successor alone. Caught at **363x** on the
    contraction and up to **3333x** on a column, with the residual at 1.0
    relative: the columns are not degraded but wrong.

* **Three test files were still written against the height coordinate, and now name
  the coordinate they mean; seven more fail for two other reasons.** Carrying the
  size distribution in birth date by default left old references stale, and the
  earlier re-blessing covered `test-strategy-ff16.R` and `test-patch.R` only. The
  ten failures separate into three causes, and only the first is a test-side fix.

  **Fixed -- a test that took its coordinate from the default.** `test-species.R`
  (45 failures), `test-node.R` (9) and `test-density-coordinate.R` (3). Every
  reference in `test-species.R` is a trapezium over node *heights*, so a
  `height_coordinate_strategy()` helper names that and the twelve construction
  sites use it. `test-node.R`'s ODE-interface reference now reads
  `control$node_density_in_birth_date` and builds the matching boundary condition
  and density rate -- `log(pr_estab)` against `log(pr_estab / g)`, mortality alone
  against mortality plus the compression term -- so it documents that the two
  coordinates are different functions instead of pinning one. Where a bare
  `Control()` stood for *the height coordinate*, it now says so.

  **Not a test bug -- seeding a patch from a size distribution is unavailable on
  the default coordinate.** `test-initial-state.R` (3) and `test-canopy-methods.R`
  (1) fail because `make_initial_state` gives every node the same birth date --
  `test-initial-state.R` asserts exactly that, `all(node_times == 0)` -- and the
  birth-date coordinate's guard rightly refuses it, since such nodes span zero
  width in its reductions. The guard's own advice is to supply
  `parameters$initial_node_times`, which the state exporter does not. **This is a
  forward-model gap, not a stale reference**, and pinning the tests to the height
  coordinate would hide it.

  **Stale numbers the coordinate moved, left un-re-blessed deliberately.**
  `test-strategy-k93.R` (2): offspring production 0.07545 against a recorded
  0.0753261, and `c(0.00254, 0.23262, 0.22014)` against
  `c(0.00254, 0.23215, 0.21944)`. **K93 is scientifically frozen**, so a 0.16
  percent move in its output is a question about whether the default belongs to a
  frozen model rather than a number to re-bless. `test-schedule-build.R` (1) 147
  schedule times against 148, `test-scm-support.R` (1) 17 auxiliary columns against
  16, `test-stochastic-patch-runner.R` (1) three stems alive against two, and
  `test-stochastic-patch.R` (3) an all-missing `min()` on an empty patch. Each needs
  a decision about what the right answer is before its reference is rewritten.

* **Rung 5's marginal-recruit limit is built, and it needed no new C++.** The check
  was skipped on the premise that nothing reports a boundary node's establishment
  probability or the production driving it. Both are already reachable --
  `new_node$individual$establishment_probability(env)` and its
  `net_mass_production_dt` aux -- so the premise was wrong.

  Its fixture assertion is derived rather than chosen. Establishment is
  `P^2/(P^2 + k^2)` times a decay in patch age with `k = a_d0 * a_0`, so the
  derivative peaks at `P = k/sqrt(3)`; at the shipped constant this fixture's
  recruits sit at **36 and 21 times** that production with establishment 0.998,
  which is the flat region where the limit would pass on the signal being small.
  Scaling `a_d0` by 30 puts them at **1.48 and 1.05** of the peak.

  From there the gradient is finite at every point of a sweep spanning **2000x** in
  establishment (0.423 down to 2.05e-04), and the establishment decay's column --
  which reaches the census through the boundary condition and nothing else -- falls
  monotonically from **1.24e-01 to 6.15e-05**, linearly in establishment. `a_d0`'s
  own column falls faster, roughly as `P^2`, which is the second-order contribution
  a recruit that cannot pay for itself is supposed to make.

* **The sweep can be split and resumed, and it is bit-identical when it is.**
  `census_trait_gradient(extra_splits)` stops and resumes at each named recorded
  step, exported as `census_trait_gradient_split_tf24(scm, splits)`. The reverse
  pass is a backward linear recursion over steps, so composition is associative and
  a split must change nothing: tolerance is exactly zero and no reference is needed,
  because it is a property the implementation either has or does not.

  What it watches is anything carried across a step boundary that is not the
  adjoint. The trait accumulator accumulates by design, but the block workspace, the
  tape, the knot adjoints and the strategy templates all live across steps, and a
  split forces a clean re-entry at the cut.

  Verified on the three-introduction stand -- 102 recorded steps, widening after
  steps 1, 24, 38 and 57 -- at an interior step, one step either side of a widening,
  and all three at once. Every one is `identical()` to the unsplit gradient.

  **Its non-vacuity guard is not decoration, and building it found the trap.**
  A split landing exactly *on* a widening is outside every segment's interior, so it
  cuts nothing and the equality holds between two identical sweeps.
  `SCM::adjoint_segments` counts the backward ranges the last gradient swept --
  12 unsplit (four segments times three metrics), 15 for one genuine cut, 21 for
  three -- and the check requires that count to rise. The boundary case is asserted
  rather than left as a trap.

* **A rejected step attempt is not a gradient question, and the ladder's framing of
  it was wrong.** The reverse pass never sees one: the adaptive run establishes the
  schedule, one state and step size are recorded per *accepted* step, and the
  replay steps by those sizes. So "require the gradient unchanged when rejections
  are excluded" tests nothing -- they are already excluded from everything the tape
  sees.

  What the design does have is a seam, which `SCM::store_trajectory` already names:
  a pinned replay does not reproduce the adaptive run, because a rejected attempt
  moves patch state that is not ODE state -- the first-same-as-last derivative
  carry, and anything cached on the patch. The check that matters is therefore
  *does the replay reproduce the run*, which is what the trajectory reference's own
  floor measures: **3.2e-12** on the four-node stand and **1.17e-12** on the
  three-introduction one. A rejected-attempt count is worth publishing only as the
  non-vacuity guard on that floor, since with no rejected attempts the equality is
  trivial.

* **Rung 5's fourth boundary probe is now buildable, and it passes.** Of the five
  single-channel probes, four test for an exact zero and were already in place; the
  fifth -- the newcomer's leaf area reaching the census through a field built
  without it -- has a **shortfall** as its dropped-channel signature rather than a
  zero, so it needed a reference and could not be written until the trajectory
  tangent existed. On the three-introduction fixture the replay reaches that
  stand's own census to 1.17e-12, and the field-borne columns agree at 4.35e-06
  (`1.k_I`), 4.98e-07 (`2.k_I`), 2.70e-04 (`1.a_l1`) and 2.21e-03 (`1.a_l2`).

  Two rung-5 items remain and both need a new entry point rather than test work:
  §8.2's marginal recruit needs something reporting a boundary node's
  establishment probability and net production, so a recruit can be asserted into
  the stiff band where the check is not vacuous; and §2.2's split-at-the-boundary
  check needs a sweep over a sub-range of a recorded trajectory, with the
  rejected-step count needing the solver's attempt count published.

* **The record-once-sweep-many permutation check cannot fail, and the injection
  above is what showed it.** With an active value deliberately held across the
  tape clear -- the exact defect that check exists for -- it passed, along with
  every other structural check in the file. Only the reference saw it.

  The cause is on the R side: `stand_gradient(scm, metrics = m)` computes
  **every** metric in C++ and subsets the result in R, so "swept alone" and "swept
  with the others" are the same sweep and comparing them is a tautology. Closing
  it needs `census_trait_gradient` to take a metric subset, so that a
  single-metric sweep is genuinely one recording; the check is left in place and
  marked, because the shape it watches for is real even though this instance of it
  is vacuous.

* **The census gradient was missing its direct term, and the trajectory reference
  is what found it.** A census reads the traits itself --
  `C = sum_k w_k n_k m(h_k, phi)` -- so its total derivative has a term at fixed
  state that no sweep produces. `SCM::census_trait_gradient` returned the
  trajectory term alone, so every metric's row was wrong for every parameter the
  metric algebra reads, and `a_l1` came back with **the wrong sign**: +0.1717
  against a true -0.2848 for leaf area.

  The suite could not see it. The check that would have -- suppressing the term
  and requiring the row to move -- was skipped for want of a switch, and the
  check that ran only established that an independently computed direct term was
  *the same size* as the reported row, not that it was in it. The gap
  `tangent - sweep` equalled that independent term on all four measurable columns
  to three significant figures.

  Fixed by `SCM::census_trait_direct<Metrics>()`, which records the census with
  the traits seeded and the state held and sweeps it once per metric -- so it is
  the metric algebra that is differentiated, and a metric added in `species.h`
  needs no edit. Adding it to the trajectory term is not double counting: the
  trajectory term is `(dC/dy)^T (dy/dphi)`, and the boundary node, which setting
  the state rebuilds through the field, is not in `y` at all.

  **This changes every census gradient `stand_gradient()` returns.** Worst
  affected on the four-node fixture: `a_l2` by a factor of 1.8, `a_l1` by a sign
  flip, `lma`'s above-ground-mass row by 13 percent. Columns the metrics do not
  read -- `k_I` among them -- are unchanged to round-off.

### Known issues

* **The scientific-surface drift guard is dark by default, and it is red.** It
  skips on CRAN, so an ordinary run has no live gate on a default that changes a
  model's output. Run with `NOT_CRAN` set it fails for all four models, on three
  changes, and only two of them are accounted for:
  * `control.node_density_in_birth_date`, `FALSE` -> `TRUE`, from the birth-date
    default. Intended.
  * `model_id`, `FF16@v1` -> `FF16@v2`. The bump was made; the snapshot was not
    re-blessed with it.
  * `control.GSS_tol_abs`, `0.001` -> `0.1`. **Unaccounted for**, and it is one of
    the four entries `gradient_control()` reports precisely because it changes the
    function being differentiated. Not blessed, for that reason.

  The third has a second problem beside its provenance: the quantity has **two
  defaults**. `Control` sets `1e-1`, `TF24_Strategy` declares `1e-3` as its own
  member, and the `Leaf` is constructed from the `Control` one -- while a comment
  in the same header reasons about "the ~GSS_tol_abs (1e-3) ceiling the leaf
  package documents", and `collar_census.h` warns that a loose value reports
  nearly every solve as pinned. So the looser value wins at the call site and two
  places in the tree assert the tighter one.

  What this does and does not put in question: every gradient measured on this
  branch was taken at `0.1`, and the interior operating point is evidently still
  well determined there -- the block's forward and reverse Jacobians agree to
  `2.4e-16` and the factorisation predicts out of sample to `2.6e-08`. That is
  consistent with the collar being obtained by solving its own first-order
  condition rather than by a golden-section search, which would make this tolerance
  largely vestigial. What it would bias is the *incidence* of points classified as
  pinned, which is a different claim and is not measured.

* **The reverse sweep's wall-clock cost is not established, and the one figure
  measured is bad enough to question the premise.** On one species at a lifetime of
  2 with the default node schedule — 81 cohorts, 171 accepted steps — the forward
  run takes **6.2 s** and one census gradient takes **1006 s**, a factor of 162.
  The sweep records and sweeps once per cohort per stage per step, so its cost is
  set by the node count, and 81 nodes is an ordinary schedule rather than a large
  one.
  * Taken at face value that inverts the project's cost argument: the adjoint
    exists because it answers 44 traits for one run plus one sweep, against 44
    re-runs, and 44 one-sided re-runs of this configuration would be about 273 s.
    A re-run difference is not actually an alternative — at production a relative
    step of `2e-7` in leaf mass per area moves a mature stand between alive and
    identically zero — so this is a statement about cost and not about method. But
    the cost premise is not currently met at this configuration.
  * **Nothing in the testing ladder measures this.** The flatness guarantees the
    design rests on are about peak memory and about the differentiation-target
    count, and rung 3 does check them (a recording of 62,896 entries at one cohort
    and at four). Time is the quantity that is not flat and no rung asserts
    anything about it, which is why a test file absorbed 30 to 40 minutes without
    anyone deciding it should.
  * One figure at one configuration, measured once. What it justifies is a cost
    assertion in the ladder, not a conclusion.

* **The introduction boundary's transpose disagrees with the trajectory tangent,
  by up to 3.4e-03 on the columns a widening carries.** This is what is left after
  the boundary node's row in both reductions was closed (below): that fix moved the
  rate-level residual by 200x and left this one untouched, so they are two defects
  and not one.

  The attribution is a prediction rather than a fit. On a stand whose trajectory
  never widens -- one cohort per species, both introduced at t = 0 -- `a_l2`'s
  residual falls from 2.9e-03 to 6.1e-06, `a_l1`'s from 6.9e-04 to 1.2e-05 and
  `k_I`'s from 1.2e-04 to 1.2e-05. So the disagreement is carried by the widening
  events, which means `Patch::introduction_adjoint` against the tangent's own
  `introduce_new_nodes`, and it is rung 5's business rather than rung 4's.

  Every column whose route does not run through a widening agrees to 1.2e-04 or
  better, and the shortlist's tightest are at 5.7e-13.

* **A dense TF24 stochastic run throws at the default ODE step cap** (#599). The
  soil water balance is stiff — the conductivity curve's exponent is
  `2*n_psi+3 ≈ 16`, so `K_sat/dz` is ~543/yr at the defaults — and at
  `Control()`'s `ode_step_size_max` of 5 yr the explicit RKCK stepper diverges
  rather than losing accuracy. A failing run reaches
  `theta = [-51.4, 52.0, 0.146, nan, nan]`, and once the rates are non-finite the
  step is not even rejected, because `adjust_step_size` derives its error ratio
  from a NaN. Downstream, the leaf's collar root-find is handed a soil potential
  no retention curve can produce and throws
  (`find_root_psi ... do not bracket the root`). Measured over seeds 1–40 at
  `patch_area = 1` (~105 individuals m⁻²): 17 of 40 runs throw on the code before
  this release's stochastic-solver work, 5 of 40 after it. Setting
  `ctrl$ode_step_size_max <- 0.05` gives 0 of 60; the relationship is not monotone
  (0.2 fails *more* often than 5 does), so that value is measured rather than
  derived. It is left as a user setting rather than clamped in the library,
  because the right fix is a per-environment bound or odelia's stiff RODAS stepper
  for environments carrying stiff state, not a global default that every model
  pays for. `test-stochastic-patch-runner.R` sets it for its own TF24 runs so the
  suite is deterministic; that is a test setting and changes nothing for callers.

* On the birth-date path, the *worst-case* node of the reported height density
  does not improve with schedule refinement, even though the typical node does
  (see the reporting note under New features). `|dh/dτ|` is a ratio of two
  differences that both shrink as the schedule is refined, so cancellation error
  sets a floor, and refinement adds nodes in the near-empty tail where that floor
  is worst. Aggregate and plotting use is sound; individual node densities far
  out in the tail are not. Note the exact `−g(H₀)` Jacobian does not help here:
  it is exact only at a node's own birth, and the worst cases are interior nodes
  in a compressed region of the size distribution (log density ~ −11), not the
  boundary. Removing this properly means evolving `∂h/∂τ` along the
  characteristic, whose rate is `∂g/∂h · ∂h/∂τ` — which is the compression term
  this change exists to avoid, so it would have to be an opt-in extra state
  wanted only for reporting.

### Minor changes & bug fixes

* **The inflow boundary's own adjoint reaches the trait accumulator.** The boundary
  node holds no ODE row, so the census's sensitivity to its density was accumulated
  into `Patch::boundary_node_adjoint` and read by nothing: written in three places
  and consumed in none. The adjoint at the boundary times the boundary condition's
  own derivative is the whole of the boundary's contribution, and the second factor
  was never applied.

  `Patch::boundary_condition_adjoint` now records the condition at an active scalar
  and takes one vector-Jacobian product, which delivers the state rows and the trait
  rows together rather than hand-writing either. Measured against a forward tangent
  of the same trajectory, on a stand carrying one cohort per species:
  `recruitment_decay` falls from **1.9e-02 to 2.9e-12** and `a_d0` from **2.0e-03
  to 2.1e-05**. Both reach the census through establishment and through nothing
  else, so their columns were the whole of the missing channel; neither had ever
  been compared against a reference. `k_I` also improves, to 1.4e-07.

  The condition is recorded in the density rather than in its logarithm, and the
  incoming adjoint divided by the density to match. At a marginal recruit the log
  density and its trait sensitivity both diverge while the density and its own
  sensitivity both tend to zero, so the logarithmic form is a vanishing adjoint
  times a diverging row.

* **A reloaded state now carries the boundary node the run carries.** A stage
  evaluates the inflow condition twice, in two different fields: once with every
  species' boundary interval left off, and again in the field rebuilt to include it.
  The reductions were built on the first; the water aggregation, an introduced node
  and a census all read the second. Every rebuild the sweep does started from a
  recorded state and stopped at the first evaluation, so the sweep linearised a
  boundary node the trajectory never carried — by **1.2 per cent** in the boundary
  density by the end of a two-year run, guarded only by a length check.

  `Patch::set_recorded_state` is the named form of the reload, and the census seed,
  the census's direct term, the replay of an introduction and the introduction's own
  transpose all take it. Within one stage the two transposes are now each linearised
  at the evaluation their own forward pass read. No forward number moves: the FF16
  bit-identity guard and its reference comparison both pass unchanged.

* **The census's direct term is reported on its own and refereed for the first
  time.** `census_trait_direct_tf24()` returns the metric's own reading of the
  traits at fixed state — the one route to the census no transpose check touches,
  and the one report 08 §6 asks for a switch on. `census_trait_difference_tf24()`
  is the same quantity differenced in plain double with the strategy moved in
  place, so the two share no path: one records the census and sweeps a tape, the
  other evaluates the census twice. They agree to **5.1e-09** over the matrix,
  which is the difference's own round-off rather than a disagreement. The switch
  test that had been skipping now runs.

  The difference perturbs the *prepared* strategy rather than rebuilding from
  `Parameters`, because a rebuild re-runs preparation and so carries the birth-size
  channel the differentiated path imposes to zero — which is exactly the channel
  the allometric constants would be refereed on.

* **An introduction's transpose derives which widened row each narrow row became,
  instead of writing it out.** `introduction_adjoint` built a per-species index map,
  copied the rows an introduction does not touch, and contracted only the
  newcomer's — which asserts that the first group is an exact identity in the state
  and carries no trait row. The recording now emits the whole widened state and the
  contraction runs over all of it, so the map is derived and the assertion is gone.
  It is the shape report 01 §5 asks for, where a narrowing written as a truncation
  of the tail is wrong for every species but the last.

  **No number moves**: every figure in the item below is bit-identical across the
  change, which is also what establishes that those rows were the identity.

* **A third reference now exists, and it says the sweep and the tangent are both
  wrong together on the traits that reach birth size — by far more than the stated
  domain records.** A whole-run central difference of the census shares no path with
  either: it runs the forward model twice. It is step-stable across `1e-05`,
  `1e-04` and `1e-03` at two- and four-year lifetimes, so it is usable at these
  sizes despite report 06 §8's warning about production stands.

  Against it, on a one-species two-cohort stand, leaf area:

  | lifetime | `a_l1` sweep/truth | `lma` sweep/truth |
  |---|---|---|
  | 0.4 yr | 2.184 | 1.507 |
  | 1 yr | 1.419 | 1.869 |
  | 2 yr | 1.059 | 2.368 |
  | 4 yr | **1.003** | **3.054** |

  Two opposite readings, and the second is the serious one. `a_l1` **converges**: its
  factor of two is a short-run artefact, where a census with almost no accumulated
  biomass is dominated by the seed and the boundary node, and it is gone by four
  years. `lma` **diverges**: its gradient is 2.4 times the truth at two years and 3.1
  at four, step-stable at every step size, and growing with run length. Report 06
  §11 records the birth-size channel as "worth about 3 per cent for leaf mass per
  area"; measured here it is 137 per cent at two years and 205 at four.

  The sweep and the tangent agree with each other throughout — `a_l1` at four years
  differs between them in the fourth digit — so this is not the sweep-versus-tangent
  residual below. It is a channel both impose to zero, and on the most-measured trait
  in the model it does not decay.

* **The two allometric constants are measured, and they are the "size unknown" entry
  in the stated domain.** Report 06 §11 lists `a_l1` and `a_l2` as *short,
  unmeasured — the reduction's contribution is dropped*, the same class as the
  extinction coefficient's documented 3.041 per cent shortfall, with the size never
  reported. It is **3.3e-02** at worst, on stem area, once a species carries more
  than one cohort. Report 05 §6.1 predicts the shape: the pair "also has a live path
  through the cohort step and through the size-space adjoint, so the reduction's
  contribution is an addition to a non-zero row rather than the whole of it" — a
  finite, correctly-signed, short row.

* **Known, and not closed: the two allometric constants disagree with a tangent by
  up to 3.4e-03 once a species carries more than one cohort.** It is per species and
  it switches on with that species' second cohort — with one cohort each `a_l2`
  measures 6.0e-05, and introducing a second cohort of one species moves that
  species to 3.4e-03 while the other stays at 1.6e-06. Further cohorts do not add to
  it.

  Five causes are now excluded by measurement rather than by argument: the inflow
  condition and the reload above (closing those moved `recruitment_decay` by ten
  orders and left this at 3.05e-03 against 3.06e-03), the census's direct term
  (exact, above), the reduction's own parameter rows (rung 3 forms them entry by
  entry on a two-cohort fixture and takes `a_l1`/`a_l2` by name), and the rows an
  introduction carries through untouched. **The whole of it is in the trajectory
  term.**

  **The introduction boundary is excluded, and that is now three checks rather than
  an argument.** The introduction map's whole Jacobian — the pre-introduction state
  and the traits in, the widened state out — agrees forward against reverse to
  **7.5e-16 over every one of its 3729 cells**, at a widening into an empty patch
  and at one into a populated patch. The rows an introduction carries through are
  bit-identical to the recorded state and the rows it writes are bit-identical to
  the boundary node the patch held, both at tolerance zero. The newcomer's own
  sensitivity to these two traits is 4.0e-07 and 1.9e-05, the seed height being
  imposed passive.

  So what is left is a species carrying **two cohorts rather than one**, with the
  introduction that necessarily accompanies it excluded. Rung 3 forms the per-stage
  Jacobian entry by entry on a two-cohort patch and finds it correct, so the
  remaining difference is the **time dimension**.

  It has a per-step character without being a clean lost term. Tightening the ODE
  tolerance from `1e-04` to `1e-10` takes the run from 98 steps to 635 and the
  residual from `3.3e-02` to `1.3e-02` — six and a half times the steps for two and
  a half times the residual, near `h^0.6`. The sweep and the tangent replay the
  **same** recorded step sizes, so a truncation error would cancel between them
  exactly; what does not cancel is a linearisation taken at a point they do not
  share.

  **It accumulates per step while a species carries two cohorts.** Moving the second
  introduction later shortens the time two cohorts coexist while widening the
  quadrature's interior interval, and the residual follows the former: `7.6e-03`
  over a duration of 0.38, `4.0e-03` over 0.20, `4.7e-04` over 0.02 — about 0.02 per
  unit time throughout. Three cohorts land where two of the same duration would, so
  it turns on at the second cohort and does not scale with the count. That excludes
  the quadrature's interior interval, which is the first thing a second cohort
  creates and which moves the opposite way.

  **The stage state and aux restore is exact**, which was the one per-step path rung
  3 does not exercise. Instrumenting that boundary — recomputing the rates after the
  restore and comparing against the stage the forward pass computed — gives a worst
  relative difference of **zero** over every stage of every step, on one cohort and
  on two.

  **And the light field's discretisation is excluded, which answers a question
  report 03 §3.3 left open.** That report asks whether the dropped knot-position
  channel shrinks with knot density, calling the convergence "the falsifier for the
  passive-position treatment", and records that it had never been checked. Built at
  33, 65 and 129 knots: `k_I`'s residual falls from `2.22e-10` to `4.58e-11`, so the
  treatment **does** converge — the falsifier comes back in its favour. `a_l2`'s is
  `1.4371e-03` at all three, identical to five figures while the forward model
  itself moves, so it is not a field quantity at all.

  **The stage recursion is excluded, and by injection rather than by reading.**
  Report 08 §7 lists a lost stage term as the one defect with no signature of its
  own, and names the fault: replace the sweep's scatter over every earlier stage
  with the immediate predecessor alone. Built into `odelia` and measured against the
  probe pair, it moves the **one-cohort control from 2.2e-07 to 7.0e-03** — a margin
  of thirty-one thousand — and saturates the longer fixtures at 1.0. So the check is
  live and hugely sensitive to that fault, and the real defect does not have its
  shape: a lost term raises the one-cohort control to the same order as the
  two-cohort case, while the defect leaves it at 2.2e-07. `odelia` is unchanged; the
  injected build was measured and reverted, and the baseline returns bit-identical.

  Two measurement corrections came with this and both had been hiding size. A
  column must be normalised **per metric row**: leaf area and above-ground mass are
  of order one while stem area is of order 1e-04, so a stem-area row wrong by three
  per cent read as 1e-06 of the column's peak. And on a stand whose species carries
  two cohorts the total is a **near-cancellation** of the direct and trajectory
  terms — the ratio reaches −16.8 on stem area — so the column's relative error is
  an amplified view of the trajectory term's, and both are now reported.

  On the birth-date coordinate a second cohort cannot exist without an introduction
  after `t = 0`, so per-introduction and per-interior-trapezium causes are not
  separable by scheduling alone. Seeding the cohorts does not separate them either:
  a run resumed from a populated state lands far outside the ladder's declared
  regime (relative reserve 0.44 to 0.73 against 0.02 to 0.30, gate slope 0.010
  against 0.4), which invalidates it as an instrument. **The resumed path is
  unrefereed and its disagreements are larger again — worth its own item.**
  `ladder_introduction_residual()` carries the measurements and the exclusions.

* **The reverse sweep's soil and trait rows are refereed for the first time, and
  they are correct.** Three references now run, each scoped to where it is valid,
  and the selection rule is the supplied row: difference where nothing is grafted,
  tangent where something is.
  * The forward tangent is exact but structurally silent about the water channel —
    the environment holds its integrated state as a `double` store and takes it
    passively, so every soil column and soil rate row of its Jacobian is exactly
    zero. That silence is now asserted rather than left to a residual, because a
    reference silent about a channel and a transpose wrong in it agree perfectly.
  * A plain-double central difference of the same right-hand side referees the
    channel the tangent cannot reach. The sweep agrees to **3.0e-06** on the state
    contraction and **4.3e-05** over the 70 trait columns such a difference can
    reach, both against the reference's own measured error rather than a literal
    tolerance. So `dsoil_K_dtheta` and `d(psi)/d(theta)` — hand-written, with
    nothing structural enforcing the pairing — are right.
  * The remaining 18 columns are the leaf's own nine traits per species. A
    perturbation of a *prepared* strategy cannot reach them: the leaf holds its own
    copy from preparation and a rate evaluation never pushes one back, so such a
    difference reads exactly zero whether the row is right, wrong or absent.
    Rebuilding the strategy from its parameters does reach them, and pays no
    birth-size penalty on exactly these columns because the seed height solves
    `mass_live(h) = omega`, which reads no leaf trait. The two schemes are
    unbiased on disjoint sets and together referee the whole row: worst
    **1.3e-03**, the leaf's own solve tolerance appearing on both sides.
  * One consequence worth knowing: because nothing on the forward path pushes a
    leaf trait into the leaf, setting one on a live strategy reaches no equation.
    The reported row is the derivative of setting it *and* re-preparing.

* **Re-blessed the references the birth-date default moved, which had been left
  stale.** Carrying the size distribution in birth date by default changed forward
  numbers -- the commit that made the change said so -- but did not update the
  references that encode them, so `test-strategy-ff16.R` and `test-patch.R` were
  failing against the height coordinate's answers. FF16 offspring production is
  `17.1720` against `16.8846`, its accepted-step count 276 against 307, and the
  two-species pair `c(12.04841, 16.59391)` against `c(11.99578, 16.47192)`.
  * Each re-blessed value now **names the coordinate it belongs to**, and the
    whole-run test asserts that the other coordinate gives a different answer
    rather than pinning a second set of literals. The two coordinates are
    different functions, not two discretisations of one, so a reference that does
    not say which one it is for cannot be interpreted.
  * `test-patch.R`'s two FF16 vectors are asserted alongside the relations that
    produce them: the density rate is exactly minus the mortality rate, which is
    the birth-date equation itself, and at a birth rate of one the log density is
    exactly minus the cumulative mortality, because the boundary condition no
    longer divides by the growth rate. Those two relations are what the stale
    literals `-0.78726` and `1.08695` were hiding -- they differ from the new
    values by `-dg/dh` and `-log(g)` exactly.
  * One of those failures was a latent test bug rather than a stale number.
    `test-patch.R` read the per-layer uptake auxiliaries *before* any rate
    evaluation at the state it had just loaded, so it was reading zeros; it passed
    only because the height coordinate's density rate solves the physiology again
    at a displaced height and populated them on the way.
  * `FF16_generate_stand_report`'s test now skips when `ggridges` is absent. It is
    a `Suggests`, so a tree without it cannot answer that check rather than
    failing it.

* **The rank-two factorisation of the marginal profit is measured out of sample,
  and its recorded blocker was stale.** It is the load-bearing claim for the whole
  water channel: the argmax channel is two scalars times closed-form vectors, and
  a compensating pair fits every row of the potential family equally well, so a
  joint residual cannot detect an error in it.
  * The blocker on record was which direction the second term differentiates. That
    is settled, in the code and in the corpus: the second intermediate is uptake's
    sensitivity to the **collar**, so the second basis vector is
    `d(dE_up/dp)/dpsi_j`. And the pair is already solved from a **cross-family**
    pair — one soil potential and one root-carbon direction — which is what a
    prediction rather than a fit requires.
  * So every layer but the solved-from one is already a prediction, and refereeing
    it needs no new machinery. The supplied uptake rows use the factorisation; a
    difference of the plain-double block re-solves the collar and so carries the
    true `dp*/dpsi_j` with none of it. Out of sample over four layers:
    **2.6e-08**, against the difference's own error of 1.2e-08. In sample at the
    solved-from layer: 4.2e-07, set by the fit's own `1e-3` step.
  * A fixture note that fell out of it: the wettest layer's potential sits on the
    root vulnerability grid's lower bound (`0.0127952` against a minimum of
    `0.012795`), so a downward perturbation at the fit's own step is infeasible and
    the leaf refuses to bracket. The declared regime covers the potential ceiling
    and not this floor.

* **The uptake-by-potential block's rank-one claim is verified, and the check that
  said otherwise was wrong.** The block is a diagonal explicit part plus a rank-one
  argmax channel. Zeroing the observed diagonal removes `D + diag(u_i v_i)`, taking
  the rank-one term's own diagonal with it and leaving a matrix of full rank — so a
  rank test on that residue fails whether or not the claim holds. Tested as an
  out-of-sample prediction of the off-diagonal entries instead: worst **4.9e-16**
  over 13 predicted entries, and all 120 off-diagonal 2x2 minors vanish to 5.7e-16.

* **A `Patch` took its disturbance regime from the constructor's argument rather
  than from the member it had just validated**, so a patch built from a
  `Parameters` whose `max_patch_lifetime` was assigned *after* that object's own
  construction ran on the default 105.32-year regime instead of the assigned one.
  Both `Patch` constructors now read `parameters.disturbance`, which
  `parameters.validate()` derives from `patch_type` and `max_patch_lifetime` one
  line earlier.

  What was wrong is every patch built by `Patch::rebind_from`, which assigns the
  lifetime onto a default-constructed `Parameters` and so never re-derived the
  regime. That is on the census gradient's path twice -- the introduction
  boundary's recorded condition reads both the patch survival and the patch
  density, and the census recording rebinds -- and it is the patch the
  forward-mode tangent reference is taken on.

  **Scope of the forward claim, stated as narrowly as it was tested.** No FF16
  result changes, and that is measured rather than argued: `test-strategy-ff16.R`
  and `test-patch.R` fail identically with this change reverted, and the FF16
  bit-identity reference comparison passes with it in. The mechanism is that
  `scm_base_parameters("FF16")` uses the same `max_patch_lifetime` as the
  `Parameters()` default, 105.32, so the argument's regime and the member's are
  built from the same number. **Whether a forward run whose lifetime differs from
  that default is affected has not been established** -- it depends on whether the
  R boundary re-derives the regime on assignment, which was not measured. Every
  fixture in the gradient ladder is such a configuration.

  Measured on a fixture at `max_patch_lifetime = 2` and `t = 1.37`: the reference's
  offspring-production rate was **93.727756x** the model's, uniformly across every
  node of every species, which is exactly
  `pr_survival(1.37 | 105.32) / pr_survival(1.37 | 2)`. Every other rate agreed to
  round-off. After the fix the reference reproduces the model's rates to 0 and
  3.9e-16 on three fixtures, and the residual disagreement between the tangent and
  a plain-double difference of the same right-hand side falls from **0.986 to
  0.120** -- with everything that survives being a soil row or column, which is a
  separate and now-isolated gap (see Known issues).

* **The scenario gateway scores numerical viability and persistence as separate
  axes, and no longer leads with a match rate** (#572). `scenario_summary()`
  classified a run as a success on `finite && total > 0`; with `birth_rate = 1`
  (what `build_scenario()` sets) `offspring_production` *is* R0, so five of the
  eight hydraulic scenarios returned R0 between 2e-15 and 6e-14 — numerically
  extinct — and were all recorded as `persisted`. All eight "succeeded", none of
  the five expected failures failed, and only one replaced itself: the gateway
  returned no signal. It now reports **numerical viability** (`n_ran`,
  `n_crashed`, `viability_rate` — did the model run, which is what the CSV's
  "Model failure" means) and **ecological persistence** (`n_persists`,
  `persistence_rate`, at R0 >= 1) as separate, labelled axes, in the summary, in
  `scripts/run_scenario_gateway.R` and in the scorecard report. `n_match` /
  `match_rate` are still reported but demoted to *agreement with the CSV's
  crash-era expectations*: with the crashes fixed (#546/#552/#554) the match rate
  mostly measures how well those predictions have aged, not model quality — which
  is why fixing the model *lowered* it from 5/8 to 3/8. `status` and `outcome`
  are unchanged and the CSV is not rewritten.

* **The scenario gateway integrates in birth date, not height** (#590), via a
  new `scenario_control()` that supplies the `Control` every entry point in the
  framework now defaults to. The package default is untouched.

  Every scenario in the gateway is TF24, and TF24 is the model the two density
  coordinates genuinely disagree on: the compression term is the total
  derivative of growth along a cohort's trajectory, which equals `dg/dh` only
  when growth is a function of size, and TF24's reserve gate (#517) breaks that.
  So the gateway had been scoring the hydraulic model through a compression term
  that is wrong for it.

  Measured on the eight scenarios at `max_patch_lifetime = 100`, the coordinate
  change raises R0 on every one — S02 2.4x, S06 8.9x, S08 15x, S05 25x, S07 47x
  — and moves **S01 across R0 = 1** (2.59e-01 to 1.49e+00). Persistence goes
  **1/8 to 2/8**. Numerical viability stays 8/8 and CSV agreement stays 3/8.
  Birth date is also about twice as fast here (36 s against 70 s).

  Note what that implies about the old guard: `observed` did not change on a
  *single* scenario across a 47x swing in R0. It tests `finite && total > 0`,
  which every scenario has satisfied since the crash fixes landed, so it was
  pinned at 8/8 and could not move. The baseline diff in
  `test-scenario-gateway.R` therefore now diffs `persists` as well as
  `observed`, and **the baseline is re-blessed** under the new coordinate.

* **Retracted: "reserve-gated growth largely excludes the slower species".**
  `test-strategy-tf24.R` recorded that as a finding about #517's NSC reserve
  gate. It is a property of the density coordinate, and #590 flagged it as
  needing re-deriving. Re-derived here, at the test's own configuration
  (`lma` 0.0825 / 0.10, `max_patch_lifetime = 5`):

  | coordinate | fast | slow | ratio |
  |---|---|---|---|
  | height | 67.32 | 2.775e-04 | 2.4e5 |
  | birth date | 287.2 | 59.53 | 4.8 |

  So the two species **coexist at comparable abundance**; they are not
  separated by five orders of magnitude. At `max_patch_lifetime = 30` the ratio
  is 2.7 (2707 against 1004, matching #590) — i.e. it *narrows* with patch
  lifetime, where progressive exclusion would widen it.

  What survives the retraction is the weaker claim: the faster species still
  leads, by roughly 3-5x. Reserve-gated growth disadvantages the slower species
  without excluding it. How much of even that gap is attributable to #517 as
  against the trait difference itself is not measured here, and should not be
  read into these numbers.

  The evidence that this is a wrong derivative rather than an under-resolved
  one is refinement. Over two halvings of the node spacing (88 → 175 → 349
  nodes) the birth-date answers are already converged — fast
  287.2/287.1/287.2, slow 59.53/59.66/59.69 — while the height answers are
  still climbing (fast 67.32/73.18/74.98) and the exclusion ratio does not
  shrink toward the birth-date one, it *grows*: 2.43e5/2.49e5/2.51e5.
  Quadrature error closes under refinement; a different derivative does not.

  The height-coordinate assertions are kept, since they pin what the package
  default still does, and a birth-date case is added alongside so the
  coordinate that is correct for TF24 is actually covered by a test rather than
  only described in a comment.

* **`run_stochastic_collect()` now reports the environment.**
  `StochasticPatch::r_get_state()` had its environment leg commented out, so
  unlike `Patch::r_get_state()` it returned only `time` and `species`. The R
  collector then read `light_env`, a name nothing had ever produced, so every
  element of that field came back `NULL`. The field is now `env`, matching what
  the patch reports and what `run_scm()`'s collected output calls it, and it is
  populated. This was inert while the stochastic solver held its environment at
  the initial state; now that the environment is integrated, the soil trajectory
  is real and worth reporting. **Breaking:** the returned list field is renamed
  from `light_env` to `env`. Nothing could have depended on its contents, since
  it was always `NULL`, but code testing for the name will need updating.

* **The stochastic arrival schedule now scales with patch area.**
  `stochastic_schedule()` passed `patch_area` into `stochastic_arrival_times()`'s
  third positional argument, which is `delta_t`, leaving `patch_area` at its
  default of 1. Two things followed. The arrival rate never scaled with area, so
  the expected number of arrivals came out as `max_patch_lifetime × birth_rate`
  whatever the patch size, and a 50 m² patch was seeded like a 1 m² one. The
  binning interval was also silently set to the area, which for a large patch
  left only a handful of intervals over which a variable birth rate was
  averaged. Arrivals now scale linearly with `patch_area` as intended, and the
  interval keeps its 0.1 yr default. `run_stochastic_collect()` is the only
  caller; its runs change accordingly, and even at the default `patch_area = 1`
  they are unchanged in expectation but not bit-identical, because the finer
  binning draws from the RNG differently. The seeded baseline in
  `test-stochastic-patch-runner.R` was re-derived and its `patch_area` reduced
  from 50 to the default 1, which keeps the ~105-individual stand that test has
  always actually run rather than the ~5300 that 50 m² now implies.
* **An empty stochastic patch no longer discards the environment's integrated
  state.** `StochasticPatch::compute_environment()` calls
  `Environment::clear_environment()` when the patch holds no individuals, which
  is right for the competition profile — nothing is casting shade — but TF24's
  override also restored the soil states and cumulative-flux accumulators to the
  values the run began with. Now that the environment is part of the ODE system
  that discarded what the solver had integrated, on every derivatives evaluation
  while the patch was empty. `clear_environment()` is now the competition
  profile alone; restoring state moved to a new `clear_state()` that only
  `Environment::clear()` calls, so resetting for a new run is unchanged. This
  matters most for a schedule whose first arrival is some way into the run,
  which is `run_stochastic_collect()`'s default.
* **The stochastic solver now integrates the environment.** `StochasticPatch`'s
  ODE system was the species alone: `ode_size()`, `set_ode_state()`,
  `ode_state()` and `ode_rates()` did not chain through the environment as
  `Patch`'s always have, and `compute_rates()` never accumulated resource
  consumption, so `Environment::compute_rates()` was never called. An
  environment carrying ODE state was therefore held at its initial value for a
  whole run. TF24 carries nine such states — five soil-moisture layers and four
  cumulative fluxes — while FF16 and K93 carry none, which is why this went
  unnoticed. Stochastic TF24 runs change: soil water recharges from 0.214 to the
  drainage equilibrium 0.3106 and is then drawn down as leaf area grows,
  tracking the SCM's own trajectory on the same drivers to ~1e-4, and total leaf
  area at patch ages 1–3 moves from 11–16% below the SCM to within 1.4% of it.
  FF16 and K93 stochastic runs are bit-identical, given the guard described in
  the next entry. The consumption vector is sized by the environment's
  `n_resources()` rather than its ODE width, so TF24's four diagnostic flux slots
  are not mistaken for resources.
* **`StochasticPatchRunner::reset()` integrates to the first arrival with error
  control.** `reset()` advanced from time zero to the first scheduled arrival in a
  single `advance_fixed` step, which for `run_stochastic_collect()`'s default
  first arrival — uniform on (0, 50) — is one step of up to fifty years. That was
  harmless while nothing was integrated over it. Once the environment's own states
  are, it drives soil water out of range and nothing establishes: fixing the state
  handling alone produces nine failures and a runner that introduces no
  individuals. The leg now uses `advance_adaptive`, as every other advance in the
  runner already does — but only when there is state to integrate. An empty patch
  whose environment carries none (FF16, K93) keeps the fixed step, and that guard
  is what preserves bit-identity for those two. `advance_adaptive` walks the
  step-size controller up to `ode_step_size_max` even on a zero-width system, and
  `step_size_last` survives into the first real step, because only
  `SolverInternal::step()` writes it and `step_to()` — which `advance_fixed`
  drives — does not. Ungated, the first step after the first arrival therefore
  begins from a rejected five-year attempt rather than
  `ode_step_size_initial`, and FF16's collected trajectory moves by up to 1e-5
  relative: inside `ode_tol_rel = 1e-4`, so a valid realisation either way, but
  not identical.
* **A stochastic recruit's strategy-specific initial states are seeded from its
  birth environment.** `StochasticSpecies::introduce_new_node(environment)`
  computed the new individual's rates without first calling
  `Individual::set_initial_states()`, which the deterministic path calls from
  `Node::compute_initial_conditions()`. A TF24 seedling was born with an empty
  carbohydrate store rather than `a_st3` of its storage capacity. FF16 and K93
  do not override `set_initial_states()`, so they are unaffected.

* `SpeciesBase::control()` called `strategy->get_control()`, which does not
  exist on any strategy. The member had never been instantiated, so the error
  had never been compiled; it now reads `strategy->control`.
* **The TF24 rainfall driver is floored at zero.** Because drivers are
  interpolated with a cubic spline, an intermittent series undershoots below
  every supplied value, and negative rainfall gave negative infiltration and an
  unphysical drying rate. That failed two different ways depending on soil
  wetness: above residual moisture the water really was removed, while at or
  below residual the guard in `compute_rates` clamped the rate, so the removal
  was recorded in `sum_rainfall` but never applied and the water budget stopped
  closing. Drylands sit at residual for much of the year, so the second case is
  the common one. Note the floor bounds the *sign* only and is not a correction
  to the interpolation: the spline conserves the integral exactly (undershoot is
  compensated by overshoot), so discarding the negative lobes raises total
  rainfall by the undershoot area — ~7% for a realistic daily series. The remedy
  is not to spline an intermittent series; see `check_driver_interpolation()`.
  Runs with constant or smooth seasonal rainfall are unaffected (the scenario
  gateway's sinusoidal driver never goes negative), so the TF24 scientific
  version does not move.
* Corrected the comment on `Leaf::assim_colimited()`, which claimed "no dark
  respiration included at the moment" while the code subtracts `R_d_`. The
  function returns a *net* rate; the comment now says so, along with how to
  recover the gross rate.
* Water-budget test coverage extended from a single layer to 1, 5 and 15 layers,
  now including root uptake in the balance (the previous check ran for 0.01 yr,
  where uptake was negligible), across the saturated-to-dry range and at and
  below the residual-moisture floor. Closure holds to round-off (relative
  residual < 1e-12) in every case. A companion test documents *why* the residual
  clamp never leaks in practice: with `n_psi = 6.57` the conductivity exponent
  is 2·n_psi+3 ≈ 16, so K(θ) collapses and ψ(θ) diverges far above θ_r —
  transport has already stopped, making the clamp a safety net rather than an
  active mass sink.

* `run_scm()` now fails with an actionable message when the SCM size-density
  (characteristic) equations run away under extreme forcing (e.g. severe
  seasonal drought in TF24): a cohort density overflowing to `+Inf`, or the
  density-weighted resource uptake driving a soil-water state non-finite. The
  guard fires each ODE step, before the non-finite value propagates into the
  competition integral or physiology, replacing the previous opaque
  `Detected non-finite contribution` / `non-finite psi_soil` errors (#550).
* `Node` now records its introduction time and patch-age density at the moment
  it is introduced, so reproduction and integration-error calculations no longer
  re-derive these from the node schedule / disturbance regime after the run.
* `SCM::run()` accumulates the per-node refinement error when `collect_errors`
  is set, exposed as `combined_node_errors`; `SCM::refine_schedule()` runs the
  adaptive node-introduction loop in C++.
* The resource spline is now floored at 0, fixing spurious negative light
  values (notably the K93 light spline at high `k_I`) (#253, #497).
* `interpolate_to_heights()` no longer silently drops individuals in the largest
  size class on a coarse grid; it appends the actual largest individual per time
  step instead (#352, #497).
* Fixed `tidyselect` deprecation warnings from `integrate_over_size_distribution()`
  by using string column names in `across()`/`rename()` selection contexts
  (#375, #501).
* Account for variable patch size in light/competition calculations (#422).
* Fixed the argument order in `net_mass_production_dt()` (#389).
* Fixed integration of density (#345).
* Fixed `scm_support` to use the correct error name `net_reproduction_ratio_errors`
  and guard the zero-offspring case (#447).
* The ODE solver now continues when the minimum step size is reached (#413);
  refined the `interpolate_to_heights()` error check (#437).
* `run_stochastic_collect()` returned empty output: it read a `state` accessor
  that `StochasticPatchRunner` no longer exposed, so every step silently
  collected `NULL` (#498, #506). It now reads `StochasticPatch$state` (a new
  accessor mirroring `Patch$state`), and each species' per-step state matrix
  carries an `is_alive` attribute, so the height/survival trajectory is
  populated again.

### Internals & performance

* **The environment carries its integrated state at the scalar the model carries.**
  `Internals<S> vars`, `std::vector<S> resource_uptake`, `compute_rates` taking an
  uptake vector at `S`, and the soil retention and conductivity curves returning
  `S`. `Environment::set_ode_state` no longer strips to a value, and the whole
  ODE-state interface moved from the untemplated base onto `TF24_Environment<S>`,
  where the scalar is known; the base keeps the no-state defaults, which is what
  FF16 and K93 already relied on -- neither references `vars`.

  This is what makes a reference possible at all. A tangent whose soil state is
  held at a value differentiates a model with frozen soil: the state's tangent is
  dropped at every stage, so the soil rate carries none, and a cohort reading a
  layer potential sees no route from a trait to the water and back. Both
  directions were dead, and a reference silent about a channel agrees perfectly
  with a transpose that is wrong in it.

  **No forward result changes.** The FF16 bit-identity reference comparison, the
  FF16 and K93 suites, `test-patch.R`, `test-scm.R` (including its pinned
  step-size replay) and `test-environment.R` all pass unchanged, and production is
  `S = double` where every added conversion is the identity.

  Two smaller repairs fell out of it. `Patch::check_finite_ode_state` read
  `environment.vars` directly and now reads the ODE interface, so it is generic
  over environments rather than reaching into one. And `TF24_Environment`'s
  potential cache is now invalidated when the state is written, not only when the
  state's *value* changes -- a cache keyed on an exact value cannot see a changed
  derivative behind an unchanged one, which is a hazard that only exists once the
  state carries a derivative.


* The scenario gateway's seasonal rainfall driver places its spline knots **per
  year** (48 yr⁻¹) rather than per run, so the realised seasonality no longer
  depends on `max_patch_lifetime`. The previous `max(200, mpl × 6)` gave 6 knots
  per annual cycle at the default `mpl = 100`. Measured, that was more accurate
  than it looked — cubic interpolation of a sine at 6 points/cycle is accurate to
  5.4e-3 on a peak of 3.0 (0.18%) and conserves the annual total to 1e-5 % — so
  this is robustness, not a bug fix; the error is 4th-order in the spacing and
  would degrade at larger `mpl`. Gateway effect: offspring production moves by up
  to 6.2% relative (demography amplifies the 0.18% driver change), with every
  success/failure classification unchanged, so the baseline needs no re-blessing.
* Each strategy now keeps its biological parameters in a value-member struct
  (`FF16_Pars`/`K93_Pars`/`TF24_Pars`) exposed to R as a nested `pars` list;
  derived/computed members (eta_c, height_0, the TF24 Leaf model, solver
  tolerances) stay plain. Hot-path access stays inlined and FF16 output is
  bit-identical (#410).
* `FF16_expand_state`/`TF24_expand_state` no longer duplicate the C++ allometry
  formulas; they call the strategy's own functions via
  `FF16/TF24_strategy_expand_allometry` (`src/strategy_expand.cpp`). Output
  columns are unchanged; `K93_expand_state` remains a no-op. Fixed a latent bug
  in the (previously dead) `mass_above_ground` (it summed root and omitted
  heartwood); it now returns leaf+bark+sapwood+heartwood, matching the rate
  version and the expand_state output (#254).
* `Individual` no longer knows any strategy's state/aux layout: it passes the
  whole `Internals` to `strategy->compute_competition(z, vars)` /
  `net_mass_production_dt(env, vars)`, and each strategy reads its own slots
  (#266, #254). A new strategy with a different layout can reuse the generic
  `Individual` unchanged.
* Hot-path optimisation of FF16/TF24 (~2.7–3.5× FF16, ~1.3× TF24) (#471), K93
  SCM perf (#493), and shared `pow(area_leaf, a_l2)` in FF16 `compute_rates`
  (#361, #494). See the `profile-plant` skill for the methodology.
* `Species` and `StochasticSpecies` now share storage, the ODE state plumbing,
  and per-element serialisation through a compile-time `SpeciesBase`; the
  stochastic species stores a `StochasticNode` element (Individual + `alive`) in
  place of a parallel `is_alive` vector, and stochastic deaths are applied to the
  solver's live system in place (#217, #506). Internal only — FF16 is
  bit-identical and shows no timing change.
* Earlier internals (since 2.0.0): default compilation set to `-O2` (#365);
  assimilation decoupled into a per-strategy `Assimilator` (#313); extrinsic
  drivers refactored into an `ExtrinsicDrivers` class (#334, #340); test-only
  reference plant relocated into the tests (#418); unused logging removed (#436);
  Makefile dev-loop speedups (#319); general `R CMD check` cleanup (#440).

### Documentation & tooling

* **CI is faster and the test suite is leaner.** The suite ran 147 s serially,
  with one file (`test-tf24-arid-corner.R`) at 24% of it, and tests were ~65% of
  each `R CMD check` job. Now 124 s (-16%), and no single file is over 19%.
  * `test-tf24-arid-corner.R` 35.5 s -> 12.7 s. It computed the same dry-start
    TF24 run four times; the completed run is now memoised and shared (the
    `short_run()` pattern from `test-density-coordinate.R`). #571's nine-point
    soil-moisture sweep is cut to the three values bracketing the
    residual-moisture floor, the full set being recorded in #571.
  * Dropped `"a failed light spline reports the patch state that caused it"`. Since
    #574 fixed the cause its `skip_if()` fired every run, so it paid for a full
    SCM run and asserted nothing. See the comment in the file for how to bring it
    back as an unconditional test.
  * Dropped `"TF24 water budget closes independently of layer count"`: it re-ran
    the same three configurations as the test above it to re-assert a bound that
    test already checks per layer.
  * Benchmarks no longer run on every push and PR — `run_plant_benchmarks()` was
    costing 3-5 min per trigger on a macOS runner while nothing compared its
    output against a baseline. Now weekly plus `workflow_dispatch`, on Linux.
  * Both workflows cancel superseded in-flight runs on the same ref, and
    `R-CMD-check` gained a `ubuntu-latest` job — commented out since CI was first
    set up (#295), so plant had never had a green Linux run.
* Narrative docs (former `vignettes/`, theory, dated posts) migrated to the
  [Overstorey](https://traitecoevo.github.io/overstorey/) site; the pkgdown
  site is now the function reference only (#496).
* The strategy-scaffolder workflow and a profiling workflow are now captured
  as the `plant-new-strategy` and `profile-plant` skills (#495, #492).
* Relicensed from GPL-2 to **AGPL-3** (#457).
* Upgraded the minimum C++ standard from C++14 — currently C++20 (#442).
* Expanded the K93 (#421) and self-thinning (#369) vignettes; added a draft
  `extrinsic_drivers` vignette (#340).

## Plant 2.0.0 release notes

v2.0.0 was released on 25/02/2021

### Major changes

* Improved templating of strategies and environments to allow for inheritance and re-use.
  See: `FF16r_strategy` for example of method overloading and `K93_environment` for environment inheritance.
* Added two new models: Kohyama 1993 (K93) and a soil water strategy
* Recovered the FF16r strategy
* Decoupled patch, environments and strategies by moving several routines to strategies, e.g. assimilation
* Hyperparameterisation now handled in R only.
* Added a strategy implementation vignette

### Minor changes

* Moved strategy defaults to header
* Moved strategy and environment specific files to `inst/include/plant/models`
* Renamed several functions, including:
  * `germination` -> `establishment`
  * `Plant` -> `Individual`
  * `area_leaf` -> `competition_effect`
  * `area_leaf_above` -> `compute_competition`
  * `vars_phys` -> `rates`
  * `scm_vars` -> `compute rates`
* Updated vignettes, documentation and tests
* Switch Ubuntu versions in Travis pipeline
* Added `r_init_interpolators` methods for environment initialisation
* A few fixes to scaffolder and tests
* Increased lenience on integration test
* Fixed compiler warnings
* Removed PlantPlus

A full account of changes from the previous version is available on GitHub: [v1.2.1...v2.0.0](https://github.com/traitecoevo/plant/compare/v1.2.1...v2.0.0)

## Plant 1.2.1  Release Notes

v1.2.1 was released on 20/09/2019

### Major Changes

### Minor Changes

- update Makefile to use `pkgbuild` instead of `devtools` for building dll ( because of upstream changes )
- switching to using `remotes` instead of `devtools` for installing from github ( because of upstream changes )

A full account of changes from the previous version is available on Github: [v1.2.0...v1.2.1](https://github.com/traitecoevo/plant/compare/v1.2.0...v1.2.1)

## Plant 1.2.0  Release Notes

v1.2.0 was released on 20/03/2018

### Major Changes

- New strategy scaffolder ensures higher level operations work across all strategies and for new strategies. This includes `lcp_whole_plant` `XXX_PlantPlus`, `grow_plant_to_size`, `grow_plant_to_height`, `grow_plant_to_time`
- website now builds via `pkgdown` package (used to use staticdocs but this is deprecated)
- simplified workflow for building website
- converted supporting materials from tex to Rmd

### Minor Changes

- start documenting notes for developers in `inst/docs/developer_notes.Rmd`
- add `CITATION` file
- Address many issues in documentation and package setup causing rcmdcheck to fail

A full account of changes from the previous version is available on Github: [v1.1.0...v1.2.0](https://github.com/traitecoevo/plant/compare/v1.1.0...v1.2.0)

## Plant 1.1.0 Release Notes

v1.1.0 was released on 2/02/2018

### Major Changes

- Now compiles and runs on Windows machines (requires R 3.3.0 or newer)
- Further details on installation
- Enable assembly_parameters to accept more named arguments

### Minor Changes

- Add Appveyor for build tests on Windows machines
- Update tests to use latest version of testthat
- Remove package traitecoevo/callr, previously used to make system calls  
- Makefile: Add roxygen and RcppR6 to compile target
- roxygen & Rcpp updates
- Added a `NEWS.md` file to track changes to the package.

A full account of changes from the previous version is available on Github: [v1.0.0...v1.1.0](https://github.com/traitecoevo/plant/compare/v1.0.0...v1.1.0)

## Plant 1.0.0 Release Notes

v1.0.0 was released on 23/02/2016

This version corresponds to the paper describing the package:

Falster, DS, RG FitzJohn, Å Brännström, U Dieckmann, M Westoby (2016) plant: A package for modelling forest trait ecology and evolution. Methods in Ecology and Evolution 7: 136-146, doi: [10.1111/2041-210X.12525](http://doi.org/10.1111/2041-210X.12525)

A full account of changes from the previous version is available on Github: [v0.2.2...v1.0.0](https://github.com/traitecoevo/plant/compare/v0.2.2...v1.0.0)


## Plant 0.2.2 Release Notes

v0.2.2 was released on 1/06/2015

Draft paper about package submitted to Methods in Ecology & Evolution.

### Major changes

First stable release of advanced implementation of plant
