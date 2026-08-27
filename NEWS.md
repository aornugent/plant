## Plant (development version)

While the reverse-mode gradient work is in progress this file carries **breaking
changes only**, and they are kept current because the `plant-update-interface`
skill reads them as its spec. Everything else that was here -- what has been
measured, what each defect excludes, and what is open -- was a development record
rather than release notes, and is in the git history rather than the tree.

### Breaking changes

These change the R-facing interface and require updating downstream code. Each
entry gives the `old -> new` migration; the `plant-update-interface` skill
(`.claude/skills/plant-update-interface/`) reads this section to migrate
products using plant.

* **A finite-difference arm that crosses the feasible boundary now refuses.**
  Requires the matching phylloptim. The leaf's environment rows are taken at a
  frozen collar, and the entry point they went through clamped that collar into
  whatever interval the PERTURBED state had -- so at a pinned operating point,
  where the collar sits about a millionth of the interval's width off the wet
  bound and a 1e-3 trait step moves that bound two orders further, one arm
  answered about a different collar. No migration: the affected metrics were
  returning a number and now report `refused` with the collar named. Measured at
  psi_soil 5.0, `dprofit/droot_b` came back as 184.699 where the same difference
  well inside the interval reads 0.0056.

* **The dry pin is reported by its arm, and the inverted interval by its own
  name.** Requires the matching phylloptim. Migration for anyone matching on the
  kind a refusal message names:
  * `"pinned-dry"`         -> `"pinned-dry-root-crit"` or `"pinned-dry-root-psi-crit"`
  * `"hydraulic-shutdown"` -> unchanged, except the inverted-interval exit, which
    is now `"infeasible-bracket"`

  The two dry arms are different functions of the inputs -- the stem's continuity
  root is a search result, the root's own critical potential is a registered
  constant -- so a consumer forming the bound's row needs to know which bound.
  Measured: at shipped defaults every dry pin is on the root-crit arm.

* **The census trait gradient says why a metric has no numbers, instead of
  raising.** A refusal used to escape to the R prompt as an error, with nothing a
  caller could inspect and no way to keep the metrics that did answer.
  Migration:
  * `census_trait_gradient_tf24(scm, m)` (a list of numeric rows)
    -> `census_trait_gradient_tf24(scm, m)$gradient` for the same rows
  * `census_trait_gradient_split_tf24(scm, splits)`
    -> `census_trait_gradient_split_tf24(scm, splits)$gradient`
  * `stand_gradient(scm)` gains `$refusal`; `$gradient` is unchanged in shape and
    meaning

  Both C++ entry points now return `list(gradient, refusal)`. `refusal` carries
  one entry per metric -- the reason, and the species it was found on -- or `NULL`
  where the metric answered.

  **A refused metric's whole gradient row is `NaN`.** Refusal is metric-level: a
  sum has no defined value with an undefined term, so no localisation within a
  metric is available. Metrics are independent of one another.

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


### Added

* **The gradient's incidence counters.** The leaf classifies its operating point
  by the branch taken and the next plant overwrites it, and a clamp that severs a
  row leaves a number indistinguishable from a true zero -- so neither was
  recoverable after a run. Additions only:
  * `census_operating_point_counts_tf24(scm)` -> per-species counts by kind
  * `census_operating_point_names_tf24()`     -> the kinds, in that order
  * `census_clamp_counts_tf24(scm)`           -> per-species counts by clamp site
  * `census_clamp_names_tf24()`               -> the sites, in that order
  * `census_clear_operating_point_counts_tf24(scm)` -> reset both, per run

  Read off the live system, so take them from the object you ran. Measured on a
  drought run that refuses: 99.71% interior, 0.29% pinned-dry -- and that 0.29%
  is what makes every metric's gradient undefined.

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
