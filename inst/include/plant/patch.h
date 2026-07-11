// -*-c++-*-
#ifndef PLANT_PLANT_PATCH_H_
#define PLANT_PLANT_PATCH_H_

#include <plant/parameters.h>
#include <plant/species.h>
#include <plant/util.h>
#include <plant/ad_value.h>
#include <odelia/ode_interface.hpp>

#include <plant/disturbance_regime.h>

#include <algorithm>
#include <limits>
#include <type_traits>

using namespace Rcpp;

namespace plant {

template <typename T, typename E> class Patch;

namespace detail {
// A strategy is AD-liftable when it provides its own scalar rebind. Detected by
// SFINAE so the double-only strategies (TF24, K93) are cleanly excluded.
template <class T, class = void>
struct strategy_is_ad : std::false_type {};
template <class T>
struct strategy_is_ad<T, std::void_t<typename T::template rebind<double>>>
    : std::true_type {};

// Carries the Patch `rebind` alias only for AD-liftable strategies, so odelia's
// rebind_or_self / has_rebind_from probes fail by SFINAE (member absent) for the
// double-only strategies instead of hard-erroring inside the alias expansion.
template <class T, class E, bool = strategy_is_ad<T>::value>
struct PatchRebind {};
template <class T, class E>
struct PatchRebind<T, E, true> {
  template <class S2>
  using rebind = Patch<typename T::template rebind<S2>,
                       typename E::template rebind<S2>>;
};
}

template <typename T, typename E>
class Patch : public detail::PatchRebind<T, E> {
public:
  using value_type = typename T::value_type;

  typedef T                 strategy_type;
  typedef E                 environment_type;
  typedef Individual<T,E>   individual_type;
  typedef Node<T,E>         node_type;
  typedef Species<T,E>      species_type;
  typedef Parameters<T,E>   parameters_type;

  Patch(parameters_type p, environment_type e, plant::Control c);
  void reset();
  size_t size() const {return species.size();}

  // odelia differentiable-System contract -----------------------------------
  // The active scalar the cohort physiology carries is value_type (above). The
  // environment stays double -- an invasion reads the competitive landscape
  // frozen -- so the rebind keeps E and lifts only the strategy. The `rebind`
  // alias itself is inherited from PatchRebind, present only for AD-liftable
  // strategies (odelia detects it by SFINAE; AD is opt-in).
  using PatchRebind_ = detail::PatchRebind<T, E>;

  // Config-only copy of this patch onto the scalar S2: the parameters (traits)
  // carry across via ad_value, the double environment and control are copied.
  // The gradient driver lifts the resident (double) patch to the active scalar
  // with this before seeding and running. Constrained to S2 != value_type: the
  // lift is always to an active scalar, and self-rebind is just a copy. Gating
  // it also keeps odelia's stiff RODAS Jacobian (which probes rebind_from at the
  // system's own scalar) on the double path, so the build stays double-only
  // until the active ODE stepper is wired.
  template <class S2, class = std::enable_if_t<!std::is_same_v<S2, value_type>>>
  auto rebind_from() const {
    return typename PatchRebind_::template rebind<S2>(
        parameters.template rebind_from<S2>(),
        environment.template rebind_from<S2>(), control);
  }

  // Handles to the differentiable trait parameters, concatenated across species
  // (species-major, each in the strategy's fixed field order). Column j of a
  // Jacobian is d(out)/d(handle j); the driver seeds a chosen subset by index.
  // Every cohort aliases its species' one shared strategy, so seeding through
  // these reaches all cohorts, including later introductions.
  std::vector<value_type*> ad_parameters() {
    std::vector<value_type*> ret;
    for (auto& s : species) {
      std::vector<value_type*> p = s.ad_parameters();
      ret.insert(ret.end(), p.begin(), p.end());
    }
    return ret;
  }

  // Registered initial-state leaves the driver seeds via DifferentiationTargets
  // ics: one active birth-rate scale per species, in species order. Seeding entry
  // s differentiates the run w.r.t. species s's birth rate (the density-feedback
  // axis). Size-distribution IC sensitivity stays a later target (Appendix A.4).
  std::vector<value_type*> ad_initial_state() {
    std::vector<value_type*> ret;
    ret.reserve(species.size());
    for (auto& s : species) ret.push_back(s.ad_birth_rate());
    return ret;
  }

  // Re-derive every species' prepare_strategy() quantities under the current
  // (seeded) parameters. This is the step that carries a seeded trait into the
  // otherwise-frozen derived quantities (eta_c, height_0, height_0_inverse,
  // area_leaf_0, canopy_shape, the bound assimilation_fn); without it those
  // differentiate to zero. reset() calls it on the active path.
  void ad_prepare() {
    for (auto& s : species) s.prepare_strategy();
  }

  //Try using pointer in place of object itself
  double time() const {return environment.time;}
  double get_area() const { return area;}
  double height_max() const;

  double compute_competition(double height) const;
  // Active-scalar total competition at height, summed over species. The resident
  // replay's canopy recompute builds the light spline from this, so a trait
  // re-shades the stand through area_leaf. Bit-identical value at S = double.
  value_type compute_competition_ad(double height) const;

  // * Lifetime fitness / offspring production
  // These are patch-level quantities: each integrates the per-node weighted
  // net reproduction over a species' node-introduction times.
  // Integrate lifetime fitness of a species' nodes, scaled per node.
  double net_reproduction_ratio_for_species(size_t species_index,
                                            std::vector<double> const& scalars) const;
  // Offspring production: fitness scaled by the birth rate over time.
  std::vector<double> offspring_production() const;
  // Overall fitness (unscaled, scalars == 1).
  std::vector<double> net_reproduction_ratios() const;
  // Sum of offspring produced across all species.
  double total_offspring_production() const;
  // Per-node reproduction integration error for each species.
  std::vector<std::vector<double>> net_reproduction_ratio_errors() const;

  // * Schedule-refinement error collection
  // Sample the competition error for each species introduced this step and
  // fold it into the running per-node max (ignoring NA, matching na.rm=TRUE in
  // R). Accumulated across the run; cleared by reset().
  void collect_competition_errors(const std::vector<size_t>& added);
  // Combine the competition error (sampled during the run) with the
  // reproduction error (computed now) into a single per-node error vector per
  // species. An all-NA node yields -Inf. Drives schedule refinement.
  std::vector<std::vector<double>> refinement_error_by_node() const;

  void introduce_new_node(size_t species_index);
  void introduce_new_nodes(const std::vector<size_t>& species_index);

  // Open to better ways to test whether nodes have been introduced
  int node_ode_size() const {
    int node_ode_size = ode_size() - environment.ode_size();
    return(node_ode_size);
  }

  const species_type& at_species(size_t species_index) const {
    return species[species_index];
  }

  // Patch disturbance
  Disturbance_Regime* survival_weighting;

  // * ODE interface
  // Caluclate size of ode system (number of equations). Is constantly changing as 
  // new nodes are introduced into size-density distreibution 
  size_t ode_size() const;
  // How many auxiallary variables are we tracking. These are being collected but
  // are not part of core ode system
  size_t aux_size() const;
  double ode_time() const;

  // Retrieve ode state from patch and save into the ode solver
  template <class It> It ode_state(It it) const;
  // Retrieve ode rates from patch and save into the ode solver
  template <class It> It ode_rates(It it) const;
  // Retrieve auxillary variables and save into the ode solver
  odelia::ode::iterator ode_aux(odelia::ode::iterator it) const;

  // Returns state in structure format as opposed to single 
  // vector as given by ode_state
  Rcpp::List r_get_state() const;

  // Two overloads: (it, time) recomputes the environment (resident); (it, stage)
  // reads the frozen environment by RK stage (mutant), via has_recorded_field().
  template <class It> It set_ode_state(It it, double time);
  template <class It> It set_ode_state(It it, int index);

  // * R interface
  // Data accessors:
  double r_density(double time) const {return survival_weighting->r_density(time);}
  double r_pr_survival(double time) const {return survival_weighting->pr_survival(time);}
  double r_disturbance_mean_interval() const {return survival_weighting->r_mean_interval();}
  double r_survival_weighting_cdf(double time) const {return survival_weighting->cdf(time);}
  double r_survival_weighting_icdf(double prob) const {return survival_weighting->icdf(prob);}

  parameters_type r_parameters() const {return parameters;}
  environment_type r_environment() const {return environment;}
  std::vector<species_type> r_species() const {return species;}
  std::vector<double> r_compute_competition_effect_error_by_node_for_species_i(size_t species_index) const;
  void r_set_time(double time);
  void r_set_state(double time,
                   const std::vector<double>& state,
                   const std::vector<size_t>& n,
                   const std::vector<double>& light_availability);
  void r_introduce_new_node(util::index species_index) {
    introduce_new_node(species_index.check_bounds(size()));
  }
  species_type r_at(util::index species_index) const {
    at(species_index.check_bounds(size()));
  }
  // These are only here because they wrap private functions.
  void r_compute_environment() {compute_environment(false);}
  void r_compute_rates() {
    environment_ptr = &environment;
    compute_rates();
  }

  // env. cache for assembly
  std::vector<double> step_history{0.0};  // always start at zero
  std::vector<std::vector<environment_type>> environment_history;
  std::vector<environment_type> environment_cache;

  // L2 slice: the light-spline knot POSITIONS per accepted step, recorded on the
  // double pass independently of the L3 field cache above. The active resident
  // replay rebuilds the canopy on these frozen positions with the active cohorts.
  std::vector<std::vector<double>> knot_history;
  std::vector<double> knot_cache;    // stage-0 knots of the step being recorded
  std::vector<double> current_knots; // this step's frozen knots, on the replay

  // odelia Replayable hooks: record_stage/record_ode_step cache the environment
  // (L3) and knot positions (L2) during a resident run; replay_step locates the
  // frozen step on a mutant replay or loads its knots on the resident replay.
  void record_stage(int step);
  void record_ode_step();
  void replay_step();

  // The has_recorded_field() query odelia's derivs dispatches on: true while a
  // mutant replays the frozen environment (set by set_mutant once it is recorded).
  // Empty on the resident replay, so derivs recomputes the canopy with the active
  // cohorts (the self-shading channel) instead of reading it frozen.
  bool has_recorded_field() const { return use_cached_environment; }

  // Record-mode flag: caches environment history during a resident run.
  bool save_RK45_cache;

  bool use_cached_environment = false;

  bool is_mutant_run = false;

  // Resident (self-shading) replay: rebuild the canopy live on the recorded knots
  // with the active cohorts. L3 stays empty (has_recorded_field() false), so the
  // trait re-shades the stand. Distinct from set_mutant (frozen-field read).
  bool replay_knots = false;
  void set_resident_replay() {
    replay_knots = true;
    save_RK45_cache = false;
    idx = 0;
  }

  void set_mutant();
  void add_strategies(std::vector<strategy_type> strategies);
  void overwrite_strategies(std::vector<strategy_type> strategies);

private:
  int idx = 0; // used to access environment cache for mutant runs
  void compute_environment(bool rescale);
  void compute_rates();

  // Seed the patch from parameters.initial_state (nodes + birth bookkeeping)
  // when present; called from reset(). Sets environment.time = initial_time.
  void set_initial_state();
  // Guard against initial conditions whose per-node log-density rates are so
  // large they would drive densities to non-finite values within a few steps.
  void check_initial_density_rates() const;
  // Guard against the SCM equations running away mid-integration: a cohort
  // density (exp(log_density)) or an environment state (e.g. TF24 soil water)
  // going non-finite. Called each derivs evaluation, before the non-finite
  // value can propagate into competition, resource uptake, or physiology and
  // surface as an opaque downstream error.
  void check_finite_ode_state() const;

  parameters_type parameters;

  double area;
  environment_type environment;
  std::vector<species_type> species;

  //TODO(#476): Move into environment?
  std::vector<double> resource_depletion;

  environment_type* environment_ptr;

  Control control;

  // Per-species running max of the competition error per node, accumulated
  // across the run via collect_competition_errors(). Entries start at -Inf and
  // ignore NA contributions, matching apply(., 2, max, na.rm=TRUE) in R.
  std::vector<std::vector<double>> competition_error_by_node;
};

template <typename T, typename E>
Patch<T,E>::Patch(parameters_type p, environment_type e, Control c)
  : parameters(p),
    area(p.patch_area),
    environment(e),
    control(c),
    environment_cache(6) {  // length of odelia::ode::Step
  
  parameters.validate();

  save_RK45_cache = control.save_RK45_cache;
  survival_weighting = p.disturbance;

  // Configure the light profile's shading model before the first
  // compute_environment() in reset(). No-op for environments without a light
  // profile (only FF16 implements alternative shading models).
  environment.set_shading_model(control.shading_model,
                                control.ppa_layer_optical_depth,
                                control.ppa_layer_smoothing);

  add_strategies(parameters.strategies);

  reset();
}

template <typename T, typename E>
void Patch<T,E>::overwrite_strategies(std::vector<strategy_type> strategies) {
  species.clear();
  add_strategies(strategies);
}

template <typename T, typename E>
void Patch<T,E>::add_strategies(std::vector<strategy_type> strategies) {
  for (auto i = 0; i < strategies.size(); ++i) {
		auto s = strategies[i];
    s.control = control; // Overwrite to take the patch control object 
    auto spec = Species<T,E>(s);
    species.push_back(spec);
  }
}

template <typename T, typename E>
void Patch<T,E>::set_mutant() {
    if (environment_history.empty()) {
       util::stop("Run a resident first to generate a competitve landscape");
    }

    is_mutant_run = true;
    save_RK45_cache = false;
    use_cached_environment = true;
  idx = 0;
}

template <typename T, typename E>
void Patch<T,E>::reset() {
  // Active path only: re-derive each strategy's prepare_strategy() quantities
  // from its own (seeded) parameters before re-initialising cohorts, so a seeded
  // trait reaches the frozen derived quantities and the new cohorts built below
  // (trap 1). Cohort state is then rebuilt from these parameters, never an
  // external double snapshot, so the seed survives reset (trap 2). The double
  // resident path is untouched.
  if constexpr (!std::is_same_v<value_type, double>) {
    ad_prepare();
  }

   for (auto& s : species) {
    s.clear();
    // allocate variables for tracking resource consumption
    s.resize_consumption_rates(environment.ode_size());
  }

  // resize to species count
  resource_depletion.reserve(environment.ode_size());

  // compute ephemeral effects like light_availability
  environment.clear();

  if (!parameters.initial_state.empty()) {
    // Seed the patch from an exported state / initial size distribution.
    // set_initial_state() does its own compute_environment(false)/compute_rates
    // (with the real node heights, no rescale), so skip the empty-patch path.
    set_initial_state();
    check_initial_density_rates();
  } else {
    compute_environment(false);

    // compute effects of resource consumption
    environment_ptr = &environment;
    compute_rates();
  }

  // clear accumulated per-node competition error
  competition_error_by_node.assign(species.size(), {});
}

// Seed the patch from parameters.initial_state. Introduces the requested number
// of nodes per species, loads the flat ODE state (nodes + environment) at the
// resume time, restores per-node birth bookkeeping (not part of the ODE state
// but feeds the rates and lifetime-fitness integrals), then computes the
// environment and rates. We load the state directly (rather than via the
// double-arg set_ode_state) so the first environment build is a full
// compute_environment(false): a rescale of the not-yet-built light spline would
// read uninitialised grid state.
template <typename T, typename E>
void Patch<T,E>::set_initial_state() {
  const size_t n_species = species.size();
  util::check_length(parameters.n_initial_cohorts.size(), n_species);

  size_t total_nodes = 0;
  for (size_t i = 0; i < n_species; ++i) {
    for (size_t j = 0; j < parameters.n_initial_cohorts[i]; ++j) {
      species[i].introduce_new_node();
    }
    total_nodes += parameters.n_initial_cohorts[i];
  }

  // Load the flat ODE state (all nodes, then environment).
  util::check_length(parameters.initial_state.size(), ode_size());
  odelia::ode::const_iterator it = parameters.initial_state.begin();
  it = odelia::ode::set_ode_state(species.begin(), species.end(), it);
  it = environment.set_ode_state(it);
  environment.time = parameters.initial_time;

  // Restore birth bookkeeping per node, sliced per species from the flat
  // parameter vectors. Skipped (left at defaults) when not supplied, e.g. a
  // from-scratch distribution seeded at patch age 0.
  if (!parameters.initial_node_times.empty()) {
    util::check_length(parameters.initial_node_times.size(), total_nodes);
    util::check_length(parameters.initial_patch_density.size(), total_nodes);
    util::check_length(parameters.initial_pr_patch_survival.size(), total_nodes);
    auto t_it = parameters.initial_node_times.begin();
    auto d_it = parameters.initial_patch_density.begin();
    auto s_it = parameters.initial_pr_patch_survival.begin();
    for (size_t i = 0; i < n_species; ++i) {
      const size_t n = parameters.n_initial_cohorts[i];
      species[i].set_birth_state(std::vector<double>(t_it, t_it + n),
                                 std::vector<double>(d_it, d_it + n),
                                 std::vector<double>(s_it, s_it + n));
      t_it += n;
      d_it += n;
      s_it += n;
    }
  }

  // Build the environment from the real node heights (full recompute, no
  // rescale) and compute rates for the seeded population.
  compute_environment(false);
  environment_ptr = &environment;
  compute_rates();
}

template <typename T, typename E>
void Patch<T,E>::check_initial_density_rates() const {
  for (const auto& s : species) {
    std::vector<double> rates = s.r_log_density_rates();
    if (std::any_of(rates.begin(), rates.end(),
                    [](double v) { return v < -100; })) {
      util::stop("Rates of initial node densities exceed ~1e43 and will likely "
                 "produce non-finite densities; provide more plausible initial "
                 "conditions (smaller sizes and/or lower densities).");
    }
  }
}

template <typename T, typename E>
void Patch<T,E>::check_finite_ode_state() const {
  // The two failure modes of the same runaway (issue #550), caught here so they
  // fail with an actionable message instead of an opaque downstream one:
  //
  //  (1) A cohort density overflowing to +Inf, which then poisons the
  //      competition integral ("Detected non-finite contribution").
  //  (2) The coupled environment state (TF24 soil water) going non-finite,
  //      because the density-weighted resource uptake diverges as density
  //      grows; the RK stage arithmetic can produce a non-finite soil state a
  //      step before density itself overflows, surfacing as a non-finite soil
  //      potential ("non-finite psi_soil").
  //
  // Density ceiling (~1e43): the same "already unphysical, will drive
  // non-finite values" bar used by check_initial_density_rates. Catching the
  // runaway at this magnitude -- before density reaches a literal +Inf --
  // widens coverage of mode (1) and heads off some instances of mode (2).
  const double log_density_ceiling = 50.0;
  for (size_t i = 0; i < species.size(); ++i) {
    for (auto it = species[i].node_begin(); it != species[i].node_end(); ++it) {
      if (!util::is_finite(it->get_density()) ||
          it->get_log_density() > log_density_ceiling) {
        // The size-density characteristic equation integrates
        //   d(log density)/dt = -d(growth)/d(height) - mortality
        // (see Node::compute_rates). When growth rate falls steeply with size
        // -- e.g. a cohort hitting an extreme drought trough -- the gradient
        // term can spike sharply positive, so log_density integrates upward
        // until density = exp(log_density) overflows to +Inf. The individual's
        // physiology (height, leaf area, per-individual competition) stays
        // finite; it is the cohort *density* that blows up. Smaller ODE steps
        // do not help (the divergence is in the equations, not the stepper), so
        // fail here with an actionable message rather than letting the +Inf
        // propagate into the competition integral or density-weighted resource
        // uptake, where it surfaces as an opaque downstream error.
        util::stop("Non-finite cohort density in the SCM size-density "
                   "(characteristic) equations: species " +
                   util::to_string(i + 1) + " has a node with density=" +
                   util::to_string(it->get_density()) + " (log_density=" +
                   util::to_string(it->get_log_density()) + ", height=" +
                   util::to_string(it->height()) + ") at time=" +
                   util::to_string(environment.time) +
                   ". The density derivative -d(growth)/d(height) - mortality "
                   "can grow without bound when growth rate falls steeply with "
                   "size under rapidly changing or extreme environmental "
                   "forcing, driving density to overflow. Try a shorter "
                   "max_patch_lifetime or less extreme environmental drivers.");
      }
    }
  }

  // Mode (2): a non-finite environment state. `vars.states` is generic across
  // environments (size 0 for light-only environments like FF16, so this loop is
  // a no-op there and cannot false-positive). For TF24 these are the per-depth
  // soil-water states.
  const Internals& env_vars = environment.vars;
  for (size_t i = 0; i < env_vars.state_size; ++i) {
    if (!util::is_finite(env_vars.states[i])) {
      util::stop("Non-finite environment state (index " + util::to_string(i) +
                 " = " + util::to_string(env_vars.states[i]) + ") at time=" +
                 util::to_string(environment.time) +
                 ". For TF24 this is a soil-water state driven non-finite by the "
                 "density-weighted resource uptake as a cohort density runs away "
                 "in the SCM size-density equations (the same failure that "
                 "otherwise overflows density to +Inf; see #550). Try a shorter "
                 "max_patch_lifetime or less extreme environmental drivers.");
    }
  }
}

template <typename T, typename E>
double Patch<T,E>::height_max() const {
  double ret = 0.0;
  for (size_t i = 0; i < species.size(); ++i) {
    if (!is_mutant_run) {
      ret = std::max(ret, species[i].height_max());
    }
  }
  return ret;
}

template <typename T, typename E>
double Patch<T,E>::compute_competition(double height) const {
  double tot = 0.0;
  for (size_t i = 0; i < species.size(); ++i) {
    if (!is_mutant_run) {
      tot += species[i].compute_competition(height) / area;
    }
  }
  return tot;
}

template <typename T, typename E>
typename Patch<T,E>::value_type
Patch<T,E>::compute_competition_ad(double height) const {
  value_type tot(0.0);
  for (size_t i = 0; i < species.size(); ++i) {
    if (!is_mutant_run) {
      tot += species[i].compute_competition_ad(height) / area;
    }
  }
  return tot;
}

template <typename T, typename E>
std::vector<double> Patch<T,E>::r_compute_competition_effect_error_by_node_for_species_i(size_t species_index) const {
  const double tot_competition_effect = compute_competition(0.0);
  return species[species_index].r_compute_competition_effect_by_nodes_error(tot_competition_effect);
}

// Integrate over lifetime fitness of individual nodes, scaled per node.
template <typename T, typename E>
double Patch<T,E>::net_reproduction_ratio_for_species(
    size_t species_index, std::vector<double> const& scalars) const {
  auto net_prod = species[species_index].net_reproduction_ratio_by_node_weighted();
  auto const times = species[species_index].node_times();
  auto net_prod_scaled = std::vector<double>(times.size());
  for (size_t i = 0; i < times.size(); ++i) {
    net_prod_scaled[i] = net_prod[i] * scalars[i];
  }
  return util::trapezium(times, net_prod_scaled);
}

// Offspring production, equal to overall fitness scaled by the birth rate.
template <typename T, typename E>
std::vector<double> Patch<T,E>::offspring_production() const {
  auto ret = std::vector<double>(species.size());
  for (size_t i = 0; i < species.size(); ++i) {
    // scale by birth rate function over time
    auto const times = species[i].node_times();
    auto scalars = std::vector<double>(times.size());
    for (size_t j = 0; j < times.size(); ++j) {
      scalars[j] = species[i].extrinsic_drivers().evaluate("birth_rate", times[j]);
    }
    ret[i] = net_reproduction_ratio_for_species(i, scalars);
  }
  return ret;
}

// Overall fitness (no scaling, ie scalars set to 1.0).
template <typename T, typename E>
std::vector<double> Patch<T,E>::net_reproduction_ratios() const {
  auto ret = std::vector<double>(species.size());
  for (size_t i = 0; i < species.size(); ++i) {
    auto scalars = std::vector<double>(species[i].size(), 1.0);
    ret[i] = net_reproduction_ratio_for_species(i, scalars);
  }
  return ret;
}

// Sum up all offspring produced.
template <typename T, typename E>
double Patch<T,E>::total_offspring_production() const {
  double total = 0.0;
  std::vector<double> offspring = offspring_production();
  for (size_t i = 0; i < species.size(); ++i) {
    total += offspring[i];
  }
  return total;
}

// Check integration errors for each species' reproduction integral.
template <typename T, typename E>
std::vector<std::vector<double>> Patch<T,E>::net_reproduction_ratio_errors() const {
  std::vector<std::vector<double>> ret;
  double total_offspring = total_offspring_production();
  for (size_t i = 0; i < species.size(); ++i) {
    ret.push_back(util::local_error_integration(
        species[i].node_times(),
        species[i].net_reproduction_ratio_by_node_weighted(),
        total_offspring));
  }
  return ret;
}

// Sample the competition error for each species introduced this step and fold
// it into the running per-node max (ignoring NA, matching na.rm=TRUE in R).
template <typename T, typename E>
void Patch<T,E>::collect_competition_errors(const std::vector<size_t>& added) {
  for (size_t idx : added) {
    std::vector<double> v =
        r_compute_competition_effect_error_by_node_for_species_i(idx);
    std::vector<double>& acc = competition_error_by_node[idx];
    if (acc.size() < v.size()) {
      acc.resize(v.size(), -std::numeric_limits<double>::infinity());
    }
    for (size_t j = 0; j < v.size(); ++j) {
      if (!ISNAN(v[j])) {
        acc[j] = std::max(acc[j], v[j]);
      }
    }
  }
}

// Combine the competition error (sampled during the run) with the reproduction
// error (computed now) into a single per-node error vector per species. An
// all-NA node yields -Inf, matching apply(rbind(...), 2, max, na.rm=TRUE) in R.
template <typename T, typename E>
std::vector<std::vector<double>> Patch<T,E>::refinement_error_by_node() const {
  std::vector<std::vector<double>> repro = net_reproduction_ratio_errors();
  std::vector<std::vector<double>> ret(species.size());
  for (size_t i = 0; i < species.size(); ++i) {
    const std::vector<double>& comp = competition_error_by_node[i];
    const std::vector<double>& rep = repro[i];
    const size_t n = species[i].size();
    std::vector<double> tot(n, -std::numeric_limits<double>::infinity());
    for (size_t j = 0; j < n; ++j) {
      if (j < comp.size() && !ISNAN(comp[j])) {
        tot[j] = std::max(tot[j], comp[j]);
      }
      if (j < rep.size() && !ISNAN(rep[j])) {
        tot[j] = std::max(tot[j], rep[j]);
      }
    }
    ret[i] = tot;
  }
  return ret;
}

// Pre-compute environment, as shaped by residents
// Creates splines of resource availability
template <typename T, typename E>
void Patch<T,E>::compute_environment(bool rescale) {

  if (!(size() > 0) || is_mutant_run) {
    return;
  }

  if constexpr (std::is_same_v<value_type, double>) {
    // Resident (double) path: fit the light spline to the double competition,
    // adaptively or by rescale -- unchanged.
    auto f = [&](double x) -> double { return compute_competition(x); };
    environment.compute_environment(f, height_max(), rescale);
  } else {
    // Active path: the competition carries the trait derivative through
    // area_leaf. On the resident replay the canopy is rebuilt on the frozen
    // recorded knots (self-shading flows without moving a node); before any step
    // is replayed (empty knots) or off the replay it falls back to the adaptive
    // build. The invasion twin never reaches here (is_mutant_run above).
    auto f = [&](double x) -> value_type { return compute_competition_ad(x); };
    if (replay_knots && !current_knots.empty()) {
      environment.compute_environment_fixed(f, current_knots);
    } else {
      environment.compute_environment(f, height_max(), rescale);
    }
  }
}


template <typename T, typename E>
void Patch<T,E>::compute_rates() {

  // Computes rates of change for the patch, including all the component species
  // While the patch has an `environment`, the rates here are calculated from
  // the env_ptr, which is a pointer to an environment object
  //  -- for the resident the pointer points to the internal environment object
  //  -- for a mutant, the pointer points to a cached environment object
  double time_ = environment_ptr->time;

  double pr_patch_survival = survival_weighting->pr_survival(time_);
  for (size_t i = 0; i < size(); ++i) {
    double birth_rate = species[i].extrinsic_drivers().evaluate("birth_rate", time_);

    // Pass the environment that pointer is tracking into compute rates.
    species[i].compute_rates(*environment_ptr, pr_patch_survival, birth_rate);
  }

  resource_depletion.reserve(environment_ptr->ode_size());
  for(size_t i = 0; i < environment_ptr->ode_size(); i++) {
    double resource_consumed = std::accumulate(species.begin(), species.end(), 0.0, [i](double r, const species_type& s) {
      return r + s.consumption_rate(i); // accumulates r from zero
    });

    resource_depletion.push_back(resource_consumed/area);
  }
  

  environment_ptr->compute_rates(resource_depletion);

  //todo do we need to clear this every step?
  resource_depletion.clear();

}

// TODO(#478): We should only be recomputing the light environment for the
// points that are below the height of the seedling -- not the entire
// light environment; probably worth just doing a rescale there?
template <typename T, typename E>
void Patch<T,E>::introduce_new_node(size_t species_index) {
  
  species[species_index].introduce_new_node();

  compute_environment(false);
}

template <typename T, typename E>
void Patch<T,E>::introduce_new_nodes(const std::vector<size_t>& species_index) {
  // Record introduction time and patch-age density on each node as it is
  // introduced, so lifetime-fitness calcs need not look these up later.
  const double t = time();
  const double patch_density = survival_weighting->density(t);
  for (size_t i : species_index) {
    species[i].introduce_new_node(t, patch_density);
  }

  compute_environment(false);
}

template <typename T, typename E>
void Patch<T,E>::r_set_time(double time) {
  environment.time = time;
}

// Arguments here are:
//   time: time
//   state: vector of ode state; we'll pass an iterator with that in
//   n: number of *individuals* of each species
template <typename T, typename E>
void Patch<T,E>::r_set_state(double time,
                           const std::vector<double>& state,
                           const std::vector<size_t>& n,
                           const std::vector<double>& light_availability) {
  const size_t n_species = species.size();
  util::check_length(n.size(), n_species);
  reset();
  for (size_t i = 0; i < n_species; ++i) {
    for (size_t j = 0; j < n[i]; ++j) {
      species[i].introduce_new_node();
    }
  }
  util::check_length(state.size(), ode_size());
  set_ode_state(state.begin(), time);
  environment.r_init_interpolators(light_availability);
}

// ODE interface
template <typename T, typename E>
size_t Patch<T,E>::ode_size() const {
  return odelia::ode::ode_size(species.begin(), species.end()) + environment.ode_size();
}

template <typename T, typename E>
size_t Patch<T,E>::aux_size() const {
  // TODO(#478): Is this useful for environment vectors?
  // no use for auxiliary environment variables (yet)
  return odelia::ode::aux_size(species.begin(), species.end());// + environment.ode_size();
}

template <typename T, typename E>
double Patch<T,E>::ode_time() const {
  return time();
}

// First set_ode_state function is for resident runs. Second is for mutant runs.
// Templated on the state iterator so the cohort hierarchy is stepped at the ODE
// vector's scalar (double resident, active under a gradient); the recursive walk
// replaces odelia's double-only free functions.
template <typename T, typename E>
template <class It>
It Patch<T,E>::set_ode_state(It it, double time) {

  // Set ode states
  for (auto& s : species) { it = s.set_ode_state(it); }
  it = environment.set_ode_state(it);

  // update time
  environment.time = time;

  // Catch a runaway size-density equation (non-finite cohort density or
  // environment state) before it feeds into competition, resource uptake, or
  // physiology below and surfaces as an opaque downstream error (issue #550).
  check_finite_ode_state();

  // Pre-compute environment, as shaped by residents
  compute_environment(true);
  environment_ptr = &environment;

  // Compute rates of change
  compute_rates();
  return it;
}

// used for mutant runs
// -- differs from above in that an index is passed in as argument
// -- environments are loaded from ODE history, instead of being calculated 
template <typename T, typename E>
template <class It>
It Patch<T,E>::set_ode_state(It it, int index) {

  for (auto& s : species) { it = s.set_ode_state(it); }

  // A cohort introduced on the final recorded step (birth >= N) reaches the last
  // ode time, one step_history entry past the last recorded RK-stage environment
  // (step_history carries the extra t=0 slot, so it is one longer). Read the final
  // frozen environment there: the boundary cohort is established at seed height and
  // never stepped past the record (the zero-height fix; without it the active crown
  // derivative dereferences a slot beyond environment_history and segfaults).
  int env_idx = idx;
  if (!environment_history.empty() &&
      static_cast<size_t>(env_idx) >= environment_history.size()) {
    env_idx = static_cast<int>(environment_history.size()) - 1;
  }
  environment_ptr = &(environment_history[env_idx][index]);
  environment.time = environment_ptr->time;

  // increment the iterator by an appropriate amount, but don't actually do anything in the env
  for (size_t i = 0; i < environment_ptr->ode_size(); i++) {*it++;}

  compute_rates();
  return it;
}

// Commit an accepted ODE step's cached environments (the 6 RK stages) to the
// history. odelia calls this per accepted step; a no-op unless recording.
template <typename T, typename E>
void Patch<T,E>::record_ode_step() {
  if(save_RK45_cache) {
    step_history.push_back(time());
    environment_history.push_back(environment_cache);
    knot_history.push_back(knot_cache);
  }
}

// Cache the environment at one RK stage (L3) and, at stage 0, the light-spline
// knot positions for this step (L2). odelia calls this per stage; a no-op unless
// recording. light_knots() is only present on light-profile environments, so the
// knot capture is guarded (a no-op for environments without one).
template <typename T, typename E>
void Patch<T,E>::record_stage(int step) {
  if(save_RK45_cache) {
    if(step == 0) {
      environment_cache.clear();
      if constexpr (requires (const environment_type& e) { e.light_knots(); }) {
        knot_cache = environment.light_knots();
      }
    }
    environment_cache.push_back(environment);
  }
}

// Locate this step's cached environment for a mutant replay. odelia calls this
// before each fixed step; the (it, stage) set_ode_state then reads by idx.
template <typename T, typename E>
void Patch<T,E>::replay_step() {
  // Both the mutant (frozen field) and resident (frozen knots) replays locate
  // this step in the recorded history; the resident additionally loads the step's
  // knot positions so the stages recompute the canopy on them.
  if (!use_cached_environment && !replay_knots) {
    return;
  }

  // Minor optimization to check the current and next index before doing a search, as the most common case is that the ODE solver is stepping through the cached environments in order. If the call sequence was not strictly sequential, we fallback to a search through the step history to find the correct environment.
  const double t = time();
  const size_t n = step_history.size();

  // Fast path: step_to() advances through ode_times in order, so this is
  // usually either the current cached step index or the next one.
  if (static_cast<size_t>(idx) < n && util::identical(step_history[idx], t)) {
    // idx already points at this step.
  } else if (static_cast<size_t>(idx + 1) < n &&
             util::identical(step_history[idx + 1], t)) {
    ++idx;
  } else {
    // Fallback to search if the call sequence was not strictly sequential.
    auto step = std::find(step_history.begin(), step_history.end(), t);
    if (step == step_history.end()) {
      util::stop("ODE time not found in step history");
    }
    idx = static_cast<int>(std::distance(step_history.begin(), step));
  }

  if (replay_knots && static_cast<size_t>(idx) < knot_history.size()) {
    current_knots = knot_history[idx];
  }
}

template <typename T, typename E>
template <class It>
It Patch<T,E>::ode_state(It it) const {
  for (const auto& s : species) { it = s.ode_state(it); }
  it = environment.ode_state(it);
  return it;
}

template <typename T, typename E>
Rcpp::List Patch<T, E>::r_get_state() const
{

  // Aseemble commkunity state, icnluding auxiallry variables
  Rcpp::List community_state;
  for (size_t i = 0; i < species.size(); ++i)
  {
    community_state.push_back(species[i].r_get_state());
  }

  return Rcpp::List::create(_["time"] = time(),
                            _["patch_density"] = r_density(time()),
                            _["species"] = community_state,
                            _["env"] = environment.r_get_state());
}

template <typename T, typename E>
template <class It>
It Patch<T,E>::ode_rates(It it) const {
  for (const auto& s : species) { it = s.ode_rates(it); }
  it = environment.ode_rates(it);
  return it;
}

template <typename T, typename E>
odelia::ode::iterator Patch<T,E>::ode_aux(odelia::ode::iterator it) const {
  it = odelia::ode::ode_aux(species.begin(), species.end(), it);
  return it;
}

}

#endif
