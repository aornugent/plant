// -*-c++-*-
#ifndef PLANT_PLANT_PATCH_H_
#define PLANT_PLANT_PATCH_H_

#include <plant/parameters.h>
#include <plant/species.h>
#include <plant/util.h>
#include <odelia/ode_interface.hpp>

#include <plant/disturbance_regime.h>

#include <algorithm>
#include <limits>

using namespace Rcpp;

namespace plant {

// Detect environments that support the exact separable competition field (only
// those with a Yokozawa deep-crown light kernel -- K93 today). For those, the
// patch assembles the field from the cohort population instead of fitting the
// light spline; others keep the spline. Structurally scopes the exact-field path.
template <class E2, class = void>
struct env_has_competition_field : std::false_type {};
template <class E2>
struct env_has_competition_field<
    E2, std::void_t<decltype(std::declval<E2&>().clear_competition_field())>>
    : std::true_type {};

template <typename T, typename E>
class Patch {
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

  //Try using pointer in place of object itself
  double time() const {return environment.time;}
  double get_area() const { return area;}
  double height_max() const;

  value_type compute_competition(double height) const;

  // * Lifetime fitness / offspring production
  // These are patch-level quantities: each integrates the per-node weighted
  // net reproduction over a species' node-introduction times.
  // Integrate lifetime fitness of a species' nodes, scaled per node.
  value_type net_reproduction_ratio_for_species(size_t species_index,
                                                std::vector<double> const& scalars) const;
  // Offspring production: fitness scaled by the birth rate over time.
  std::vector<value_type> offspring_production() const;
  // Overall fitness (unscaled, scalars == 1).
  std::vector<value_type> net_reproduction_ratios() const;
  // Sum of offspring produced across all species.
  value_type total_offspring_production() const;
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
  template <typename It> It ode_state(It it) const;
  // Retrieve ode rates from patch and save into the ode solver
  template <typename It> It ode_rates(It it) const;
  // Retrieve auxillary variables and save into the ode solver
  template <typename It> It ode_aux(It it) const;

  // Returns state in structure format as opposed to single 
  // vector as given by ode_state
  Rcpp::List r_get_state() const;

  // Set patch state from the solver's estimate. Two overloads, chosen by
  // odelia::ode::derivs on has_recorded_field():
  //   - (it, time):  resident -- recompute the environment at this time.
  //   - (it, index): mutant   -- read the recorded environment for this RK stage.
  template <typename It> It set_ode_state(It it, double time);
  template <typename It> It set_ode_state(It it, int index);

  // The differentiable inputs the gradient driver seeds a subset of (§8.1):
  // species-major, each species' low-level strategy parameters. A resident
  // patch's cohort birth states are taped intermediates of the population (set
  // by compute_initial_conditions each step), not independent inputs, so it
  // seeds no initial state of its own.
  std::vector<value_type*> ad_parameters() {
    std::vector<value_type*> ptrs;
    for (auto& sp : species) {
      auto block = sp.ad_parameters();
      ptrs.insert(ptrs.end(), block.begin(), block.end());
    }
    return ptrs;
  }
  std::vector<value_type*> ad_initial_state() { return {}; }

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

  // Recorded environment for mutant replay. The resident (recording) pass stores
  // the accepted step times and, per step, a whole-environment snapshot at each of
  // the six RK stages; the mutant (replay) pass reads them back as fixed doubles.
  // step_history[0] == 0.
  std::vector<double> step_history{0.0};
  std::vector<std::vector<environment_type>> environment_history;  // [step][RK stage]
  std::vector<environment_type> environment_cache;                 // this step's stages

  // The odelia::ode::Replayable hooks the Solver drives: record_stage per RK
  // stage, record_ode_step per accepted step, replay_step per step on the replay
  // pass. record_* append only while recording; replay_step only when a recording
  // is present.
  void record_stage(int stage);
  void record_ode_step();
  void replay_step();

  // Is the recorded environment populated for replay? True on a mutant pass (a
  // recording is present and this pass is not the one producing it), false while a
  // resident records or when nothing was recorded. odelia::ode::derivs reads this
  // to route between the recompute and replay set_ode_state overloads.
  bool has_recorded_field() const {
    return !recording && !environment_history.empty();
  }

  // Set while a resident produces the recording: enabled up front for a resident
  // that will be replayed, cleared by set_mutant() so the mutant reads it back.
  bool recording;

  // A mutant (invasion) run: the focal strategy reads the recorded background and
  // does not shape the field (no self-shading).
  bool is_mutant_run = false;

  void set_mutant();
  void add_strategies(std::vector<strategy_type> strategies);
  void overwrite_strategies(std::vector<strategy_type> strategies);

private:
  int idx = 0; // used to access environment cache for mutant runs
  void compute_environment(bool rescale);
  void assemble_competition_field();  // exact separable field (deep-crown kernels)
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
  // Per-ODE-channel resource consumption handed to environment.compute_rates.
  // Carries value_type so the resident soil coupling (uptake -> soil state)
  // differentiates; double on the production path.
  std::vector<value_type> resource_depletion;

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

  recording = control.save_RK45_cache;
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
    recording = false;  // stop recording so has_recorded_field() reads it back
  idx = 0;
}

template <typename T, typename E>
void Patch<T,E>::reset() {
   for (auto& s : species) {
    // Re-derive each strategy's precomputed quantities from its current
    // parameters before the run: the gradient driver seeds the parameters then
    // reset()s, so the birth size / canopy shape / eta_c must be recomputed here
    // to carry their parameter derivative (cf. IndividualRunner::reset()). On the
    // double path this is idempotent.
    s.prepare_strategy();
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
      // A coincident cohort on the chart has log_density = -inf (density 0) from
      // odelia's zero-spacing convention, which passes both checks below
      // (is_finite(0), -inf < ceiling) -- a zero-density degenerate cohort is
      // safe, the opposite of the #550 overflow this guards.
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
                   util::to_string(xad::value(it->get_density())) + " (log_density=" +
                   util::to_string(xad::value(it->get_log_density())) + ", height=" +
                   util::to_string(xad::value(it->height())) + ") at time=" +
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
  const auto& env_vars = environment.vars;  // Internals_<value_type>
  for (size_t i = 0; i < env_vars.state_size; ++i) {
    if (!util::is_finite(env_vars.states[i])) {  // finiteness guard (Kind A)
      util::stop("Non-finite environment state (index " + util::to_string(i) +
                 " = " + util::to_string(xad::value(env_vars.states[i])) + ") at time=" +
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
typename Patch<T,E>::value_type
Patch<T,E>::compute_competition(double height) const {
  // Resident self-shading: the stand competition carries value_type so a trait
  // re-shades the environment. (At the double instantiation this is double, so
  // the R-facing signature is unchanged.)
  value_type tot = 0.0;
  for (size_t i = 0; i < species.size(); ++i) {
    if (!is_mutant_run) {
      tot += species[i].compute_competition(height) / area;
    }
  }
  return tot;
}

template <typename T, typename E>
std::vector<double> Patch<T,E>::r_compute_competition_effect_error_by_node_for_species_i(size_t species_index) const {
  const double tot_competition_effect = xad::value(compute_competition(0.0));
  return species[species_index].r_compute_competition_effect_by_nodes_error(tot_competition_effect);
}

// Integrate over lifetime fitness of individual nodes, scaled per node.
template <typename T, typename E>
typename Patch<T,E>::value_type Patch<T,E>::net_reproduction_ratio_for_species(
    size_t species_index, std::vector<double> const& scalars) const {
  auto net_prod = species[species_index].net_reproduction_ratio_by_node_weighted();
  auto const times = species[species_index].node_times();
  auto net_prod_scaled = std::vector<value_type>(times.size());
  for (size_t i = 0; i < times.size(); ++i) {
    net_prod_scaled[i] = net_prod[i] * scalars[i];
  }
  return util::trapezium(times, net_prod_scaled);
}

// Offspring production, equal to overall fitness scaled by the birth rate.
template <typename T, typename E>
std::vector<typename Patch<T,E>::value_type> Patch<T,E>::offspring_production() const {
  auto ret = std::vector<value_type>(species.size());
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
std::vector<typename Patch<T,E>::value_type> Patch<T,E>::net_reproduction_ratios() const {
  auto ret = std::vector<value_type>(species.size());
  for (size_t i = 0; i < species.size(); ++i) {
    auto scalars = std::vector<double>(species[i].size(), 1.0);
    ret[i] = net_reproduction_ratio_for_species(i, scalars);
  }
  return ret;
}

// Sum up all offspring produced.
template <typename T, typename E>
typename Patch<T,E>::value_type Patch<T,E>::total_offspring_production() const {
  value_type total = 0.0;
  std::vector<value_type> offspring = offspring_production();
  for (size_t i = 0; i < species.size(); ++i) {
    total += offspring[i];
  }
  return total;
}

// Check integration errors for each species' reproduction integral.
template <typename T, typename E>
std::vector<std::vector<double>> Patch<T,E>::net_reproduction_ratio_errors() const {
  std::vector<std::vector<double>> ret;
  // Schedule-refinement error is a double-only diagnostic; drop any active
  // derivative here with xad::value (identity on the double production path).
  double total_offspring = xad::value(total_offspring_production());
  for (size_t i = 0; i < species.size(); ++i) {
    auto weighted = species[i].net_reproduction_ratio_by_node_weighted();
    std::vector<double> weighted_d(weighted.size());
    for (size_t j = 0; j < weighted.size(); ++j) weighted_d[j] = xad::value(weighted[j]);
    ret.push_back(util::local_error_integration(
        species[i].node_times(), weighted_d, total_offspring));
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
  if (size() == 0 || is_mutant_run) return;

  // Transport-log-mass chart: reconstruct the log_density/density view from the
  // transported lambda before anything reads density (the competition field and
  // the fitted spline both do). No-op for strategies off the chart.
  for (auto& s : species) s.reconstruct_densities();

  // Assemble the exact separable field the cohorts read (deep-crown environments).
  if constexpr (env_has_competition_field<E>::value) {
    assemble_competition_field();
    // K93's field fully serves the read, so the fitted spline is dead weight and
    // is skipped (height_max() reads the species not the spline, the slope surface
    // has no callers, fixed environments build their own spline). FF16 keeps the
    // spline (field_supersedes_spline=false) for its not-yet-assembled reads.
    if constexpr (E::field_supersedes_spline) return;
  }

  // The environment-creation function returns value_type: the light spline is
  // built from active stand competition, so on a resident gradient pass the field
  // knots carry the self-shading derivative.
  auto f = [&](double x) -> value_type { return compute_competition(x); };
  environment.compute_environment(f, height_max(), rescale);
}

// Assemble the exact separable competition field from the whole cohort population.
// A(z) = sum_{H_j >= z} amp_j * M_j * Q(z/H_j), with amp_j = density_j * wpc_j (=
// node.compute_competition(0), since Q(0)=1), M_j the trapezium measure over the
// (descending) cohort heights, and Q's rank-3 factors from the canopy. Sources are
// all cohorts across species merged in descending height; the query factors carry
// the active query-height derivative. (Single shared canopy; per-species eta would
// need per-species fields -- noted, not exercised here.)
template <typename T, typename E>
void Patch<T,E>::assemble_competition_field() {
  // Build one source per cohort: source_weight[p] = amp * M * b_p(H). The
  // trapezium measure M and the source factors b_p use SAME-SPECIES neighbours
  // and the species' own canopy (compute_competition integrates each species
  // separately), so this is computed per species; the sources are then merged in
  // descending height for a single cumulative field.
  std::vector<value_type> H;                          // per-cohort height
  std::array<std::vector<value_type>, E::comp_rank> sw;  // per-cohort source weight
  for (size_t i = 0; i < species.size(); ++i) {
    const size_t m = species[i].size();
    if (m == 0) continue;
    std::vector<value_type> h(m), amp(m);
    size_t k = 0;
    // Copy the (small) canopy by value: r_get_strategy() returns a temporary, so
    // a pointer into it would dangle. Same shape for every cohort of the species.
    CanopyShape<value_type> canopy = species[i].node_begin()->individual.r_get_strategy().canopy_shape;
    for (auto it = species[i].node_begin(); it != species[i].node_end(); ++it, ++k) {
      h[k]   = it->height();
      // density*wpc (Q(0)=1), per patch area to match Patch::compute_competition.
      // On the chart a coincident cohort has density 0 (odelia's -inf spacing
      // convention), so amp = 0 there rather than a NaN from exp(lambda)/spacing.
      amp[k] = it->compute_competition(0.0) / area;
    }
    for (size_t j = 0; j < m; ++j) {
      // Same-species trapezium measure: half the span to the two neighbours
      // (one-sided at the crown top / smallest cohort). (Skips the new_node
      // boundary term -- a small, documented re-baseline of the field value.)
      const value_type hi = (j == 0)     ? h[0]     : h[j - 1];
      const value_type lo = (j == m - 1) ? h[m - 1] : h[j + 1];
      const value_type M = 0.5 * (hi - lo);
      const auto b = canopy.template shading_source_factors<value_type>(h[j]);
      H.push_back(h[j]);
      for (size_t p = 0; p < E::comp_rank; ++p) sw[p].push_back(amp[j] * M * b[p]);
    }
  }
  const size_t n = H.size();
  if (n == 0) { environment.clear_competition_field(); return; }

  // Merge all species' sources in descending height for the cumulative field.
  std::vector<size_t> ord(n);
  for (size_t j = 0; j < n; ++j) ord[j] = j;
  std::sort(ord.begin(), ord.end(), [&](size_t a, size_t b) {
    return odelia::util::to_passive(H[a]) > odelia::util::to_passive(H[b]);
  });
  std::array<std::vector<value_type>, E::comp_rank> sw_sorted;
  for (size_t p = 0; p < E::comp_rank; ++p) sw_sorted[p].resize(n);
  std::vector<double> heights_d(n);
  for (size_t j = 0; j < n; ++j) {
    heights_d[j] = odelia::util::to_passive(H[ord[j]]);
    for (size_t p = 0; p < E::comp_rank; ++p) sw_sorted[p][j] = sw[p][ord[j]];
  }
  // Query factors come from any cohort's canopy (shared shape); use species 0's.
  // Value copy: r_get_strategy() is a temporary, so a reference would dangle.
  CanopyShape<value_type> canopy = species[0].node_begin()->individual.r_get_strategy().canopy_shape;
  environment.assemble_competition_field(sw_sorted, heights_d, canopy);
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
    value_type resource_consumed = std::accumulate(species.begin(), species.end(), value_type(0.0), [i](value_type r, const species_type& s) {
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

// First set_ode_state function is for resident runs. Second is for mutant runs
template <typename T, typename E>
template <typename It>
It Patch<T,E>::set_ode_state(It it, double time) {

  // Set ode states
  it = odelia::ode::set_ode_state(species.begin(), species.end(), it);
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
template <typename It>
It Patch<T,E>::set_ode_state(It it, int index) {

  it = odelia::ode::set_ode_state(species.begin(), species.end(), it);

  // using a pointer here to avoid copying environment object
  // just point the pointer, used inside compute rates to get env, to relevant env object
  environment_ptr = &(environment_history[idx][index]);
  environment.time = environment_ptr->time;

  // increment the iterator by an appropriate amount, but don't actually do anything in the env
  for (size_t i = 0; i < environment_ptr->ode_size(); i++) {*it++;}
 
  compute_rates();
  return it;
}

// Per accepted ODE step: commit this step's stage snapshots and its time to the
// recording (resident pass only).
template <typename T, typename E>
void Patch<T,E>::record_ode_step() {
  if (recording) {
    step_history.push_back(time());
    environment_history.push_back(environment_cache);
  }
}

// Per RK stage: snapshot the environment into this step's stage buffer; stage 0
// starts a fresh buffer (resident pass only).
template <typename T, typename E>
void Patch<T,E>::record_stage(int stage) {
  if (recording) {
    if (stage == 0) {
      environment_cache.clear();
    }
    environment_cache.push_back(environment);
  }
}

// Per step on the replay pass: advance the recording cursor to this step's time
// so the mutant set_ode_state reads the right stage snapshots.
template <typename T, typename E>
void Patch<T,E>::replay_step() {
  if (has_recorded_field())
  {
    // Minor optimization to check the current and next index before doing a search, as the most common case is that the ODE solver is stepping through the cached environments in order. If the call sequence was not strictly sequential, we fallback to a search through the step history to find the correct environment.

    const double t = time();
    const size_t n = step_history.size();

    // Fast path: step_to() advances through ode_times in order, so this is
    // usually either the current cached step index or the next one.
    if (static_cast<size_t>(idx) < n && util::identical(step_history[idx], t)) {
      return;
    }
    if (static_cast<size_t>(idx + 1) < n &&
        util::identical(step_history[idx + 1], t)) {
      ++idx;
      return;
    }

    // Fallback to search if the call sequence was not strictly sequential.
    auto step = std::find(step_history.begin(), step_history.end(), t);
    if (step == step_history.end()) {
      util::stop("ODE time not found in step history");
    }
    idx = static_cast<int>(std::distance(step_history.begin(), step));
  }
}

template <typename T, typename E>
template <typename It>
It Patch<T,E>::ode_state(It it) const {
  it = odelia::ode::ode_state(species.begin(), species.end(), it);
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
template <typename It>
It Patch<T,E>::ode_rates(It it) const {
  it = odelia::ode::ode_rates(species.begin(), species.end(), it);
  it = environment.ode_rates(it);
  return it;
}

template <typename T, typename E>
template <typename It>
It Patch<T,E>::ode_aux(It it) const {
  it = odelia::ode::ode_aux(species.begin(), species.end(), it);
  return it;
}

}

#endif
