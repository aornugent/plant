// -*-c++-*-
#ifndef PLANT_PLANT_PATCH_H_
#define PLANT_PLANT_PATCH_H_

#include <plant/parameters.h>
#include <plant/species.h>
#include <plant/util.h>
#include <plant/adaptive_interpolator.h> // interpolator::refinement_failure
#include <odelia/ode_interface.hpp>
#include <odelia/gradient.hpp>

#include <plant/disturbance_regime.h>

#include <algorithm>
#include <cmath>
#include <limits>
#include <string>

using namespace Rcpp;

namespace plant {

// A strategy or environment named at scalar U. Reached through U so the lookup
// of the rebind alias waits until a rebind is actually asked for.
template <typename U>
struct at_scalar {
  template <typename X> using apply = typename X::template rebind<U>;
};

// One accepted step. The state widens at an introduction, so the record is ragged.
struct ode_step_record { double time; double step_size; std::vector<double> state; };

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
  typedef typename strategy_type::ptr strategy_type_ptr;

  Patch(parameters_type p, environment_type e, plant::Control c);
  void reset();

  // This patch at scalar U: strategies (already prepared), environment, node
  // structure, ODE state, and the birth stamps that divide the fecundity rate.
  template <class U,
            class T2 = typename at_scalar<U>::template apply<T>,
            class E2 = typename at_scalar<U>::template apply<E>>
  Patch<T2,E2> rebind_from() const;

  // Every scalar's Patch is one class, so a rebind reaches the rebound patch's
  // members.
  template <typename, typename> friend class Patch;
  size_t size() const {return species.size();}

  //Try using pointer in place of object itself
  double time() const {return environment.time;}
  double get_area() const { return area;}
  value_type height_max() const;

  value_type compute_competition(double height) const;

  // The competition profile and its vertical derivative at z, from one pass over
  // the species. The first entry equals compute_competition(z) bit for bit.
  std::pair<value_type, value_type>
  compute_competition_and_slope(double z) const;
  // The same pair for R, which takes only doubles: value first, then slope.
  std::vector<double> r_compute_competition_and_slope(double z) const {
    const std::pair<value_type, value_type> fs =
      compute_competition_and_slope(z);
    return {odelia::util::to_passive(fs.first),
            odelia::util::to_passive(fs.second)};
  }

  // Describe the size distribution around `height`, for error messages. The
  // light spline is built from the cohort heights, so when its refinement fails
  // the useful thing to report is what the cohorts are doing there.
  std::string describe_nodes_near(double height) const;

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

  // One entry per species gaining a node. This is what the solver carries back to
  // the patch when it walks a recording it widened, and it never reads it.
  using widening = std::vector<size_t>;

  // Widen the state by one node per species named, and bring what the width
  // changed up to date: the field, and the rates that describe neither once it
  // moves. narrow() is the exact inverse, dropping each species' NEWEST node --
  // so a sequence of widenings must be undone in reverse, which nothing in the
  // signature can say and only a round trip can check.
  void widen(const widening& species_index);
  void narrow(const widening& species_index);

  // The widening as a map: the state before it in, the whole widened state out,
  // and nothing rebuilt. This is the one the sweep transposes, so it runs at
  // whatever scalar it is called on and loads the state itself rather than
  // asking the caller to. It leaves this patch holding what it added.
  template <typename It>
  void widened_state(const widening& species_index, double time, It x,
                     std::vector<value_type>& y);

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
  // Compute rates of change at the state currently loaded and save them into the
  // ode solver. Computing here rather than in set_ode_state is what makes the
  // rates always those of that state, however the caller arrived at it.
  template <typename It> It ode_rates(It it);
  // Retrieve auxillary variables and save into the ode solver
  template <typename It> It ode_aux(It it) const;
  // Hand them back, in the order ode_aux wrote them
  template <typename It> It set_ode_aux(It it);

  // The mirror of ode_rates: the adjoints of dydt in, the adjoints of the state
  // out through the iterator ode_state would write to.
  template <class ItIn, class ItOut>
  ItOut ode_rates_adjoint(ItIn lambda_dydt, ItOut lambda_y);

  // The same transpose over several metrics at once. Everything that depends on
  // the STATE -- the two boundary evaluations, the field rebuild, and every
  // block recording -- runs once; only the seeding, the sweeps and the linear
  // maps run per metric. That is report 05's record-once-sweep-many, and here
  // the record is 99 per cent of the call.
  void ode_rates_adjoint_batched(
      const std::vector<std::vector<double>>& lambda_dydt,
      std::vector<std::vector<double>>& lambda_y);

  size_t node_count() const;

  // How often the boundary's own term was asked for, and how often it carried
  // anything. The boundary node stands at the seed's height for a whole run while
  // its condition is re-evaluated at every stage of every step, so this row acts
  // once per step where the quantity it belongs to changes once per plant. Its
  // value being right is therefore not the same claim as its being multiplied by
  // the right number of evaluations, and only the second is a count.
  //
  // The two agree now: one recording answers for every channel, so there is no
  // seed that leaves it with nothing to carry. They are kept apart because the
  // count is the claim and an early return would have to show up here.
  size_t boundary_condition_asked = 0;
  size_t boundary_condition_carried = 0;

  // The trait adjoints, species-major in each strategy's ad_parameters() order.
  // One trait is one input every cohort at every step reads, so this accumulates.
  //
  // One ROW per census metric being swept. A sweep records a block once and
  // sweeps it per metric, so the rows fill together and the accumulator has to
  // be as wide as the batch; a single-metric caller reads row zero and is
  // otherwise unaffected.
  std::vector<std::vector<double>> trait_adjoint;

  // Every species' differentiable parameters, species-major: the order a trait
  // row is indexed in, answered once rather than walked by each caller.
  std::vector<typename T::value_type*> ad_parameters();

  size_t trait_adjoint_size() const;
  // The same order, named. Each name carries its species index, because
  // concatenating the strategies' own names repeats every one of them per
  // species and character indexing then resolves each to species one's column,
  // which an unknown-name check cannot see.
  std::vector<std::string> trait_adjoint_names() const;
  void clear_trait_adjoint(size_t n_metrics = 1);

  // How large the stage's recording grew, and how many blocks the one tape has
  // carried. The stage holds every unit, so the size is linear in the unit count
  // rather than flat in it; what a leak looks like here is a slope that climbs.
  size_t block_recording_size = 0;
  size_t block_sweeps = 0;

  // The widening map's whole Jacobian, by forward tangent: one row per widened state entry
  // and one column per input, the state's entries first and the traits after.
  // Forming it entirely is what localises a disagreement to a cell, which a
  // contraction against the transpose cannot do.
  std::vector<std::vector<double>>
  introduction_jacobian(const std::vector<size_t>& species_index,
                        const std::vector<double>& state_before,
                        double time_before);

  // Returns state in structure format as opposed to single 
  // vector as given by ode_state
  Rcpp::List r_get_state() const;

  // Set state of patch, based on estimate of future state estimated by the solver.
  // The first is for resident runs and computes the environment as it goes.
  //
  // The second is for mutant runs, which read a cached environment by RK stage
  // rather than computing one. **The solver has no route to it.** It is reached
  // only through the stage-indexed derivs, which asks a System whether it has
  // recorded field values, and this one does not record them -- so the stage index
  // it writes is never written, and rate_environment() reads stage zero's
  // environment for all six stages of a mutant step. Finishing the mutant path
  // means recording the field, not calling this from somewhere else.
  template <typename It> It set_ode_state(It it, double time);
  template <typename It> It set_ode_state(It it, int index);

  // The state and the field: what set_ode_state establishes before it computes
  // rates. ode_rates_adjoint is taken here, where the block recordings carry the
  // rate chain and a rate evaluation in double would repeat all of it.
  template <typename It> It set_ode_state_and_field(It it, double time);

  // A recorded state loaded as the run itself carries it. set_ode_state_and_field
  // evaluates the inflow condition in the field that leaves the boundary interval
  // off, then rebuilds the field including it; the run then rates the nodes and
  // evaluates the condition a second time, in that second field. It is the second
  // value an introduced node inherits and the census reads, so reloading a state
  // without it linearises a boundary node the trajectory never carried.
  template <typename It> It set_recorded_state(It it, double time);

  // The inflow condition alone, in the field as it now stands. Public because the
  // two evaluations above have to be taken one at a time to be told apart.
  void compute_boundary_nodes();

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
    return species[species_index.check_bounds(size())];
  }
  // This is only here because it wraps a private function.
  void r_compute_environment() {compute_environment(false);}

  // env. cache for assembly
  std::vector<double> step_history{0.0};  // always start at zero
  std::vector<std::vector<environment_type>> environment_history;
  std::vector<environment_type> environment_cache;

  void cache_ode_step();
  void cache_RK45_step(int step);
  void load_ode_step();

  // Per-accepted-step recording. Called by the stepper once a step is committed, so
  // the states it keeps are the run's own. The field is not kept with them: it is
  // rebuilt at the recorded positions, which is what lets its derivative flow.
  void record_ode_step();

  bool record_steps = false;
  // One record per accepted step, the first being the initial state. The step size is
  // the solver's, filled in beside its own state once the run is over: the stepper
  // does not pass it here, and it cannot be differenced out of the times.
  std::vector<ode_step_record> trajectory;
  
  // used cache_ode_step for mutant runs
  bool save_RK45_cache;

  // used in load_ode_step for mutant runs
  bool use_cached_environment = false;

  bool is_mutant_run = false;

  void set_mutant();
  void add_strategies(std::vector<strategy_type> strategies);
  void overwrite_strategies(std::vector<strategy_type> strategies);

private:
  // A patch whose species take the prepared strategies given, rather than
  // preparing the ones in the parameters. rebind_from's only route in.
  Patch(parameters_type p, environment_type e, plant::Control c,
        const std::vector<strategy_type_ptr>& prepared);

  // Step (b)'s tape and templates, held so one tape spans the cohort loop.
  // Type-erased: a model declaring no rebind has no active type to form here.
  // What each species' reduction had accumulated at each knot before its
  // closing trapezium, kept from the field built without the boundary interval
  // so the field built with it costs one trapezium per knot rather than a second
  // walk over every node.
  mutable std::vector<std::vector<typename species_type::competition_split>>
    competition_capture;
  mutable size_t capture_at = 0;
  void compute_environment_excl_capturing(bool rescale);
  void compute_environment_closing(bool rescale);

  int idx = 0; // used to access environment cache for mutant runs
  void compute_environment(bool rescale);
  // One field build, with every species' inflow boundary interval included or not.
  void compute_environment_once(bool rescale, bool include_boundary);
  // The competition profile with every species' boundary interval left off: a
  // function of the ODE state alone.
  value_type compute_competition_excl_boundary(double height) const;
  // That same profile and its vertical derivative, from one pass over the species.
  std::pair<value_type, value_type>
  compute_competition_and_slope_excl_boundary(double height) const;
  void compute_rates();

  // The environment the rates are computed against: the patch's own on a
  // resident run, and on a mutant step the one recorded for that step, which the
  // mutant experiences rather than shapes. Derived from what the patch owns, so
  // it still refers into the patch after the patch is copied.
  environment_type& rate_environment() {
    return use_cached_environment
      ? environment_history.at(idx).at(cached_environment_index)
      : environment;
  }

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
  // Guard the birth-date coordinate's quadrature grid: it is the per-node
  // introduction times, so nodes sharing one give zero-width intervals that drop
  // silently out of the integral. Only reachable for a patch whose nodes were
  // created outside the schedule (seeded or imported without per-node times).
  void check_birth_dates_distinct() const;

  parameters_type parameters;

  double area;
  environment_type environment;
  std::vector<species_type> species;

  //TODO(#476): Move into environment?
  std::vector<value_type> resource_depletion;

  // Which recorded environment a mutant step is evaluated against; see
  // rate_environment(). Unused on a resident run.
  int cached_environment_index = 0;

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
  // The validated member, not the argument: validate() derives the regime from
  // patch_type and max_patch_lifetime, so an argument whose lifetime was
  // assigned after its own construction still carries the regime of the lifetime
  // it was constructed with.
  survival_weighting = parameters.disturbance;

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
Patch<T,E>::Patch(parameters_type p, environment_type e, Control c,
                  const std::vector<strategy_type_ptr>& prepared)
  : parameters(p),
    area(p.patch_area),
    environment(e),
    control(c),
    environment_cache(6) {

  parameters.validate();

  save_RK45_cache = control.save_RK45_cache;
  survival_weighting = parameters.disturbance;

  environment.set_shading_model(control.shading_model,
                                control.ppa_layer_optical_depth,
                                control.ppa_layer_smoothing);

  for (const strategy_type_ptr& s : prepared) {
    species.push_back(Species<T,E>(s));
  }

  reset();
}

template <typename T, typename E>
template <class U, class T2, class E2>
Patch<T2,E2> Patch<T,E>::rebind_from() const {
  Parameters<T2,E2> p2;
  p2.patch_area = parameters.patch_area;
  p2.n_patches = parameters.n_patches;
  p2.patch_type = parameters.patch_type;
  p2.max_patch_lifetime = parameters.max_patch_lifetime;
  p2.node_schedule_times_default = parameters.node_schedule_times_default;
  p2.node_schedule_times = parameters.node_schedule_times;
  p2.ode_times = parameters.ode_times;
  p2.initial_time = parameters.initial_time;
  p2.strategy_default = parameters.strategy_default.template rebind_from<U>();

  // The strategies the species run, not the ones in the parameters: those are
  // the prepared copies, and preparing again is refused at an active scalar.
  std::vector<typename T2::ptr> prepared;
  for (const species_type& s : species) {
    prepared.push_back(std::make_shared<T2>(
      s.strategy_ptr()->template rebind_from<U>()));
    p2.strategies.push_back(*prepared.back());
  }

  E2 env = environment.template rebind_from<U>();
  Patch<T2,E2> out(p2, env, control, prepared);

  for (size_t i = 0; i < species.size(); ++i) {
    for (size_t j = 0; j < species[i].size(); ++j) {
      out.species[i].introduce_new_node();
    }
    // Not in the ODE state, and pr_patch_survival_at_birth divides the
    // fecundity rate: without these the rebound rates differ.
    out.species[i].set_birth_state(species[i].node_times(),
                                   species[i].r_patch_densities(),
                                   species[i].r_pr_patch_survival_at_birth());
  }

  std::vector<U> node_state(node_ode_size());
  odelia::ode::ode_state(species.begin(), species.end(), node_state.begin());
  odelia::ode::set_ode_state(out.species.begin(), out.species.end(),
                             node_state.begin());

  // reset() in the constructor cleared the environment back to its initial
  // soil state, so restore the current one before the spline is rebuilt.
  out.environment = env;
  out.environment.set_shading_model(control.shading_model,
                                    control.ppa_layer_optical_depth,
                                    control.ppa_layer_smoothing);
  // The field and the boundary node are left to the caller. Every caller sets a
  // state through set_ode_state_and_field or set_recorded_state before reading
  // either, and both rebuild the field and re-evaluate the inflow condition, so
  // computing them here solves the boundary leaf twice per right-hand side and
  // then discards it. A caller that reads before setting gets an unbuilt field.
  return out;
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
   for (auto& s : species) {
    s.clear();
    // allocate variables for tracking resource consumption
    s.resize_consumption_rates(environment.n_resources());
  }

  // resize to species count
  resource_depletion.reserve(environment.n_resources());

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
    compute_rates();
  }

  // clear accumulated per-node competition error
  competition_error_by_node.assign(species.size(), {});

  trajectory.clear();
  // No step reached the initial time, so it records no size.
  record_ode_step();
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
  auto it = parameters.initial_state.begin();
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

  check_birth_dates_distinct();

  // Build the environment from the real node heights (full recompute, no
  // rescale) and compute rates for the seeded population.
  compute_environment(false);
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
void Patch<T,E>::check_birth_dates_distinct() const {
  for (size_t i = 0; i < species.size(); ++i) {
    // Only the birth-date coordinate integrates over these times. On the height
    // path they feed the lifetime-fitness integral after the run, where a
    // repeated time has always been tolerated.
    if (!species[i].density_in_birth_date()) {
      continue;
    }
    if (!species[i].birth_dates_are_distinct()) {
      util::stop("Species " + util::to_string(i + 1) + " has nodes sharing an "
                 "introduction time, which the birth-date size-density "
                 "coordinate integrates over: the repeated nodes span zero width "
                 "and drop out of the competition and resource integrals. Supply "
                 "per-node introduction times (parameters$initial_node_times) "
                 "with the initial state, or run with "
                 "control$node_density_in_birth_date = FALSE.");
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

  // Mode (2): a non-finite environment state, read through the ODE interface so
  // an environment with no integrated state contributes an empty loop rather
  // than a special case. For TF24 these are the per-depth soil-water states.
  std::vector<value_type> env_state(environment.ode_size());
  environment.ode_state(env_state.begin());
  for (size_t i = 0; i < env_state.size(); ++i) {
    const double state_i = odelia::util::to_passive(env_state[i]);
    if (!util::is_finite(state_i)) {
      util::stop("Non-finite environment state (index " + util::to_string(i) +
                 " = " + util::to_string(state_i) + ") at time=" +
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
typename Patch<T,E>::value_type Patch<T,E>::height_max() const {
  value_type ret = 0.0;
  for (size_t i = 0; i < species.size(); ++i) {
    if (!is_mutant_run) {
      const value_type h = species[i].height_max();
      if (h > ret) {
        ret = h;
      }
    }
  }
  return ret;
}

template <typename T, typename E>
typename Patch<T,E>::value_type
Patch<T,E>::compute_competition(double height) const {
  value_type tot = 0.0;
  for (size_t i = 0; i < species.size(); ++i) {
    if (!is_mutant_run) {
      tot += species[i].compute_competition(height) / area;
    }
  }
  return tot;
}

template <typename T, typename E>
std::pair<typename Patch<T,E>::value_type, typename Patch<T,E>::value_type>
Patch<T,E>::compute_competition_and_slope(double z) const {
  value_type tot = 0.0, tot_slope = 0.0;
  for (size_t i = 0; i < species.size(); ++i) {
    if (!is_mutant_run) {
      const std::pair<value_type, value_type> fs =
        species[i].compute_competition_and_slope(z);
      tot       += fs.first / area;
      tot_slope += fs.second / area;
    }
  }
  return {tot, tot_slope};
}

template <typename T, typename E>
std::vector<double> Patch<T,E>::r_compute_competition_effect_error_by_node_for_species_i(size_t species_index) const {
  const double tot_competition_effect =
    odelia::util::to_passive(compute_competition(0.0));
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
// The competition profile with every species' inflow boundary interval left off.
template <typename T, typename E>
typename Patch<T,E>::value_type
Patch<T,E>::compute_competition_excl_boundary(double height) const {
  value_type tot = 0.0;
  for (size_t i = 0; i < size(); ++i) {
    tot += species[i].compute_competition_excl_boundary(height) / area;
  }
  return tot;
}

template <typename T, typename E>
std::pair<typename Patch<T,E>::value_type, typename Patch<T,E>::value_type>
Patch<T,E>::compute_competition_and_slope_excl_boundary(double height) const {
  value_type tot = 0.0, tot_slope = 0.0;
  for (size_t i = 0; i < size(); ++i) {
    const std::pair<value_type, value_type> fs =
      species[i].compute_competition_and_slope_excl_boundary(height);
    tot       += fs.first / area;
    tot_slope += fs.second / area;
  }
  return {tot, tot_slope};
}

// Evaluate every species' inflow boundary condition in the field as it currently
// stands. Owned by the field build rather than by compute_rates(), so that the
// field reads a boundary density derived from this state instead of one carried
// from the previous evaluation.
// A mutant experiences a recorded environment rather than shaping one, and this
// reads the patch's own, so it is not the mutant's condition to evaluate.
template <typename T, typename E>
void Patch<T,E>::compute_boundary_nodes() {
  if (is_mutant_run) {
    return;
  }
  const double time_ = environment.time;
  const double pr_patch_survival = survival_weighting->pr_survival(time_);
  for (size_t i = 0; i < size(); ++i) {
    const double birth_rate =
      species[i].extrinsic_drivers().evaluate("birth_rate", time_);
    species[i].compute_boundary_node(environment, pr_patch_survival, birth_rate);
  }
}

// Creates splines of resource availability.
//
// The reduction's closing trapezium is the inflow boundary condition
// n_b = birth_rate * pr_estab / g, and n_b needs a field to be evaluated in --
// so the field and the boundary condition are mutually dependent. Ordering the
// build removes the cycle rather than iterating it:
//
//   A0  the reduction excluding the boundary interval   (the ODE state alone)
//   n_b the boundary condition evaluated in A0          (the ODE state alone)
//   A   A0 plus the boundary interval formed from n_b   (the ODE state alone)
//
// so a stage is a function of (y, t) and nothing else. Closing the fixed point by
// iteration instead only attenuates the carried dependence by the contraction
// modulus, which is ~1e-3.
template <typename T, typename E>
void Patch<T,E>::compute_environment(bool rescale) {
  if (!(size() > 0 && !is_mutant_run)) {
    return;
  }
  compute_environment_excl_capturing(rescale);
  compute_boundary_nodes();
  compute_environment_closing(rescale);
}

// The field without the boundary interval, keeping each species' reduction at
// the point its closing trapezium would be added.
template <typename T, typename E>
void Patch<T,E>::compute_environment_excl_capturing(bool rescale) {
  for (auto& s : species) {
    s.set_new_node_birth_date(environment.time);
  }
  competition_capture.assign(size(), {});
  capture_at = 0;
  auto f = [&](double x) -> std::pair<value_type, value_type> {
    value_type tot = 0.0, tot_slope = 0.0;
    for (size_t i = 0; i < size(); ++i) {
      const typename species_type::competition_split c =
        species[i].compute_competition_and_slope_split(x);
      competition_capture[i].push_back(c);
      tot       += c.excl.first / area;
      tot_slope += c.excl.second / area;
    }
    return {tot, tot_slope};
  };
  if (size() > 0 & !is_mutant_run) {
    try {
      environment.compute_environment(f, height_max(), rescale);
    } catch (const interpolator::refinement_failure& e) {
      util::stop(std::string(e.what()) + " " + describe_nodes_near(e.report.x_lo));
    }
  }
}

// The same field with the boundary interval closed, from the kept reductions and
// the boundary node the call between the two established.
template <typename T, typename E>
void Patch<T,E>::compute_environment_closing(bool rescale) {
  for (auto& s : species) {
    s.set_new_node_birth_date(environment.time);
  }
  capture_at = 0;
  auto f = [&](double x) -> std::pair<value_type, value_type> {
    value_type tot = 0.0, tot_slope = 0.0;
    for (size_t i = 0; i < species.size(); ++i) {
      if (!is_mutant_run) {
        if (capture_at >= competition_capture[i].size()) {
          util::stop("compute_environment_closing: the field without the "
                     "boundary interval was built over a different knot set");
        }
        const std::pair<value_type, value_type> fs =
          species[i].close_competition_and_slope(competition_capture[i][capture_at], x);
        tot       += fs.first / area;
        tot_slope += fs.second / area;
      }
    }
    ++capture_at;
    return {tot, tot_slope};
  };
  if (size() > 0 & !is_mutant_run) {
    try {
      environment.compute_environment(f, height_max(), rescale);
    } catch (const interpolator::refinement_failure& e) {
      util::stop(std::string(e.what()) + " " + describe_nodes_near(e.report.x_lo));
    }
  }
}

template <typename T, typename E>
void Patch<T,E>::compute_environment_once(bool rescale, bool include_boundary) {

  // The boundary node is one end of the birth-date quadrature and its birth date
  // is the current time, so refresh it before the profile is built. Its own
  // stamp is set in compute_rates(), which the stepper calls *after* the
  // set_ode_state() that brings us here, so reading that stamp would use the
  // previous derivs call's time and shorten the boundary interval by a
  // Runge-Kutta stage. Measured effect on FF16 offspring production is below
  // 1e-6 -- the boundary node carries almost no leaf area, so the segment it
  // ends contributes little -- but the interval is then a function of the step
  // size, which the spatial quadrature has no business depending on, and the
  // whole integral *is* that one segment while a species has a single node.
  // No-op for the height coordinate, where this abscissa is the constant
  // initial height.
  for (auto& s : species) {
    s.set_new_node_birth_date(environment.time);
  }

  // The competition profile and its vertical derivative at x. The field carries a
  // slope at every knot, so the build asks for the pair.
  // Written as std::pair<double, double> this still compiles, taking the value
  // of an active profile, and every knot value and slope in the field would then
  // be a constant with nothing raised to say the cohorts had stopped reaching it.
  auto f = [&](double x) -> std::pair<value_type, value_type> {
    return include_boundary ? compute_competition_and_slope(x)
                            : compute_competition_and_slope_excl_boundary(x);
  };

  if (size() > 0 & !is_mutant_run) {
    try {
      environment.compute_environment(f, height_max(), rescale);
    } catch (const interpolator::refinement_failure& e) {
      // The refiner can only say that the profile has a feature it cannot
      // resolve, at some height. We know what is at that height, and it is
      // usually the actual problem: cohorts stacked at one size put a step in
      // the competition profile (#571). Say so rather than making the reader
      // reconstruct the patch state by hand.
      util::stop(std::string(e.what()) + " " +
                 describe_nodes_near(e.report.x_lo));
    }
  }
}

// Report what the size distribution is doing around `height`: how many cohorts
// sit within a narrow window of it, and -- the usual culprit -- whether the node
// list is still ordered by decreasing height. Species::compute_competition() and
// Species::height_max() both take that ordering as given (see the invariant noted
// on both), so once it is violated the competition profile they build has
// fictitious steps in it, and the refiner fails on one of them (#571).
//
// A window rather than a single point because the unresolved interval is ~1e-5 m
// wide, while what matters is whether cohorts share effectively the same size.
template <typename T, typename E>
std::string Patch<T,E>::describe_nodes_near(double height) const {
  const double window = 1e-3; // m

  std::string ret = "Patch state at that height (time " +
                    util::format_double(environment.time) + "):";

  for (size_t i = 0; i < species.size(); ++i) {
    size_t n_near = 0, n_inversions = 0, n_inversions_live = 0, n_zero_density = 0;
    double h_near_min = std::numeric_limits<double>::infinity(),
           h_near_max = -std::numeric_limits<double>::infinity(),
           h_max = -std::numeric_limits<double>::infinity(),
           h_front = NA_REAL, h_prev = NA_REAL;
    bool prev_live = false;

    for (auto it = species[i].node_begin(); it != species[i].node_end(); ++it) {
      // Everything in this scan is counted, compared or formatted into the
      // message, so the height is read at its value.
      const double h = odelia::util::to_passive(it->height());
      const bool live = it->get_density() > 0.0;
      if (!live) {
        n_zero_density++;
      }
      if (!util::is_finite(h_front)) {
        h_front = h;
      } else if (h > h_prev) {
        n_inversions++;
        // Whether the *live* cohorts are still ordered is the question that
        // matters: the method of characteristics guarantees it for them (growth
        // trajectories sharing an environment cannot cross), so inversions among
        // zero-density nodes are bookkeeping debris in the quadrature grid,
        // whereas an inversion between two live cohorts would mean the
        // characteristics themselves had crossed.
        if (live && prev_live) {
          n_inversions_live++;
        }
      }
      h_prev = h;
      prev_live = live;
      h_max = std::max(h_max, h);
      if (fabs(h - height) <= window) {
        n_near++;
        h_near_min = std::min(h_near_min, h);
        h_near_max = std::max(h_near_max, h);
      }
    }

    ret += " species " + util::to_string(i + 1) + ": " +
           util::to_string(n_near) + " of " +
           util::to_string(species[i].size()) + " cohorts within " +
           util::format_double(window) + " m";
    if (n_near > 0) {
      ret += ", spanning [" + util::format_double(h_near_min) + ", " +
             util::format_double(h_near_max) + "]";
      // Several cohorts at one size means growth has stalled and the schedule is
      // introducing new cohorts into a patch with nothing to separate them.
      if (n_near > 1 && (h_near_max - h_near_min) < window / 100) {
        ret += " -- effectively a single size, so the size distribution has"
               " collapsed here";
      }
    }
    if (n_inversions > 0) {
      ret += "; node heights are NOT decreasing (" +
             util::to_string(n_inversions) + " inversions, " +
             util::to_string(n_inversions_live) +
             " of them between two cohorts of non-zero density; " +
             util::to_string(n_zero_density) + " of " +
             util::to_string(species[i].size()) +
             " nodes have zero density), which breaks the ordering that"
             " Species::compute_competition() and height_max() assume: the"
             " tallest cohort is " + util::format_double(h_max) +
             " but height_max() reports " + util::format_double(h_front) +
             " (the front node), so the competition profile is wrong and its"
             " steps are an artefact rather than a feature of the model";
      if (n_inversions_live == 0) {
        ret += " (the live cohorts are still correctly ordered, so this is the"
               " quadrature grid being scrambled by zero-density nodes rather"
               " than growth trajectories crossing)";
      }
    }
    ret += ";";
  }

  return ret;
}


template <typename T, typename E>
void Patch<T,E>::compute_rates() {

  // Computes rates of change for the patch, including all the component species,
  // against the environment the patch experiences (see rate_environment()).
  environment_type& env = rate_environment();
  double time_ = env.time;

  double pr_patch_survival = survival_weighting->pr_survival(time_);
  for (size_t i = 0; i < size(); ++i) {
    double birth_rate = species[i].extrinsic_drivers().evaluate("birth_rate", time_);

    species[i].compute_rates(env, pr_patch_survival, birth_rate);
  }

  resource_depletion.reserve(env.n_resources());
  for(size_t i = 0; i < env.n_resources(); i++) {
    value_type resource_consumed = std::accumulate(species.begin(), species.end(), value_type(0.0), [i](const value_type& r, const species_type& s) -> value_type {
      return r + s.consumption_rate(i); // accumulates r from zero
    });

    resource_depletion.push_back(resource_consumed / area);
  }
  

  env.compute_rates(resource_depletion);

  //todo do we need to clear this every step?
  resource_depletion.clear();

}

// The whole light environment is rebuilt here, where only the knots below the
// seedling's height change. The knot fractions are fixed, so a narrower rebuild
// would write a subrange of the same positions.
template <typename T, typename E>
void Patch<T,E>::introduce_new_node(size_t species_index) {
  
  species[species_index].introduce_new_node();

  compute_environment(false);
}

template <typename T, typename E>
void Patch<T,E>::widen(const widening& species_index) {
  // Record introduction time and patch-age density on each node as it is
  // introduced, so lifetime-fitness calcs need not look these up later.
  const double t = time();
  const double patch_density = survival_weighting->density(t);
  for (size_t i : species_index) {
    species[i].introduce_new_node(t, patch_density);
  }

  compute_environment(false);

  // New nodes have just changed the state and the light field, so the stored
  // rates now describe neither. The solver reads them next without checking, so
  // they have to be brought up to date here.
  compute_rates();
}

template <typename T, typename E>
void Patch<T,E>::narrow(const widening& species_index) {
  for (size_t i : species_index) {
    species[i].remove_newest_node();
  }
  compute_environment(false);
  compute_rates();
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
  // Every node here is a copy of the boundary node, so they all carry the same
  // birth date and the ODE state does not restore it (it is bookkeeping, not a
  // state variable). Fine for the height coordinate; fatal for the birth-date
  // one, which uses these times as its quadrature grid.
  check_birth_dates_distinct();
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
  return odelia::ode::aux_size(species.begin(), species.end()) +
    environment.aux_size();
}

template <typename T, typename E>
double Patch<T,E>::ode_time() const {
  return time();
}

// First set_ode_state function is for resident runs. Second is for mutant runs
template <typename T, typename E>
template <typename It>
It Patch<T,E>::set_ode_state_and_field(It it, double time) {

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

  return it;
}

template <typename T, typename E>
template <typename It>
It Patch<T,E>::set_ode_state(It it, double time) {
  it = set_ode_state_and_field(it, time);

  return it;
}

// The second evaluation of the inflow condition, in the field the first one was
// folded into. See the declaration for why a reloaded state needs it.
template <typename T, typename E>
template <typename It>
It Patch<T,E>::set_recorded_state(It it, double time) {
  it = set_ode_state_and_field(it, time);
  compute_boundary_nodes();
  return it;
}

// used for mutant runs
// -- differs from above in that an index is passed in as argument
// -- environments are loaded from ODE history, instead of being calculated 
template <typename T, typename E>
template <typename It>
It Patch<T,E>::set_ode_state(It it, int index) {

  it = odelia::ode::set_ode_state(species.begin(), species.end(), it);

  // Record which of this step's cached environments the rates are evaluated
  // against; rate_environment() reads it back without copying the object.
  cached_environment_index = index;
  const environment_type& cached = rate_environment();
  environment.time = cached.time;

  // increment the iterator by an appropriate amount, but don't actually do anything in the env
  for (size_t i = 0; i < cached.ode_size(); i++) {*it++;}

  return it;
}

// called from ode_solver->cache
// saves cached set of environments(6) from each ODE step to the step history
template <typename T, typename E>
void Patch<T,E>::cache_ode_step() {
  if(save_RK45_cache) { 
    step_history.push_back(time());
    environment_history.push_back(environment_cache);
  }
}

// called from ode_step->cache
// saves environment at each RK45 step to the environment cache
template <typename T, typename E>
void Patch<T,E>::cache_RK45_step(int step) {
  if(save_RK45_cache) {  
    if(step == 0) {
      environment_cache.clear();
    }
    environment_cache.push_back(environment);
  }
}

// called from ode_solver->load, only gets called for mutant runs
template <typename T, typename E>
void Patch<T,E>::load_ode_step() {
  if (use_cached_environment)
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
void Patch<T,E>::record_ode_step() {
  if (!record_steps) {
    return;
  }
  ode_step_record record{time(), std::numeric_limits<double>::quiet_NaN(),
                         std::vector<double>(ode_size())};
  ode_state(record.state.begin());
  trajectory.push_back(std::move(record));
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
It Patch<T,E>::ode_rates(It it) {
  compute_rates();
  it = odelia::ode::ode_rates(species.begin(), species.end(), it);
  it = environment.ode_rates(it);
  return it;
}

template <typename T, typename E>
template <typename It>
It Patch<T,E>::ode_aux(It it) const {
  it = odelia::ode::ode_aux(species.begin(), species.end(), it);
  it = environment.ode_aux(it);
  return it;
}

template <typename T, typename E>
template <typename It>
It Patch<T,E>::set_ode_aux(It it) {
  it = odelia::ode::set_ode_aux(species.begin(), species.end(), it);
  it = environment.set_ode_aux(it);
  return it;
}

template <typename T, typename E>
size_t Patch<T,E>::node_count() const {
  size_t n = 0;
  for (size_t i = 0; i < species.size(); ++i) {
    n += species[i].size();
  }
  return n;
}

// The differentiable parameters of every species, species-major, which is the
// order every row indexed by a trait is in. A patch answers for this rather than
// each caller walking the species, because the order is the layout of a gradient
// and a caller that walks it itself is free to walk it differently.
template <typename T, typename E>
std::vector<typename T::value_type*> Patch<T,E>::ad_parameters() {
  std::vector<typename T::value_type*> ret;
  for (size_t i = 0; i < species.size(); ++i) {
    for (typename T::value_type* p :
         species[i].strategy_ptr()->ad_parameters()) {
      ret.push_back(p);
    }
  }
  return ret;
}

template <typename T, typename E>
size_t Patch<T,E>::trait_adjoint_size() const {
  size_t n = 0;
  for (size_t i = 0; i < species.size(); ++i) {
    n += species[i].strategy_ptr()->ad_parameters().size();
  }
  return n;
}

template <typename T, typename E>
std::vector<std::string> Patch<T,E>::trait_adjoint_names() const {
  std::vector<std::string> ret;
  ret.reserve(trait_adjoint_size());
  for (size_t i = 0; i < species.size(); ++i) {
    const std::string species_index =
      util::to_string(static_cast<int>(i + 1)) + ".";
    for (const std::string& n :
         species[i].strategy_ptr()->ad_parameter_names()) {
      ret.push_back(species_index + n);
    }
  }
  return ret;
}

template <typename T, typename E>
void Patch<T,E>::clear_trait_adjoint(size_t n_metrics) {
  if (n_metrics == 0) {
    util::stop("clear_trait_adjoint: a sweep accumulates for at least one metric");
  }
  trait_adjoint.assign(n_metrics, std::vector<double>(trait_adjoint_size(), 0.0));
}

template <typename T, typename E>
template <typename It>
void Patch<T,E>::widened_state(const widening& species_index, double time, It x,
                               std::vector<value_type>& y) {
  set_recorded_state(x, time);
  for (size_t i : species_index) {
    species[i].introduce_new_node();
  }
  util::check_length(y.size(), ode_size());
  ode_state(y.begin());
}

// One tangent seed per input column, over the same map. The node the map pushes is
// removed between columns, so every column is taken from the same width.
template <typename T, typename E>
std::vector<std::vector<double>>
Patch<T,E>::introduction_jacobian(const std::vector<size_t>& species_index,
                                  const std::vector<double>& state_before,
                                  double time_before) {
  using tangent = xad::fwd<double>::active_type;
  const size_t n_state = ode_size();
  const size_t n_trait = trait_adjoint_size();
  const size_t n_out = n_state + species_index.size() * node_type::ode_size();
  util::check_length(state_before.size(), n_state);

  std::vector<double> in(state_before);
  in.reserve(n_state + n_trait);
  for (const typename T::value_type* p : ad_parameters()) {
    in.push_back(odelia::util::to_passive(*p));
  }
  util::check_length(in.size(), n_state + n_trait);

  auto active = this->template rebind_from<tangent>();
  std::vector<std::vector<double>> ret(n_out,
                                       std::vector<double>(in.size(), 0.0));
  for (size_t c = 0; c < in.size(); ++c) {
    std::vector<tangent> x(in.size());
    for (size_t j = 0; j < in.size(); ++j) {
      x[j] = in[j];
      xad::derivative(x[j]) = (j == c) ? 1.0 : 0.0;
    }
    // The traits first, for the reason the transpose writes them first: a
    // quantity the state determines reads them while deriving it.
    size_t at = n_state;
    for (tangent* p : active.ad_parameters()) {
      *p = x[at++];
    }
    util::check_length(at, x.size());
    std::vector<tangent> y(n_out);
    active.widened_state(species_index, time_before, x.begin(), y);
    for (size_t r = 0; r < n_out; ++r) {
      ret[r][c] = xad::derivative(y[r]);
    }
    for (size_t i : species_index) {
      active.species[i].remove_newest_node();
    }
  }
  return ret;
}


template <typename T, typename E>
void Patch<T,E>::ode_rates_adjoint_batched(
    const std::vector<std::vector<double>>& lambda_dydt,
    std::vector<std::vector<double>>& lambda_y) {
  using scalar = odelia::ode::active_scalar<double>;
  using active_strategy = typename at_scalar<scalar>::template apply<T>;
  const size_t n_state = ode_size();
  const size_t n_seed = lambda_dydt.size();
  if (n_seed == 0) {
    util::stop("ode_rates_adjoint: needs at least one rate adjoint");
  }
  for (size_t m = 0; m < n_seed; ++m) {
    util::check_length(lambda_dydt[m].size(), n_state);
  }
  if (trait_adjoint.size() < n_seed) {
    clear_trait_adjoint(n_seed);
  }

  std::vector<value_type> current(n_state);
  ode_state(current.begin());
  std::vector<double> state(n_state);
  for (size_t j = 0; j < n_state; ++j) {
    state[j] = odelia::util::to_passive(current[j]);
  }

  // Built with no tape active, so it holds no slot the recording could alias.
  Patch<active_strategy, typename at_scalar<scalar>::template apply<E>> active =
    this->template rebind_from<scalar>();

  const double time_ = environment.time;
  // The stage, recorded whole. Everything between the state and the rates is an
  // intermediate of it -- the field and both its reductions, the inflow
  // condition, every individual's physiology, the water aggregation and the soil
  // -- so none of them needs a transpose written for it, and none can drift from
  // the forward it transposes because there is one function here.
  auto rates = [&](typename std::vector<scalar>::const_iterator x,
                   std::vector<scalar>& y) -> void {
    active.set_ode_state_and_field(x, time_);
    active.ode_rates(y.begin());
  };

  // The inflow condition is inside this recording, so its row enters once per
  // stage per step per metric and this is still where that count is taken. It
  // is no longer skippable: one recording answers for every channel, so there is
  // no seed that leaves it with nothing to carry.
  boundary_condition_asked += n_seed;
  boundary_condition_carried += n_seed;

  typename scalar::tape_type tape(false);
  block_recording_size = odelia::ode::state_and_parameter_adjoints(
    tape, active, state, lambda_dydt, rates, lambda_y, trait_adjoint);
  block_sweeps = n_seed;
}

// One metric's transpose, accumulating into this patch's own trait adjoint. The
// batched form is the real one; every caller outside the sweep takes this.
template <typename T, typename E>
template <class ItIn, class ItOut>
ItOut Patch<T,E>::ode_rates_adjoint(ItIn lambda_dydt, ItOut lambda_y) {
  const size_t n = ode_size();
  std::vector<std::vector<double>> lambda_in(1, std::vector<double>(n));
  for (size_t i = 0; i < n; ++i) {
    lambda_in[0][i] = *lambda_dydt++;
  }
  if (trait_adjoint.size() != 1 ||
      trait_adjoint[0].size() != trait_adjoint_size()) {
    clear_trait_adjoint();
  }
  std::vector<std::vector<double>> out;
  ode_rates_adjoint_batched(lambda_in, out);
  for (size_t i = 0; i < n; ++i) {
    *lambda_y++ = out[0][i];
  }
  return lambda_y;
}

}

#endif
