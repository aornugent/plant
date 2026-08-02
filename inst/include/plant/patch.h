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
  void introduce_new_nodes(const std::vector<size_t>& species_index);
  // The mirror of introduce_new_nodes: drop each species' newest node and
  // rebuild the field and the rates at the narrower width.
  void remove_new_nodes(const std::vector<size_t>& species_index);

  // The species whose newest node was introduced at `time`, in species order.
  // Read off the nodes' own birth times, which is what the run did; the schedule
  // is what it was asked to do, and it is consumed by the time a sweep runs.
  std::vector<size_t> nodes_introduced_at(double time) const;

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
  // Hand them back, in the order ode_aux wrote them
  template <typename It> It set_ode_aux(It it);

  // Block output adjoints the closed-form steps form and the per-cohort sweep
  // consumes, in the order ode_state visits the nodes.
  struct block_seeds {
    std::vector<double> rate;
    std::vector<double> transport;
    std::vector<double> uptake;
  };

  // Adjoints of the light field's two data vectors, one entry per knot.
  struct light_knot_adjoints {
    std::vector<double> value;
    std::vector<double> slope;
  };

  // The mirror of ode_rates: the adjoints of dydt in, the adjoints of the state
  // out through the iterator ode_state would write to.
  template <class ItIn, class ItOut>
  ItOut ode_rates_adjoint(ItIn lambda_dydt, ItOut lambda_y);

  size_t node_count() const;

  // The trait adjoints, species-major in each strategy's ad_parameters() order.
  // One trait is one input every cohort at every step reads, so this accumulates.
  std::vector<double> trait_adjoint;
  size_t trait_adjoint_size() const;
  void clear_trait_adjoint();

  // How large the recording grew on the last block, and how many blocks the one
  // tape has carried. Either climbing with the cohort count is the tape leaking.
  size_t block_recording_size = 0;
  size_t block_sweeps = 0;

  // One recording and one sweep per cohort, on a tape held across the loop,
  // seeded from the block output adjoints the closed-form steps left in `seeds`.
  void cohort_block_adjoint(const block_seeds& seeds,
                            std::vector<double>& lambda_state,
                            light_knot_adjoints& lambda_knot);

  // Carry lambda back across a node introduction. Called with this patch already
  // narrowed to the width `state_before` has, so it records the inflow boundary
  // condition n_b = birth_rate * pr_estab / g in the field `state_before`
  // builds, and reads the newcomers' rows of `lambda_after` back into
  // `lambda_before` and into trait_adjoint.
  //
  // The newcomers' rows must be contracted here and never dropped: pr_estab and
  // g are a full rate evaluation at the seed size, so they carry every trait the
  // seedling's physiology reads and every other cohort's state through the
  // field. Dropping them narrows the width just as well and returns a gradient
  // that is finite, correctly signed and wrong.
  void introduction_adjoint(const std::vector<size_t>& species_index,
                            const std::vector<double>& state_before,
                            double time_before,
                            const std::vector<double>& lambda_after,
                            std::vector<double>& lambda_before);

  // The soil drainage cascade and the water aggregation above it, which
  // between them write the soil state adjoints and the per-cohort uptake.
  void soil_adjoint(const std::vector<double>& lambda_dydt,
                    std::vector<double>& lambda_state,
                    block_seeds& seeds) const;

  // offspring_produced_survival_weighted_dt, which reaches the fecundity rate
  // and, through exp(-mortality), the mortality state directly.
  void offspring_adjoint(const std::vector<double>& lambda_dydt,
                         std::vector<double>& lambda_state,
                         block_seeds& seeds) const;

  // The light knots pulled back to each cohort's leaf area, height and log
  // density. `out` holds one entry per node in ode_state's order.
  void light_knot_adjoint(const light_knot_adjoints& lambda_knot,
                          std::vector<node_size_adjoints>& out) const;

  // The leaf-area adjoints folded onto height, and the whole of `sizes`
  // scattered into the state adjoints.
  void allometry_adjoint(const std::vector<node_size_adjoints>& sizes,
                         std::vector<double>& lambda_state) const;

  // Returns state in structure format as opposed to single 
  // vector as given by ode_state
  Rcpp::List r_get_state() const;

  // Set state of patch, based on estimate of future state estimated by the solver
  // There are two implementations.
  //   - first function is for resident runs.
  //   - second is for mutant runs.
  // The second does not calculate environment when states are updated, as mutants only experience the environment
  // The decision which to use is determined by `use_cached_environment` below
  template <typename It> It set_ode_state(It it, double time);
  template <typename It> It set_ode_state(It it, int index);

  // The state and the field: what set_ode_state establishes before it computes
  // rates. ode_rates_adjoint is taken here, where the block recordings carry the
  // rate chain and a rate evaluation in double would repeat all of it.
  template <typename It> It set_ode_state_and_field(It it, double time);

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

  void cache_ode_step();
  void cache_RK45_step(int step);
  void load_ode_step();

  // Per-accepted-step recording. record_ode_step() is called by the stepper once a
  // step is committed, so the states it keeps are the run's own. The other three are
  // what the solver's concept asks for and no more; has_recorded_field() is false, so
  // a derivative evaluation takes the same path it takes for a patch that keeps none
  // of this.
  void record_ode_step();
  void record_stage(int) {}
  void replay_step() {}
  bool has_recorded_field() const { return false; }

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
  mutable std::shared_ptr<void> block_workspace;

  int idx = 0; // used to access environment cache for mutant runs
  void compute_environment(bool rescale);
  // One field build, with every species' inflow boundary interval included or not.
  void compute_environment_once(bool rescale, bool include_boundary);
  // Evaluate every species' inflow boundary condition in the field as it stands.
  void compute_boundary_nodes();
  // The competition profile with every species' boundary interval left off: a
  // function of the ODE state alone.
  value_type compute_competition_excl_boundary(double height) const;
  // That same profile and its vertical derivative, from one pass over the species.
  std::pair<value_type, value_type>
  compute_competition_and_slope_excl_boundary(double height) const;
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
Patch<T,E>::Patch(parameters_type p, environment_type e, Control c,
                  const std::vector<strategy_type_ptr>& prepared)
  : parameters(p),
    area(p.patch_area),
    environment(e),
    control(c),
    environment_cache(6) {

  parameters.validate();

  save_RK45_cache = control.save_RK45_cache;
  survival_weighting = p.disturbance;

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
  out.compute_environment(false);
  out.environment_ptr = &out.environment;
  out.compute_rates();
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
    environment_ptr = &environment;
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
  const Internals<double>& env_vars = environment.vars;
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
template <typename T, typename E>
void Patch<T,E>::compute_boundary_nodes() {
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
  compute_environment_once(rescale, false);
  compute_boundary_nodes();
  compute_environment_once(rescale, true);
}

template <typename T, typename E>
void Patch<T,E>::compute_environment_once(bool rescale, bool include_boundary) {

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

  resource_depletion.reserve(environment_ptr->n_resources());
  for(size_t i = 0; i < environment_ptr->n_resources(); i++) {
    value_type resource_consumed = std::accumulate(species.begin(), species.end(), value_type(0.0), [i](const value_type& r, const species_type& s) -> value_type {
      return r + s.consumption_rate(i); // accumulates r from zero
    });

    // The environment's own store is Internals<double>, so the uptake is read at
    // its value here; templating that store is what would carry it further.
    resource_depletion.push_back(
      odelia::util::to_passive(resource_consumed / area));
  }
  

  environment_ptr->compute_rates(resource_depletion);

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
void Patch<T,E>::introduce_new_nodes(const std::vector<size_t>& species_index) {
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
void Patch<T,E>::remove_new_nodes(const std::vector<size_t>& species_index) {
  for (size_t i : species_index) {
    species[i].remove_newest_node();
  }
  compute_environment(false);
  compute_rates();
}

template <typename T, typename E>
std::vector<size_t> Patch<T,E>::nodes_introduced_at(double time_) const {
  std::vector<size_t> ret;
  for (size_t i = 0; i < species.size(); ++i) {
    const size_t n = species[i].size();
    if (n > 0 &&
        util::identical(species[i].node_at(n - 1).introduction_time(), time_)) {
      ret.push_back(i);
    }
  }
  return ret;
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
  environment_ptr = &environment;
  return it;
}

template <typename T, typename E>
template <typename It>
It Patch<T,E>::set_ode_state(It it, double time) {
  it = set_ode_state_and_field(it, time);

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
It Patch<T,E>::ode_rates(It it) const {
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

// The drainage cascade is lower bidiagonal and the environment transposes it;
// the aggregation above it is one trapezium per species per resource.
template <typename T, typename E>
void Patch<T,E>::soil_adjoint(const std::vector<double>& lambda_dydt,
                              std::vector<double>& lambda_state,
                              block_seeds& seeds) const {
  const size_t n_env = environment.ode_size();
  const size_t env_offset = ode_size() - n_env;
  const size_t n_resource = environment.n_resources();
  util::check_length(lambda_dydt.size(), ode_size());
  util::check_length(lambda_state.size(), ode_size());
  util::check_length(seeds.uptake.size(), node_count() * n_resource);

  const std::vector<double> lambda_env_rate(
    lambda_dydt.begin() + env_offset, lambda_dydt.end());
  std::vector<double> lambda_env_state(n_env, 0.0);
  std::vector<double> lambda_uptake(n_resource, 0.0);
  environment.compute_rates_adjoint(lambda_env_rate, lambda_env_state,
                                    lambda_uptake);
  for (size_t i = 0; i < n_env; ++i) {
    lambda_state[env_offset + i] += lambda_env_state[i];
  }

  const size_t node_stride = node_type::ode_size();
  std::vector<node_uptake_adjoints> per_node(node_count(),
                                             node_uptake_adjoints{0, 0, 0});
  for (size_t r = 0; r < n_resource; ++r) {
    if (lambda_uptake[r] == 0.0) {
      continue;
    }
    for (size_t k = 0; k < node_count(); ++k) {
      per_node[k] = node_uptake_adjoints{0, 0, 0};
    }
    size_t at = 0;
    for (size_t i = 0; i < species.size(); ++i) {
      species[i].consumption_rate_adjoint(static_cast<int>(r),
                                          lambda_uptake[r] / area,
                                          per_node.data() + at);
      at += species[i].size();
    }
    for (size_t k = 0; k < node_count(); ++k) {
      seeds.uptake[k * n_resource + r] += per_node[k].uptake;
      lambda_state[k * node_stride + HEIGHT_INDEX] += per_node[k].height;
      lambda_state[k * node_stride + T::state_size() + 1] +=
        per_node[k].log_density;
    }
  }
}

template <typename T, typename E>
void Patch<T,E>::offspring_adjoint(const std::vector<double>& lambda_dydt,
                                   std::vector<double>& lambda_state,
                                   block_seeds& seeds) const {
  const size_t node_stride = node_type::ode_size();
  const double pr_patch_survival = survival_weighting->pr_survival(time());
  size_t k = 0;
  for (size_t i = 0; i < species.size(); ++i) {
    for (size_t j = 0; j < species[i].size(); ++j, ++k) {
      const node_type& node = species[i].node_at(j);
      const double lambda_offspring =
        lambda_dydt[k * node_stride + T::state_size()];
      if (lambda_offspring == 0.0) {
        continue;
      }
      const double weight = odelia::util::to_passive(
        node.offspring_dt_dfecundity_rate(pr_patch_survival));
      seeds.rate[k * T::state_size() + FECUNDITY_INDEX] +=
        lambda_offspring * weight;
      // exp(-mortality) puts the offspring rate on a state as well as a rate,
      // and the squashed non-finite case carries a zero survival with it.
      const double fecundity_rate = odelia::util::to_passive(
        node.individual.rate(FECUNDITY_INDEX));
      lambda_state[k * node_stride + MORTALITY_INDEX] -=
        lambda_offspring * fecundity_rate * weight;
    }
  }
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
void Patch<T,E>::clear_trait_adjoint() {
  trait_adjoint.assign(trait_adjoint_size(), 0.0);
}

// The block: Individual::compute_rates at the active scalar, recorded once per
// cohort and swept once, with the leaf held constant at its declared boundary.
template <typename T, typename E>
void Patch<T,E>::cohort_block_adjoint(const block_seeds& seeds,
                                      std::vector<double>& lambda_state,
                                      light_knot_adjoints& lambda_knot) {
  using scalar = odelia::ode::active_scalar<double>;
  using active_strategy = typename at_scalar<scalar>::template apply<T>;
  using active_environment = typename at_scalar<scalar>::template apply<E>;
  using active_individual = Individual<active_strategy, active_environment>;

  // The tape and the buffers persist. No active value does: clearAll() returns
  // the tape's slot counter to zero, so a value outliving a recording aliases.

  // The templates are the exception, and never carry a slot: they are built
  // from doubles with no tape active, and are only ever copied from.
  struct block_state {
    typename scalar::tape_type tape{false};
    std::vector<double> in, out_adjoint, in_adjoint;
    std::vector<active_strategy> strategy_template;
    active_environment environment_template;
  };

  const size_t n_state = T::state_size();
  const size_t n_resource = environment.n_resources();
  const size_t n_knot = environment.light_availability.knot_count();
  const size_t n_layer = static_cast<size_t>(environment.get_soil_number_of_depths());
  const size_t node_stride = node_type::ode_size();
  const size_t env_offset = ode_size() - environment.ode_size();
  util::check_length(lambda_state.size(), ode_size());
  util::check_length(seeds.rate.size(), node_count() * n_state);
  util::check_length(seeds.transport.size(), node_count());
  util::check_length(seeds.uptake.size(), node_count() * n_resource);
  util::check_length(lambda_knot.value.size(), n_knot);
  util::check_length(lambda_knot.slope.size(), n_knot);
  if (trait_adjoint.size() != trait_adjoint_size()) {
    clear_trait_adjoint();
  }

  if (!block_workspace) {
    std::shared_ptr<block_state> fresh = std::make_shared<block_state>();
    for (size_t i = 0; i < species.size(); ++i) {
      fresh->strategy_template.push_back(
        species[i].strategy_ptr()->template rebind_from<scalar>());
    }
    block_workspace = fresh;
  }
  block_state& ws = *std::static_pointer_cast<block_state>(block_workspace);

  // This stage's soil state, time and drivers. The knot positions come across
  // because set_cohort_reads writes the field's data and not its grid.
  ws.environment_template = environment.template rebind_from<scalar>();
  ws.environment_template.light_availability.spline.set_nodes(
    environment.light_availability.spline.knots());

  size_t k = 0;
  size_t trait_at = 0;
  for (size_t i = 0; i < species.size(); ++i) {
    const size_t n_trait = species[i].strategy_ptr()->ad_parameters().size();
    for (size_t j = 0; j < species[i].size(); ++j, ++k) {
      typename active_strategy::ptr strategy =
        std::make_shared<active_strategy>(ws.strategy_template[i]);
      active_individual individual(strategy);
      active_environment block_environment = ws.environment_template;

      ws.in.resize(individual.block_input_size(block_environment));
      species[i].node_at(j).individual.block_inputs(ws.in.begin(), environment);

      ws.out_adjoint.assign(n_state + 1 + n_resource, 0.0);
      for (size_t s = 0; s < n_state; ++s) {
        ws.out_adjoint[s] = seeds.rate[k * n_state + s];
      }
      ws.out_adjoint[n_state] = seeds.transport[k];
      for (size_t r = 0; r < n_resource; ++r) {
        ws.out_adjoint[n_state + 1 + r] = seeds.uptake[k * n_resource + r];
      }

      // The transport term inside block_outputs evaluates the cohort a second
      // time at a displaced height, so the recording holds both evaluations and
      // the quotient over them.
      auto block = [&](const std::vector<scalar>& x,
                       std::vector<scalar>& y) -> void {
        individual.set_block_inputs(x.begin(), block_environment);
        individual.compute_rates(block_environment);
        individual.block_outputs(y.begin(), block_environment);
      };
      block_recording_size = odelia::ode::vector_jacobian_product(
        ws.tape, ws.in, ws.out_adjoint, block, ws.in_adjoint);
      ++block_sweeps;

      size_t at = 0;
      for (size_t s = 0; s < n_state; ++s) {
        lambda_state[k * node_stride + s] += ws.in_adjoint[at++];
      }
      for (size_t c = 0; c < n_knot; ++c) {
        lambda_knot.value[c] += ws.in_adjoint[at++];
      }
      for (size_t c = 0; c < n_knot; ++c) {
        lambda_knot.slope[c] += ws.in_adjoint[at++];
      }
      for (size_t layer = 0; layer < n_layer; ++layer) {
        lambda_state[env_offset + layer] +=
          ws.in_adjoint[at++] *
          environment.dpsi_from_soil_moist_dtheta(environment.vars.state(layer),
                                                  layer);
      }
      for (size_t p = 0; p < n_trait; ++p) {
        trait_adjoint[trait_at + p] += ws.in_adjoint[at++];
      }
      util::check_length(at, ws.in.size());
    }
    trait_at += n_trait;
  }
}

template <typename T, typename E>
void Patch<T,E>::introduction_adjoint(const std::vector<size_t>& species_index,
                                      const std::vector<double>& state_before,
                                      double time_before,
                                      const std::vector<double>& lambda_after,
                                      std::vector<double>& lambda_before) {
  using scalar = odelia::ode::active_scalar<double>;
  using active_strategy = typename at_scalar<scalar>::template apply<T>;
  const size_t node_stride = node_type::ode_size();
  const size_t n_trait = trait_adjoint_size();
  util::check_length(state_before.size(), ode_size());
  if (trait_adjoint.size() != n_trait) {
    clear_trait_adjoint();
  }

  // Each newcomer's rows sit at the end of its own species' node block, so every
  // later species and the environment sit one node higher after the widening.
  std::vector<size_t> add(species.size(), 0);
  for (size_t i : species_index) {
    ++add[i];
  }
  std::vector<size_t> carried;   // widened row each narrow row is the same as
  std::vector<size_t> newcomer;  // widened rows the introduction wrote
  size_t post = 0;
  for (size_t i = 0; i < species.size(); ++i) {
    for (size_t r = 0; r < species[i].size() * node_stride; ++r) {
      carried.push_back(post++);
    }
    for (size_t r = 0; r < add[i] * node_stride; ++r) {
      newcomer.push_back(post++);
    }
  }
  for (size_t r = 0; r < environment.ode_size(); ++r) {
    carried.push_back(post++);
  }
  util::check_length(post, lambda_after.size());
  util::check_length(carried.size(), ode_size());

  lambda_before.assign(ode_size(), 0.0);
  for (size_t j = 0; j < carried.size(); ++j) {
    lambda_before[j] = lambda_after[carried[j]];
  }

  // The twin is built with no tape active, so it holds no derivative slot that
  // the recording below could alias.
  Patch<active_strategy, typename at_scalar<scalar>::template apply<E>> active =
    this->template rebind_from<scalar>();

  std::vector<double> in(state_before);
  in.reserve(ode_size() + n_trait);
  for (size_t i = 0; i < species.size(); ++i) {
    for (const typename T::value_type* p :
         species[i].strategy_ptr()->ad_parameters()) {
      in.push_back(odelia::util::to_passive(*p));
    }
  }
  util::check_length(in.size(), ode_size() + n_trait);

  // The traits go in before the state: area_leaf(height) reads lma, so a state
  // set first is derived at the previous value. set_ode_state_and_field is what
  // computes the boundary node, and introduce_new_node only pushes it.
  const size_t n_state = ode_size();
  auto introduce = [&](const std::vector<scalar>& x,
                       std::vector<scalar>& y) -> void {
    size_t at = n_state;
    for (size_t i = 0; i < active.species.size(); ++i) {
      for (scalar* p : active.species[i].strategy_ptr()->ad_parameters()) {
        *p = x[at++];
      }
    }
    active.set_ode_state_and_field(x.begin(), time_before);
    for (size_t i : species_index) {
      active.species[i].introduce_new_node();
    }
    std::vector<scalar> widened(active.ode_size());
    active.ode_state(widened.begin());
    for (size_t j = 0; j < y.size(); ++j) {
      y[j] = widened[newcomer[j]];
    }
  };

  std::vector<double> out_adjoint(newcomer.size());
  for (size_t j = 0; j < newcomer.size(); ++j) {
    out_adjoint[j] = lambda_after[newcomer[j]];
  }
  std::vector<double> in_adjoint;
  typename scalar::tape_type tape(false);
  odelia::ode::vector_jacobian_product(tape, in, out_adjoint, introduce,
                                       in_adjoint);

  for (size_t j = 0; j < n_state; ++j) {
    lambda_before[j] += in_adjoint[j];
  }
  for (size_t p = 0; p < n_trait; ++p) {
    trait_adjoint[p] += in_adjoint[n_state + p];
  }
}

// The knots hold L = exp(-A), so the value adjoint chains through dL/dA = -L
// and the slope adjoint reaches both data vectors before either is distributed.
template <typename T, typename E>
void Patch<T,E>::light_knot_adjoint(const light_knot_adjoints& lambda_knot,
                                    std::vector<node_size_adjoints>& out) const {
  const std::vector<double>& knots = environment.light_availability.spline.knots();
  util::check_length(lambda_knot.value.size(), knots.size());
  util::check_length(lambda_knot.slope.size(), knots.size());
  util::check_length(out.size(), node_count());
  if (is_mutant_run) {
    return;
  }
  for (size_t k = 0; k < knots.size(); ++k) {
    if (lambda_knot.value[k] == 0.0 && lambda_knot.slope[k] == 0.0) {
      continue;
    }
    const double z = knots[k];
    const std::pair<value_type, value_type> as =
      compute_competition_and_slope(z);
    const double competition = odelia::util::to_passive(as.first);
    const double competition_slope = odelia::util::to_passive(as.second);
    const double transmittance = std::exp(-competition);
    // m = -L A', so the slope data carries the value with it.
    const double lambda_competition_slope =
      -(lambda_knot.slope[k] * transmittance);
    const double lambda_transmittance =
      lambda_knot.value[k] - lambda_knot.slope[k] * competition_slope;
    const double lambda_competition = -(lambda_transmittance * transmittance);
    size_t at = 0;
    for (size_t i = 0; i < species.size(); ++i) {
      species[i].compute_competition_and_slope_adjoint(
        z, lambda_competition / area, lambda_competition_slope / area,
        out.data() + at);
      at += species[i].size();
    }
  }
}

template <typename T, typename E>
void Patch<T,E>::allometry_adjoint(const std::vector<node_size_adjoints>& sizes,
                                   std::vector<double>& lambda_state) const {
  util::check_length(sizes.size(), node_count());
  util::check_length(lambda_state.size(), ode_size());
  const size_t node_stride = node_type::ode_size();
  size_t k = 0;
  for (size_t i = 0; i < species.size(); ++i) {
    for (size_t j = 0; j < species[i].size(); ++j, ++k) {
      const double darea_leaf_dheight = odelia::util::to_passive(
        species[i].node_at(j).individual.darea_leaf_dheight());
      lambda_state[k * node_stride + HEIGHT_INDEX] +=
        sizes[k].height + sizes[k].area_leaf * darea_leaf_dheight;
      lambda_state[k * node_stride + T::state_size() + 1] +=
        sizes[k].log_density;
    }
  }
}

template <typename T, typename E>
template <class ItIn, class ItOut>
ItOut Patch<T,E>::ode_rates_adjoint(ItIn lambda_dydt, ItOut lambda_y) {
  const size_t n = ode_size();
  const size_t n_resource = environment.n_resources();
  std::vector<double> lambda_in(n);
  for (size_t i = 0; i < n; ++i) {
    lambda_in[i] = *lambda_dydt++;
  }
  std::vector<double> lambda_state(n, 0.0);
  block_seeds seeds{std::vector<double>(node_count() * T::state_size(), 0.0),
                    std::vector<double>(node_count(), 0.0),
                    std::vector<double>(node_count() * n_resource, 0.0)};
  // The strategy rate adjoints the stage recursion supplies, and beside them the
  // transport term's, which is a block output rather than a closed-form seed.
  for (size_t k = 0; k < node_count(); ++k) {
    for (size_t s = 0; s < T::state_size(); ++s) {
      seeds.rate[k * T::state_size() + s] =
        lambda_in[k * node_type::ode_size() + s];
    }
    seeds.transport[k] =
      lambda_in[k * node_type::ode_size() + T::state_size() + 1];
  }

  soil_adjoint(lambda_in, lambda_state, seeds);
  offspring_adjoint(lambda_in, lambda_state, seeds);

  // Per cohort: record the block, seed it from `seeds`, sweep, and read the
  // state, knot and trait adjoints back.
  light_knot_adjoints lambda_knot{
    std::vector<double>(environment.light_availability.spline.knots().size(), 0.0),
    std::vector<double>(environment.light_availability.spline.knots().size(), 0.0)};
  cohort_block_adjoint(seeds, lambda_state, lambda_knot);

  std::vector<node_size_adjoints> sizes(node_count(),
                                        node_size_adjoints{0, 0, 0});
  light_knot_adjoint(lambda_knot, sizes);
  allometry_adjoint(sizes, lambda_state);

  for (size_t i = 0; i < n; ++i) {
    *lambda_y++ = lambda_state[i];
  }
  return lambda_y;
}

}

#endif
