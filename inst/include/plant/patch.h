// -*-c++-*-
#ifndef PLANT_PLANT_PATCH_H_
#define PLANT_PLANT_PATCH_H_

#include <plant/parameters.h>
#include <plant/species.h>
#include <plant/util.h>
#include <plant/clamp_sites.h>
#include <plant/gradient_status.h>
#include <odelia/ode_interface.hpp>
#include <odelia/gradient.hpp>

#include <plant/disturbance_regime.h>

#include <algorithm>
#include <cmath>
#include <limits>
#include <string>

using namespace Rcpp;

namespace plant {

// A strategy whose inner solve makes a choice the state leaves open, and which keeps
// what it chose so a pass re-running the model over the same states places it rather
// than searching for it again. The patch hands down the address of the rate
// evaluation now running; a strategy with no such solve declares neither member and
// the forwarding compiles away.
template <typename T>
concept KeepsSolvedChoices =
  requires(T& s, odelia::ode::recorded_stage at, bool keeping) {
    s.begin_stage(at, keeping);
    s.end_stage();
  };

// A strategy or environment named at scalar U: the type its own rebind returns.
// Named from the factory rather than from a second alias beside it, so there is
// one answer to "what is this at another scalar" and a type that cannot rebind
// says so here rather than further in.
template <typename X, typename U>
using at_scalar = decltype(std::declval<const X&>().template rebind_from<U>());

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
            class T2 = at_scalar<T, U>, class E2 = at_scalar<E, U>>
  Patch<T2,E2> rebind_from() const;

  // The values rebind_from copies, written into a patch that already exists: the
  // strategies, the environment, the node structure and its birth stamps, by the
  // same calls in the same order, so what is left is what a rebind would have
  // returned. An active patch holding what a recording wrote carries that recording's tape
  // slots into the next one, and the sweep comes back wrong with nothing raised,
  // so it is assigned before every recording rather than once per step.
  template <typename T1, typename E1>
  void assign_from(const Patch<T1,E1>& src);

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

  // One entry per species gaining a node. This is what the solver carries back to
  // the patch when it walks a recording it widened, and it never reads it.
  using widening = std::vector<size_t>;

  // The one insertion. One node per species named, stamped with the time it is
  // introduced at -- the time is an argument because it is the schedule's, and a
  // patch that reads it off its own clock inserts whatever the last load left
  // there. Brings the field and the rates up to date, because a node changes both.
  void introduce_nodes(const widening& species_index, double time);

  // One species, for a caller building a patch by hand. The time is an argument
  // for the same reason: routed through the clock instead, every node of a
  // species shares a date and the grid the birth-date coordinate integrates over
  // is tied.
  void r_introduce_new_node(util::index species_index, double time) {
    introduce_nodes(widening{species_index.check_bounds(size())}, time);
  }

  // The insertion as a map: the state before it in, the whole widened state out,
  // and nothing rebuilt. This is the one the sweep transposes, so it runs at
  // whatever scalar it is called on and loads the state itself rather than
  // asking the caller to. It leaves this patch holding what it added.
  template <typename Widening, typename It>
  void widened_state(const odelia::ode::recorded_insertion<Widening>& insertion,
                     It x, std::vector<value_type>& y);

  // Open to better ways to test whether nodes have been introduced
  int node_ode_size() const {
    int node_ode_size = ode_size() - environment.ode_size();
    return(node_ode_size);
  }

  const species_type& at_species(size_t species_index) const {
    return species[species_index];
  }

  // Patch disturbance
  std::shared_ptr<Disturbance_Regime> survival_weighting;

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

  size_t node_count() const;

  // Every species' differentiable parameters, species-major: the order a trait
  // row is indexed in, answered once rather than walked by each caller.
  std::vector<typename T::value_type*> ad_parameters();

  size_t trait_adjoint_size() const;
  // The same order, named. Each name carries its species index, because
  // concatenating the strategies' own names repeats every one of them per
  // species and character indexing then resolves each to species one's column,
  // which an unknown-name check cannot see.
  std::vector<std::string> trait_adjoint_names() const;
  // What an exactly-zero entry in each column would mean, same order and width.
  std::vector<gradient_status::Kind> trait_adjoint_zero_classes() const;

  // The environment's own clamp tally. The site list is shared with the
  // strategy's, so a caller adds the two rather than reading them apart; this is
  // the one route to it, because the environment itself is not the caller's.
  clamp_counter& environment_clamps() const { return environment.clamps; }

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

  // Set state of patch, based on estimate of future state estimated by the solver,
  // computing the environment as it goes.
  template <typename It> It set_ode_state(It it, double time);

  // The state, and the address of the choices the run made at the rate evaluation
  // this one stands in for. The third of the family: nothing completes the state on
  // the plain loader, the inflow condition's second evaluation on the recorded one,
  // and here what an inner solve chose that the state leaves open.
  //
  // The address is opened HERE and not at the rates, because the inflow condition's
  // own leaf solves happen inside the field build below; it is closed where the rates
  // are read, so one rate evaluation is one address and a solve outside a step
  // neither keeps nor places.
  template <typename It>
  It set_ode_state(It it, double time, odelia::ode::recorded_stage at);

  // A recorded state loaded as the run itself carries it. set_ode_state evaluates
  // the inflow condition in the field that leaves the boundary interval off, then
  // rebuilds the field including it; the run then rates the nodes and evaluates
  // the condition a second time, in that second field. It is the second value an
  // introduced node inherits and the census reads, so reloading a state without it
  // linearises a boundary node the trajectory never carried.
  template <typename It> It set_recorded_state(It it, double time);

  // The same, at the node structure the record describes: whatever the patch was
  // seeded with, plus one node per species named by each of the first `applied`
  // insertions, each stamped with the recorded time it was inserted at.
  //
  // A reconciliation rather than a sequence of insertions and removals, so it is
  // idempotent: being at a recorded step twice is being there once, and a walk
  // can be run again over a recording it has already walked. Every node the
  // structure gains carries three numbers that are not ODE state -- its birth
  // date, and the patch density and survival there -- and all three are functions
  // of the time it was inserted at, which is why the record needs to hold only
  // WHICH species gained a node after which step.
  template <typename Widening>
  void set_recorded_state(
      const std::vector<value_type>& y, double time,
      const std::vector<odelia::ode::recorded_insertion<Widening>>& insertions,
      size_t applied);

  // The inflow condition alone, in the field as it now stands. Public because the
  // two evaluations above have to be taken one at a time to be told apart.
  void compute_boundary_nodes();

  // Forget what every species' inner solve chose. The run calls it where it starts,
  // because the record outlives this patch: it is shared with every copy and every
  // rebound version of it, which is what lets a sweep read what the run wrote.
  void clear_solved_choices() {
    if constexpr (KeepsSolvedChoices<strategy_type>) {
      for (species_type& s : species) {
        s.strategy_ptr()->leaf_points->clear();
      }
    }
  }

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
  species_type r_at(util::index species_index) const {
    return species[species_index.check_bounds(size())];
  }
  // This is only here because it wraps a private function.
  void r_compute_environment() {compute_environment(false);}

  // Whether this run is the one whose choices a later pass reads back. The states
  // are the SOLVER's record now, beside the times and the sizes; what is left here
  // is the one bit the loader needs, and the run sets it.
  bool recording = false;

  void add_strategies(std::vector<strategy_type> strategies);
  void overwrite_strategies(std::vector<strategy_type> strategies);

private:
  // A patch whose species take the prepared strategies given, rather than
  // preparing the ones in the parameters. rebind_from's only route in.
  Patch(parameters_type p, environment_type e, plant::Control c,
        const std::vector<strategy_type_ptr>& prepared);

  // What each species' reduction had accumulated at each knot before its
  // closing trapezium, kept from the field built without the boundary interval
  // so the field built with it costs one trapezium per knot rather than a second
  // walk over every node.
  std::vector<std::vector<typename species_type::competition_split>>
    competition_capture;
  size_t capture_at = 0;
  void compute_environment_excl_capturing(bool rescale);
  void compute_environment_closing(bool rescale);

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
  // Guard the birth-date coordinate's quadrature grid: it is the per-node
  // introduction times, so nodes sharing one give zero-width intervals that drop
  // silently out of the integral. Reached by a schedule carrying a repeated time
  // and by a patch whose nodes were seeded or imported without per-node times.
  void check_birth_dates_distinct() const;

  // One node per species named, stamped from the time alone. The insertion and
  // the reconciling loader share it, so a node the run made and a node a walk
  // rebuilt are stamped by the same expression.
  void push_nodes(const widening& species_index, double time);

  parameters_type parameters;

  double area;
  environment_type environment;
  std::vector<species_type> species;

  //TODO(#476): Move into environment?
  std::vector<value_type> resource_depletion;

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
    control(c) {

  parameters.validate();

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
    control(c) {

  parameters.validate();

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
  p2.ode_step_sizes = parameters.ode_step_sizes;
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
  // soil state, so restore the current one before the spline is rebuilt. Moved,
  // not copied: this is the last read of env, and the object carries the light
  // spline and a vector per soil layer.
  out.environment = std::move(env);
  out.environment.set_shading_model(control.shading_model,
                                    control.ppa_layer_optical_depth,
                                    control.ppa_layer_smoothing);
  // The field and the boundary node are left to the caller. Every caller sets a
  // state through set_ode_state or set_recorded_state before reading
  // either, and both rebuild the field and re-evaluate the inflow condition, so
  // computing them here solves the boundary leaf twice per right-hand side and
  // then discards it. A caller that reads before setting gets an unbuilt field.
  return out;
}

template <typename T, typename E>
template <typename T1, typename E1>
void Patch<T,E>::assign_from(const Patch<T1,E1>& src) {
  using U = value_type;
  if (species.size() != src.species.size()) {
    util::stop("assign_from: this patch runs a different number of species from "
               "the one it is assigned from");
  }

  environment = src.environment.template rebind_from<U>();
  environment.set_shading_model(control.shading_model,
                                control.ppa_layer_optical_depth,
                                control.ppa_layer_smoothing);

  for (size_t i = 0; i < species.size(); ++i) {
    *species[i].strategy_ptr() =
      src.species[i].strategy_ptr()->template rebind_from<U>();
    // Back to the state reset() leaves a species in, so the nodes pushed below
    // are the same copies of the same boundary node a rebind pushes.
    species[i].clear();
    species[i].resize_consumption_rates(environment.n_resources());
    for (size_t j = 0; j < src.species[i].size(); ++j) {
      species[i].introduce_new_node();
    }
    species[i].set_birth_state(src.species[i].node_times(),
                               src.species[i].r_patch_densities(),
                               src.species[i].r_pr_patch_survival_at_birth());
  }

  std::vector<U> node_state(src.node_ode_size());
  odelia::ode::ode_state(src.species.begin(), src.species.end(),
                         node_state.begin());
  odelia::ode::set_ode_state(species.begin(), species.end(),
                             node_state.begin());
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
      const std::vector<double> times = species[i].node_times();
      double repeated = times.empty() ? 0.0 : times.front();
      for (size_t k = 1; k < times.size(); ++k) {
        if (times[k] == times[k - 1]) {
          repeated = times[k];
          break;
        }
      }
      util::stop("Species " + util::to_string(i + 1) + " has nodes sharing an "
                 "introduction time (" + util::to_string(repeated) + "), which "
                 "the birth-date size-density coordinate integrates over: the "
                 "repeated nodes span zero width and drop out of the competition "
                 "and resource integrals. Remove the repeat from the node "
                 "schedule, supply per-node introduction times "
                 "(parameters$initial_node_times) with an initial state, or run "
                 "with control$node_density_in_birth_date = FALSE.");
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
  // Double, not the working scalar: this reads the state to test it for
  // finiteness and never differentiates it, and the iterator write converts on
  // the way out. Held at the working scalar it was a slot per entry per stage
  // for a value flattened on the next line.
  std::vector<double> env_state(environment.ode_size());
  environment.ode_state(env_state.begin());
  for (size_t i = 0; i < env_state.size(); ++i) {
    const double state_i = env_state[i];
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
    const value_type h = species[i].height_max();
    if (h > ret) {
      ret = h;
    }
  }
  return ret;
}

template <typename T, typename E>
typename Patch<T,E>::value_type
Patch<T,E>::compute_competition(double height) const {
  value_type tot = 0.0;
  for (size_t i = 0; i < species.size(); ++i) {
    tot += species[i].compute_competition(height) / area;
  }
  return tot;
}

template <typename T, typename E>
std::pair<typename Patch<T,E>::value_type, typename Patch<T,E>::value_type>
Patch<T,E>::compute_competition_and_slope(double z) const {
  value_type tot = 0.0, tot_slope = 0.0;
  for (size_t i = 0; i < species.size(); ++i) {
    const std::pair<value_type, value_type> fs =
      species[i].compute_competition_and_slope(z);
    tot       += fs.first / area;
    tot_slope += fs.second / area;
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

// Evaluate every species' inflow boundary condition in the field as it currently
// stands. Owned by the field build rather than by compute_rates(), so that the
// field reads a boundary density derived from this state instead of one carried
// from the previous evaluation.
// A mutant experiences a recorded environment rather than shaping one, and this
// reads the patch's own, so it is not the mutant's condition to evaluate.
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
  if (size() == 0) {
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
  // The boundary node is one end of the birth-date quadrature and its birth date
  // is the current time, so refresh it before the profile is built. Its own stamp
  // is set in compute_rates(), which the stepper calls *after* the set_ode_state()
  // that brings us here, so reading that stamp would use the previous derivs
  // call's time and shorten the boundary interval by a Runge-Kutta stage. The
  // measured effect is below 1e-6, but the interval is then a function of the step
  // size, which the spatial quadrature has no business depending on, and the whole
  // integral *is* that one segment while a species has a single node. No-op for
  // the height coordinate, where this abscissa is the constant initial height.
  for (auto& s : species) {
    s.set_new_node_birth_date(environment.time);
  }
  // Cleared without dropping capacity: assign() would, and the reallocation that
  // follows copies active scalars, which each cost a slot to register and another
  // to release.
  competition_capture.resize(size());
  for (auto& c : competition_capture) {
    c.clear();
  }
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
  if (size() > 0) {
    environment.compute_environment(f, height_max(), rescale);
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
      if (capture_at >= competition_capture[i].size()) {
        util::stop("compute_environment_closing: the field without the "
                   "boundary interval was built over a different knot set");
      }
      const std::pair<value_type, value_type> fs =
        species[i].close_competition_and_slope(competition_capture[i][capture_at], x);
      tot       += fs.first / area;
      tot_slope += fs.second / area;
    }
    ++capture_at;
    return {tot, tot_slope};
  };
  if (size() > 0) {
    environment.compute_environment(f, height_max(), rescale);
  }
}



template <typename T, typename E>
void Patch<T,E>::compute_rates() {

  // Computes rates of change for the patch, including all the component species,
  // against the environment the patch experiences.
  environment_type& env = environment;
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
void Patch<T,E>::introduce_nodes(const widening& species_index, double time) {
  push_nodes(species_index, time);

  // A schedule carrying the same time twice for one species stamps two nodes
  // with it, and the grid both reductions integrate over is those times.
  check_birth_dates_distinct();

  compute_environment(false);

  // New nodes have just changed the state and the light field, so the stored
  // rates now describe neither. The solver reads them next without checking, so
  // they have to be brought up to date here.
  compute_rates();
}

// The three numbers a node carries that the ODE state does not: its birth date,
// and the patch density and survival there. All three are functions of the time
// it is introduced at, which is what lets a record hold only the time.
template <typename T, typename E>
void Patch<T,E>::push_nodes(const widening& species_index, double time) {
  const double patch_density = survival_weighting->density(time);
  const double pr_survival = survival_weighting->pr_survival(time);
  for (size_t i : species_index) {
    if (i >= species.size()) {
      util::stop("introduce_nodes: species " +
                 util::to_string(static_cast<int>(i)) + " of a patch holding " +
                 util::to_string(static_cast<int>(species.size())));
    }
    species[i].introduce_new_node(time, patch_density, pr_survival);
  }
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

  // Build the field the rates will be taken in.
  compute_environment(true);

  return it;
}

// The same load, with the choices the run made at this rate evaluation reachable by
// every species that keeps any. See the declaration for where it is closed.
template <typename T, typename E>
template <typename It>
It Patch<T,E>::set_ode_state(It it, double time, odelia::ode::recorded_stage at) {
  if constexpr (KeepsSolvedChoices<strategy_type>) {
    for (species_type& s : species) {
      s.strategy_ptr()->begin_stage(at, recording);
    }
  }
  return set_ode_state(it, time);
}

// The second evaluation of the inflow condition, in the field the first one was
// folded into. See the declaration for why a reloaded state needs it.
template <typename T, typename E>
template <typename It>
It Patch<T,E>::set_recorded_state(It it, double time) {
  it = set_ode_state(it, time);
  compute_boundary_nodes();
  return it;
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
It Patch<T,E>::ode_rates(It it) {
  compute_rates();
  it = odelia::ode::ode_rates(species.begin(), species.end(), it);
  it = environment.ode_rates(it);
  // One rate evaluation is one address: what a solve outside a step chooses belongs
  // to no recorded stage, so nothing is open to keep it in or place it from.
  if constexpr (KeepsSolvedChoices<strategy_type>) {
    for (species_type& s : species) {
      s.strategy_ptr()->end_stage();
    }
  }
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
std::vector<gradient_status::Kind>
Patch<T,E>::trait_adjoint_zero_classes() const {
  std::vector<gradient_status::Kind> ret;
  ret.reserve(trait_adjoint_size());
  for (size_t i = 0; i < species.size(); ++i) {
    for (const gradient_status::Kind k :
         species[i].strategy_ptr()->ad_parameter_zero_classes()) {
      ret.push_back(k);
    }
  }
  return ret;
}

template <typename T, typename E>
template <typename Widening, typename It>
void Patch<T,E>::widened_state(
    const odelia::ode::recorded_insertion<Widening>& insertion, It x,
    std::vector<value_type>& y) {
  set_recorded_state(x, insertion.time);
  push_nodes(insertion.what, insertion.time);
  y.assign(ode_size(), value_type(0.0));
  ode_state(y.begin());
}

template <typename T, typename E>
template <typename Widening>
void Patch<T,E>::set_recorded_state(
    const std::vector<value_type>& y, double time,
    const std::vector<odelia::ode::recorded_insertion<Widening>>& insertions,
    size_t applied) {
  // What the patch was seeded with is the structure no insertion accounts for;
  // everything above it is one node per species named by the insertions so far.
  const size_t n_species = species.size();
  const bool seeded = !parameters.initial_state.empty();
  std::vector<std::vector<double>> when(n_species);
  for (size_t j = 0; j < applied; ++j) {
    for (size_t i : insertions.at(j).what) {
      if (i >= n_species) {
        util::stop("set_recorded_state: insertion " +
                   util::to_string(static_cast<int>(j)) + " names species " +
                   util::to_string(static_cast<int>(i)) +
                   " of a patch holding " +
                   util::to_string(static_cast<int>(n_species)));
      }
      when[i].push_back(insertions[j].time);
    }
  }

  bool moved = false;
  for (size_t i = 0; i < n_species; ++i) {
    const size_t base = seeded ? parameters.n_initial_cohorts.at(i) : 0;
    const size_t target = base + when[i].size();
    while (species[i].size() > target) {
      species[i].remove_newest_node();
      moved = true;
    }
    if (species[i].size() < target) {
      // The values a pushed node carries are the state's and arrive with the
      // load below; only its stamps are its own, and those are set here.
      widening one{i};
      while (species[i].size() < target) {
        push_nodes(one, when[i][species[i].size() - base]);
        moved = true;
      }
    }
  }
  if (moved) {
    check_birth_dates_distinct();
    compute_environment(false);
    compute_rates();
  }
  util::check_length(y.size(), ode_size());
  set_recorded_state(y.begin(), time);
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
    active.widened_state(
        odelia::ode::recorded_insertion<widening>{species_index, 0, time_before},
        x.begin(), y);
    for (size_t r = 0; r < n_out; ++r) {
      ret[r][c] = xad::derivative(y[r]);
    }
    for (size_t i : species_index) {
      active.species[i].remove_newest_node();
    }
  }
  return ret;
}


}

#endif
