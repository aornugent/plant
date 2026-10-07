// -*-c++-*-
#ifndef PLANT_PLANT_PATCH_H_
#define PLANT_PLANT_PATCH_H_

#include <plant/parameters.h>
#include <plant/species.h>
#include <plant/util.h>
#include <plant/clamp_sites.h>
#include <odelia/implicit_node.hpp>
#include <odelia/ode_interface.hpp>
#include <odelia/ode_util.hpp> // odelia::util::stop_domain

#include <plant/disturbance_regime.h>
#include <plant/with_slope.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <concepts>
#include <limits>
#include <memory>
#include <sstream>
#include <string>
#include <type_traits>
#include <utility> // std::pair

using namespace Rcpp;

namespace plant {

// A strategy or environment named at scalar U: the type its own rebind returns.
// Named from the factory rather than from a second alias beside it, so there is
// one answer to "what is this at another scalar" and a type that cannot rebind
// says so here rather than further in.
template <typename X, typename U>
using at_scalar = decltype(std::declval<const X&>().template rebind_from<U>());

// The scalar a forward tangent runs on, and the two things done to one: seed a
// direction, read the derivative it produced. Brought in here once, so no header
// in this package names the AD library at all -- and so that seeding an adjoint
// scalar by mistake, which the library spells identically and answers with
// nothing, does not compile.
using tangent = odelia::ode::tangent_scalar<double>;
using odelia::ode::seed_direction;
using odelia::ode::derivative_along;

// One accepted step. The state widens at an introduction, so the record is ragged.
// ⚠️ DOUBLES, NOT AN ENVIRONMENT. A recording is taken at `double` and handed to
// whatever scalar a later pass runs at -- the sweep lifts the system and replays
// the same rows -- so a recorded number that carried the working scalar could not
// be handed over at all. It is the same reason `leaf_solved_point` is a collar
// and an arm rather than a Leaf. And it is the cheap form besides: what an
// environment costs is its interpolant's adaptive builder and band-solve
// workspace, which no replay reads.
struct recorded_field {
  // What Environment::set_knot_data() reads: the light field as its interpolant
  // held it.
  std::vector<double> knot_data;
  // The environment's own ODE state -- for TF24 the soil water, which is half of
  // what a competitor competes for.
  std::vector<double> state;
  double time = 0.0;
};

// A strategy whose rates change form where one of its auxiliaries changes sign.
template <typename T>
concept NamesSignValue = requires(const T& s) {
  { s.sign_value_aux() } -> std::same_as<int>;
};

// An environment whose soil a step can take alone: its first n_resources()
// states, which read the plants only through the uptake from each.
template <typename E>
concept SoilStepsAlone =
  requires(const E& e, std::vector<typename E::value_type>& u) {
  e.alone_inputs(u);
  { e.uptake_share() } -> std::convertible_to<double>;
};

// One rate evaluation's record: what each species solved for, the field it was
// taken in, and the field it built. Templated on the strategy's own recorded type
// rather than nested in Patch, so that a patch and the same patch lifted to an
// active scalar name ONE type -- see the alias in Patch for why that is
// load-bearing.
template <typename StrategySolved>
struct patch_solved_values {
  std::vector<StrategySolved> strategies;
  // Null where the evaluation built its own. Shared, because a walk copies every
  // row it is seeded with.
  std::shared_ptr<const recorded_field> field;
  // Null where no run kept it. An invasion of the run moves it into `field`.
  std::shared_ptr<const recorded_field> built;
};

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
  with_slope<value_type>
  compute_competition_and_slope(double z) const;
  // The same pair for R, which takes only doubles: value first, then slope.
  std::vector<double> r_compute_competition_and_slope(double z) const {
    const with_slope<value_type> fs =
      compute_competition_and_slope(z);
    return {odelia::util::to_passive(fs.value),
            odelia::util::to_passive(fs.slope)};
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
  // Accumulated across the run and cleared by reset(); see the definition.
  void collect_competition_errors(const std::vector<size_t>& added);
  // What drives schedule refinement; see the definition.
  std::vector<std::vector<double>> refinement_error_by_node() const;

  // The species gaining a node at one introduction. Named here because it is
  // this model's, not the solver's: a solver that walks a recording of this patch
  // never learns what an inserted entry is.
  using introduction = std::vector<size_t>;

  // One node per species named, stamped with the time it is introduced at, for a
  // caller building a patch by hand -- a run applies its schedule's entries with
  // apply_insertion(). The time is an argument because it is the schedule's.
  // Brings the field and the rates up to date, because a node changes both.
  void introduce_nodes(const introduction& species_index, double time);

  // One species, for a caller building a patch by hand. The time is an argument
  // for the same reason: routed through the clock instead, every node of a
  // species shares a date and the grid the birth-date coordinate integrates over
  // is tied.
  void r_introduce_new_node(util::index species_index, double time) {
    introduce_nodes(introduction{species_index.check_bounds(size())}, time);
  }

  // Apply one non-introduction scheduled event (issue #628). Called between
  // solver legs, so it may change state and ODE size but must not move the
  // clock -- the caller re-reads both with set_state_from_system(). An action
  // is free to reach the answer however it likes, including by integrating its
  // own fast sub-model over the event's nominal duration with the patch's
  // demography frozen; to this solver that is still one instantaneous jump.
  EventRecord apply_event(const NodeScheduleEvent& event);

  // Scale every selected node's density by phi, and report how many nodes were
  // touched. Shared by harvest and climate extremes, which differ only in how phi
  // is chosen per node -- so there is one place where a discrete removal meets
  // the state, and one place to get it right.
  //
  // Both accountings have to move together. On the birth-date coordinate
  // log_density is log(birth_rate) - mortality(t), derived rather than a state,
  // and moving both keeps it current. Moving
  // only log_density leaves fecundity undiscounted (it is weighted by
  // exp(-mortality) independently); moving only mortality leaves the standing
  // density untouched. Either way the two views of the same cohort disagree.
  // Returns {nodes_affected, density_removed} -- the count is bookkeeping, the
  // density is the quantity actually taken out.
  template <typename Select>
  std::pair<size_t, double> scale_node_densities(size_t species_index,
                                                 Select select);

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


  // Every species' differentiable parameters, species-major: the order a trait
  // row is indexed in, answered once rather than walked by each caller.
  std::vector<typename T::value_type*> ad_parameters();

  // Every active value this patch holds, wherever it lives. A walk that keeps a
  // active patch across recordings hands each of these back before it clears the
  // tape, and the tape then says whether any was missed.
  //
  // `area` and the disturbance regime are not here: they are double. Nor are the
  // per-node error tallies, for the same reason.
  template <class F>
  void for_each_active(F&& f) {
    odelia::ode::visit_active(f, parameters, environment, species,
                              competition_capture);
  }

  size_t trait_adjoint_size() const;
  // The same order, named. Each name carries its species index, because
  // concatenating the strategies' own names repeats every one of them per
  // species and character indexing then resolves each to species one's column,
  // which an unknown-name check cannot see.
  std::vector<std::string> trait_adjoint_names() const;

  // The environment's own clamp tally. The site list is shared with the
  // strategy's, so a caller adds the two rather than reading them apart; this is
  // the one route to it, because the environment itself is not the caller's.
  clamp_counter& environment_clamps() const { return environment.clamps; }

  // The insertion map's whole Jacobian, by forward tangent: one row per wider state entry
  // and one column per input, the state's entries first and the traits after.
  // Forming it entirely is what localises a disagreement to a cell, which a
  // contraction against the transpose cannot do.
  std::vector<std::vector<double>>
  introduction_jacobian(const std::vector<size_t>& species_index,
                        const std::vector<double>& state_before,
                        double time_before);

  // The strategy's own bound on its states, read by the solver: a step landing
  // on a state it refuses is rejected and retried smaller, rather than
  // committed and clamped by whoever reads it next. Walks the state vector in
  // the order ode_state writes it.
  bool ode_state_valid(const std::vector<double>& y) const;

  // What the adaptive stepper multiplies each state's error level by on a step
  // starting at `time` (see Control's ode_weight_soil, ode_weight_times and
  // ode_weight_max).
  void error_weights(double time, std::vector<double>& w) const;
  std::vector<double> r_error_weights(double time) const {
    std::vector<double> w;
    error_weights(time, w);
    return w;
  }

  // The environment's soil, which a step takes alone where the plants' uptake is
  // under ode_soil_alone_share of the water moving through it. Its layers open
  // the environment's states, which close the patch's.
  std::pair<size_t, size_t> alone_block() const requires SoilStepsAlone<E> {
    return {ode_size() - environment.ode_size(), environment.n_resources()};
  }
  bool steps_alone() const requires SoilStepsAlone<E> {
    return environment.uptake_share() < control.ode_soil_alone_share;
  }
  void alone_inputs(std::vector<value_type>& u) const
    requires SoilStepsAlone<E> {
    environment.alone_inputs(u);
  }
  template <class U>
  void alone_rates(double t, const std::vector<U>& y, const std::vector<U>& u,
                   std::vector<U>& out) const requires SoilStepsAlone<E> {
    environment.alone_rates(t, y, u, out);
  }

  // Returns state in structure format as opposed to single 
  // vector as given by ode_state
  Rcpp::List r_get_state() const;

  // Set state of patch, based on estimate of future state estimated by the solver,
  // computing the environment as it goes.
  template <typename It> It set_ode_state(It it, double time);

  // What one of this patch's rate evaluations solves for: each species' own inner
  // solves, the field the evaluation was taken in, and the field it built where a
  // run keeps it.
  //
  // ⚠️ A RUN KEEPS ITS FIELD IN `built`, WHICH NO EVALUATION READS. A run's own
  // field is recomputable from its state, so a sweep taken in a recorded copy would
  // drop the field's derivative with every number finite. An invasion moves the
  // run's fields into `field`: an invader's field is not its own, and a recorded
  // copy gives the derivative it has, zero.
  // ⚠️ AN ALIAS AND NOT A NESTED STRUCT: a nested struct of a class template is a
  // distinct type per instantiation, and the sweep hands a recording taken on the
  // double patch to the patch at an active scalar.
  using solved_values =
      patch_solved_values<odelia::ode::solved_values_t<strategy_type>>;

  // The walk opens this extent before the state is loaded, when the field is not
  // built yet, so compute_environment() reads the open slot's field there.
  void store_solved(solved_values& into) {
    into.strategies.resize(species.size());
    for (size_t i = 0; i < species.size(); ++i) {
      odelia::ode::store_solved(*species[i].strategy_ptr(), into.strategies[i]);
    }
    storing = &into;
  }
  // The species count is checked: a recorded list of another length would
  // otherwise read as a species that solved for nothing.
  void load_solved(const solved_values& from) {
    util::check_length(from.strategies.size(), species.size());
    for (size_t i = 0; i < species.size(); ++i) {
      odelia::ode::load_solved(*species[i].strategy_ptr(), from.strategies[i]);
    }
    loading = &from;
  }
  void end_solved() {
    for (species_type& s : species) {
      odelia::ode::end_solved(*s.strategy_ptr());
    }
    storing = nullptr;
    loading = nullptr;
  }

  // Keep the field each evaluation builds; set by a run that keeps its states (see
  // `solved_values`).
  void set_keep_field(bool keep) { keep_field = keep; }
  bool keeps_field() const { return keep_field; }

  // What a step records of each node it split, and of the field it sampled.
  using split_blocks = std::vector<odelia::ode::split_block<solved_values>>;
  using split_samples = std::vector<solved_values>;
  // Each node's net production after an evaluation, in the order ode_state writes
  // the nodes.
  void sign_values(std::vector<double>& out) const requires NamesSignValue<T>;
  // Integrates each node whose net production changed sign in the step in pieces
  // between its sign changes; true if it evaluated anything, split or not.
  template <class Step>
  bool split_sign_changes(const Step& step, split_samples& samples,
                          split_blocks& record) requires NamesSignValue<T>;
  // The sweep's split: each recorded sign change moves as the zero of its node's net
  // production, and the node is integrated again loading what the run solved for.
  template <class Step>
  void split_as_recorded(const Step& step, const split_samples& samples,
                         const split_blocks& recorded) requires NamesSignValue<T>;
  // A walk integrates each invader's node born when a node the run split was born
  // in pieces at the run's sign changes, held, in the field the run sampled.
  template <class Step, class Row>
  void take_recorded_splits(const Step& step, const Row& recorded,
                            split_samples& samples, split_blocks& record)
    requires NamesSignValue<T>;
  // The birth dates of the walked run's nodes, in the order its state holds them at
  // its end, by which a walk finds an invader's copy of the run's node.
  void set_run_birth_dates(std::vector<double> dates) {
    run_birth_dates = std::move(dates);
  }
  const std::vector<double>& get_run_birth_dates() const { return run_birth_dates; }

  // The entries this patch is run on, which a walk applies and works out the
  // shape at a step from. Set where a run begins, since the schedule can change
  // between runs.
  void set_schedule(std::shared_ptr<const std::vector<schedule_entry>> entries) {
    schedule = std::move(entries);
  }
  // How many nodes species `i` holds at a recorded `time`: what it was seeded
  // with, plus the introductions strictly below. Strictly, because an
  // introduction at `time` follows the step recorded there.
  size_t nodes_at(size_t i, double time) const;
  // Become that shape, stamping what it gains with the date the schedule gives.
  // Idempotent, so a walk can be run again over a recording it has walked.
  void reshape_to(double time);

  // The entry at `time`, applied to the state held: its events in schedule order,
  // then one node per species it introduces. What each event did is returned.
  std::vector<EventRecord> apply_insertion(double time);
  // One node per species named, stamped with `time`, each a copy of the boundary
  // node as the last rate evaluation left it.
  void push_nodes(const introduction& species_index, double time);
  // The same introductions as a map from the state below them, `x`, at any
  // scalar: that state's evaluation, then the nodes. `y` is the wider state.
  template <typename It>
  void introduce_from(const introduction& species_index, double time, It x,
                      std::vector<value_type>& y);

  // The inflow condition alone, in the field as it now stands. Public because the
  // two evaluations of it have to be taken one at a time to be told apart.
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
  species_type r_at(util::index species_index) const {
    return species[species_index.check_bounds(size())];
  }
  // This is only here because it wraps a private function.
  void r_compute_environment() {compute_environment();}

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
  void compute_environment_excl_capturing();
  void compute_environment_closing();

  void compute_environment();
  void compute_rates();

  // The field at each of the step's five sample fractions, kept in each sample's
  // slot as an evaluation keeps its own: the light's knot data, then the soil's.
  using field_samples = std::array<std::vector<value_type>, 5>;
  template <class Step, class Samples>
  void sample_field(const Step& step, field_samples& field, Samples& samples);
  // How many components node `node` holds, counting nodes as ode_state() does.
  size_t node_width(size_t node) const;
  // A node's net production at fraction u of the step, its own components read
  // from the dense output and its field from the samples, evaluated at `time`.
  template <class Step, class U>
  value_type node_value_at(const Step& step, const field_samples& field,
                           size_t node, size_t first, size_t width, const U& u,
                           double time)
    requires NamesSignValue<T>;
  // A node integrated again over the step in pieces meeting at `split_at`, each
  // evaluation inside the extent `solved()` opens; its end goes into the step's.
  template <class Step, class U, class Solved>
  void split_node(const Step& step, const field_samples& field, size_t node,
                  size_t first, size_t width, const std::vector<U>& split_at,
                  Solved&& solved)
    requires NamesSignValue<T>;
  // One node's rates and net production, with its state `own` and the patch's field
  // `field` installed; the patch is left off its state.
  value_type node_rates_in_field(size_t node, const std::vector<value_type>& own,
                                 const std::vector<value_type>& field,
                                 double time, std::vector<value_type>& rates)
    requires NamesSignValue<T>;

  // Set where a run begins, read by compute_environment().
  bool keep_field = false;
  // The slot the rate evaluation now running stores into or loads from.
  solved_values* storing = nullptr;
  const solved_values* loading = nullptr;
  // See set_schedule(). Empty for a patch built by hand, which is given none.
  std::shared_ptr<const std::vector<schedule_entry>> schedule =
      std::make_shared<const std::vector<schedule_entry>>();
  // See set_run_birth_dates().
  std::vector<double> run_birth_dates;

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

  parameters_type parameters;

  double area;
  environment_type environment;
  std::vector<species_type> species;

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
  validate_ode_weights(control);

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
  // What a rebound patch reads, and nothing else. patch_type and
  // max_patch_lifetime are read by validate(), which the constructor below runs
  // to rebuild the disturbance regime: without them the rebound patch gets a
  // different one, pr_patch_survival moves, and the fecundity rates differ with
  // nothing raised. node_schedule_times is length-checked there too.
  //
  // Left behind: ode_times and ode_step_sizes, which are a record of the last
  // run rather than configuration and are read only by make_node_schedule, and
  // n_patches, which nothing reads. A rebind happens once per width a sweep
  // crosses, so copying two vectors of the whole trajectory into it was a copy
  // per width per gradient, into an object that never looked at them.
  Parameters<T2,E2> p2;
  p2.patch_area = parameters.patch_area;
  p2.patch_type = parameters.patch_type;
  p2.max_patch_lifetime = parameters.max_patch_lifetime;
  p2.node_schedule_times_default = parameters.node_schedule_times_default;
  p2.node_schedule_times = parameters.node_schedule_times;
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
  // The entries the rebound patch applies, which the constructor does not take.
  out.schedule = schedule;
  // The field and the boundary node are left to the caller. Every caller
  // evaluates before reading either, which rebuilds the field and re-evaluates
  // the inflow condition, so computing them here would solve the boundary leaf
  // twice and discard it. A caller that reads before evaluating gets an unbuilt
  // field.
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
void Patch<T,E>::reset() {
   for (auto& s : species) {
    s.clear();
    // allocate variables for tracking resource consumption
    s.resize_consumption_rates(environment.n_resources());
  }

  // compute ephemeral effects like light_availability
  environment.clear();

  if (!parameters.initial_state.empty()) {
    // Seed the patch from an exported state / initial size distribution.
    // set_initial_state() does its own compute_environment()/compute_rates at
    // the real node heights, so skip the empty-patch path.
    set_initial_state();
    check_initial_density_rates();
  } else {
    compute_environment();

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
// environment and rates. We load the state directly, rather than via the
// double-arg set_ode_state, so the first environment build reads node heights
// that exist.
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

  // Build the environment from the real node heights and compute rates for the
  // seeded population.
  compute_environment();
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
                 "introduction time (" + util::to_string(repeated) + "), but the "
                 "birth-date size-density coordinate needs each node's own birth "
                 "date. Remove the repeat from the node schedule, supply per-node "
                 "introduction times (parameters$initial_node_times) with an "
                 "initial state, or run with "
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
  // Double, not the working scalar: this reads the state to test it for
  // finiteness and never differentiates it, and the iterator write converts on
  // the way out. Held at the working scalar it was a slot per entry per stage
  // for a value flattened on the next line.
  std::vector<double> env_state(environment.ode_size());
  environment.ode_state(env_state.begin());
  for (size_t i = 0; i < env_state.size(); ++i) {
    const double state_i = env_state[i];
    if (!util::is_finite(state_i)) {
      // A rejection, not a failure: measured, every observed TF24 soil
      // excursion is explicit-integrator overshoot across a discrete change and
      // not a defect in the balance -- hand it to the stepper to shrink and
      // retry rather than failing the run. Mode (1) above stays fatal.
      odelia::util::stop_domain("Non-finite environment state (index " + util::to_string(i) +
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

// The pair's first entry, which is what the pair is built to be. Nothing on the
// field's path arrives here -- that takes the split reduction and closes it -- so
// this is an accessor, and an accessor paying one extra crown evaluation per node
// is better than a second loop that can disagree with the one below it.
template <typename T, typename E>
typename Patch<T,E>::value_type
Patch<T,E>::compute_competition(double height) const {
  return compute_competition_and_slope(height).value;
}

template <typename T, typename E>
with_slope<typename Patch<T,E>::value_type>
Patch<T,E>::compute_competition_and_slope(double z) const {
  with_slope<value_type> sum{0.0, 0.0};
  for (size_t i = 0; i < species.size(); ++i) {
    const with_slope<value_type> fs =
      species[i].compute_competition_and_slope(z);
    sum.value += fs.value / area;
    sum.slope += fs.slope / area;
  }
  return sum;
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
  return odelia::util::to_passive(
      species[species_index].net_reproduction_ratio(scalars));
}

// Offspring production, equal to overall fitness scaled by the birth rate.
template <typename T, typename E>
std::vector<double> Patch<T,E>::offspring_production() const {
  auto ret = std::vector<double>(species.size());
  for (size_t i = 0; i < species.size(); ++i) {
    ret[i] = odelia::util::to_passive(species[i].offspring_production());
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
// An invader replaying a recorded field reaches this with that field already
// installed, which is what makes its own inflow condition right to evaluate
// here: the condition is what a species does IN the field, not to it.
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
void Patch<T,E>::compute_environment() {
  if (size() == 0) {
    return;
  }
  // The boundary node is one end of the birth-date quadrature and its birth date
  // is the current time, so it is stamped before either profile is built. Its own
  // stamp is set in compute_rates(), which the stepper calls *after* the
  // set_ode_state() that brings us here, so reading that stamp would use the
  // previous derivs call's time and shorten the boundary interval by a
  // Runge-Kutta stage. The measured effect is below 1e-6, but the interval is
  // then a function of the step size, which the spatial quadrature has no
  // business depending on, and the whole integral *is* that one interval while a
  // species has a single node. No-op for the height coordinate, where this
  // abscissa is the constant initial height.
  for (auto& s : species) {
    s.set_new_node_birth_date(environment.time);
  }

  // An invader is evaluated in a field it did not build: the one its slot holds.
  // It adds nothing to the reduction below, and only the inflow condition is its
  // own.
  //
  // ⚠️ THE INSTANT IS CHECKED, to every digit: a recording paired one step out
  // would hand every evaluation its neighbour's field with every number finite.
  const recorded_field* field =
      storing != nullptr ? storing->field.get()
      : loading != nullptr ? loading->field.get() : nullptr;
  if (field != nullptr) {
    if (!util::identical(field->time, environment.time)) {
      std::ostringstream m;
      m.precision(17);
      m << "The recorded field is for t=" << field->time
        << " and this evaluation is at t=" << environment.time
        << ": the recording and the run it is replayed through disagree about "
           "which rate evaluation this is";
      util::stop(m.str());
    }
    environment.set_knot_data(field->knot_data);
    environment.set_ode_state(field->state.begin());
    compute_boundary_nodes();
    return;
  }

  compute_environment_excl_capturing();
  compute_boundary_nodes();
  compute_environment_closing();

  // Kept here and nowhere else, because this is where the field comes to exist.
  if (storing != nullptr && keep_field) {
    recorded_field kept;
    kept.knot_data = environment.knot_data();
    kept.state.assign(environment.ode_size(), 0.0);
    environment.ode_state(kept.state.begin());
    kept.time = environment.time;
    storing->built = std::make_shared<const recorded_field>(std::move(kept));
  }
}

// The field without the boundary interval, keeping each species' reduction at
// the point its closing trapezium would be added.
template <typename T, typename E>
void Patch<T,E>::compute_environment_excl_capturing() {
  // The reduction fills one entry per knot, so the only thing needed here is one
  // vector per species; field_splits sizes each of them.
  competition_capture.resize(size());
  auto f = [&](const std::vector<value_type>& x, std::vector<value_type>& y,
               std::vector<value_type>& m) -> void {
    for (size_t k = 0; k < x.size(); ++k) {
      y[k] = 0.0;
      m[k] = 0.0;
    }
    // Species-major, so each knot's sum still runs over the species in index
    // order: the reduction is one pass over a species' nodes for the whole knot
    // set, and asking it per knot is what made the build quadratic.
    for (size_t i = 0; i < size(); ++i) {
      species[i].field_splits(x, competition_capture[i]);
      for (size_t k = 0; k < x.size(); ++k) {
        const with_slope<value_type> open =
          competition_capture[i][k].without_boundary;
        y[k] += open.value / area;
        m[k] += open.slope / area;
      }
    }
  };
  if (size() > 0) {
    environment.compute_environment(f, height_max());
  }
}

// The same field with the boundary interval closed, from the kept reductions and
// the boundary node the call between the two established.
template <typename T, typename E>
void Patch<T,E>::compute_environment_closing() {
  auto f = [&](const std::vector<value_type>& x, std::vector<value_type>& y,
               std::vector<value_type>& m) -> void {
    for (size_t k = 0; k < x.size(); ++k) {
      y[k] = 0.0;
      m[k] = 0.0;
    }
    for (size_t i = 0; i < species.size(); ++i) {
      if (competition_capture[i].size() != x.size()) {
        util::stop("compute_environment_closing: the field without the "
                   "boundary interval was built over a different knot set");
      }
      for (size_t k = 0; k < x.size(); ++k) {
        const with_slope<value_type> fs =
          species[i].close_competition_and_slope(competition_capture[i][k], x[k]);
        y[k] += fs.value / area;
        m[k] += fs.slope / area;
      }
    }
  };
  if (size() > 0) {
    environment.compute_environment(f, height_max());
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

  // Produced here and drained by the call below, so it lives here. As a member
  // it needed clearing by hand, and a refusal thrown between the fill and the
  // clear left active values on a Patch across a tape reset.
  std::vector<value_type> resource_depletion;
  resource_depletion.reserve(env.n_resources());
  for (size_t i = 0; i < env.n_resources(); i++) {
    value_type resource_consumed = std::accumulate(
      species.begin(), species.end(), value_type(0.0),
      [i](const value_type& r, const species_type& s) -> value_type {
        return r + s.consumption_rate(i);  // accumulates r from zero
      });
    resource_depletion.push_back(resource_consumed / area);
  }
  env.compute_rates(resource_depletion);
}

// The whole light environment is rebuilt here, where only the knots below the
// seedling's height change. The knot fractions are fixed, so a narrower rebuild
// would write a subrange of the same positions.
template <typename T, typename E>
void Patch<T,E>::introduce_nodes(const introduction& species_index, double time) {
  push_nodes(species_index, time);

  // A schedule carrying the same time twice for one species stamps two nodes
  // with it, and the grid both reductions integrate over is those times.
  check_birth_dates_distinct();

  compute_environment();

  // New nodes have just changed the state and the light field, so the stored
  // rates now describe neither. The solver reads them next without checking, so
  // they have to be brought up to date here.
  compute_rates();
}

// The numbers a node carries that the ODE state does not: its birth date, and
// the patch density, survival and birth rate there. All are functions of the
// time it is introduced at, which is what lets a record hold only the time.
template <typename T, typename E>
void Patch<T,E>::push_nodes(const introduction& species_index, double time) {
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
template <typename Select>
std::pair<size_t, double>
Patch<T,E>::scale_node_densities(size_t species_index, Select select) {
  species_type& sp = species[species_index];
  size_t n_affected = 0;
  double density_removed = 0.0;
  for (auto n = sp.node_begin(); n != sp.node_end(); ++n) {
    const double phi = select(odelia::util::to_passive(n->height()));
    if (phi >= 1.0) {
      continue;
    }
    if (!(phi > 0.0)) {
      util::stop("An event must leave a positive fraction of each cohort: "
                 "removing all of one would take its density to -Inf, which "
                 "the density transport cannot carry back.");
    }
    const double log_phi = std::log(phi);
    const double before = odelia::util::to_passive(n->get_density());
    n->set_log_density(n->get_log_density() + log_phi);
    if (util::is_finite(before)) {
      density_removed += before - odelia::util::to_passive(n->get_density());
    }
    n->individual.set_state("mortality",
                            n->individual.state(MORTALITY_INDEX) - log_phi);
    ++n_affected;
  }
  return {n_affected, density_removed};
}

template <typename T, typename E>
EventRecord Patch<T,E>::apply_event(const NodeScheduleEvent& event) {
  EventRecord rec;
  rec.time = time();
  rec.type = event.type;
  rec.target = event.target;
  rec.target_index = event.target_index;
  rec.requested = event.params;

  // Which species the event reaches: one, or all of them.
  std::vector<size_t> targets;
  if (event.target == EventTarget::Species) {
    targets.push_back(event.target_index);
  } else {
    for (size_t i = 0; i < species.size(); ++i) {
      targets.push_back(i);
    }
  }

  switch (event.type) {
  case EventType::NodeIntroduction:
    // Introductions are batched by the caller so that a run of them costs one
    // environment recompute rather than one each; they never arrive here.
    util::stop("Node introductions are applied via introduce_new_nodes()");
    break;

  case EventType::ResourcePulse: {
    // The environment only: no individual changes, so the competition profile
    // is untouched and the nodes keep their state. The rates are stale
    // afterwards, but the solver recomputes them when it re-reads the system
    // (Patch::ode_rates computes), so there is nothing to do here.
    if (event.target_index >= environment.n_resources()) {
      util::stop("Resource " + util::to_string(event.target_index + 1) +
                 " does not exist: this environment has " +
                 util::to_string(environment.n_resources()) + " resources");
    }
    rec.applied = environment.add_resource_pulse(event.target_index,
                                                 event.params.at(0));
    break;
  }

  case EventType::Harvest: {
    // Remove a fraction of the individuals in a size band. One rule covers the
    // cases that matter: a whole-population removal (the default band), taking
    // everything above a size, and the size-selective removal the restoration
    // work asks for (#627).
    //
    // Individuals are removed rather than shrunk. Changing sizes instead would
    // have to keep them in decreasing order, which Species' competition
    // integral relies on, and the unordered fallback exists only on the size
    // coordinate. Removing part of an individual rather than the whole of it
    // is a real and different thing, and needs that question answered first.
    const double fraction = event.params.at(0);
    const double size_min = event.params.at(1);
    const double size_max = event.params.at(2);
    if (!(fraction >= 0.0) || fraction >= 1.0) {
      util::stop("Harvest fraction must be in [0, 1)");
    }
    const double phi = 1.0 - fraction;
    size_t n_affected = 0;
    double removed = 0.0;
    for (size_t i : targets) {
      const auto r = scale_node_densities(i, [&](double size) {
        return (size >= size_min && size <= size_max) ? phi : 1.0;
      });
      n_affected += r.first;
      removed += r.second;
    }
    rec.applied = {fraction, static_cast<double>(n_affected), removed};
    break;
  }

  case EventType::ClimateExtreme: {
    // An episode of extreme conditions -- heat, cold, salinity, whatever the
    // model's `intensity` means -- that kills in proportion to the dose
    // accumulated above a threshold.
    //
    // This is a worked example of the sub-integration pattern rather than
    // defensible physiology. The event has a nominal duration in the world;
    // the solver's clock does not move across it, so the action steps its own
    // damage variable over that duration under a simple daily cycle and hands
    // back a single survival fraction. The hook is the part worth keeping: a
    // model with a real damage state (#566 for TF24's leaf) plugs in here, and
    // until one exists a dose-dependent mortality is the honest crude stand-in.
    const double intensity = event.params.at(0);
    const double duration = event.params.at(1);
    const double threshold = event.params.at(2);
    const double sensitivity = event.params.at(3);
    if (!util::is_finite(duration) || duration < 0.0) {
      util::stop("Climate extreme duration must be finite and non-negative");
    }

    // Half-hourly sub-steps through a sinusoidal daily cycle peaking at
    // `intensity`, with amplitude set by its excess over `threshold`. Dose
    // accrues only above the threshold, so a mild episode accrues none at all.
    const double dt_days = 0.5 / 24.0;
    const size_t n_steps =
      static_cast<size_t>(std::max(1.0, std::ceil(duration * 365.0 / dt_days)));
    const double amplitude = std::max(0.0, intensity - threshold);
    double damage = 0.0;
    for (size_t k = 0; k < n_steps; ++k) {
      const double day_fraction =
        std::fmod(static_cast<double>(k) * dt_days, 1.0);
      // Peak mid-cycle; the trough sits one amplitude below the threshold.
      const double level =
        threshold + amplitude * (2.0 * std::sin(M_PI * day_fraction) - 1.0);
      const double excess = std::max(0.0, level - threshold);
      damage += sensitivity * excess * dt_days / 365.0;
    }

    const double phi = std::exp(-damage);
    size_t n_affected = 0;
    double removed = 0.0;
    for (size_t i : targets) {
      const auto r = scale_node_densities(i, [&](double) { return phi; });
      n_affected += r.first;
      removed += r.second;
    }
    rec.applied = {1.0 - phi, static_cast<double>(n_affected), removed};
    break;
  }
  }
  return rec;
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
  compute_environment();

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

// Every species in a patch carries the same strategy type, so this is a count
// per species rather than a walk that builds a vector of pointers to measure its
// length. It is asked several times per gradient.
template <typename T, typename E>
size_t Patch<T,E>::trait_adjoint_size() const {
  return species.size() * T::ad_column_count;
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
size_t Patch<T,E>::nodes_at(size_t i, double time) const {
  size_t n = parameters.initial_state.empty() ? 0
                                              : parameters.n_initial_cohorts.at(i);
  for (const schedule_entry& e : *schedule) {
    if (e.time < time &&
        std::find(e.species.begin(), e.species.end(), i) != e.species.end()) {
      ++n;
    }
  }
  return n;
}

// Derived rather than replayed: the schedule this patch is run from says which
// species gain a node when, so a walk that jumps into the middle of a recording
// works out the shape there instead of being handed a list of what happened.
template <typename T, typename E>
void Patch<T,E>::reshape_to(double time) {
  bool moved = false;
  for (size_t i = 0; i < species.size(); ++i) {
    const size_t target = nodes_at(i, time);
    while (species[i].size() > target) {
      species[i].remove_newest_node();
      moved = true;
    }
    if (species[i].size() < target) {
      const size_t base = parameters.initial_state.empty()
                              ? 0
                              : parameters.n_initial_cohorts.at(i);
      // The values a pushed node carries are the state's and arrive with the
      // load; only its stamps are its own, and those are the schedule's date.
      // Species i's j-th introduction is its (base + j)-th node.
      size_t node = base;
      const introduction one{i};
      for (const schedule_entry& e : *schedule) {
        if (species[i].size() == target) {
          break;
        }
        if (std::find(e.species.begin(), e.species.end(), i) == e.species.end()) {
          continue;
        }
        if (node++ == species[i].size()) {
          push_nodes(one, e.time);
          moved = true;
        }
      }
    }
  }
  // The field is built by the evaluation that follows, at the recorded time.
  if (moved) {
    check_birth_dates_distinct();
  }
}

template <typename T, typename E>
std::vector<EventRecord> Patch<T,E>::apply_insertion(double time) {
  std::vector<EventRecord> ret;
  const auto at = std::find_if(
      schedule->begin(), schedule->end(),
      [time](const schedule_entry& e) { return util::identical(e.time, time); });
  if (at == schedule->end()) {
    return ret;
  }
  for (const NodeScheduleEvent& a : at->actions) {
    ret.push_back(apply_event(a));
  }
  push_nodes(at->species, time);
  // A schedule carrying the same time twice for one species stamps two nodes
  // with it, and the grid both reductions integrate over is those times.
  check_birth_dates_distinct();
  return ret;
}

template <typename T, typename E>
template <typename It>
void Patch<T,E>::introduce_from(const introduction& species_index, double time,
                                It x, std::vector<value_type>& y) {
  std::vector<value_type> below(x, x + static_cast<std::ptrdiff_t>(ode_size()));
  std::vector<value_type> rates(below.size());
  odelia::ode::derivs(*this, below, rates, time);
  push_nodes(species_index, time);
  y.assign(ode_size(), value_type(0.0));
  ode_state(y.begin());
}

// One tangent seed per input column, over the same map. The node the map pushes is
// removed between columns, so every column is taken from the same width.
template <typename T, typename E>
std::vector<std::vector<double>>
Patch<T,E>::introduction_jacobian(const std::vector<size_t>& species_index,
                                  const std::vector<double>& state_before,
                                  double time_before) {
  const size_t n_state = ode_size();
  const size_t n_trait = trait_adjoint_size();
  size_t n_out = n_state;
  for (size_t i : species_index) {
    n_out += species.at(i).r_new_node().ode_size();
  }
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
      seed_direction(x[j], (j == c) ? 1.0 : 0.0);
    }
    // The traits first, for the reason the transpose writes them first: a
    // quantity the state determines reads them while deriving it.
    size_t at = n_state;
    for (tangent* p : active.ad_parameters()) {
      *p = x[at++];
    }
    util::check_length(at, x.size());
    std::vector<tangent> y(n_out);
    active.introduce_from(species_index, time_before, x.begin(), y);
    for (size_t r = 0; r < n_out; ++r) {
      ret[r][c] = derivative_along(y[r]);
    }
    for (size_t i : species_index) {
      active.species[i].remove_newest_node();
    }
  }
  return ret;
}

template <typename T, typename E>
bool Patch<T,E>::ode_state_valid(const std::vector<double>& y) const {
  // The environment block first, and separately: it is the trailing part of the
  // ODE vector, it is where integrator overshoot shows up (TF24's soil water,
  // driven past its bounds by a step inherited across a discrete change), and
  // the completed step's state never reaches check_finite_ode_state(), which
  // runs per stage. Finiteness only -- a node's log_density is legitimately
  // -Inf, so this must not be extended over the species block.
  const size_t n_env = environment.ode_size();
  if (n_env > 0 && y.size() >= n_env) {
    for (size_t i = y.size() - n_env; i < y.size(); ++i) {
      if (!util::is_finite(y[i])) {
        return false;
      }
    }
  }

  // Resolved on the concrete strategy, so a model that appends or inserts a
  // state gets its own positions rather than its base's. Named rather than
  // indexed because a strategy that declared the index would be declaring
  // something its own state_names() already says.
  static const std::vector<size_t> bounded = [] {
    const std::vector<std::string> names = T::state_names();
    std::vector<size_t> ret;
    for (const std::string& name : T::non_negative_states()) {
      const auto found = std::find(names.begin(), names.end(), name);
      if (found == names.end()) {
        util::stop("Strategy declares '" + name + "' as a non-negative state "
                   "and does not carry it");
      }
      ret.push_back(static_cast<size_t>(found - names.begin()));
    }
    return ret;
  }();
  if (bounded.empty()) {
    return true;
  }
  size_t at = 0;
  for (const auto& sp : species) {
    for (auto n = sp.node_begin(); n != sp.node_end(); at += n->ode_size(), ++n) {
      for (const size_t i : bounded) {
        if (y[at + i] < 0.0) {
          return false;
        }
      }
    }
  }
  return true;
}

template <typename T, typename E>
void Patch<T,E>::sign_values(std::vector<double>& out) const
  requires NamesSignValue<T> {
  out.clear();
  for (const species_type& s : species) {
    const int at = s.strategy_ptr()->sign_value_aux();
    for (auto n = s.node_begin(); n != s.node_end(); ++n) {
      out.push_back(n->individual.aux(at));
    }
  }
}

template <typename T, typename E>
template <class Step>
bool Patch<T,E>::split_sign_changes(const Step& step, split_samples& samples,
                                    split_blocks& record)
  requires NamesSignValue<T> {
  if (!control.ode_split_sign_changes) {
    return false;
  }
  field_samples field;
  bool sampled = false;
  size_t node = 0, first = 0;
  for (species_type& s : species) {
    for (size_t j = 0; j < s.size(); ++j, ++node) {
      const size_t width =
        (s.node_begin() + static_cast<std::ptrdiff_t>(j))->ode_size();
      auto changes = step.template sign_changes<solved_values>(
        node, [&](double u, solved_values* into) -> double {
          if (!sampled) {
            sample_field(step, field, samples);
            sampled = true;
          }
          const double time = step.time + u * step.h;
          if (into == nullptr) {
            return node_value_at(step, field, node, first, width, u, time);
          }
          const odelia::ode::solved_scope<Patch, solved_values> extent{*this,
                                                                       *into};
          return node_value_at(step, field, node, first, width, u, time);
        });
      if (!changes.empty()) {
        odelia::ode::split_block<solved_values>& block = record.emplace_back();
        block.block = node;
        block.first = first;
        std::vector<double> split_at;
        for (const auto& change : changes) {
          split_at.push_back(change.u);
        }
        block.sign_changes = std::move(changes);
        split_node(step, field, node, first, width, split_at, [&] {
          return odelia::ode::solved_scope<Patch, solved_values>{
            *this, block.solved.emplace_back()};
        });
      }
      first += width;
    }
  }
  return sampled;
}

template <typename T, typename E>
template <class Step>
void Patch<T,E>::split_as_recorded(const Step& step,
                                   const split_samples& samples,
                                   const split_blocks& recorded)
  requires NamesSignValue<T> {
  field_samples field;
  sample_field(step, field, samples);
  for (const auto& block : recorded) {
    const size_t width = node_width(block.block);
    std::vector<value_type> split_at;
    for (const auto& change : block.sign_changes) {
      // A zero slope gives no derivative: the sign change stays where the run found it.
      if (change.slope == 0.0) {
        split_at.emplace_back(change.u);
        continue;
      }
      split_at.push_back(odelia::implicit_value<value_type>(
        change.u, change.slope, [&](const value_type& u) -> value_type {
          const odelia::ode::solved_scope<Patch, const solved_values> extent{
            *this, change.solved};
          return node_value_at(step, field, block.block, block.first, width, u,
                               step.time + change.u * step.h);
        }));
    }
    size_t next = 0;
    split_node(step, field, block.block, block.first, width, split_at, [&] {
      return odelia::ode::solved_scope<Patch, const solved_values>{
        *this, block.solved.at(next++)};
    });
    util::check_length(next, block.solved.size());
  }
}

template <typename T, typename E>
template <class Step, class Row>
void Patch<T,E>::take_recorded_splits(const Step& step, const Row& recorded,
                                      split_samples& samples,
                                      split_blocks& record)
  requires NamesSignValue<T> {
  field_samples field;
  bool sampled = false;
  size_t node = 0, first = 0;
  for (species_type& s : species) {
    const size_t width = s.size() > 0 ? s.node_begin()->ode_size() : 0;
    for (const auto& run_block : recorded.solved.split_blocks) {
      // A block is the run's node, of its one species: only such a run keeps its
      // splits for a walk.
      const double born = run_birth_dates.at(run_block.block);
      const auto copy =
        std::find_if(s.node_begin(), s.node_end(), [&](const auto& n) {
          return util::identical(n.introduction_time(), born);
        });
      if (copy == s.node_end()) {
        continue;
      }
      if (!sampled) {
        sample_field(step, field, samples);
        sampled = true;
      }
      const size_t k = static_cast<size_t>(copy - s.node_begin());
      odelia::ode::split_block<solved_values>& block = record.emplace_back();
      block.block = node + k;
      block.first = first + k * width;
      std::vector<double> split_at;
      for (const auto& change : run_block.sign_changes) {
        // A zero slope holds the sign change where the run found it.
        block.sign_changes.push_back({change.u, 0.0, {}});
        split_at.push_back(change.u);
      }
      split_node(step, field, block.block, block.first, width, split_at, [&] {
        return odelia::ode::solved_scope<Patch, solved_values>{
          *this, block.solved.emplace_back()};
      });
    }
    node += s.size();
    first += s.size() * width;
  }
}

template <typename T, typename E>
template <class Step, class Samples>
void Patch<T,E>::sample_field(const Step& step, field_samples& field,
                              Samples& samples) {
  if constexpr (!std::is_const_v<Samples>) {
    samples.resize(field.size());
  }
  util::check_length(samples.size(), field.size());
  std::vector<value_type> state(ode_size());
  for (size_t m = 0; m < field.size(); ++m) {
    const double u = step.sample_fractions[m];
    step.dense_state(u, 0, state);
    {
      const odelia::ode::solved_scope<
        Patch, std::remove_reference_t<decltype(samples[m])>> extent{*this,
                                                                     samples[m]};
      set_ode_state(state.begin(), step.time + u * step.h);
    }
    field[m].resize(environment.light_availability.knot_data_size() +
                    environment.ode_size());
    environment.ode_state(
      environment.light_availability.knot_data(field[m].begin()));
  }
}

template <typename T, typename E>
size_t Patch<T,E>::node_width(size_t node) const {
  for (const species_type& s : species) {
    if (node < s.size()) {
      return (s.node_begin() + static_cast<std::ptrdiff_t>(node))->ode_size();
    }
    node -= s.size();
  }
  util::stop("node_width: the patch holds fewer nodes than the one named");
}

template <typename T, typename E>
template <class Step, class U>
typename Patch<T,E>::value_type
Patch<T,E>::node_value_at(const Step& step, const field_samples& field,
                          size_t node, size_t first, size_t width, const U& u,
                          double time)
  requires NamesSignValue<T> {
  std::vector<value_type> own(width), field_at, rates(width);
  step.dense_state(u, first, own);
  step.sample_at(u, field, field_at);
  return node_rates_in_field(node, own, field_at, time, rates);
}

template <typename T, typename E>
template <class Step, class U, class Solved>
void Patch<T,E>::split_node(const Step& step, const field_samples& field,
                            size_t node, size_t first, size_t width,
                            const std::vector<U>& split_at, Solved&& solved)
  requires NamesSignValue<T> {
  std::vector<value_type> field_at, own(width);
  step.integrate_pieces(
    first, split_at,
    [&](const U& u, const std::vector<value_type>& at,
        std::vector<value_type>& out) {
      step.sample_at(u, field, field_at);
      const auto extent = solved();
      // ⚠️ THE TIME IS HELD WHERE THE RUN EVALUATED: rates reading the time itself
      // miss their derivative in time times each sign change's move.
      node_rates_in_field(node, at, field_at,
                          step.time + odelia::util::to_passive(u) * step.h, out);
    },
    own);
  std::copy(own.begin(), own.end(),
            step.y_end.begin() + static_cast<std::ptrdiff_t>(first));
}

template <typename T, typename E>
typename Patch<T,E>::value_type
Patch<T,E>::node_rates_in_field(size_t node, const std::vector<value_type>& own,
                                const std::vector<value_type>& field,
                                double time, std::vector<value_type>& rates)
  requires NamesSignValue<T> {
  util::check_length(field.size(),
                     environment.light_availability.knot_data_size() +
                       environment.ode_size());
  environment.time = time;
  environment.set_ode_state(
    environment.light_availability.set_knot_data(field.begin()));
  for (species_type& s : species) {
    s.set_new_node_birth_date(time);
  }
  for (species_type& s : species) {
    if (node >= s.size()) {
      node -= s.size();
      continue;
    }
    const auto it = s.node_begin() + static_cast<std::ptrdiff_t>(node);
    util::check_length(own.size(), it->ode_size());
    util::check_length(rates.size(), it->ode_size());
    s.set_node_ode_state(node, own.begin());
    s.compute_node_rates(node, environment, survival_weighting->pr_survival(time),
                         s.extrinsic_drivers().evaluate("birth_rate", time));
    it->ode_rates(rates.begin());
    return it->individual.aux(s.strategy_ptr()->sign_value_aux());
  }
  util::stop("node_rates_in_field: the patch holds fewer nodes than the one named");
}

template <typename T, typename E>
void Patch<T,E>::error_weights(double time, std::vector<double>& w) const {
  // The environment's states close the patch's, as ode_state writes them.
  w.assign(ode_size(), 1.0);
  environment.error_weights(
    control, w.end() - static_cast<std::ptrdiff_t>(environment.ode_size()));
  const std::vector<double>& times = control.ode_weight_times;
  double factor = 1.0;
  if (!times.empty()) {
    const auto after = std::upper_bound(times.begin(), times.end(), time);
    factor = control.ode_weight_factors.at(
      static_cast<size_t>(after - times.begin()) - 1);
  }
  for (double& x : w) {
    x = std::min(x * factor, control.ode_weight_max);
  }
}

}

#endif
