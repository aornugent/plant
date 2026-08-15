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
  // Compute rates of change at the state currently loaded and save them into the
  // ode solver. Computing here rather than in set_ode_state is what makes the
  // rates always those of that state, however the caller arrived at it.
  template <typename It> It ode_rates(It it);
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

  // Everything one census metric's sweep accumulates. Several ride through the
  // transpose together so that a block is RECORDED once and swept once per
  // metric: the recording is a model evaluation and the sweep is arithmetic, so
  // where a metric costs a recording it costs a hundred times what it needs to.
  //
  // It exists as a bundle rather than four parallel arguments because the
  // batched forms would otherwise carry four vectors that must stay the same
  // length and in the same order, and a transposed pair there is silent.
  struct sweep_adjoints {
    std::vector<double> state;
    light_knot_adjoints knot;
    std::vector<boundary_node_adjoints> boundary_node;
  };

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

  // Every grid point both reductions integrate over: the introduced nodes and,
  // last within each species, its boundary node. The reductions run from the
  // boundary node up, so this is the index their transposes scatter over, and it
  // is wider than the ODE state by one node per species.
  size_t reduction_node_count() const;

  // The boundary node's own size and density adjoints, one entry per species. No
  // ODE row holds them, so boundary_condition_adjoint() pulls each back through
  // the condition that sets it. The height and leaf-area slots reach the seed
  // height, which solves its own condition, so they carry a trait row rather
  // than terminating.
  std::vector<boundary_node_adjoints> boundary_node_adjoint;

  // How often the boundary's own term was asked for, and how often it carried
  // anything. The boundary node stands at the seed's height for a whole run while
  // its condition is re-evaluated at every stage of every step, so this row acts
  // once per step where the quantity it belongs to changes once per plant. Its
  // value being right is therefore not the same claim as its being multiplied by
  // the right number of evaluations, and only the second is a count.
  //
  // `asked` counts every call and `carried` only the calls that recorded, which
  // differ whenever no boundary adjoint was seeded. A caller resets both and reads
  // them back around one sweep.
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

  // The knot adjoints the cohort blocks produced on the last sweep, kept so a
  // trait row can be read back as its block half and its reduction half rather
  // than only as their sum.
  light_knot_adjoints last_knot_adjoint{std::vector<double>(), std::vector<double>()};
  size_t trait_adjoint_size() const;
  // The same order, named. Each name carries its species index, because
  // concatenating the strategies' own names repeats every one of them per
  // species and character indexing then resolves each to species one's column,
  // which an unknown-name check cannot see.
  std::vector<std::string> trait_adjoint_names() const;
  void clear_trait_adjoint(size_t n_metrics = 1);

  // How large the recording grew on the last block, and how many blocks the one
  // tape has carried. Either climbing with the cohort count is the tape leaking.
  size_t block_recording_size = 0;
  size_t block_sweeps = 0;

  // One recording per cohort, on a tape held across the loop, swept once per
  // seed set. `seeds[m]` and `out[m]` are one census metric's, and the cost of a
  // metric past the first is its sweeps alone.
  void cohort_block_adjoint(const std::vector<block_seeds>& seeds,
                            std::vector<sweep_adjoints>& out);

  // The same for one metric, accumulating into this patch's own trait and
  // boundary-node adjoints. Every caller that is not the batched transpose takes
  // this form.
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
  // `metric` selects the row of trait_adjoint the newcomer's trait rows land
  // in, for allometry_adjoint's reason.
  void introduction_adjoint(const std::vector<size_t>& species_index,
                            const std::vector<double>& state_before,
                            double time_before,
                            const std::vector<double>& lambda_after,
                            std::vector<double>& lambda_before,
                            size_t metric = 0);

  // The soil drainage cascade and the water aggregation above it, which
  // between them write the soil state adjoints and the per-cohort uptake.
  // `boundary_out` takes the boundary node's own rows, which no ODE row holds.
  // It is a parameter rather than this patch's member because a batched sweep
  // fills one set per census metric.
  void soil_adjoint(const std::vector<double>& lambda_dydt,
                    std::vector<double>& lambda_state,
                    block_seeds& seeds,
                    std::vector<boundary_node_adjoints>& boundary_out);

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
  // `metric` selects which row of trait_adjoint the reduction's own trait rows
  // accumulate into, because a batched sweep fills one per census metric.
  void allometry_adjoint(const std::vector<node_size_adjoints>& sizes,
                         std::vector<double>& lambda_state,
                         std::vector<boundary_node_adjoints>& boundary_out,
                         size_t metric = 0);

  // The boundary node's density adjoints pulled back through the condition that
  // sets it, which is the whole of the inflow boundary's contribution: the adjoint
  // at the boundary times the boundary condition's own derivative. Both
  // evaluations are recorded, so each accumulator contracts against the one its
  // own forward pass read.
  //
  // Without this the accumulator is written and never read, and every parameter
  // reaching the census only through establishment comes back short while staying
  // finite and correctly signed.
  //
  // The two densities are the ones the accumulators were taken at, read at the
  // points the stage's own two evaluations were restored.
  //
  // Batched over metrics for cohort_block_adjoint's reason: the condition is a
  // physiology evaluation at the seed's size, so recording it once and sweeping
  // it per metric is what keeps a metric past the first cheap.
  void boundary_condition_adjoint(const std::vector<double>& density_in_field,
                                  const std::vector<double>& density_in_uptake,
                                  std::vector<sweep_adjoints>& out);

  // Each species' inflow boundary density, as the patch now holds it.
  std::vector<double> boundary_density() const;

  // The introduction as a map: the pre-introduction state and the traits in, the
  // whole widened state out. introduction_adjoint records it and
  // introduction_jacobian evaluates it at a tangent, so the transpose and its
  // reference differentiate one function rather than two spellings of one.
  //
  // `active` is left holding the node it pushed; a caller evaluating the map more
  // than once removes it between calls.
  template <class Active, class S>
  static void introduce_over(Active& active,
                             const std::vector<size_t>& species_index,
                             double time_before, size_t n_state,
                             const std::vector<S>& x, std::vector<S>& y);

  // That map's whole Jacobian, by forward tangent: one row per widened state entry
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
  out.compute_environment(false);
  // What a rebound patch owes its callers is the boundary node: each of them
  // reads r_new_node() before setting a state of its own, and none reads a
  // cohort's rates. Rating every cohort to reach the same new_node evaluates one
  // leaf per node, and this runs once per adjoint of the right-hand side.
  out.compute_boundary_nodes();
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
  compute_environment_once(rescale, false);
  compute_boundary_nodes();
  compute_environment_once(rescale, true);
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

template <typename T, typename E>
size_t Patch<T,E>::reduction_node_count() const {
  return node_count() + species.size();
}

// The drainage cascade is lower bidiagonal and the environment transposes it;
// the aggregation above it is one trapezium per species per resource.
template <typename T, typename E>
void Patch<T,E>::soil_adjoint(const std::vector<double>& lambda_dydt,
                              std::vector<double>& lambda_state,
                              block_seeds& seeds,
                              std::vector<boundary_node_adjoints>& boundary_out) {
  const size_t n_env = environment.ode_size();
  const size_t env_offset = ode_size() - n_env;
  const size_t n_resource = environment.n_resources();
  util::check_length(lambda_dydt.size(), ode_size());
  util::check_length(lambda_state.size(), ode_size());
  util::check_length(seeds.uptake.size(), reduction_node_count() * n_resource);

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
  std::vector<node_uptake_adjoints> per_node(reduction_node_count(),
                                             node_uptake_adjoints{0, 0, 0});
  if (boundary_node_adjoint.size() != species.size()) {
    boundary_node_adjoint.assign(species.size(),
                                 boundary_node_adjoints{0, 0, 0, 0, 0});
  }
  for (size_t r = 0; r < n_resource; ++r) {
    if (lambda_uptake[r] == 0.0) {
      continue;
    }
    for (size_t k = 0; k < per_node.size(); ++k) {
      per_node[k] = node_uptake_adjoints{0, 0, 0};
    }
    size_t at = 0;
    for (size_t i = 0; i < species.size(); ++i) {
      species[i].consumption_rate_adjoint(static_cast<int>(r),
                                          lambda_uptake[r] / area,
                                          per_node.data() + at);
      at += species[i].size() + 1;
    }
    // The boundary node's slot is last within its species and owns no ODE row,
    // so its uptake seeds a block sweep and its two size rows are read out.
    size_t slot = 0, state_at = 0;
    for (size_t i = 0; i < species.size(); ++i) {
      for (size_t j = 0; j <= species[i].size(); ++j, ++slot) {
        seeds.uptake[slot * n_resource + r] += per_node[slot].uptake;
        if (j == species[i].size()) {
          boundary_out[i].height += per_node[slot].height;
          boundary_out[i].density_in_uptake += per_node[slot].log_density;
          continue;
        }
        lambda_state[state_at * node_stride + HEIGHT_INDEX] +=
          per_node[slot].height;
        lambda_state[state_at * node_stride + T::state_size() + 1] +=
          per_node[slot].log_density;
        ++state_at;
      }
    }
  }
}

template <typename T, typename E>
void Patch<T,E>::offspring_adjoint(const std::vector<double>& lambda_dydt,
                                   std::vector<double>& lambda_state,
                                   block_seeds& seeds) const {
  const size_t node_stride = node_type::ode_size();
  const double pr_patch_survival = survival_weighting->pr_survival(time());
  // Two index spaces, and they are not the same: `slot` runs over the grid points
  // the block loop sweeps, `state_at` over the ODE rows. A boundary node has the
  // first and not the second.
  size_t slot = 0, state_at = 0;
  for (size_t i = 0; i < species.size(); ++i) {
    for (size_t j = 0; j <= species[i].size(); ++j, ++slot) {
      if (j == species[i].size()) {
        continue;
      }
      const node_type& node = species[i].node_at(j);
      const double lambda_offspring =
        lambda_dydt[state_at * node_stride + T::state_size()];
      if (lambda_offspring != 0.0) {
        const double weight = odelia::util::to_passive(
          node.offspring_dt_dfecundity_rate(pr_patch_survival));
        seeds.rate[slot * T::state_size() + FECUNDITY_INDEX] +=
          lambda_offspring * weight;
        // exp(-mortality) puts the offspring rate on a state as well as a rate,
        // and the squashed non-finite case carries a zero survival with it.
        const double fecundity_rate = odelia::util::to_passive(
          node.individual.rate(FECUNDITY_INDEX));
        lambda_state[state_at * node_stride + MORTALITY_INDEX] -=
          lambda_offspring * fecundity_rate * weight;
      }
      ++state_at;
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

// The block: Individual::compute_rates at the active scalar, recorded once per
// cohort and swept once, with the leaf held constant at its declared boundary.
template <typename T, typename E>
void Patch<T,E>::cohort_block_adjoint(const std::vector<block_seeds>& seeds,
                                      std::vector<sweep_adjoints>& out) {
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
    std::vector<double> in;
    std::vector<std::vector<double>> out_adjoint, in_adjoint;
    std::vector<active_strategy> strategy_template;
    active_environment environment_template;
  };

  const size_t n_state = T::state_size();
  const size_t n_resource = environment.n_resources();
  const size_t n_knot = environment.light_availability.knot_count();
  const size_t n_layer = static_cast<size_t>(environment.get_soil_number_of_depths());
  const size_t node_stride = node_type::ode_size();
  const size_t env_offset = ode_size() - environment.ode_size();
  const size_t n_seed = seeds.size();
  if (n_seed == 0) {
    util::stop("cohort_block_adjoint: needs at least one seed set");
  }
  util::check_length(out.size(), n_seed);
  for (size_t m = 0; m < n_seed; ++m) {
    util::check_length(out[m].state.size(), ode_size());
    util::check_length(seeds[m].rate.size(), reduction_node_count() * n_state);
    util::check_length(seeds[m].transport.size(), reduction_node_count());
    util::check_length(seeds[m].uptake.size(), reduction_node_count() * n_resource);
    util::check_length(out[m].knot.value.size(), n_knot);
    util::check_length(out[m].knot.slope.size(), n_knot);
    if (out[m].boundary_node.size() != species.size()) {
      out[m].boundary_node.assign(species.size(),
                                  boundary_node_adjoints{0, 0, 0, 0, 0});
    }
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

  // One block per grid point of the reductions, so the boundary node's own draw
  // and its own radiation reach the soil state and the traits. Its rate seeds are
  // zero -- its rates are not ODE rates -- and its state adjoints have no row.
  size_t k = 0;
  size_t state_at = 0;
  size_t trait_at = 0;
  for (size_t i = 0; i < species.size(); ++i) {
    const size_t n_trait = species[i].strategy_ptr()->ad_parameters().size();
    for (size_t j = 0; j <= species[i].size(); ++j, ++k) {
      const bool boundary = j == species[i].size();
      const node_type& node =
        boundary ? species[i].r_new_node() : species[i].node_at(j);
      typename active_strategy::ptr strategy =
        std::make_shared<active_strategy>(ws.strategy_template[i]);
      active_individual individual(strategy);
      active_environment block_environment = ws.environment_template;

      ws.in.resize(individual.block_input_size(block_environment));
      node.individual.block_inputs(ws.in.begin(), environment);

      ws.out_adjoint.assign(n_seed,
                            std::vector<double>(n_state + 1 + n_resource, 0.0));
      for (size_t m = 0; m < n_seed; ++m) {
        for (size_t s = 0; s < n_state; ++s) {
          ws.out_adjoint[m][s] = seeds[m].rate[k * n_state + s];
        }
        ws.out_adjoint[m][n_state] = seeds[m].transport[k];
        for (size_t r = 0; r < n_resource; ++r) {
          ws.out_adjoint[m][n_state + 1 + r] = seeds[m].uptake[k * n_resource + r];
        }
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
      // ONE recording of the block, swept once per metric. The block is a model
      // evaluation and a sweep is arithmetic, so a metric past the first is
      // nearly free -- which is the whole of this function's economy.
      block_recording_size = odelia::ode::vector_jacobian_products(
        ws.tape, ws.in, ws.out_adjoint, block, ws.in_adjoint);
      ++block_sweeps;

      for (size_t m = 0; m < n_seed; ++m) {
        const std::vector<double>& row_in = ws.in_adjoint[m];
        size_t at = 0;
        for (size_t s = 0; s < n_state; ++s) {
          const double row = row_in[at++];
          if (boundary) {
            // No ODE row holds the boundary node's states: its height solves the
            // seed's own condition and the rest are the inflow condition's.
            if (s == HEIGHT_INDEX) {
              out[m].boundary_node[i].height += row;
            }
            continue;
          }
          out[m].state[state_at * node_stride + s] += row;
        }
        for (size_t c = 0; c < n_knot; ++c) {
          out[m].knot.value[c] += row_in[at++];
        }
        for (size_t c = 0; c < n_knot; ++c) {
          out[m].knot.slope[c] += row_in[at++];
        }
        for (size_t layer = 0; layer < n_layer; ++layer) {
          out[m].state[env_offset + layer] +=
            row_in[at++] *
            environment.dpsi_from_soil_moist_dtheta(
              odelia::util::to_passive(environment.vars.state(layer)), layer);
        }
        for (size_t p = 0; p < n_trait; ++p) {
          trait_adjoint[m][trait_at + p] += row_in[at++];
        }
        util::check_length(at, ws.in.size());
      }
      if (!boundary) {
        ++state_at;
      }
    }
    trait_at += n_trait;
  }
}

// One metric's sweep, accumulating into this patch's own trait and boundary-node
// adjoints. The batched form is the real one; this packs the members into a single
// seed set and unpacks them again, so a caller outside the transpose is unchanged.
template <typename T, typename E>
void Patch<T,E>::cohort_block_adjoint(const block_seeds& seeds,
                                      std::vector<double>& lambda_state,
                                      light_knot_adjoints& lambda_knot) {
  if (trait_adjoint.size() != 1 ||
      trait_adjoint[0].size() != trait_adjoint_size()) {
    clear_trait_adjoint();
  }
  if (boundary_node_adjoint.size() != species.size()) {
    boundary_node_adjoint.assign(species.size(),
                                 boundary_node_adjoints{0, 0, 0, 0, 0});
  }
  std::vector<sweep_adjoints> out(1);
  out[0].state = std::move(lambda_state);
  out[0].knot = std::move(lambda_knot);
  out[0].boundary_node = std::move(boundary_node_adjoint);
  cohort_block_adjoint(std::vector<block_seeds>{seeds}, out);
  lambda_state = std::move(out[0].state);
  lambda_knot = std::move(out[0].knot);
  boundary_node_adjoint = std::move(out[0].boundary_node);
}

// The traits go in before the state: area_leaf(height) reads lma, so a state set
// first is derived at the previous value. set_recorded_state is what computes the
// boundary node the run introduced, and introduce_new_node only pushes it.
template <typename T, typename E>
template <class Active, class S>
void Patch<T,E>::introduce_over(Active& active,
                                const std::vector<size_t>& species_index,
                                double time_before, size_t n_state,
                                const std::vector<S>& x, std::vector<S>& y) {
  size_t at = n_state;
  for (size_t i = 0; i < active.species.size(); ++i) {
    for (S* p : active.species[i].strategy_ptr()->ad_parameters()) {
      *p = x[at++];
    }
  }
  util::check_length(at, x.size());
  active.set_recorded_state(x.begin(), time_before);
  for (size_t i : species_index) {
    active.species[i].introduce_new_node();
  }
  util::check_length(y.size(), active.ode_size());
  active.ode_state(y.begin());
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
  for (size_t i = 0; i < species.size(); ++i) {
    for (const typename T::value_type* p :
         species[i].strategy_ptr()->ad_parameters()) {
      in.push_back(odelia::util::to_passive(*p));
    }
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
    std::vector<tangent> y(n_out);
    introduce_over(active, species_index, time_before, n_state, x, y);
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
void Patch<T,E>::introduction_adjoint(const std::vector<size_t>& species_index,
                                      const std::vector<double>& state_before,
                                      double time_before,
                                      const std::vector<double>& lambda_after,
                                      std::vector<double>& lambda_before,
                                      size_t metric) {
  using scalar = odelia::ode::active_scalar<double>;
  using active_strategy = typename at_scalar<scalar>::template apply<T>;
  const size_t n_trait = trait_adjoint_size();
  util::check_length(state_before.size(), ode_size());
  if (trait_adjoint.size() <= metric ||
      trait_adjoint[metric].size() != n_trait) {
    clear_trait_adjoint(metric + 1);
  }

  // The recording below emits the WHOLE widened state, so which widened row each
  // narrow row became is derived rather than written out. Copying the rows an
  // introduction does not touch and contracting only the newcomer's asserts that
  // the first group is an exact identity in the state and carries no trait row --
  // true of the nodes themselves, and an assertion about every other quantity a
  // state load rebuilds. Seeding the whole vector costs one more seed per row and
  // asserts nothing.
  lambda_before.assign(ode_size(), 0.0);

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

  const size_t n_state = ode_size();
  auto introduce = [&](const std::vector<scalar>& x,
                       std::vector<scalar>& y) -> void {
    introduce_over(active, species_index, time_before, n_state, x, y);
  };

  std::vector<double> in_adjoint;
  typename scalar::tape_type tape(false);
  odelia::ode::vector_jacobian_product(tape, in, lambda_after, introduce,
                                       in_adjoint);

  for (size_t j = 0; j < n_state; ++j) {
    lambda_before[j] = in_adjoint[j];
  }
  for (size_t p = 0; p < n_trait; ++p) {
    trait_adjoint[metric][p] += in_adjoint[n_state + p];
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
  util::check_length(out.size(), reduction_node_count());
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
      at += species[i].size() + 1;
    }
  }
}

// The size-space adjoint pulled back: to height through the allometry, and to
// the traits the reduction reads that no cohort step writes. A trait is one
// input read by every node of its species, so its row is a sum over them.
template <typename T, typename E>
void Patch<T,E>::allometry_adjoint(const std::vector<node_size_adjoints>& sizes,
                                   std::vector<double>& lambda_state,
                                   std::vector<boundary_node_adjoints>& boundary_out,
                                   size_t metric) {
  util::check_length(sizes.size(), reduction_node_count());
  util::check_length(lambda_state.size(), ode_size());
  if (trait_adjoint.size() <= metric ||
      trait_adjoint[metric].size() != trait_adjoint_size()) {
    clear_trait_adjoint(metric + 1);
  }
  if (boundary_node_adjoint.size() != species.size()) {
    boundary_node_adjoint.assign(species.size(),
                                 boundary_node_adjoints{0, 0, 0, 0, 0});
  }
  const size_t node_stride = node_type::ode_size();
  size_t k = 0;
  size_t state_at = 0;
  size_t trait_at = 0;
  for (size_t i = 0; i < species.size(); ++i) {
    typename T::ptr strategy = species[i].strategy_ptr();
    const size_t n_trait = strategy->ad_parameters().size();
    const typename T::light_reduction_slots slots =
      strategy->light_reduction_trait_slots();
    for (size_t j = 0; j <= species[i].size(); ++j, ++k) {
      const bool boundary = j == species[i].size();
      const individual_type& individual =
        (boundary ? species[i].r_new_node() : species[i].node_at(j)).individual;
      // The trait rows arrive for every grid point, because a trait is one input
      // read by every node of its species and the boundary node is one of them.
      strategy->light_reduction_trait_adjoint(
        slots, individual.area_leaf(), sizes[k].area_leaf,
        sizes[k].extinction, trait_adjoint[metric].data() + trait_at);
      if (boundary) {
        boundary_out[i].area_leaf += sizes[k].area_leaf;
        boundary_out[i].height += sizes[k].height;
        boundary_out[i].density_in_field += sizes[k].log_density;
        boundary_out[i].extinction += sizes[k].extinction;
        continue;
      }
      const double darea_leaf_dheight =
        odelia::util::to_passive(individual.darea_leaf_dheight());
      lambda_state[state_at * node_stride + HEIGHT_INDEX] +=
        sizes[k].height + sizes[k].area_leaf * darea_leaf_dheight;
      lambda_state[state_at * node_stride + T::state_size() + 1] +=
        sizes[k].log_density;
      ++state_at;
    }
    trait_at += n_trait;
  }
}

// The inflow boundary's two density adjoints, each contracted against the
// evaluation of the boundary condition its own forward pass read. The condition
// runs a whole physiology at the seed size, so it depends on the field, the soil
// and the traits; recording it and taking one vector-Jacobian product delivers
// every one of those rows rather than hand-writing them.
//
// The condition is recorded in the density and not in its logarithm, and the
// incoming adjoints are divided by the density to match. At a marginal recruit the
// log density and its trait sensitivity both diverge while the density and its own
// sensitivity both tend to zero, so the logarithmic form is a vanishing adjoint
// times a diverging row and the census gradient through a recruit that cannot pay
// for itself comes back as nothing rather than as the zero it is.
template <typename T, typename E>
std::vector<double> Patch<T,E>::boundary_density() const {
  std::vector<double> ret;
  ret.reserve(species.size());
  for (size_t i = 0; i < species.size(); ++i) {
    ret.push_back(
      odelia::util::to_passive(species[i].r_new_node().get_density()));
  }
  return ret;
}

template <typename T, typename E>
void Patch<T,E>::boundary_condition_adjoint(
    const std::vector<double>& density_in_field,
    const std::vector<double>& density_in_uptake,
    std::vector<sweep_adjoints>& out) {
  using scalar = odelia::ode::active_scalar<double>;
  using active_strategy = typename at_scalar<scalar>::template apply<T>;
  const size_t n_state = ode_size();
  const size_t n_trait = trait_adjoint_size();
  const size_t n_species = species.size();
  const size_t n_seed = out.size();
  if (n_seed == 0) {
    util::stop("boundary_condition_adjoint: needs at least one seed set");
  }
  for (size_t m = 0; m < n_seed; ++m) {
    util::check_length(out[m].state.size(), n_state);
    util::check_length(out[m].boundary_node.size(), n_species);
  }
  // The recording is the expensive part and a seed of zeros buys nothing. Every
  // channel the product carries has to be tested here, or a stand whose only
  // boundary sensitivity is the seed's size returns before recording anything.
  //
  // Batched, the test is over ALL the metrics: one recording serves them, so it
  // is worth making if any of them is seeded, and a metric that is not seeded
  // still gets its sweep -- a row of zeros, which is what its accumulator
  // expects. Testing per metric and recording per metric would give back the
  // cost this batching exists to remove.
  bool seeded = false;
  for (size_t m = 0; m < n_seed && !seeded; ++m) {
    for (const boundary_node_adjoints& b : out[m].boundary_node) {
      seeded = seeded || b.density_in_field != 0.0 || b.density_in_uptake != 0.0 ||
               b.height != 0.0 || b.area_leaf != 0.0;
    }
  }
  // Counted per METRIC, not per recording. One recording now serves every
  // metric, but the term still enters each metric's adjoint once per stage per
  // step, and it is that count a row is multiplied by -- which is what these
  // exist to assert (report 08 s4.8).
  boundary_condition_asked += n_seed;
  if (!seeded) {
    return;
  }
  boundary_condition_carried += n_seed;

  std::vector<value_type> current(n_state);
  ode_state(current.begin());
  std::vector<double> in(n_state);
  for (size_t j = 0; j < n_state; ++j) {
    in[j] = odelia::util::to_passive(current[j]);
  }
  in.reserve(n_state + n_trait);
  for (size_t i = 0; i < n_species; ++i) {
    for (const typename T::value_type* p :
         species[i].strategy_ptr()->ad_parameters()) {
      in.push_back(odelia::util::to_passive(*p));
    }
  }
  util::check_length(in.size(), n_state + n_trait);

  // Built with no tape active, so it holds no slot the recording could alias.
  Patch<active_strategy, typename at_scalar<scalar>::template apply<E>> active =
    this->template rebind_from<scalar>();

  const double time_ = environment.time;
  auto condition = [&](const std::vector<scalar>& x,
                       std::vector<scalar>& y) -> void {
    size_t at = n_state;
    for (size_t i = 0; i < active.species.size(); ++i) {
      for (scalar* p : active.species[i].strategy_ptr()->ad_parameters()) {
        *p = x[at++];
      }
    }
    // The traits go in before the state: the seed's leaf area reads them, so a
    // state set first derives it at the previous value.
    active.set_ode_state_and_field(x.begin(), time_);
    for (size_t i = 0; i < n_species; ++i) {
      y[i] = active.species[i].r_new_node().get_density();
    }
    active.compute_boundary_nodes();
    for (size_t i = 0; i < n_species; ++i) {
      y[n_species + i] = active.species[i].r_new_node().get_density();
      // The seed's height, which solves its own condition here rather than
      // arriving as a value, so this output is what carries dh_0/dtrait.
      y[2 * n_species + i] = active.species[i].r_new_node().height();
    }
  };

  // The reductions differentiate the log density, so each adjoint is divided by
  // the density it was taken at to seed the recording above. A zero density
  // carries a zero adjoint with it, and the quotient is the contribution's limit.
  util::check_length(density_in_field.size(), n_species);
  util::check_length(density_in_uptake.size(), n_species);
  std::vector<std::vector<double>> out_adjoint(
    n_seed, std::vector<double>(3 * n_species, 0.0));
  for (size_t m = 0; m < n_seed; ++m) {
  for (size_t i = 0; i < n_species; ++i) {
    if (density_in_field[i] > 0.0) {
      out_adjoint[m][i] =
        out[m].boundary_node[i].density_in_field / density_in_field[i];
    }
    if (density_in_uptake[i] > 0.0) {
      out_adjoint[m][n_species + i] =
        out[m].boundary_node[i].density_in_uptake / density_in_uptake[i];
    }
    // The leaf-area adjoint is converted to a height one and summed with it,
    // exactly as an interior node's is, because the seed's leaf area is the
    // allometry at the seed's height. Seeding leaf area as a second output
    // instead would deliver its partials at fixed height a second time, and the
    // reduction has already taken those.
    const double darea_leaf_dheight = odelia::util::to_passive(
      species[i].r_new_node().individual.darea_leaf_dheight());
    out_adjoint[m][2 * n_species + i] =
      out[m].boundary_node[i].height +
      out[m].boundary_node[i].area_leaf * darea_leaf_dheight;
  }
  }
  std::vector<std::vector<double>> in_adjoint;
  typename scalar::tape_type tape(false);
  odelia::ode::vector_jacobian_products(tape, in, out_adjoint, condition,
                                        in_adjoint);

  for (size_t m = 0; m < n_seed; ++m) {
    for (size_t j = 0; j < n_state; ++j) {
      out[m].state[j] += in_adjoint[m][j];
    }
    for (size_t p = 0; p < n_trait; ++p) {
      trait_adjoint[m][p] += in_adjoint[m][n_state + p];
    }
  }
}

template <typename T, typename E>
void Patch<T,E>::ode_rates_adjoint_batched(
    const std::vector<std::vector<double>>& lambda_dydt,
    std::vector<std::vector<double>>& lambda_y) {
  const size_t n = ode_size();
  const size_t n_resource = environment.n_resources();
  const size_t n_seed = lambda_dydt.size();
  if (n_seed == 0) {
    util::stop("ode_rates_adjoint: needs at least one rate adjoint");
  }
  std::vector<sweep_adjoints> out(n_seed);
  const size_t n_slot = reduction_node_count();
  const size_t n_knot = environment.light_availability.spline.knots().size();

  std::vector<block_seeds> seeds(
    n_seed, block_seeds{std::vector<double>(n_slot * T::state_size(), 0.0),
                        std::vector<double>(n_slot, 0.0),
                        std::vector<double>(n_slot * n_resource, 0.0)});
  for (size_t m = 0; m < n_seed; ++m) {
    util::check_length(lambda_dydt[m].size(), n);
    out[m].state.assign(n, 0.0);
    out[m].knot.value.assign(n_knot, 0.0);
    out[m].knot.slope.assign(n_knot, 0.0);
    out[m].boundary_node.assign(species.size(),
                                boundary_node_adjoints{0, 0, 0, 0, 0});
  }

  // The strategy rate adjoints the stage recursion supplies, and beside them the
  // transport term's, which is a block output rather than a closed-form seed.
  // A boundary node's slot takes none of them: its rates are not ODE rates.
  for (size_t m = 0; m < n_seed; ++m) {
    size_t slot = 0, state_at = 0;
    for (size_t i = 0; i < species.size(); ++i) {
      for (size_t j = 0; j <= species[i].size(); ++j, ++slot) {
        if (j == species[i].size()) {
          continue;
        }
        for (size_t s = 0; s < T::state_size(); ++s) {
          seeds[m].rate[slot * T::state_size() + s] =
            lambda_dydt[m][state_at * node_type::ode_size() + s];
        }
        seeds[m].transport[slot] =
          lambda_dydt[m][state_at * node_type::ode_size() + T::state_size() + 1];
        ++state_at;
      }
    }
  }

  // The stage evaluated the inflow condition twice and the two reductions read
  // different ones: the field was built with the first, and the water aggregation
  // ran after the second. So each transpose is linearised at the boundary node its
  // own forward pass saw -- the second here, the first restored below.
  //
  // Both of these are the STATE's, not a metric's, so they run once however many
  // metrics ride along.
  compute_boundary_nodes();
  const std::vector<double> density_in_uptake = boundary_density();
  for (size_t m = 0; m < n_seed; ++m) {
    soil_adjoint(lambda_dydt[m], out[m].state, seeds[m], out[m].boundary_node);
    offspring_adjoint(lambda_dydt[m], out[m].state, seeds[m]);
  }
  // Rebuilt on the grid the blocks below index their knot adjoints against, so
  // restoring the first evaluation cannot move a knot.
  compute_environment(false);
  const std::vector<double> density_in_field = boundary_density();

  cohort_block_adjoint(seeds, out);
  last_knot_adjoint = out.front().knot;

  for (size_t m = 0; m < n_seed; ++m) {
    std::vector<node_size_adjoints> sizes(reduction_node_count(),
                                          node_size_adjoints{0, 0, 0, 0});
    light_knot_adjoint(out[m].knot, sizes);
    allometry_adjoint(sizes, out[m].state, out[m].boundary_node, m);
  }
  boundary_condition_adjoint(density_in_field, density_in_uptake, out);

  lambda_y.resize(n_seed);
  for (size_t m = 0; m < n_seed; ++m) {
    lambda_y[m] = std::move(out[m].state);
  }
  // Kept so the single-seed entry point can hand this patch's own accumulator
  // back to a caller that reads it between steps.
  boundary_node_adjoint = std::move(out.front().boundary_node);
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
