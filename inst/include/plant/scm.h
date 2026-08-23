// -*-c++-*-
#ifndef PLANT_PLANT_SCM_H_
#define PLANT_PLANT_SCM_H_

#include <plant/node_schedule.h>
#include <plant/gradient_status.h>
#include <odelia/ode_solver.hpp>
#include <plant/patch.h>
#include <plant/scm_utils.h>

#include <odelia/gradient.hpp>

#include <algorithm>
#include <cmath>
#include <limits>
#include <tuple>

using namespace Rcpp;

namespace plant {

// A census differentiated on both of its input sets, one row per metric each.
// Kept as a pair because they come out of one recording: handing them back
// separately is what let two callers take two recordings of one function.
struct census_rows {
  odelia::ode::row_batch state;  // one column per ODE state entry
  odelia::ode::row_batch trait;  // one column per registered trait
};

// SCM: the "Solver for Characteristics Method" driver.
//
// Owns a Patch (the population being integrated), a NodeSchedule (when each
// species' nodes are introduced), and an ODE Solver. It steps the patch
// forward by repeatedly: introducing all nodes due at the current time, then
// integrating the patch ODE system up to the next introduction time.
//
// The patch owns all of the ecology (fitness, offspring, competition, error
// computations); the SCM is the time-stepping/scheduling layer on top of it.
// Most r_* members are thin facades that expose the C++ API to R via RcppR6.
template <typename T, typename E> class SCM {
public:
  // ---- Type aliases ------------------------------------------------------
  typedef T                strategy_type;
  typedef E                environment_type;
  typedef Individual<T, E> individual_type;
  typedef Node<T, E>       node_type;
  typedef Species<T, E>    species_type;
  typedef Patch<T, E>      patch_type;
  typedef Parameters<T, E> parameters_type;

  // ---- Construction ------------------------------------------------------
  SCM(parameters_type p, environment_type e, plant::Control c);

  // ---- Simulation lifecycle ----------------------------------------------

  // Run the whole schedule from t = 0 to completion.
  void run();

  // Advance one step: introduce every node due at the current time, then
  // integrate the patch forward to the next introduction (or over the fixed
  // ode times). Returns the species indices introduced this step.
  std::vector<size_t> run_next();

  // Replay the resident's saved environment for a mutant strategy: swap in
  // mutant parameters, reuse the cached environment/ode times, and run.
  void run_mutant(parameters_type p);

  // Run, keeping the state at each accepted step, and return one record per step.
  //
  // A run pinned to this run's times and sizes DOES reproduce these states, bit for
  // bit over 3381 steps, and pays 15% less because it attempts no step it will
  // reject. This said the opposite -- that a rejected attempt moves patch state
  // which is not ODE state, and a pinned run makes none -- and the states are
  // measurably identical, so whatever a rejected attempt leaves behind is rebuilt
  // from the state before anything reads it.
  std::vector<ode_step_record> store_trajectory();

  // Set before run() to keep the state at every accepted step. The reverse pass
  // needs those states and cannot recover them from a finished run, so a run that
  // did not keep them has to be repeated -- one whole forward integration. Off by
  // default because the store is one double per state entry per step and a forward
  // run has no use for it.
  bool record_trajectory = false;

  // Adaptively refine the node-introduction schedule entirely in C++:
  // repeatedly run, flag nodes whose combined error exceeds schedule_eps,
  // and bisect the interval below each flagged node (upwind scheme), up to
  // schedule_nsteps times. Replaces the R build_schedule loop.
  void refine_schedule();

  // Return patch, schedule and solver to their t = 0 state; clear history.
  void reset();

  // True once every scheduled node introduction has been consumed.
  bool complete() const;

  // Current patch time.
  double time() const;

  // ---- Outputs -----------------------------------------------------------
  // Total (not per-capita) offspring. These delegate to the patch, which owns
  // the fitness/offspring computations.
  std::vector<double> net_reproduction_ratios() const { return patch.net_reproduction_ratios(); }
  std::vector<double> offspring_production() const { return patch.offspring_production(); }

  // Every metric the strategy declares, summed over the species, in the order it
  // declares them. The codomain is the list's length.
  std::vector<double> census() const;

  // What a census of this model reads, at whatever scalar `P` carries. The one
  // place a strategy is asked, and the one place the question is refused for a
  // strategy that does not answer it: a census over a list nothing declared
  // would otherwise fail inside whichever loop reached for it.
  template <class P>
  static std::vector<census_metric<typename P::strategy_type>> metrics_of() {
    static_assert(Censusable<typename P::strategy_type>,
                  "a census of this model is being taken, so its strategy must "
                  "declare census_metrics()");
    return P::strategy_type::census_metrics();
  }

  // One metric summed over every species of `p`. Templated on the patch so the
  // value and its derivative are the same reduction at two scalars.
  template <class P>
  static typename P::value_type census_over(
      const P& p, const census_metric<typename P::strategy_type>& metric) {
    typename P::value_type tot = 0.0;
    for (size_t i = 0; i < p.size(); ++i) {
      tot += p.at_species(i).census(metric);
    }
    return tot;
  }

  // The reverse pass runs on the birth-date coordinate only, and refuses the
  // other one here rather than answering it. On the height coordinate the
  // abscissa is state, so the quadrature weights carry a derivative nothing
  // supplies and the density rate carries a compression term the recorded step
  // does not compute: the sweep is then the transpose of a function the forward
  // model is not evaluating. Nothing about the arithmetic
  // complains, and the two coordinates are different functions rather than two
  // discretisations of one -- one census metric's trait sensitivity changes
  // sign between them -- so the answer would be finite, plausible and wrong.
  void require_birth_date_coordinate(const char* entry) const {
    if (!control.node_density_in_birth_date) {
      util::stop(std::string(entry) +
                 ": the reverse-mode gradient runs on the birth-date "
                 "size-density coordinate only. Set "
                 "control$node_density_in_birth_date = TRUE and re-run.");
    }
  }

  // d(census)/d(ODE state) and d(census)/d(trait) at the current time, one row
  // per metric each, from one recording. The state half is what the reverse pass
  // is seeded with; the trait half is what no sweep produces, because a metric
  // reads the traits itself and the boundary node's own quantities are rebuilt
  // when the state is set. Columns as ode_state writes them, and as
  // census_trait_gradient reports them.
  census_rows census_state_and_trait_rows() const;

  // d(census)/d(trait), one row per metric and one column per trait in each
  // strategy's ad_parameters() order, species-major. Requires an adaptive run to
  // have resolved the schedule this replays.
  //
  // `extra_splits` names recorded steps at which the sweep stops and resumes. The
  // adjoint recursion is linear in the step, so composition over steps is
  // associative and any split must give the same numbers bit for bit; a
  // difference is something carried across a step boundary that is not the
  // adjoint. Splits outside a segment's interior are ignored, so a caller may
  // pass a boundary index without special-casing it.
  //
  // Returns the numbers AND what each of them is. The two travel together
  // because a row of doubles cannot say whether it is an answer: a refusal
  // anywhere in a metric's sweep makes that metric's whole gradient undefined,
  // and an exact zero is more often a slot nothing reached than a sensitivity
  // the model means.
  census_gradient
  // `which_metrics` names the rows to sweep, empty meaning every one. A metric
  // not asked for is not seeded and not swept, so asking for one costs one --
  // which is what a caller differentiating a single census wants and what
  // computing all of them and subsetting the answer does not give.
  census_trait_gradient(const std::vector<size_t>& extra_splits = {},
                        const std::vector<std::string>& which_metrics = {});

  // The four references the ladder checks census_trait_gradient against.

  // The same quantity differenced, by moving the prepared strategy exactly where
  // the recording seeds it. It referees the trait half above while sharing none of
  // it: that one records the census and sweeps a tape, this one evaluates the
  // census twice. A difference that rebuilt from Parameters would re-run
  // preparation and carry the birth-size channel the differentiated path imposes
  // to zero, so this one perturbs in place.
  std::vector<std::vector<double>> census_trait_difference(double rel);

  // One exact directional derivative of the census, by a forward tangent of the
  // same trajectory stepped at the sizes the run recorded. `direction` carries
  // one weight per trait, species-major in each strategy's ad_parameters()
  // order; a coordinate direction gives one Jacobian column and a mixed one a
  // contraction. Returns one tangent per metric, and writes the metrics the
  // replay itself reached: a reference whose value disagrees with the model is a
  // reference to a different function, and the gap is this check's own floor.
  std::vector<double> census_trait_tangent(const std::vector<double>& direction,
                                           std::vector<double>& value);

  // One exact directional derivative of the census with respect to the first
  // recorded state, by a forward tangent stepped at the sizes the run recorded.
  // `direction` carries one weight per component of that state.
  //
  // No trait, no derived quantity and no census direct term is on this path, so
  // it isolates how the trajectory carries a perturbation to a state a cohort
  // starts at, and census_initial_state_replay is what refereees it.
  // `segment` picks where the seeding happens: 0 is the first recorded state,
  // which on this coordinate holds the environment and no cohort, and `j` is the
  // state the run reached just after the jth introduction. A cohort's own birth
  // height is seedable only from `j >= 1`, because that is where it first exists.
  std::vector<double>
  census_initial_state_tangent(const std::vector<double>& direction,
                               std::vector<double>& value,
                               size_t segment = 0);

  // The census a plain-double replay of the recorded steps reaches from
  // `state0`. Differencing it moves the state the tangent above seeds, through
  // the same steps and the same introductions, so the two differentiate one
  // function and a disagreement is the propagation's own.
  std::vector<double>
  census_initial_state_replay(const std::vector<double>& state0,
                              size_t segment = 0);

  // The state a segment's first step ran from, which is what the two calls above
  // index their arguments against.
  std::vector<double> segment_base_state(size_t segment);

  // The replay both entry points above run. `seed` fills the scalar state the
  // replay starts from, given the recorded one.
  //
  // Storing a trajectory runs the model, so the seeding is handed in rather than
  // applied by the caller: a caller reading the recorded state for itself would
  // store twice and run twice.
  template <class Scalar, class Seed>
  std::vector<Scalar> replay_initial_state(size_t from_segment, Seed seed);

  // The Control entries that move the trajectory or move which states answer,
  // and so move the gradient, in the order stand_gradient() compares them. The
  // curvature floor is here for the second reason rather than the first: it
  // changes no forward number and still decides which rows exist.
  std::vector<double> gradient_control() const {
    return {control.GSS_tol_abs, control.ci_abs_tol, control.node_gradient_eps,
            control.schedule_eps, control.gradient_curvature_floor};
  }

  // ---- R interface -------------------------------------------------------

  // Run / parameters / state access
  parameters_type r_parameters() const { return parameters; }
  const patch_type &r_patch() const { return patch; }

  // The classification tally the FORWARD run built, one row per species and one
  // column per operating-point kind. Read off the live system rather than
  // r_patch(), which is a snapshot the run copies out and whose counters are
  // whatever they were when it was taken.
  // One environment serves every species, so its tally has no species of its
  // own; it lands on the first row rather than being dropped or duplicated.
  static void add_environment_clamps(std::vector<std::vector<size_t>>& ret,
                                     const std::vector<size_t>& env) {
    if (ret.empty()) {
      ret.push_back(env);
      return;
    }
    for (size_t s = 0; s < env.size() && s < ret[0].size(); ++s) {
      ret[0][s] += env[s];
    }
  }

  std::vector<std::vector<size_t>> operating_point_counts() {
    const patch_type& live = solver.get_system_ref();
    std::vector<std::vector<size_t>> ret;
    ret.reserve(live.size());
    for (size_t i = 0; i < live.size(); ++i) {
      ret.push_back(live.at_species(i).strategy_ptr()->operating_point_counts);
    }
    return ret;
  }
  // The same for the clamp sites, which are counted for the same reason: a
  // clamped row and a true zero are the same number.
  //
  // The sites live in two objects and the list is one, so the environment's
  // tally is ADDED to the species' rather than reported beside it -- a caller
  // asking how often a site fired is asking about the site, not about which
  // object happened to reach it. The environment is one, so its counts land on
  // the first row.
  std::vector<std::vector<size_t>> clamp_counts() {
    const patch_type& live = solver.get_system_ref();
    std::vector<std::vector<size_t>> ret;
    ret.reserve(live.size());
    for (size_t i = 0; i < live.size(); ++i) {
      const auto s = live.at_species(i).strategy_ptr();
      std::vector<size_t> row = s->clamps.forward;
      // The leaf's own sites keep ONE tally across both paths, because the leaf
      // solves in double on both and every rebound copy shares its storage. So the
      // forward share is the total less what the sweep was measured to take.
      const std::vector<std::size_t> leaf_total = s->leaf.clamp_counts();
      const std::vector<size_t>& swept = *s->clamps.differentiated;
      for (std::size_t k = 0; k < leaf_total.size(); ++k) {
        const std::size_t at = CLAMP_LEAF_FIRST + k;
        if (at >= CLAMP_SITE_COUNT) { break; }
        row[at] += leaf_total[k] > swept[at] ? leaf_total[k] - swept[at] : 0;
      }
      ret.push_back(row);
    }
    add_environment_clamps(ret, live.environment_clamps().forward);
    return ret;
  }
  // And the same sites counted where the sweep runs, which is the only path a
  // clamp severs a row on. The forward tally cannot stand in for this one: the
  // sweep visits the recorded steps rather than every solve.
  std::vector<std::vector<size_t>> clamp_counts_differentiated() {
    const patch_type& live = solver.get_system_ref();
    std::vector<std::vector<size_t>> ret;
    ret.reserve(live.size());
    for (size_t i = 0; i < live.size(); ++i) {
      ret.push_back(*live.at_species(i).strategy_ptr()->clamps.differentiated);
    }
    add_environment_clamps(ret, *live.environment_clamps().differentiated);
    return ret;
  }
  // The smallest profit curvature the differentiated path met, one per species.
  // A guard that held and a guard nothing reached report the same green, so the
  // distance to the floor is carried out beside the refusals it did not raise.
  std::vector<double> curvature_margins() {
    const patch_type& live = solver.get_system_ref();
    std::vector<double> ret;
    ret.reserve(live.size());
    for (size_t i = 0; i < live.size(); ++i) {
      ret.push_back(*live.at_species(i).strategy_ptr()->curvature_margin);
    }
    return ret;
  }

  void clear_operating_point_counts() {
    const patch_type& live = solver.get_system_ref();
    for (size_t i = 0; i < live.size(); ++i) {
      std::vector<size_t>& c =
        live.at_species(i).strategy_ptr()->operating_point_counts;
      c.assign(c.size(), 0);
      live.at_species(i).strategy_ptr()->clamps.clear();
      live.at_species(i).strategy_ptr()->leaf.clear_clamp_counts();
      *live.at_species(i).strategy_ptr()->curvature_margin = -1.0;
    }
    live.environment_clamps().clear();
  }
  const std::vector<patch_type> &r_history() const { return history; }

  // How many times the inflow boundary was evaluated over one sweep. It is
  // evaluated once per rate evaluation, and a step is recorded once and swept
  // per metric, so the count is six per step and does not scale with the metrics
  // asked for. A row that acts once per stage is multiplied by that count, so
  // the count is part of the row and belongs beside its value.
  size_t boundary_condition_evaluations() { return solver.recorded_rates(); }
  void clear_boundary_condition_evaluations() { solver.clear_recorded_rates(); }
  Rcpp::List r_get_state() const { return patch.r_get_state(); };

  // Fitness / reproduction
  double r_net_reproduction_ratio_for_species(util::index species_index) const;
  std::vector<std::vector<double>> r_net_reproduction_ratio_errors() const;

  // Schedule-refinement error signals.
  // Per-node refinement error: element-wise max of the competition error
  // (sampled during the run) and the reproduction error (computed at the end).
  // This is the signal that drives refine_schedule().
  std::vector<std::vector<double>> refinement_error_by_node() const;
  std::vector<double>
  r_compute_competition_effect_error_by_node_for_species_i(util::index species_index) const;

  // Node schedule access
  NodeSchedule r_node_schedule() const { return node_schedule; }
  void r_set_node_schedule(NodeSchedule x);
  void r_set_node_schedule_times(std::vector<std::vector<double>> x);

  // ODE times: the step times the solver actually used on the last run.
  // (Whether to *pin* integration to a fixed set of times is controlled on the
  // NodeSchedule via its own use_ode_times flag.)
  std::vector<double> r_ode_times() const;

  // The size of the step that reached each of r_ode_times(), NaN first. Pin
  // these on the NodeSchedule alongside the times to replay a run faithfully.
  std::vector<double> r_ode_step_sizes() const;

  // The trajectory as a list of records, each a time, the step size that reached it,
  // and the state there.
  Rcpp::List r_store_trajectory();

  // ---- Public state ------------------------------------------------------
  // The two toggles are exposed to R directly (access: field), so they need
  // no getter/setter wrappers.
  // How many backward ranges the last census_trait_gradient swept, summed over
  // metrics. One per segment with no splits requested, and one more per split
  // that fell inside a segment -- which is what says a requested split cut.
  size_t adjoint_segments = 0;

  // The adjoint the last census_trait_gradient's walk ended holding: one row per
  // metric swept, one column per entry of the first recorded state. It is
  // d(census)/d(that state), which census_initial_state_tangent computes forwards
  // from the same state and over the same steps.
  std::vector<std::vector<double>> adjoint_at_first_state;

  // The time of each recorded step, which the walks over a widened state index
  // their states against.
  std::vector<double> recorded_times(
      const std::vector<ode_step_record>& trajectory) const {
    std::vector<double> ret;
    ret.reserve(trajectory.size());
    for (const ode_step_record& record : trajectory) {
      ret.push_back(record.time);
    }
    return ret;
  }

  // Where the run widened the state, and by what. Filled as the run takes each
  // widening and cleared with the rest of the run's state.
  std::vector<odelia::ode::recorded_widening<typename patch_type::widening>>
      widenings;

  bool collect;                    // record a patch snapshot after each step
  bool collect_refinement_errors;  // accumulate competition errors during run
  std::vector<patch_type> history; // per-step patch snapshots when collect

private:
  // Upwind bisection: insert the midpoint of the interval below each flagged
  // node. Mirrors split_times() in build_schedule.R.
  static std::vector<double> bisect_flagged_intervals(const std::vector<double>& times,
                                                      const std::vector<bool>& split);

  // Uniform grid for fixed-step forward-Euler integration (control.fixed_time_step).
  static std::vector<double> uniform_euler_times(double t0, double t1, double dt);

  // Shared implementation of run_next(). The solver owns the patch system
  // (odelia::ode::Solver), so the live state lives in solver.get_system_ref();
  // sync_patch controls whether the `patch` member snapshot is refreshed from
  // it on return (skipped inside the run() loop to avoid per-step copies).
  std::vector<size_t> run_next_impl(bool sync_patch);

  parameters_type parameters;
  Control control;
  patch_type patch;
  NodeSchedule node_schedule;
  odelia::ode::Solver<patch_type> solver;
};

// ---- Construction --------------------------------------------------------

template <typename T, typename E>
SCM<T, E>::SCM(parameters_type p, environment_type e, Control c)
    : parameters(p), control(c), patch(parameters, e, c),
      node_schedule(make_node_schedule(parameters)),
      solver(patch, make_ode_control(c)) {

  parameters.validate();

  collect = false;
  collect_refinement_errors = false;
  solver.set_collect(false);

  if (!util::identical(parameters.patch_area, 1.0)) {
    util::warning("We recommened keeping patch_area = 1 for the SCM, as need to check units for all other sizes");
  }

}

// Build a uniform grid {t0, t0 + dt, ..., t1} with spacing dt, starting exactly
// at t0 and ending exactly at t1 (the final interval may be shorter than dt).
// Used to drive forward-Euler integration between schedule events.
template <typename T, typename E>
std::vector<double> SCM<T, E>::uniform_euler_times(double t0, double t1,
                                                   double dt) {
  std::vector<double> times;
  times.push_back(t0);
  if (t1 <= t0) {
    return times;
  }
  // Number of (mostly dt-sized) intervals; the small tolerance avoids spawning
  // a spurious tiny final interval when (t1 - t0) is an FP-near multiple of dt.
  const size_t n =
      static_cast<size_t>(std::ceil((t1 - t0) / dt - 1e-10));
  for (size_t i = 1; i < n; ++i) {
    times.push_back(t0 + static_cast<double>(i) * dt);
  }
  times.push_back(t1); // exact endpoint
  return times;
}

// ---- Simulation lifecycle ------------------------------------------------

template <typename T, typename E> void SCM<T, E>::run() {
  // Before reset(), which records the initial state and then copies the patch
  // into the solver: set after it, that first record is missed and the store is
  // one short of the step sizes it is read beside.
  patch.recording = record_trajectory;
  // Set before reset(), which records the state the run starts from.
  solver.set_keep_states(record_trajectory);
  // The choices this run's rate evaluations make are the same recording as its
  // states, so they start over together. Cleared HERE and not in the patch's own
  // reset, which a rebound patch runs on the shared strategies mid-sweep.
  patch.clear_solved_choices();
  reset();
  // The solver owns the live patch system; operate on it directly during the
  // run and avoid per-step copies into the `patch` member.
  if (collect) {
    history.push_back(solver.get_system_ref());
  }

  while (!complete()) {
    std::vector<size_t> added = run_next_impl(false);
    if (collect_refinement_errors) {
      solver.get_system_ref().collect_competition_errors(added);
    }
    if (collect) {
      history.push_back(solver.get_system_ref());
    }
  }

  // Expose the final state through the `patch` accessor after the loop.
  patch = solver.get_system_ref();
}

template <typename T, typename E> std::vector<size_t> SCM<T, E>::run_next() {
  return run_next_impl(true);
}

template <typename T, typename E>
std::vector<size_t> SCM<T, E>::run_next_impl(bool sync_patch) {
  std::vector<size_t> ret;
  const double t0 = time();
  // The live patch system is owned by the solver; mutate it in place.
  auto &sys = solver.get_system_ref();

  NodeSchedule::Event e = node_schedule.next_event();

  // Resume support: if the next scheduled introduction is in the future,
  // integrate the gap up to it without introducing any node. This happens on
  // the first step of a run resumed from an exported state -- the patch is
  // already populated (in reset()) and starts at parameters.initial_time, which
  // falls before the first residual schedule entry. It never happens for an
  // empty patch, whose schedule always starts at t0 = 0, so the normal path
  // below is unchanged. The next call will then introduce at e's time.
  if (e.time_introduction() > t0) {
    solver.set_state_from_system();
    if (node_schedule.using_ode_steps()) {
      util::stop("Resuming from an initial state is not supported for "
                 "replaying a recorded run");
    } else if (control.fixed_time_step > 0.0) {
      solver.advance_euler(
          uniform_euler_times(t0, e.time_introduction(), control.fixed_time_step));
    } else {
      solver.advance_adaptive({solver.time(), e.time_introduction()});
    }
    if (sync_patch) {
      patch = sys;
    }
    return ret; // empty: nothing introduced this step
  }

  // Consume every event scheduled at the current time t0: each contributes a
  // species to introduce. Stop once the next event ends later than t0 (i.e. it
  // belongs to a later introduction) or the schedule is exhausted.
  while (true) {
    if (!util::identical(t0, e.time_introduction())) {
      util::stop("Start time not what was expected");
    }
    ret.push_back(e.species_index);
    node_schedule.pop();
    if (e.time_end() > t0 || complete()) {
      break;
    } else {
      e = node_schedule.next_event();
    }
  }

  sys.introduce_nodes(ret, e.time_introduction());
  // Declared where it happens. Read back off the state's width afterwards it
  // cannot be wrong, because the width is what the reading is derived from.
  if (record_trajectory) {
    widenings.push_back({solver.recorded_steps() - 1, ret});
  }
  solver.set_state_from_system();

  // Three integration modes:
  //  - pinned ode times (resident replay for a mutant): step exactly to the
  //    cached times via the full RKCK stepper, by their recorded step sizes
  //    when the schedule carries them;
  //  - fixed-step forward Euler (control.fixed_time_step > 0): walk a uniform
  //    sub-grid between this event and the next introduction;
  //  - otherwise: adaptive, error-controlled RKCK to the next event time.
  if (node_schedule.using_ode_steps()) {
    if (control.fixed_time_step > 0.0) {
      // The replay relies on the RK sub-step environment cache, which forward
      // Euler does not populate. Refuse rather than mis-integrate.
      util::stop("fixed_time_step (forward Euler) is not supported for a pinned "
                 "ODE schedule");
    }
    // Each recorded step carries the size it took and the time it reached, so a
    // replay lands where the run landed rather than a rounding short of it.
    // Stepping to the times instead would take different steps, because a size
    // differenced back out of two recorded times is not the size that was taken.
    //
    // A size is NaN where the schedule is a grid rather than a recording, and
    // the step is taken TO that time instead. The schedule holds no step at an
    // interval's own end -- a run stops there by clamping its last step to the
    // boundary, and the step below does the same arithmetic.
    // An interval with no step inside it is one the schedule crosses in a single
    // step, so there is nothing to replay before the step to its end. Deciding
    // that per interval rather than per schedule is what silently integrated
    // those intervals adaptively.
    if (!e.steps.empty()) {
      solver.advance_recorded(e.steps);
    }
    solver.advance_fixed({solver.time(), e.times.back()});
  } else if (control.fixed_time_step > 0.0) {
    solver.advance_euler(
        uniform_euler_times(t0, e.time_end(), control.fixed_time_step));
  } else {
    solver.advance_adaptive({solver.time(), e.time_end()});
  }

  if (sync_patch) {
    patch = sys;
  }

  return ret;
}

// An invader integrates against a field it does not move, so the run it needs is
// one that replays a resident's field rather than rebuilding it. The recorder
// that would supply it never ran: it was reached through hooks the solver stopped
// calling, so the field it filled stayed empty and every call arrived here and
// stopped. It is deleted rather than left standing, and what replaces it is
// odelia's ReplaysField, against which this is to be written.
template <typename T, typename E>
void SCM<T, E>::run_mutant(parameters_type /* p */) {
  util::stop("run_mutant needs a recorded resident field to integrate against, "
             "and nothing records one");
}

template <typename T, typename E>
std::vector<ode_step_record> SCM<T, E>::store_trajectory() {
  // What the sweep does with the store: every walk over the widenings takes
  // the patch through these, and each was reached by name until it was said here.
  static_assert(odelia::ode::WidensState<patch_type>,
                "Patch must satisfy WidensState or the segment walk cannot "
                "put it on a recorded step or transpose an insertion");
  // Run only if this run kept no states to read. Reading does not consume them,
  // and a walk puts the patch back on the last recorded step when it is done, so
  // a second consumer reads the same record rather than repeating the run.
  if (!solver.keeps_states()) {
    record_trajectory = true;
    run();
  }

  // Read in place, because the state, the time it was reached at and the size that
  // reached it are now one record on the solver. They used to be two stores, and
  // pairing them was only honest while the states were emptied as they were read --
  // which is why every consumer after the first repeated the whole run. One record
  // cannot be mispaired, so nothing is emptied and nothing is repeated.
  const std::vector<double> times = solver.times();
  const std::vector<double> sizes = solver.step_sizes();
  util::check_length(sizes.size(), times.size());
  std::vector<ode_step_record> ret;
  ret.reserve(times.size());
  for (size_t i = 0; i < times.size(); ++i) {
    ret.push_back({times[i], sizes[i], solver.recorded_state(i)});
  }
  return ret;
}

// Upwind bisection of flagged intervals. For each flagged node j (j >= 1; the
// first and last nodes are never flagged), insert the midpoint of the interval
// (t[j-1], t[j]). Equivalent to sort(c(times, times[i] - dt[i-1]/2)) in R.
template <typename T, typename E>
std::vector<double> SCM<T, E>::bisect_flagged_intervals(const std::vector<double>& times,
                                                        const std::vector<bool>& split) {
  std::vector<double> ret = times;
  for (size_t j = 1; j < split.size(); ++j) {
    if (split[j]) {
      ret.push_back(0.5 * (times[j] + times[j - 1]));
    }
  }
  std::sort(ret.begin(), ret.end());
  return ret;
}

template <typename T, typename E>
void SCM<T, E>::refine_schedule() {
  collect_refinement_errors = true;
  const double eps = control.schedule_eps;

  for (size_t step = 0; step < control.schedule_nsteps; ++step) {
    run(); // resets, then runs with collect_refinement_errors set

    std::vector<std::vector<double>> node_error = refinement_error_by_node();

    // Flag nodes whose refinement error exceeds the threshold.
    std::vector<std::vector<bool>> split(node_error.size());
    bool any = false;
    for (size_t i = 0; i < node_error.size(); ++i) {
      split[i].assign(node_error[i].size(), false);
      for (size_t j = 0; j < node_error[i].size(); ++j) {
        if (node_error[i][j] > eps) {
          split[i][j] = true;
          any = true;
        }
      }
    }
    if (!any) {
      break; // converged: no interval needs refining
    }

    // Bisect flagged intervals and install the denser schedule.
    std::vector<std::vector<double>> times = node_schedule.get_times();
    for (size_t i = 0; i < times.size(); ++i) {
      times[i] = bisect_flagged_intervals(times[i], split[i]);
    }
    node_schedule.set_times(times);
  }

  // Leave Parameters self-describing: record the refined schedule and the
  // ode times from the final run (mirrors build_schedule.R).
  parameters.node_schedule_times = node_schedule.get_times();
  parameters.ode_times = r_ode_times();
  parameters.ode_step_sizes = solver.step_sizes();
}

// NOTE: solver.reset() sets the solver's internal time to zero. There is
// currently no other way to set that time; it might be cleaner to add an
// odelia::ode::Solver::set_time and call set_time(0) explicitly here.
template <typename T, typename E> void SCM<T, E>::reset() {
  patch.reset();
  node_schedule.reset();
  // Seed the solver's owned system from the freshly reset patch, then reset
  // the solver's time/step state and sync the snapshot back.
  solver.get_system_ref() = patch;
  solver.reset();
  patch = solver.get_system_ref();
  history.clear();
  widenings.clear();
  adjoint_segments = 0;
  adjoint_at_first_state.clear();
}


template <typename T, typename E> bool SCM<T, E>::complete() const {
  return node_schedule.remaining() == 0;
}

template <typename T, typename E> double SCM<T, E>::time() const {
  return solver.time();
}

// ---- R interface ---------------------------------------------------------
//
// The fitness/offspring and per-node error computations live on the patch
// (patch.h); the SCM methods below are thin facades that preserve the R API.
//
// Several of these are diagnostic/inspection hooks rather than part of the
// production run path: outside the C++ refinement loop they are only called
// from the test suite and the node_spacing vignette (noted per method below).

// Per-species fitness: the net reproduction ratio (expected offspring per seed)
// for one species. A genuine biological quantity, not just a diagnostic.
template <typename T, typename E>
double SCM<T, E>::r_net_reproduction_ratio_for_species(
    util::index species_index) const {
  const size_t idx = species_index.check_bounds(patch.size());
  auto scalars = std::vector<double>(patch.at_species(idx).size(), 1.0);
  return patch.net_reproduction_ratio_for_species(idx, scalars);
}

// Diagnostic: per-node discretisation error in the reproduction integral. One
// of the two components of refinement_error_by_node; exposed for inspection
// and validation (tests / vignette).
template <typename T, typename E>
std::vector<std::vector<double>>
SCM<T, E>::r_net_reproduction_ratio_errors() const {
  return patch.net_reproduction_ratio_errors();
}

// The combined per-node refinement error. Used internally by refine_schedule();
// also exposed to R so tests / the vignette can inspect the signal that drives
// schedule refinement.
template <typename T, typename E>
std::vector<std::vector<double>> SCM<T, E>::refinement_error_by_node() const {
  return patch.refinement_error_by_node();
}

// Diagnostic probe: the per-node competition (light) error for one species --
// the per-step sample that collect_competition_errors() accumulates. Exposed
// mainly so tests / the vignette can reconstruct the error signal by hand.
template <typename T, typename E>
std::vector<double>
SCM<T, E>::r_compute_competition_effect_error_by_node_for_species_i(util::index species_index) const {
  // The per-node error is scaled by the total competition effect inside the
  // patch-level call below: it computes compute_competition(0.0) -- which
  // already divides by patch area -- and passes it through as the scaling
  // argument. The live schedule-refinement collector
  // (Patch::collect_competition_errors) reconstructs the signal via this same
  // path, so no extra area scaling is needed here to keep them consistent
  // (resolves the scaling question in #478).
  const size_t idx = species_index.check_bounds(patch.size());
  return patch.r_compute_competition_effect_error_by_node_for_species_i(idx);
}

template <typename T, typename E>
void SCM<T, E>::r_set_node_schedule(NodeSchedule x) {
  if (patch.node_ode_size() > 0) {
    util::stop("Cannot set schedule without resetting first");
  }
  util::check_length(x.get_n_species(), patch.size());
  node_schedule = x;

  // Update here so that extracting Parameters reflects the new schedule,
  // keeping Parameters self-sufficient.
  parameters.node_schedule_times = node_schedule.get_times();
}

template <typename T, typename E>
void SCM<T, E>::r_set_node_schedule_times(
    std::vector<std::vector<double>> x) {
  if (patch.node_ode_size() > 0) {
    util::stop("Cannot set schedule without resetting first");
  }
  node_schedule.set_times(x);
  parameters.node_schedule_times = x;
}

template <typename T, typename E>
std::vector<double> SCM<T, E>::r_ode_times() const {
  return solver.times();
}

template <typename T, typename E>
std::vector<double> SCM<T, E>::r_ode_step_sizes() const {
  return solver.step_sizes();
}

template <typename T, typename E>
Rcpp::List SCM<T, E>::r_store_trajectory() {
  const std::vector<ode_step_record> steps = store_trajectory();
  Rcpp::List ret(steps.size());
  for (size_t i = 0; i < steps.size(); ++i) {
    ret[i] = Rcpp::List::create(Rcpp::_["time"] = steps[i].time,
                                Rcpp::_["step_size"] = steps[i].step_size,
                                Rcpp::_["state"] = steps[i].state);
  }
  return ret;
}

template <typename T, typename E>
std::vector<double> SCM<T, E>::census() const {
  const std::vector<census_metric<T>> metrics = metrics_of<patch_type>();
  std::vector<double> ret;
  ret.reserve(metrics.size());
  for (const census_metric<T>& metric : metrics) {
    ret.push_back(odelia::util::to_passive(census_over(patch, metric)));
  }
  return ret;
}

// The census differentiated with respect to the state AND the traits, from ONE
// recording of one metric algebra. The two used to be two functions over the same
// arithmetic, which put the seam between the halves in two places; here it is the
// solver's own, written once for every transpose in the tree.
//
// The order the halves are written in is the reason they are one recording rather
// than two calls: the traits are seated first and the state loaded after, so a
// quantity the state determines is derived at the traits the recording registered.
// Loading the state first derives it at the values they had before.
//
// set_recorded_state rebuilds the environment and the boundary node from the state
// it is given. Both are on the census's path -- the boundary node is the
// reduction's lower grid point and is not ODE state -- so the recording must carry
// that rebuild. Loading the state without it leaves the boundary node at the values
// it was copied with, and its whole contribution to the seed is then exactly zero
// with nothing thrown. Loading it with set_ode_state alone leaves the condition at
// its first evaluation, which is not the one census() reads.
//
// One patch, one tape, one recording, and a seed per metric. The recording does
// not depend on which metric is being asked for -- it writes every metric into y
// and only the seed picks one out -- so a recording per metric was a recording
// repeated.
//
// What made that repetition look necessary is real and is worth stating, because
// it is the trap next door. Clearing a tape returns its derivative-slot counter
// to zero, so an active value built outside a sweep loop and read inside it
// refers, after the first clear, to a slot that now belongs to something else.
// Measured when that was live: the second and third metrics' seeds were wrong by
// three orders and their heartwood columns read exactly zero -- the first metric
// correct and lending its credibility to the rest. The answer is not a patch per
// metric, it is to clear once and record once, which is what sweeping a batch
// does: the clear happens before the recording, and between sweeps only the
// derivative slots are returned to zero.
template <typename T, typename E>
census_rows SCM<T, E>::census_state_and_trait_rows() const {
  require_birth_date_coordinate("census_state_and_trait_rows");
  using scalar = odelia::ode::active_scalar<double>;

  std::vector<double> state(patch.ode_size());
  patch.ode_state(state.begin());

  const size_t n_metric = metrics_of<patch_type>().size();
  census_rows ret;
  ret.trait.assign(n_metric, patch.trait_adjoint_size());

  typename scalar::tape_type tape(false);
  auto reduce = [&](auto& active, typename std::vector<scalar>::const_iterator x,
                    std::vector<scalar>& y) -> void {
    // The traits are already written from the other half of the recorded inputs.
    active.set_recorded_state(x, time());
    const auto metrics = metrics_of<std::decay_t<decltype(active)>>();
    for (size_t m = 0; m < metrics.size(); ++m) {
      y[m] = census_over(active, metrics[m]);
    }
  };
  odelia::ode::state_and_parameter_adjoints(
      tape, patch, state, odelia::ode::row_batch::all_rows(n_metric), reduce,
      ret.state, ret.trait);
  return ret;
}

// The census twice per trait, at the state held, with the strategy moved in place.
// See the declaration for why it perturbs rather than rebuilds.
template <typename T, typename E>
std::vector<std::vector<double>>
SCM<T, E>::census_trait_difference(double rel) {
  require_birth_date_coordinate("census_trait_difference");
  const std::vector<census_metric<T>> metrics = metrics_of<patch_type>();
  const size_t n_metric = metrics.size();
  const size_t n_state = patch.ode_size();

  std::vector<double> state(n_state);
  patch.ode_state(state.begin());
  const double time_ = time();

  // The patch answers for this order; walking the species here would be free to
  // walk it differently.
  const std::vector<typename T::value_type*> pars = patch.ad_parameters();

  // The state is re-set on every evaluation, which is what makes the moved trait
  // reach the quantities a state determines -- the boundary node among them.
  auto census_at = [&](std::vector<double>& out) -> void {
    patch.set_recorded_state(state.begin(), time_);
    out.clear();
    for (const census_metric<T>& metric : metrics) {
      out.push_back(odelia::util::to_passive(census_over(patch, metric)));
    }
  };

  std::vector<std::vector<double>> ret(n_metric,
                                       std::vector<double>(pars.size(), 0.0));
  std::vector<double> up, dn;
  for (size_t c = 0; c < pars.size(); ++c) {
    const double base = odelia::util::to_passive(*pars[c]);
    const double h = std::max(std::abs(base) * rel, rel);
    *pars[c] = base + h;
    census_at(up);
    *pars[c] = base - h;
    census_at(dn);
    *pars[c] = base;
    for (size_t m = 0; m < n_metric; ++m) {
      ret[m][c] = (up[m] - dn[m]) / (2.0 * h);
    }
  }
  // Leave the patch where it was found, so this call is repeatable beside the
  // recording that shares its state.
  census_at(up);
  return ret;
}

// Seed lambda on the states the census reads at T, then run the reverse pass
// back over the recorded steps. The trait adjoints accumulate across every
// cohort and every step, so the accumulator is cleared once per metric and read
// once the sweep is done. It is the solver's system that accumulates: `patch` is
// a snapshot the run copies out, and reading its accumulator gives zeros.
template <typename T, typename E>
census_gradient
SCM<T, E>::census_trait_gradient(const std::vector<size_t>& extra_splits,
                                 const std::vector<std::string>& which_metrics) {
  require_birth_date_coordinate("census_trait_gradient");
  // Which rows to sweep, resolved against the strategy's own list before
  // anything runs, so the shape of the answer is known on the refusal path too.
  // Named rather than positional: a caller indexing by position gets a different
  // metric's gradient when the list changes, and nothing says so.
  const std::vector<census_metric<T>> metrics = metrics_of<patch_type>();
  std::vector<size_t> rows;
  if (which_metrics.empty()) {
    rows.resize(metrics.size());
    for (size_t m = 0; m < rows.size(); ++m) {
      rows[m] = m;
    }
  } else {
    for (const std::string& want : which_metrics) {
      size_t at = metrics.size();
      for (size_t m = 0; m < metrics.size(); ++m) {
        if (want == metrics[m].name) {
          at = m;
          break;
        }
      }
      if (at == metrics.size()) {
        std::string known;
        for (const census_metric<T>& metric : metrics) {
          known += known.empty() ? "" : ", ";
          known += metric.name;
        }
        util::stop("census_trait_gradient: this model has no census metric `" +
                   want + "`; it has " + known);
      }
      rows.push_back(at);
    }
  }
  // The sweep needs the state at every accepted step. store_trajectory() repeats
  // the run to get them unless record_trajectory kept them the first time, and
  // either way it may run, so the seeds below are taken after it.
  const std::vector<ode_step_record> trajectory = store_trajectory();
  std::vector<std::vector<double>> states;
  states.reserve(trajectory.size());
  for (const ode_step_record& record : trajectory) {
    states.push_back(record.state);
  }

  patch_type& live = solver.get_system_ref();

  adjoint_segments = 0;
  adjoint_at_first_state.clear();
  // The missing-row flag latches, and it is shared storage rather than a member
  // of the copy that sets it -- so it is cleared where the call that reads it
  // starts, and a previous call's degeneracy cannot refuse this one.
  for (size_t i = 0; i < live.size(); ++i) {
    const auto s = live.at_species(i).strategy_ptr();
    *s->uptake_rows_unavailable = false;
    s->uptake_rows_reason->clear();
  }
  // The seeds and the direct term are on the gradient path too, and a refusal
  // raised while forming them has no sweep to be caught in -- so it escaped as
  // an error while every refusal from the sweep itself came back as a status.
  // Caught here, where the shape of the answer is still known from what the
  // caller asked for and the patch's trait width rather than from the seeds
  // themselves. Nothing is restored: the walk that borrows the width runs below.
  odelia::ode::row_batch all_seeds, all_direct;
  try {
    const census_rows both = census_state_and_trait_rows();
    all_seeds = std::move(both.state);
    all_direct = std::move(both.trait);
  } catch (gradient_refusal& e) {
    const size_t width = live.trait_adjoint_size();
    census_gradient ret;
    for (size_t m = 0; m < rows.size(); ++m) {
      ret.gradient.push_back(std::vector<double>(
          width, std::numeric_limits<double>::quiet_NaN()));
      ret.status.push_back(std::vector<gradient_status>(width, e.status));
    }
    return ret;
  }
  const odelia::ode::row_batch seeds = all_seeds.select(rows);
  const odelia::ode::row_batch direct = all_direct.select(rows);
  // Every metric's sweep visits the same trajectory and differs only in its
  // seed, so they are carried TOGETHER: a block is recorded once and swept once
  // per metric, where the loop this replaces recorded it once per metric. The
  // recording is a model evaluation and a sweep is arithmetic, so the second and
  // third metrics were costing what the first did and now cost almost nothing.
  const size_t n_metric = seeds.rows();

  // The accumulator the sweep adds into, owned here for the length of the sweep
  // rather than kept on the patch between calls. Every writer reaches it through
  // the driver, so a row that does not match the seeds is a length mismatch.
  odelia::ode::row_batch trait_adjoint(n_metric, live.trait_adjoint_size());
  odelia::ode::row_batch lambda = seeds;
  // One segment per width, highest first, narrowing across each widening and
  // transposing the map that took it. The solver owns that walk: what is left
  // here is the census the sweep is seeded from and the direct term below.
  //
  // A refusal costs every metric, and the grain is forced rather than chosen: the
  // row that could not be supplied is an intermediate of a recording spanning six
  // stages and every cohort in them, so no seed carries a component to attribute
  // it to. Nothing is unwound here -- the walk restores the width it borrowed at
  // the tail either way.
  bool refused = false;
  gradient_status refusal;
  try {
    adjoint_segments =
        n_metric * odelia::ode::solve_adjoint_over_widenings(
                       solver, states, widenings, lambda, trait_adjoint,
                       extra_splits);
    adjoint_at_first_state = lambda.to_rows();
  } catch (gradient_refusal& e) {
    refusal = e.status;
    refused = true;
    adjoint_segments = 0;
    adjoint_at_first_state.clear();
  }
  // A row the water outputs need can go missing without any refusal being
  // thrown: the recording carries the values on and leaves those rows off the
  // tape, which makes them exactly zero, and a zero row and an absent one are
  // the same number. So the leaf records the loss instead of throwing, and it is
  // read here -- the one place in this call where a refusal has somewhere to go.
  //
  // The strategy that measured it is a per-unit copy the sweep discards, so the
  // flag lives in storage the run still owns, shared through the rebind exactly
  // as the clamp counts are.
  if (!refused) {
    for (size_t i = 0; i < live.size(); ++i) {
      const auto s = live.at_species(i).strategy_ptr();
      if (!*s->uptake_rows_unavailable) {
        continue;
      }
      refusal.kind = gradient_status::Kind::refused;
      refusal.reason = *s->uptake_rows_reason;
      refusal.species = static_cast<int>(i + 1);
      refused = true;
      break;
    }
  }

  census_gradient ret;
  ret.gradient.reserve(n_metric);
  ret.status.reserve(n_metric);
  const std::vector<gradient_status::Kind> zero_classes =
      live.trait_adjoint_zero_classes();
  for (size_t m = 0; m < n_metric; ++m) {
    if (refused) {
      // A sum has no defined value with an undefined term, so the whole metric
      // goes -- not the cohort's column and not the parameter's entry. The
      // numbers are not-a-number rather than absent so that a caller indexing by
      // position still finds the shape it expects.
      ret.gradient.push_back(std::vector<double>(
          direct.width(), std::numeric_limits<double>::quiet_NaN()));
      ret.status.push_back(std::vector<gradient_status>(direct.width(), refusal));
      continue;
    }
    std::vector<double> row(trait_adjoint[m].begin(), trait_adjoint[m].end());
    util::check_length(row.size(), direct.width());
    for (size_t p = 0; p < row.size(); ++p) {
      row[p] += direct[m][p];
    }
    // An exact zero is read against what the column declares one would mean,
    // and a column that declares nothing says so rather than passing the number
    // off as an answer. Only exact zeros are classified: a parameter that is
    // live is not zero, so no declaration is consulted for it.
    std::vector<gradient_status> row_status(row.size());
    for (size_t p = 0; p < row.size(); ++p) {
      if (row[p] == 0.0) {
        row_status[p].kind = zero_classes[p];
      }
    }
    ret.status.push_back(row_status);
    ret.gradient.push_back(row);
  }

  // Leave the system at the width the run left it, so this call is repeatable.
  const std::vector<double> times = recorded_times(trajectory);
  odelia::ode::be_at_step(live,
                          odelia::ode::insertions_of(widenings, times), states,
                          times, states.size() - 1);
  return ret;
}


// The reference the trajectory sweep is checked against. It is a tangent of the
// same forward source: exact, with no step size of its own and no truncation,
// and it traverses both reductions and the introduction boundary while none of
// the transposes under test are on its path.
//
// The recorded step sizes are replayed rather than the times, and rather than a
// controller of its own. A tangent run left to choose its own steps
// differentiates the controller, which the model does not contain -- and a size
// differenced back out of two recorded times is not the size that was taken,
// since fl(fl(t + h) - t) != h.
template <typename T, typename E>
std::vector<double>
SCM<T, E>::census_trait_tangent(const std::vector<double>& direction,
                                std::vector<double>& value) {
  require_birth_date_coordinate("census_trait_tangent");

  const std::vector<ode_step_record> trajectory = store_trajectory();
  std::vector<std::vector<double>> states;
  states.reserve(trajectory.size());
  for (const ode_step_record& record : trajectory) {
    states.push_back(record.state);
  }
  const std::vector<double> times = recorded_times(trajectory);
  const auto insertions = odelia::ode::insertions_of(widenings, times);
  patch_type& live = solver.get_system_ref();
  odelia::ode::be_at_step(live, insertions, states, times, 0);
  auto active = live.template rebind_from<tangent>();

  // Seeded before the state is set: the quantities a state determines read the
  // parameters, and would otherwise be derived at the unseeded values.
  size_t at = 0;
  for (tangent* p : active.ad_parameters()) {
    if (at >= direction.size()) {
      util::stop("census_trait_tangent: one weight per trait, species-major");
    }
    seed_direction(*p, direction[at++]);
  }
  util::check_length(direction.size(), at);
  std::vector<tangent> x0(states[0].size());
  for (size_t i = 0; i < x0.size(); ++i) {
    x0[i] = states[0][i];
  }
  active.set_ode_state(x0.begin(), trajectory[0].time);

  odelia::ode::Solver<decltype(active)> forward(active, make_ode_control(control));
  forward.set_collect(false);
  forward.set_state_from_system();

  std::vector<odelia::ode::recorded_step> steps;
  steps.reserve(trajectory.size());
  for (const ode_step_record& record : trajectory) {
    steps.push_back({record.time, record.step_size});
  }
  odelia::ode::advance_over_widenings(forward, insertions, steps, 0, 0);

  // Leave the double system where the run left it, so this call is repeatable
  // beside the sweep that shares its trajectory.
  odelia::ode::be_at_step(solver.get_system_ref(), insertions, states, times,
                          states.size() - 1);

  const auto& reached = forward.get_system_ref();
  const auto metrics = metrics_of<std::decay_t<decltype(reached)>>();
  std::vector<double> ret;
  ret.reserve(metrics.size());
  value.clear();
  value.reserve(metrics.size());
  // Reduced ONCE per metric, and its value read off the same tangent as its
  // derivative. Two reductions would be two evaluations of one function, which is
  // the shape this whole check exists to catch elsewhere.
  for (const auto& metric : metrics) {
    const tangent reached_metric = census_over(reached, metric);
    ret.push_back(derivative_along(reached_metric));
    value.push_back(odelia::util::to_passive(reached_metric));
  }
  return ret;
}


template <typename T, typename E>
std::vector<double> SCM<T, E>::segment_base_state(size_t segment) {
  const std::vector<ode_step_record> trajectory = store_trajectory();
  std::vector<std::vector<double>> states;
  states.reserve(trajectory.size());
  for (const ode_step_record& record : trajectory) {
    states.push_back(record.state);
  }
  patch_type& live = solver.get_system_ref();
  const std::vector<double> times = recorded_times(trajectory);
  const auto insertions = odelia::ode::insertions_of(widenings, times);

  std::vector<double> base;
  size_t start = 0;
  odelia::ode::state_at_segment(live, insertions, states, times, segment, base,
                                start);
  odelia::ode::be_at_step(live, insertions, states, times, states.size() - 1);
  return base;
}

template <typename T, typename E>
template <class Scalar, class Seed>
std::vector<Scalar> SCM<T, E>::replay_initial_state(size_t from_segment,
                                                    Seed seed) {
  const std::vector<ode_step_record> trajectory = store_trajectory();
  std::vector<std::vector<double>> states;
  states.reserve(trajectory.size());
  for (const ode_step_record& record : trajectory) {
    states.push_back(record.state);
  }
  patch_type& live = solver.get_system_ref();
  const std::vector<double> times = recorded_times(trajectory);
  const auto insertions = odelia::ode::insertions_of(widenings, times);

  std::vector<double> base;
  size_t start = 0;
  const double t0 = odelia::ode::state_at_segment(live, insertions, states, times,
                                                  from_segment, base, start);

  auto active = live.template rebind_from<Scalar>();

  std::vector<Scalar> x0(base.size());
  seed(x0, base);
  active.set_ode_state(x0.begin(), t0);

  odelia::ode::Solver<decltype(active)> forward(active, make_ode_control(control));
  forward.set_collect(false);
  forward.set_state_from_system();

  std::vector<odelia::ode::recorded_step> steps;
  steps.reserve(trajectory.size());
  for (const ode_step_record& record : trajectory) {
    steps.push_back({record.time, record.step_size});
  }
  odelia::ode::advance_over_widenings(forward, insertions, steps, from_segment,
                                     start);

  // Leave the double system where the run left it, so this call is repeatable
  // beside the sweep that shares its trajectory.
  odelia::ode::be_at_step(live, insertions, states, times, states.size() - 1);

  const auto& reached = forward.get_system_ref();
  const auto metrics = metrics_of<std::decay_t<decltype(reached)>>();
  std::vector<Scalar> out;
  out.reserve(metrics.size());
  for (const auto& metric : metrics) {
    out.push_back(census_over(reached, metric));
  }
  return out;
}

template <typename T, typename E>
std::vector<double>
SCM<T, E>::census_initial_state_tangent(const std::vector<double>& direction,
                                        std::vector<double>& value,
                                        size_t segment) {
  require_birth_date_coordinate("census_initial_state_tangent");

  const std::vector<tangent> reached =
    replay_initial_state<tangent>(segment,
      [&](std::vector<tangent>& x0,
          const std::vector<double>& base) -> void {
        util::check_length(direction.size(), base.size());
        for (size_t i = 0; i < x0.size(); ++i) {
          x0[i] = base[i];
          seed_direction(x0[i], direction[i]);
        }
      });

  std::vector<double> ret;
  ret.reserve(reached.size());
  value.clear();
  value.reserve(reached.size());
  for (const tangent& metric : reached) {
    ret.push_back(derivative_along(metric));
    value.push_back(odelia::util::to_passive(metric));
  }
  return ret;
}

template <typename T, typename E>
std::vector<double>
SCM<T, E>::census_initial_state_replay(const std::vector<double>& state0,
                                       size_t segment) {
  require_birth_date_coordinate("census_initial_state_replay");
  return replay_initial_state<double>(segment,
    [&](std::vector<double>& x0, const std::vector<double>& base) -> void {
      util::check_length(state0.size(), base.size());
      x0 = state0;
    });
}


} // namespace plant

#endif
