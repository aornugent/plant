// -*-c++-*-
#ifndef PLANT_PLANT_SCM_H_
#define PLANT_PLANT_SCM_H_

#include <plant/node_schedule.h>
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
  // The states are the run's own, recorded as it goes: a run pinned to this run's
  // times and sizes does not reproduce them, because a rejected step attempt moves
  // patch state that is not part of the ODE state and a pinned run makes no such
  // attempts.
  std::vector<ode_step_record> store_trajectory();

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

  // Each metric of `Metrics`, summed over the species, in tuple order. The
  // codomain is the tuple's size.
  template <class Metrics> std::vector<double> census() const;

  // One metric summed over every species of `p`. Templated on the patch so the
  // value and its derivative are the same reduction at two scalars.
  template <class P, class Psi>
  static typename P::value_type census_over(const P& p, Psi psi) {
    typename P::value_type tot = 0.0;
    for (size_t i = 0; i < p.size(); ++i) {
      tot += p.at_species(i).census(psi);
    }
    return tot;
  }

  // The reverse pass runs on the birth-date coordinate only, and refuses the
  // other one here rather than answering it. On the height coordinate the
  // abscissa is state, so the quadrature weights carry a derivative the
  // reduction transposes omit and the density rate carries a compression term
  // the recorded step does not compute: the sweep is then the transpose of a
  // function the forward model is not evaluating. Nothing about the arithmetic
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

  // d(census)/d(ODE state) at the current time, one row per metric and one
  // column per ODE state entry. This is what the reverse pass is seeded with.
  template <class Metrics>
  std::vector<std::vector<double>> census_state_adjoint() const;

  // d(census)/d(trait) at the state held, which no sweep produces: a metric
  // reads the traits itself, and the boundary node's own quantities are rebuilt
  // when the state is set. One row per metric, columns as census_trait_gradient.
  template <class Metrics>
  std::vector<std::vector<double>> census_trait_direct();

  // The same quantity differenced, by moving the prepared strategy exactly where
  // the recording seeds it. It referees census_trait_direct while sharing none of
  // it: that one records the census and sweeps a tape, this one evaluates the
  // census twice. A difference that rebuilt from Parameters would re-run
  // preparation and carry the birth-size channel the differentiated path imposes
  // to zero, so this one perturbs in place.
  template <class Metrics>
  std::vector<std::vector<double>> census_trait_difference(double rel);

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
  template <class Metrics>
  std::vector<std::vector<double>>
  // `which_metrics` names the rows to sweep, empty meaning every one. A metric
  // not asked for is not seeded and not swept, so asking for one costs one --
  // which is what a caller differentiating a single census wants and what
  // computing all of them and subsetting the answer does not give.
  census_trait_gradient(const std::vector<size_t>& extra_splits = {},
                        const std::vector<size_t>& which_metrics = {});

  // One exact directional derivative of the census, by a forward tangent of the
  // same trajectory stepped at the sizes the run recorded. `direction` carries
  // one weight per trait, species-major in each strategy's ad_parameters()
  // order; a coordinate direction gives one Jacobian column and a mixed one a
  // contraction. Returns one tangent per metric, and writes the metrics the
  // replay itself reached: a reference whose value disagrees with the model is a
  // reference to a different function, and the gap is this check's own floor.
  template <class Metrics>
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
  template <class Metrics>
  std::vector<double>
  census_initial_state_tangent(const std::vector<double>& direction,
                               std::vector<double>& value,
                               size_t segment = 0);

  // The census a plain-double replay of the recorded steps reaches from
  // `state0`. Differencing it moves the state the tangent above seeds, through
  // the same steps and the same introductions, so the two differentiate one
  // function and a disagreement is the propagation's own.
  template <class Metrics>
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
  template <class Metrics, class Scalar, class Seed>
  std::vector<Scalar> replay_initial_state(size_t from_segment, Seed seed);

  // The state and time a segment's first step ran from, with the solver's system
  // left carrying that segment's width. Reached by replaying the introductions
  // before it, since a widened state is not what any record holds.
  double narrow_to_segment(const std::vector<std::vector<double>>& states,
                           const std::vector<ode_step_record>& trajectory,
                           const std::vector<size_t>& boundary,
                           const std::vector<std::vector<size_t>>& introduced,
                           size_t from_segment,
                           std::vector<double>& base,
                           size_t& start);

  // Where the recorded trajectory widens and which species widened it. Returns
  // with the solver's system narrowed to the first segment's width.
  void narrow_over_introductions(
      const std::vector<std::vector<double>>& states,
      const std::vector<ode_step_record>& trajectory,
      std::vector<size_t>& boundary,
      std::vector<std::vector<size_t>>& introduced);

  // Re-introduce every node the sweep's narrowing removed, recovering the state
  // each block's first step ran from. See the definition.
  void widen_over_introductions(
      const std::vector<std::vector<double>>& states,
      const std::vector<ode_step_record>& trajectory,
      const std::vector<size_t>& boundary,
      const std::vector<std::vector<size_t>>& introduced,
      std::vector<std::vector<double>>& sweep_states);

  // The four Control entries that move the trajectory and so move the gradient,
  // in the order stand_gradient() compares them.
  std::vector<double> gradient_control() const {
    return {control.GSS_tol_abs, control.ci_abs_tol, control.node_gradient_eps,
            control.schedule_eps};
  }

  // ---- R interface -------------------------------------------------------

  // Run / parameters / state access
  parameters_type r_parameters() const { return parameters; }
  const patch_type &r_patch() const { return patch; }
  const std::vector<patch_type> &r_history() const { return history; }

  // How many times the inflow boundary's own term entered the trait adjoint over
  // one sweep, read from the system the sweep runs on rather than from the copy
  // r_patch() hands out. A row that acts once per step is multiplied by that
  // count, so the count is part of the row and belongs beside its value.
  std::vector<size_t> boundary_condition_evaluations() {
    const patch_type& live = solver.get_system_ref();
    return {live.boundary_condition_asked, live.boundary_condition_carried};
  }
  void clear_boundary_condition_evaluations() {
    patch_type& live = solver.get_system_ref();
    live.boundary_condition_asked = 0;
    live.boundary_condition_carried = 0;
  }
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

  // Forward-Euler integration (control.fixed_time_step > 0) has no analogue for
  // the RK sub-step environment cache used to replay residents for mutants, so
  // refuse the combination up front rather than produce a wrong fitness.
  if (control.fixed_time_step > 0.0 && control.save_RK45_cache) {
    util::stop("fixed_time_step (forward Euler) is incompatible with "
               "save_RK45_cache / the mutant-fitness replay path");
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
    if (node_schedule.using_ode_times()) {
      util::stop("Resuming from an initial state is not supported for "
                 "ode-time replay / mutant runs");
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

  sys.introduce_new_nodes(ret);
  solver.set_state_from_system();

  // Three integration modes:
  //  - pinned ode times (resident replay for a mutant): step exactly to the
  //    cached times via the full RKCK stepper, by their recorded step sizes
  //    when the schedule carries them;
  //  - fixed-step forward Euler (control.fixed_time_step > 0): walk a uniform
  //    sub-grid between this event and the next introduction;
  //  - otherwise: adaptive, error-controlled RKCK to the next event time.
  if (node_schedule.using_ode_times()) {
    if (control.fixed_time_step > 0.0) {
      // The mutant replay path relies on the RK sub-step environment cache,
      // which forward Euler does not populate. Refuse rather than mis-integrate.
      util::stop("fixed_time_step (forward Euler) is not supported for "
                 "ode-time replay / mutant runs");
    }
    if (e.step_sizes.empty()) {
      solver.advance_fixed(e.times);
    } else {
      // fl(fl(t + h) - t) != h, so a size differenced back out of the recorded
      // times is not the size that was taken; step by the recorded sizes and
      // difference only the last step, which is where the free run itself
      // landed on the interval end rather than accumulating.
      solver.advance_fixed_steps(e.step_sizes);
      solver.advance_fixed({solver.time(), e.times.back()});
    }
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

template <typename T, typename E>
void SCM<T, E>::run_mutant(parameters_type p) {

  // Switch the patch to its cached (resident) environment.
  patch.set_mutant();

  // Destructive: overwrite the resident parameters with the mutant's.
  parameters = p;

  // Swap in the mutant strategies.
  patch.overwrite_strategies(parameters.strategies);

  // Rebuild the schedule for the new parameters, then pin its integration
  // points to the resident's step history so the mutant sees the same
  // environment trajectory.
  node_schedule = make_node_schedule(parameters);
  node_schedule.r_set_ode_times(patch.step_history);
  node_schedule.r_set_use_ode_times(true);
  node_schedule.reset();

  // Re-initialise solver/patch and run.
  reset();
  run();
}

template <typename T, typename E>
std::vector<ode_step_record> SCM<T, E>::store_trajectory() {
  // The stepper calls record_ode_step() only for a System it recognises as
  // Replayable, so this is what makes the store's mechanism present rather than
  // silently absent -- the four members are satisfied by name, not by declaration.
  static_assert(odelia::ode::Replayable<patch_type>,
                "Patch must satisfy Replayable or the stepper records nothing");
  patch.record_steps = true;
  run();
  patch.record_steps = false;

  // The solver holds the size of the step that reached each state, one per record
  // including the initial NaN, so each size lands beside the state it reached.
  std::vector<ode_step_record> ret = std::move(patch.trajectory);
  patch.trajectory.clear();
  const std::vector<double> sizes = solver.step_sizes();
  util::check_length(sizes.size(), ret.size());
  for (size_t i = 0; i < ret.size(); ++i) {
    ret[i].step_size = sizes[i];
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
template <class Metrics>
std::vector<double> SCM<T, E>::census() const {
  std::vector<double> ret;
  ret.reserve(std::tuple_size<Metrics>::value);
  std::apply(
      [&](auto... psi) -> void {
        (ret.push_back(odelia::util::to_passive(census_over(patch, psi))), ...);
      },
      Metrics{});
  return ret;
}

// One recording of the census over the whole patch at the active scalar, swept
// once per metric. The inputs are the ODE state, so the rows come back in the
// order ode_state writes and can be handed straight to the reverse pass.
//
// set_recorded_state rebuilds the environment and the boundary node from the state
// it is given. Both are on the census's path -- the boundary node is the
// reduction's lower grid point and is not ODE state -- so the recording must carry
// that rebuild. Loading the state without it leaves the boundary node at the values
// it was copied with, and its whole contribution to the seed is then exactly zero
// with nothing thrown. Loading it with set_ode_state alone leaves the condition at
// its first evaluation, which is not the one census() reads.
template <typename T, typename E>
template <class Metrics>
std::vector<std::vector<double>> SCM<T, E>::census_state_adjoint() const {
  require_birth_date_coordinate("census_state_adjoint");
  using scalar = odelia::ode::active_scalar<double>;
  const size_t n_metric = std::tuple_size<Metrics>::value;
  const size_t n_state = patch.ode_size();

  std::vector<double> in(n_state);
  patch.ode_state(in.begin());

  std::vector<std::vector<double>> ret(n_metric,
                                       std::vector<double>(n_state, 0.0));
  for (size_t m = 0; m < n_metric; ++m) {
    // The active patch and the tape are built inside the loop, one per metric,
    // and that is a correctness requirement rather than tidiness.
    //
    // Clearing a tape returns its derivative-slot counter to zero, so an active
    // value constructed OUTSIDE the sweep loop and read inside it refers, after
    // the first clear, to a slot that now belongs to something else. Hoisting
    // the patch out of the loop is the tempting economy -- rebinding it is not
    // free -- and it is the worst failure shape available here: the first metric
    // comes back correct and lends its credibility to the rest, while every
    // later one reads unrelated storage. Measured before this was closed, the
    // second and third metrics' seeds were wrong by three orders and their
    // heartwood columns read exactly zero.
    //
    // Recording once and sweeping many is what a shared tape would buy, and this
    // reduction cannot express it: the states are re-set inside the recording,
    // so each metric's recording is a different one. The saving is not available
    // and pretending otherwise is what produced the defect.
    typename scalar::tape_type tape(false);
    auto active = patch.template rebind_from<scalar>();
    auto reduce = [&](const std::vector<scalar>& x,
                      std::vector<scalar>& y) -> void {
      active.set_recorded_state(x.begin(), time());
      size_t at = 0;
      std::apply(
          [&](auto... psi) -> void {
            ((y[at++] = census_over(active, psi)), ...);
          },
          Metrics{});
    };
    std::vector<double> seed(n_metric, 0.0);
    seed[m] = 1.0;
    odelia::ode::vector_jacobian_product(tape, in, seed, reduce, ret[m]);
  }
  return ret;
}

// The census's own reading of the traits, with the state held. This is not a
// sensitivity of the state, so the sweep below cannot produce it and adding it is
// not double counting: the trajectory term is (dC/dy)^T (dy/dphi), and the
// boundary node -- which a set of the state rebuilds, through the field -- is not
// in y at all.
//
// Recorded rather than written out, so it is the metric algebra that is
// differentiated and a metric added in species.h needs no edit here. One patch
// and one tape per metric, for the reason census_state_adjoint gives.
template <typename T, typename E>
template <class Metrics>
std::vector<std::vector<double>> SCM<T, E>::census_trait_direct() {
  require_birth_date_coordinate("census_trait_direct");
  using scalar = odelia::ode::active_scalar<double>;
  const size_t n_metric = std::tuple_size<Metrics>::value;

  std::vector<double> state(patch.ode_size());
  patch.ode_state(state.begin());

  std::vector<double> in;
  for (size_t i = 0; i < patch.size(); ++i) {
    for (const double* p : patch.at_species(i).strategy_ptr()->ad_parameters()) {
      in.push_back(*p);
    }
  }

  std::vector<std::vector<double>> ret(n_metric,
                                       std::vector<double>(in.size(), 0.0));
  for (size_t m = 0; m < n_metric; ++m) {
    typename scalar::tape_type tape(false);
    auto active = patch.template rebind_from<scalar>();
    auto reduce = [&](const std::vector<scalar>& x,
                      std::vector<scalar>& y) -> void {
      size_t at = 0;
      for (size_t i = 0; i < active.size(); ++i) {
        std::vector<scalar*> pars =
          active.at_species(i).strategy_ptr()->ad_parameters();
        for (size_t p = 0; p < pars.size(); ++p) {
          *pars[p] = x[at++];
        }
      }
      util::check_length(at, x.size());
      // The state is handed in at its value, which is what holds it fixed.
      active.set_recorded_state(state.begin(), time());
      at = 0;
      std::apply(
          [&](auto... psi) -> void {
            ((y[at++] = census_over(active, psi)), ...);
          },
          Metrics{});
    };
    std::vector<double> seed(n_metric, 0.0);
    seed[m] = 1.0;
    odelia::ode::vector_jacobian_product(tape, in, seed, reduce, ret[m]);
  }
  return ret;
}

// The census twice per trait, at the state held, with the strategy moved in place.
// See the declaration for why it perturbs rather than rebuilds.
template <typename T, typename E>
template <class Metrics>
std::vector<std::vector<double>>
SCM<T, E>::census_trait_difference(double rel) {
  require_birth_date_coordinate("census_trait_difference");
  const size_t n_metric = std::tuple_size<Metrics>::value;
  const size_t n_state = patch.ode_size();

  std::vector<double> state(n_state);
  patch.ode_state(state.begin());
  const double time_ = time();

  std::vector<typename T::value_type*> pars;
  for (size_t i = 0; i < patch.size(); ++i) {
    for (typename T::value_type* p :
         patch.at_species(i).strategy_ptr()->ad_parameters()) {
      pars.push_back(p);
    }
  }

  // The state is re-set on every evaluation, which is what makes the moved trait
  // reach the quantities a state determines -- the boundary node among them.
  auto census_at = [&](std::vector<double>& out) -> void {
    patch.set_recorded_state(state.begin(), time_);
    out.clear();
    std::apply(
        [&](auto... psi) -> void {
          (out.push_back(odelia::util::to_passive(census_over(patch, psi))), ...);
        },
        Metrics{});
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
template <class Metrics>
std::vector<std::vector<double>>
SCM<T, E>::census_trait_gradient(const std::vector<size_t>& extra_splits,
                                 const std::vector<size_t>& which_metrics) {
  require_birth_date_coordinate("census_trait_gradient");
  // The sweep needs the state at every accepted step, and store_trajectory()
  // re-runs to record them, so the seeds below are taken after it.
  const std::vector<ode_step_record> trajectory = store_trajectory();
  std::vector<std::vector<double>> states;
  states.reserve(trajectory.size());
  for (const ode_step_record& record : trajectory) {
    states.push_back(record.state);
  }

  std::vector<size_t> boundary;
  std::vector<std::vector<size_t>> introduced;
  narrow_over_introductions(states, trajectory, boundary, introduced);
  patch_type& live = solver.get_system_ref();

  adjoint_segments = 0;
  const std::vector<std::vector<double>> all_seeds =
      census_state_adjoint<Metrics>();
  const std::vector<std::vector<double>> all_direct =
      census_trait_direct<Metrics>();
  // Which rows to sweep. Empty is every one, and a named row outside the census
  // is refused rather than clamped: a caller indexing by position would
  // otherwise get a different metric's gradient back.
  std::vector<size_t> rows = which_metrics;
  if (rows.empty()) {
    rows.resize(all_seeds.size());
    for (size_t m = 0; m < rows.size(); ++m) {
      rows[m] = m;
    }
  }
  std::vector<std::vector<double>> seeds, direct;
  for (const size_t r : rows) {
    if (r >= all_seeds.size()) {
      util::stop("census_trait_gradient: metric " + util::to_string(static_cast<int>(r)) +
                 " is outside the census");
    }
    seeds.push_back(all_seeds[r]);
    direct.push_back(all_direct[r]);
  }
  std::vector<std::vector<double>> sweep_states = states;

  // Every metric's sweep visits the same trajectory and differs only in its
  // seed, so they are carried TOGETHER: a block is recorded once and swept once
  // per metric, where the loop this replaces recorded it once per metric. The
  // recording is a model evaluation and a sweep is arithmetic, so the second and
  // third metrics were costing what the first did and now cost almost nothing.
  const size_t n_metric = seeds.size();
  widen_over_introductions(states, trajectory, boundary, introduced,
                           sweep_states);

  live.clear_trait_adjoint(n_metric);
  std::vector<std::vector<double>> lambda = seeds;
  // One segment per width, highest first: the decomposition the tangent runs
  // forwards, so there is one more segment than there are widenings and the
  // lowest of them runs down to the initial state. Every segment but that one
  // has an introduction at its foot, transposed once the sweep reaches it.
  for (size_t j = boundary.size() + 1; j-- > 0;) {
    const size_t first = j > 0 ? boundary[j - 1] : 0;
    const size_t last = j < boundary.size() ? boundary[j] : states.size() - 1;
    // Stopped and resumed at each requested step inside this segment, highest
    // first, so the pieces compose in the order the whole sweep would take
    // them. With none requested this is the single call it replaces.
    std::vector<size_t> cuts;
    for (const size_t s : extra_splits) {
      if (s > first && s < last) {
        cuts.push_back(s);
      }
    }
    std::sort(cuts.begin(), cuts.end());
    size_t upper = last;
    for (size_t c = cuts.size(); c-- > 0;) {
      solver.solve_adjoint_batched(sweep_states, lambda, cuts[c], upper);
      upper = cuts[c];
      adjoint_segments += n_metric;
    }
    // A widening at the first recorded step leaves the lowest segment with no
    // step in it, which is what a run from bare ground gives.
    if (first < upper) {
      solver.solve_adjoint_batched(sweep_states, lambda, first, upper);
      adjoint_segments += n_metric;
    }
    if (j > 0) {
      live.remove_new_nodes(introduced[j - 1]);
      std::vector<std::vector<double>> narrowed;
      live.introduction_adjoint(introduced[j - 1], states[first],
                                trajectory[first].time, lambda, narrowed);
      lambda = std::move(narrowed);
    }
  }

  std::vector<std::vector<double>> ret;
  ret.reserve(n_metric);
  for (size_t m = 0; m < n_metric; ++m) {
    std::vector<double> row = live.trait_adjoint[m];
    util::check_length(row.size(), direct[m].size());
    for (size_t p = 0; p < row.size(); ++p) {
      row[p] += direct[m][p];
    }
    ret.push_back(row);
  }

  // Leave the system at the width the run left it, so this call is repeatable.
  widen_over_introductions(states, trajectory, boundary, introduced,
                           sweep_states);
  return ret;
}

// A node introduction widens the state between two steps, so a trajectory runs
// one segment per width. states[b] is the state before the introduction at its
// own time: the last step of the segment ended exactly there, since
// advance_adaptive lands its final step on the event time, and the introduction
// followed. So the width changes across b, and step b + 1 starts from the
// widened state, which no record holds.
//
// The species are read newest-first with the system narrowed as it goes, so each
// read sees the nodes that were the newest then -- which leaves the system at the
// first segment's width, where both the sweep and the tangent start.
template <typename T, typename E>
void SCM<T, E>::narrow_over_introductions(
    const std::vector<std::vector<double>>& states,
    const std::vector<ode_step_record>& trajectory,
    std::vector<size_t>& boundary,
    std::vector<std::vector<size_t>>& introduced) {
  if (states.size() < 2) {
    util::stop("no recorded steps to sweep");
  }
  boundary.clear();
  for (size_t k = 1; k < states.size(); ++k) {
    if (states[k].size() != states[k - 1].size()) {
      boundary.push_back(k - 1);
    }
  }
  introduced.assign(boundary.size(), {});
  patch_type& live = solver.get_system_ref();
  for (size_t j = boundary.size(); j-- > 0;) {
    const size_t b = boundary[j];
    introduced[j] = live.nodes_introduced_at(trajectory[b].time);
    live.remove_new_nodes(introduced[j]);
    util::check_length(live.ode_size(), states[b].size());
  }
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
template <class Metrics>
std::vector<double>
SCM<T, E>::census_trait_tangent(const std::vector<double>& direction,
                                std::vector<double>& value) {
  require_birth_date_coordinate("census_trait_tangent");
  using tangent = xad::fwd<double>::active_type;

  const std::vector<ode_step_record> trajectory = store_trajectory();
  std::vector<std::vector<double>> states;
  states.reserve(trajectory.size());
  for (const ode_step_record& record : trajectory) {
    states.push_back(record.state);
  }
  std::vector<size_t> boundary;
  std::vector<std::vector<size_t>> introduced;
  narrow_over_introductions(states, trajectory, boundary, introduced);

  patch_type& live = solver.get_system_ref();
  live.set_ode_state_and_field(states[0].begin(), trajectory[0].time);
  auto active = live.template rebind_from<tangent>();

  // Seeded before the state is set: the quantities a state determines read the
  // parameters, and would otherwise be derived at the unseeded values.
  size_t at = 0;
  for (tangent* p : active.ad_parameters()) {
    if (at >= direction.size()) {
      util::stop("census_trait_tangent: one weight per trait, species-major");
    }
    xad::derivative(*p) = direction[at++];
  }
  util::check_length(direction.size(), at);
  std::vector<tangent> x0(states[0].size());
  for (size_t i = 0; i < x0.size(); ++i) {
    x0[i] = states[0][i];
  }
  active.set_ode_state_and_field(x0.begin(), trajectory[0].time);

  odelia::ode::Solver<decltype(active)> forward(active, make_ode_control(control));
  forward.set_collect(false);
  forward.set_state_from_system();

  size_t first = 0;
  for (size_t j = 0; j <= boundary.size(); ++j) {
    const size_t last =
      j < boundary.size() ? boundary[j] : states.size() - 1;
    // The first entry is the size no step reached, which is how
    // advance_fixed_steps reads a recorded run.
    std::vector<double> sizes(1, std::numeric_limits<double>::quiet_NaN());
    for (size_t k = first + 1; k <= last; ++k) {
      sizes.push_back(trajectory[k].step_size);
    }
    if (sizes.size() > 1) {
      forward.advance_fixed_steps(sizes);
    }
    if (j < boundary.size()) {
      forward.get_system_ref().introduce_new_nodes(introduced[j]);
      forward.set_state_from_system();
      first = last;
    }
  }

  // Leave the double system where the run left it, so this call is repeatable
  // beside the sweep that shares its trajectory.
  std::vector<std::vector<double>> sweep_states = states;
  widen_over_introductions(states, trajectory, boundary, introduced,
                           sweep_states);

  std::vector<double> ret;
  ret.reserve(std::tuple_size<Metrics>::value);
  value.clear();
  value.reserve(std::tuple_size<Metrics>::value);
  const auto& reached = forward.get_system_ref();
  std::apply(
      [&](auto... psi) -> void {
        ((ret.push_back(xad::derivative(census_over(reached, psi))),
          value.push_back(xad::value(census_over(reached, psi)))), ...);
      },
      Metrics{});
  return ret;
}

template <typename T, typename E>
double SCM<T, E>::narrow_to_segment(
    const std::vector<std::vector<double>>& states,
    const std::vector<ode_step_record>& trajectory,
    const std::vector<size_t>& boundary,
    const std::vector<std::vector<size_t>>& introduced,
    size_t from_segment,
    std::vector<double>& base,
    size_t& start) {
  patch_type& live = solver.get_system_ref();
  live.set_ode_state_and_field(states[0].begin(), trajectory[0].time);
  start = 0;
  for (size_t j = 0; j < from_segment; ++j) {
    const size_t b = boundary[j];
    live.set_recorded_state(states[b].begin(), trajectory[b].time);
    live.introduce_new_nodes(introduced[j]);
    start = b;
  }
  base.assign(live.ode_size(), 0.0);
  live.ode_state(base.begin());
  return trajectory[start].time;
}

template <typename T, typename E>
std::vector<double> SCM<T, E>::segment_base_state(size_t segment) {
  const std::vector<ode_step_record> trajectory = store_trajectory();
  std::vector<std::vector<double>> states;
  states.reserve(trajectory.size());
  for (const ode_step_record& record : trajectory) {
    states.push_back(record.state);
  }
  std::vector<size_t> boundary;
  std::vector<std::vector<size_t>> introduced;
  narrow_over_introductions(states, trajectory, boundary, introduced);
  if (segment > boundary.size()) {
    util::stop("segment_base_state: the recording has no such segment");
  }
  std::vector<double> base;
  size_t start = 0;
  narrow_to_segment(states, trajectory, boundary, introduced, segment, base,
                    start);
  // narrow_to_segment introduced its way forward to that segment's width; undo
  // exactly those introductions, because the restore below replays from the
  // first boundary and starts at the narrowest width.
  for (size_t j = segment; j-- > 0;) {
    solver.get_system_ref().remove_new_nodes(introduced[j]);
  }
  std::vector<std::vector<double>> sweep_states = states;
  widen_over_introductions(states, trajectory, boundary, introduced,
                           sweep_states);
  return base;
}

template <typename T, typename E>
template <class Metrics, class Scalar, class Seed>
std::vector<Scalar> SCM<T, E>::replay_initial_state(size_t from_segment,
                                                    Seed seed) {
  const std::vector<ode_step_record> trajectory = store_trajectory();
  std::vector<std::vector<double>> states;
  states.reserve(trajectory.size());
  for (const ode_step_record& record : trajectory) {
    states.push_back(record.state);
  }
  std::vector<size_t> boundary;
  std::vector<std::vector<size_t>> introduced;
  narrow_over_introductions(states, trajectory, boundary, introduced);
  if (from_segment > boundary.size()) {
    util::stop("replay_initial_state: the recording has no such segment");
  }

  std::vector<double> base;
  size_t start = 0;
  const double t0 = narrow_to_segment(states, trajectory, boundary, introduced,
                                      from_segment, base, start);

  patch_type& live = solver.get_system_ref();
  auto active = live.template rebind_from<Scalar>();

  std::vector<Scalar> x0(base.size());
  seed(x0, base);
  active.set_ode_state_and_field(x0.begin(), t0);

  odelia::ode::Solver<decltype(active)> forward(active, make_ode_control(control));
  forward.set_collect(false);
  forward.set_state_from_system();

  size_t first = start;
  for (size_t j = from_segment; j <= boundary.size(); ++j) {
    const size_t last =
      j < boundary.size() ? boundary[j] : states.size() - 1;
    // The first entry is the size no step reached, which is how
    // advance_fixed_steps reads a recorded run.
    std::vector<double> sizes(1, std::numeric_limits<double>::quiet_NaN());
    for (size_t k = first + 1; k <= last; ++k) {
      sizes.push_back(trajectory[k].step_size);
    }
    if (sizes.size() > 1) {
      forward.advance_fixed_steps(sizes);
    }
    if (j < boundary.size()) {
      forward.get_system_ref().introduce_new_nodes(introduced[j]);
      forward.set_state_from_system();
      first = last;
    }
  }

  // Leave the double system where the run left it, so this call is repeatable
  // beside the sweep that shares its trajectory. narrow_to_segment left it at
  // this segment's width and the restore replays from the first boundary, so it
  // needs the narrowest width back first, so undo exactly the introductions
  // narrow_to_segment made on the way forward.
  for (size_t j = from_segment; j-- > 0;) {
    live.remove_new_nodes(introduced[j]);
  }
  std::vector<std::vector<double>> sweep_states = states;
  widen_over_introductions(states, trajectory, boundary, introduced,
                           sweep_states);

  std::vector<Scalar> out;
  out.reserve(std::tuple_size<Metrics>::value);
  const auto& reached = forward.get_system_ref();
  std::apply(
      [&](auto... psi) -> void {
        ((out.push_back(census_over(reached, psi))), ...);
      },
      Metrics{});
  return out;
}

template <typename T, typename E>
template <class Metrics>
std::vector<double>
SCM<T, E>::census_initial_state_tangent(const std::vector<double>& direction,
                                        std::vector<double>& value,
                                        size_t segment) {
  require_birth_date_coordinate("census_initial_state_tangent");
  using tangent = xad::fwd<double>::active_type;

  const std::vector<tangent> reached =
    replay_initial_state<Metrics, tangent>(segment,
      [&](std::vector<tangent>& x0,
          const std::vector<double>& base) -> void {
        util::check_length(direction.size(), base.size());
        for (size_t i = 0; i < x0.size(); ++i) {
          x0[i] = base[i];
          xad::derivative(x0[i]) = direction[i];
        }
      });

  std::vector<double> ret;
  ret.reserve(reached.size());
  value.clear();
  value.reserve(reached.size());
  for (const tangent& metric : reached) {
    ret.push_back(xad::derivative(metric));
    value.push_back(xad::value(metric));
  }
  return ret;
}

template <typename T, typename E>
template <class Metrics>
std::vector<double>
SCM<T, E>::census_initial_state_replay(const std::vector<double>& state0,
                                       size_t segment) {
  require_birth_date_coordinate("census_initial_state_replay");
  return replay_initial_state<Metrics, double>(segment,
    [&](std::vector<double>& x0, const std::vector<double>& base) -> void {
      util::check_length(state0.size(), base.size());
      x0 = state0;
    });
}

// Replay each introduction on the solver's system, which must be narrowed to the
// first boundary's width, and write the state it produced into sweep_states at
// that boundary: the state the block's first step ran from.
template <typename T, typename E>
void SCM<T, E>::widen_over_introductions(
    const std::vector<std::vector<double>>& states,
    const std::vector<ode_step_record>& trajectory,
    const std::vector<size_t>& boundary,
    const std::vector<std::vector<size_t>>& introduced,
    std::vector<std::vector<double>>& sweep_states) {
  patch_type& live = solver.get_system_ref();
  for (size_t j = 0; j < boundary.size(); ++j) {
    const size_t b = boundary[j];
    util::check_length(live.ode_size(), states[b].size());
    live.set_recorded_state(states[b].begin(), trajectory[b].time);
    live.introduce_new_nodes(introduced[j]);
    sweep_states[b].assign(live.ode_size(), 0.0);
    live.ode_state(sweep_states[b].begin());
    util::check_length(sweep_states[b].size(), states[b + 1].size());
  }
}

} // namespace plant

#endif
