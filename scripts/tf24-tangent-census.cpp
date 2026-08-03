// A whole-run TANGENT (forward-mode) reference for the TF24 census: one trait
// seeded, the recorded ODE schedule replayed at the tangent scalar, and the
// census read through the same reduction the double run uses. This is the
// reference the adjoint's trait gradient is judged against.
//
//   Rscript scripts/tangent-reference-driver.R <plant worktree>
//
// Two limits the reference shares with the adjoint, so a pass over them is
// not evidence:
//
//   The leaf's supplied partials are common to both. graft_leaf_outputs fires
//   for FReal and AReal alike and both read the same double rows out of
//   Leaf::input_adjoints, so a wrong row cancels between the two and this
//   comparison cannot see it.
//
//   Both drop d(height_0)/d(trait). rebind_from carries height_0,
//   area_leaf_0 and eta_c across as values, because prepare_strategy and
//   height_seed refuse an active scalar. This is why the tangent and a
//   converged central difference differ by the same amount at every step size
//   -- 2.6% for lma -- rather than converging: the difference carries that
//   channel and neither derivative does.
//
// SCM<TF24_Strategy<FReal>, ...> is not constructible: its constructor runs
// add_strategies -> prepare_strategy, whose static assertions refuse an active
// scalar. So the prepared double patch is taken through r_patch(), returned to
// t = 0, rebound, and the solver is driven here.
#include <plant.h>
#include <odelia/ode_solver.hpp>
#include <XAD/XAD.hpp>
#include <string>
#include <vector>

using dstrat = plant::TF24_Strategy<double>;
using denv = plant::TF24_Environment<double>;
using dscm = plant::SCM<dstrat, denv>;
using dpatch = plant::Patch<dstrat, denv>;

using F = typename xad::fwd<double>::active_type;
using fstrat = plant::TF24_Strategy<F>;
using fenv = plant::TF24_Environment<F>;
using fpatch = plant::Patch<fstrat, fenv>;

// [[Rcpp::export]]
Rcpp::List tangent_census_tf24(plant::RcppR6::RcppR6<dscm> obj_,
                               std::string trait,
                               Rcpp::List control_list) {

  // rebind_from must start from a prepared double patch, and at t = 0.
  dpatch dp = obj_->r_patch();
  dp.reset();

  // The trait is seeded on the same ad_parameters() address the adjoint
  // accumulates into, so the two are transposes of one map.
  fpatch fp = dp.template rebind_from<F>();
  int seeded = 0;
  for (size_t i = 0; i < fp.size(); ++i) {
    fstrat* s = fp.at_species(i).strategy_ptr().get();
    const std::vector<std::string> nm = s->ad_parameter_names();
    std::vector<F*> ad = s->ad_parameters();
    for (size_t j = 0; j < nm.size(); ++j) {
      if (nm[j] == trait) { xad::derivative(*ad[j]) = 1.0; ++seeded; }
    }
  }
  if (seeded == 0) Rcpp::stop("trait '" + trait + "' is not an ad_parameter");

  // The replayed schedule: this run's own introduction times, ODE times and
  // accepted step sizes. Without the sizes the tangent run takes its own steps
  // and the comparison is between two different discretisations.
  plant::NodeSchedule sched = obj_->r_node_schedule();
  sched.r_set_ode_times(obj_->r_ode_times());
  sched.r_set_ode_step_sizes(obj_->r_ode_step_sizes());
  sched.r_set_use_ode_times(true);
  sched.reset();
  if (!sched.using_ode_times()) Rcpp::stop("schedule did not pin");

  plant::Control ctrl = Rcpp::as<plant::Control>(control_list);
  odelia::ode::Solver<fpatch> solver(fp, plant::make_ode_control(ctrl));
  solver.set_collect(true);
  solver.reset();

  int n_intro = 0;
  while (sched.remaining() > 0) {
    const double t0 = solver.time();
    plant::NodeSchedule::Event e = sched.next_event();
    std::vector<size_t> add;
    while (true) {
      if (!plant::util::identical(t0, e.time_introduction()))
        Rcpp::stop("start time not what was expected: t0=" +
                   std::to_string(t0) + " event=" +
                   std::to_string(e.time_introduction()));
      add.push_back(e.species_index);
      sched.pop();
      if (e.time_end() > t0 || sched.remaining() == 0) break;
      e = sched.next_event();
    }
    solver.get_system_ref().introduce_new_nodes(add);
    n_intro += static_cast<int>(add.size());
    solver.set_state_from_system();
    if (e.step_sizes.empty()) {
      solver.advance_fixed(e.times);
    } else {
      solver.advance_fixed_steps(e.step_sizes);
      solver.advance_fixed({solver.time(), e.times.back()});
    }
  }

  // SCM::census_over is templated on the patch, so the tangent run is reduced
  // through the same code the double run is.
  const fpatch& live = solver.get_system_ref();
  std::vector<double> value, deriv;
  auto one = [&](auto psi) -> void {
    const F c = dscm::census_over(live, psi);
    value.push_back(xad::value(c));
    deriv.push_back(xad::derivative(c));
  };
  std::apply([&](auto... psi) -> void { (one(psi), ...); },
             plant::tf24_census{});

  return Rcpp::List::create(
      Rcpp::_["value"] = value,
      Rcpp::_["tangent"] = deriv,
      Rcpp::_["time"] = solver.time(),
      Rcpp::_["ode_size"] = static_cast<int>(live.ode_size()),
      Rcpp::_["n_introduced"] = n_intro,
      Rcpp::_["accepted_steps"] = static_cast<int>(solver.times().size()));
}
