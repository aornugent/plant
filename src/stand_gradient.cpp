#include <plant.h>
#include <plant/emergent_functional.h>
#include <odelia/gradient.hpp>
#include <algorithm>

namespace plant {
namespace gradient {

// Resolve trait names to differentiation-target columns in the AD column order:
// species-major, each species contributing the strategy's field order. The seed
// value of each target is the resident's current parameter value.
odelia::ode::DifferentiationTargets
resolve_targets(const std::vector<std::string>& traits,
                const std::vector<std::string>& field_names,
                size_t species_index,
                const std::vector<double>& param_values) {
  odelia::ode::DifferentiationTargets targets;
  const size_t n_fields = field_names.size();
  for (const auto& trait : traits) {
    auto it = std::find(field_names.begin(), field_names.end(), trait);
    if (it == field_names.end()) {
      util::stop("unknown trait for gradient: " + trait);
    }
    const size_t f = static_cast<size_t>(std::distance(field_names.begin(), it));
    const size_t col = species_index * n_fields + f;
    if (col >= param_values.size()) {
      util::stop("trait column out of range; is the species index valid?");
    }
    targets.params.push_back(static_cast<int>(col));
    targets.values.push_back(param_values[col]);
  }
  return targets;
}

// Pack the metric values and the m x n Jacobian (rows = metrics, cols = traits)
// into an R list, carrying the names so the R wrapper stays thin.
Rcpp::List pack(const std::vector<double>& values,
                const std::vector<std::vector<double>>& jacobian,
                const std::vector<std::string>& metrics,
                const std::vector<std::string>& traits) {
  const size_t m = jacobian.size();
  const size_t n = traits.size();
  Rcpp::NumericMatrix grad(m, n);
  for (size_t i = 0; i < m; ++i) {
    for (size_t j = 0; j < n && j < jacobian[i].size(); ++j) {
      grad(i, j) = jacobian[i][j];
    }
  }
  Rcpp::CharacterVector row_names(metrics.begin(), metrics.end());
  Rcpp::CharacterVector col_names(traits.begin(), traits.end());
  grad.attr("dimnames") = Rcpp::List::create(row_names, col_names);
  Rcpp::NumericVector value(values.begin(), values.end());
  value.attr("names") = row_names;
  return Rcpp::List::create(Rcpp::_["value"] = value, Rcpp::_["gradient"] = grad);
}

// Reverse-mode Jacobian of the requested emergent metrics w.r.t. the requested
// traits: unwrap the resident SCM, lift it (strategy AND environment) to the
// active scalar, hand over the recording, and drive compute_jacobian. The two
// workflows differ by exactly one bit -- whether the canopy is read frozen
// (invasion) or recomputed live on the recorded knots (resident).
template <class Strategy, class Env>
Rcpp::List gradient_impl(SEXP scm_sexp,
                         const std::vector<std::string>& metrics,
                         const std::vector<std::string>& traits,
                         size_t species_index,
                         const std::vector<std::string>& field_names,
                         bool invasion) {
  using SCM_d = plant::SCM<Strategy, Env>;
  using patch_type = plant::Patch<Strategy, Env>;
  using active_scalar = typename odelia::ode::Solver<patch_type>::active_scalar;

  plant::RcppR6::RcppR6<SCM_d> handle(scm_sexp);
  SCM_d& resident = *handle;
  patch_type& resident_patch = resident.get_system_ref();

  std::vector<double> param_values;
  {
    std::vector<double*> ptrs = resident_patch.ad_parameters();
    param_values.reserve(ptrs.size());
    for (double* p : ptrs) param_values.push_back(*p);
  }
  odelia::ode::DifferentiationTargets targets =
      resolve_targets(traits, field_names, species_index, param_values);

  // Lift to the active twin (strategy + environment) and pin the schedule to the
  // recorded ode times -- both replays need a fixed schedule.
  auto active = resident.template rebind_from<active_scalar>();
  auto& active_patch = active.get_system_ref();
  active_patch.step_history = resident_patch.step_history;

  NodeSchedule sched = active.r_node_schedule();
  sched.r_set_ode_times(resident_patch.step_history);
  sched.r_set_use_ode_times(true);
  sched.reset();
  active.r_set_node_schedule(sched);

  if (invasion) {
    // Frozen-canopy replay: lift the recorded resident environments to the active
    // scalar as passive constants (derivative through the canopy is zero) and
    // switch to mutant mode. L3 populated -> derivs reads the field frozen.
    active_patch.environment_history.reserve(
        resident_patch.environment_history.size());
    for (const auto& row : resident_patch.environment_history) {
      typename std::decay_t<decltype(active_patch.environment_history)>::value_type
          arow;
      arow.reserve(row.size());
      for (const auto& e : row) {
        arow.push_back(e.template rebind_from<active_scalar>());
      }
      active_patch.environment_history.push_back(std::move(arow));
    }
    active_patch.set_mutant();
  } else {
    // Resident (self-shading) replay: hand over the recorded knot POSITIONS and
    // rebuild the canopy on them with the active cohorts. L3 stays empty
    // (has_recorded_field() false), so a trait re-shades the stand through
    // area_leaf and the self-shading cross term flows.
    active_patch.knot_history = resident_patch.knot_history;
    active_patch.set_resident_replay();
  }

  EmergentFunctional functional(metrics);
  auto result = odelia::ode::compute_jacobian(active, targets, functional);
  return pack(result.first, result.second, metrics, traits);
}

// Reverse-mode gradient of the requested metrics w.r.t. a focal species' birth
// rate, on the coupled resident replay (AD-8): the birth-rate scale is a single
// registered initial-state leaf, seeded to `scale`, that multiplies the extrinsic
// birth rate at every node introduction. It flows through the cohort density and,
// because the canopy is recomputed live on the recorded knots, through the whole
// stand -- so the feedback axis is present and can flip the sign of biomass.
//
// The seed multiplies the base rate, so d/d(scale) = base * d/d(birth_rate); the
// reverse-mode row is divided by the base rate to report the derivative w.r.t. the
// birth rate itself. `value` is the metric at the seeded scale, so a caller can
// finite-difference the same frozen-knot replay by perturbing `scale`.
template <class Strategy, class Env>
Rcpp::List birth_rate_gradient_impl(SEXP scm_sexp,
                                    const std::vector<std::string>& metrics,
                                    size_t species_index,
                                    double scale) {
  using SCM_d = plant::SCM<Strategy, Env>;
  using patch_type = plant::Patch<Strategy, Env>;
  using active_scalar = typename odelia::ode::Solver<patch_type>::active_scalar;

  plant::RcppR6::RcppR6<SCM_d> handle(scm_sexp);
  SCM_d& resident = *handle;
  patch_type& resident_patch = resident.get_system_ref();

  if (species_index >= resident_patch.size()) {
    util::stop("species index out of range for the birth-rate gradient");
  }
  const double base = resident_patch.at_species(species_index)
                          .extrinsic_drivers()
                          .evaluate("birth_rate", 0.0);
  if (!(base > 0.0)) {
    util::stop("birth-rate gradient needs a positive base birth rate");
  }

  odelia::ode::DifferentiationTargets targets;
  targets.ics.push_back(static_cast<int>(species_index));
  targets.values.push_back(scale);

  auto active = resident.template rebind_from<active_scalar>();
  auto& active_patch = active.get_system_ref();
  active_patch.step_history = resident_patch.step_history;

  NodeSchedule sched = active.r_node_schedule();
  sched.r_set_ode_times(resident_patch.step_history);
  sched.r_set_use_ode_times(true);
  sched.reset();
  active.r_set_node_schedule(sched);

  active_patch.knot_history = resident_patch.knot_history;
  active_patch.set_resident_replay();

  EmergentFunctional functional(metrics);
  auto result = odelia::ode::compute_jacobian(active, targets, functional);

  const auto& values = result.first;
  const auto& jac = result.second;
  Rcpp::NumericVector value(values.begin(), values.end());
  Rcpp::NumericVector grad(metrics.size());
  for (size_t i = 0; i < metrics.size(); ++i) {
    grad[i] = (i < jac.size() && !jac[i].empty()) ? jac[i][0] / base : NA_REAL;
  }
  Rcpp::CharacterVector names(metrics.begin(), metrics.end());
  value.attr("names") = names;
  grad.attr("names") = names;
  return Rcpp::List::create(Rcpp::_["value"] = value,
                            Rcpp::_["gradient"] = grad,
                            Rcpp::_["birth_rate"] = base);
}

} // namespace gradient
} // namespace plant

// [[Rcpp::export]]
Rcpp::List invasion_gradient_cpp(SEXP scm,
                                 std::vector<std::string> metrics,
                                 std::vector<std::string> traits,
                                 int species,
                                 std::string strategy) {
  const size_t sp = species > 0 ? static_cast<size_t>(species - 1) : 0;
  if (strategy == "FF16") {
    return plant::gradient::gradient_impl<plant::FF16_Strategy, plant::FF16_Environment>(
        scm, metrics, traits, sp, plant::FF16_Pars::field_names(), true);
  }
  plant::util::stop("stand gradients are only available for the FF16 strategy; got: " +
                    strategy);
  return Rcpp::List();
}

// [[Rcpp::export]]
Rcpp::List birth_rate_gradient_cpp(SEXP scm,
                                   std::vector<std::string> metrics,
                                   int species,
                                   std::string strategy,
                                   double scale) {
  const size_t sp = species > 0 ? static_cast<size_t>(species - 1) : 0;
  if (strategy == "FF16") {
    return plant::gradient::birth_rate_gradient_impl<plant::FF16_Strategy,
                                                     plant::FF16_Environment>(
        scm, metrics, sp, scale);
  }
  plant::util::stop("birth-rate gradients are only available for the FF16 strategy; got: " +
                    strategy);
  return Rcpp::List();
}

// [[Rcpp::export]]
Rcpp::List stand_gradient_cpp(SEXP scm,
                              std::vector<std::string> metrics,
                              std::vector<std::string> traits,
                              int species,
                              std::string strategy) {
  const size_t sp = species > 0 ? static_cast<size_t>(species - 1) : 0;
  if (strategy == "FF16") {
    return plant::gradient::gradient_impl<plant::FF16_Strategy, plant::FF16_Environment>(
        scm, metrics, traits, sp, plant::FF16_Pars::field_names(), false);
  }
  plant::util::stop("stand gradients are only available for the FF16 strategy; got: " +
                    strategy);
  return Rcpp::List();
}
