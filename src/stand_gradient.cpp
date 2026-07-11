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
// traits: unwrap the resident SCM, lift it to the active scalar, hand over the
// recording for the frozen-canopy (invasion) replay, and drive compute_jacobian.
// Only the invasion workflow is wired; the resident recompute needs the active
// environment.
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

  if (!invasion) {
    util::stop("the resident (self-shading feedback) stand gradient is not "
               "available yet; use invasion_gradient for the frozen-canopy gradient");
  }

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

  // Lift to the active twin and hand over the recording: the resident's step
  // schedule and frozen environment, patch switched to mutant mode, schedule
  // pinned to the recorded ode times (the active replay needs a fixed schedule).
  auto active = resident.template rebind_from<active_scalar>();
  auto& active_patch = active.get_system_ref();
  active_patch.step_history = resident_patch.step_history;
  active_patch.environment_history = resident_patch.environment_history;
  active_patch.set_mutant();

  NodeSchedule sched = active.r_node_schedule();
  sched.r_set_ode_times(resident_patch.step_history);
  sched.r_set_use_ode_times(true);
  sched.reset();
  active.r_set_node_schedule(sched);

  EmergentFunctional functional(metrics);
  auto result = odelia::ode::compute_jacobian(active, targets, functional);
  return pack(result.first, result.second, metrics, traits);
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
