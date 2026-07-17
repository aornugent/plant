// [[Rcpp::plugins(cpp20)]]
// FF16 SCM R0 (offspring) gradient: the exact reverse-mode trait gradient of
// lifetime offspring production over the whole FF16 method-of-characteristics
// solve, verified by the JVP=VJP dot-product oracle.
//
// This reuses the generic v2 machinery proven on K93 with no strategy-specific
// additions: the SCM duck-types the odelia gradient-driver contract, the active
// pass replays the L1 ode-time schedule recorded on a double run (resident L2, no
// set_mutant), and offspring flows through the value_type reproduction chain.
//
// FF16's density transport is the FD upwind stencil (the geometric mass chart is
// numerically unstable for FF16 -- it drives cohort density to overflow -- so
// FF16 does not declare the geometric_transport marker and ignores the Control
// flag). The stencil's transport term carries no parameter-derivative on the
// active pass, but offspring does not route through it: the fitness integral
// weights by birth-time density (patch_density_at_birth) and accumulated
// survival-weighted fecundity, not the transported cohort density. So R0 is
// exact here; a density-weighted CENSUS would need a differentiable transport and
// is deferred.
#include <Rcpp.h>
#include <plant/models/ff16_strategy.h>
#include <plant/models/ff16_environment.h>
#include <plant/individual.h>
#include <plant/parameters.h>
#include <plant/control.h>
#include <plant/scm.h>
#include <odelia/gradient.hpp>
#include <string>
#include <vector>

using namespace plant;

// Set one named FF16 field on a strategy of any scalar S (value in, narrowed).
template <class S>
static void set_field(FF16_Strategy_<S>& s, const std::string& name, double v) {
  auto names = FF16_Strategy_<S>::field_names();
  auto ptrs  = s.field_ptrs();
  for (std::size_t i = 0; i < names.size(); ++i)
    if (names[i] == name) { *ptrs[i] = S(v); return; }
  Rcpp::stop("unknown FF16 field: " + name);
}

// Target traits for the offspring gradient (FF16_AD_FIELDS order): lma (leaf mass
// per area, 0), a_l1 (leaf-area allometry, 6), k_l (leaf turnover, 16). All feed
// the fitness integral through growth/allocation.
static const std::vector<int> FF16_R0_TARGET_IDX = {0, 6, 16};

// [[Rcpp::export]]
Rcpp::List ff16_scm_offspring_gradient(double lma = 0.1978791, double birth_rate = 20.0,
                                       double max_patch_lifetime = 105.32) {
  using RevS = xad::adj<double>::active_type;  // AReal: compute_gradient
  using FwdS = xad::fwd<double>::active_type;  // FReal: compute_jvp

  auto make_params = [&](auto strat_tag) {
    using Strat = std::decay_t<decltype(strat_tag)>;
    using Env   = typename Strat::environment_type;
    Strat strat;
    set_field(strat, "lma", lma);
    strat.is_variable_birth_rate = false;
    strat.birth_rate_y = {birth_rate};
    Parameters<Strat, Env> p;
    p.strategies.push_back(strat);
    p.max_patch_lifetime = max_patch_lifetime;
    p.validate();
    return p;
  };

  // Recording Control (double run) records the L1 ode-time schedule; the active
  // runs replay it and recompute the light field at the active scalar.
  Control ctrl_record;
  ctrl_record.save_RK45_cache = true;
  Control ctrl_replay;

  std::vector<double> step_history;
  double offspring_double = 0.0;
  {
    FF16_Environment_<double> env;
    SCM<FF16_Strategy_<double>, FF16_Environment_<double>> scm(
        make_params(FF16_Strategy_<double>()), env, ctrl_record);
    scm.run();
    step_history     = scm.get_system_ref().step_history;
    offspring_double = scm.get_system_ref().offspring_production()[0];
  }

  auto pin_replay = [&](auto& scm) {
    scm.reset();
    NodeSchedule ns = scm.r_node_schedule();
    ns.r_set_ode_times(step_history);
    ns.r_set_use_ode_times(true);
    scm.r_set_node_schedule(ns);
  };
  auto targets_for = [&](auto& scm) {
    odelia::ode::DifferentiationTargets t;
    t.params = FF16_R0_TARGET_IDX;
    auto ptrs = scm.get_system_ref().ad_parameters();
    for (int i : FF16_R0_TARGET_IDX) t.values.push_back(xad::value(*ptrs[i]));
    return t;
  };

  // Reverse leg: value + gradient g = Jᵀ·1 over an AReal SCM.
  double value = 0.0;
  std::vector<double> g;
  {
    FF16_Environment_<RevS> env;
    SCM<FF16_Strategy_<RevS>, FF16_Environment_<RevS>> scm(
        make_params(FF16_Strategy_<RevS>()), env, ctrl_replay);
    pin_replay(scm);
    auto functional = [](decltype(scm)& s) -> RevS {
      return s.get_system_ref().offspring_production()[0];
    };
    auto res = odelia::ode::compute_gradient(scm, targets_for(scm), functional);
    value = res.first;
    g     = res.second;
  }

  // Forward leg: directional derivative J v; the oracle wants it to equal dot(g, v).
  std::vector<double> v(FF16_R0_TARGET_IDX.size());
  for (std::size_t i = 0; i < v.size(); ++i) v[i] = 0.3 + 0.15 * double(i);
  double jvp = 0.0;
  {
    FF16_Environment_<FwdS> env;
    SCM<FF16_Strategy_<FwdS>, FF16_Environment_<FwdS>> scm(
        make_params(FF16_Strategy_<FwdS>()), env, ctrl_replay);
    pin_replay(scm);
    auto functional = [](decltype(scm)& s) -> FwdS {
      return s.get_system_ref().offspring_production()[0];
    };
    auto dd = odelia::ode::compute_directional_derivative(scm, targets_for(scm), v,
                                                          functional);
    jvp = dd.second;
  }

  double dot = 0.0;
  for (std::size_t i = 0; i < g.size(); ++i) dot += g[i] * v[i];

  return Rcpp::List::create(
      Rcpp::Named("value")            = value,
      Rcpp::Named("offspring_double") = offspring_double,
      Rcpp::Named("grad")             = Rcpp::wrap(g),
      Rcpp::Named("jvp")              = jvp,
      Rcpp::Named("dot")              = dot,
      Rcpp::Named("n_steps")          = double(step_history.size()));
}
