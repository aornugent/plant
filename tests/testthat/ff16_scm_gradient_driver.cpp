// [[Rcpp::plugins(cpp20)]]
// FF16 SCM R0 (offspring) gradient: the reverse-mode trait gradient of lifetime
// offspring production over the whole FF16 method-of-characteristics solve.
//
// Reuses the generic v2 machinery proven on K93 with no strategy-specific
// additions: the SCM duck-types the odelia gradient-driver contract, the active
// pass replays the L1 ode-time schedule recorded on a double run (resident L2, no
// set_mutant), and offspring flows through the value_type reproduction chain.
//
// The replay grid is the SOLVER-OWNED L1 schedule (SCM::r_ode_times() ==
// solver.times()), the single source odelia's AD workflow keeps consistent -- NOT
// patch.step_history (the save_RK45_cache / run_mutant L3 legacy record, a different
// grid). This is the correct grid regardless; it does NOT by itself make FF16's R0
// gradient correct.
//
// OPEN (2026-07-19): on this grid the reverse AD (+729) and the pinned-schedule FD
// on the SAME grid (-255) disagree, while the fully-adaptive real-model FD (ground
// truth, measured at R level) is +4.2. Two facts follow: (1) frozen-schedule FD
// (-255) != adaptive FD (+4.2), so NO frozen-replay gradient reaches the true value
// for FF16 (unlike K93, where all three agree); (2) AD != FD on the IDENTICAL frozen
// grid, a derivative discrepancy independent of schedule. It reproduces on metric=2
// (pure growth, sum of heights), isolating it to the coupled growth trajectory, not
// reproduction/census. Under investigation. `freeze_query` / `metric` / `fd_rel`
// (FD step, sweep for the plateau) are the channel-isolation diagnostics.
#include <Rcpp.h>
#include <plant/models/ff16_strategy.h>
#include <plant/models/ff16_environment.h>
#include <plant/individual.h>
#include <plant/parameters.h>
#include <plant/control.h>
#include <plant/scm.h>
#include <odelia/gradient.hpp>
#include <cmath>
#include <string>
#include <vector>

using namespace plant;

// Reduce the positioned patch to a scalar by metric:
//   0 = offspring (R0, through reproduction); 1 = census (compute_competition(0),
//   density-weighted); 2 = sum of cohort heights (PURE GROWTH -- no density, no
//   reproduction, no field-at-0, so it isolates the growth trajectory / RK
//   integration from everything downstream).
template <class Patch>
static auto reduce_metric(Patch& p, int metric)
    -> std::decay_t<decltype(p.compute_competition(0.0))> {
  using V = std::decay_t<decltype(p.compute_competition(0.0))>;
  if (metric == 1) return p.compute_competition(0.0);
  if (metric == 2) {
    V h = 0.0;
    for (std::size_t i = 0; i < p.size(); ++i) {
      auto& sp = p.at_species(i);
      for (auto it = sp.node_begin(); it != sp.node_end(); ++it) h += it->height();
    }
    return h;
  }
  return p.offspring_production()[0];
}

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
                                       double max_patch_lifetime = 105.32,
                                       bool freeze_query = false, int metric = 0,
                                       double fd_rel = 3e-4, bool freeze_field = false) {
  using RevS = xad::adj<double>::active_type;  // AReal: compute_gradient
  using FwdS = xad::fwd<double>::active_type;  // FReal: compute_jvp

  // Channel isolation: freeze the field's query-height derivative on the reverse
  // type only, so the reverse gradient carries the field's SOURCE self-shading but
  // not its query-height feedback (what the fitted spline also does). Comparing to
  // the full field gradient splits the two channels.
  FF16_Environment_<RevS>::freeze_query_derivative = freeze_query;
  FF16_Environment_<RevS>::freeze_field_derivative = freeze_field;

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

  std::vector<double> ode_times;
  double offspring_double = 0.0;
  {
    FF16_Environment_<double> env;
    SCM<FF16_Strategy_<double>, FF16_Environment_<double>> scm(
        make_params(FF16_Strategy_<double>()), env, ctrl_record);
    scm.run();
    // The replay grid is the SOLVER-OWNED L1 schedule (solver.times()), the single
    // source odelia's AD workflow guarantees consistent -- NOT patch.step_history
    // (the save_RK45_cache / run_mutant L3 legacy record, a separate grid).
    ode_times        = scm.r_ode_times();
    offspring_double = scm.get_system_ref().offspring_production()[0];
  }

  auto pin_replay = [&](auto& scm) {
    scm.reset();
    NodeSchedule ns = scm.r_node_schedule();
    ns.r_set_ode_times(ode_times);
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
    auto functional = [metric](decltype(scm)& s) -> RevS {
      return reduce_metric(s.get_system_ref(), metric);
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
    auto functional = [metric](decltype(scm)& s) -> FwdS {
      return reduce_metric(s.get_system_ref(), metric);
    };
    auto dd = odelia::ode::compute_directional_derivative(scm, targets_for(scm), v,
                                                          functional);
    jvp = dd.second;
  }

  double dot = 0.0;
  for (std::size_t i = 0; i < g.size(); ++i) dot += g[i] * v[i];

  // The trustworthy gate: a pinned-schedule central FD of offspring over the
  // DOUBLE SCM (perturb each target on the SAME recorded ode schedule -- no
  // adaptive-schedule noise). Unlike the JVP=VJP oracle, this is a check against
  // the model, so it catches the shared-representation-bias class the oracle hid
  // for K93. The caller sweeps delta for the stable plateau.
  auto offspring_pinned = [&](int field_idx, double d) -> double {
    FF16_Strategy_<double> s;
    set_field(s, "lma", lma);
    s.is_variable_birth_rate = false; s.birth_rate_y = {birth_rate};
    *s.field_ptrs()[field_idx] += d;
    Parameters<FF16_Strategy_<double>, FF16_Environment_<double>> p;
    p.strategies.push_back(s); p.max_patch_lifetime = max_patch_lifetime; p.validate();
    FF16_Environment_<double> env;
    SCM<FF16_Strategy_<double>, FF16_Environment_<double>> scm(p, env, ctrl_replay);
    pin_replay(scm); scm.run();
    return reduce_metric(scm.get_system_ref(), metric);
  };
  std::vector<double> fd_grad(FF16_R0_TARGET_IDX.size());
  {
    std::vector<double> tv;
    { FF16_Strategy_<double> s0; set_field(s0, "lma", lma);
      for (int i : FF16_R0_TARGET_IDX) tv.push_back(*s0.field_ptrs()[i]); }
    for (std::size_t i = 0; i < FF16_R0_TARGET_IDX.size(); ++i) {
      const double dl = fd_rel * (std::abs(tv[i]) + 1e-3);  // caller sweeps for the plateau
      fd_grad[i] = (offspring_pinned(FF16_R0_TARGET_IDX[i], dl) -
                    offspring_pinned(FF16_R0_TARGET_IDX[i], -dl)) / (2.0 * dl);
    }
  }

  return Rcpp::List::create(
      Rcpp::Named("value")            = value,
      Rcpp::Named("offspring_double") = offspring_double,
      Rcpp::Named("grad")             = Rcpp::wrap(g),
      Rcpp::Named("jvp")              = jvp,
      Rcpp::Named("dot")              = dot,
      Rcpp::Named("fd_grad")          = Rcpp::wrap(fd_grad),
      Rcpp::Named("n_steps")          = double(ode_times.size()));
}
