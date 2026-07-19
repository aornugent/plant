// [[Rcpp::plugins(cpp20)]]
// FF16 SCM R0 (offspring) gradient: the reverse-mode trait gradient of lifetime
// offspring production over the whole FF16 method-of-characteristics solve.
//
// Reuses the generic v2 machinery proven on K93 with no strategy-specific
// additions: the SCM duck-types the odelia gradient-driver contract, the active
// pass replays the L1 ode-time schedule recorded on a double run (resident L2, no
// set_mutant), and offspring flows through the value_type reproduction chain.
//
// The replay uses the RESOLVED schedule (L0 node_schedule_times + L1 ode_times)
// from a refined base run, passed in from R (run_scm(refine_schedule=TRUE)), so it
// replays EXACTLY what run_scm(use_ode_times=TRUE) does. Two settled facts:
//
//   (1) NO schedule sensitivity. The frozen replay on the resolved schedule matches
//       the fully-adaptive model: d(offspring)/d(lma) FD = +4.2 both ways (measured
//       in double at R level). The earlier "schedule sensitivity" and "r_ode_times
//       alone fixes it" claims are RETIRED -- the fix is the resolved L0+L1, and the
//       old drivers were wrong because they pinned only L1 onto the default L0.
//   (2) OPEN adjoint bug: reverse AD != FD on the IDENTICAL resolved schedule
//       (delta-independent, both AD modes agree with each other yet disagree with
//       FD). It reproduces on metric=2 (pure growth), isolating it to FF16's coupled
//       self-shading growth trajectory (NOT reproduction/census, NOT the query
//       channel). So the code computes an analytically wrong derivative FF16's rate/
//       field path drops; FD catches it. This is the remaining work for a correct
//       FF16 gradient. `freeze_query` / `freeze_field` / `metric` / `fd_rel` are the
//       channel-isolation diagnostics.
//
// Tape memory: reverse AD fits to ~life 40; life 50 exceeds memory (checkpointing at
// the node-introduction boundary is deferred).
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
Rcpp::List ff16_scm_offspring_gradient(Rcpp::List node_sched, std::vector<double> ode_times_in,
                                       double lma = 0.1978791, double birth_rate = 20.0,
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

  // The RESOLVED schedule (L0 node-introduction times + L1 ode times) is computed
  // by the R harness via run_scm(refine_schedule=TRUE) and passed in, so the active
  // and FD runs replay EXACTLY what run_scm(use_ode_times=TRUE) does. Pinning only
  // L1 onto the DEFAULT (unrefined) L0 is an inconsistent schedule that gives a
  // value-correct but derivative-wrong trajectory.
  std::vector<std::vector<double>> resolved_node_times;
  for (R_xlen_t i = 0; i < node_sched.size(); ++i)
    resolved_node_times.push_back(Rcpp::as<std::vector<double>>(node_sched[i]));
  const std::vector<double> resolved_ode_times = ode_times_in;

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
    p.node_schedule_times = resolved_node_times;  // L0
    p.ode_times           = resolved_ode_times;   // L1
    p.validate();
    return p;
  };

  Control ctrl_replay;

  // Double replay on the resolved schedule -> the reference offspring value.
  double offspring_double = 0.0;
  {
    FF16_Environment_<double> env;
    SCM<FF16_Strategy_<double>, FF16_Environment_<double>> scm(
        make_params(FF16_Strategy_<double>()), env, ctrl_replay);
    NodeSchedule ns = scm.r_node_schedule();
    ns.r_set_use_ode_times(true);
    scm.r_set_node_schedule(ns);
    scm.run();
    offspring_double = scm.get_system_ref().offspring_production()[0];
  }

  // The schedule is already loaded from parameters; reset and flip on ode-time replay.
  auto pin_replay = [&](auto& scm) {
    scm.reset();
    NodeSchedule ns = scm.r_node_schedule();
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
    p.strategies.push_back(s); p.max_patch_lifetime = max_patch_lifetime;
    p.node_schedule_times = resolved_node_times;  // resolved L0
    p.ode_times           = resolved_ode_times;    // resolved L1
    p.validate();
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
      Rcpp::Named("n_steps")          = double(resolved_ode_times.size()));
}

// Per-cohort localisation: forward-mode d(height_i)/d(lma) (AD tangent) vs a
// per-cohort central FD, both on the resolved schedule. Which cohorts diverge
// (all uniformly vs deep-shade only) pins where the coupled-growth derivative is
// dropped.
// [[Rcpp::export]]
Rcpp::List ff16_cohort_height_tangents(Rcpp::List node_sched,
                                       std::vector<double> ode_times_in,
                                       double lma = 0.1978791, double birth_rate = 20.0,
                                       double max_patch_lifetime = 40.0) {
  using FwdS = xad::fwd<double>::active_type;
  std::vector<std::vector<double>> nst;
  for (R_xlen_t i = 0; i < node_sched.size(); ++i)
    nst.push_back(Rcpp::as<std::vector<double>>(node_sched[i]));

  auto make_params = [&](auto tag, double lma_v) {
    using Strat = std::decay_t<decltype(tag)>; using Env = typename Strat::environment_type;
    Strat strat; set_field(strat, "lma", lma_v);
    strat.is_variable_birth_rate = false; strat.birth_rate_y = {birth_rate};
    Parameters<Strat, Env> p; p.strategies.push_back(strat);
    p.max_patch_lifetime = max_patch_lifetime;
    p.node_schedule_times = nst; p.ode_times = ode_times_in; p.validate();
    return p;
  };
  Control ctrl;
  auto pin = [&](auto& scm){ scm.reset(); NodeSchedule ns=scm.r_node_schedule();
                             ns.r_set_use_ode_times(true); scm.r_set_node_schedule(ns); };

  // Forward AD: seed lma tangent, run, read each cohort's height value + tangent.
  std::vector<double> h_val, h_tan, grate;
  {
    FF16_Environment_<FwdS> env;
    SCM<FF16_Strategy_<FwdS>, FF16_Environment_<FwdS>> scm(make_params(FF16_Strategy_<FwdS>(), lma), env, ctrl);
    auto ptrs = scm.get_system_ref().ad_parameters();
    xad::derivative(*ptrs[0]) = 1.0;   // lma is FF16 field index 0
    pin(scm); scm.run();
    auto& p = scm.get_system_ref();
    for (std::size_t s = 0; s < p.size(); ++s) {
      auto& sp = p.at_species(s);
      for (auto it = sp.node_begin(); it != sp.node_end(); ++it) {
        h_val.push_back(xad::value(it->height()));
        h_tan.push_back(xad::derivative(it->height()));
        // Growth rate: net<=0 forces this to 0 (the net_mass_production_dt>0
        // branch), so ~0 flags a cohort at/below the kink.
        grate.push_back(xad::value(it->individual.rate(HEIGHT_INDEX)));
      }
    }
  }
  // Per-cohort central FD on the same schedule (same introductions => same order).
  auto heights_at = [&](double lma_v) {
    FF16_Environment_<double> env;
    SCM<FF16_Strategy_<double>, FF16_Environment_<double>> scm(make_params(FF16_Strategy_<double>(), lma_v), env, ctrl);
    pin(scm); scm.run();
    std::vector<double> h; auto& p = scm.get_system_ref();
    for (std::size_t s = 0; s < p.size(); ++s) { auto& sp = p.at_species(s);
      for (auto it = sp.node_begin(); it != sp.node_end(); ++it) h.push_back(xad::value(it->height())); }
    return h;
  };
  const double d = 3e-4 * (std::abs(lma) + 1e-3);
  auto hp = heights_at(lma + d), hm = heights_at(lma - d);
  std::vector<double> h_fd(hp.size());
  for (std::size_t i = 0; i < hp.size() && i < hm.size(); ++i) h_fd[i] = (hp[i] - hm[i]) / (2.0 * d);

  return Rcpp::List::create(
      Rcpp::Named("height")   = Rcpp::wrap(h_val),
      Rcpp::Named("ad_tan")   = Rcpp::wrap(h_tan),
      Rcpp::Named("grate")    = Rcpp::wrap(grate),
      Rcpp::Named("fd")       = Rcpp::wrap(h_fd));
}
