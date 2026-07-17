// [[Rcpp::plugins(cpp20)]]
// Layer (a) of the K93 SCM census gradient (CD-G): build and run a double K93
// SCM entirely in C++ (there is no C++ SCM-construction path today; params are
// set up R-side), then reduce a census metric. Verified against R's run_scm so
// the C++ construction is trustworthy before the active gradient is layered on.
#include <Rcpp.h>
#include <plant/models/k93_strategy.h>
#include <plant/models/k93_environment.h>
#include <plant/individual.h>
#include <plant/parameters.h>
#include <plant/control.h>
#include <plant/scm.h>
#include <odelia/gradient.hpp>
#include <vector>

using namespace plant;

// [[Rcpp::export]]
Rcpp::List k93_scm_census_double(double b_0 = 0.059, double birth_rate = 20.0,
                                 double max_patch_lifetime = 35.10667,
                                 bool geometric = true) {
  using Strat = K93_Strategy;
  using Env   = K93_Environment;

  Strat strat;
  strat.pars.b_0 = b_0;
  strat.is_variable_birth_rate = false;
  strat.birth_rate_y = {birth_rate};

  Parameters<Strat, Env> p;
  p.strategies.push_back(strat);
  p.max_patch_lifetime = max_patch_lifetime;
  p.validate();

  Env env;
  Control ctrl;
  ctrl.node_geometric_compression = geometric;
  ctrl.save_RK45_cache = true;

  SCM<Strat, Env> scm(p, env, ctrl);
  scm.run();

  std::vector<double> offspring = scm.offspring_production();
  double census = scm.get_system_ref().compute_competition(0.0);

  return Rcpp::List::create(
      Rcpp::Named("offspring") = Rcpp::wrap(offspring),
      Rcpp::Named("census")    = census);
}

// Layers (c)/(d) of the K93 SCM census gradient (CD-G): the exact reverse-mode
// trait gradient of an emergent SCM census metric.
//
// The SCM itself is the runnable the gradient is taken over -- it already
// duck-types the odelia gradient-driver contract (value_type, get_system_ref,
// ad_parameters, ad_initial_state, reset, run, an assignable `tape`). No
// wrapper: the driver seeds the strategy params on the patch, reset()s, run()s,
// and reduces the positioned patch through a census functional.
//
// The active scalar's adaptive stepper is compiled out (OdeControl is double
// only), so the active SCM replays the L1 ode-time schedule recorded on a double
// run. Crucially this is the RESIDENT path, not run_mutant: we pin the recorded
// step times through the public node-schedule surface but do NOT set_mutant, so
// the environment field is recomputed at the active scalar (L2) rather than
// frozen (L3). That recompute is what makes the gradient exact.
//
// Correctness is checked by the JVP=VJP dot-product oracle <J v, u> = <v, Jᵀ u>:
// for the scalar census functional the reverse gradient g is Jᵀ·1, the forward
// directional derivative along v is J v = g·v, so compute_directional_derivative
// must equal dot(g, v) to machine precision -- an FD-free correctness gate.
//
// Differentiate the census metric w.r.t. a representative set of K93 traits:
// b_0, b_1 (growth) and d_0 (mortality). Field order is K93_AD_FIELDS:
// {height_0, b_0, b_1, b_2, c_0, c_1, d_0, d_1, S_D, eta, k_I}.
static const std::vector<int> K93_CENSUS_TARGET_IDX = {1, 2, 6};

// [[Rcpp::export]]
Rcpp::List k93_scm_census_gradient(double b_0 = 0.059, double birth_rate = 20.0,
                                   double max_patch_lifetime = 35.10667) {
  using DStrat  = K93_Strategy_<double>;
  using DEnv    = K93_Environment_<double>;
  using RevS    = xad::adj<double>::active_type;  // AReal: compute_gradient
  using FwdS    = xad::fwd<double>::active_type;  // FReal: compute_jvp

  // Build validated K93 parameters for a strategy scalar S from the shared double
  // inputs. birth_rate_y is a plain double vector (the birth-rate schedule is not
  // itself a differentiated trait), so it takes the double directly.
  auto make_params = [&](auto strat_tag) {
    using Strat = std::decay_t<decltype(strat_tag)>;
    using Env   = typename Strat::environment_type;
    Strat strat;
    strat.pars.b_0 = b_0;
    strat.is_variable_birth_rate = false;
    strat.birth_rate_y = {birth_rate};
    Parameters<Strat, Env> p;
    p.strategies.push_back(strat);
    p.max_patch_lifetime = max_patch_lifetime;
    p.validate();
    return p;
  };

  // Recording Control (double run): save_RK45_cache turns on step recording so
  // the accepted ode-step times land in patch.step_history. Replay Control
  // (active runs): recording OFF, so the active SCM recomputes its environment
  // field at the active scalar (resident L2) rather than freezing a recorded one.
  Control ctrl_record;
  ctrl_record.node_geometric_compression = true;
  ctrl_record.save_RK45_cache = true;
  Control ctrl_replay;
  ctrl_replay.node_geometric_compression = true;

  // (1) Double run: record the L1 ode-time schedule the adaptive stepper used.
  std::vector<double> step_history;
  double census_double = 0.0;
  {
    DEnv env;
    SCM<DStrat, DEnv> scm(make_params(DStrat()), env, ctrl_record);
    scm.run();
    step_history  = scm.get_system_ref().step_history;
    census_double = scm.get_system_ref().compute_competition(0.0);
  }

  // Configure an active SCM for the RESIDENT replay: pin the recorded ode times
  // through the public node-schedule surface (reset first so the guard passes),
  // leaving the environment free to recompute at the active scalar -- resident
  // L2, not a frozen mutant (run_mutant). The census functional reads the total
  // stand competition (a basal-area moment) straight off the positioned patch.
  auto pin_replay = [&](auto& scm) {
    scm.reset();
    NodeSchedule ns = scm.r_node_schedule();
    ns.r_set_ode_times(step_history);
    ns.r_set_use_ode_times(true);
    scm.r_set_node_schedule(ns);
  };
  auto targets_for = [&](auto& scm) {
    odelia::ode::DifferentiationTargets t;
    t.params = K93_CENSUS_TARGET_IDX;
    auto ptrs = scm.get_system_ref().ad_parameters();
    for (int i : K93_CENSUS_TARGET_IDX) t.values.push_back(xad::value(*ptrs[i]));
    return t;
  };

  // Reverse leg: value + gradient g = Jᵀ·1 over an AReal SCM.
  double value = 0.0;
  std::vector<double> g;
  {
    K93_Environment_<RevS> env;
    SCM<K93_Strategy_<RevS>, K93_Environment_<RevS>> scm(
        make_params(K93_Strategy_<RevS>()), env, ctrl_replay);
    pin_replay(scm);
    auto functional = [](decltype(scm)& s) -> RevS {
      return s.get_system_ref().compute_competition(0.0);
    };
    auto res = odelia::ode::compute_gradient(scm, targets_for(scm), functional);
    value = res.first;
    g     = res.second;
  }

  // Forward leg: directional derivative J v along a deterministic v over an FReal
  // SCM. The JVP=VJP oracle wants it to equal dot(g, v) to machine precision.
  std::vector<double> v(K93_CENSUS_TARGET_IDX.size());
  for (std::size_t i = 0; i < v.size(); ++i) v[i] = 0.3 + 0.15 * double(i);
  double jvp = 0.0;
  {
    K93_Environment_<FwdS> env;
    SCM<K93_Strategy_<FwdS>, K93_Environment_<FwdS>> scm(
        make_params(K93_Strategy_<FwdS>()), env, ctrl_replay);
    pin_replay(scm);
    auto functional = [](decltype(scm)& s) -> FwdS {
      return s.get_system_ref().compute_competition(0.0);
    };
    auto dd = odelia::ode::compute_directional_derivative(scm, targets_for(scm), v,
                                                          functional);
    jvp = dd.second;
  }

  double dot = 0.0;
  for (std::size_t i = 0; i < g.size(); ++i) dot += g[i] * v[i];

  return Rcpp::List::create(
      Rcpp::Named("value")         = value,
      Rcpp::Named("census_double") = census_double,
      Rcpp::Named("grad")          = Rcpp::wrap(g),
      Rcpp::Named("jvp")           = jvp,
      Rcpp::Named("dot")           = dot,
      Rcpp::Named("n_steps")       = double(step_history.size()));
}
