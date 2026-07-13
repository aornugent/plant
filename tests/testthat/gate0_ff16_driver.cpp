// [[Rcpp::plugins(cpp20)]]
#include <Rcpp.h>
#include <plant/models/ff16_strategy.h>
#include <plant/individual.h>
#include <plant/individual_runner.h>
#include <odelia/ode_solver.hpp>
#include <odelia/gradient.hpp>
#include <string>
#include <vector>

using namespace plant;

// Set one named low-level FF16 parameter on a strategy of any scalar (value in,
// narrowed from double).
template <class S>
static void set_param(FF16_Strategy_<S>& s, const std::string& name, double v) {
  auto names = FF16_Strategy_<S>::field_names();
  auto ptrs  = s.field_ptrs();
  for (std::size_t i = 0; i < names.size(); ++i)
    if (names[i] == name) { *ptrs[i] = S(v); return; }
  Rcpp::stop("unknown FF16 parameter: " + name);
}

// Final plant height after replaying `sched` at S=double with the parameter set.
static double metric_double(const std::string& name, double v,
                            const std::vector<double>& sched) {
  FF16_Strategy_<double> s;
  set_param(s, name, v);
  auto sp = make_strategy_ptr(s);
  Individual<FF16_Strategy_<double>, FF16_Environment_<double>> ind(sp);
  FF16_Environment_<double> env;
  tools::IndividualRunner<FF16_Strategy_<double>, FF16_Environment_<double>> runner(ind, env);
  odelia::ode::Solver<decltype(runner)> solver(runner, odelia::ode::OdeControl());
  solver.set_schedule(sched);
  solver.reset();
  solver.run();
  std::vector<double> st(solver.get_system_ref().ode_size());
  solver.get_system_ref().ode_state(st.begin());
  return st[0];  // height
}

// [[Rcpp::export]]
Rcpp::List ff16_gate0_fd_check(std::string param, double t_end = 5.0,
                               double delta = 1e-5) {
  using ActiveS = xad::adj<double>::active_type;
  using DStrat = FF16_Strategy_<double>;
  using DEnv   = FF16_Environment_<double>;
  using AStrat = FF16_Strategy_<ActiveS>;
  using AEnv   = FF16_Environment_<ActiveS>;

  // Baseline value of the chosen parameter.
  DStrat base_s;
  double p0 = 0.0;
  {
    auto names = DStrat::field_names();
    auto ptrs  = base_s.field_ptrs();
    bool found = false;
    for (std::size_t i = 0; i < names.size(); ++i)
      if (names[i] == param) { p0 = xad::value(*ptrs[i]); found = true; }
    if (!found) Rcpp::stop("unknown FF16 parameter: " + param);
  }

  // Record the fixed RK schedule on a double adaptive pass [0, t_end].
  std::vector<double> sched;
  {
    auto sp = make_strategy_ptr(DStrat());
    Individual<DStrat, DEnv> ind(sp);
    DEnv env;
    tools::IndividualRunner<DStrat, DEnv> runner(ind, env);
    odelia::ode::Solver<decltype(runner)> solver(runner, odelia::ode::OdeControl());
    solver.advance_adaptive({0.0, t_end});
    sched = solver.times();
  }

  // AD gradient: active runner, replay the recorded schedule, d(height)/d(param).
  double value = 0.0, ad_grad = 0.0;
  {
    auto sp = make_strategy_ptr(AStrat());
    Individual<AStrat, AEnv> ind(sp);
    AEnv env;
    tools::IndividualRunner<AStrat, AEnv> runner(ind, env);
    odelia::ode::Solver<decltype(runner)> solver(runner, odelia::ode::OdeControl());
    solver.set_schedule(sched);

    // Seed the chosen parameter by its index in field order.
    auto names = AStrat::field_names();
    int idx = -1;
    for (std::size_t i = 0; i < names.size(); ++i) if (names[i] == param) idx = (int)i;
    odelia::ode::DifferentiationTargets targets;
    targets.params = {idx};
    targets.values = {p0};

    auto functional = [](odelia::ode::Solver<decltype(runner)>& sv) -> ActiveS {
      std::vector<ActiveS> st(sv.get_system_ref().ode_size());
      sv.get_system_ref().ode_state(st.begin());
      return st[0];  // final height
    };
    auto res = odelia::ode::compute_gradient(solver, targets, functional);
    value = res.first;
    ad_grad = res.second[0];
  }

  // Central finite difference of the same replayed metric.
  double fd_grad = (metric_double(param, p0 + delta, sched) -
                    metric_double(param, p0 - delta, sched)) / (2.0 * delta);

  return Rcpp::List::create(
      Rcpp::Named("param")   = param,
      Rcpp::Named("p0")      = p0,
      Rcpp::Named("value")   = value,
      Rcpp::Named("ad_grad") = ad_grad,
      Rcpp::Named("fd_grad") = fd_grad,
      Rcpp::Named("abs_err") = std::abs(ad_grad - fd_grad),
      Rcpp::Named("n_steps") = (double)sched.size());
}
