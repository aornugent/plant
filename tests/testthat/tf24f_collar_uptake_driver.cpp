// [[Rcpp::plugins(cpp20)]]
// Verify the TF24f collar-psi channel (Stage D): the leaf runs at a tracked
// collar-psi ODE state (off the optimum), so d(profit)/d(psi) != 0 and the seam
// must inject it onto the tracked state. Seed the tracked state active, compute
// the growth rate through the leaf, and FD-check d(growth)/d(tracked_psi).
#include <Rcpp.h>
#include <plant/models/tf24f_strategy.h>
#include <plant/individual.h>
#include <XAD/XAD.hpp>

using namespace plant;

static double uptake_double(double tracked_psi, double height, int layer) {
  TF24f_Strategy_<double> s;
  auto sp = make_strategy_ptr(s);
  Individual<TF24f_Strategy_<double>, TF24_Environment_<double>> ind(sp);
  TF24_Environment_<double> env;
  env.set_fixed_environment(1.0);
  ind.set_state(HEIGHT_INDEX, height);
  // tracked collar-psi state is the 6th state (index 5)
  const int idx = static_cast<int>(TF24_Strategy_<double>::state_size());
  ind.set_state(idx, tracked_psi);
  ind.compute_rates(env);
  return ind.consumption_rate(layer);
}

// [[Rcpp::export]]
Rcpp::List tf24f_collar_uptake_check(double tracked_psi = 1.5, double height = 1.0, int layer = 0, double delta = 1e-6) {
  using ActiveS = xad::adj<double>::active_type;
  using Tape    = xad::adj<double>::tape_type;
  using AStrat  = TF24f_Strategy_<ActiveS>;
  using AEnv    = TF24_Environment_<ActiveS>;

  Tape tape;  // ctor auto-activates

  AStrat s;
  auto sp = make_strategy_ptr(s);
  Individual<AStrat, AEnv> ind(sp);
  AEnv env;
  env.set_fixed_environment(1.0);

  ActiveS psi = tracked_psi;
  tape.registerInput(psi);
  tape.newRecording();

  ind.set_state(HEIGHT_INDEX, ActiveS(height));
  // tracked collar-psi state is appended after TF24's 5 states (index 5)
  const int idx = static_cast<int>(TF24_Strategy_<ActiveS>::state_size());
  ind.set_state(idx, psi);
  ind.compute_rates(env);

  ActiveS out = ind.consumption_rate(layer);
  tape.registerOutput(out);
  xad::derivative(out) = 1.0;
  tape.computeAdjoints();
  const double ad_grad = xad::derivative(psi);
  const double value = xad::value(out);
  tape.deactivate();

  const double fd_grad = (uptake_double(tracked_psi + delta, height, layer) -
                          uptake_double(tracked_psi - delta, height, layer)) / (2.0 * delta);

  return Rcpp::List::create(
      Rcpp::Named("value")   = value,
      Rcpp::Named("ad_grad") = ad_grad,
      Rcpp::Named("fd_grad") = fd_grad,
      Rcpp::Named("abs_err") = std::abs(ad_grad - fd_grad));
}
