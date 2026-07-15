// [[Rcpp::plugins(cpp20)]]
// Verify the resident soil->leaf channel of Stage C: seed a soil-moisture layer
// active, compute the individual's growth rate through the TF24 leaf (which reads
// the active soil psi via the supplied_derivative seam), and FD-check
// d(growth_rate)/d(soil_moisture) at one compute_rates call.
#include <Rcpp.h>
#include <plant/models/tf24_strategy.h>
#include <plant/individual.h>
#include <XAD/XAD.hpp>

using namespace plant;

static double growth_double(double theta0, double height) {
  TF24_Strategy_<double> s;
  auto sp = make_strategy_ptr(s);
  Individual<TF24_Strategy_<double>, TF24_Environment_<double>> ind(sp);
  TF24_Environment_<double> env;
  env.set_fixed_environment(1.0);
  env.vars.set_state(0, theta0);
  env.psi_soil_cache_valid_ = false;
  ind.set_state(HEIGHT_INDEX, height);
  ind.compute_rates(env);
  return ind.rate(HEIGHT_INDEX);
}

// [[Rcpp::export]]
Rcpp::List tf24_soil_channel_check(double height = 1.0, double delta = 1e-6) {
  using ActiveS = xad::adj<double>::active_type;
  using Tape    = xad::adj<double>::tape_type;
  using AStrat  = TF24_Strategy_<ActiveS>;
  using AEnv    = TF24_Environment_<ActiveS>;

  const double theta0 = 0.428 * 0.5;  // default initial soil moisture (layer 0)

  Tape tape;  // ctor auto-activates

  AStrat s;
  auto sp = make_strategy_ptr(s);
  Individual<AStrat, AEnv> ind(sp);
  AEnv env;
  env.set_fixed_environment(1.0);

  ActiveS theta = theta0;
  tape.registerInput(theta);
  tape.newRecording();

  env.vars.set_state(0, theta);
  env.psi_soil_cache_valid_ = false;
  ind.set_state(HEIGHT_INDEX, ActiveS(height));
  ind.compute_rates(env);

  ActiveS growth = ind.rate(HEIGHT_INDEX);
  tape.registerOutput(growth);
  xad::derivative(growth) = 1.0;
  tape.computeAdjoints();
  const double ad_grad = xad::derivative(theta);
  const double value = xad::value(growth);
  tape.deactivate();

  const double fd_grad = (growth_double(theta0 + delta, height) -
                          growth_double(theta0 - delta, height)) / (2.0 * delta);

  return Rcpp::List::create(
      Rcpp::Named("value")   = value,
      Rcpp::Named("ad_grad") = ad_grad,
      Rcpp::Named("fd_grad") = fd_grad,
      Rcpp::Named("abs_err") = std::abs(ad_grad - fd_grad));
}
