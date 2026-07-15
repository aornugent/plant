// [[Rcpp::plugins(cpp20)]]
// Verify the resident light->leaf channel (Stage E, light channel): seed the
// resident light (canopy openness) active, compute the individual growth rate
// through the leaf (which reads the active light via the seam's light channel),
// and FD-check d(growth)/d(light) at one compute_rates. Tests both single-solve
// shading models.
#include <Rcpp.h>
#include <plant/models/tf24_strategy.h>
#include <plant/individual.h>
#include <XAD/XAD.hpp>

using namespace plant;

static double growth_double(double light, double height, const std::string& sm) {
  TF24_Strategy_<double> s;
  s.control.shading_model = sm;
  auto sp = make_strategy_ptr(s);
  Individual<TF24_Strategy_<double>, TF24_Environment_<double>> ind(sp);
  TF24_Environment_<double> env;
  env.set_fixed_environment(light);
  ind.set_state(HEIGHT_INDEX, height);
  ind.compute_rates(env);
  return ind.rate(HEIGHT_INDEX);
}

// [[Rcpp::export]]
Rcpp::List tf24_light_channel_check(std::string shading_model = "",
                                    double light = 0.6, double height = 1.0,
                                    double delta = 1e-6) {
  using ActiveS = xad::adj<double>::active_type;
  using Tape    = xad::adj<double>::tape_type;
  using AStrat  = TF24_Strategy_<ActiveS>;
  using AEnv    = TF24_Environment_<ActiveS>;

  Tape tape;  // ctor auto-activates

  AStrat s;
  s.control.shading_model = shading_model;
  auto sp = make_strategy_ptr(s);
  Individual<AStrat, AEnv> ind(sp);
  AEnv env;

  ActiveS L = light;
  tape.registerInput(L);
  tape.newRecording();

  env.set_fixed_environment_scalar(L, 150.0);
  ind.set_state(HEIGHT_INDEX, ActiveS(height));
  ind.compute_rates(env);

  ActiveS growth = ind.rate(HEIGHT_INDEX);
  tape.registerOutput(growth);
  xad::derivative(growth) = 1.0;
  tape.computeAdjoints();
  const double ad_grad = xad::derivative(L);
  const double value = xad::value(growth);
  tape.deactivate();

  const double fd_grad = (growth_double(light + delta, height, shading_model) -
                          growth_double(light - delta, height, shading_model)) / (2.0 * delta);

  return Rcpp::List::create(
      Rcpp::Named("value")   = value,
      Rcpp::Named("ad_grad") = ad_grad,
      Rcpp::Named("fd_grad") = fd_grad,
      Rcpp::Named("abs_err") = std::abs(ad_grad - fd_grad));
}
