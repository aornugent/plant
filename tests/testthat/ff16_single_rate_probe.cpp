// [[Rcpp::plugins(cpp20)]]
// The simplest possible check: a SINGLE FF16 individual in a FIXED light
// environment (no field, no feedback, no trajectory), compute its growth (height)
// rate, and compare d(growth)/d(lma) reverse-AD vs central FD. If this is off,
// the FF16 detached edge is a plain bug in the single-plant rate path's lma
// derivative -- not feedback, not schedule, not the field.
#include <Rcpp.h>
#include <plant/models/ff16_strategy.h>
#include <plant/models/ff16_environment.h>
#include <plant/individual.h>
#include <odelia/gradient.hpp>
#include <vector>

using namespace plant;

// Growth (height) rate of one FF16 individual of the given (double) height in a
// fixed-light environment, at strategy scalar S with lma active.
template <class S>
static S single_growth_rate(S lma, double height, double light) {
  FF16_Strategy_<S> strat;
  strat.pars.lma = lma;
  auto s = std::make_shared<FF16_Strategy_<S>>(strat);
  FF16_Environment_<S> env;
  env.set_fixed_environment(light);          // constant light: no field, no feedback
  Individual<FF16_Strategy_<S>, FF16_Environment_<S>> ind(s);
  ind.prepare_strategy();
  ind.set_state("height", S(height));
  ind.compute_rates(env);
  return ind.rate(0);                         // HEIGHT_INDEX rate
}

// [[Rcpp::export]]
Rcpp::List ff16_single_rate_probe(double lma = 0.1978791, double height = 5.0,
                                  double light = 0.5, double fd_rel = 1e-4) {
  using RevS = xad::adj<double>::active_type;
  double grad_ad = 0.0, value = 0.0;
  {
    xad::Tape<double> tape;
    RevS L = lma; tape.registerInput(L); tape.newRecording();
    RevS g = single_growth_rate<RevS>(L, height, light);
    tape.registerOutput(g); xad::derivative(g) = 1.0; tape.computeAdjoints();
    value = xad::value(g); grad_ad = xad::derivative(L);
  }
  double d = fd_rel * (lma + 1e-3);
  double fd = (single_growth_rate<double>(lma + d, height, light) -
               single_growth_rate<double>(lma - d, height, light)) / (2.0 * d);
  return Rcpp::List::create(Rcpp::Named("value") = value,
                            Rcpp::Named("grad_ad") = grad_ad,
                            Rcpp::Named("fd") = fd);
}
