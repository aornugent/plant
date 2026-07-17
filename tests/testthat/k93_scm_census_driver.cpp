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
