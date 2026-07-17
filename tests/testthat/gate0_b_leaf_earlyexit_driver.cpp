// [[Rcpp::plugins(cpp20)]]
// Gate-0 check B (#4): are the TF24 leaf shut-down early-exits C0 (continuous)
// in the driving variable, so they classify as decide() predicates rather than
// breakpoint kinks? Sweep soil moisture theta from wet to the residual bound
// (which drives psi_soil past psi_crit -> the shut-down early-exits fire) and
// record the leaf profit at each theta. A decide() classification requires
// profit(theta) to have NO jump across the transition; a value jump would force
// a breakpoint node instead.
#include <Rcpp.h>
#include <plant/models/tf24_strategy.h>
#include <XAD/XAD.hpp>

using namespace plant;

// [[Rcpp::export]]
Rcpp::List tf24_leaf_earlyexit_sweep(double height = 5.0,
                                     double theta_lo = 0.011,
                                     double theta_hi = 0.214,
                                     int n = 400) {
  TF24_Strategy_<double> s;
  s.prepare_strategy();

  TF24_Environment_<double> env;   // 5 soil layers by default
  env.set_fixed_environment(1.0);  // full light openness

  const int nsoil = env.get_soil_number_of_depths();
  const double area_leaf_ = s.area_leaf(height);
  const double height_inverse = 1.0 / height;

  Rcpp::NumericVector theta(n), profit(n), net(n), psi_stem(n), collar(n);
  for (int i = 0; i < n; ++i) {
    const double th = theta_lo + (theta_hi - theta_lo) * i / (n - 1);
    // Set every soil layer to the same moisture so the whole column crosses the
    // shut-down thresholds together (the sharpest test of continuity).
    for (int L = 0; L < nsoil; ++L) env.vars.set_state(L, th);
    env.psi_soil_cache_valid_ = false;

    const double np = s.net_mass_production_dt(env, height, area_leaf_, height_inverse);
    theta[i]    = th;
    profit[i]   = s.leaf.profit_;
    net[i]      = np;
    psi_stem[i] = s.leaf.opt_psi_stem_;
    collar[i]   = s.leaf.root_collar_psi_;
  }
  return Rcpp::List::create(
      Rcpp::Named("theta") = theta, Rcpp::Named("profit") = profit,
      Rcpp::Named("net") = net, Rcpp::Named("psi_stem") = psi_stem,
      Rcpp::Named("collar") = collar);
}
