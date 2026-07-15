// [[Rcpp::plugins(cpp20)]]
#include <Rcpp.h>
#include <plant/models/tf24_strategy.h>
#include <plant/individual.h>
#include <plant/leaf_model.h>
#include <string>
#include <vector>
using namespace plant;

// Directly verify Leaf::dsoil_consumption_dpsi_collar_perlayer (the #47 analytic
// per-layer d(soil_consumption_)/d(P_x_r)) against a pure-double FD of
// E_from_Soil_to_Root_Collar w.r.t. the collar potential P_x_r. No tape, no seam,
// no census metric -- isolates the analytic derivative's correctness (math + sign).
// [[Rcpp::export]]
Rcpp::List tf24_collar_deriv_check(double height = 5.0, double delta = 1e-6) {
  TF24_Strategy_<double> s;
  auto sp = make_strategy_ptr(s);
  Individual<TF24_Strategy_<double>, TF24_Environment_<double>> ind(sp);
  ind.set_state("height", height);
  TF24_Environment_<double> env; env.set_fixed_environment(1.0);
  ind.compute_rates(env);            // sets up sp->leaf at its operating point

  Leaf& L = sp->leaf;
  const double P_x_r = L.root_collar_psi_;
  const std::vector<double> psi = L.psi_soil_inverted_;

  std::vector<double> analytic;
  L.dsoil_consumption_dpsi_collar_perlayer(P_x_r, psi, analytic);

  // central FD of E_from_Soil_to_Root_Collar's per-layer soil_consumption_.
  L.E_from_Soil_to_Root_Collar(P_x_r + delta, psi);
  const std::vector<double> plus = L.soil_consumption_;
  L.E_from_Soil_to_Root_Collar(P_x_r - delta, psi);
  const std::vector<double> minus = L.soil_consumption_;
  L.E_from_Soil_to_Root_Collar(P_x_r, psi);  // restore

  const int n = (int)analytic.size();
  Rcpp::NumericVector an(n), fd(n), rel(n);
  for (int i = 0; i < n; ++i) {
    an[i] = analytic[i];
    fd[i] = (plus[i] - minus[i]) / (2.0 * delta);
    rel[i] = std::abs(an[i] - fd[i]) / (std::abs(fd[i]) + 1e-30);
  }
  return Rcpp::List::create(
      Rcpp::Named("P_x_r") = P_x_r, Rcpp::Named("analytic") = an,
      Rcpp::Named("fd") = fd, Rcpp::Named("rel") = rel);
}
