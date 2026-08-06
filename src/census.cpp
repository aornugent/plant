#include <plant.h>
#include <plant/census.h>

// R entry points for the stand census and its direct (fixed-state) sensitivity
// term. See inst/include/plant/census.h for the algebra and for why this is one
// of two terms.
//
// The state arrives as plain vectors, never as a live patch. That is deliberate:
// it is what makes "at fixed state" true by construction rather than by
// convention. Nothing prepare_strategy() or a parameter assignment does can
// reach back and move a node's height, density or quadrature weight, because the
// grid has already left the object.

namespace plant {
namespace {

template <typename S>
Rcpp::List census_direct_impl(S s,
                              const Rcpp::NumericVector& height,
                              const Rcpp::NumericVector& density,
                              const Rcpp::NumericVector& area_heartwood,
                              const Rcpp::NumericVector& mass_heartwood) {
  // Resolve eta_c and friends: a strategy arriving from R is unprepared, and
  // eta_c is what d/d(eta) is a derivative of.
  s.prepare_strategy();

  const R_xlen_t n = height.size();
  if (density.size() != n || area_heartwood.size() != n ||
      mass_heartwood.size() != n) {
    Rcpp::stop("height, density, area_heartwood and mass_heartwood must have "
               "equal length");
  }

  const std::vector<double> h(height.begin(), height.end());
  const std::vector<double> d(density.begin(), density.end());
  const std::vector<double> ah(area_heartwood.begin(), area_heartwood.end());
  const std::vector<double> mh(mass_heartwood.begin(), mass_heartwood.end());

  const std::size_t M = census::n_metrics(), P = census::n_parameters();
  std::vector<double> total(M, 0.0), grad(M * P, 0.0);
  census::accumulate(s, h, d, ah, mh, total.data(), grad.data());

  const std::vector<std::string> mn = census::metric_names();
  const std::vector<std::string> pn = census::parameter_names();

  Rcpp::NumericVector census(M);
  for (std::size_t i = 0; i < M; ++i) {
    census[i] = total[i];
  }
  census.attr("names") = Rcpp::wrap(mn);

  Rcpp::NumericMatrix gradient(static_cast<int>(M), static_cast<int>(P));
  Rcpp::LogicalMatrix support(static_cast<int>(M), static_cast<int>(P));
  for (std::size_t i = 0; i < M; ++i) {
    for (std::size_t j = 0; j < P; ++j) {
      gradient(i, j) = grad[i * P + j];
      support(i, j) = census::reaches(i, j);
    }
  }
  Rcpp::List dn = Rcpp::List::create(Rcpp::wrap(mn), Rcpp::wrap(pn));
  gradient.attr("dimnames") = dn;
  support.attr("dimnames") = dn;

  return Rcpp::List::create(
    Rcpp::_["census"] = census,
    Rcpp::_["gradient"] = gradient,
    Rcpp::_["support"] = support,
    Rcpp::_["grid_is_monotone"] = census::grid_is_monotone(h));
}

}
}

// [[Rcpp::export]]
std::vector<std::string> census_metric_names() {
  return plant::census::metric_names();
}

// [[Rcpp::export]]
std::vector<std::string> census_parameter_names() {
  return plant::census::parameter_names();
}

// [[Rcpp::export]]
Rcpp::List FF16_census_direct(plant::FF16_Strategy s,
                              Rcpp::NumericVector height,
                              Rcpp::NumericVector density,
                              Rcpp::NumericVector area_heartwood,
                              Rcpp::NumericVector mass_heartwood) {
  return plant::census_direct_impl(s, height, density, area_heartwood,
                                   mass_heartwood);
}

// [[Rcpp::export]]
Rcpp::List TF24_census_direct(plant::TF24_Strategy s,
                              Rcpp::NumericVector height,
                              Rcpp::NumericVector density,
                              Rcpp::NumericVector area_heartwood,
                              Rcpp::NumericVector mass_heartwood) {
  return plant::census_direct_impl(s, height, density, area_heartwood,
                                   mass_heartwood);
}

// [[Rcpp::export]]
Rcpp::List TF24f_census_direct(plant::TF24f_Strategy s,
                               Rcpp::NumericVector height,
                               Rcpp::NumericVector density,
                               Rcpp::NumericVector area_heartwood,
                               Rcpp::NumericVector mass_heartwood) {
  return plant::census_direct_impl(s, height, density, area_heartwood,
                                   mass_heartwood);
}
