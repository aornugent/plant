// [[Rcpp::plugins(cpp20)]]
// Full-SCM reverse-mode trait gradient for TF24 and TF24f, through the SAME
// run-shaped entry (plant/scm_gradient.h) FF16/K93 use: the caller passes trait
// values + target field names + a metric; the entry owns the adaptive double
// refine, the double->active rebind (SCM::rebind_from, carrying TF24's soil
// config), and the resolved-schedule replay. This exercises the growing-dimension
// SCM gradient path a user actually calls for R0 (offspring) and census, which the
// bespoke leaf-coupling drivers do NOT. FD is the entry's own adaptive
// reoptimising central difference (the R test perturbs the trait and re-refines).
#include <Rcpp.h>
#include <plant/models/tf24_strategy.h>
#include <plant/models/tf24_environment.h>
#include <plant/models/tf24f_strategy.h>
#include <plant/individual.h>
#include <plant/parameters.h>
#include <plant/control.h>
#include <plant/scm.h>
#include <plant/scm_gradient.h>
#include <odelia/gradient.hpp>
#include <string>
#include <vector>

using namespace plant;

// Apply optional {name -> value} trait overrides to a strategy, so the R side can
// perturb a trait for a finite-difference reference (the entry re-refines adaptively
// each call). Empty override_names is a no-op (defaults).
template <class StratD>
static void apply_overrides(StratD& s, const std::vector<std::string>& override_names,
                            const std::vector<double>& override_values) {
  if (override_names.empty()) return;
  auto names = StratD::field_names();
  auto ptrs  = s.field_ptrs();
  for (std::size_t k = 0; k < override_names.size(); ++k) {
    int idx = -1;
    for (std::size_t i = 0; i < names.size(); ++i)
      if (names[i] == override_names[k]) idx = (int)i;
    if (idx < 0) Rcpp::stop("unknown TF24 override field: " + override_names[k]);
    *ptrs[idx] = override_values[k];
  }
}

// Resolve trait names to field indices on the double strategy, reading the seed
// values off a configured strategy (birth_rate held fixed, not differentiated).
template <template <class> class StratT>
static std::pair<double, std::vector<double>>
tf24_entry(std::vector<std::string> targets, double birth_rate,
           double max_patch_lifetime, int metric,
           std::vector<std::string> override_names = {},
           std::vector<double> override_values = {}) {
  using StratD = StratT<double>;
  using EnvD   = typename StratD::environment_type;

  StratD s;
  apply_overrides(s, override_names, override_values);
  s.is_variable_birth_rate = false;
  s.birth_rate_y = {birth_rate};
  Parameters<StratD, EnvD> p;
  p.strategies.push_back(s);
  p.max_patch_lifetime = max_patch_lifetime;
  p.validate();

  Control ctrl;

  auto names = StratD::field_names();
  odelia::ode::DifferentiationTargets t;
  {
    StratD s0;
    apply_overrides(s0, override_names, override_values);
    auto ptrs = s0.field_ptrs();
    for (auto const& nm : targets) {
      int idx = -1;
      for (std::size_t i = 0; i < names.size(); ++i) if (names[i] == nm) idx = (int)i;
      if (idx < 0) Rcpp::stop("unknown TF24 field: " + nm);
      t.params.push_back(idx);
      t.values.push_back(*ptrs[idx]);
    }
  }

  return (metric == 1) ? scm_gradient(p, ctrl, t, census_metric{})
                       : scm_gradient(p, ctrl, t, offspring_metric{});
}

template <template <class> class StratT>
static Rcpp::List tf24_census_vec(std::vector<std::string> targets, double birth_rate,
                                  double max_patch_lifetime,
                                  std::vector<std::string> override_names = {},
                                  std::vector<double> override_values = {}) {
  using StratD = StratT<double>;
  using EnvD   = typename StratD::environment_type;

  StratD s;
  apply_overrides(s, override_names, override_values);
  s.is_variable_birth_rate = false;
  s.birth_rate_y = {birth_rate};
  Parameters<StratD, EnvD> p;
  p.strategies.push_back(s);
  p.max_patch_lifetime = max_patch_lifetime;
  p.validate();

  Control ctrl;
  auto names = StratD::field_names();
  odelia::ode::DifferentiationTargets t;
  {
    StratD s0;
    apply_overrides(s0, override_names, override_values);
    auto ptrs = s0.field_ptrs();
    for (auto const& nm : targets) {
      int idx = -1;
      for (std::size_t i = 0; i < names.size(); ++i) if (names[i] == nm) idx = (int)i;
      if (idx < 0) Rcpp::stop("unknown TF24 field: " + nm);
      t.params.push_back(idx);
      t.values.push_back(*ptrs[idx]);
    }
  }

  auto out = scm_jacobian(p, ctrl, t, census_vector{});
  Rcpp::NumericMatrix J(out.second.size(),
                        out.second.empty() ? 0 : out.second[0].size());
  for (std::size_t r = 0; r < out.second.size(); ++r)
    for (std::size_t c = 0; c < out.second[r].size(); ++c) J(r, c) = out.second[r][c];
  return Rcpp::List::create(Rcpp::Named("values")   = Rcpp::wrap(out.first),
                            Rcpp::Named("jacobian") = J);
}

// [[Rcpp::export]]
Rcpp::List tf24_scm_gradient(std::vector<std::string> targets, double birth_rate,
                             double max_patch_lifetime, int metric,
                             std::vector<std::string> override_names,
                             std::vector<double> override_values) {
  auto r = tf24_entry<TF24_Strategy_>(targets, birth_rate, max_patch_lifetime, metric,
                                      override_names, override_values);
  TF24_Strategy_<double> s0; apply_overrides(s0, override_names, override_values);
  auto nm = TF24_Strategy_<double>::field_names(); auto pt = s0.field_ptrs();
  std::vector<double> seed;
  for (auto const& t : targets) for (std::size_t i=0;i<nm.size();++i) if(nm[i]==t) seed.push_back(*pt[i]);
  return Rcpp::List::create(Rcpp::Named("value") = r.first,
                            Rcpp::Named("grad")  = Rcpp::wrap(r.second),
                            Rcpp::Named("seed")  = Rcpp::wrap(seed));
}

// [[Rcpp::export]]
Rcpp::List tf24f_scm_gradient(std::vector<std::string> targets, double birth_rate,
                              double max_patch_lifetime, int metric,
                              std::vector<std::string> override_names,
                              std::vector<double> override_values) {
  auto r = tf24_entry<TF24f_Strategy_>(targets, birth_rate, max_patch_lifetime, metric,
                                       override_names, override_values);
  return Rcpp::List::create(Rcpp::Named("value") = r.first,
                            Rcpp::Named("grad")  = Rcpp::wrap(r.second));
}

// [[Rcpp::export]]
Rcpp::List tf24_scm_census_vector(std::vector<std::string> targets, double birth_rate,
                                  double max_patch_lifetime,
                                  std::vector<std::string> override_names,
                                  std::vector<double> override_values) {
  return tf24_census_vec<TF24_Strategy_>(targets, birth_rate, max_patch_lifetime,
                                         override_names, override_values);
}

// [[Rcpp::export]]
Rcpp::List tf24f_scm_census_vector(std::vector<std::string> targets, double birth_rate,
                                   double max_patch_lifetime,
                                   std::vector<std::string> override_names,
                                   std::vector<double> override_values) {
  return tf24_census_vec<TF24f_Strategy_>(targets, birth_rate, max_patch_lifetime,
                                          override_names, override_values);
}
