// [[Rcpp::plugins(cpp20)]]
// Exercises the run-shaped gradient entry (plant/scm_gradient.h). The caller passes
// ONLY trait values + target field indices + a metric -- NO schedule. The entry owns
// the adaptive double refine, the double->active rebind (SCM::rebind_from), and the
// resolved-schedule replay, so this driver cannot express a wrong replay grid. The R
// test gates the returned reverse gradient against a reoptimising (fully adaptive)
// run_scm FD -- the correctness reference.
#include <Rcpp.h>
#include <plant/models/ff16_strategy.h>
#include <plant/models/ff16_environment.h>
#include <plant/models/k93_strategy.h>
#include <plant/models/k93_environment.h>
#include <plant/individual.h>
#include <plant/parameters.h>
#include <plant/control.h>
#include <plant/scm.h>
#include <plant/scm_gradient.h>
#include <odelia/gradient.hpp>
#include <string>
#include <vector>

using namespace plant;

template <template <class> class StratT, class Setup>
static Rcpp::List run_entry(std::vector<int> target_idx, double birth_rate,
                            double max_patch_lifetime, int metric, Setup setup) {
  using StratD = StratT<double>;
  using EnvD   = typename StratD::environment_type;

  // Build a configured double Parameters -- the ordinary user object. No schedule.
  StratD s;
  setup(s);
  s.is_variable_birth_rate = false;
  s.birth_rate_y = {birth_rate};
  Parameters<StratD, EnvD> p;
  p.strategies.push_back(s);
  p.max_patch_lifetime = max_patch_lifetime;
  p.validate();

  Control ctrl;

  // Seed values, in field-index order, read off the configured strategy.
  odelia::ode::DifferentiationTargets t;
  t.params = target_idx;
  {
    StratD s0; setup(s0);
    auto ptrs = s0.field_ptrs();
    for (int i : target_idx) t.values.push_back(*ptrs[i]);
  }

  // scm_gradient returns odelia's {value, gradient}; the value-reproduces-double
  // R5 check is a structural assert inside the entry (it stops on a config-crossing
  // gap), so the driver need not surface a separate double value.
  std::pair<double, std::vector<double>> r =
      (metric == 1) ? scm_gradient(p, ctrl, t, census_metric{})
                    : scm_gradient(p, ctrl, t, offspring_metric{});
  return Rcpp::List::create(
      Rcpp::Named("value") = r.first,
      Rcpp::Named("grad")  = Rcpp::wrap(r.second),
      Rcpp::Named("names") = StratD::field_names());
}

// [[Rcpp::export]]
Rcpp::List ff16_entry_gradient(std::vector<int> target_idx, double lma = 0.1978791,
                               double birth_rate = 20.0,
                               double max_patch_lifetime = 40.0, int metric = 0) {
  auto setup = [&](FF16_Strategy_<double>& s) {
    auto p = s.field_ptrs(); auto nm = FF16_Strategy_<double>::field_names();
    for (std::size_t i = 0; i < nm.size(); ++i) if (nm[i] == "lma") *p[i] = lma;
  };
  return run_entry<FF16_Strategy_>(target_idx, birth_rate, max_patch_lifetime, metric, setup);
}

// The FF16 standing-stock census VECTOR gradient (codomain 3: LAI, above-ground
// biomass, basal area). Exercises scm_jacobian's multi-output path -- one adaptive
// recording, one active replay, three reverse sweeps -- and the mass-weighted
// Patch::census operator. Returns the 3 values and the 3 x n Jacobian.
// [[Rcpp::export]]
Rcpp::List ff16_census_vector_gradient(std::vector<int> target_idx,
                                       double lma = 0.1978791, double birth_rate = 20.0,
                                       double max_patch_lifetime = 40.0) {
  using StratD = FF16_Strategy_<double>;
  using EnvD   = StratD::environment_type;
  auto setup = [&](StratD& s) {
    auto p = s.field_ptrs(); auto nm = StratD::field_names();
    for (std::size_t i = 0; i < nm.size(); ++i) if (nm[i] == "lma") *p[i] = lma;
  };
  StratD s; setup(s);
  s.is_variable_birth_rate = false; s.birth_rate_y = {birth_rate};
  Parameters<StratD, EnvD> p;
  p.strategies.push_back(s);
  p.max_patch_lifetime = max_patch_lifetime; p.validate();
  Control ctrl;
  odelia::ode::DifferentiationTargets t; t.params = target_idx;
  { StratD s0; setup(s0); auto ptrs = s0.field_ptrs();
    for (int i : target_idx) t.values.push_back(*ptrs[i]); }

  auto out = scm_jacobian(p, ctrl, t, census_vector{});  // {values(3), jacobian(3 x n)}
  Rcpp::NumericMatrix J(out.second.size(),
                        out.second.empty() ? 0 : out.second[0].size());
  for (std::size_t r = 0; r < out.second.size(); ++r)
    for (std::size_t c = 0; c < out.second[r].size(); ++c) J(r, c) = out.second[r][c];
  return Rcpp::List::create(
      Rcpp::Named("values")   = Rcpp::wrap(out.first),
      Rcpp::Named("jacobian") = J,
      Rcpp::Named("names")    = StratD::field_names());
}

// [[Rcpp::export]]
Rcpp::List k93_entry_gradient(std::vector<int> target_idx, double b_0 = 0.059,
                              double birth_rate = 20.0,
                              double max_patch_lifetime = 40.0, int metric = 0) {
  auto setup = [&](K93_Strategy_<double>& s) { s.pars.b_0 = b_0; };
  return run_entry<K93_Strategy_>(target_idx, birth_rate, max_patch_lifetime, metric, setup);
}
