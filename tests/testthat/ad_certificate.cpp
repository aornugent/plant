// [[Rcpp::plugins(cpp20)]]
// AD completeness certificate: reverse-mode gradient over EVERY registered
// AD_FIELDS leaf vs a per-field pinned-schedule central FD, for FF16 and K93.
// Per leaf: intact (AD==FD, nonzero) | structural-zero (both 0) | SEVERED
// (AD==0, FD!=0) | KINK/partial (both nonzero, ratio off). Certifies every
// registered leaf's pathway in one sweep -- what per-site reading cannot.
#include <Rcpp.h>
#include <plant/models/ff16_strategy.h>
#include <plant/models/ff16_environment.h>
#include <plant/models/k93_strategy.h>
#include <plant/models/k93_environment.h>
#include <plant/models/tf24_strategy.h>
#include <plant/models/tf24_environment.h>
#include <plant/individual.h>
#include <plant/parameters.h>
#include <plant/control.h>
#include <plant/scm.h>
#include <odelia/gradient.hpp>
#include <cmath>
#include <string>
#include <vector>

using namespace plant;

// One reducer for both strategies.
template <class Patch>
static auto reduce_metric(Patch& p, int metric)
    -> std::decay_t<decltype(p.compute_competition(0.0))> {
  using V = std::decay_t<decltype(p.compute_competition(0.0))>;
  if (metric == 2) {
    V h = 0.0;
    for (std::size_t s = 0; s < p.size(); ++s) { auto& sp = p.at_species(s);
      for (auto it = sp.node_begin(); it != sp.node_end(); ++it) h += it->height(); }
    return h;
  }
  return p.offspring_production()[0];
}

// do_ad: run the one reverse-AD sweep (all fields at once). fd_fields: the
// (0-based) field indices to central-difference; empty => none. Splitting the
// two halves lets a caller parallelise the embarrassingly-parallel FD sweep
// (one independent pinned SCM run per field per sign) across processes.
template <template <class> class StratT, class Setup>
static Rcpp::List sweep(Rcpp::List node_sched, std::vector<double> ode_times_in,
                        double max_patch_lifetime, int metric, double fd_rel,
                        double birth_rate, Setup setup,
                        bool do_ad, std::vector<int> fd_fields) {
  using RevS = xad::adj<double>::active_type;
  std::vector<std::vector<double>> nst;
  for (R_xlen_t i = 0; i < node_sched.size(); ++i)
    nst.push_back(Rcpp::as<std::vector<double>>(node_sched[i]));

  const int N = (int)StratT<double>().field_ptrs().size();
  std::vector<std::string> names = StratT<double>::field_names();
  std::vector<int> idx(N); for (int i = 0; i < N; ++i) idx[i] = i;

  auto make_params = [&](auto tag, int perturb, double dv) {
    using S = std::decay_t<decltype(tag)>;
    using Env = typename StratT<S>::environment_type;
    StratT<S> s; setup(s);
    s.is_variable_birth_rate = false; s.birth_rate_y = {birth_rate};
    if (perturb >= 0) *s.field_ptrs()[perturb] += S(dv);
    Parameters<StratT<S>, Env> p;
    p.strategies.push_back(s); p.max_patch_lifetime = max_patch_lifetime;
    p.node_schedule_times = nst; p.ode_times = ode_times_in; p.validate();
    return p;
  };
  Control ctrl;
  auto pin = [&](auto& scm){ scm.reset(); NodeSchedule ns=scm.r_node_schedule();
                             ns.r_set_use_ode_times(true); scm.r_set_node_schedule(ns); };

  // Reverse AD: all leaves in one sweep.
  std::vector<double> g;
  if (do_ad) {
    using Env = typename StratT<RevS>::environment_type;
    Env env;
    SCM<StratT<RevS>, Env> scm(make_params(RevS(), -1, 0.0), env, ctrl);
    pin(scm);
    odelia::ode::DifferentiationTargets t; t.params = idx;
    auto ptrs = scm.get_system_ref().ad_parameters();
    for (int i : idx) t.values.push_back(xad::value(*ptrs[i]));
    auto functional = [metric](decltype(scm)& s) -> RevS { return reduce_metric(s.get_system_ref(), metric); };
    g = odelia::ode::compute_gradient(scm, t, functional).second;
  }

  // Per-field pinned-schedule central FD (double), over the requested subset.
  auto run_double = [&](int perturb, double dv) -> double {
    using Env = typename StratT<double>::environment_type;
    Env env;
    SCM<StratT<double>, Env> scm(make_params(double(), perturb, dv), env, ctrl);
    pin(scm); scm.run();
    return reduce_metric(scm.get_system_ref(), metric);
  };
  std::vector<double> base(N);
  { StratT<double> s0; setup(s0); for (int i=0;i<N;++i) base[i] = *s0.field_ptrs()[i]; }
  std::vector<double> fd(fd_fields.size());
  for (std::size_t k = 0; k < fd_fields.size(); ++k) {
    int i = fd_fields[k];
    double d = fd_rel * (std::abs(base[i]) + 1e-3);
    fd[k] = (run_double(i, d) - run_double(i, -d)) / (2.0 * d);
  }
  return Rcpp::List::create(Rcpp::Named("names")=names, Rcpp::Named("base")=base,
                            Rcpp::Named("ad")=g, Rcpp::Named("fd")=fd,
                            Rcpp::Named("fd_fields")=fd_fields);
}

// Build the full 0..N-1 field-index list for a strategy (all fields).
template <template <class> class StratT>
static std::vector<int> all_fields() {
  const int N = (int)StratT<double>().field_ptrs().size();
  std::vector<int> v(N); for (int i = 0; i < N; ++i) v[i] = i; return v;
}

// [[Rcpp::export]]
Rcpp::List ff16_allfield(Rcpp::List node_sched, std::vector<double> ode_times_in,
                         double lma = 0.1978791, double birth_rate = 20.0,
                         double max_patch_lifetime = 40.0, int metric = 0, double fd_rel = 3e-4) {
  auto setup = [&](auto& s){ auto p=s.field_ptrs(); auto nm=std::decay_t<decltype(s)>::field_names();
    for (std::size_t i=0;i<nm.size();++i) if (nm[i]=="lma") *p[i]=lma; };
  return sweep<FF16_Strategy_>(node_sched, ode_times_in, max_patch_lifetime, metric, fd_rel, birth_rate, setup,
                               true, all_fields<FF16_Strategy_>());
}

// [[Rcpp::export]]
Rcpp::List k93_allfield(Rcpp::List node_sched, std::vector<double> ode_times_in,
                        double b_0 = 0.059, double birth_rate = 20.0,
                        double max_patch_lifetime = 40.0, int metric = 0, double fd_rel = 3e-4) {
  auto setup = [&](auto& s){ s.pars.b_0 = b_0; };
  return sweep<K93_Strategy_>(node_sched, ode_times_in, max_patch_lifetime, metric, fd_rel, birth_rate, setup,
                              true, all_fields<K93_Strategy_>());
}

// [[Rcpp::export]]
std::vector<std::string> tf24_field_names() { return TF24_Strategy_<double>::field_names(); }

// TF24: the reoptimising FD is just the standard double SCM run (the leaf re-solves
// its collar-psi optimum every compute_rates), so comparing reverse-AD (which uses
// the frozen-p* supplied_derivative seam) to this FD directly tests plant#60.
// field_values (length N, field-name order) injects the R hyperpar-resolved params
// -- TF24's Leaf hydraulics are only valid at the resolved operating point.
// [[Rcpp::export]]
Rcpp::List tf24_allfield(Rcpp::List node_sched, std::vector<double> ode_times_in,
                         std::vector<double> field_values, double birth_rate = 20.0,
                         double max_patch_lifetime = 10.0, int metric = 2, double fd_rel = 3e-4) {
  auto setup = [&](auto& s){ auto p=s.field_ptrs();
    for (std::size_t i=0;i<p.size() && i<field_values.size();++i)
      *p[i] = std::decay_t<decltype(*p[i])>(field_values[i]); };
  return sweep<TF24_Strategy_>(node_sched, ode_times_in, max_patch_lifetime, metric, fd_rel, birth_rate, setup,
                               true, all_fields<TF24_Strategy_>());
}

// Split TF24 entries so the caller can parallelise the FD sweep across processes.
// tf24_ad: the single reverse-AD sweep (all fields, no FD).
// tf24_fd: the pinned-schedule central FD for the requested (0-based) field subset.
// Together they reproduce tf24_allfield exactly; run tf24_fd chunks under mclapply.
static auto tf24_setup(std::vector<double> const& field_values) {
  return [field_values](auto& s){ auto p=s.field_ptrs();
    for (std::size_t i=0;i<p.size() && i<field_values.size();++i)
      *p[i] = std::decay_t<decltype(*p[i])>(field_values[i]); };
}

// [[Rcpp::export]]
Rcpp::List tf24_ad(Rcpp::List node_sched, std::vector<double> ode_times_in,
                   std::vector<double> field_values, double birth_rate = 20.0,
                   double max_patch_lifetime = 10.0, int metric = 2) {
  return sweep<TF24_Strategy_>(node_sched, ode_times_in, max_patch_lifetime, metric,
                               3e-4, birth_rate, tf24_setup(field_values), true, {});
}

// [[Rcpp::export]]
Rcpp::List tf24_fd(Rcpp::List node_sched, std::vector<double> ode_times_in,
                   std::vector<double> field_values, std::vector<int> fields,
                   double birth_rate = 20.0, double max_patch_lifetime = 10.0,
                   int metric = 2, double fd_rel = 3e-4) {
  return sweep<TF24_Strategy_>(node_sched, ode_times_in, max_patch_lifetime, metric,
                               fd_rel, birth_rate, tf24_setup(field_values), false, fields);
}
