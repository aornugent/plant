// [[Rcpp::plugins(cpp20)]]
// FF16 is forward-mode instantiable: its whole rate path compiles and runs at
// the active reverse scalar AND at the nested forward-over-reverse type that the
// rebind alias enables. Instantiating Species<FF16_active>::compute_rates pulls
// in both the Species-level geometric-compression branch and, through
// Node::compute_rates -> growth_rate_gradient, the nested forward-over-reverse
// block (Fwd = tangent_of<active>). This is the compile-time guard for the FF16
// census-gradient port: the heavy assimilation quadrature (QK), the light
// spline, and the birth-height lift all have to strip every AD layer where they
// narrow to double, or this translation unit fails to build.
#include <Rcpp.h>
#include <plant/models/ff16_strategy.h>
#include <plant/individual.h>
#include <plant/node.h>
#include <plant/species.h>
#include <odelia/gradient.hpp>

using namespace plant;

template <class S>
static double build_and_compute(bool geometric) {
  using Strat = FF16_Strategy_<S>;
  using Env   = FF16_Environment_<S>;

  Strat s;
  s.control.node_geometric_compression = geometric;

  Species<Strat, Env> sp(s);

  Env env;
  env.set_fixed_environment(1.0);  // full light, fixed field

  // Two cohorts so the geometric neighbour-difference has neighbours.
  sp.introduce_new_node(0.0, 1.0);
  sp.introduce_new_node(0.05, 1.0);

  sp.compute_rates(env, 1.0, 1.0);
  return odelia::util::to_passive(sp.height_max());
}

// [[Rcpp::export]]
Rcpp::List ff16_forward_instantiable_smoke() {
  using ActiveS = xad::adj<double>::active_type;

  // Reverse scalar, geometric off -> instantiates the nested forward-over-reverse
  // block (Fwd = tangent_of<ActiveS>) inside Node::growth_rate_gradient.
  const double a = build_and_compute<ActiveS>(false);
  // Reverse scalar, geometric on -> instantiates the Species geometric branch.
  const double b = build_and_compute<ActiveS>(true);
  // Double path (sanity; must match the active values, which carry the same
  // primal since no parameter is seeded here).
  const double c = build_and_compute<double>(false);
  const double d = build_and_compute<double>(true);

  return Rcpp::List::create(
      Rcpp::Named("height_max_active_stencil")   = a,
      Rcpp::Named("height_max_active_geometric") = b,
      Rcpp::Named("height_max_double_stencil")   = c,
      Rcpp::Named("height_max_double_geometric") = d);
}
