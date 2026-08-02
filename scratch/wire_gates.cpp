// Gates for the leaf Jacobian wired into the cohort block: the block's adjoint
// against a finite difference of the same block, the knot accumulation, and the
// decomposition against one whole-patch recording.
#include <Rcpp.h>
#include <plant/models/tf24_strategy.h>
#include <plant/patch.h>
#include <odelia/gradient.hpp>
#include <cmath>
#include <vector>

using strategy_type = plant::TF24_Strategy<double>;
using environment_type = plant::TF24_Environment<double>;
using patch_type = plant::Patch<strategy_type, environment_type>;
using scalar = odelia::ode::active_scalar<double>;
using active_strategy = plant::TF24_Strategy<scalar>;
using active_environment = plant::TF24_Environment<scalar>;
using active_individual = plant::Individual<active_strategy, active_environment>;
using active_patch = plant::Patch<active_strategy, active_environment>;

static patch_type* as_patch(SEXP obj) {
  Rcpp::Environment r6(obj);
  return Rcpp::as<Rcpp::XPtr<patch_type> >(r6.get(".ptr"));
}

// One cohort's block at the active scalar, built the way cohort_block_adjoint
// builds it: a per-block copy from a template no recording has touched.
struct block_at {
  strategy_type tmpl;
  environment_type env_tmpl;
  std::vector<double> in;
  size_t n_out;

  block_at(patch_type& p, size_t species_index, size_t node_index)
      : tmpl(*p.r_species()[species_index].strategy_ptr()),
        env_tmpl(p.r_environment()) {
    typename strategy_type::ptr s = std::make_shared<strategy_type>(tmpl);
    plant::Individual<strategy_type, environment_type> ind(s);
    in.resize(ind.block_input_size(env_tmpl));
    n_out = ind.block_output_size(env_tmpl);
    p.r_species()[species_index].node_at(node_index).individual.block_inputs(
        in.begin(), env_tmpl);
  }

  void run_active(const std::vector<scalar>& x, std::vector<scalar>& y) {
    typename active_strategy::ptr s =
        std::make_shared<active_strategy>(tmpl.rebind_from<scalar>());
    active_individual ind(s);
    active_environment e = env_tmpl.rebind_from<scalar>();
    e.light_availability.spline.set_nodes(
        env_tmpl.light_availability.spline.knots());
    ind.set_block_inputs(x.begin(), e);
    ind.compute_rates(e);
    ind.block_outputs(y.begin());
  }

  void run_double(const std::vector<double>& x, std::vector<double>& y) {
    typename strategy_type::ptr s = std::make_shared<strategy_type>(tmpl);
    plant::Individual<strategy_type, environment_type> ind(s);
    environment_type e = env_tmpl;
    ind.set_block_inputs(x.begin(), e);
    ind.compute_rates(e);
    ind.block_outputs(y.begin());
  }
};

// [[Rcpp::export]]
Rcpp::List block_vjp(SEXP obj, int species_index, int node_index,
                     Rcpp::NumericVector out_adjoint) {
  patch_type& p = *as_patch(obj);
  block_at b(p, species_index, node_index);
  std::vector<double> lam(out_adjoint.begin(), out_adjoint.end());
  std::vector<double> adj;
  odelia::ode::vector_jacobian_product(
      b.in, lam, [&](const std::vector<scalar>& x, std::vector<scalar>& y) -> void {
        b.run_active(x, y);
      }, adj);
  return Rcpp::List::create(Rcpp::_["inputs"] = b.in, Rcpp::_["adjoint"] = adj,
                            Rcpp::_["n_out"] = (int)b.n_out);
}

// The same block at double, so a difference reads the identical code path.
// [[Rcpp::export]]
Rcpp::NumericVector block_value(SEXP obj, int species_index, int node_index,
                                Rcpp::NumericVector x_in) {
  patch_type& p = *as_patch(obj);
  block_at b(p, species_index, node_index);
  std::vector<double> x(x_in.begin(), x_in.end()), y(b.n_out, 0.0);
  b.run_double(x, y);
  return Rcpp::wrap(y);
}

// [[Rcpp::export]]
Rcpp::List patch_adjoint(SEXP obj, Rcpp::NumericVector lambda_dydt) {
  patch_type& p = *as_patch(obj);
  std::vector<double> lam(lambda_dydt.begin(), lambda_dydt.end());
  std::vector<double> out(p.ode_size(), 0.0);
  p.clear_trait_adjoint();
  p.ode_rates_adjoint(lam.begin(), out.begin());
  return Rcpp::List::create(Rcpp::_["lambda_y"] = out,
                            Rcpp::_["trait"] = p.trait_adjoint,
                            Rcpp::_["sweeps"] = (double)p.block_sweeps,
                            Rcpp::_["recording"] = (double)p.block_recording_size);
}

// The decomposition one contribution at a time, in the order ode_rates_adjoint
// applies them. `upto` counts contributions: 1 soil, 2 offspring, 3 the cohort
// blocks, 4 the knot pullback, 5 the allometry.
// [[Rcpp::export]]
Rcpp::List patch_adjoint_partial(SEXP obj, Rcpp::NumericVector lambda_dydt,
                                 int upto) {
  patch_type& p = *as_patch(obj);
  const size_t n = p.ode_size();
  const size_t n_resource = p.r_environment().n_resources();
  std::vector<double> lam(lambda_dydt.begin(), lambda_dydt.end());
  std::vector<double> lambda_state(n, 0.0);
  typename patch_type::block_seeds seeds{
      std::vector<double>(p.node_count() * strategy_type::state_size(), 0.0),
      std::vector<double>(p.node_count() * n_resource, 0.0)};
  for (size_t k = 0; k < p.node_count(); ++k) {
    for (size_t s = 0; s < strategy_type::state_size(); ++s) {
      seeds.rate[k * strategy_type::state_size() + s] =
          lam[k * patch_type::node_type::ode_size() + s];
    }
  }
  if (upto >= 1) p.soil_adjoint(lam, lambda_state, seeds);
  if (upto >= 2) p.offspring_adjoint(lam, lambda_state, seeds);
  typename patch_type::light_knot_adjoints knot{
      std::vector<double>(p.r_environment().light_availability.spline.knots().size(), 0.0),
      std::vector<double>(p.r_environment().light_availability.spline.knots().size(), 0.0)};
  p.clear_trait_adjoint();
  if (upto >= 3) p.cohort_block_adjoint(seeds, lambda_state, knot);
  std::vector<plant::node_size_adjoints> sizes(
      p.node_count(), plant::node_size_adjoints{0, 0, 0});
  if (upto >= 4) p.light_knot_adjoint(knot, sizes);
  if (upto >= 5) p.allometry_adjoint(sizes, lambda_state);
  return Rcpp::List::create(Rcpp::_["lambda_y"] = lambda_state,
                            Rcpp::_["knot_value"] = knot.value,
                            Rcpp::_["knot_slope"] = knot.slope,
                            Rcpp::_["trait"] = p.trait_adjoint);
}

// One whole-patch recording at the current state: the reference the
// decomposition is compared against.
// [[Rcpp::export]]
Rcpp::NumericVector patch_recording_vjp(SEXP obj, Rcpp::NumericVector lambda_dydt) {
  patch_type& p = *as_patch(obj);
  std::vector<double> y(p.ode_size());
  p.ode_state(y.begin());
  const double time = p.time();
  std::vector<double> lam(lambda_dydt.begin(), lambda_dydt.end()), adj;
  active_patch ap = p.rebind_from<scalar>();
  odelia::ode::vector_jacobian_product(
      y, lam, [&](const std::vector<scalar>& x, std::vector<scalar>& out) -> void {
        ap.set_ode_state(x.begin(), time);
        ap.ode_rates(out.begin());
      }, adj);
  return Rcpp::wrap(adj);
}

// The per-cohort knot-value contributions the accumulation sums, taken one
// cohort at a time so the sum can be asserted against it.
// [[Rcpp::export]]
Rcpp::NumericMatrix knot_contributions(SEXP obj, Rcpp::NumericVector lambda_dydt) {
  patch_type& p = *as_patch(obj);
  const size_t n = p.ode_size();
  const size_t n_resource = p.r_environment().n_resources();
  const size_t n_state = strategy_type::state_size();
  const size_t n_knot = p.r_environment().light_availability.knot_count();
  std::vector<double> lam(lambda_dydt.begin(), lambda_dydt.end());
  std::vector<double> lambda_state(n, 0.0);
  typename patch_type::block_seeds seeds{
      std::vector<double>(p.node_count() * n_state, 0.0),
      std::vector<double>(p.node_count() * n_resource, 0.0)};
  for (size_t k = 0; k < p.node_count(); ++k) {
    for (size_t s = 0; s < n_state; ++s) {
      seeds.rate[k * n_state + s] = lam[k * patch_type::node_type::ode_size() + s];
    }
  }
  p.soil_adjoint(lam, lambda_state, seeds);
  p.offspring_adjoint(lam, lambda_state, seeds);

  Rcpp::NumericMatrix out(p.node_count(), n_knot);
  size_t k = 0;
  std::vector<patch_type::species_type> sp = p.r_species();
  for (size_t i = 0; i < sp.size(); ++i) {
    for (size_t j = 0; j < sp[i].size(); ++j, ++k) {
      block_at b(p, i, j);
      std::vector<double> lam_out(b.n_out, 0.0);
      for (size_t s = 0; s < n_state; ++s) lam_out[s] = seeds.rate[k * n_state + s];
      for (size_t r = 0; r < n_resource; ++r) {
        lam_out[n_state + r] = seeds.uptake[k * n_resource + r];
      }
      std::vector<double> adj;
      odelia::ode::vector_jacobian_product(
          b.in, lam_out,
          [&](const std::vector<scalar>& x, std::vector<scalar>& y) -> void {
            b.run_active(x, y);
          }, adj);
      for (size_t c = 0; c < n_knot; ++c) out(k, c) = adj[n_state + c];
    }
  }
  return out;
}

// Where the collar solve of one block lands: the residual at the point it left,
// the recorded curvature, and the pinned flag. Taken at the three points a
// central difference of one input visits.
// [[Rcpp::export]]
Rcpp::DataFrame collar_residuals(SEXP obj, int species_index, int node_index,
                                 int input, double h) {
  patch_type& p = *as_patch(obj);
  block_at b(p, species_index, node_index);
  Rcpp::NumericVector at(3), collar(3), resid(3), curv(3);
  Rcpp::LogicalVector pinned(3);
  const double base = b.in[input];
  for (int k = 0; k < 3; ++k) {
    std::vector<double> x = b.in, y(b.n_out, 0.0);
    x[input] = base + (k == 0 ? 0.0 : (k == 1 ? h : -h));
    typename strategy_type::ptr s = std::make_shared<strategy_type>(b.tmpl);
    plant::Individual<strategy_type, environment_type> ind(s);
    environment_type e = b.env_tmpl;
    ind.set_block_inputs(x.begin(), e);
    ind.compute_rates(e);
    plant::Leaf& l = s->leaf;
    at[k] = x[input];
    collar[k] = l.root_collar_psi_;
    curv[k] = l.dR_dcollar_;
    pinned[k] = l.collar_pinned_;
    resid[k] = l.dprofit_droot_collar_psi(-l.root_collar_psi_);
  }
  return Rcpp::DataFrame::create(Rcpp::_["x"] = at, Rcpp::_["collar"] = collar,
                                 Rcpp::_["residual"] = resid,
                                 Rcpp::_["dR_dcollar"] = curv,
                                 Rcpp::_["pinned"] = pinned);
}

// The collar curvature the reverse pass divides by, taken over a range of
// difference steps, beside the value the solve recorded.
// [[Rcpp::export]]
Rcpp::DataFrame collar_curvature_sweep(SEXP obj, int species_index, int node_index,
                                       Rcpp::NumericVector hs) {
  patch_type& p = *as_patch(obj);
  block_at b(p, species_index, node_index);
  std::vector<double> x = b.in, y(b.n_out, 0.0);
  typename strategy_type::ptr s = std::make_shared<strategy_type>(b.tmpl);
  plant::Individual<strategy_type, environment_type> ind(s);
  environment_type e = b.env_tmpl;
  ind.set_block_inputs(x.begin(), e);
  ind.compute_rates(e);
  plant::Leaf& l = s->leaf;
  const double pmag = -l.root_collar_psi_;
  const double recorded = l.dR_dcollar_;
  Rcpp::NumericVector out(hs.size()), rec(hs.size());
  for (int i = 0; i < hs.size(); ++i) {
    out[i] = l.dR_dcollar_at(pmag, hs[i]);
    rec[i] = recorded;
  }
  return Rcpp::DataFrame::create(Rcpp::_["h"] = hs, Rcpp::_["dR_dcollar"] = out,
                                 Rcpp::_["recorded"] = rec);
}
