// [[Rcpp::plugins(cpp20)]]
// Isolate the separable_field's query derivative in FF16's CROWN pattern:
// a focal plant reads the field at z = node * H_focal (query height tied to its
// own active height), summed over crown nodes. Compare the separable field's AD
// d/dtheta against a DIRECT O(N^2) sum A(z)=sum_j w_j Q(z/H_j) and against FD.
// If separable-AD != direct-AD, the separable factoring is wrong in this usage.
#include <Rcpp.h>
#include <plant/canopy_shape.h>
#include <odelia/gradient.hpp>
#include <odelia/separable_field.hpp>
#include <array>
#include <vector>
#include <algorithm>

using namespace plant;

// One "crown assimilation"-like reduction: sum over crown nodes of the optical
// depth A(z=node*H_focal), reading shading from ALL sources taller than z. The
// sources are the stand cohorts; every height and weight scales with theta (as
// a single-species stand's cohorts all move together when a shared trait moves).
template <class S>
static S reduce_separable(S theta, const std::vector<double>& h0,
                          const std::vector<double>& w0, std::size_t foc,
                          const CanopyShape& canopy) {
  const std::size_t n = h0.size();
  // Active heights + source weights (descending in height already).
  std::vector<S> H(n);
  std::array<std::vector<S>, 3> sw;
  for (auto& v : sw) v.resize(n);
  std::vector<double> heights_d(n);
  for (std::size_t j = 0; j < n; ++j) {
    H[j] = S(h0[j]) * theta;
    heights_d[j] = h0[j] * xad::value(theta);  // passive positions for the rank
    const auto b = canopy.template shading_source_factors<S>(H[j]);
    for (int p = 0; p < 3; ++p) sw[p][j] = S(w0[j]) * theta * b[p];
  }
  odelia::separable_field<S, 3> field;
  field.assemble(sw);

  const S Hf = H[foc];
  S acc = 0.0;
  for (double node : {0.2, 0.4, 0.6, 0.8}) {
    const S z = S(node) * Hf;                 // query tied to focal height
    const double zp = xad::value(z);
    // rank = #sources with height >= z
    std::size_t k = 0;
    for (std::size_t j = 0; j < n; ++j) if (heights_d[j] >= zp) ++k;
    if (k == 0) continue;
    const S A = field.at(canopy.template shading_query_factors<S>(z), k - 1);
    acc += A;
  }
  return acc;
}

template <class S>
static S reduce_direct(S theta, const std::vector<double>& h0,
                       const std::vector<double>& w0, std::size_t foc,
                       const CanopyShape& canopy) {
  const std::size_t n = h0.size();
  std::vector<S> H(n), W(n);
  for (std::size_t j = 0; j < n; ++j) { H[j] = S(h0[j]) * theta; W[j] = S(w0[j]) * theta; }
  const S Hf = H[foc];
  S acc = 0.0;
  for (double node : {0.2, 0.4, 0.6, 0.8}) {
    const S z = S(node) * Hf;
    for (std::size_t j = 0; j < n; ++j) {
      // shade only from sources at least as tall as z
      if (xad::value(H[j]) >= xad::value(z))
        acc += W[j] * canopy.template Q<S>(z / H[j]);
    }
  }
  return acc;
}

// [[Rcpp::export]]
Rcpp::List field_crown_probe() {
  using RevS = xad::adj<double>::active_type;
  CanopyShape canopy; canopy.initialise(12.0);  // deep-crown, eta=12

  // A little descending-height stand; focal is a mid cohort.
  std::vector<double> h0 = {3.0, 2.2, 1.5, 1.0, 0.6};
  std::vector<double> w0 = {0.4, 0.5, 0.6, 0.5, 0.3};
  std::size_t foc = 2;  // height 1.5

  double val_sep=0, grad_sep=0, val_dir=0, grad_dir=0;
  {
    xad::Tape<double> tape;
    RevS th = 1.0; tape.registerInput(th); tape.newRecording();
    RevS out = reduce_separable<RevS>(th, h0, w0, foc, canopy);
    tape.registerOutput(out); xad::derivative(out) = 1.0; tape.computeAdjoints();
    val_sep = xad::value(out); grad_sep = xad::derivative(th);
  }
  {
    xad::Tape<double> tape;
    RevS th = 1.0; tape.registerInput(th); tape.newRecording();
    RevS out = reduce_direct<RevS>(th, h0, w0, foc, canopy);
    tape.registerOutput(out); xad::derivative(out) = 1.0; tape.computeAdjoints();
    val_dir = xad::value(out); grad_dir = xad::derivative(th);
  }

  // FD on the direct double form.
  auto val_d = [&](double th, auto reduce){ return reduce(th, h0, w0, foc, canopy); };
  double d = 1e-6;
  double fd_sep = (val_d(1.0+d, reduce_separable<double>) - val_d(1.0-d, reduce_separable<double>))/(2*d);
  double fd_dir = (val_d(1.0+d, reduce_direct<double>)   - val_d(1.0-d, reduce_direct<double>))  /(2*d);

  return Rcpp::List::create(
    Rcpp::Named("val_sep")=val_sep, Rcpp::Named("grad_sep")=grad_sep, Rcpp::Named("fd_sep")=fd_sep,
    Rcpp::Named("val_dir")=val_dir, Rcpp::Named("grad_dir")=grad_dir, Rcpp::Named("fd_dir")=fd_dir);
}
