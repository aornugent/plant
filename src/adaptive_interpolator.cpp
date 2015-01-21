#include <tree2/adaptive_interpolator.h>
#include <Rcpp.h>
#include <tree2/util_post_rcpp.h>

namespace tree2 {
namespace interpolator {

// Given our current set of x/y points, construct the interpolated
// fucntion.  This is the only function that actually modifies the
// interpolator object.
void AdaptiveInterpolator::update_spline() {
  std::vector<double> x(xx.begin(), xx.end());
  std::vector<double> y(yy.begin(), yy.end());
  interpolator.init(x, y);
}

void AdaptiveInterpolator::check_bounds(double a, double b) {
  if (a >= b) {
    Rcpp::stop("Impossible bounds");
  }
  if (!util::is_finite(a) || !util::is_finite(b)) {
    Rcpp::stop("Infinite bounds");
  }
}

// Determines if difference between predicted and true values falls
// within error bounds.
bool AdaptiveInterpolator::check_err(double y_true, double y_pred) const {
  const double err_abs = fabs(y_true - y_pred);
  const double err_rel = fabs(1 - y_pred / y_true);
  return err_abs < atol || err_rel < rtol;
}

}
}

// [[Rcpp::export]]
tree2::interpolator::Interpolator
test_adaptive_interpolator(Rcpp::Function f, double a, double b) {
  tree2::util::RFunctionWrapper fw(f);
  const double atol = 1e-6, rtol = 1e-6;
  const int nbase = 17, max_depth = 16;
  tree2::interpolator::AdaptiveInterpolator
    generator(atol, rtol, nbase, max_depth);
  return generator.construct(fw, a, b);
}
