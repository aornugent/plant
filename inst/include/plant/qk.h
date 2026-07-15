// -*-c++-*-

// The class `QK` provides methods to use Gauss-Kronrod quadrature to numerically integrate a function.
// This is code is ported from GSL
// The integration has several "rules", defined in qk_rules.cpp. These include QK15, QK21, QK31, QK41, QK51. These allow for different numbers of points in the integration.
// General background on the method is available at: https://en.wikipedia.org/wiki/ or
// or https://www.gnu.org/software/gsl/doc/html/integration.html
// These methods can also be accessed via the class QAG, which offers an adaptive extension of the method.
// Implemented by Rich FitzJohn

#ifndef PLANT_PLANT_QK_H_
#define PLANT_PLANT_QK_H_

#include <vector>
#include <cmath> // std::abs, std::isnan
#include <RcppCommon.h> // SEXP
#include <XAD/XAD.hpp> // xad::value (narrow the error estimate to double)
#include <odelia/ode_util.hpp> // odelia::util::to_passive (full AD strip, nesting-safe)
#include <plant/util.h> // util::stop

namespace plant {
namespace quadrature {

class QK {
public:
  QK();
  QK(size_t rule);
  // Fixed Gauss-Kronrod rule templated on the scalar S of the integrand/bounds.
  // A fixed rule places its nodes as a deterministic affine function of the
  // bounds, so an active bound tapes exactly (the integral value carries S,
  // including the half-length Jacobian). The error estimate is diagnostic and
  // off the differentiated path, so it is computed in double via xad::value.
  // S = double is the production path and is bit-identical to the pre-template
  // code. (Cluster 4 / plan §4.4.)
  template <typename S, typename Function>
  S integrate(Function f, S a, S b);

  // These two provide very low level access to the integration
  // routines.
  std::vector<double> integrate_vector_x(double a, double b) const;
  double integrate_vector(const std::vector<double>& y,
                          double a, double b);

  double get_last_area()     const {return last_result;    }
  double get_last_error()    const {return last_error;     }
  double get_last_area_abs() const {return last_result_abs;}
  double get_last_area_asc() const {return last_result_asc;}

  // * R interface
  double r_integrate(SEXP f, double a, double b);

private:
  void initialise(size_t rule);
  static double rescale_error(double err, double result_abs,
                              double result_asc);

  double last_result;
  double last_result_abs;
  double last_result_asc;
  double last_error;

  size_t n;
  std::vector<double> fv1;
  std::vector<double> fv2;
  std::vector<double> xgk, wg, wgk;
};

// Skipping dealing with Functors by just requiring something
// callable using templates.
//
// The integral accumulators and the (bound-dependent) nodes carry S so an
// active bound differentiates through exactly; the abs/error-estimate machinery
// (result_abs/result_asc/last_error, the fv1/fv2 scratch) stays double via
// util::to_passive, both because it is diagnostic and to keep abs()'s kink off
// the tape. to_passive strips every AD layer, so this also compiles when S is a
// nested forward-over-reverse type. At S = double it is the identity, so this is
// bit-identical to the previous double-only body.
template <typename S, typename Function>
S QK::integrate(Function f, S a, S b) {
  const S      center          = 0.5 * (a + b);
  const S      half_length     = 0.5 * (b - a);
  const double abs_half_length = std::abs(odelia::util::to_passive(half_length));
  const S      f_center        = f(center);
  if (std::isnan(odelia::util::to_passive(f_center))) {
    util::stop("Integrand returned NaN at x=" +
               std::to_string(odelia::util::to_passive(center)));
  }

  S result_gauss = 0;
  S result_kronrod = f_center * wgk[n - 1];
  double result_abs = std::abs(odelia::util::to_passive(result_kronrod));

  if (n % 2 == 0) {
    result_gauss = f_center * wg[n / 2 - 1];
  }

  for (size_t j = 0; j < (n - 1) / 2; j++) {
    const size_t jtw = j * 2 + 1;  /* in original fortran j=1,2,3 jtw=2,4,6 */
    const S abscissa = half_length * xgk[jtw];
    const S fval1 = f(center - abscissa);
    const S fval2 = f(center + abscissa);
    const S fsum = fval1 + fval2;
    fv1[jtw] = odelia::util::to_passive(fval1);
    fv2[jtw] = odelia::util::to_passive(fval2);
    result_gauss   += wg[j] * fsum;
    result_kronrod += wgk[jtw] * fsum;
    result_abs     += wgk[jtw] * (std::abs(odelia::util::to_passive(fval1)) +
                                  std::abs(odelia::util::to_passive(fval2)));
  }

  for (size_t j = 0; j < n / 2; j++) {
    size_t jtwm1 = j * 2;
    const S abscissa = half_length * xgk[jtwm1];
    const S fval1 = f(center - abscissa);
    const S fval2 = f(center + abscissa);
    fv1[jtwm1] = odelia::util::to_passive(fval1);
    fv2[jtwm1] = odelia::util::to_passive(fval2);
    result_kronrod += wgk[jtwm1] * (fval1 + fval2);
    result_abs     += wgk[jtwm1] * (std::abs(odelia::util::to_passive(fval1)) +
                                    std::abs(odelia::util::to_passive(fval2)));
  };

  const double mean = odelia::util::to_passive(result_kronrod) * 0.5;

  double result_asc = wgk[n - 1] * std::abs(odelia::util::to_passive(f_center) - mean);

  for (size_t j = 0; j < n - 1; j++) {
    result_asc += wgk[j] * (std::abs(fv1[j] - mean) +
                            std::abs(fv2[j] - mean));
  }

  /* scale by the width of the integration region */

  const S err = (result_kronrod - result_gauss) * half_length;

  result_kronrod *= half_length;
  result_abs *= abs_half_length;
  result_asc *= abs_half_length;

  last_result     = odelia::util::to_passive(result_kronrod);
  last_result_abs = result_abs;
  last_result_asc = result_asc;
  last_error      = rescale_error(odelia::util::to_passive(err), result_abs, result_asc);

  if (std::isnan(last_result)) {
    util::stop("Integrand produced NaN result over [" +
               std::to_string(odelia::util::to_passive(a)) + ", " +
               std::to_string(odelia::util::to_passive(b)) + "]");
  }

  return result_kronrod;
}

}
}

#endif
