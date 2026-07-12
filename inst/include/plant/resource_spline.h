// -*-c++-*-

#ifndef PLANT_PLANT_RESOURCE_SPLINE_H_
#define PLANT_PLANT_RESOURCE_SPLINE_H_

#include <odelia/interpolator.hpp>
#include <odelia/ode_interface.hpp>
#include <plant/util.h>
#include <algorithm> // std::max, for the resource-availability floor (#253)

using namespace Rcpp;

namespace plant {

// Interpolated resource availability (light) as a function of size. Templated on
// the scalar S of the knot VALUES; knot POSITIONS stay double. S = double is the
// production path (the `ResourceSpline` alias below); an active S makes the
// interpolated light differentiable w.r.t. the knot values, for the resident
// self-shading gradient. The point read also accepts an active query height, so
// d(light)/d(height) flows even when the values are frozen double (mutant /
// fixed environment). The odelia interpolator owns the adaptive refiner.
template <class S = double>
class ResourceSpline_ {
public:

  // Constructors
  ResourceSpline_() {
    setup(1e-6, 17, 16, false);
  }

  ResourceSpline_(double tol, size_t nbase, size_t max_depth, bool rescale_usually) {
    setup(tol, nbase, max_depth, rescale_usually);
  }

  void setup(double tol, size_t nbase, size_t max_depth, bool rescale_usually) {
    // Adaptive-refinement controls, handed to the odelia interpolator's
    // construct() from compute_environment (the interpolator owns the refiner).
    spline_tol = tol;
    spline_nbase = nbase;
    spline_max_depth = max_depth;
    // Initialise to a flat, fully-open spline (light = 1 everywhere over [0,1]),
    // matching an empty environment; over-written by the first
    // compute_environment() before any real light read. A constant target
    // reproduces the previous dummy (which read the not-yet-built spline and so
    // saw the open value 1.0) without reading the spline as it is built.
    spline.construct([](double) { return 1.0; }, 0.0, 1.0,
                     spline_tol, spline_tol, spline_nbase, spline_max_depth);
    spline_rescale_usually = rescale_usually;
  };

  template <typename Function>
  void compute_environment(Function f_compute_competition, double height_max, bool rescale) {
    if (rescale & spline_rescale_usually) {
      rescale_spline(f_compute_competition, height_max);
    } else {
      construct_spline(f_compute_competition, height_max);
    }
  };

  void set_fixed_value(double value, double height_max) {
    std::vector<double> x = {0, height_max/2.0, height_max};
    std::vector<S> y = {static_cast<S>(value), static_cast<S>(value),
                        static_cast<S>(value)};
    clear();
    spline.init(x, y);
  }

  void clear() {
    spline.clear();
  }

  // Highest height covered by the spline; above this get_value_at_height()
  // returns the hard-coded open value (1.0). Hoist this out of hot per-point
  // loops with get_value_at_height(height, cap).
  double max_height() const { return spline.max(); }

  S get_value_at_height(S height) const {
    return get_value_at_height(height, static_cast<S>(spline.max()));
  }

  // Variant taking a pre-fetched cap (= max_height()) so callers integrating
  // over many points pay the spline.max() lookup once rather than per point.
  S get_value_at_height(S height, S cap) const {
    // TODO(#385): change maximum - here hard-coded to 1.0
    // `cap` already guards the upper bound and the crown integral keeps
    // height >= 0 = spline.min(), so use the unchecked operator() rather than
    // eval() to avoid re-running check_active()/bound checks per quadrature
    // point. Same underlying tk_spline(height) call.
    //
    // Floor the result at 0 (#253): the cubic spline can undershoot below
    // zero between knots (notably the K93 light spline at high k_I), which is
    // non-physical for a resource availability. The clamp is a no-op for the
    // usual positive case (so values stay bit-identical there) and only bites
    // on the spurious negative undershoot. This is the single chokepoint for
    // FF16/K93/TF24, and belongs here in plant rather than in the
    // general-purpose interpolator (which is migrating to odelia).
    return height <= cap ? std::max(static_cast<S>(0.0), spline(height))
                         : static_cast<S>(1.0);
  }

  virtual void r_init_interpolators(const std::vector<double>& state) {
    // See issue #144; this is important as we have to at least refine
    // the light environment, but doing this is better because it means
    // that if rescale_usually is on we do get the same light
    // environment as before.
    if (state.size() % 2 != 0) {
      util::stop("Expected even number of elements in light environment");
    }
    const size_t state_n = state.size() / 2;
    auto it = state.begin();
    std::vector<double> state_x;
    std::vector<S> state_y;
    std::copy_n(it,         state_n, std::back_inserter(state_x));
    std::copy_n(it + state_n, state_n, std::back_inserter(state_y));
    spline.init(state_x, state_y);
  }

  // Stores the interpolator spline of resource availability as a function of
  // size. The odelia interpolator owns the adaptive refiner (construct); the
  // controls below are handed to it. Positions double, values S.
  odelia::interpolator::basic_interpolator<S> spline;

  // Adaptive-refinement controls for spline.construct().
  double spline_tol;
  size_t spline_nbase;
  size_t spline_max_depth;

  // flag, do we try to rescale the spline when possible? this is quicker
  bool spline_rescale_usually;

  Rcpp::NumericMatrix r_get_state() const
  {

    // format spline as Matrix; only double crosses to R, so narrow the (possibly
    // active) knot values via xad::value (a no-op on the double path).
    std::vector<std::vector<double>> xy;
    xy.push_back(spline.get_x());
    std::vector<double> y_values;
    for (auto const& v : spline.get_y()) {
      y_values.push_back(xad::value(v));
    }
    xy.push_back(y_values);

    const size_t n = xy.size();
    Rcpp::NumericMatrix ret(static_cast<int>(xy.begin()->size()),
                            static_cast<int>(n));
    Rcpp::NumericMatrix::iterator it = ret.begin();
    for (size_t i = 0; i < n; ++i)
    {
      it = std::copy(xy[i].begin(), xy[i].end(), it);
    }

    // Add colnames
    ret.attr("dimnames") = Rcpp::List::create(R_NilValue, Rcpp::CharacterVector::create("height", "light_availability"));
    return ret;
  }

private:

  template <typename Function>
  void construct_spline(Function f_compute_competition, double height_max)
  {
    const double lower_bound = 0.0;
    double upper_bound = height_max;

    spline.construct(f_compute_competition, lower_bound, upper_bound,
                     spline_tol, spline_tol, spline_nbase, spline_max_depth);
  }

  template <typename Function>
  void rescale_spline(Function f_compute_competition, double height_max) {
    std::vector<double> h = spline.get_x();
    const double min = spline.min(), // 0.0?
      height_max_old = spline.max();

    util::rescale(h.begin(), h.end(), min, height_max_old, min, height_max);
    h.back() = height_max; // Avoid round-off error.

    spline.clear();
    for (auto hi : h) {
      spline.add_point(hi, f_compute_competition(hi));
    }
    spline.initialise();
  }

  };

// The double instantiation is the production path bound by RcppR6 as
// `plant::ResourceSpline` (ResourceSpline getters/setters, the three
// environments' light_availability member).
using ResourceSpline = ResourceSpline_<double>;

} // plant namespace

#endif
