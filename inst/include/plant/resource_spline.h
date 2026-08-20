// -*-c++-*-

#ifndef PLANT_PLANT_RESOURCE_SPLINE_H_
#define PLANT_PLANT_RESOURCE_SPLINE_H_

#include <odelia/interpolator.hpp>
#include <odelia/ode_util.hpp>
#include <odelia/ode_interface.hpp>
#include <plant/util.h>
#include <algorithm> // std::max, for the resource-availability floor (#253)
#include <utility>   // std::pair, the value and slope a build supplies

using namespace Rcpp;

namespace plant {

// Templated on the scalar S the resource values carry; the knot positions stay
// double. S = double is production.
//
// The knot values and slopes live on the interpolant and are read back from it.
// They were mirrored here once, on the ground that a slope recovered from a span
// is not bit-identical -- true of a span, and not of what set_data was handed,
// which the interpolant keeps unchanged.
template <typename S = double>
class ResourceSpline {
public:
  using value_type = S;

  // Constructors. The refinement arguments no longer select anything: the knot
  // positions are the fixed fractions times height_max. They are still taken
  // because the R constructor and the model environments pass them.
  ResourceSpline() {
    setup(1e-6, 17, 16, false);
  }

  ResourceSpline(double tol, size_t nbase, size_t max_depth, bool rescale_usually) {
    setup(tol, nbase, max_depth, rescale_usually);
  }

  void setup(double, size_t, size_t, bool) {
    // Uniform, and fixed for the run: every build places its knots at
    // u_k * height_max, so the positions and the count depend on height_max and
    // on nothing else in the state. 1/64 is exact, so u_k is too.
    knot_fractions_ = util::seq_len(0.0, 1.0, knot_count_);

    // A field to answer queries with until the first build.
    set_fixed_value(S(1.0), S(1.0));
  };

  // f returns the field's value and its vertical derivative at a height.
  template <typename Function>
  void compute_environment(Function f_value_and_slope, S height_max, bool) {
    rebuild_spline(f_value_and_slope, height_max);
  };

  void set_fixed_value(S value, S height_max) {
    // Knot positions are the interpolant's grid and stay double, so a position
    // built from an active height_max is read at its value.
    const double top = odelia::util::to_passive(height_max);
    std::vector<double> x = {0, top/2.0, top};
    std::vector<S> y = {value, value, value};
    std::vector<S> m = {S(0.0), S(0.0), S(0.0)};
    spline.clear();
    spline.init(x, y, m);
  }

  // Restores the open field rather than leaving no field at all: every query
  // reads the interpolant's bounds, and an emptied interpolant has none.
  void clear() {
    set_fixed_value(S(1.0), S(1.0));
  }

  // Highest height covered by the spline; above this get_value_at_height()
  // returns the hard-coded open value (1.0). Hoist this out of hot per-point
  // loops with get_value_at_height(height, cap).
  double max_height() const { return spline.max(); }

  S get_value_at_height(S height) const {
    return get_value_at_height(height, spline.max());
  }

  // Variant taking a pre-fetched cap (= max_height()) so callers integrating
  // over many points pay the spline.max() lookup once rather than per point.
  S get_value_at_height(S height, double cap) const {
    // TODO(#385): change maximum - here hard-coded to 1.0
    // `cap` already guards the upper bound and the crown integral keeps
    // height >= 0 = spline.min(), so use the unchecked operator() rather than
    // eval() to avoid re-running check_active()/bound checks per quadrature
    // point. Same underlying tk_spline(height) call.
    //
    // Floor the result at 0 (#253): an interpolated resource availability must
    // not be negative. The clamp is a no-op for the usual positive case, so
    // values stay bit-identical there. This is the single chokepoint for
    // FF16/K93/TF24.
    return height <= cap ? std::max(S(0.0), spline(height)) : S(1.0);
  }

  // The heights, then the values, then the slopes: the columns of r_get_state()
  // laid end to end.
  virtual void r_init_interpolators(const std::vector<double>& state) {
    if (state.size() % 3 != 0) {
      util::stop("Expected a height, a value and a slope for every knot");
    }
    const size_t state_n = state.size() / 3;
    auto it = state.begin();
    std::vector<double> state_x(it, it + state_n);
    std::vector<S> state_y(it + state_n, it + 2 * state_n);
    std::vector<S> state_m(it + 2 * state_n, state.end());
    spline.init(state_x, state_y, state_m);
  }

  // Knots the run places, fixed by the fractions and not by any build.
  size_t knot_count() const { return knot_fractions_.size(); }

  // The data the field was last built from, as the interpolant was handed it.
  const std::vector<S>& knot_values() const { return spline.values(); }
  const std::vector<S>& knot_slopes() const { return spline.slopes(); }

  // Rebuild the spans from supplied data, leaving the knot positions alone.
  // set_data length-checks, so injecting into an unbuilt field throws.
  void set_knot_data(const std::vector<S>& y, const std::vector<S>& m) {
    spline.set_data(y, m);
  }

  // Resource availability as a function of size, carrying a value and a slope at
  // every knot: what a caller reads as the slope is the derivative of what it
  // reads as the value.
  odelia::interpolator::hermite_interpolator<S> spline;

  // Knot positions in units of the canopy top, u_k = x_k / height_max, uniform and
  // fixed for the run. Nothing may reassign them: a rebuild places knots at
  // u_k * height_max, and that is what makes the positions run-constant.
  std::vector<double> knot_fractions_;

  // Knot heights, values and slopes, read back out of the spline that holds them.
  Rcpp::NumericMatrix r_get_state() const
  {
    const int n = spline.is_initialised()
                    ? static_cast<int>(spline.size()) : 0;
    Rcpp::NumericMatrix ret(n, 3);
    for (int i = 0; i < n; ++i) {
      const double x = spline.knots()[static_cast<size_t>(i)];
      S value, slope;
      spline.value_and_slope(x, value, slope);
      ret(i, 0) = x;
      ret(i, 1) = odelia::util::to_passive(value);
      ret(i, 2) = odelia::util::to_passive(slope);
    }
    ret.attr("dimnames") = Rcpp::List::create(
      R_NilValue,
      Rcpp::CharacterVector::create("height", "light_availability", "slope"));
    return ret;
  }

private:

  // Refill at the held fractions. The first and last are exactly 0 and 1, so the
  // domain is exactly [0, height_max]. The positions are u_k * height_max and
  // nothing else, so they are laid out again only when the canopy top moves;
  // every build refreshes the values and slopes.
  template <typename Function>
  void rebuild_spline(Function f_value_and_slope, S height_max) {
    if (spline.size() != knot_fractions_.size() || spline.max() != height_max) {
      // The grid stays double for the reason set_fixed_value() gives: a position
      // built from an active canopy top is laid out at its value, and the
      // field's dependence on the cohorts travels in the values and slopes.
      const double top = odelia::util::to_passive(height_max);
      std::vector<double> x(knot_fractions_.size());
      for (size_t k = 0; k < x.size(); ++k) {
        x[k] = knot_fractions_[k] * top;
      }
      spline.set_nodes(x);
    }
    const std::vector<double>& x = spline.knots();
    std::vector<S> y(x.size()), m(x.size());
    for (size_t k = 0; k < x.size(); ++k) {
      const std::pair<S, S> vs = f_value_and_slope(x[k]);
      y[k] = vs.first;
      m[k] = vs.second;
    }
    spline.set_data(y, m);
  }

  // Chosen from the re-blessing tolerance: the crown-mean light shift against an
  // adaptive fit is 1.7e-03 at worst here, and halving the spacing divides it by
  // about five.
  static constexpr size_t knot_count_ = 65;

  };


} // plant namespace

#endif
