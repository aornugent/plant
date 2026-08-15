// -*-c++-*-

#ifndef PLANT_PLANT_RESOURCE_SPLINE_H_
#define PLANT_PLANT_RESOURCE_SPLINE_H_

#include <odelia/hermite_interpolator.hpp>
#include <odelia/ode_interface.hpp>
#include <plant/util.h>
#include <algorithm> // std::max, for the resource-availability floor (#253)
#include <utility>   // std::pair, the value and slope a build supplies
#include <vector>

using namespace Rcpp;

namespace plant {

class ResourceSpline {
public:

  // Knots sit at k * spacing, so no cohort height reaches a position.
  //
  // No default: a spacing sized to the wrong model is silent. K93's stand is a
  // quarter the height of the others and reads 1.8e-03 off the refined answer at
  // theirs, with nothing raised, so every environment states its own.
  explicit ResourceSpline(double knot_spacing) {
    setup(knot_spacing);
  }

  void setup(double knot_spacing) {
    if (!(knot_spacing > 0.0)) {
      util::stop("ResourceSpline: knot spacing must be positive");
    }
    knot_spacing_ = knot_spacing;
    set_fixed_value(1.0, 1.0);
  };

  // f fills the value and dA/dz at every knot of the grid it is handed. Whole
  // grid, because the reduction behind it costs cohorts once for a grid and
  // cohorts per height otherwise.
  template <typename Function>
  void compute_environment(Function f_value_and_slope, double height_max) {
    rebuild_spline(f_value_and_slope, height_max);
  };

  void set_fixed_value(double value, double height_max) {
    std::vector<double> x = {0, height_max / 2.0, height_max};
    std::vector<double> y = {value, value, value};
    std::vector<double> m = {0.0, 0.0, 0.0};
    spline.clear();
    spline.init(x, y, m);
  }

  // An open field, which is what a patch with nothing in it casts. Every query
  // reads the interpolant's bounds, so leaving it empty would answer none.
  void clear() {
    set_fixed_value(1.0, 1.0);
  }

  // Highest height covered by the spline; above this get_value_at_height()
  // returns the hard-coded open value (1.0). Hoist this out of hot per-point
  // loops with get_value_at_height(height, cap).
  double max_height() const { return spline.max(); }

  double get_value_at_height(double height) const {
    return get_value_at_height(height, spline.max());
  }

  // Variant taking a pre-fetched cap (= max_height()) so callers integrating
  // over many points pay the spline.max() lookup once rather than per point.
  double get_value_at_height(double height, double cap) const {
    // TODO(#385): change maximum - here hard-coded to 1.0
    //
    // Floor the result at 0 (#253): an interpolated resource availability must
    // not be negative. The clamp is a no-op for the usual positive case, so
    // values stay bit-identical there. This is the single chokepoint for
    // FF16/K93/TF24.
    return height <= cap ? std::max(0.0, spline(height)) : 1.0;
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
    std::vector<double> state_y(it + state_n, it + 2 * state_n);
    std::vector<double> state_m(it + 2 * state_n, state.end());
    spline.init(state_x, state_y, state_m);
  }

  // Read off the grid rather than stored beside it, so a count cannot drift from
  // the vector a build fills.
  size_t knot_count() const { return spline.size(); }

  // Resource availability as a function of size, carrying a value and a slope at
  // every knot: what a caller reads as the slope is the derivative of what it
  // reads as the value.
  odelia::interpolator::hermite_interpolator<double> spline;

  // Knot heights, values and slopes, as the field was built from them.
  Rcpp::NumericMatrix r_get_state() const
  {
    const int n = static_cast<int>(spline.values().size());
    Rcpp::NumericMatrix ret(n, 3);
    for (int i = 0; i < n; ++i) {
      const size_t k = static_cast<size_t>(i);
      ret(i, 0) = spline.knots()[k];
      ret(i, 1) = spline.values()[k];
      ret(i, 2) = spline.slopes()[k];
    }
    ret.attr("dimnames") = Rcpp::List::create(
      R_NilValue,
      Rcpp::CharacterVector::create("height", "light_availability", "slope"));
    return ret;
  }

private:

  // Reach the canopy, then refill. The canopy sets how many knots there are and
  // nothing else about them, and the lattice only ever grows: an extension adds
  // knots above every span the old grid covered, so every query at or below the
  // canopy reads the same number after a stand has grown.
  template <typename Function>
  void rebuild_spline(Function f_value_and_slope, double height_max) {
    // The interpolant lays whatever lattice it is asked for; only this class can
    // say what ran away, so the refusal names the height.
    if (!(height_max >= 0.0) || height_max > knot_spacing_ * max_knots_) {
      util::stop("ResourceSpline: canopy height " +
                 util::format_double(height_max) + " m needs more than " +
                 util::to_string(max_knots_) + " knots at a spacing of " +
                 util::format_double(knot_spacing_) +
                 " m; the size-density equation has run away");
    }
    spline.ensure_lattice(knot_spacing_,
                          spline.lattice_size(knot_spacing_, height_max));
    const std::vector<double>& x = spline.knots();
    std::vector<double> y(x.size()), m(x.size());
    f_value_and_slope(x, y, m);
    spline.set_data(y, m);
  }

  // Metres between knots. Every environment states its own; see setup().
  double knot_spacing_;

  // 1250 m of canopy at the 0.05 m spacing FF16 and TF24 use, against a
  // tallest-tree record near 130 m: it bounds a runaway without a stand
  // reaching it.
  static constexpr size_t max_knots_ = 25000;

  };


} // plant namespace

#endif
