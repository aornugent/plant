// -*-c++-*-

#ifndef PLANT_PLANT_RESOURCE_SPLINE_H_
#define PLANT_PLANT_RESOURCE_SPLINE_H_

#include <odelia/hermite_interpolator.hpp>
#include <odelia/ode_util.hpp>
#include <odelia/ode_interface.hpp>
#include <plant/util.h>
#include <algorithm> // std::max, for the resource-availability floor (#253)
#include <cmath>     // std::ceil, for how far the grid has to reach
#include <utility>   // std::pair, the value and slope a build supplies

using namespace Rcpp;

namespace plant {

// Templated on the scalar S the resource values carry; the knot positions stay
// double. S = double is production.
template <typename S = double>
class ResourceSpline {
public:
  using value_type = S;

  // The spacing is the only thing a caller chooses. Knots sit at k * spacing,
  // which is a set of constants: no cohort height reaches a position, so the
  // chain from the canopy top into every knot position does not exist.
  //
  // No default. A spacing sized to the wrong model is silent -- K93's stand is
  // a quarter the height of the others and reads 1.8e-03 off the refined answer
  // at theirs, with nothing raised -- so every environment states its own.
  explicit ResourceSpline(double knot_spacing) {
    setup(knot_spacing);
  }

  void setup(double knot_spacing) {
    if (!(knot_spacing > 0.0)) {
      util::stop("ResourceSpline: knot spacing must be positive");
    }
    knot_spacing_ = knot_spacing;

    // A field to answer queries with until the first build.
    set_fixed_value(S(1.0), S(1.0));
  };

  // f fills the field's value and its vertical derivative at every knot of the
  // grid it is handed. Whole-grid rather than per height, because the reduction
  // behind it costs cohorts once for the grid and cohorts per height otherwise.
  template <typename Function>
  void compute_environment(Function f_value_and_slope, S height_max) {
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
    knot_values_ = y;
    knot_slopes_ = m;
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
    knot_values_ = state_y;
    knot_slopes_ = state_m;
  }

  // Read off the grid rather than stored beside it. set_cohort_reads() fills a
  // vector of this length and hands it to set_data(), which length-checks
  // against the grid, so a stored count that drifted from the grid would be a
  // buffer overrun rather than an error.
  size_t knot_count() const { return spline.size(); }

  // The data the field was last built from, held rather than read back out of
  // the interpolant: a knot slope recovered from a span is not bit-identical.
  const std::vector<S>& knot_values() const { return knot_values_; }
  const std::vector<S>& knot_slopes() const { return knot_slopes_; }

  // Rebuild the spans from supplied data, leaving the knot positions alone.
  // set_data length-checks, so injecting into an unbuilt field throws.
  void set_knot_data(const std::vector<S>& y, const std::vector<S>& m) {
    spline.set_data(y, m);
    knot_values_ = y;
    knot_slopes_ = m;
  }

  // Resource availability as a function of size, carrying a value and a slope at
  // every knot: what a caller reads as the slope is the derivative of what it
  // reads as the value.
  odelia::interpolator::hermite_interpolator<S> spline;

  // Mirrors of the interpolant's data, kept as the pack's source of truth.
  std::vector<S> knot_values_;
  std::vector<S> knot_slopes_;

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

  // Reach the canopy, then refill. The canopy decides how many knots there are
  // and nothing else about them: knot k sits at k * spacing whatever the stand
  // is doing, so no cohort height reaches a position and there is no position
  // derivative to drop.
  //
  // The lattice only ever grows. Above the canopy the profile is exactly flat --
  // the crown shape's value and slope both vanish at a cohort's own top -- so
  // the knots an extension adds carry (1, 0) with a derivative that is exactly
  // zero, and no query reaches the span they open. Extending is therefore
  // bit-identical for every query at or below the canopy, which is what lets
  // the count depend on the state without the count carrying a derivative.
  template <typename Function>
  void rebuild_spline(Function f_value_and_slope, S height_max) {
    const double top = odelia::util::to_passive(height_max);
    // The canopy sizes the grid, where a grid spanning the canopy's own height
    // did not, so a runaway size-density equation now reaches the allocation
    // rather than only making the field coarse. Refuse it here, naming the
    // height, rather than reaching for the memory: nothing in this model grows a
    // canopy of that size, so this is a runaway and not a forest.
    if (!(top >= 0.0) || top > knot_spacing_ * max_knots_) {
      util::stop("ResourceSpline: canopy height " + util::format_double(top) +
                 " m needs more than " + util::to_string(max_knots_) +
                 " knots at a spacing of " + util::format_double(knot_spacing_) +
                 " m; the size-density equation has run away");
    }
    // One knot clear of the canopy, so a query at exactly height_max is inside
    // the grid rather than on its last node.
    const size_t wanted =
      static_cast<size_t>(std::ceil(top / knot_spacing_)) + 2;
    // Long enough is not the whole test: set_fixed_value() and a restored state
    // both leave a grid this class did not lay out, and one of those can be long
    // enough while sitting somewhere else entirely. Read whether the grid is the
    // lattice off the grid rather than remembering it, so the two cannot drift.
    const size_t held = spline.size();
    const bool on_lattice =
      held >= 2 && spline.knots()[1] == knot_spacing_ &&
      spline.max() == static_cast<double>(held - 1) * knot_spacing_;
    if (!on_lattice || held < wanted) {
      std::vector<double> x(wanted);
      for (size_t k = 0; k < wanted; ++k) {
        x[k] = static_cast<double>(k) * knot_spacing_;
      }
      spline.set_nodes(x);
    }
    const std::vector<double>& x = spline.knots();
    std::vector<S> y(x.size()), m(x.size());
    f_value_and_slope(x, y, m);
    spline.set_data(y, m);
    knot_values_ = y;
    knot_slopes_ = m;
  }

  // Metres between knots. Chosen against the crown-mean light a cohort reads:
  // 0.1 holds that to about 1e-6 at every canopy height from a seedling to 35 m,
  // where knots tied to the canopy top run from 7e-8 to 6e-4 over the same range.
  double knot_spacing_ = 0.1;

  // Far above any canopy this model grows -- 2.5 km at the coarsest spacing in
  // use, against a tallest-tree record near 130 m -- so it bounds a runaway
  // without being reachable by a stand.
  static constexpr size_t max_knots_ = 25000;

  };


} // plant namespace

#endif
