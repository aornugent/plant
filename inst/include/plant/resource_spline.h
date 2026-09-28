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

  ResourceSpline() { setup(); }

  void setup() {
    // Uniform and fixed for the run, and they ARE the interpolant's grid: the
    // field is held against u = height / height_max, so a build refreshes values
    // and slopes and never moves a knot. 1/64 is exact, so u_k is too.
    //
    // ⚠️ DO NOT PUT THE KNOTS BACK IN HEIGHT. They were at u_k * height_max,
    // laid out at to_passive(height_max) on the ground that the field's
    // dependence on the cohorts travels in the values and slopes. It does not:
    // the canopy top is the tallest cohort's height, so 64 of 65 positions moved
    // with it and the field at a fixed height moved 2.4 per cent for a 0.1 per
    // cent change in the top at k_I = 100 -- a channel no reader of this class
    // could carry, because a position was an input nowhere.
    knot_fractions_ = util::seq_len(0.0, 1.0, knot_count_);

    // A field to answer queries with until the first build.
    set_fixed_value(S(1.0), S(1.0));
  };

  // f fills the field's value and its vertical derivative at EVERY knot from one
  // call. The whole knot set rather than one knot at a time, because the reduction
  // that produces the values is linear in the nodes plus the knots and quadratic
  // only if it is asked per knot.
  template <typename Function>
  void compute_environment(Function f_all_knots, S height_max) {
    lay_out_knots();
    height_max_ = height_max;
    // The reduction is asked at the knots' HEIGHTS, and they are active: a knot
    // sits at u_k * height_max, so what it samples moves with the canopy top and
    // the value it returns carries that.
    const size_t n = knot_fractions_.size();
    std::vector<S> x(n), y(n), m(n);
    knot_heights_.resize(n);
    for (size_t k = 0; k < n; ++k) {
      x[k] = knot_fractions_[k] * height_max;
      knot_heights_[k] = odelia::util::to_passive(x[k]);
    }
    f_all_knots(x, y, m);
    // The reduction's slope is per unit height and the interpolant's abscissa is
    // u, so the chain factor is the canopy top itself.
    for (size_t k = 0; k < n; ++k) { m[k] = m[k] * height_max; }
    spline.set_data(y, m);
  };

  void set_fixed_value(S value, S height_max) {
    // In u, so the grid is the same three points whatever the canopy top is; the
    // top is held beside it and is what a query divides by.
    height_max_ = height_max;
    const double top = odelia::util::to_passive(height_max);
    knot_heights_ = {0.0, top / 2.0, top};
    std::vector<double> x = {0.0, 0.5, 1.0};
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
  double max_height() const { return odelia::util::to_passive(height_max_); }

  // The canopy top the field is indexed against. A reader that queries the field
  // needs it, because the abscissa is height / height_max and that division is
  // where the field's dependence on the canopy top lives.
  S height_max() const { return height_max_; }

  S get_value_at_height(S height) const {
    return get_value_at_height(height, max_height());
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
    return height <= cap ? std::max(S(0.0), spline(height / height_max_))
                         : S(1.0);
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
    // Heights in, u held: the caller's grid names the canopy top as its last
    // knot, and the slopes it carries are per unit height.
    const double top = state_x.empty() ? 1.0 : state_x.back();
    height_max_ = S(top);
    knot_heights_ = state_x;
    for (double& xk : state_x) { xk /= top; }
    for (S& mk : state_m) { mk = mk * top; }
    spline.init(state_x, state_y, state_m);
  }

  // Knots the run places, fixed by the fractions and not by any build.
  size_t knot_count() const { return knot_fractions_.size(); }

  // The knot values, the knot slopes, then the canopy top: the field as the
  // interpolant holds it, which set_knot_data() restores to every digit. A record
  // in heights does not, because u_k * top / top need not be u_k.
  //
  // ⚠️ THE FIELD AND NOT THE BUILDER. What an interpolant costs is dominated by
  // the builder and band-solve workspace it drags, which no replay reads: copying
  // those whole is what put an invasion at 6.8 GB and OOM past ~10 yr.
  size_t knot_data_size() const { return 2 * knot_count() + 1; }
  template <typename It> It knot_data(It it) const {
    // A field set from heights may hold another count than knot_data_size() says.
    util::check_length(spline.size(), knot_count());
    for (const S& v : spline.values()) { util::write_iterator_scalar(it, v); }
    for (const S& v : spline.slopes()) { util::write_iterator_scalar(it, v); }
    util::write_iterator_scalar(it, height_max_);
    return it;
  }
  template <typename It> It set_knot_data(It it) {
    const size_t n = knot_count();
    std::vector<S> y(n), m(n);
    for (S& v : y) { v = *it++; }
    for (S& v : m) { v = *it++; }
    lay_out_knots();
    spline.set_data(y, m);
    height_max_ = *it++;
    const double top = odelia::util::to_passive(height_max_);
    knot_heights_.resize(n);
    for (size_t k = 0; k < n; ++k) { knot_heights_[k] = knot_fractions_[k] * top; }
    return it;
  }

  // Resource availability as a function of size, carrying a value and a slope at
  // every knot: what a caller reads as the slope is the derivative of what it
  // reads as the value.
  odelia::interpolator::hermite_interpolator<S> spline;

  // The canopy top the knots are fractions of. Active, and walked with the
  // spline: a query divides by it, so it carries the field's dependence on the
  // tallest cohort and a walk that skipped it would drop that channel again.
  S height_max_ = S(1.0);

  // The heights the field was built at, kept for reporting so that a round trip
  // through r_init_interpolators() returns what it was handed. Passive and read
  // by nothing on the rate path: u_k * height_max is what the model uses.
  std::vector<double> knot_heights_ = {0.0, 0.5, 1.0};

  template <class F>
  void for_each_active(F&& f) {
    odelia::ode::visit_active(f, spline, height_max_);
  }

  // The interpolant's abscissae, u_k = height / height_max, uniform and fixed for
  // the run. Nothing may reassign them: they are the grid itself, which is what
  // keeps a knot from carrying a derivative nobody reads.
  std::vector<double> knot_fractions_;

  // Knot heights, values and slopes, read back out of the spline that holds them.
  Rcpp::NumericMatrix r_get_state() const
  {
    const int n = spline.is_initialised()
                    ? static_cast<int>(spline.size()) : 0;
    Rcpp::NumericMatrix ret(n, 3);
    const double top = odelia::util::to_passive(height_max_);
    for (int i = 0; i < n; ++i) {
      const double u = spline.knots()[static_cast<size_t>(i)];
      S value, slope;
      spline.value_and_slope(u, value, slope);
      ret(i, 0) = knot_heights_[static_cast<size_t>(i)];
      ret(i, 1) = odelia::util::to_passive(value);
      ret(i, 2) = odelia::util::to_passive(slope) / top;
    }
    ret.attr("dimnames") = Rcpp::List::create(
      R_NilValue,
      Rcpp::CharacterVector::create("height", "light_availability", "slope"));
    return ret;
  }

private:

  // Place the knots at the held fractions, once. The first and last are exactly
  // 0 and 1, so the domain is exactly [0, 1] in u and [0, height_max] in height.
  // Nothing moves them afterwards; every build refreshes the values and slopes.
  void lay_out_knots() {
    if (spline.size() == knot_fractions_.size()) {
      return;
    }
    spline.set_nodes(knot_fractions_);
  }

  // Chosen from the re-blessing tolerance: the crown-mean light shift against an
  // adaptive fit is 1.7e-03 at worst here, and halving the spacing divides it by
  // about five.
  static constexpr size_t knot_count_ = 65;

  };


// Beer's law on a cumulative extinction profile: a resource left at a height is
// E = exp(-A) where A is the amount intercepted above it, so dE/dz = -A' exp(-A).
// One place, because the arithmetic is the same whatever else an environment holds
// and it was written out once per environment.
//
// The lambda declares its return type: an active E returned through a deduced one
// hands the field an expression template referencing operands that die here.
template <typename Field, typename Function, typename S>
void build_extinction_field(Field& field, Function f_competition_all_knots,
                            S height_max) {
  using value_type = typename Field::value_type;
  field.compute_environment(
    [&](const std::vector<value_type>& x, std::vector<value_type>& y,
        std::vector<value_type>& m) -> void {
      f_competition_all_knots(x, y, m);
      for (std::size_t k = 0; k < x.size(); ++k) {
        const value_type E = exp(-y[k]);
        m[k] = -(m[k] * E);
        y[k] = E;
      }
    },
    height_max);
}

} // plant namespace

#endif
