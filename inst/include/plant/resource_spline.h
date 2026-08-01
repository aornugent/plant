// -*-c++-*-

#ifndef PLANT_PLANT_RESOURCE_SPLINE_H_
#define PLANT_PLANT_RESOURCE_SPLINE_H_

#include <odelia/interpolator.hpp>
#include <odelia/ode_interface.hpp>
#include <plant/util.h>
#include <algorithm> // std::max, for the resource-availability floor (#253)

using namespace Rcpp;

namespace plant {

// Templated on the scalar S the resource values carry; the knot positions stay
// double. S = double is production.
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

  template <typename Function>
  void compute_environment(Function f_compute_competition, S height_max, bool) {
    rebuild_spline(f_compute_competition, height_max);
  };

  void set_fixed_value(S value, S height_max) {
    std::vector<double> x = {0, height_max/2.0, height_max};
    std::vector<S> y = {value, value, value};
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
    // Floor the result at 0 (#253): the cubic spline can undershoot below
    // zero between knots (notably the K93 light spline at high k_I), which is
    // non-physical for a resource availability. The clamp is a no-op for the
    // usual positive case (so values stay bit-identical there) and only bites
    // on the spurious negative undershoot. This is the single chokepoint for
    // FF16/K93/TF24, and belongs here in plant rather than in the
    // general-purpose interpolator (which is migrating to odelia).
    return height <= cap ? std::max(S(0.0), spline(height)) : S(1.0);
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

  // This object will store an interpolator spline of 
  // resource availability as a function of size
  odelia::interpolator::basic_interpolator<S> spline;

  // Knot positions in units of the canopy top, u_k = x_k / height_max, uniform and
  // fixed for the run. Nothing may reassign them: a rebuild places knots at
  // u_k * height_max, and that is what makes the positions run-constant.
  std::vector<double> knot_fractions_;

  Rcpp::NumericMatrix r_get_state() const
  {

    // format spline as Matrix
    std::vector<std::vector<double>> xy;
    xy.push_back(spline.get_x());
    xy.push_back(spline.get_y());

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

  // Refit at the held fractions. The first and last are exactly 0 and 1, so the
  // rebuilt domain is exactly [0, height_max].
  template <typename Function>
  void rebuild_spline(Function f_compute_competition, S height_max) {
    spline.clear();
    for (double u : knot_fractions_) {
      const double x = u * height_max;
      spline.add_point(x, f_compute_competition(x));
    }
    spline.initialise();
  }

  // Chosen from the re-blessing tolerance: the crown-mean light shift against an
  // adaptive fit is 1.7e-03 at worst here, and halving the spacing divides it by
  // about five.
  static constexpr size_t knot_count_ = 65;

  };


} // plant namespace

#endif
