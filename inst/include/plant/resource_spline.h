// -*-c++-*-

#ifndef PLANT_PLANT_RESOURCE_SPLINE_H_
#define PLANT_PLANT_RESOURCE_SPLINE_H_

#include <odelia/interpolator.hpp>
#include <plant/adaptive_interpolator.h>
#include <odelia/ode_interface.hpp>
#include <plant/ad_value.h>
#include <plant/util.h>
#include <algorithm> // std::max, for the resource-availability floor (#253)
#include <type_traits>

using namespace Rcpp;

namespace plant {

// Templated on the scalar S of the light VALUES; knot positions stay double. At
// S = double this is the resident numerics (bit-identical, the RcppR6 alias
// below). At an active S the interpolated light differentiates w.r.t. the
// competition that shaped it -- the self-shading channel of the resident
// gradient. The adaptive build stays double (plant's AdaptiveInterpolator); an
// active spline is only ever built on frozen knots (compute_environment_fixed).
template <typename S = double>
class ResourceSpline_ {
  // Sibling instantiations lift each other's private build parameters.
  template <typename> friend class ResourceSpline_;
public:

  // Constructors
  ResourceSpline_() {
    setup(1e-6, 17, 16, false);
  }

  ResourceSpline_(double tol, size_t nbase, size_t max_depth, bool rescale_usually) {
    setup(tol, nbase, max_depth, rescale_usually);
  }

  void setup(double tol, size_t nbase, size_t max_depth, bool rescale_usually) {
    tol_ = tol;
    nbase_ = nbase;
    max_depth_ = max_depth;

    // Initialise adaptive interpolator. This object can create an interpolator spline
    spline_construction =
        interpolator::AdaptiveInterpolator(tol, tol, nbase, max_depth);
    // Create an actual spline, For initalisation
    // Provide a dummy function and construct
    // This will be over-written later with actual function
    construct_spline([&](double height) { return get_value_at_height(height); }, 1.0);

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

  // Rebuild the light spline on a set of frozen knot positions, recomputing the
  // values from the (active) competition. The resident replay's L2 path: the
  // positions the double pass chose stay fixed, the values re-shade with the
  // trait so the self-shading derivative flows through.
  template <typename Function>
  void compute_environment_fixed(Function f_compute_competition,
                                 const std::vector<double>& knots) {
    std::vector<S> y;
    y.reserve(knots.size());
    for (double xi : knots) {
      y.push_back(f_compute_competition(xi));
    }
    spline.clear();
    spline.init(knots, y);
  }

  void set_fixed_value(double value, double height_max) {
    std::vector<double> x = {0, height_max/2.0, height_max};
    std::vector<S> y = {S(value), S(value), S(value)};
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

  // Knot positions of the current spline (always double). Recorded per accepted
  // step so the resident replay rebuilds on them.
  std::vector<double> get_knots() const { return spline.get_x(); }

  S get_value_at_height(double height) const {
    return get_value_at_height(height, spline.max());
  }

  // Variant taking a pre-fetched cap (= max_height()) so callers integrating
  // over many points pay the spline.max() lookup once rather than per point.
  S get_value_at_height(double height, double cap) const {
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

  // Analytic slope of the (frozen) spline as a passive double, matching
  // get_value_at_height's clamps: zero above the cap or where the value is
  // floored to zero, else the interpolator's exact derivative. Feeds the crown
  // integral's height-linearisation on the gradient path; kept passive so only
  // the light value (not its slope) carries the self-shading derivative.
  double get_deriv_at_height(double height, double cap) const {
    if (height > cap) {
      return 0.0;
    }
    return ad_value(std::max(S(0.0), spline(height))) > 0.0
               ? ad_value(spline.deriv(height))
               : 0.0;
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
    std::copy_n(it, state_n, std::back_inserter(state_x));
    for (auto yi = it + state_n; yi != it + 2 * state_n; ++yi) {
      state_y.push_back(S(*yi));
    }
    spline.init(state_x, state_y);
  }

  // This object will store an interpolator spline of
  // resource availability as a function of size
  odelia::interpolator::basic_interpolator<S> spline;

  // This object can create an interpolator spline via adaptive refinement
  // (double only; an active spline is built on frozen knots, not refined).
  interpolator::AdaptiveInterpolator spline_construction;

  // flag, do we try to rescale the spline when possible? this is quicker
  bool spline_rescale_usually;

  // Config-only lift onto another scalar S2: positions carry across as-is, the
  // recorded light values become passive constants (via ad_value). Used to lift
  // a frozen resident environment for the invasion (frozen-canopy) replay.
  template <class S2>
  ResourceSpline_<S2> rebind_from() const {
    ResourceSpline_<S2> out;
    out.tol_ = tol_;
    out.nbase_ = nbase_;
    out.max_depth_ = max_depth_;
    out.spline_construction = spline_construction;
    out.spline_rescale_usually = spline_rescale_usually;
    std::vector<double> x = spline.get_x();
    if (x.size() >= 3) {
      std::vector<S> y = spline.get_y();
      std::vector<S2> y2;
      y2.reserve(y.size());
      for (const auto& v : y) {
        y2.push_back(S2(ad_value(v)));
      }
      out.spline.init(x, y2);
    }
    return out;
  }

  Rcpp::NumericMatrix r_get_state() const
  {

    // format spline as Matrix
    std::vector<std::vector<double>> xy;
    xy.push_back(spline.get_x());
    std::vector<double> yy;
    for (const auto& v : spline.get_y()) {
      yy.push_back(ad_value(v));
    }
    xy.push_back(yy);

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

  double tol_;
  size_t nbase_, max_depth_;

  template <typename Function>
  void construct_spline(Function f_compute_competition, double height_max)
  {
    const double lower_bound = 0.0;
    double upper_bound = height_max;

    if constexpr (std::is_same_v<S, double>) {
      spline =
        spline_construction.construct(f_compute_competition, lower_bound, upper_bound);
    } else {
      // Active spline never crosses the adaptive path in a resident replay (it is
      // built on frozen knots). This branch only fires for the empty-patch build
      // at reset; odelia's construct is scalar-templated, so it compiles and runs
      // there without dragging the double-only AdaptiveInterpolator to S.
      spline.construct(f_compute_competition, lower_bound, upper_bound,
                       tol_, tol_, nbase_, max_depth_);
    }
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

// S = double crosses the R boundary (the RcppR6 name_cpp is this alias) and is
// the type K93/TF24 environments hold; every existing caller is unchanged.
using ResourceSpline = ResourceSpline_<double>;

} // plant namespace

#endif
