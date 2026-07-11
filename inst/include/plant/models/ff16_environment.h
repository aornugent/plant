// -*-c++-*-
#ifndef PLANT_PLANT_FF16_ENVIRONMENT_H_
#define PLANT_PLANT_FF16_ENVIRONMENT_H_

#include <plant/environment.h>
#include <plant/resource_spline.h>
#include <plant/ad_value.h>
#include <odelia/interpolator.hpp>
#include <plant/canopy_shape.h> // ShadingModel, shading_model_from_string
#include <cmath>                // std::floor (PPA stepping)

using namespace Rcpp;

namespace plant {

// Templated on the scalar S of the light values (knot positions stay double). At
// S = double this is the resident numerics (bit-identical, the RcppR6 alias
// below). At an active S the light a plant reads carries the derivative of the
// competition that shaded it -- the self-shading channel the resident gradient
// needs; the invasion path lifts a frozen (passive) environment through the same
// template so a mutant reads a constant canopy.
template <class S = double>
class FF16_Environment_ : public Environment {
public:
  // constructor for R interface - default settings can be modified
  // except for light_availability_spline_rescale_usually
  // which are only updated on construction
  FF16_Environment_() {
    time = 0.0;

    // Shading defaults have lower tolerance which are overwritten for speed
    light_availability = ResourceSpline_<S>(
        1e-4, // light_availability_spline_tol,
        17,   // light_availability_spline_nbase,
        16,   // light_availability_spline_max_depth,
        true  // light_availability_spline_rescale_usually)
    );

  };

  // A ResourceSpline used for storing light availbility (0-1)
  ResourceSpline_<S> light_availability;

  // PPA: when true, the light a plant experiences is the stepped (layered)
  // profile rather than the smooth one stored in light_availability. The
  // underlying spline is still fitted to the smooth optical depth (which
  // refines cleanly); the discretisation is applied at read time below.
  bool light_profile_stepped = false;
  // Thickness of one canopy layer in optical-depth units (Control::ppa_layer_optical_depth).
  double layer_optical_depth = 0.5;
  // Smoothing fraction of each layer boundary (Control::ppa_layer_smoothing).
  double layer_smoothing = 0.3;

  // Called once from the Patch constructor. Selects the stepped profile for PPA;
  // deep-crown and crown-centre keep the smooth profile.
  void set_shading_model(const std::string& model,
                         double layer_optical_depth_,
                         double layer_smoothing_) override {
    // Only PPA builds a stepped profile; every other model (including the ""
    // default and mean-light) keeps the smooth profile. Compare the string
    // directly so the "" default does not hit the throwing parser.
    light_profile_stepped = (model == "ppa");
    layer_optical_depth = layer_optical_depth_;
    layer_smoothing = layer_smoothing_;
  }

  // Ability to prescribe a fixed value
  // TODO(#476): add setting to set other variables like water
  void set_fixed_environment(double value, double height_max) {
    light_availability.set_fixed_value(value, height_max);
  }

  void set_fixed_environment(double value) {
    double height_max = 150.0;
    set_fixed_environment(value, height_max);
  }

  S get_environment_at_height(double height) const {
    return step_light(light_availability.get_value_at_height(height));
  }

  // Highest height covered by the light spline; hoist out of hot per-point
  // loops and feed back into the capped get_environment_at_height() overload.
  double max_environment_height() const {
    return light_availability.max_height();
  }

  S get_environment_at_height(double height, double cap) const {
    return step_light(light_availability.get_value_at_height(height, cap));
  }

  // Slope of the frozen light profile at height as a passive double, for the
  // crown integral's active-height linearisation. step_light is the identity for
  // the smooth (default) shading models, so this is the spline slope; the PPA
  // staircase's extra chain-rule factor is out of that path.
  double get_environment_deriv_at_height(double height, double cap) const {
    return light_availability.get_deriv_at_height(height, cap);
  }

  // Knot positions of the current light spline (always double). Recorded per
  // accepted step so the resident replay rebuilds the canopy on them.
  std::vector<double> light_knots() const {
    return light_availability.get_knots();
  }

  // Discretise a smooth light value into PPA canopy layers. For the smooth
  // models this is a single predicted branch returning the input unchanged, so
  // it adds no measurable cost to deep-crown/crown-centre. For PPA it maps the
  // optical depth tau = -log(E) onto a smoothed integer number of layers of
  // thickness layer_optical_depth and back-transforms:
  //   E_step = exp(-d * smooth_floor(tau / d)).
  S step_light(S E) const {
    if (!light_profile_stepped || E >= 1.0) {
      return E;
    }
    // Guard the log: the smooth spline can undershoot to <= 0, which would make
    // tau non-finite. Such a point is fully shaded, so return 0.
    if (E <= 0.0) {
      return S(0.0);
    }
    const S tau = -log(E);
    return exp(-layer_optical_depth * smooth_floor(tau / layer_optical_depth));
  }

  // Monotone, C1-continuous smooth staircase. Each layer is flat over its lower
  // (1 - layer_smoothing) and ramps to the next integer via a cubic smoothstep
  // over its top layer_smoothing fraction. C1 at the joins because smoothstep
  // has zero slope at both ends, so the resulting light profile is smooth enough
  // for the adaptive ODE solver. With layer_smoothing -> 0 this recovers the
  // hard floor (and its instability). floor is a passive integer boundary (zero
  // derivative); the trait derivative flows through the smoothstep ramp.
  S smooth_floor(S u) const {
    const S n = S(std::floor(ad_value(u)));
    const double w = layer_smoothing;
    if (w <= 0.0) {
      return n; // hard step
    }
    const S f = u - n; // fractional position within the layer, [0, 1)
    if (f <= 1.0 - w) {
      return n; // flat lower part of the layer
    }
    const S t = (f - (1.0 - w)) / w; // [0, 1] across the transition
    return n + t * t * (3.0 - 2.0 * t);    // cubic smoothstep
  }

  virtual void r_init_interpolators(const std::vector<double> &state)
  {
    light_availability.r_init_interpolators(state);
  }

  virtual void compute_rates(std::vector<double> const& resource_depletion) {

  }

  virtual Rcpp::List r_get_state() const {
    return Rcpp::List::create(
              _["light_availability"] = light_availability.r_get_state()
            );
  }

  // Pre-compute resources available in the environment, as a function of height
  template <typename Function>
  void compute_environment(Function f_compute_competition, double height_max, bool rescale) {

    // Define an anonymous function to use in creation of light_availability spline
    // Note: extinction coefficient was already applied in strategy, so
    // f_compute_competition gives sum of projected leaf area (k L) across species. Just need to apply Beer's law, E = exp(- (k L))
    auto f_light_availability = [&](double height) -> S
    { return exp(-f_compute_competition(height)); };

    // Calculates the light_availability spline, by fitting to the function
    // `f_compute_competition` as a function of height
    light_availability.compute_environment(f_light_availability, height_max, rescale);
  }

  // Rebuild the canopy on frozen knot positions with the (active) competition:
  // the resident replay's L2 recompute. The positions the double pass chose stay
  // fixed; only the values re-shade with the trait, so the self-shading
  // derivative flows through area_leaf without moving any node.
  template <typename Function>
  void compute_environment_fixed(Function f_compute_competition,
                                 const std::vector<double>& knots) {
    auto f_light_availability = [&](double height) -> S
    { return exp(-f_compute_competition(height)); };
    light_availability.compute_environment_fixed(f_light_availability, knots);
  }

  virtual void clear_environment() {
    light_availability.clear();
  }

  // double -> active mould (so the Patch/SCM rebind lifts the environment too),
  // and a config-only copy onto S2 carrying the recorded canopy across as passive
  // constants for the invasion (frozen-canopy) replay.
  template <class S2> using rebind = FF16_Environment_<S2>;

  template <class S2>
  FF16_Environment_<S2> rebind_from() const {
    FF16_Environment_<S2> out;
    out.time = time;
    out.vars = vars;
    out.extrinsic_drivers = extrinsic_drivers;
    out.species_arriving_index = species_arriving_index;
    out.light_profile_stepped = light_profile_stepped;
    out.layer_optical_depth = layer_optical_depth;
    out.layer_smoothing = layer_smoothing;
    out.light_availability = light_availability.template rebind_from<S2>();
    return out;
  }
};

// S = double crosses the R boundary (the RcppR6 name_cpp is this alias) and is
// the invasion (frozen-canopy) environment; the active twin is instantiated only
// where a resident gradient recomputes the canopy.
using FF16_Environment = FF16_Environment_<double>;


}

#endif
