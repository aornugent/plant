// -*-c++-*-
#ifndef PLANT_PLANT_FF16_ENVIRONMENT_H_
#define PLANT_PLANT_FF16_ENVIRONMENT_H_

#include <plant/environment.h>
#include <plant/resource_spline.h>
#include <odelia/interpolator.hpp>
#include <plant/canopy_shape.h> // ShadingModel, shading_model_from_string
#include <cmath>                // std::log, std::exp, std::floor (PPA stepping)

using namespace Rcpp;

namespace plant {

// FF16's environment: a single light field, no ODE state (ode_size() == 0).
// Templated on the scalar S of the light knot values; S = double is the
// production path (the `FF16_Environment` alias below). Base members are reached
// through this-> because Environment_<S> is a dependent base.
template <class S = double>
class FF16_Environment_ : public Environment_<S> {
public:
  // Same environment at a different scalar U (nested forward-over-reverse for
  // dg/dh); paired with FF16_Strategy_'s rebind so generic code can lift the
  // whole rate path to a forward tangent type.
  template <class U> using rebind = FF16_Environment_<U>;

  // constructor for R interface - default settings can be modified
  // except for light_availability_spline_rescale_usually
  // which are only updated on construction
  FF16_Environment_() {
    this->time = 0.0;

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

  // Fix the light field to a scalar value that may carry AD derivatives (the
  // resident competition a cohort sees). Unlike set_fixed_environment(double),
  // which strips those derivatives, this keeps them. Used to build the frozen
  // environment for the forward-over-reverse dg/dh evaluation: the query is
  // frozen, so a field fixed at the cohort's current competition reproduces the
  // read while preserving d(competition)/dtheta.
  void set_fixed_environment_scalar(S value, double height_max) {
    light_availability.set_fixed_value_scalar(value, height_max);
  }

  // d(light)/d(height) by secant -- the coupling channel for the density-transport
  // dg/dh. The field owns its slope; the dg/dh seam reads this rather than
  // hand-rolling a secant, with step+direction from the same Control the
  // production stencil uses.
  S get_environment_slope_at_height(S height, double step, int direction) const {
    return light_availability.slope_at_height(height, step, direction);
  }

  S get_environment_at_height(S height) const {
    return step_light(light_availability.get_value_at_height(height));
  }

  // Highest height covered by the light spline; hoist out of hot per-point
  // loops and feed back into the capped get_environment_at_height() overload.
  double max_environment_height() const {
    return light_availability.max_height();
  }

  S get_environment_at_height(S height, S cap) const {
    return step_light(light_availability.get_value_at_height(height, cap));
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
      return static_cast<S>(0.0);
    }
    const S tau = -log(E);
    return exp(-layer_optical_depth * smooth_floor(tau / layer_optical_depth));
  }

  // Monotone, C1-continuous smooth staircase. Each layer is flat over its lower
  // (1 - layer_smoothing) and ramps to the next integer via a cubic smoothstep
  // over its top layer_smoothing fraction. C1 at the joins because smoothstep
  // has zero slope at both ends, so the resulting light profile is smooth enough
  // for the adaptive ODE solver. With layer_smoothing -> 0 this recovers the
  // hard floor (and its instability). The integer floor is a piecewise-constant
  // selector (derivative zero), so it is taken on the double value.
  S smooth_floor(S u) const {
    const double n = std::floor(odelia::util::to_passive(u));
    const double w = layer_smoothing;
    if (w <= 0.0) {
      return static_cast<S>(n); // hard step
    }
    const S f = u - n; // fractional position within the layer, [0, 1)
    if (xad::value(f) <= 1.0 - w) {
      return static_cast<S>(n); // flat lower part of the layer
    }
    const S t = (f - (1.0 - w)) / w; // [0, 1] across the transition
    return n + t * t * (3.0 - 2.0 * t);    // cubic smoothstep
  }

  virtual void r_init_interpolators(const std::vector<double> &state)
  {
    light_availability.r_init_interpolators(state);
  }

  virtual void compute_rates(std::vector<S> const& resource_depletion) {

  }

  virtual Rcpp::List r_get_state() const {
    return Rcpp::List::create(
              _["light_availability"] = light_availability.r_get_state()
            );
  }

  // Pre-compute resources available in the environment, as a function of height.
  // The field CONSTRUCTION stays double (adaptive knot placement + a double
  // competition function); an active resident light field is rebuilt on the
  // recorded knots later (the deferred L2 path).
  template <typename Function>
  void compute_environment(Function f_compute_competition, double height_max, bool rescale) {

    // Define an anonymous function to use in creation of light_availability spline
    // Note: extinction coefficient was already applied in strategy, so
    // f_compute_competition gives sum of projected leaf area (k L) across species. Just need to apply Beer's law, E = exp(- (k L))
    // Build the light field at the environment's scalar S (mirrors K93): the
    // knot heights are chosen from double values, but the knot VALUES carry the
    // resident self-shading derivative through Beer's law -- the active resident
    // field (L2 recompute), not a double-only construction.
    auto f_light_availability = [&](double height) -> S
    { return exp(-f_compute_competition(height)); };

    // Calculates the light_availability spline, by fitting to the function
    // `f_compute_competition` as a function of height
    light_availability.compute_environment(f_light_availability, height_max, rescale);
  }

  virtual void clear_environment() {
    light_availability.clear();
  }
};

// The double instantiation is the production path bound by RcppR6 as
// `plant::FF16_Environment`.
using FF16_Environment = FF16_Environment_<double>;

}

#endif
