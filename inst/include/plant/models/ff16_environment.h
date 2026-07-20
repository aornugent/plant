// -*-c++-*-
#ifndef PLANT_PLANT_FF16_ENVIRONMENT_H_
#define PLANT_PLANT_FF16_ENVIRONMENT_H_

#include <plant/environment.h>
#include <plant/resource_spline.h>
#include <odelia/interpolator.hpp>
#include <odelia/separable_field.hpp>
#include <plant/canopy_shape.h> // ShadingModel, shading_model_from_string
#include <algorithm>
#include <array>
#include <functional>
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

  // Exact separable competition field (same rank-3 Yokozawa field K93 uses),
  // assembled from the cohorts by the patch. FF16 reads its deep-crown light
  // through this so the query-height self-shading feedback flows (the derivative
  // the fitted spline freezes). The spline is still built (field_supersedes_spline
  // = false) for the not-yet-assembled / fixed-environment reads.
  static constexpr std::size_t comp_rank = CanopyShape<S>::shading_rank;  // 3
  static constexpr bool field_supersedes_spline = false;
  odelia::separable_field<S, comp_rank> competition_field;
  std::vector<double> competition_source_heights;  // descending
  CanopyShape<S> competition_canopy;
  bool competition_field_ready = false;

  // DEBUG (channel isolation): when true, optical_depth uses PASSIVE query factors
  // a_p(z), so the field carries the SOURCE self-shading derivative but not the
  // query-height one -- reproducing what the fitted spline carries. Lets a driver
  // split the field's source channel from its query channel. Off in production.
  static inline bool freeze_query_derivative = false;

  void assemble_competition_field(
      const std::array<std::vector<S>, comp_rank>& source_weight,
      const std::vector<double>& source_heights, const CanopyShape<S>& canopy) {
    competition_field.assemble(source_weight);
    competition_source_heights = source_heights;
    competition_canopy = canopy;
    competition_field_ready = true;
  }
  void clear_competition_field() { competition_field_ready = false; }

  std::size_t n_sources_at_least(double z) const {
    auto it = std::upper_bound(competition_source_heights.begin(),
                               competition_source_heights.end(), z,
                               std::greater<double>());
    return static_cast<std::size_t>(it - competition_source_heights.begin());
  }

  // DEBUG (channel isolation): when true, the whole optical depth is stripped to
  // passive, so the self-shading light feedback derivative d(light)/dtheta is zero
  // (as if the field were a recorded double). Off in production.
  static inline bool freeze_field_derivative = false;

  // Optical depth A(z) at the query; light = exp(-A). freeze_query_derivative
  // strips the query-height derivative by evaluating a_p at passive z.
  S field_optical_depth(S height) const {
    const double z = odelia::util::to_passive(height);
    const std::size_t k = n_sources_at_least(z);
    if (k == 0) return S(0.0);
    const S q = freeze_query_derivative ? S(z) : height;
    const S A = competition_field.at(competition_canopy.template shading_query_factors<S>(q), k - 1);
    return freeze_field_derivative ? S(odelia::util::to_passive(A)) : A;
  }
  double competition_max_source_height() const {
    return competition_source_heights.empty() ? 0.0 : competition_source_heights.front();
  }

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
    if (competition_field_ready)
      return step_light(exp(-field_optical_depth(height)));
    return step_light(light_availability.get_value_at_height(height));
  }

  // Highest height covered by the light field (or spline): the crown-integral
  // upper bound. Above it the light is full.
  double max_environment_height() const {
    if (competition_field_ready) return competition_max_source_height();
    return light_availability.max_height();
  }

  // The cap only clamps spline extrapolation; the exact field needs none.
  S get_environment_at_height(S height, S cap) const {
    if (competition_field_ready)
      return step_light(exp(-field_optical_depth(height)));
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
