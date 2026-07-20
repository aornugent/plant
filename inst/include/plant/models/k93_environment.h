// -*-c++-*-
#ifndef PLANT_PLANT_K93_ENVIRONMENT_H_
#define PLANT_PLANT_K93_ENVIRONMENT_H_

#include <plant/environment.h>
#include <plant/resource_spline.h>
#include <plant/canopy_shape.h>
#include <odelia/separable_field.hpp>
#include <algorithm>
#include <array>
#include <functional>
#include <vector>

using namespace Rcpp;

namespace plant {

// K93's environment: a single light field, no ODE state (ode_size() == 0).
// Templated on the scalar S of the light knot values; S = double is the
// production path (the `K93_Environment` alias below). Base members are reached
// through this-> because Environment_<S> is a dependent base.
template <class S = double>
class K93_Environment_ : public Environment_<S> {
public:
  // Same environment at a different scalar U (nested forward-over-reverse for
  // dg/dh, plant#39).
  template <class U> using rebind = K93_Environment_<U>;

  K93_Environment_() {
    this->time = 0.0;
    // Match FF16: loosen the light-availability spline tolerance from the
    // ResourceSpline default (1e-6) to 1e-4 for speed. The spline is rebuilt
    // every ODE step, so its construction dominates K93 runtime; 1e-6 was 100x
    // tighter than FF16 for no comparable accuracy need.
    light_availability = ResourceSpline_<S>(
        1e-4, // light_availability_spline_tol
        17,   // light_availability_spline_nbase
        16,   // light_availability_spline_max_depth
        true  // light_availability_spline_rescale_usually
    );
  };

  // Light interface
  ResourceSpline_<S> light_availability;

  // Exact competition field (odelia::separable_field) replacing the fitted spline
  // for the resident coupling read: assembled from the cohort population by the
  // patch each step, it makes a cohort's light read carry BOTH the source
  // self-shading AND the query-height feedback (the derivative the spline drops).
  // Active until the patch assembles it, get_environment_at_height keeps the spline
  // path, so this is behaviour-preserving on the double path until wired.
  static constexpr std::size_t comp_rank = CanopyShape<S>::shading_rank;  // 3
  // K93's exact field fully serves the light read, so the patch drops the fitted
  // spline build (FF16 keeps it -- see FF16_Environment).
  static constexpr bool field_supersedes_spline = true;
  odelia::separable_field<S, comp_rank> competition_field;
  std::vector<double> competition_source_heights;  // descending, for the rank search
  CanopyShape<S> competition_canopy;                // supplies the query factors a_p(z)
  bool competition_field_ready = false;

  // Assemble the exact field from per-cohort source weights (descending height).
  // source_weight[p][j] = amplitude_j * measure_j * b_p(H_j); heights are the
  // (passive) source positions used for the query-rank search.
  void assemble_competition_field(
      const std::array<std::vector<S>, comp_rank>& source_weight,
      const std::vector<double>& source_heights, const CanopyShape<S>& canopy) {
    competition_field.assemble(source_weight);
    competition_source_heights = source_heights;
    competition_canopy = canopy;
    competition_field_ready = true;
  }
  void clear_competition_field() { competition_field_ready = false; }

  // Count of source cohorts at least as tall as z (competition_source_heights is
  // descending) -- the query's rank for separable_field.at(a(z), rank-1).
  std::size_t n_sources_at_least(double z) const {
    auto it = std::upper_bound(competition_source_heights.begin(),
                               competition_source_heights.end(), z,
                               std::greater<double>());
    return static_cast<std::size_t>(it - competition_source_heights.begin());
  }

  void set_fixed_environment(double value, double height_max) {
    light_availability.set_fixed_value(value, height_max);
  }

  void set_fixed_environment(double value) {
    double height_max = 150.0;
    set_fixed_environment(value, height_max);
  }

  // Fix the light field to a scalar value that may carry AD derivatives (e.g. the
  // resident-feedback derivative of the competition a cohort sees). Unlike
  // set_fixed_environment(double), which would strip those derivatives, this keeps
  // them. Used to build the frozen environment for the forward-over-reverse dg/dh
  // evaluation (plant#39): the query is frozen, so a field fixed at the cohort's
  // current competition reproduces the read while preserving d(competition)/dtheta.
  void set_fixed_environment_scalar(S value, double height_max) {
    light_availability.set_fixed_value_scalar(value, height_max);
  }

  S get_environment_at_height(S height) const {
    if (competition_field_ready) {
      // Exact field: light = exp(-A(z)). The query factors a_p(z) carry the active
      // query-height derivative (the self-shading feedback the spline drops); the
      // rank restricts the sum to cohorts at least as tall (descending heights).
      const double z = odelia::util::to_passive(height);
      const std::size_t k = n_sources_at_least(z);  // cohorts with H >= z
      if (k == 0) return S(1.0);                     // nothing taller -> full light
      const S A = competition_field.at(
          competition_canopy.template shading_query_factors<S>(height), k - 1);
      return exp(-A);
    }
    // Fallback (field not assembled): the fitted spline with the query-height
    // derivative frozen (§15 / plant#39 -- an unreliable interpolant tangent).
    return light_availability.get_value_at_height_frozen_query(height);
  }

  // d(light)/d(height) by secant -- the coupling channel for the density-transport
  // dg/dh (plant#39 / design B). The field owns its slope; the dg/dh seam reads this
  // rather than hand-rolling a secant, with step+direction from the same Control the
  // production stencil uses.
  S get_environment_slope_at_height(S height, double step, int direction) const {
    return light_availability.slope_at_height(height, step, direction);
  }

  virtual void r_init_interpolators(const std::vector<double> &state)
  {
    light_availability.r_init_interpolators(state);
  }

  virtual Rcpp::List r_get_state() const
  {
    return Rcpp::List::create(_["light_availability"] = this->time); //      light_availability);
  }

  // Core functions
  template <typename Function>
  void compute_environment(Function f_compute_competition, double height_max, bool rescale) {

    // Define an anonymous function to use in creation of light_availability spline
    // Note: extinction coefficient was already applied in strategy, so
    // f_compute_competition gives sum of projected leaf area (k L) across species. Just need to apply Beer's law, E = exp(- (k L))
    // Returns S: the competition is active on a resident gradient pass, so the
    // light-availability knots carry the self-shading derivative (Beer's law).
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
// `plant::K93_Environment`.
using K93_Environment = K93_Environment_<double>;

}

#endif
