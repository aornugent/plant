// -*-c++-*-
#ifndef PLANT_PLANT_K93_ENVIRONMENT_H_
#define PLANT_PLANT_K93_ENVIRONMENT_H_

#include <plant/environment.h>
#include <plant/resource_spline.h>

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
    // Freeze the query-height derivative on the rate path (§15 Gate 1 finding):
    // the query height is a cohort's evolving ODE state, and the interpolant's
    // analytic tangent w.r.t. it is an unreliable slope that compounds across the
    // replay. Knot-value derivatives (resident self-shading) are still carried.
    // Also nested-type safe, so the forward-over-reverse dg/dh evaluation reads
    // the same frozen field.
    return light_availability.get_value_at_height_frozen_query(height);
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
