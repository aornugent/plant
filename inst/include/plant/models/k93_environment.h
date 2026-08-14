// -*-c++-*-
#ifndef PLANT_PLANT_K93_ENVIRONMENT_H_
#define PLANT_PLANT_K93_ENVIRONMENT_H_

#include <plant/environment.h>
#include <plant/resource_spline.h>

using namespace Rcpp;

namespace plant {

class K93_Environment : public Environment {
public:
  K93_Environment() {
    time = 0.0;
  };

  // Metres between the light field's knots. Four times finer than FF16's and
  // TF24's because K93's stand is four times shorter -- it starts at 2 m and
  // reaches about 8.5, where theirs reach 18 -- and a grid of constants gives a
  // stand resolution in proportion to its height. At 0.1 K93's offspring
  // production sits 1.8e-03 from the refined answer; at 0.025 it sits 7.7e-05,
  // which is where knots tied to the canopy top had it.
  constexpr static double light_knot_spacing = 0.025;

  // Light interface
  ResourceSpline<double> light_availability{light_knot_spacing};

  void set_fixed_environment(double value, double height_max) {
    light_availability.set_fixed_value(value, height_max);
  }

  void set_fixed_environment(double value) {
    double height_max = 150.0;
    set_fixed_environment(value, height_max);
  }

  double get_environment_at_height(double height) const {
    return light_availability.get_value_at_height(height);
  }

  virtual void r_init_interpolators(const std::vector<double> &state)
  {
    light_availability.r_init_interpolators(state);
  }

  virtual Rcpp::List r_get_state() const
  {
    return Rcpp::List::create(_["light_availability"] = time); //      light_availability);
  }

  // Core functions
  template <typename Function>
  void compute_environment(Function f_compute_competition_and_slope, double height_max) {

    // Beer's law on the competition profile A, whose extinction coefficient the
    // strategy has already applied: E = exp(-A) and dE/dz = -A' exp(-A).
    auto f_light_availability = [&](double height) -> std::pair<double, double>
    {
      const std::pair<double, double> as = f_compute_competition_and_slope(height);
      const double E = exp(-as.first);
      return {E, -(as.second * E)};
    };

    light_availability.compute_environment(f_light_availability, height_max);
  }

  virtual void clear_environment() {
    light_availability.clear();
  }

};

}

#endif
