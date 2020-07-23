// -*-c++-*-
#ifndef PLANT_PLANT_ENVIRONMENT_H_
#define PLANT_PLANT_ENVIRONMENT_H_

#include <plant/control.h>
#include <plant/disturbance.h>
#include <plant/interpolator.h>
#include <plant/adaptive_interpolator.h>
#include <plant/util.h>
#include <plant/ode_interface.h>

using namespace Rcpp;

namespace plant {

class Environment {
public:

  Environment() :
    time(0.0),
    disturbance_regime(30),
    seed_rain({}),
    seed_rain_index(0),
    environment_generator(interpolator::AdaptiveInterpolator()) {
  };

  template <typename Function>
  void compute_environment(Function f, double height_max);
  template <typename Function>
  void rescale_environment(Function f, double height_max);

  void set_fixed_environment(double value, double height_max) {
    std::vector<double> x = {0, height_max/2.0, height_max};
    std::vector<double> y = {value, value, value};
    clear_environment();
    environment_interpolator.init(x, y);
  }

  void set_fixed_environment(double value) {
    double height_max = 150.0;
    set_fixed_environment(value, height_max);
  }

  void compute_rates() {}

  // Dummy iterators: do nothing if the environment has no state. 
  ode::iterator ode_state(ode::iterator it) const { return it; }
  ode::const_iterator set_ode_state(ode::const_iterator it) { return it; }
  ode::iterator ode_rates(ode::iterator it) const { return it; }
  int ode_size() const { return 0; }

  // Computes the probability of survival from 0 to time.
  double patch_survival() const {
    return disturbance_regime.pr_survival(time);
  }
  // Computes the probability of survival from time_at_birth to time, by
  // conditioning survival over [0,time] on survival over
  // [0,time_at_birth].
  double patch_survival_conditional(double time_at_birth) const {
    return disturbance_regime.pr_survival_conditional(time, time_at_birth);
  }
  // Reset the environment.
  void clear() {
    time = 0.0;
    clear_environment();
  }
  void clear_environment() {
    environment_interpolator.clear();
  }
  double seed_rain_dt() const {
    if (seed_rain.empty()) {
      Rcpp::stop("Cannot get seed rain for empty environment");
    }
    return seed_rain[seed_rain_index];
  }
  void set_seed_rain_index(size_t x) {
    seed_rain_index = x;
  }
  // * R interface
  void r_set_seed_rain_index(util::index x) {
    set_seed_rain_index(x.check_bounds(seed_rain.size()));
  }


  double get_environment_at_height(double height) const {
    const bool within = height <= environment_interpolator.max();
    // TODO: change maximum - here hard-coded to 1.0
    return within ? environment_interpolator.eval(height) : 1.0;
  }

  double time;
  Disturbance disturbance_regime;
  interpolator::Interpolator environment_interpolator;
  std::vector<double> seed_rain;
  size_t seed_rain_index;
  interpolator::AdaptiveInterpolator environment_generator;
};

template <typename Function>
void Environment::compute_environment(Function f, double height_max) {
  environment_interpolator = environment_generator.construct(f, 0, height_max);
}

template <typename Function>
void Environment::rescale_environment(Function f, double height_max) {
  std::vector<double> h = environment_interpolator.get_x();
  const double min = environment_interpolator.min(), // 0.0?
    height_max_old = environment_interpolator.max();

  util::rescale(h.begin(), h.end(), min, height_max_old, min, height_max);
  h.back() = height_max; // Avoid round-off error.

  environment_interpolator.clear();
  for (auto hi : h) {
    environment_interpolator.add_point(hi, f(hi));
  }
  environment_interpolator.initialise();
}

inline interpolator::AdaptiveInterpolator
make_interpolator(const Control& control) {
  using namespace interpolator;
  return AdaptiveInterpolator(control.environment_light_tol,
                              control.environment_light_tol,
                              control.environment_light_nbase,
                              control.environment_light_max_depth);
}

}

#endif
