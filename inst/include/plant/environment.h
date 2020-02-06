// -*-c++-*-
#ifndef PLANT_PLANT_ENVIRONMENT_H_
#define PLANT_PLANT_ENVIRONMENT_H_

#include <plant/control.h>
#include <plant/disturbance.h>
#include <plant/interpolator.h>
#include <plant/adaptive_interpolator.h>
#include <plant/util.h>

using namespace Rcpp;

namespace plant {

class Environment {
public:
  double canopy_openness(double height) const;
  template <typename Function>
  void compute_environment(Function f_canopy_openness, double height_max);
  template <typename Function>
  void rescale_environment(Function f_canopy_openness, double height_max);
  double patch_survival() const;
  double patch_survival_conditional(double time_at_birth) const;
  void clear();
  void clear_environment();

  // NOTE: Interface here will change
  double seed_rain_dt() const;
  void set_seed_rain_index(size_t x);

  // * R interface
  void r_set_seed_rain_index(util::index x);

  double time;
  Disturbance disturbance_regime;
  interpolator::Interpolator environment_interpolator;

private:
  std::vector<double> seed_rain;
  size_t seed_rain_index;
  interpolator::AdaptiveInterpolator environment_generator;
};

template <typename Function>
void Environment::compute_environment(Function f_canopy_openness,
                                            double height_max) {
  environment_interpolator =
    environment_generator.construct(f_canopy_openness, 0, height_max);
}

template <typename Function>
void Environment::rescale_environment(Function f_canopy_openness,
                                            double height_max) {
  std::vector<double> h = environment_interpolator.get_x();
  const double min = environment_interpolator.min(), // 0.0?
    height_max_old = environment_interpolator.max();

  util::rescale(h.begin(), h.end(), min, height_max_old, min, height_max);
  h.back() = height_max; // Avoid round-off error.

  environment_interpolator.clear();
  for (auto hi : h) {
    environment_interpolator.add_point(hi, f_canopy_openness(hi));
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

template <typename T>
Environment make_environment(T p) {
  return Environment(p.disturbance_mean_interval, p.seed_rain, p.control);
}

}

#endif
