// -*-c++-*-
#ifndef PLANT_PLANT_ENVIRONMENT_H_
#define PLANT_PLANT_ENVIRONMENT_H_

#include <plant/control.h>
#include <odelia/interpolator.hpp>
#include <odelia/ode_interface.hpp>
#include <plant/internals.h>
#include <plant/util.h>
#include <unordered_map>
#include <Rcpp.h>
#include <plant/extrinsic_drivers.h>

using namespace Rcpp;

namespace plant {

// The environment a plant couples to (light, and for TF24 soil water). Templated
// on the scalar S carried by its ODE state / knot values; S = double is the
// production path (the `Environment` alias below). FF16/K93 carry no environment
// ODE state (ode_size() == 0), so their seam methods are pass-throughs; TF24's
// soil state rides these on the double alias. The ODE-serialisation seam is
// templated on the iterator so the state flows at whatever scalar the solver
// drives.
template <class S = double>
class Environment_ {
public:
  template <typename Function>
  void compute_environment(Function f, double height_max, bool rescale);

  void set_fixed_environment(double value, double height_max);
  void set_fixed_environment(double value);

  // Configure the crown shading model for the light profile. Default: no-op;
  // only FF16_Environment builds an alternative (stepped) profile. Called once
  // from the Patch constructor with the run's Control settings.
  virtual void set_shading_model(const std::string& /*model*/,
                                 double /*layer_optical_depth*/,
                                 double /*layer_smoothing*/) {}

  // ODE interface: do nothing if the environment has no state.
  size_t ode_size() const { return vars.state_size; }
  // resource_depletion carries S so the resident soil coupling (plant water
  // uptake -> soil state) differentiates; FF16/K93 have no soil state and ignore
  // it. At S = double this is the previous std::vector<double> signature.
  virtual void compute_rates(std::vector<S> const& resource_depletion){};

  template <typename It>
  It set_ode_state(It it) {
    for (size_t i = 0; i < vars.state_size; i++) {
      vars.states[i] = *it++;
    }
    return it;
  }

  template <typename It>
  It ode_state(It it) const {
    for (size_t i = 0; i < vars.state_size; i++) {
      *it++ = vars.states[i];
    }
    return it;
  }

  template <typename It>
  It ode_rates(It it) const {
    for (size_t i = 0; i < vars.state_size; i++) {
      *it++ = vars.rates[i];
    }
    return it;
  }

  virtual Rcpp::List r_get_state() const
  {
    return Rcpp::List::create(_["time"] = time);
  }

  // Reset the environment
  void clear() {
    time = 0.0;
    clear_environment();
  }

  virtual void clear_environment() {}

  virtual void r_init_interpolators(const std::vector<double>& state) {}

  S get_environment_at_height(S height) const { return static_cast<S>(0.0); };

  virtual ~Environment_() = default;

  double time;

  size_t species_arriving_index;

  Internals_<S> vars;
  ExtrinsicDrivers extrinsic_drivers;

  // The
  std::vector<std::string> extrinsic_drivers_get_names() const
  {
    return  extrinsic_drivers.get_names();
  }

  void extrinsic_drivers_set_constant(std::string driver_name, double value)
  {
    extrinsic_drivers.set_constant(driver_name, value);
  }

  void extrinsic_drivers_set_variable(std::string driver_name, std::vector<double> const &x, std::vector<double> const &y)
  {
    extrinsic_drivers.set_variable(driver_name, x, y);
  }

  double extrinsic_drivers_evaluate(std::string driver_name, double x) const
  {
    return extrinsic_drivers.evaluate(driver_name, x);
  }

  std::vector<double> extrinsic_drivers_evaluate_range(std::string driver_name, std::vector<double> const &x) const
  {
    return extrinsic_drivers.evaluate_range(driver_name, x);
  }

};

// The double instantiation is the production path and what every current caller
// (TF24_Environment, and the FF16/K93 environments' base) means by `Environment`.
using Environment = Environment_<double>;
}
#endif
