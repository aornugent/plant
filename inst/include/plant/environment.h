// -*-c++-*-
#ifndef PLANT_PLANT_ENVIRONMENT_H_
#define PLANT_PLANT_ENVIRONMENT_H_

#include <plant/control.h>
#include <odelia/interpolator.hpp>
#include <plant/adaptive_interpolator.h>
#include <odelia/ode_interface.hpp>
#include <plant/internals.h>
#include <plant/util.h>
#include <unordered_map>
#include <Rcpp.h>
#include <plant/extrinsic_drivers.h>

using namespace Rcpp;

namespace plant {

class Environment {
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

  // How many entries of resource_depletion compute_rates reads, and so how
  // long each individual's consumption vector must be.
  virtual size_t n_resources() const { return 0; }

  // What a cohort reads out of the shared environment, so a tape can supply
  // those values as active inputs instead of reading them as constants.
  virtual size_t n_cohort_reads() const { return 0; }

  template <typename It> It cohort_reads(It it) const { return it; }

  template <typename It> It set_cohort_reads(It it) { return it; }

  virtual void compute_rates(std::vector<double> const& resource_depletion){};

  // One aux slot per resource, holding the uptake the individuals supplied to
  // compute_rates. The soil rates subtract it and keep no record of it, so
  // without this slot the consumption is unrecoverable from the state.
  size_t aux_size() const { return n_resources(); }

  // The soil state is passive by declaration: its parameter sensitivity travels
  // by the adjoint ODE, and set_cohort_reads is where an active value arrives.
  template <typename It> It set_ode_state(It it) {
    for (size_t i = 0; i < vars.state_size; i++) {
      vars.states[i] = odelia::util::to_passive(*it++);
    }
    return it;
  }

  template <typename It> It ode_state(It it) const {
    for (size_t i = 0; i < vars.state_size; i++) {
      *it++ = vars.states[i];
    }
    return it;
  }

  template <typename It> It ode_rates(It it) const {
    for (size_t i = 0; i < vars.state_size; i++) {
      *it++ = vars.rates[i];
    }
    return it;
  }

  template <typename It> It ode_aux(It it) const {
    util::check_length(resource_uptake.size(), aux_size());
    for (size_t i = 0; i < aux_size(); i++) {
      *it++ = resource_uptake[i];
    }
    return it;
  }

  // Passive: aux carries a linearisation point and a branch condition, and both
  // want a value. A metric reading aux is a block output and overturns this.
  template <typename It> It set_ode_aux(It it) {
    util::check_length(resource_uptake.size(), aux_size());
    for (size_t i = 0; i < aux_size(); i++) {
      resource_uptake[i] = odelia::util::to_passive(*it++);
    }
    return it;
  }

  // n_resources() is the count and this is the buffer it sizes, so an environment
  // that changes its resource count calls this and the two cannot disagree.
  void resize_resource_uptake() { resource_uptake.assign(n_resources(), 0.0); }

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

  double get_environment_at_height(double height) const { return 0.0; };

  virtual ~Environment() = default;

  double time;

  size_t species_arriving_index;

  Internals<double> vars;
  ExtrinsicDrivers extrinsic_drivers;

  // Uptake per resource, as compute_rates received it. Sized to n_resources()
  // by whichever derived environment defines the resources.
  std::vector<double> resource_uptake;

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
}
#endif
