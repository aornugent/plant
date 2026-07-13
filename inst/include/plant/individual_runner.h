// -*-c++-*-
#ifndef PLANT_RUNNER
#define PLANT_RUNNER

#include <plant/individual.h>
#include <plant/environment.h>
#include <odelia/ode_solver.hpp>

namespace plant {
namespace tools {

// A fixed-dimension single-plant ODE System: one Individual growing in a given
// environment. The scalar it carries is the Individual's -- `double` on the
// production path, an active type when a strategy parameter is seeded for a
// gradient. The ODE-serialization seam is templated on the iterator so the same
// body drives both.
template <typename T, typename E>
class IndividualRunner {
public:
  using value_type = typename Individual<T,E>::value_type;

  IndividualRunner(Individual<T,E> individual_, E environment_)
    : individual(individual_), environment(environment_),
      t0(environment_.time) {
    individual.compute_rates(environment);
  }

  static size_t ode_size() {return Individual<T,E>::ode_size();}

  double ode_time() const {return environment.time;}
  double ode_t0() const {return t0;}

  template <typename It>
  It set_ode_state(It it, double time) {
    it = individual.set_ode_state(it);
    environment.time = time;
    individual.compute_rates(environment);
    return it;
  }
  template <typename It>
  It ode_state(It it) const {
    return individual.ode_state(it);
  }
  template <typename It>
  It ode_rates(It it) const {
    return individual.ode_rates(it);
  }

  // Re-derive the birth state from the current (seeded) parameters and recompute
  // rates, so the gradient driver replays the same growth from freshly-seeded
  // parameters (§8.3 run()). prepare_strategy() must run here, not just at
  // construction, or a parameter whose effect is mediated by a precomputed
  // quantity (area_leaf_0, eta_c) loses that part of its derivative.
  void reset() {
    individual.prepare_strategy();
    individual.set_birth_state();
    environment.time = t0;
    individual.compute_rates(environment);
  }

  // The differentiable inputs the driver seeds a subset of: the strategy's
  // low-level parameters (§8.1). A single plant's birth ODE state is derived from
  // those parameters (via prepare_strategy / height_seed), not an independent
  // input, so it seeds no initial state of its own.
  std::vector<value_type*> ad_parameters() { return individual.ad_parameters(); }
  std::vector<value_type*> ad_initial_state() { return {}; }

  Individual<T,E> individual;
  E environment;

private:
  double t0;
};

}
}

#endif /* PLANT_RUNNER */
