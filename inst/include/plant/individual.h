// -*-c++-*-
#ifndef PLANT_PLANT_PLANT_MINIMAL_H_
#define PLANT_PLANT_PLANT_MINIMAL_H_

#include <memory> // std::shared_ptr
#include <odelia/ode_interface.hpp>
#include <vector>
#include <plant/internals.h>
#include <plant/uniroot.h>


namespace plant {

template <typename T, typename E> class Individual {
public:
  typedef T strategy_type;
  typedef E environment_type;
  typedef typename strategy_type::ptr strategy_type_ptr;
  // The scalar this individual's state and rates carry, taken from the strategy
  // (double on the production path; active when a trait is seeded).
  using value_type = typename strategy_type::value_type;
  // for the time being...
  Individual(strategy_type_ptr s) : strategy(s) {
    if (strategy->aux_index.size() != s->aux_size()) {
      strategy->refresh_indices();
    }
    // The strategy owns the knowledge of which state/aux slots its rate
    // functions need: Individual passes the whole Internals to the strategy
    // (see compute_competition / net_mass_production_dt below) rather than
    // resolving strategy-specific aux slots here (#266). The strategy still
    // reads them by cached integer index, so the #466 hot-path optimisation is
    // preserved.
    vars.resize(strategy_type::state_size(), s->aux_size()); // = Internals(strategy_type::state_size());
    set_state("height", strategy->initial_height());
  }

  // useage: state(HEIGHT_INDEX)
  value_type state(std::string name) const {
    return vars.state(strategy->state_index.at(name));
  }
  value_type state(int i) const { return vars.state(i); }

  // useage:_rate("area_heartwood")
  value_type rate(std::string name) const {
    return vars.rate(strategy->state_index.at(name));
  }
  value_type rate(int i) const { return vars.rate(i); }

  // useage: set_state("height", 2.0)
  void set_state(std::string name, value_type v) {
    int i = strategy->state_index.at(name);
    vars.set_state(i, v);
    strategy->update_dependent_aux(i, vars);
  }
  void set_state(int i, value_type v) {
    vars.set_state(i, v);
    strategy->update_dependent_aux(i, vars);
  }

  // aux vars by name and index
  value_type aux(std::string name) const {
    return vars.aux(strategy->aux_index.at(name));
  }
  value_type aux(int i) const { return vars.aux(i); }

  // set # consumable resources based on env. variables
  void resize_consumption_rates(int i) {
    vars.resize_consumption_rates(i);
  }
  value_type consumption_rate(int i) const { return vars.consumption_rate(i); }

  value_type compute_competition(double z) const {
    return strategy->compute_competition(z, vars);
  }

  // Per-individual census quantities (the Psi of a population census); delegate to
  // the strategy, which reads this individual's multivariate state (#266). Only
  // instantiated for strategies that define them (e.g. FF16's census gradient).
  value_type census_leaf_area() const { return strategy->census_leaf_area(vars); }
  value_type census_mass() const { return strategy->census_mass(vars); }
  value_type census_basal_area() const { return strategy->census_basal_area(vars); }

  // Seed strategy-specific initial ODE states (e.g. an acclimating tracked
  // state) given the birth environment. No-op for strategies that don't need it.
  void set_initial_states(const environment_type& environment) {
    strategy->set_initial_states(environment, vars);
  }

  void compute_rates(const environment_type& environment) {
    if (vars.resource_size != environment.ode_size()) {
      // handles when Individual hasn't been instantiated in a Patch (ie with an environment)
      vars.resize_consumption_rates(environment.ode_size());
    }
    strategy->compute_rates(environment, vars);
  }

  value_type establishment_probability(const environment_type &environment) {
    return strategy->establishment_probability(environment);
  }

  value_type net_mass_production_dt(const environment_type &environment) {
    // TODO(#483):  maybe reuse intervals? default false
    return strategy->net_mass_production_dt(environment, vars);
  }

  // The low-level strategy parameters this individual's rates depend on, as
  // pointers for the gradient driver to seed active (§8.1). Forwarded from the
  // strategy and valid for the individual's lifetime. Only instantiated for a
  // strategy whose IndividualRunner is actually differentiated, so a strategy
  // without field_ptrs() is not forced to grow one.
  std::vector<value_type*> ad_parameters() { return strategy->field_ptrs(); }

  // Re-run the strategy's one-time precomputation so its derived quantities are
  // recomputed from the current parameters. Needed on a gradient pass, where the
  // parameters are seeded active *after* construction: without this the cached
  // derived quantities (e.g. area_leaf_0) stay constant w.r.t. the seed and their
  // contribution to the gradient is silently dropped.
  void prepare_strategy() { strategy->prepare_strategy(); }

  // Restore the birth ODE state: the strategy-derived initial height, every other
  // state zero. Mirrors construction, but explicitly zeros first (resize only
  // clears newly-added slots) so it also re-inits an already-grown individual for
  // a fresh replay from a re-prepared strategy.
  void set_birth_state() {
    vars.resize(strategy_type::state_size(), strategy->aux_size());
    for (size_t i = 0; i < strategy_type::state_size(); i++) {
      set_state(static_cast<int>(i), value_type(0.0));
    }
    set_state(HEIGHT_INDEX, strategy->initial_height());
  }

  // * ODE interface
  static size_t ode_size() { return strategy_type::state_size(); }
  static std::vector<std::string> ode_names() { return strategy_type::state_names(); }

  // why doesn't strategy->aux_size() also need to be const-qualified? is it because strategy is a pointer?
  size_t aux_size() const { return strategy->aux_size(); }
  std::vector<std::string> aux_names() { return strategy->aux_names(); }

  // The ODE-serialization seam is templated on the iterator so the state flows
  // at whatever scalar the solver drives -- vector<double>::iterator on the
  // production pass, vector<active>::iterator on a gradient pass.
  template <typename It>
  It set_ode_state(It it) {
    for (size_t i = 0; i < vars.state_size; i++) {
      vars.states[i] = *it++;
      strategy->update_dependent_aux(i, vars);
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

  template <typename It>
  It ode_aux(It it) const {
    for (size_t i = 0; i < vars.aux_size; i++) {
      *it++ = vars.auxs[i];
    }
    return it;
  }

  // Single individual methods

  // Used in the stochastic model:
  value_type mortality_probability() const { return 1 - exp(-state(MORTALITY_INDEX)); }

  void reset_mortality() { set_state("mortality", 0.0); }

  value_type growth_rate_given_height(value_type height, const environment_type& environment) {
    // Called repeatedly from the finite-difference gradient (Node::
    // growth_rate_gradient), so address height by integer slot rather than the
    // "height" string-map lookup (see #466). `height` carries value_type so the
    // active pass can perturb around an ACTIVE abscissa (the cohort's current
    // height, which depends on the seeded parameters); the double path passes a
    // plain double as before.
    set_state(HEIGHT_INDEX, height);
    compute_rates(environment);
    return rate(HEIGHT_INDEX);
  }

  double resource_compensation_point() {
    environment_type env = environment_type();

    auto target = [&] (double x) mutable -> double {
      env.set_fixed_environment(x, 100);
      compute_rates(env);
      return net_mass_production_dt(env);
    };

    const double f1 = target(1.0);
    if (f1 < 0.0) {
      return NA_REAL;
    } else {
      const double tol = control().offspring_production_tol;
      const size_t max_iterations = control().offspring_production_iterations;
      return util::uniroot(target, 0.0, 1.0, tol, max_iterations);
    }
  }



  std::string strategy_name() const { return strategy->name; }

  // * R interface
  strategy_type r_get_strategy() const { return *strategy.get(); }

  // ! External R code depends on knowing r internals for like growing plant to
  // ! height or something
  Internals_<value_type> r_internals() const { return vars; }
  const Control &control() const { return strategy->control; }

private:
  strategy_type_ptr strategy;
  Internals_<value_type> vars;
};

template <typename T, typename E> Individual<T,E> make_individual(T s) {
  return Individual<T,E>(make_strategy_ptr(s));
}

} // namespace plant

#endif
