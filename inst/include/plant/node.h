// -*-c++-*-
#ifndef NODE
#define NODE

#include <plant/environment.h>
#include <plant/individual.h>
#include <odelia/ode_interface.hpp>
#include <limits> // std::numeric_limits

namespace plant {

template <typename T, typename E>
class Node {
public:
  using value_type = typename T::value_type;

  typedef T        strategy_type;
  typedef E        environment_type;
  typedef Individual<T,E> individual_type;
  typedef typename strategy_type::ptr strategy_type_ptr;
  Node(strategy_type_ptr s);

  void compute_rates(const environment_type& environment, double pr_patch_survival);
  void compute_initial_conditions(const environment_type& environment, double pr_patch_survival, double birth_rate);

  double height() const {return individual.state(HEIGHT_INDEX);}
  double compute_competition(double z) const;
  // The node's competition contribution and its vertical derivative, both
  // weighted by density. The first entry equals compute_competition(z) exactly.
  std::pair<double, double> compute_competition_and_slope(double z) const;
  double fecundity() const {return offspring_produced_survival_weighted;}

  // The two rates the transport term is built from, and the write for its
  // result. Species::growth_rate_gradient differences growth across neighbouring
  // nodes, so it reads these rather than perturbing a copy.
  double growth_rate() const {return individual.rate(HEIGHT_INDEX);}
  double mortality_rate() const {return individual.rate(MORTALITY_INDEX);}
  void set_log_density_rate(double rate) {log_density_dt = rate;}

  // Bookkeeping recorded at the moment the node is introduced, so that
  // lifetime-fitness calculations need not look these up after the run.
  // patch_density_at_birth is the (unnormalised) probability density of a
  // patch having the node's introduction age, i.e. survival_weighting->density.
  void set_introduction(double time, double patch_density) {
    node_introduction_time = time;
    patch_density_at_birth = patch_density;
  }
  double introduction_time() const {return node_introduction_time;}
  double patch_density() const {return patch_density_at_birth;}
  double get_pr_patch_survival_at_birth() const {return pr_patch_survival_at_birth;}
  double get_log_density_rate() const {return log_density_dt;}

  // Restore birth bookkeeping for a node imported from an exported patch state,
  // without re-running compute_initial_conditions (which would overwrite the
  // loaded ODE state). pr_patch_survival_at_birth feeds the fecundity rate;
  // node_introduction_time and patch_density_at_birth feed lifetime-fitness
  // integrals. Required for a resumed run to reproduce the original trajectory.
  void set_birth_state(double time, double patch_density_in,
                       double pr_patch_survival) {
    node_introduction_time = time;
    patch_density_at_birth = patch_density_in;
    pr_patch_survival_at_birth = pr_patch_survival;
  }

  // Lifetime offspring of this node, weighted by the probability of
  // landing in a patch of the node's age and by survival during dispersal.
  double weighted_fecundity(double S_D) const {
    return offspring_produced_survival_weighted * patch_density_at_birth * S_D;
  }

  // Unfortunate, but need a get_ here because of name shadowing...
  double get_log_density() const {return log_density;}
  // exp(log_density); can overflow to +Inf when the SCM density equation runs
  // away (see Patch::check_finite_node_densities).
  double get_density() const {return density;}
  void set_log_density(double x) {
    log_density = x;
    density = exp(log_density);
  }

  // ODE interface.
  //
  // NOTE: We are a time-independent model here so no need to pass
  // time in as an argument.  All the bits involving time are taken
  // care of by Environment for us.
  // +2 for log_density and offspring_production_dt
  static size_t ode_size() { return strategy_type::state_size() + 2; }
  size_t aux_size() const { return individual.aux_size(); }
  template <typename It> It set_ode_state(It it);
  template <typename It> It ode_state(It it) const;
  template <typename It> It ode_rates(It it) const;
  template <typename It> It ode_aux(It it) const;
  template <typename It> It set_ode_aux(It it);

  static std::vector<std::string> ode_names() {
    std::vector<std::string> names = strategy_type::state_names();
    names.push_back("offspring_produced_survival_weighted");
    names.push_back("log_density");
    return names;
  }

  void resize_consumption_rates(int i) {
    individual.resize_consumption_rates(i);
  }

  double consumption_rate(int i) const {
    return individual.consumption_rate(i) * density;
  }

  individual_type individual;

private:
  double log_density;
  double log_density_dt;
  double density; // hmm...
  double offspring_produced_survival_weighted;
  double offspring_produced_survival_weighted_dt;
  double pr_patch_survival_at_birth;

  // Recorded at introduction (see set_introduction).
  double node_introduction_time;
  double patch_density_at_birth;
};

template <typename T, typename E>
Node<T,E>::Node(strategy_type_ptr s)
  : individual(s),
    log_density(-std::numeric_limits<double>::infinity()),
    log_density_dt(0),
    density(0),
    offspring_produced_survival_weighted(0),
    offspring_produced_survival_weighted_dt(0),
    node_introduction_time(0),
    patch_density_at_birth(0) {
}

template <typename T, typename E>
void Node<T,E>::compute_rates(const environment_type& environment,
                                double pr_patch_survival) {
  individual.compute_rates(environment);

  // log_density_dt is not set here: it needs the transport term, which is a
  // property of the node's place among its neighbours rather than of the node.
  // Species::compute_rates writes it.

  // survival_individual: converts from the mean of the poisson process (on
  // [0,Inf)) to a probability (on [0,1]).
  double survival_individual = exp(-individual.state(MORTALITY_INDEX));
  if (!util::is_finite(survival_individual)) {
    // This is caused by NaN values in plant.mortality and log
    // density; this should only be an issue when density is so low
    // that we can throw these away.  I think that with smaller step
    // sizes this is better behaved too?
    survival_individual = 0.0;
  }

  offspring_produced_survival_weighted_dt =
    individual.rate(FECUNDITY_INDEX) * survival_individual *
    pr_patch_survival / pr_patch_survival_at_birth;
}

// NOTE: There will be a discussion of why the mortality rate initial
// condition is -log(establishment_probability) in the documentation
// that Daniel is working out.
//
// NOTE: The initial condition for log_density is also a bit tricky, and
// defined on p 7 at the moment.
template <typename T, typename E>
void Node<T,E>::compute_initial_conditions(const environment_type& environment,
                                             double pr_patch_survival, double birth_rate) {
  pr_patch_survival_at_birth = pr_patch_survival;
  // Seed strategy-specific initial states (e.g. TF24f's tracked psi at its
  // optimum) before the first rates evaluation, so the birth growth rate uses
  // the initialised operating point rather than a default.
  individual.set_initial_states(environment);
  compute_rates(environment, pr_patch_survival);
  // The transport term needs an interval below, and the inflow boundary node has
  // none. Its log_density rate is never integrated: a node introduced from this
  // one is rated by Species::compute_rates before the solver reads it.
  log_density_dt = 0.0;

  const double pr_estab =
    individual.establishment_probability_of_newborn(environment);
  individual.set_state("mortality", -log(pr_estab));
  const double g = individual.rate(HEIGHT_INDEX);
  // NOTE: log(0.0) -> -Inf, which should behave fine.
  set_log_density(g > 0 ? log(birth_rate * pr_estab / g) : log(0.0));

  // NOTE: It's *possible* here that we need to set
  // individual.vars.mortality_dt to zero here, but I don't see that's
  // likely.
}

template <typename T, typename E>
double Node<T,E>::compute_competition(double height_) const {
  return density * individual.compute_competition(height_);
}

template <typename T, typename E>
std::pair<double, double>
Node<T,E>::compute_competition_and_slope(double height_) const {
  const std::pair<double, double> fs =
    individual.compute_competition_and_slope(height_);
  return {density * fs.first, density * fs.second};
}

// ODE interface -- note that the don't care about time in the node;
// only Patch and above does.
template <typename T, typename E>
template <typename It>
It Node<T,E>::set_ode_state(It it) {
  for (size_t i = 0; i < individual.ode_size(); i++) {
    individual.set_state(i, *it++);
  }
  offspring_produced_survival_weighted = *it++;
  set_log_density(*it++);
  return it;
}
template <typename T, typename E>
template <typename It>
It Node<T,E>::ode_state(It it) const {
  for (size_t i = 0; i < individual.ode_size(); i++) {
    *it++ = individual.state(i);
  }
  *it++ = offspring_produced_survival_weighted;
  *it++ = log_density;
  return it;
}
template <typename T, typename E>
template <typename It>
It Node<T,E>::ode_rates(It it) const {
  for (size_t i = 0; i < individual.ode_size(); i++) {
    *it++ = individual.rate(i);
  }
  *it++ = offspring_produced_survival_weighted_dt;
  *it++ = log_density_dt;
  return it;
}

template <typename T, typename E>
template <typename It>
It Node<T,E>::ode_aux(It it) const {
  for (size_t i = 0; i < individual.aux_size(); i++) {
    *it++ = individual.aux(i);
  }
  return it;
}

template <typename T, typename E>
template <typename It>
It Node<T,E>::set_ode_aux(It it) {
  return individual.set_ode_aux(it);
}


template <typename T, typename E>
Node<T,E> make_node(typename Node<T,E>::strategy_type s) {
  return Node<T,E>(make_strategy_ptr(s));
}

}

#endif /* NODE */
