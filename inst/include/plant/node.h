// -*-c++-*-
#ifndef NODE
#define NODE

#include <plant/environment.h>
#include <plant/gradient.h>
#include <plant/individual.h>
#include <odelia/ode_interface.hpp>
#include <limits> // std::numeric_limits

namespace plant {

// One node's adjoints from the water aggregation: the individual's uptake, which
// is a block output, and the density the quadrature weights it by. A height
// reaches total uptake through the individual's own rate, which the block
// carries.
struct node_uptake_adjoints {
  double uptake;
  double log_density;
};

// The inflow boundary node's own adjoints. It holds no ODE row, so each of these
// is pulled back through the condition that sets it rather than written to a
// state slot.
//
// The density has two slots because a stage evaluates the boundary condition
// twice, in two different fields: the light field is built with the first and the
// water aggregation runs after the second. They are sensitivities to different
// quantities, so one accumulator would transpose one derivative through the
// other's argument.
//
// `height` and `area_leaf` are pulled back together, through the seed's height.
struct boundary_node_adjoints {
  double area_leaf;
  double height;
  double density_in_field;
  double density_in_uptake;
};

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

  // Wrapper to growth_rate_gradient for testing
  value_type r_growth_rate_gradient(const environment_type& environment);

  value_type height() const {return individual.state(HEIGHT_INDEX);}
  value_type compute_competition(const value_type& z) const;
  // The node's competition contribution and its vertical derivative, both
  // weighted by density. The first entry equals compute_competition(z) exactly.
  std::pair<value_type, value_type>
  compute_competition_and_slope(const value_type& z) const;
  value_type fecundity() const {return offspring_produced_survival_weighted;}

  // The survival factor offspring_produced_survival_weighted_dt multiplies,
  // zero where compute_rates squashed a non-finite one.
  value_type survival_individual() const {
    const value_type s = exp(-individual.state(MORTALITY_INDEX));
    return util::is_finite(s) ? s : value_type(0.0);
  }

  // d(offspring_produced_survival_weighted_dt)/d(fecundity rate).
  value_type offspring_dt_dfecundity_rate(double pr_patch_survival) const {
    return survival_individual() * pr_patch_survival /
      pr_patch_survival_at_birth;
  }

  // The two rates the transport term is built from, and the write for its
  // result. Species::growth_rate_gradient differences growth across neighbouring
  // nodes, so it reads these rather than perturbing a copy.
  value_type growth_rate() const {return individual.rate(HEIGHT_INDEX);}
  value_type mortality_rate() const {return individual.rate(MORTALITY_INDEX);}

  // Bookkeeping recorded at the moment the node is introduced, so that
  // lifetime-fitness calculations need not look these up after the run.
  // patch_density_at_birth is the (unnormalised) probability density of a
  // patch having the node's introduction age, i.e. survival_weighting->density.
  void set_introduction(double time, double patch_density) {
    node_introduction_time = time;
    patch_density_at_birth = patch_density;
  }
  // The height growth rate this node was born at, i.e. |dh/dtau| exactly, at
  // birth. Current only for the boundary node (re-evaluated every step); frozen
  // at its own birth for an introduced one. Zero for a node loaded from an
  // exported state, which never ran compute_initial_conditions().
  double growth_rate_at_birth() const {return birth_growth_rate;}

  // Refresh only the birth date, leaving the rest of the bookkeeping alone.
  // Used for the not-yet-introduced boundary node, whose birth date is the
  // current time and so moves with every step; see
  // Species::set_new_node_birth_date().
  void set_introduction_time(double time) {node_introduction_time = time;}
  double introduction_time() const {return node_introduction_time;}
  double patch_density() const {return patch_density_at_birth;}
  double get_pr_patch_survival_at_birth() const {return pr_patch_survival_at_birth;}
  value_type get_log_density_rate() const {return log_density_dt;}

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
  value_type weighted_fecundity(const value_type& S_D) const {
    return offspring_produced_survival_weighted * patch_density_at_birth * S_D;
  }

  // Unfortunate, but need a get_ here because of name shadowing...
  value_type get_log_density() const {return log_density;}
  // exp(log_density); can overflow to +Inf when the SCM density equation runs
  // away (see Patch::check_finite_node_densities).
  value_type get_density() const {return density;}
  void set_log_density(const value_type& x) {
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

  value_type consumption_rate(int i) const {
    return individual.consumption_rate(i) * density;
  }

  individual_type individual;

private:
  // This is the gradient of growth rate with respect to height:
  value_type growth_rate_gradient(const environment_type& environment) const;

  value_type log_density;
  value_type log_density_dt;
  value_type density; // hmm...
  value_type offspring_produced_survival_weighted;
  value_type offspring_produced_survival_weighted_dt;
  double pr_patch_survival_at_birth;

  // Recorded at introduction (see set_introduction).
  double node_introduction_time;
  double patch_density_at_birth;
  // |dh/dtau| at birth; see growth_rate_at_birth().
  double birth_growth_rate;
};

template <typename T, typename E>
Node<T,E>::Node(strategy_type_ptr s)
  : individual(s),
    log_density(value_type(-std::numeric_limits<double>::infinity())),
    log_density_dt(0),
    density(0),
    offspring_produced_survival_weighted(0),
    offspring_produced_survival_weighted_dt(0),
    node_introduction_time(0),
    patch_density_at_birth(0),
    birth_growth_rate(0) {
}

template <typename T, typename E>
void Node<T,E>::compute_rates(const environment_type& environment,
                                double pr_patch_survival) {
  individual.compute_rates(environment);

  // NOTE: This must be called *after* compute_rates, but given we
  // need mortality_dt() that's always going to be the case.
  //
  // The coordinate branch lives in Individual::log_density_rate, which is one of
  // the cohort block's outputs, so the recorded block and this path cannot
  // disagree about which coordinate they are on.
  log_density_dt = individual.log_density_rate(environment);
  // survival_individual: converts from the mean of the poisson process (on
  // [0,Inf)) to a probability (on [0,1]).
  value_type survival_individual = exp(-individual.state(MORTALITY_INDEX));
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

  const value_type pr_estab =
    individual.establishment_probability_of_newborn(environment);
  individual.set_state("mortality", -log(pr_estab));
  // The birth-date axis of the node about to be introduced; Patch re-stamps
  // this with the exact introduction time as the node is pushed.
  node_introduction_time = environment.time;
  // dh/dtau = -g(H_0) at birth, so this is |dh/dtau| exactly, with no
  // differencing: characteristics are labelled by birth date, and a cohort born
  // an instant later starts an instant's growth behind. Recorded for every node
  // (a cheap read -- compute_rates() above has just set it, and the height
  // branch below needs it anyway), but only current for the *boundary* node,
  // which compute_initial_conditions() re-evaluates every step. For an
  // introduced node it is frozen at its own birth and so cannot serve as its
  // present-day Jacobian; see Species::height_jacobian().
  const value_type g = individual.rate(HEIGHT_INDEX);
  // A double diagnostic read back over the R boundary, so it takes the value.
  birth_growth_rate = odelia::util::to_passive(g);
  if (individual.control().node_density_in_birth_date) {
    set_log_density(log(birth_rate * pr_estab));
  } else {
    // NOTE: log(0.0) -> -Inf, which should behave fine.
    set_log_density(g > 0 ? log(birth_rate * pr_estab / g)
                          : value_type(log(0.0)));
  }

  // Need to check that the rates are valid after setting the
  // mortality value here (can go to -Inf and that requires squashing
  // the rate to zero).
  if (!util::is_finite(log_density)) {
    // Can do this at the same time that we do set_log_density, I think.
    log_density_dt = 0.0;
  }
  // NOTE: It's *possible* here that we need to set
  // individual.vars.mortality_dt to zero here, but I don't see that's
  // likely.
}

template <typename T, typename E>
typename Node<T,E>::value_type
Node<T,E>::growth_rate_gradient(const environment_type& environment) const {
  return individual.growth_rate_gradient(environment);
}

// Wrapper to growth_rate_gradient for testing
template <typename T, typename E>
typename Node<T,E>::value_type
Node<T,E>::r_growth_rate_gradient(const environment_type& environment) {
  // We need to compute the physiological variables here, first, so
  // that reusing intervals works as expected.  This would ordinarily
  // be taken care of because of the calling order of
  // compute_rates / growth_rate_gradient.
  individual.compute_rates(environment);
  return growth_rate_gradient(environment);
}

template <typename T, typename E>
typename Node<T,E>::value_type
Node<T,E>::compute_competition(const value_type& height_) const {
  return density * individual.compute_competition(height_);
}

template <typename T, typename E>
std::pair<typename Node<T,E>::value_type, typename Node<T,E>::value_type>
Node<T,E>::compute_competition_and_slope(const value_type& height_) const {
  const std::pair<value_type, value_type> fs =
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
    *it++ = util::as_iterator_scalar<It>(individual.state(i));
  }
  *it++ = util::as_iterator_scalar<It>(offspring_produced_survival_weighted);
  *it++ = util::as_iterator_scalar<It>(log_density);
  return it;
}
template <typename T, typename E>
template <typename It>
It Node<T,E>::ode_rates(It it) const {
  for (size_t i = 0; i < individual.ode_size(); i++) {
    *it++ = util::as_iterator_scalar<It>(individual.rate(i));
  }
  *it++ = util::as_iterator_scalar<It>(offspring_produced_survival_weighted_dt);
  *it++ = util::as_iterator_scalar<It>(log_density_dt);
  return it;
}

template <typename T, typename E>
template <typename It>
It Node<T,E>::ode_aux(It it) const {
  for (size_t i = 0; i < individual.aux_size(); i++) {
    *it++ = util::as_iterator_scalar<It>(individual.aux(i));
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
