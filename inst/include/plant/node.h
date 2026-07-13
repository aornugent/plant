// -*-c++-*-
#ifndef NODE
#define NODE

#include <plant/environment.h>
#include <plant/gradient.h>
#include <odelia/ode_interface.hpp>
#include <optional>
#include <limits> // std::numeric_limits
#include <type_traits> // std::is_same_v

namespace plant {

template <typename T, typename E>
class Node {
public:
  typedef T        strategy_type;
  typedef E        environment_type;
  typedef Individual<T,E> individual_type;
  typedef typename strategy_type::ptr strategy_type_ptr;
  // The scalar the node's demographic state (log_density, offspring) carries,
  // taken from the individual it wraps.
  using value_type = typename individual_type::value_type;
  Node(strategy_type_ptr s);

  void compute_rates(const environment_type& environment, double pr_patch_survival);
  void compute_initial_conditions(const environment_type& environment, double pr_patch_survival, double birth_rate);

  // Wrapper to growth_rate_gradient for testing
  double r_growth_rate_gradient(const environment_type& environment);

  value_type height() const {return individual.state(HEIGHT_INDEX);}
  value_type compute_competition(double z) const;
  double fecundity() const {return offspring_produced_survival_weighted;}

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
  double weighted_fecundity(double S_D) const {
    return offspring_produced_survival_weighted * patch_density_at_birth * S_D;
  }

  // Unfortunate, but need a get_ here because of name shadowing...
  value_type get_log_density() const {return log_density;}
  // exp(log_density); can overflow to +Inf when the SCM density equation runs
  // away (see Patch::check_finite_node_densities).
  value_type get_density() const {return density;}
  void set_log_density(value_type x) {
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
  // This is the gradient of growth rate with respect to height. The result
  // carries value_type: the finite-difference stencil (abscissa/step) is a
  // double primitive, but the growth rate is evaluated on the active parameters,
  // so its parameter-derivative flows into log_density_dt (the density-transport
  // term) instead of being dropped (§0.5).
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

  // NOTE: This must be called *after* compute_rates, but given we
  // need mortality_dt() that's always going to be the case.
  log_density_dt =
    - growth_rate_gradient(environment)
    - individual.rate(MORTALITY_INDEX);
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

  const value_type pr_estab = individual.establishment_probability(environment);
  individual.set_state("mortality", -log(pr_estab));
  const value_type g = individual.rate(HEIGHT_INDEX);
  // NOTE: log(0.0) -> -Inf, which should behave fine.
  // Collapse both arms to value_type: at an active scalar the two log(...)
  // expressions have different XAD expression types, which a raw ?: cannot
  // reconcile. This is the log-density birth kink (recorded in the manifest).
  set_log_density(g > 0 ? value_type(log(birth_rate * pr_estab / g))
                        : value_type(log(value_type(0.0))));

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
  const Control& control = individual.control();
  const double eps = control.node_gradient_eps;

  if constexpr (std::is_same_v<value_type, double>) {
    // Production / R-facing path: unchanged. Finite-difference the growth rate
    // on a thread-local scratch (copy assignment reuses the vector storage, so
    // steady-state calls don't allocate), reusing the already-computed rate as
    // fx and honouring the direction / Richardson controls.
    thread_local std::optional<individual_type> scratch;
    if (scratch.has_value()) { *scratch = individual; }
    else                     { scratch.emplace(individual); }
    individual_type& p = *scratch;
    auto fun = [&] (double h) -> double {
      return p.growth_rate_given_height(h, environment);
    };
    const double h0 = individual.state(HEIGHT_INDEX);
    if (control.node_gradient_richardson) {
      return util::gradient_richardson(fun, h0, eps,
                                       control.node_gradient_richardson_depth);
    } else {
      return util::gradient_fd(fun, h0, eps, individual.rate(HEIGHT_INDEX),
                               control.node_gradient_direction);
    }
  } else {
    // Active pass: replicate the double path's stencil EXACTLY (same direction /
    // Richardson controls, same eps, same frozen abscissa h0 = value(height)) so
    // the growth-gradient VALUE is bit-identical to the double replay -- the
    // active forward pass must reproduce the double trajectory, and dg/dh feeds
    // log_density -> density -> competition, so a different stencil here silently
    // forks the trajectory once cohorts shade each other. The ONE change: source
    // fx from the scratch (fun(h0)) rather than the real node's cached rate. The
    // value is the same (p is a copy, so p.growth_rate_given_height(h0) equals
    // individual.rate(HEIGHT_INDEX)), but both stencil points now flow through
    // the SAME tape subgraph, so their difference is a clean parameter-derivative.
    // Reusing the real node's rate as fx instead mixes two subgraphs whose O(eps)
    // value match hides an O(1) derivative mismatch that /eps amplifies to garbage.
    // Fresh copy-CONSTRUCT records `p = individual` onto the live tape so the
    // perturbed rates link back to the real active parameters.
    // NOTE (TF24 retrofit trigger): this differentiates the growth rate on the
    // outer tape, which is only valid while g is retapeable. When TF24 makes g
    // re-run a non-retapeable leaf optimiser at the perturbed height, dg/dh must
    // instead be computed off-tape and injected via odelia::supplied_derivative
    // (Kind B). See docs/ad-implementation.md.
    individual_type p = individual;
    auto fun = [&] (double h) -> value_type {
      return p.growth_rate_given_height(value_type(h), environment);
    };
    const double h0 = xad::value(individual.state(HEIGHT_INDEX));
    value_type dgdh = control.node_gradient_richardson
      ? util::gradient_richardson(fun, h0, eps, control.node_gradient_richardson_depth)
      : util::gradient_fd(fun, h0, eps, fun(h0), control.node_gradient_direction);
    // Deferred (candidate B): drop the parameter-derivative of dg/dh, keeping its
    // value. Differentiating the FD stencil on the outer tape is unreliable -- for
    // a suppressed cohort the backward stencil straddles the size_dt growth clamp,
    // so d(dg/dh)/dtheta = (dg/dtheta(h0) - dg/dtheta(h0-eps))/eps picks up the
    // kink and /eps amplifies it (measured: the two-cohort resident gradient blows
    // to ~2.3x FD, while dropping the term lands a clean, consistent 0.96x across
    // every seeded parameter -- the residual IS this dropped d2g/dh.dtheta). The
    // correct fix is to compute dg/dh off-tape (double) and inject its parameter
    // partials via odelia::supplied_derivative (Kind B, ad-implementation §15) --
    // the same seam TF24 needs when g re-runs a non-retapeable leaf optimiser.
    // Until then the single-cohort gradient is exact and the multi-cohort one
    // carries this bounded ~4% underestimate.
    return value_type(xad::value(dgdh));
  }
}

// Wrapper to growth_rate_gradient for testing
template <typename T, typename E>
double Node<T,E>::r_growth_rate_gradient(const environment_type& environment) {
  // We need to compute the physiological variables here, first, so
  // that reusing intervals works as expected.  This would ordinarily
  // be taken care of because of the calling order of
  // compute_rates / growth_rate_gradient.
  individual.compute_rates(environment);
  return xad::value(growth_rate_gradient(environment));  // R-facing: double only
}

template <typename T, typename E>
typename Node<T,E>::value_type
Node<T,E>::compute_competition(double height_) const {
  return density * individual.compute_competition(height_);
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
Node<T,E> make_node(typename Node<T,E>::strategy_type s) {
  return Node<T,E>(make_strategy_ptr(s));
}

}

#endif /* NODE */
