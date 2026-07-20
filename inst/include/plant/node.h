// -*-c++-*-
#ifndef NODE
#define NODE

#include <plant/environment.h>
#include <plant/gradient.h>
#include <odelia/ode_interface.hpp>
#include <odelia/ode_util.hpp> // odelia::util::diagnostic (intent-named AD strip)
#include <optional>
#include <limits> // std::numeric_limits
#include <type_traits> // std::is_same_v

namespace plant {

// Detect a strategy that SUPPORTS the geometric mass chart (neighbour-secant
// transport, odelia::log_density_rate), via a nested `geometric_transport` marker
// type. The mass chart is the stable, differentiable transport for K93; the FD
// secant is numerically unstable for FF16/TF24 (it drives cohort density to
// overflow), so they do not declare it and never take that arm regardless of the
// Control flag. For supporting strategies the flag node_geometric_compression
// (default on) selects it, and can still turn it off per run.
template <typename S2, typename = void>
struct strategy_supports_geometric_transport : std::false_type {};
template <typename S2>
struct strategy_supports_geometric_transport<S2, std::void_t<typename S2::geometric_transport>>
    : std::true_type {};

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
  // Growth rate dh/dt (the characteristic velocity g). Exposed so Species can
  // build the geometric compression -∂ₓg from neighbouring cohorts (Tier 1).
  value_type growth_rate() const {return individual.rate(HEIGHT_INDEX);}
  value_type compute_competition(double z) const;

  value_type mortality_rate() const {return individual.rate(MORTALITY_INDEX);}
  // Tier 1: Species computes the log-density transport rate from the geometric
  // compression across neighbouring cohorts (it needs the neighbour list) and
  // sets it here. See Species::compute_rates.
  void set_log_density_dt(value_type v) {log_density_dt = v;}
  value_type fecundity() const {return offspring_produced_survival_weighted;}

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
  value_type weighted_fecundity(value_type S_D) const {
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

  // Transport-log-mass chart (geometric strategies with the flag on). The
  // ODE-integrated slot carries lambda = log_density + log(cohort_spacing),
  // whose rate is just -mortality (the compression cancels; see
  // odelia::log_mass_from_log_density). log_density/density are a read-side
  // view Species reconstructs from lambda and the shared spacing. Whether this
  // node is on the chart is fixed by the strategy marker and the control flag.
  bool on_mass_chart() const {
    if constexpr (strategy_supports_geometric_transport<strategy_type>::value)
      return individual.control().node_geometric_compression;
    else
      return false;
  }
  value_type get_log_mass() const {return log_mass_;}
  void set_log_mass(value_type x) {log_mass_ = x;}

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
  value_type log_mass_;  // transported quantity on the mass chart (see on_mass_chart)
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
    log_mass_(-std::numeric_limits<double>::infinity()),
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
  bool geometric = false;
  if constexpr (strategy_supports_geometric_transport<strategy_type>::value)
    geometric = individual.control().node_geometric_compression;
  if (geometric) {
    // Species::compute_rates sets -dg/dx - mortality from cohort neighbour
    // spacings for every node with a neighbour. Seed the mortality term here so
    // a lone cohort (which Species leaves untouched) still decays; Species
    // overwrites this for nodes with >= 2 cohorts.
    log_density_dt = - individual.rate(MORTALITY_INDEX);
  } else {
    // Default path (flag off, and always for FF16/TF24/TF24f): the
    // finite-difference stencil. Those strategies' differentiated metrics do
    // not route through ∂ₓg (docs §15), so Tier 1 offers them nothing; with the
    // flag off the published models are bit-for-bit unchanged.
    log_density_dt =
      - growth_rate_gradient(environment)
      - individual.rate(MORTALITY_INDEX);
  }
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
  // Default the transported log mass to the birth log density (the lone-cohort
  // dx=1 convention). Species::seed_newborn_log_mass overwrites this with the
  // spacing-adjusted lambda once the node has a neighbour; off the mass chart
  // log_mass_ is unused. This keeps a bare node's ode_state == log_density.
  log_mass_ = log_density;

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
  // dg/dh (the McKendrick density-transport term). Its VALUE always comes from the
  // finite-difference stencil, on both the double and active passes. Two reasons:
  //   1. STABILITY. The one-sided FD stencil is the UPWIND discretisation of the
  //      advection (transport) term d(log_density)/dt = -dg/dh - mortality. Upwinding
  //      is the standard stabilisation for hyperbolic transport on a coarse grid; the
  //      exact analytic dg/dh is the centred/exact scheme, which is numerically
  //      unstable here -- fed into log_density_dt it drives the density transport out
  //      of bounds once cohorts shade each other (verified: the analytic trajectory
  //      only stays bounded with a growth clamp smoothed to eps ~ 5e-2, a ~6% change
  //      to the K93 demography). So the FD stencil DEFINES the production trajectory
  //      and it is not merely a clamp-corner workaround (plant#39).
  //   2. CONSISTENCY. The active forward pass must reproduce the double trajectory
  //      (a different dg/dh value forks it once cohorts shade each other), so the
  //      active pass keeps the same FD value.
  const Control& control = individual.control();
  const double eps = control.node_gradient_eps;

  if constexpr (std::is_same_v<value_type, double>) {
    // Production / R-facing: FD on a reused thread_local scratch (no per-call alloc).
    thread_local std::optional<individual_type> scratch;
    if (scratch.has_value()) { *scratch = individual; }
    else                     { scratch.emplace(individual); }
    individual_type& p = *scratch;
    auto fun = [&] (double h) -> double { return p.growth_rate_given_height(h, environment); };
    const double h0 = individual.state(HEIGHT_INDEX);
    if (control.node_gradient_richardson)
      return util::gradient_richardson(fun, h0, eps, control.node_gradient_richardson_depth);
    return util::gradient_fd(fun, h0, eps, individual.rate(HEIGHT_INDEX),
                             control.node_gradient_direction);
  } else {
    // Active pass. The FD VALUE (matching the double trajectory) is computed in
    // double on a fresh copy (value only -- its stencil derivative is discarded).
    // The parameter-derivative of the transport term is deliberately dropped here:
    // strategies whose census gradient needs it (K93) take the geometric mass
    // chart instead (Species::compute_rates), whose derivative is exact and well
    // conditioned; the other strategies' differentiated metrics (mutant fitness)
    // carry no density factor, so this term is off that graph (docs §15). This is
    // why the old forward-over-reverse dg/dh injection (plant#39) is gone -- the
    // mass chart superseded it.
    individual_type p = individual;
    auto fun_d = [&] (double h) -> double {
      return xad::value(p.growth_rate_given_height(value_type(h), environment));
    };
    const double h0 = xad::value(individual.state(HEIGHT_INDEX));
    const double fd_value = control.node_gradient_richardson
      ? util::gradient_richardson(fun_d, h0, eps, control.node_gradient_richardson_depth)
      : util::gradient_fd(fun_d, h0, eps, fun_d(h0), control.node_gradient_direction);
    return value_type(fd_value);
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
  return odelia::util::diagnostic(growth_rate_gradient(environment));  // R-facing: double only
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
  // Mass chart: the integrated slot is lambda (log mass); log_density/density
  // are reconstructed later (Species::reconstruct_densities, before the field is
  // read). Off the chart, the slot is log_density directly.
  if (on_mass_chart()) { log_mass_ = *it++; }
  else                 { set_log_density(*it++); }
  return it;
}
template <typename T, typename E>
template <typename It>
It Node<T,E>::ode_state(It it) const {
  for (size_t i = 0; i < individual.ode_size(); i++) {
    *it++ = individual.state(i);
  }
  *it++ = offspring_produced_survival_weighted;
  *it++ = on_mass_chart() ? log_mass_ : log_density;
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
