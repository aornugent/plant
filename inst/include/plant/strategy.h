// -*-c++-*-
#ifndef PLANT_PLANT_STRATEGY_H_
#define PLANT_PLANT_STRATEGY_H_

#include <memory>
#include <plant/control.h>
#include <plant/internals.h>
#include <plant/uniroot.h>
#include <RcppCommon.h> // NA_REAL
#include <plant/uniroot.h>
#include <plant/extrinsic_drivers.h>


namespace plant {

template <typename E>
class Strategy {
public:
  typedef E             environment_type;
  typedef std::shared_ptr<Strategy> ptr;

  // The scalar the strategy's physiology runs at. double is the production path;
  // a strategy templated on S overrides this so a trait derivative flows through
  // Individual/Node/Species/Patch, which read it as value_type.
  using value_type = double;

  // update this when the length of state_names changes
  static size_t state_size ();
  // update this when the length of aux_names changes
  size_t aux_size ();

  static std::vector<std::string> state_names();

  std::vector<std::string> aux_names();

  // TODO(#483) : expose this so can access state_names directly
  // In previous attempt couldn't get it to run
  // static std::vector<std::string> state_names() { return strategy_type::state_names(); }
  // the index of variables in the internals extra vector
  std::map<std::string, int> state_index; 
  std::map<std::string, int> aux_index;

  // birth rate spline control points for each species
  // default is constant birth_rate of 1.0
  std::vector<double> birth_rate_x;
  std::vector<double> birth_rate_y = {1.0};
  // whether the spline for each species should be constant fn or not (extrapolation on/off)
  bool is_variable_birth_rate = false;

  bool collect_all_auxiliary;

  // Copy the scalar-independent configuration (birth-rate spline, flags, control,
  // name, drivers, size_0) from another strategy, whatever scalar it runs at. The
  // differentiable parameters and precomputed state are the concrete strategy's own
  // concern -- its rebind_from() carries the former via field_ptrs() and leaves the
  // latter for prepare_strategy(). Used to lift a configured double strategy onto an
  // active scalar for a gradient (the odelia System rebind_from contract).
  template <class Other>
  void copy_config_from(const Other& o) {
    birth_rate_x           = o.birth_rate_x;
    birth_rate_y           = o.birth_rate_y;
    is_variable_birth_rate = o.is_variable_birth_rate;
    collect_all_auxiliary  = o.collect_all_auxiliary;
    control                = o.control;
    name                   = o.name;
    extrinsic_drivers      = o.extrinsic_drivers;
    size_0                 = o.size_0;
  }

  void refresh_indices();

  double competition_effect(double size) const;

  double competition_effect_state(Internals& vars);

  void compute_rates(const environment_type& environment, Internals& vars);

  void update_dependent_aux(const int index, Internals& vars);

  // Seed strategy-specific initial ODE states for a newly introduced individual,
  // given its birth environment (called once from Node::compute_initial_conditions
  // before the first compute_rates). Default no-op; strategies that carry an
  // acclimating/tracked state (e.g. TF24f, #525) override this to initialise it at
  // its optimum so there is no birth transient. Resolved on the concrete strategy
  // type by Individual<T,E>, so overriding it here is not required to be virtual.
  // Templated on the Internals scalar so the concrete strategy's value_type
  // (double or active) flows through this inherited no-op; a strategy that
  // carries a tracked state overrides it at its own scalar.
  template <typename V>
  void set_initial_states(const environment_type& environment, Internals_<V>& vars) {
    (void)environment;
    (void)vars;
  }

  double net_mass_production_dt(const environment_type& environment,
                                double size, double competition_effect_);

  double establishment_probability(const environment_type& environment);

  double fecundity_dt(double net_mass_production_dt,
                      double fraction_allocation_reproduction) const;

  double mortality_dt(double productivity_area, double cumulative_mortality) const;

  double compute_competition(double z, double size) const;

  double initial_size(void) const;

  double size_0;

  // Every Strategy needs a set of Control objects -- these govern
  // things to do with how numerical calculations are performed,
  // rather than the biological control that this class has.
  Control control;

  std::string name;

  ExtrinsicDrivers extrinsic_drivers;
};

// Prepare a strategy and hand back a shared pointer. Generic over the concrete
// strategy type (double or scalar-templated), so each strategy needs no bespoke
// overload -- Individual/Node/Species construct through this one entry point.
template <class Strat>
typename Strat::ptr make_strategy_ptr(Strat s) {
  s.prepare_strategy();
  return std::make_shared<Strat>(s);
}

// The whole double->active rebind mechanic, in one place: a fresh strategy of type
// To carrying `from`'s scalar-independent config (copy_config_from) and its
// differentiable parameters widened through field_ptrs() (a new AD field is carried
// automatically). Precomputed state is left for prepare_strategy() to rebuild from
// the seeded parameters (the reset-timing contract). Each concrete strategy's
// rebind_from() -- the odelia System hook -- is a one-line call to this, so the
// boilerplate lives here, not once per model.
template <class To, class From>
To rebind_strategy_fields(From& from) {
  To a;
  a.copy_config_from(from);
  auto dp = from.field_ptrs();
  auto ap = a.field_ptrs();
  for (std::size_t i = 0; i < dp.size(); ++i)
    *ap[i] = typename To::value_type(*dp[i]);
  return a;
}

// One token that makes a scalar-templated strategy differentiable: the two odelia
// System hooks every strategy needs identically -- the `rebind` alias (name the same
// strategy at another scalar) and `rebind_from` (copy this configured strategy onto
// that scalar, via rebind_strategy_fields: scalar-independent config + the field_ptrs()
// parameters, precomputed state rebuilt later by prepare_strategy()). Invoke in the
// class body with the strategy's own template name; it composes with the AD_FIELDS
// X-macro (which defines field_ptrs). A strategy carrying scalar-independent config
// beyond field_ptrs (e.g. TF24f's acclimation settings) writes rebind_from by hand.
#define PLANT_DIFFERENTIABLE(STRATEGY_TMPL)                                    \
  template <class U> using rebind = STRATEGY_TMPL<U>;                          \
  template <class S2> STRATEGY_TMPL<S2> rebind_from() {                        \
    return plant::rebind_strategy_fields<STRATEGY_TMPL<S2>>(*this);            \
  }

}

#endif
