// -*-c++-*-
#ifndef PLANT_PLANT_TF24F_STRATEGY_H_
#define PLANT_PLANT_TF24F_STRATEGY_H_

#include <plant/models/tf24_strategy.h>

namespace plant {

// TF24f ("f" for fast / forecasting): a variant of TF24 that lets the optimal leaf
// hydraulic state chase its optimum via an extra ODE state (gradient-ascent),
// instead of re-solving the nested leaf optimisation from scratch at every step
// (issue #525). It inherits TF24_Strategy_<S> and reuses TF24_Environment + TF24_Pars;
// only the leaf-solve hook and the extra state are overridden.
//
// Adds one extra ODE state, opt_root_psi_state, holding the tracked optimal
// root-collar water potential; the state relaxes toward its optimum by gradient
// ascent on carbon profit. The tracked leaf state is a plant-local double, so the
// state's rate is a passive constant on the reverse tape -- TF24f's coupled resident
// gradient at long horizons is out of v1 (stiffness, design; gated by the entry).
template <class S = double>
class TF24f_Strategy_ : public TF24_Strategy_<S> {
public:
  using value_type = S;
  typedef std::shared_ptr<TF24f_Strategy_> ptr;
  TF24f_Strategy_();

  // TF24's five states + the tracked root-collar psi (appended last so the
  // inherited state indices 0..4 are unchanged).
  static size_t state_size() { return TF24_Strategy_<S>::state_size() + 1; }
  static std::vector<std::string> state_names() {
    std::vector<std::string> ret = TF24_Strategy_<S>::state_names();
    ret.push_back("opt_root_psi_state");
    return ret;
  }

  void refresh_indices();

  void compute_rates(const TF24_Environment_<S>& environment, Internals_<S>& vars);

  // Override the leaf solve: evaluate the leaf at the tracked state and leave the
  // profit gradient in dprofit_dpsi_ for compute_rates to turn into the state rate.
  void solve_leaf();

  // Seed the tracked state at its optimum for a newly introduced individual.
  void set_initial_states(const TF24_Environment_<S>& environment, Internals_<S>& vars);

  // odelia differentiable-System handles ------------------------------------
  template <class S2> using rebind = TF24f_Strategy_<S2>;

  template <class S2>
  TF24f_Strategy_<S2> rebind_from() const {
    TF24f_Strategy_<S2> out;
    out.pars = this->pars.template rebind<S2>();
    out.control = this->control;
    out.name = this->name;
    out.birth_rate_x = this->birth_rate_x;
    out.birth_rate_y = this->birth_rate_y;
    out.is_variable_birth_rate = this->is_variable_birth_rate;
    out.collect_all_auxiliary = this->collect_all_auxiliary;
    out.k_acclim = k_acclim;
    out.psi_fd_step = psi_fd_step;
    out.use_ad_gradient = use_ad_gradient;
    out.refresh_indices();
    return out;
  }

  // Acclimation gain k in  dpsi/dt = k * d(profit)/d(psi). Exposed to R (#525).
  double k_acclim = 1.0;
  // Finite-difference step (MPa) for d(profit)/d(psi) when use_ad_gradient is false.
  double psi_fd_step = 1e-3;
  // Use the exact AD/IFT gradient (Leaf::dprofit_droot_collar_psi) instead of FD (#527).
  bool use_ad_gradient = true;

  // Cached slot for the tracked state, resolved in refresh_indices().
  int state_idx_opt_root_psi_state = -1;

  // Channel between solve_leaf() (writes) and compute_rates() (reads).
  double tracked_root_psi_ = 0.0;
  double dprofit_dpsi_ = 0.0;
  // When true, solve_leaf() runs the base optimiser (find_root_collar_psi) rather
  // than the tracked evaluation; used by set_initial_states to read the optimum.
  bool initializing_ = false;
};

using TF24f_Strategy = TF24f_Strategy_<double>;

TF24f_Strategy::ptr make_strategy_ptr(TF24f_Strategy s);

template <class S>
typename TF24f_Strategy_<S>::ptr make_strategy_ptr(TF24f_Strategy_<S> s) {
  s.prepare_strategy();
  return std::make_shared<TF24f_Strategy_<S>>(s);
}

}

#endif
