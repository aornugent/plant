// -*-c++-*-
#ifndef PLANT_PLANT_TF24F_STRATEGY_H_
#define PLANT_PLANT_TF24F_STRATEGY_H_

#include <plant/models/tf24_strategy.h>

namespace plant {

// TF24f ("f" for fast / forecasting): a variant of TF24 that will let the
// optimal leaf hydraulic state chase its optimum via an extra ODE state
// (gradient-ascent), instead of re-solving the nested leaf optimisation /
// root-collar root-find from scratch at every step (issue #525). It inherits
// TF24_Strategy_<S> and reuses TF24_Environment_<S> + TF24_Pars_<S>; only the
// leaf-solve hook and the extra state are overridden. Templated on the same S
// as its base; S = double is the production path (the `TF24f_Strategy` alias
// below). The tracked state rides the double Leaf, so its parameter sensitivity
// (like TF24's) reaches the tape via supplied_derivative, not by templating.
//
// Adds one extra ODE state, opt_root_psi_state, holding the tracked optimal
// root-collar water potential. Instead of TF24's per-step golden-section
// optimisation, solve_leaf() evaluates the leaf at the tracked state and the
// state relaxes toward its optimum by gradient ascent on carbon profit
// (dpsi/dt = k_acclim * d(profit)/d(psi)); the state is seeded at its optimum at
// birth. The five states shared with TF24 are computed by the inherited
// compute_rates, so TF24f does NOT reproduce TF24 bit-for-bit -- it tracks the
// optimum with a (gain-dependent) lag and is faster.
template <class S = double>
class TF24f_Strategy_ : public TF24_Strategy_<S> {
public:
  using environment_type = TF24_Environment_<S>;
  using value_type = S;
  typedef std::shared_ptr<TF24f_Strategy_<S>> ptr;
  TF24f_Strategy_();

  // TF24's five states + the tracked root-collar psi (appended last so the
  // inherited state indices 0..4 are unchanged). state_size()/state_names() are
  // static and resolved on the concrete type by Individual<TF24f, ...>.
  static size_t state_size() { return TF24_Strategy_<S>::state_size() + 1; }
  static std::vector<std::string> state_names() {
    std::vector<std::string> ret = TF24_Strategy_<S>::state_names();
    ret.push_back("opt_root_psi_state");
    return ret;
  }

  // Re-register indices including the new state slot. Base refresh_indices()
  // calls the *static* (base) state_names(), so it would otherwise miss the
  // appended state; we call it then add the extra slot.
  void refresh_indices();

  void compute_rates(const environment_type& environment, Internals_<S>& vars);

  // Override the leaf solve: instead of optimising the root-collar psi, evaluate
  // the leaf at the tracked state and finite-difference the profit gradient
  // (left in dprofit_dpsi_ for compute_rates to turn into the state's rate).
  void solve_leaf();

  // Leaf-seam collar-psi channel (Stage D). The leaf operates at the tracked
  // collar-psi ODE state, so the seam must add d(profit)/d(psi) * d(psi)/d(theta):
  // the input is the ACTIVE tracked state (its derivative is on the tape) and the
  // partial is d(profit)/d(psi) (= dprofit_dpsi_, the acclimation gradient).
  S* seam_collar_psi_input() override { return &tracked_root_psi_active_; }
  double seam_collar_psi_partial() const override { return dprofit_dpsi_; }

  // Seed the tracked state at its optimum for a newly introduced individual, so
  // gradient ascent starts at the optimum (no birth transient / no climb from 0,
  // which otherwise lets shaded recruits escape suppression). Runs the base
  // optimiser once via the initializing_ flag below.
  void set_initial_states(const environment_type& environment, Internals_<S>& vars);

  // Acclimation gain k in  dpsi/dt = k * d(profit)/d(psi). Exposed to R so the
  // stiffness / accuracy-vs-speed k-sweep (#525) can be driven without a rebuild;
  // large k recovers the quasi-steady-state (TF24) optimum.
  double k_acclim = 1.0;
  // Finite-difference step (positive magnitude, MPa) for d(profit)/d(psi); used
  // only when use_ad_gradient is false.
  double psi_fd_step = 1e-3;
  // Use the exact AD/IFT gradient (Leaf::dprofit_droot_collar_psi) instead of the
  // finite difference. Exposed to R for A/B comparison (#527).
  bool use_ad_gradient = true;

  // Cached slot for the tracked state, resolved in refresh_indices().
  int state_idx_opt_root_psi_state = -1;

  // Channel between solve_leaf() (writes) and compute_rates() (reads). Default is
  // finite so solve_leaf is safe before the first compute_rates (e.g. during
  // establishment_probability), where the tracked state has not been read yet.
  // Off the S rate path (the Leaf is double), so held as double.
  double tracked_root_psi_ = 0.0;
  double dprofit_dpsi_ = 0.0;
  // The tracked collar-psi state kept as an ACTIVE scalar (value == the double
  // tracked_root_psi_) so the leaf seam can inject d(profit)/d(psi) onto its
  // tape slot. Set each compute_rates before the base rate computation runs the
  // seam. Default 0 (passive) so a stray seam read before the first set drops out.
  S tracked_root_psi_active_ = 0.0;
  // When true, solve_leaf() runs the base optimiser (find_root_collar_psi) rather
  // than the tracked evaluation; used by set_initial_states to read the optimum.
  bool initializing_ = false;
};

using TF24f_Strategy = TF24f_Strategy_<double>;

}

#endif
