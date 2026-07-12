#include <plant/models/tf24f_strategy.h>
#include <plant/individual.h>
#include <plant/patch.h>
#include <odelia/ode_solver.hpp>
#include <XAD/XAD.hpp>

namespace plant {

// The base constructor sets collect_all_auxiliary and runs the base
// refresh_indices(); we re-run our own so the appended state slot is registered.
template <class S>
TF24f_Strategy_<S>::TF24f_Strategy_() {
  this->name = "TF24f";
  refresh_indices();
}

// Build on the base index maps, then register the extra tracked state.
template <class S>
void TF24f_Strategy_<S>::refresh_indices() {
  TF24_Strategy_<S>::refresh_indices();
  const int idx = static_cast<int>(TF24_Strategy_<S>::state_size());
  this->state_index["opt_root_psi_state"] = idx;
  state_idx_opt_root_psi_state = idx;
}

// Reuse TF24's rates for the five shared states; the tracked-state value feeds
// solve_leaf() via tracked_root_psi_ (a plant-local double), and the resulting
// profit gradient comes back in dprofit_dpsi_, which becomes the tracked state's
// rate (gradient ascent). That rate is a passive constant on the reverse tape.
template <class S>
void TF24f_Strategy_<S>::compute_rates(const TF24_Environment_<S>& environment,
                                       Internals_<S>& vars) {
  if (!util::is_finite(k_acclim) || k_acclim < 0.0) {
    util::stop("TF24f: k_acclim must be finite and >= 0 (got " +
               util::to_string(k_acclim) + ")");
  }
  tracked_root_psi_ = ad_value(vars.state(state_idx_opt_root_psi_state));
  TF24_Strategy_<S>::compute_rates(environment, vars);
  vars.set_rate(state_idx_opt_root_psi_state, S(k_acclim * dprofit_dpsi_));
}

// Track instead of optimise: evaluate the leaf at the tracked collar psi and supply
// d(profit)/d(psi) for the gradient-ascent rate.
template <class S>
void TF24f_Strategy_<S>::solve_leaf() {
  if (initializing_) {
    this->leaf.find_root_collar_psi();
    return;
  }
  if (use_ad_gradient) {
    // Exact gradient (default, #527): forward-mode AD + IFT + analytic spline slopes.
    this->leaf.evaluate_root_collar_psi(tracked_root_psi_);
    const double used = -this->leaf.root_collar_psi_;
    dprofit_dpsi_ = this->leaf.dprofit_droot_collar_psi(used);
    this->leaf.evaluate_root_collar_psi(used);  // restore operating-point outputs
  } else {
    // Centred finite-difference fallback (#526), sharing one prepare_collar_solve.
    const double h = psi_fd_step;
    if (!util::is_finite(h) || h <= 0.0) {
      util::stop("TF24f: psi_fd_step must be finite and > 0 (got " +
                 util::to_string(h) + ")");
    }
    double bound_a, bound_b;
    if (!this->leaf.prepare_collar_solve(bound_a, bound_b)) {
      dprofit_dpsi_ = 0.0;
      return;
    }
    const double used = std::min(std::max(tracked_root_psi_, bound_a), bound_b);
    const double p_plus  = this->leaf.profit_at_collar_psi(used + h, bound_a, bound_b);
    const double p_minus = this->leaf.profit_at_collar_psi(used - h, bound_a, bound_b);
    dprofit_dpsi_ = (p_plus - p_minus) / (2.0 * h);
    this->leaf.profit_at_collar_psi(used, bound_a, bound_b);  // restore operating point
  }
}

// Seed the tracked state at its optimum: run the base optimiser once (via the
// initializing_ flag) and store the resulting collar psi as the initial state.
template <class S>
void TF24f_Strategy_<S>::set_initial_states(const TF24_Environment_<S>& environment,
                                            Internals_<S>& vars) {
  initializing_ = true;
  this->net_mass_production_dt(environment, vars);
  initializing_ = false;
  vars.set_state(state_idx_opt_root_psi_state, S(-this->leaf.root_collar_psi_));
}

TF24f_Strategy::ptr make_strategy_ptr(TF24f_Strategy s) {
  s.prepare_strategy();
  return std::make_shared<TF24f_Strategy>(s);
}

template class TF24f_Strategy_<double>;

using ad_reverse =
    odelia::ode::Solver<Patch<TF24f_Strategy_<double>, TF24_Environment>>::active_scalar;
template class TF24f_Strategy_<ad_reverse>;
}
