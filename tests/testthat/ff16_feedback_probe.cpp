// [[Rcpp::plugins(cpp20)]]
// Isolate ONE compute_rates + field assembly (the instantaneous self-shading
// feedback) from the SCM trajectory. Procedure: run a short FF16 SCM once in
// double to get a realistic multi-cohort ODE state; freeze that state to
// trait-INDEPENDENT doubles; then build an SCM at scalar S with lma ACTIVE, load
// the SAME frozen state, recompute the environment (field) + rates, and sum the
// cohort growth (height) rates. Because the state is fixed, d(sum g)/d(lma) is
// PURELY the instantaneous feedback derivative (through area_leaf and the field),
// with the trajectory removed. Compare reverse-AD to a central FD on the same
// frozen state.
//
// RESULT (2026-07-19): AD == FD to 6 digits (Node::ode_size()==7 per cohort, so
// height rates stride by 7 -- an interim stride bug of 5 produced spurious
// mismatches). So ONE compute_rates + field assembly is faithful, and there is NO
// detached edge: the reverse AD correctly computes the frozen-schedule gradient.
// The FF16 R0 gap vs the adaptive model is purely L1 schedule sensitivity.
#include <Rcpp.h>
#include <plant/models/ff16_strategy.h>
#include <plant/models/ff16_environment.h>
#include <plant/individual.h>
#include <plant/parameters.h>
#include <plant/control.h>
#include <plant/scm.h>
#include <odelia/gradient.hpp>
#include <vector>

using namespace plant;

// Build a patch at scalar S with lma active, with `n_cohorts` introduced, load the
// frozen state, recompute env+rates, return sum of per-cohort HEIGHT rates.
template <class S>
static S instant_rate(S lma, std::size_t n_cohorts,
                      const std::vector<double>& frozen_state, double t, double L) {
  FF16_Strategy_<S> strat;
  strat.pars.lma = lma;                    // active input flows through the params
  strat.is_variable_birth_rate = false;
  strat.birth_rate_y = {20.0};
  Parameters<FF16_Strategy_<S>, FF16_Environment_<S>> p;
  p.strategies.push_back(strat);
  p.max_patch_lifetime = L;
  p.validate();
  FF16_Environment_<S> env;
  Control ctrl;
  SCM<FF16_Strategy_<S>, FF16_Environment_<S>> scm(p, env, ctrl);
  scm.reset();
  auto& patch = scm.get_system_ref();
  for (std::size_t k = 0; k < n_cohorts; ++k) patch.introduce_new_node(0);

  // Load the frozen (double) state as S, then set it -- this recomputes the
  // environment (assembles the field from the frozen heights + active area_leaf)
  // and the rates.
  std::vector<S> st(frozen_state.size());
  for (std::size_t i = 0; i < frozen_state.size(); ++i) st[i] = S(frozen_state[i]);
  patch.set_ode_state(st.begin(), t);

  // Sum the HEIGHT rate of each cohort. Node::ode_size() = state_size()+2 = 7 per
  // cohort (5 individual states + log_density + offspring), so the flat rate vector
  // strides by 7; height is slot 0.
  std::vector<S> rates(patch.ode_size());
  patch.ode_rates(rates.begin());
  S g = 0.0;
  for (std::size_t k = 0; k < n_cohorts; ++k) g += rates[k * 7 + 0];
  return g;
}

// [[Rcpp::export]]
Rcpp::List ff16_feedback_probe(double lma = 0.1978791, double fd_rel = 3e-4, double L = 10.0) {
  using RevS = xad::adj<double>::active_type;

  // 1. Short double run to get a realistic multi-cohort state; freeze it.
  std::vector<double> frozen_state; std::size_t n_cohorts = 0; double t = 0.0;
  {
    FF16_Strategy_<double> strat; strat.pars.lma = lma;
    strat.is_variable_birth_rate = false; strat.birth_rate_y = {20.0};
    Parameters<FF16_Strategy_<double>, FF16_Environment_<double>> p;
    p.strategies.push_back(strat); p.max_patch_lifetime = L; p.validate();
    FF16_Environment_<double> env; Control ctrl;
    SCM<FF16_Strategy_<double>, FF16_Environment_<double>> scm(p, env, ctrl);
    scm.run();
    auto& patch = scm.get_system_ref();
    n_cohorts = patch.at_species(0).size();
    frozen_state.resize(patch.ode_size());
    patch.ode_state(frozen_state.begin());
    t = patch.ode_time();
  }

  // 2. Reverse-AD of the instantaneous rate on the frozen state.
  double grad_ad = 0.0, value = 0.0;
  {
    xad::Tape<double> tape;
    RevS L_in = lma; tape.registerInput(L_in); tape.newRecording();
    RevS out = instant_rate<RevS>(L_in, n_cohorts, frozen_state, t, L);
    tape.registerOutput(out); xad::derivative(out) = 1.0; tape.computeAdjoints();
    value = xad::value(out); grad_ad = xad::derivative(L_in);
  }

  // 3. Central FD on the SAME frozen state (double).
  double d = fd_rel * (lma + 1e-3);
  double fd = (instant_rate<double>(lma + d, n_cohorts, frozen_state, t, L) -
               instant_rate<double>(lma - d, n_cohorts, frozen_state, t, L)) / (2.0 * d);
  (void)value;

  return Rcpp::List::create(
      Rcpp::Named("n_cohorts") = double(n_cohorts),
      Rcpp::Named("value") = value,
      Rcpp::Named("grad_ad") = grad_ad,
      Rcpp::Named("fd") = fd);
}
