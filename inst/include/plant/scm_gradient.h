// -*-c++-*-
#ifndef PLANT_PLANT_SCM_GRADIENT_H_
#define PLANT_PLANT_SCM_GRADIENT_H_

// Run-shaped reverse-mode trait gradient of an SCM reduction.
//
// One entry point, scm_gradient/scm_jacobian, replaces the per-gradient bespoke
// driver. The caller hands in a double Parameters, a set of target trait indices,
// and a functional (the scalar reduction to differentiate); it gets back the value
// and its gradient. The caller NEVER hands in a replay schedule: the entry runs the
// adaptive double pass itself (refine_schedule), so the L1 grid it replays is the
// solver-recorded r_ode_times() of that pass -- there is no schedule argument, hence
// no way to replay a wrong one. This is the guarantee the standalone drivers lacked
// (they pinned patch.step_history, a separate mutant-path record, giving a 60x-wrong
// gradient); here it is kept true by the signature, not by convention.
//
// The double->active crossing is SCM::rebind_from<S2>() (the odelia System contract):
// it carries the same configuration the double run used -- differentiable parameters,
// scalar-independent config, and the resolved schedule -- and leaves precomputed
// per-strategy state for prepare_strategy() to rebuild after seeding (the reset-timing
// contract in odelia/AUTODIFF.md). A dropped config member would make the active value
// disagree with the double reference, which the entry returns alongside for that check.

#include <plant/scm.h>
#include <odelia/gradient.hpp>
#include <cmath>
#include <cstddef>
#include <tuple>
#include <utility>
#include <vector>

namespace plant {

// Reductions of the positioned patch to the output(s) to differentiate. Each reads
// the active patch directly (get_system_ref()), never the SCM's double R facades,
// and returns a vector so it doubles as a compute_jacobian functional (codomain()
// outputs). A scalar reduction reports codomain() == 1 and returns a one-element
// vector.
struct offspring_metric {
  std::size_t codomain() const { return 1; }
  template <class Runner>
  std::vector<typename Runner::value_type> operator()(Runner& s) const {
    return { s.get_system_ref().offspring_production()[0] };
  }
};

// Total competition at ground level -- the density-weighted census reduction.
struct census_metric {
  std::size_t codomain() const { return 1; }
  template <class Runner>
  std::vector<typename Runner::value_type> operator()(Runner& s) const {
    return { s.get_system_ref().compute_competition(0.0) };
  }
};

// Pin an SCM to its own recorded ode times (the resident L1 replay). The schedule is
// already loaded from Parameters (make_node_schedule); this flips the NodeSchedule
// onto it after a reset.
template <class Runner>
static void use_recorded_ode_times(Runner& scm) {
  scm.reset();
  NodeSchedule ns = scm.r_node_schedule();
  ns.r_set_use_ode_times(true);
  scm.r_set_node_schedule(ns);
}

// The Jacobian (m outputs x n targets) of `functional` over the SCM solve, seeding
// the named parameter targets. The plant analogue of odelia's jacobian_on_double
// (plant's SCM is not an ode::Solver -- it HAS-A one and adds node scheduling -- so
// it cannot use that entry): it owns the record->replay odelia's Solver owns for a
// bare ODE. An adaptive double run discovers the resolved schedule; the active twin
// (SCM::rebind_from, the odelia System hook) replays exactly it; odelia's
// compute_jacobian does the seed/tape/sweep. Returns odelia's {values, jacobian}.
// `targets.params` indexes field_ptrs()/field_names(); ics are unused here.
template <class StratD, class EnvD, class Functional>
std::pair<std::vector<double>, std::vector<std::vector<double>>>
scm_jacobian(Parameters<StratD, EnvD> p, Control control,
             const odelia::ode::DifferentiationTargets& targets,
             Functional functional) {
  using RevS = xad::adj<double>::active_type;

  // Adaptive double pass: discover and record the resolved L0+L1 schedule.
  EnvD env;
  SCM<StratD, EnvD> scm(p, env, control);
  scm.refine_schedule();

  // Double reference value, replayed on the resolved schedule through the same
  // rebind_from path the active pass uses -- the R5 check below.
  std::vector<double> value_double;
  {
    auto ref = scm.template rebind_from<double>();
    use_recorded_ode_times(ref);
    ref.run();
    for (auto const& v : functional(ref)) value_double.push_back(xad::value(v));
  }

  // Active pass: rebind config to the reverse scalar, replay the SAME schedule, one
  // reverse sweep for the whole Jacobian.
  auto active = scm.template rebind_from<RevS>();
  use_recorded_ode_times(active);
  auto out = odelia::ode::compute_jacobian(active, targets, functional);

  // R5 as structure, not convention: the active value must reproduce the double
  // reference. It does (bit-identically) when every configuration member crosses
  // double->active; a member dropped by rebind_from shifts the value by O(1), so a
  // config-incomplete gradient fails loudly here instead of returning a plausible
  // wrong number. (TF24's env soil config does not cross yet -- this trips if a TF24
  // gradient is attempted; b1.)
  for (std::size_t i = 0; i < out.first.size(); ++i) {
    if (std::abs(out.first[i] - value_double[i]) >
        1e-8 * (std::abs(value_double[i]) + 1e-8)) {
      util::stop("scm_jacobian: active value does not reproduce the double "
                 "reference -- a configuration member was not carried across "
                 "double->active (rebind_from is incomplete for this System)");
    }
  }
  return out;
}

// Scalar convenience: the gradient (one row) of a scalar functional, as odelia's
// {value, gradient}.
template <class StratD, class EnvD, class Functional>
std::pair<double, std::vector<double>>
scm_gradient(Parameters<StratD, EnvD> p, Control control,
             const odelia::ode::DifferentiationTargets& targets,
             Functional functional) {
  auto [values, jacobian] = scm_jacobian(std::move(p), std::move(control), targets,
                                         std::move(functional));
  return {values[0], jacobian[0]};
}

} // namespace plant

#endif
