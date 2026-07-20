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

struct scm_gradient_result {
  std::vector<double>              value;         // active values (m outputs)
  std::vector<std::vector<double>> jacobian;      // m x n
  std::vector<double>              value_double;  // double reference (m)
};

// The Jacobian (m outputs x n targets) of `functional` over the SCM solve, seeding
// the named parameter targets. Owns record->replay: an adaptive double run discovers
// the resolved schedule, the active run (SCM::rebind_from) replays exactly it. Returns
// {values, jacobian, value_double}; value vs value_double is the config-crossing
// check. `targets.params` indexes field_ptrs()/field_names(); ics are unused here.
template <class StratD, class EnvD, class Functional>
scm_gradient_result
scm_jacobian(Parameters<StratD, EnvD> p, Control control,
             const odelia::ode::DifferentiationTargets& targets,
             Functional functional) {
  using RevS = xad::adj<double>::active_type;

  // 1. Adaptive double pass: discover and record the resolved L0+L1 schedule into
  //    the SCM's Parameters.
  EnvD env;
  SCM<StratD, EnvD> scm(p, env, control);
  scm.refine_schedule();

  // 2. Double reference: rebind to double (an identity copy through the same
  //    contract path, so construction matches the active pass exactly), replay the
  //    resolved schedule, reduce through the functional.
  scm_gradient_result out;
  {
    auto ref = scm.template rebind_from<double>();
    use_recorded_ode_times(ref);
    ref.run();
    auto vs = functional(ref);
    out.value_double.reserve(vs.size());
    for (auto const& v : vs) out.value_double.push_back(xad::value(v));
  }

  // 3. Active pass: rebind config to the reverse scalar, replay the SAME resolved
  //    schedule, one reverse sweep for the whole Jacobian.
  {
    auto active = scm.template rebind_from<RevS>();
    use_recorded_ode_times(active);
    auto [values, jacobian] =
        odelia::ode::compute_jacobian(active, targets, functional);
    out.value    = std::move(values);
    out.jacobian = std::move(jacobian);
  }

  // R5 as structure, not convention: the active value must reproduce the double
  // reference. It does (bit-identically) when every configuration member crosses
  // double->active; a member dropped by rebind_from shifts the value by O(1). So a
  // silently config-incomplete gradient becomes a loud failure here rather than a
  // plausible wrong number the caller has to catch. (TF24's env soil config does not
  // cross yet -- this is what would trip if a TF24 gradient were attempted; b1.)
  for (std::size_t i = 0; i < out.value.size(); ++i) {
    const double a = out.value[i], d = out.value_double[i];
    if (std::abs(a - d) > 1e-8 * (std::abs(d) + 1e-8)) {
      util::stop("scm_gradient: active value does not reproduce the double "
                 "reference -- a configuration member was not carried across "
                 "double->active (rebind_from is incomplete for this System)");
    }
  }
  return out;
}

// Scalar convenience: the gradient (one row) of a scalar functional. Returns
// {value, gradient, value_double}.
template <class StratD, class EnvD, class Functional>
std::tuple<double, std::vector<double>, double>
scm_gradient(Parameters<StratD, EnvD> p, Control control,
             const odelia::ode::DifferentiationTargets& targets,
             Functional functional) {
  auto r = scm_jacobian(std::move(p), std::move(control), targets,
                        std::move(functional));
  return {r.value[0], r.jacobian[0], r.value_double[0]};
}

} // namespace plant

#endif
