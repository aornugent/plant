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
#include <cstdio>
#include <cstdlib>
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

// FF16 standing-stock census vector -- LAI, above-ground biomass, and stem basal
// area, each per patch area (codomain 3). One SCM recording feeds three reverse
// sweeps via scm_jacobian. Every metric is Patch::census of a per-individual
// quantity: the mass-weighted population reduction Sum_i n_i * psi(state_i), so a
// trait's effect flows through the whole re-shaded stand, not one cohort.
struct census_vector {
  std::size_t codomain() const { return 3; }
  template <class Runner>
  std::vector<typename Runner::value_type> operator()(Runner& s) const {
    auto& patch = s.get_system_ref();
    return {
      patch.census([](auto const& ind) { return ind.census_leaf_area(); }),
      patch.census([](auto const& ind) { return ind.census_mass(); }),
      patch.census([](auto const& ind) { return ind.census_basal_area(); }),
    };
  }
};

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

  // Adaptive double pass: discover and record the resolved L1 schedule. Its
  // recorded_steps() is the single source of the replay grid for both passes below.
  EnvD env;
  SCM<StratD, EnvD> scm(p, env, control);
  scm.refine_schedule();
  const std::vector<double> schedule = scm.recorded_steps();

  // Double reference value, replayed on the recorded schedule through the same
  // rebind_from path the active pass uses -- the R5 check below.
  std::vector<double> value_double;
  {
    auto ref = scm.template rebind_from<double>();
    ref.set_schedule(schedule);
    ref.run();
    for (auto const& v : functional(ref)) value_double.push_back(xad::value(v));
  }

  // Active pass: rebind config to the reverse scalar, replay the SAME recorded
  // schedule, one reverse sweep for the whole Jacobian.
  auto active = scm.template rebind_from<RevS>();
  active.set_schedule(schedule);
  auto out = odelia::ode::compute_jacobian(active, targets, functional);

  // Diagnostic (PLANT_TAPE_STATS): the reverse tape is still populated here (the
  // driver's guard deactivates but does not free it). Report its peak memory and
  // operation count alongside the run shape -- steps (L1 grid) and final node ODE
  // width -- so the OOM growth curve is a measurement, not a guess. Env-gated,
  // stderr, no effect on the returned value.
  if (std::getenv("PLANT_TAPE_STATS") && active.tape) {
    std::fprintf(stderr,
                 "TAPE_STATS steps=%zu node_ode_final=%d mem_bytes=%zu ops=%zu stmts=%zu "
                 "deriv=%zu chkpt=%zu maxderiv=%zu\n",
                 schedule.size(),
                 active.get_system_ref().node_ode_size(),
                 active.tape->getMemory(),
                 static_cast<std::size_t>(active.tape->getNumOperations()),
                 static_cast<std::size_t>(active.tape->getNumStatements()),
                 static_cast<std::size_t>(active.tape->diagNumDerivatives()),
                 static_cast<std::size_t>(active.tape->diagNumCheckpoints()),
                 static_cast<std::size_t>(active.tape->diagMaxDerivative()));
    std::fflush(stderr);
  }

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
