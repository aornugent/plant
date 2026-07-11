// -*-c++-*-
#ifndef PLANT_PLANT_PARAMETERS_H_
#define PLANT_PLANT_PARAMETERS_H_

#include <vector>
#include <RcppCommon.h> // SEXP

#include <plant/control.h>
#include <plant/models/ff16_strategy.h>
#include <plant/node_schedule.h>
#include <plant/scm_utils.h> // Unfortunately needed for setup_node_schedule

#include <plant/disturbance_regime.h>
#include <plant/disturbances/no_disturbance.h>
#include <plant/disturbances/weibull_disturbance.h>

// TODO(#483): I will possibly move out the "Patch" parameters out into
// their own simple list class at some point, to make this a bit more
// coherent.
//
// TODO(#483): Will require some free functions on the R side:
//   * add_strategy (with flag for mutant/non mutant)

namespace plant {

template <typename T, typename E>
struct Parameters {
  typedef T strategy_type;
  typedef E environment_type;

  Parameters() :
    patch_area(1.0),
    n_patches(1),
    patch_type("meta-population"),
    max_patch_lifetime(105.32), // designed to agree with Daniel's implementation
    initial_time(0.0)
  {
    validate();
  }

  // Data -- public for now (see github issue #17).
  double patch_area; // Size of the patch (m^2)
  size_t n_patches;  // Number of patches in the metacommunity
  std::string patch_type;
  double max_patch_lifetime; // Disturbance interval (years)
  std::vector<strategy_type> strategies;

  Disturbance_Regime* disturbance;

  // Default strategy.
  strategy_type strategy_default;

  // Node information.
  std::vector<double> node_schedule_times_default;
  std::vector<std::vector<double> > node_schedule_times;
  std::vector<double> ode_times;

  // Initial patch state. When initial_state is non-empty the patch is seeded
  // with these nodes at reset() instead of starting empty -- used to resume an
  // exported patch run or to seed an arbitrary initial size distribution.
  // Carried here so a run stays fully self-describing/serialisable (see
  // agents.md), and so the seeding survives the reset() at the start of every
  // SCM::run()/refine_schedule().
  //   initial_state             flat ODE state (all nodes, then environment),
  //                             in the order Patch::set_ode_state expects.
  //   n_initial_cohorts         number of nodes per species.
  //   initial_node_times        per-node introduction time (flat across species).
  //   initial_patch_density     per-node patch-age density at birth (flat).
  //   initial_pr_patch_survival per-node pr_patch_survival at birth (flat).
  //   initial_time              patch age to resume at (0 for a fresh seed).
  // Length consistency is checked in Patch::set_initial_state(), where the
  // node/ode sizes are known.
  std::vector<double> initial_state;
  std::vector<size_t> n_initial_cohorts;
  std::vector<double> initial_node_times;
  std::vector<double> initial_patch_density;
  std::vector<double> initial_pr_patch_survival;
  double initial_time;

  // Config-only copy onto another scalar S2 (the strategies carry their trait
  // values across via ad_value; the schedule/patch fields are plain doubles).
  // Used by Patch::rebind_from to lift the whole configuration to the active
  // scalar. The returned Parameters is validated (disturbance + node schedule).
  // Deduced return type: the body (which needs T::rebind) is instantiated only
  // when used, so strategies without a rebind still form a valid Parameters.
  template <class S2>
  auto rebind_from() const {
    Parameters<typename T::template rebind<S2>, E> out;
    out.patch_area = patch_area;
    out.n_patches = n_patches;
    out.patch_type = patch_type;
    out.max_patch_lifetime = max_patch_lifetime;
    for (const auto& s : strategies) {
      out.strategies.push_back(s.template rebind_from<S2>());
    }
    out.strategy_default = strategy_default.template rebind_from<S2>();
    out.node_schedule_times_default = node_schedule_times_default;
    out.node_schedule_times = node_schedule_times;
    out.ode_times = ode_times;
    out.initial_state = initial_state;
    out.n_initial_cohorts = n_initial_cohorts;
    out.initial_node_times = initial_node_times;
    out.initial_patch_density = initial_patch_density;
    out.initial_pr_patch_survival = initial_pr_patch_survival;
    out.initial_time = initial_time;
    out.validate();
    return out;
  }

  // Some little query functions for use on the C side:
  size_t size() const;
  void validate();

private:
  void setup_node_schedule();
};

template <typename T, typename E>
size_t Parameters<T,E>::size() const {
  return strategies.size();
}

// NOTE: this will be called *every time* that the object is passed in
// from R -> C++.  That's unlikely to be that often, but it does incur
// a penalty.  So don't put anything too stupidly heavy in here.
template <typename T, typename E>
void Parameters<T,E>::validate() {
  const size_t n_spp = size();

  setup_node_schedule();
  if (node_schedule_times.size() != n_spp) {
    util::stop("Incorrect length node_schedule_times");
  }

  // Disturbances used to describe evolution of a metapopulation of patches
  // when calculating fitness, otherwise defaults to fixed-duration run without
  // disturbance
  if(patch_type == "meta-population") {
    disturbance = new Weibull_Disturbance_Regime(max_patch_lifetime);
  }
  else {
    disturbance = new No_Disturbance();
  }
}

// Separating this out just because it's a bit crap:
// TODO(#483): Consider adding this to scm_utils.h perhaps?
template <typename T, typename E>
void Parameters<T,E>::setup_node_schedule() {
  node_schedule_times_default =
      plant::node_schedule_times_default(max_patch_lifetime);

  if ((node_schedule_times.empty() && size() > 0)) {
    node_schedule_times.clear();
    for (size_t i = 0; i < size(); ++i) {
      node_schedule_times.push_back(node_schedule_times_default);
    }
  }
}
}

#endif
