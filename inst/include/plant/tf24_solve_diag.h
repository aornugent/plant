#ifndef PLANT_TF24_SOLVE_DIAG_H
#define PLANT_TF24_SOLVE_DIAG_H

// Stage-1 event classifier: per-cohort branch-signature sink for the TF24
// root-collar solve (off by default; zero cost when off -- one bool test).
//
// The Leaf hydraulic solve (Leaf::prepare_collar_solve) exits through one of a
// small set of branches -- shutdown (wettest layer drier than psi_crit /
// E_column<0 / root_crit), the assim<0 zero-transpiration point, the collapsed
// feasible interval, or the full golden-section optimisation. Which branch a
// cohort takes is a discrete event surface in the RHS: when a cohort crosses
// from one branch to another between steps, the derivative kinks. The Leaf is
// reused across cohorts, so the per-cohort branch code cannot be read back after
// the fact -- it is appended here as each cohort solves. compute_species_rates
// resets the sink at the top of every full RHS sweep, so after any evaluation
// the sink holds exactly that sweep's cohorts. Patch::step_monitor aggregates it
// (branch-code histogram + closest shutdown/interval margins) into the odelia
// step-monitor record. Not part of any model or gradient path.

#include <vector>

namespace plant {
namespace tf24_solve_diag {

extern bool enabled;
extern std::vector<int>    branch;    // per-cohort branch code (see below)
extern std::vector<double> shutdown;  // psi_crit + wettest (<=0 => in shutdown)
extern std::vector<double> interval;  // feasible interval width bound_b-bound_a
                                      // (only meaningful for the GSS branches)

// Branch codes for Leaf::prepare_collar_solve exits.
enum Branch {
  BRANCH_SHUTDOWN_WETTEST = 0,  // -wettest >= psi_crit
  BRANCH_SHUTDOWN_ECOLUMN = 1,  // E_column(-psi_crit) < 0
  BRANCH_SHUTDOWN_ROOTCRIT = 2, // -root_crit >= psi_crit
  BRANCH_ZERO_E = 3,            // assim_max_ < 0 (zero-transpiration point)
  BRANCH_COLLAPSED = 4,         // |bound_b - bound_a| <= GSS_tol
  BRANCH_GSS = 5,               // full golden-section optimisation
  N_BRANCH = 6
};

inline void reset() {
  branch.clear();
  shutdown.clear();
  interval.clear();
}

inline void record(int code, double shutdown_margin, double interval_width) {
  if (enabled) {
    branch.push_back(code);
    shutdown.push_back(shutdown_margin);
    interval.push_back(interval_width);
  }
}

}  // namespace tf24_solve_diag
}  // namespace plant

#endif  // PLANT_TF24_SOLVE_DIAG_H
