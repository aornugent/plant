// Storage + R interface for the Stage-1 event classifier per-cohort
// branch-signature sink (plant/tf24_solve_diag.h). Enabled alongside odelia's
// step monitor; read by Patch::step_monitor each accepted step.
#include <Rcpp.h>
#include <plant/tf24_solve_diag.h>

namespace plant {
namespace tf24_solve_diag {
bool enabled = false;
std::vector<int>    branch;
std::vector<double> shutdown;
std::vector<double> interval;
}  // namespace tf24_solve_diag
}  // namespace plant

// [[Rcpp::export]]
void tf24_solve_diag_enable(bool on) {
  plant::tf24_solve_diag::enabled = on;
}
