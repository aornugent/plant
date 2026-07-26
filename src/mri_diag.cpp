// Diagnostic counter for the multirate (method="mri") fast sub-cycle.
//
// Each Patch::fast_rates call evaluates the soil (fast-block) coupling once,
// which costs O(M) cohort physiology solves -- the fast-block cost driver. This
// counter accumulates those calls across a run so we can read the per-macro-step
// fast-evaluation count (Lever-1 headroom: how many micro evaluations the stiff
// drainage currently forces) and compare schemes. Not part of the model; used
// only by the multirate cost measurements.
#include <Rcpp.h>

namespace plant {
long mri_fast_rate_calls = 0;
// Full patch RHS evaluations (Patch::compute_rates), i.e. one evaluation of the
// whole coupled derivative -- the cost/step proxy shared by every global stepper
// (rkck, rodas). Lets us compare accepted-step economics across methods and
// diagnose whether the coupled system is accuracy- or stability-limited.
long patch_rhs_calls = 0;
// Expensive coupling snapshots under method="mri_uptake": each Patch::refresh_anchor
// re-captures a0 = fast_block_uptake() (the O(M) cohort sum) + the uptake Jacobian.
// This is the scarce resource the uptake arbitrage rations -- comparing it against
// mri_fast_rate_calls (the cheap frozen-coupling residual evals) measures the
// realised cohort-sum reduction. Not part of the model.
long mri_coupling_evals = 0;
}

// [[Rcpp::export]]
double mri_fast_rate_calls_get() {
  return static_cast<double>(plant::mri_fast_rate_calls);
}

// [[Rcpp::export]]
void mri_fast_rate_calls_reset() {
  plant::mri_fast_rate_calls = 0;
}

// [[Rcpp::export]]
double patch_rhs_calls_get() {
  return static_cast<double>(plant::patch_rhs_calls);
}

// [[Rcpp::export]]
void patch_rhs_calls_reset() {
  plant::patch_rhs_calls = 0;
}

// [[Rcpp::export]]
double mri_coupling_evals_get() {
  return static_cast<double>(plant::mri_coupling_evals);
}

// [[Rcpp::export]]
void mri_coupling_evals_reset() {
  plant::mri_coupling_evals = 0;
}

// TEMPORARY (anchor-fusion diagnostic, remove once the skip is settled): the
// largest relative disagreement seen between the anchor published by slow_rates
// and the anchor refresh_anchor would have computed by sweeping at the same
// (x, u). Zero means the skip is exact; nonzero localises the drift.
namespace plant {
double anchor_skip_max_reldiff = 0.0;
long   anchor_skip_checks = 0;
}

// [[Rcpp::export]]
double anchor_skip_max_reldiff_get() {
  return plant::anchor_skip_max_reldiff;
}

// [[Rcpp::export]]
double anchor_skip_checks_get() {
  return static_cast<double>(plant::anchor_skip_checks);
}

// [[Rcpp::export]]
void anchor_skip_diag_reset() {
  plant::anchor_skip_max_reldiff = 0.0;
  plant::anchor_skip_checks = 0;
}

namespace plant { bool anchor_skip_diag = false; }

// [[Rcpp::export]]
void anchor_skip_diag_set(bool on) { plant::anchor_skip_diag = on; }
