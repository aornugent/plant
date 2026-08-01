// Census of how the collar-potential polish terminates: at an interior stationary
// point of profit, or against a bracket end. Every entry point is behind
// PLANT_COLLAR_CENSUS and the class is never constructed without it.
//
// The class is whichever guard in the polish broke its loop. A distance from the
// returned potential to the bracket end would instead be tolerance-dependent,
// reporting nearly every solve as pinned at a loose GSS_tol_abs.
//
// Accumulates counts and extremes rather than rows: one record per leaf solve, of
// order 4 million on a production lifetime.
#ifndef PLANT_PLANT_COLLAR_CENSUS_H_
#define PLANT_PLANT_COLLAR_CENSUS_H_

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <string>

namespace plant {
namespace internals {

inline bool collar_census_active() {
  static const bool active = getenv("PLANT_COLLAR_CENSUS") != nullptr;
  return active;
}

// The polish's termination branches. INTERIOR and EXHAUSTED are not pinned; the
// four BOUND_ classes are, and name which guard stopped Newton.
enum collar_class {
  COLLAR_INTERIOR = 0,     // |R| <= R_tol
  COLLAR_R_NONFINITE,      // R not finite, so no step is defined
  COLLAR_BOUND_A,          // psi - h <= bound_a: no room for the probe
  COLLAR_BOUND_B,          // psi + h >= bound_b: no room for the probe
  COLLAR_BOUND_CURVATURE,  // dR/d(collar) at psi is not negative
  COLLAR_BOUND_STEP,       // the Newton step from psi lands outside the bracket
  COLLAR_EXHAUSTED,        // five steps taken, none of them closing
  COLLAR_N_CLASS
};

// Running min/max/mean of one quantity.
struct census_span {
  double lo = 0.0, hi = 0.0, sum = 0.0;
  long long n = 0;

  void offer(double v) {
    if (n == 0) { lo = v; hi = v; }
    else { if (v < lo) lo = v; if (v > hi) hi = v; }
    sum += v;
    n++;
  }
};

class collar_census {
public:
  // psi_wet is the wettest layer's potential as a positive magnitude, radiation
  // the value the strategy passed as PPFD, area_leaf the solving leaf's area.
  void add(int cls, double abs_R, double psi_wet, double radiation,
           double area_leaf) {
    n++;
    if (cls < 0 || cls >= COLLAR_N_CLASS) return;
    n_class[cls]++;
    if (std::isfinite(abs_R)) R_span[cls].offer(abs_R);
    else n_R_nonfinite[cls]++;
    if (cls >= COLLAR_BOUND_A && cls <= COLLAR_BOUND_STEP) {
      pinned_psi.offer(psi_wet);
      pinned_radiation.offer(radiation);
      pinned_area_leaf.offer(area_leaf);
    } else {
      interior_psi.offer(psi_wet);
      interior_radiation.offer(radiation);
      interior_area_leaf.offer(area_leaf);
    }
  }

  void add_early_exit() { n_early_exit++; }

  void add_bracket(double width) {
    if (n_bracket == 0 || width < min_bracket) min_bracket = width;
    n_bracket++;
  }

  ~collar_census() {
    if (n == 0 && n_early_exit == 0) return;
    const char* path = getenv("PLANT_COLLAR_CENSUS");
    FILE* f = (path && path[0] && std::string(path) != "1") ? fopen(path, "w")
                                                           : stderr;
    fprintf(f, "# collar census: polish termination class per leaf solve\n");
    fprintf(f, "polished_solves %lld\n", n);
    fprintf(f, "early_exit_solves %lld\n", n_early_exit);
    fprintf(f, "min_bracket %.17g\n", n_bracket ? min_bracket : 0.0);
    long long pinned = 0;
    for (int c = COLLAR_BOUND_A; c <= COLLAR_BOUND_STEP; ++c) pinned += n_class[c];
    fprintf(f, "pinned_total %lld\n", pinned);
    for (int c = 0; c < COLLAR_N_CLASS; ++c) {
      fprintf(f, "class %s count %lld R_nonfinite %lld", name(c), n_class[c],
              n_R_nonfinite[c]);
      const census_span& s = R_span[c];
      if (s.n == 0) fprintf(f, " absR none\n");
      else fprintf(f, " absR_min %.17g absR_max %.17g absR_mean %.17g\n", s.lo,
                   s.hi, s.sum / static_cast<double>(s.n));
    }
    report(f, "pinned", pinned_psi, pinned_radiation, pinned_area_leaf);
    report(f, "not_pinned", interior_psi, interior_radiation, interior_area_leaf);
    if (f != stderr) fclose(f);
  }

private:
  static const char* name(int c) {
    switch (c) {
    case COLLAR_INTERIOR: return "interior";
    case COLLAR_R_NONFINITE: return "R_nonfinite";
    case COLLAR_BOUND_A: return "pinned_bound_a";
    case COLLAR_BOUND_B: return "pinned_bound_b";
    case COLLAR_BOUND_CURVATURE: return "pinned_curvature";
    case COLLAR_BOUND_STEP: return "pinned_step_outside";
    case COLLAR_EXHAUSTED: return "exhausted";
    default: return "unknown";
    }
  }

  static void report(FILE* f, const char* tag, const census_span& psi,
                     const census_span& rad, const census_span& area) {
    if (psi.n == 0) { fprintf(f, "%s_states none\n", tag); return; }
    const double d = static_cast<double>(psi.n);
    fprintf(f,
            "%s_states n %lld psi_wet %.17g .. %.17g mean %.17g"
            " radiation %.17g .. %.17g mean %.17g"
            " area_leaf %.17g .. %.17g mean %.17g\n",
            tag, psi.n, psi.lo, psi.hi, psi.sum / d, rad.lo, rad.hi, rad.sum / d,
            area.lo, area.hi, area.sum / d);
  }

  long long n = 0, n_early_exit = 0, n_bracket = 0;
  long long n_class[COLLAR_N_CLASS] = {0};
  long long n_R_nonfinite[COLLAR_N_CLASS] = {0};
  double min_bracket = 0.0;
  census_span R_span[COLLAR_N_CLASS];
  census_span pinned_psi, pinned_radiation, pinned_area_leaf;
  census_span interior_psi, interior_radiation, interior_area_leaf;
};

inline collar_census& the_collar_census() {
  static collar_census census;
  return census;
}

}  // namespace internals
}  // namespace plant

#endif
