// Transport-term census: the cohort-grid stencil against the sub-grid probe, on a
// production run. Instrumentation, not model code -- every entry point is behind
// PLANT_TRANSPORT_CENSUS and the class is never constructed without it.
//
// Two questions, one run. How far log_density_dt moves when the stencil changes
// grids, and how large the cohort-grid value gets where the cohort spacing is
// small enough for the divisor to matter.
//
// Accumulates moments and extremes rather than rows: the census is one record per
// node per stage, of order 4 million on a production lifetime.
#ifndef PLANT_PLANT_TRANSPORT_CENSUS_H_
#define PLANT_PLANT_TRANSPORT_CENSUS_H_

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <string>

namespace plant {
namespace internals {

inline bool transport_census_active() {
  static const bool active = getenv("PLANT_TRANSPORT_CENSUS") != nullptr;
  return active;
}

// One extreme, with enough context to find the state it came from.
struct census_peak {
  double value = 0.0;
  double time = 0.0, height = 0.0, dh = 0.0, g = 0.0, g_below = 0.0;
  double subgrid = 0.0, cohort = 0.0;
  bool seen = false;

  void offer(double magnitude, double t, double h, double d, double gi, double gb,
             double sub, double coh) {
    if (seen && !(magnitude > value)) return;
    seen = true;
    value = magnitude; time = t; height = h; dh = d;
    g = gi; g_below = gb; subgrid = sub; cohort = coh;
  }
};

class transport_census {
public:
  void add(double t, double h, double dh, double gi, double gb, double sub,
           double coh) {
    n++;
    const double d = coh - sub;
    sum_d += d; sum_d2 += d * d;
    if (!std::isfinite(d)) n_nonfinite++;
    peak_diff.offer(std::fabs(d), t, h, dh, gi, gb, sub, coh);
    peak_cohort.offer(std::fabs(coh), t, h, dh, gi, gb, sub, coh);
    peak_subgrid.offer(std::fabs(sub), t, h, dh, gi, gb, sub, coh);
    if (dh == 0.0) n_dh_zero++;
    if (dh < 1e-4) {
      n_dh_tight++;
      peak_cohort_tight.offer(std::fabs(coh), t, h, dh, gi, gb, sub, coh);
      peak_diff_tight.offer(std::fabs(d), t, h, dh, gi, gb, sub, coh);
    }
    if (!(dh > min_dh)) min_dh = dh;
    // The stencil's own magnitude bracketed, so the tail is visible without rows.
    const double a = std::fabs(coh);
    if (a > 1e5) n_coh_gt_1e5++;
    else if (a > 1e3) n_coh_gt_1e3++;
    else if (a > 1e1) n_coh_gt_1e1++;
  }

  ~transport_census() {
    if (n == 0) return;
    const char* path = getenv("PLANT_TRANSPORT_CENSUS");
    FILE* f = (path && path[0] && std::string(path) != "1") ? fopen(path, "w")
                                                           : stderr;
    const double mean = sum_d / static_cast<double>(n);
    const double var = sum_d2 / static_cast<double>(n) - mean * mean;
    fprintf(f, "# transport census: cohort-grid stencil vs sub-grid probe\n");
    fprintf(f, "records %lld\n", n);
    fprintf(f, "nonfinite_diff %lld\n", n_nonfinite);
    fprintf(f, "diff_mean %.17g\ndiff_sd %.17g\n", mean,
            var > 0 ? std::sqrt(var) : 0.0);
    fprintf(f, "dh_min %.17g\ndh_zero %lld\ndh_below_1e-4 %lld\n", min_dh,
            n_dh_zero, n_dh_tight);
    fprintf(f, "cohort_abs_gt_1e1 %lld\ncohort_abs_gt_1e3 %lld\n"
               "cohort_abs_gt_1e5 %lld\n",
            n_coh_gt_1e1, n_coh_gt_1e3, n_coh_gt_1e5);
    report(f, "peak_diff", peak_diff);
    report(f, "peak_diff_tight", peak_diff_tight);
    report(f, "peak_cohort", peak_cohort);
    report(f, "peak_cohort_tight", peak_cohort_tight);
    report(f, "peak_subgrid", peak_subgrid);
    if (f != stderr) fclose(f);
  }

private:
  static void report(FILE* f, const char* name, const census_peak& p) {
    if (!p.seen) { fprintf(f, "%s none\n", name); return; }
    fprintf(f,
            "%s value %.17g time %.17g height %.17g dh %.17g g %.17g g_below %.17g"
            " subgrid %.17g cohort %.17g\n",
            name, p.value, p.time, p.height, p.dh, p.g, p.g_below, p.subgrid,
            p.cohort);
  }

  long long n = 0, n_nonfinite = 0, n_dh_zero = 0, n_dh_tight = 0;
  long long n_coh_gt_1e1 = 0, n_coh_gt_1e3 = 0, n_coh_gt_1e5 = 0;
  double sum_d = 0.0, sum_d2 = 0.0;
  double min_dh = 1.0 / 0.0;
  census_peak peak_diff, peak_diff_tight, peak_cohort, peak_cohort_tight,
      peak_subgrid;
};

inline transport_census& the_transport_census() {
  static transport_census census;
  return census;
}

}  // namespace internals
}  // namespace plant

#endif
