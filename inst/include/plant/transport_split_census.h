// Splits the cohort-grid transport stencil by the pair it was taken over: the
// boundary pair (lowest node against the inflow node) against interior pairs.
// Instrumentation only; every entry point is behind PLANT_SPLIT_CENSUS and
// nothing is accumulated when the variable is unset.
#ifndef PLANT_PLANT_TRANSPORT_SPLIT_CENSUS_H_
#define PLANT_PLANT_TRANSPORT_SPLIT_CENSUS_H_

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>

namespace plant {
namespace internals {

inline bool split_census_active() {
  static const bool active = getenv("PLANT_SPLIT_CENSUS") != nullptr;
  return active;
}

struct split_record {
  double value = 0.0, time = 0.0, height = 0.0, dh = 0.0, g = 0.0, g_below = 0.0;
};

// One arm of the split: count, mean/max of |stencil|, and the largest few
// records kept in full so the extremes can be traced back to a state.
struct split_arm {
  long long n = 0, n_guarded = 0;
  double sum_abs = 0.0, max_abs = 0.0, min_dh = 1.0 / 0.0;
  std::vector<split_record> top;

  void add(double t, double h, double dh, double gi, double gb, double v) {
    n++;
    const double a = std::fabs(v);
    sum_abs += a;
    if (a > max_abs) max_abs = a;
    if (dh < min_dh) min_dh = dh;
    if (!(dh > 0.0)) n_guarded++;
    if (top.size() < 8 || a > top.back().value) {
      top.push_back(split_record{a, t, h, dh, gi, gb});
      std::sort(top.begin(), top.end(),
                [](const split_record& x, const split_record& y) {
                  return x.value > y.value;
                });
      if (top.size() > 8) top.resize(8);
    }
  }

  void report(FILE* f, const char* name) const {
    fprintf(f, "%s n %lld guarded %lld mean_abs %.17g max_abs %.17g min_dh %.17g\n",
            name, n, n_guarded, n ? sum_abs / static_cast<double>(n) : 0.0,
            max_abs, min_dh);
    for (const auto& r : top) {
      fprintf(f, "%s_top value %.17g time %.17g height %.17g dh %.17g g %.17g"
                 " g_below %.17g\n",
              name, r.value, r.time, r.height, r.dh, r.g, r.g_below);
    }
  }
};

class split_census {
public:
  void add(bool boundary, double t, double h, double dh, double gi, double gb,
           double v) {
    (boundary ? bnd : intr).add(t, h, dh, gi, gb, v);
  }

  ~split_census() {
    if (bnd.n == 0 && intr.n == 0) return;
    const char* path = getenv("PLANT_SPLIT_CENSUS");
    FILE* f = (path && path[0] && std::string(path) != "1") ? fopen(path, "w")
                                                            : stderr;
    fprintf(f, "# transport stencil split by pair\n");
    bnd.report(f, "boundary");
    intr.report(f, "interior");
    if (f != stderr) fclose(f);
  }

private:
  split_arm bnd, intr;
};

inline split_census& the_split_census() {
  static split_census census;
  return census;
}

}  // namespace internals
}  // namespace plant

#endif
