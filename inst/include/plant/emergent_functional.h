// -*-c++-*-
#ifndef PLANT_PLANT_EMERGENT_FUNCTIONAL_H_
#define PLANT_PLANT_EMERGENT_FUNCTIONAL_H_

#include <plant/util.h>
#include <cstddef>
#include <string>
#include <vector>

namespace plant {

// A pure reduction of a replayed stand run into emergent metrics, matching
// odelia's functional shape (codomain() + operator()(runnable)). It reads the
// runnable's final system through get_system_ref() and drives nothing -- the
// gradient driver owns the replay. Each requested metric name selects one
// reduction kernel; codomain() is the number of metrics, so several stand
// summaries come back as one Jacobian.
class EmergentFunctional {
public:
  enum class Metric { offspring_production, LAI, biomass, basal_area };

  explicit EmergentFunctional(const std::vector<std::string>& names) {
    metrics.reserve(names.size());
    for (const auto& n : names) {
      metrics.push_back(resolve(n));
    }
  }

  std::size_t codomain() const { return metrics.size(); }

  template <class Runnable>
  std::vector<typename Runnable::value_type> operator()(Runnable& runnable) const {
    using S = typename Runnable::value_type;
    const auto& patch = runnable.get_system_ref();
    std::vector<S> out;
    out.reserve(metrics.size());
    for (Metric m : metrics) {
      switch (m) {
      case Metric::offspring_production:
        out.push_back(offspring_production<S>(patch));
        break;
      case Metric::LAI:
      case Metric::biomass:
      case Metric::basal_area:
        util::stop("census metrics (LAI/biomass/basal_area) require the native "
                   "crown integral and are not available yet");
      }
    }
    return out;
  }

private:
  std::vector<Metric> metrics;

  static Metric resolve(const std::string& name) {
    if (name == "offspring_production" || name == "offspring") {
      return Metric::offspring_production;
    }
    if (name == "LAI") return Metric::LAI;
    if (name == "biomass") return Metric::biomass;
    if (name == "basal_area") return Metric::basal_area;
    util::stop("unknown emergent metric: " + name);
    return Metric::offspring_production; // unreachable
  }

  // Total survival-weighted offspring across the stand: the trapezium of each
  // species' per-node weighted fecundity, scaled by the birth rate, over the
  // frozen introduction times (constant weights). Reuses the model's own
  // node-level reduction so the gradient is structural, not a bit-copy.
  template <class S, class Patch>
  static S offspring_production(const Patch& patch) {
    S total(0.0);
    for (std::size_t s = 0; s < patch.size(); ++s) {
      const auto& species = patch.at_species(s);
      const std::vector<double> times = species.node_times();
      if (times.size() < 2) continue;
      const auto weighted = species.net_reproduction_ratio_by_node_weighted();
      const auto drivers = species.extrinsic_drivers();
      for (std::size_t i = 1; i < times.size(); ++i) {
        const double t0 = times[i - 1], t1 = times[i];
        const double b0 = drivers.evaluate("birth_rate", t0);
        const double b1 = drivers.evaluate("birth_rate", t1);
        total += 0.5 * (t1 - t0) * (weighted[i] * b1 + weighted[i - 1] * b0);
      }
    }
    return total;
  }
};

}

#endif
