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
  enum class Metric { offspring_production, net_reproduction_ratio, LAI, biomass, basal_area };

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
      case Metric::net_reproduction_ratio:
        out.push_back(net_reproduction_ratio<S>(patch));
        break;
      case Metric::LAI:
      case Metric::biomass:
      case Metric::basal_area:
        out.push_back(census<S>(patch, m));
        break;
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
    if (name == "net_reproduction_ratio" || name == "R0") {
      return Metric::net_reproduction_ratio;
    }
    if (name == "LAI") return Metric::LAI;
    if (name == "biomass") return Metric::biomass;
    if (name == "basal_area") return Metric::basal_area;
    util::stop("unknown emergent metric: " + name);
    return Metric::offspring_production; // unreachable
  }

  // Stand census metric: the size-distribution integral summed over species, i.e.
  // the descending-height trapezium of density * per-plant psi over the active
  // cohorts (Species::census). Each metric selects a per-plant kernel reusing the
  // strategy's own reductions -- LAI = k_I * leaf area (compute_competition),
  // biomass = live + heartwood mass, basal_area = stem area.
  template <class S, class Patch>
  static S census(const Patch& patch, Metric m) {
    S total(0.0);
    for (std::size_t s = 0; s < patch.size(); ++s) {
      const auto& sp = patch.at_species(s);
      switch (m) {
      case Metric::LAI:
        total += sp.census([](const auto& ind) { return ind.compute_competition(0.0); });
        break;
      case Metric::biomass:
        total += sp.census([](const auto& ind) { return ind.census_biomass(); });
        break;
      case Metric::basal_area:
        total += sp.census([](const auto& ind) { return ind.census_basal_area(); });
        break;
      default:
        break;
      }
    }
    return total;
  }

  // Total survival-weighted offspring across the stand: the trapezium of each
  // species' per-node weighted fecundity, scaled by the birth rate, over the
  // frozen introduction times (constant weights). The birth rate carries the
  // active scale, so the birth-rate gradient sees both this direct scaling (the
  // identity metric/birth_rate) and the density feedback through the recomputed
  // canopy. At S = double / an unseeded scale the value is bit-identical.
  template <class S, class Patch>
  static S offspring_production(const Patch& patch) {
    S total(0.0);
    for (std::size_t s = 0; s < patch.size(); ++s) {
      const auto& species = patch.at_species(s);
      const std::vector<double> times = species.node_times();
      if (times.size() < 2) continue;
      const auto weighted = species.net_reproduction_ratio_by_node_weighted();
      const auto drivers = species.extrinsic_drivers();
      const S scale = species.birth_rate_scale;
      for (std::size_t i = 1; i < times.size(); ++i) {
        const double t0 = times[i - 1], t1 = times[i];
        const S b0 = scale * drivers.evaluate("birth_rate", t0);
        const S b1 = scale * drivers.evaluate("birth_rate", t1);
        total += 0.5 * (t1 - t0) * (weighted[i] * b1 + weighted[i - 1] * b0);
      }
    }
    return total;
  }

  // Net reproduction ratio R0 (expected offspring per seed), summed over species:
  // the trapezium of per-node weighted fecundity over the introduction times, with
  // NO birth-rate factor (R0 is per seed). It depends on the birth rate only
  // through the density-feedback canopy on the resident replay, so d R0/d birth_rate
  // is purely that feedback -- the plant-side derivative for the R0 = 1 Newton solve.
  // Mirrors Patch::net_reproduction_ratio_for_species with unit scalars.
  template <class S, class Patch>
  static S net_reproduction_ratio(const Patch& patch) {
    S total(0.0);
    for (std::size_t s = 0; s < patch.size(); ++s) {
      const auto& species = patch.at_species(s);
      const std::vector<double> times = species.node_times();
      if (times.size() < 2) continue;
      const auto weighted = species.net_reproduction_ratio_by_node_weighted();
      for (std::size_t i = 1; i < times.size(); ++i) {
        total += 0.5 * (times[i] - times[i - 1]) * (weighted[i] + weighted[i - 1]);
      }
    }
    return total;
  }
};

}

#endif
