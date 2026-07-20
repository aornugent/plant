#ifndef PLANT_EXTRINSIC_DRIVERS_H
#define PLANT_EXTRINSIC_DRIVERS_H

#include <plant/util.h>
#include <odelia/interpolator.hpp>
#include <Rcpp.h>
#include <algorithm>
#include <limits>
#include <vector>

namespace plant {
namespace {

class Function {
public:
  Function() = default;

  Function(std::vector<double> const &x, std::vector<double> const &y) {
    variable.init(x, y);
    variable.set_extrapolate(false);
    is_variable = true;
    // Feature times for the forcing-kink clip: knots where the driver value
    // changes (onset/cessation of an event). A run of equal values (e.g. a dry
    // spell of zero rain) contributes none, so the clip leaves smooth stretches
    // alone and only lands steps on genuine forcing features -- not on every
    // daily knot, which would force sub-daily steps everywhere.
    for (size_t i = 1; i < x.size(); ++i) {
      if (y[i] != y[i - 1]) {
        feature_x_.push_back(x[i]);
      }
    }
  }

  Function(double k) {
    constant = k;
    is_variable = false;
  }

  double evaluate(double u) const {
    if (is_variable) {
      return variable.eval(u);
    } else {
      return constant;
    }
  }

  std::vector<double> evaluate_range(const std::vector<double> &u) const {
    if (is_variable) {
      return variable.r_eval(u);
    } else {
      return std::vector<double>(u.size(), constant);
    }
  }

  void set_extrapolate(bool extrapolate) {
    variable.set_extrapolate(extrapolate);
  }

  // Next forcing feature time strictly after t (+inf for a constant driver or
  // past the last feature). Used by the forcing-kink clip. Features are the
  // value-change knots precomputed at construction, not every spline knot.
  double next_node_after(double t) const {
    if (!is_variable) return std::numeric_limits<double>::infinity();
    auto it = std::upper_bound(feature_x_.begin(), feature_x_.end(), t);
    return it != feature_x_.end() ? *it
                                  : std::numeric_limits<double>::infinity();
  }

private:
    odelia::interpolator::Interpolator variable;
    std::vector<double> feature_x_;
    double constant;
    bool is_variable;
};
}

class ExtrinsicDrivers {
  
public:
  // this will override any previously defined drivers with the same name
  void set_constant(const std::string &driver_name, double k) {

    if (drivers.find(driver_name) != drivers.end())
    {
      drivers.erase(driver_name);
    }
    drivers.insert({driver_name, Function(k)});
  }

  // initialise spline of driver with x, y control points
  void set_variable(const std::string &driver_name, std::vector<double> const &x, std::vector<double> const &y) {
    if (drivers.find(driver_name) != drivers.end())
    {
      drivers.erase(driver_name);
    }
    drivers.insert({driver_name, Function(x, y)});
  }

  void set_extrapolate(const std::string &driver_name, bool extrapolate) {
    drivers.at(driver_name).set_extrapolate(extrapolate);
  }

  // evaluate/query interpolated spline for driver at point u, return s(x), where s is interpolated function
  // (taken by const ref to avoid a per-call std::string copy on the hot path)
  double evaluate(const std::string &driver_name, double x) const {

    return drivers.at(driver_name).evaluate(x);
  }


  // evaluate/query interpolated spline for driver at vector of points, return vector of values
  std::vector<double> evaluate_range(const std::string &driver_name, std::vector<double> x) const {
    return drivers.at(driver_name).evaluate_range(x);
  }

  // returns the name of each active driver - useful for R output
  std::vector<std::string> get_names() const {
    auto ret = std::vector<std::string>();
    for (auto const &driver: drivers) {
      ret.push_back(driver.first);
    }
    return ret;
  }

  void clear() {
    drivers.clear();
  }

  // Smallest feature time strictly after t across all variable drivers, or +inf
  // if none (all constant, or t past every node). The forcing-kink clip lands
  // trial steps on these times. Cheap: one binary search per variable driver.
  double next_node_after(double t) const {
    double best = std::numeric_limits<double>::infinity();
    for (auto const &driver : drivers) {
      best = std::min(best, driver.second.next_node_after(t));
    }
    return best;
  }

private:
  std::unordered_map <std::string, Function> drivers;
};

}

#endif //PLANT_EXTRINSIC_DRIVERS_H
