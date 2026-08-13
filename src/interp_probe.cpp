// Scaffolding for the light-field placement study: selects the knot placement
// at run time so one build can time and referee every candidate, and reports
// the grid a build actually laid out.

#include <plant.h>
#include <plant/resource_spline.h>
#include <Rcpp.h>
#include <string>

// [[Rcpp::export]]
void interp_policy_set(std::string mode, double delta, int count, int pad,
                       double ceiling = 0.0) {
  if (mode == "canopy") {
    plant::ResourceGridPolicy::mode =
      plant::ResourceGridPolicy::Mode::CanopyUniform;
  } else if (mode == "fixed") {
    plant::ResourceGridPolicy::mode =
      plant::ResourceGridPolicy::Mode::FixedAbsolute;
  } else {
    plant::util::stop("interp_policy_set: mode must be 'canopy' or 'fixed'");
  }
  if (!(delta > 0.0)) plant::util::stop("interp_policy_set: delta must be > 0");
  if (count < 2)      plant::util::stop("interp_policy_set: count must be >= 2");
  if (pad < 0)        plant::util::stop("interp_policy_set: pad must be >= 0");
  plant::ResourceGridPolicy::delta = delta;
  plant::ResourceGridPolicy::count = static_cast<size_t>(count);
  plant::ResourceGridPolicy::pad   = static_cast<size_t>(pad);
  plant::ResourceGridPolicy::ceiling = ceiling;
}

// [[Rcpp::export]]
Rcpp::List interp_policy_get() {
  const bool fixed = plant::ResourceGridPolicy::mode ==
                     plant::ResourceGridPolicy::Mode::FixedAbsolute;
  return Rcpp::List::create(
    Rcpp::_["mode"]  = std::string(fixed ? "fixed" : "canopy"),
    Rcpp::_["delta"] = plant::ResourceGridPolicy::delta,
    Rcpp::_["count"] = static_cast<int>(plant::ResourceGridPolicy::count),
    Rcpp::_["pad"]   = static_cast<int>(plant::ResourceGridPolicy::pad),
    Rcpp::_["ceiling"] = plant::ResourceGridPolicy::ceiling);
}
