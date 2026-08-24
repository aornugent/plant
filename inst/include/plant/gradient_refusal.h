// -*-c++-*-
#ifndef PLANT_PLANT_GRADIENT_REFUSAL_H_
#define PLANT_PLANT_GRADIENT_REFUSAL_H_

#include <stdexcept>
#include <string>
#include <vector>

namespace plant {

// Why a metric has no gradient, and where that was found. A refusal is the only
// thing a number needs said beside it: a parameter no gradient exists for cannot
// be asked for, so it never reaches an answer, and every column that does reach
// one carries a number the sweep computed.
//
// The location is known in pieces -- the leaf knows the reason, the block loop
// knows which plant, the sweep knows which steps -- so it is filled in by the
// frames that know and left at -1 by the frames that do not.
struct refusal {
  std::string reason;
  int species = -1;
  // 0 for the inflow boundary, 1-based for an introduced node.
  int node = -1;
  int step_first = -1;
  int step_last = -1;

  bool happened() const { return !reason.empty(); }
};

// Thrown where a gradient row does not exist, so the reason reaches the metric
// boundary as a refusal rather than the R prompt as an error.
class gradient_refusal : public std::runtime_error {
public:
  explicit gradient_refusal(const std::string& reason)
    : std::runtime_error(reason) {
    site.reason = reason;
  }
  refusal site;
};

// One census metric's gradient. A refused metric's whole row is not-a-number and
// `why` says what happened: a sum has no defined value with an undefined term, so
// refusal is metric-level and carries no localisation within a metric.
struct census_gradient {
  // metric-major, one column per registered parameter
  std::vector<std::vector<double>> gradient;
  std::vector<refusal> why;
};

}  // namespace plant

#endif
