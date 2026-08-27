// -*-c++-*-
#ifndef PLANT_PLANT_GRADIENT_REFUSAL_H_
#define PLANT_PLANT_GRADIENT_REFUSAL_H_

#include <stdexcept>
#include <string>
#include <vector>

namespace plant {

// Why a metric has no gradient, and which species it was found on. A refusal is
// the only thing a number needs said beside it: a parameter no gradient exists
// for cannot be asked for, so it never reaches an answer, and every column that
// does reach one carries a number the sweep computed.
//
// The species is as fine as the grain gets: a row that could not be supplied is
// an intermediate of a recording spanning six stages and every cohort in them, so
// nothing below the species has a component to attribute it to.
struct refusal {
  std::string reason;
  // 1-based, and -1 where the frame that found it knew no species.
  int species = -1;

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
