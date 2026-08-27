// -*-c++-*-
#ifndef PLANT_PLANT_CENSUS_GRADIENT_H_
#define PLANT_PLANT_CENSUS_GRADIENT_H_

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
  // 1-based. A refusal is recorded on the species whose strategy holds it, so one
  // that happened always names one.
  int species = -1;

  bool happened() const { return !reason.empty(); }
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
