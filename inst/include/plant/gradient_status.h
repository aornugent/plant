// -*-c++-*-
#ifndef PLANT_PLANT_GRADIENT_STATUS_H_
#define PLANT_PLANT_GRADIENT_STATUS_H_

#include <stdexcept>
#include <string>
#include <vector>

namespace plant {

// What a gradient entry is, beside its number. A bare double cannot say whether
// it is an answer, a zero the model means, or a slot no sweep reached -- and an
// exact zero is the signature of the last of those far more often than of the
// second, so a reader shown only numbers reads a missing accumulator as an
// ecological finding.
struct gradient_status {
  enum class Kind {
    answered,
    // Live parameter, and the model's own answer here is zero: complementary
    // slackness at an interior optimum, a constraint the operating point is not
    // sitting on. Becomes non-zero at a state where the constraint binds.
    zero_slack,
    // The parameter reaches no equation this metric reads. Zero for the whole
    // run, whatever the state, and zero on any other trajectory too.
    zero_structural,
    // Exactly zero, and nothing declares why. This is a finding rather than an
    // answer: a missing accumulator and a true insensitivity are the same
    // number, so the entry is marked instead of being read as either.
    zero_undeclared,
    // No gradient exists at some state the trajectory visited. The number
    // beside this is not an answer and must not be read as one.
    refused
  };

  Kind kind = Kind::answered;
  // Set only when refused. Names the quantity and the constraint, from the
  // frame that found it.
  std::string reason;
  // Where it was found, filled in by the frames that know: 1-based species,
  // node 0 for the inflow boundary and 1-based for an introduced node, and the
  // range of recorded steps the sweep was covering. -1 is "not known here".
  int species = -1;
  int node = -1;
  int step_first = -1;
  int step_last = -1;

  static const char* kind_name(Kind k) {
    switch (k) {
    case Kind::answered:        return "answered";
    case Kind::zero_slack:      return "zero-slack";
    case Kind::zero_structural: return "zero-structural";
    case Kind::zero_undeclared: return "zero-undeclared";
    case Kind::refused:         return "refused";
    }
    return "unknown";
  }
};

// Thrown where a gradient row does not exist, so the reason reaches the metric
// boundary as a status rather than the R prompt as an error. It carries a
// status rather than only a message because the location is known in pieces:
// the leaf knows the reason, the block loop knows which plant, and the sweep
// knows which steps.
class gradient_refusal : public std::runtime_error {
public:
  explicit gradient_refusal(const std::string& reason)
    : std::runtime_error(reason) {
    status.kind = gradient_status::Kind::refused;
    status.reason = reason;
  }
  gradient_status status;
};

// One census metric's gradient and what each entry of it is. The two travel
// together so that a caller cannot obtain a row of numbers without the reading
// that says whether they are answers.
struct census_gradient {
  // metric-major, one column per registered parameter
  std::vector<std::vector<double>> gradient;
  std::vector<std::vector<gradient_status>> status;
};

}  // namespace plant

#endif
