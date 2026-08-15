// -*-c++-*-
#ifndef PLANT_PLANT_CANOPY_SHAPE_H_
#define PLANT_PLANT_CANOPY_SHAPE_H_

#include <cmath>
#include <string>
#include <stdexcept>
#include <utility>

namespace plant {

// How the crown intercepts light. Resolved once per strategy in
// prepare_strategy() (string -> enum), never compared on the hot path.
//
// All except FlatTopBox share the same per-plant competition contribution (the
// smooth Yokozawa leaf-area profile Q); they differ in how a plant's own
// assimilation is computed and how the patch light profile is built:
//   DeepCrown    - assimilation integrated over crown depth against the smooth
//                  light profile: the leaf-area-weighted mean of the (concave)
//                  photosynthetic rate. The original plant behaviour.
//   MeanLight    - integrate the *light* over crown depth to a leaf-area-
//                  weighted mean, then a single photosynthesis evaluation on
//                  that mean light. Partway between DeepCrown and CrownCentre: it
//                  captures the mean light exactly but ignores the curvature of
//                  photosynthesis across the within-crown light distribution.
//   CrownCentre  - identical light profile to DeepCrown, but assimilation is a
//                  single evaluation of the light at the crown centre
//                  (z = H*eta_c) rather than any integral over depth.
//   FlatTopBox   - like CrownCentre for assimilation, but the plant's *competition*
//                  contribution is also collapsed into the thin crown-centre
//                  layer (a hard step: full shade below z = H*eta_c, none above)
//                  instead of the smooth Yokozawa profile. A deliberately naive
//                  variant: it casts shade *incorrectly*, so the patch light
//                  profile is discontinuous and the light-environment spline
//                  cannot be built -- the model does not run. See the vignette.
//   FlatTopSoftBox - a runnable version of FlatTopBox: the step competition is
//                  smoothed into a continuous C1 drop concentrated near the
//                  crown centre (so the light environment can be built), but the
//                  shape is still wrong (box-like, not the gradual Yokozawa
//                  taper). It runs but gives a biased fitness landscape -- the
//                  point being that a *wrong* competition profile yields wrong
//                  evolutionary predictions even when it is numerically fine.
//   PPA          - perfect-plasticity approximation: the patch light profile is
//                  built as a *stepped* function (cumulative leaf area floored
//                  into discrete canopy layers); assimilation then reads that
//                  stepped profile at the crown centre, as CrownCentre does. See
//                  FF16_Environment::compute_environment.
enum class ShadingModel {
  DeepCrown, MeanLight, CrownCentre, FlatTopBox, FlatTopSoftBox, PPA
};

inline ShadingModel shading_model_from_string(const std::string& name) {
  if (name == "deep-crown") {
    return ShadingModel::DeepCrown;
  } else if (name == "mean-light") {
    return ShadingModel::MeanLight;
  } else if (name == "crown-centre") {
    return ShadingModel::CrownCentre;
  } else if (name == "flat-top-box") {
    return ShadingModel::FlatTopBox;
  } else if (name == "flat-top-soft-box") {
    return ShadingModel::FlatTopSoftBox;
  } else if (name == "ppa") {
    return ShadingModel::PPA;
  }
  throw std::invalid_argument("Unknown shading model: " + name);
}

// As above, but an empty string selects the supplied per-strategy default
// (the shared Control default is "" so each strategy picks its own).
inline ShadingModel shading_model_from_string(const std::string& name,
                                              ShadingModel fallback) {
  if (name.empty()) {
    return fallback;
  }
  return shading_model_from_string(name);
}

// Canopy profile used by the FF16/TF24/K93 strategies. The equations follow
// the Yokozawa (1995) foliage-profile model, written in terms of the
// height-normalised coordinate u = z / H:
//
//   q(z, H) = 2 eta (1 - u^eta) u^eta / z
//   Q(z, H) = (1 - u^eta)^2
//   Qp(x, H) = (1 - sqrt(x))^(1 / eta) H
//
// In the solvers these functions are called very frequently from two hot
// paths: crown assimilation quadrature and competition/environment rebuilds.
// The primary q() and Q() methods therefore take the height-normalised ratio
// u = z / H directly, so callers can hoist the z / H division out of inner
// loops. The *_from_height() helpers keep the full-height form available for
// less performance-sensitive code and for reading the original equations.
// initialise() selects an eta-specialised multiplication chain once, and
// caches 1 / eta for Qp(), avoiding repeated generic pow setup where possible
// while preserving the original q(z, H), Q(z, H), and Qp(x, H) semantics. When a
// model can choose eta without changing its intended biology, prefer one of
// the specialised values below (1, 2, 4, 8, 10, 12); other eta values still
// work, but fall back to std::pow().

class CanopyShape {
public:
  CanopyShape()
    : eta_(12.0), eta_inverse_(1.0 / 12.0), eta_c_(eta_c(12.0)),
      pow_eta_(&pow_eta_12), leaf_above_(&leaf_above_deep),
      leaf_above_q_(&leaf_above_q_deep) {
  }

  explicit CanopyShape(double eta) {
    initialise(eta);
  }

  void initialise(double eta, ShadingModel shading_model = ShadingModel::DeepCrown) {
    eta_ = eta;
    eta_inverse_ = 1.0 / eta;
    eta_c_ = eta_c(eta);
    pow_eta_ = select_pow_eta(eta);
    // Most models cast shade via the smooth Yokozawa Q (leaf_area_above == Q).
    // FlatTopBox collapses it to a hard step; FlatTopSoftBox to a smoothed step.
    switch (shading_model) {
    case ShadingModel::FlatTopBox:
      leaf_above_ = &leaf_above_box;
      leaf_above_q_ = &leaf_above_q_box;
      break;
    case ShadingModel::FlatTopSoftBox:
      leaf_above_ = &leaf_above_softbox;
      leaf_above_q_ = &leaf_above_q_softbox;
      break;
    default:
      leaf_above_ = &leaf_above_deep;
      leaf_above_q_ = &leaf_above_q_deep;
      break;
    }
  }

  // [eqn 11] Fraction of projected leaf area above the height-normalised
  // coordinate u = z / H -- the shading a plant casts at u. Smooth Yokozawa Q
  // for every model except FlatTopBox, which uses a step at the crown centre.
  // Bound once in initialise(), so the competition hot path makes one predicted
  // indirect call with no branch.
  double leaf_area_above(double z_over_height) const {
    return leaf_above_(*this, z_over_height);
  }

  // Undefined at the crown base, where z is 0: use q_from_height, which carries
  // the height the limit there needs.
  double q(double z_over_height, double z) const {
    const double u_eta = pow_eta_(z_over_height, eta_);
    return 2.0 * eta_ * (1.0 - u_eta) * u_eta / z;
  }

  double q_from_height(double z, double height) const {
    // The 1 / z above is 0 / 0 at the crown base, so take the limit: 0 for every
    // eta above 1, and 2 / height at eta = 1. The light field's lowest knot asks
    // for exactly this, and the crown integral never does.
    if (z <= 0.0) {
      return eta_ == 1.0 ? 2.0 / height : 0.0;
    }
    return q(z / height, z);
  }

  double Q(double z_over_height) const {
    if (z_over_height > 1.0) {
      return 0.0;
    }
    const double tmp = 1.0 - pow_eta_(z_over_height, eta_);
    return tmp * tmp;
  }

  double Q_from_height(double z, double height) const {
    if (z > height) {
      return 0.0;
    }
    return Q(z / height);
  }

  // Q and q from the one u^eta both need, for the callers that read the smooth
  // profile directly rather than through leaf_area_above.
  std::pair<double, double> Q_and_q(double z_over_height,
                                    double height_inverse) const {
    return {Q(z_over_height), leaf_above_q_deep(*this, z_over_height, height_inverse)};
  }

  // The shading a crown casts at u = z / H and its q there: the pair the light
  // field reads at every knot. Defined at the crown base, where q_from_height
  // takes a limit.
  std::pair<double, double> leaf_area_above_and_q(double z_over_height,
                                                  double height_inverse) const {
    return {leaf_above_(*this, z_over_height),
            leaf_above_q_(*this, z_over_height, height_inverse)};
  }

  // Q(u) = 1 - 2 u^eta + u^(2 eta), so a reduction over cohorts at many heights
  // needs three running sums rather than one Q per cohort per height.
  //
  // The sums hold the cohorts whose crowns reach `at`, scaled by `at` itself.
  // Carrying them unscaled loses the sum's own small terms: the factor h^-eta at
  // eta 12 spans nineteen orders over one stand's heights. They carry `at` so a
  // caller cannot descend by a ratio taken from some other height -- the result
  // would be a plausible field with no symptom.
  struct ProfileSums {
    double flat = 0.0;
    double first = 0.0;
    double second = 0.0;
    double at = 0.0;
  };

  // Whether Q is a polynomial in u^eta, which is what the sums above rest on.
  // The box profiles are not, and the hard box has no vertical slope at all.
  bool profile_expands_in_u_eta() const {
    return leaf_above_ == &leaf_above_deep;
  }

  // Whether a light field can be built from this profile. Named profiles only:
  // a grid samples a step at the crown centre as a ramp one span wide and
  // returns numbers, so a profile that has not been checked for continuity must
  // not reach a field by default.
  bool profile_builds_a_field() const {
    return leaf_above_ == &leaf_above_deep || leaf_above_ == &leaf_above_softbox;
  }

  // Rescale the sums to a lower height. The ratio is taken from what they hold,
  // so it cannot come from anywhere else; every admitted cohort's u falls by it.
  void descend_to(ProfileSums& s, double z) const {
    if (s.at > 0.0) {
      const double r = pow_eta_(z / s.at, eta_);
      s.first *= r;
      s.second *= r * r;
    }
    s.at = z;
  }

  // Admit one cohort of amplitude `scale` whose crown top is `height`, at the
  // height the sums now hold. Every term is at most `scale`, because a cohort is
  // admitted only once that height has fallen to its own crown top.
  void admit(ProfileSums& s, double scale, double height) const {
    const double u_eta = pow_eta_(s.at / height, eta_);
    s.flat += scale;
    s.first += scale * u_eta;
    s.second += scale * u_eta * u_eta;
  }

  // The pair leaf_area_above_and_q returns, summed over everything admitted.
  // Undefined at the crown base for the reason q is.
  std::pair<double, double> Q_and_q_from_sums(const ProfileSums& s) const {
    const double gain = 2.0 * eta_ / s.at;
    return {s.flat - 2.0 * s.first + s.second, (s.first - s.second) * gain};
  }

  double Qp(double x, double height) const {
    return std::pow(1.0 - std::sqrt(x), eta_inverse_) * height;
  }

  // [eqn 12] Crown-centre coordinate u = z / H. Static because the strategies
  // need the same number for their sapwood and conductance terms, and one
  // formula is better than three.
  static double eta_c(double eta) {
    return 1.0 - 2.0 / (1.0 + eta) + 1.0 / (1.0 + 2.0 * eta);
  }

private:
  typedef double (*pow_eta_fn)(double, double);
  typedef double (*leaf_above_fn)(const CanopyShape&, double);
  typedef double (*leaf_above_q_fn)(const CanopyShape&, double, double);

  // Smooth Yokozawa profile -- the correct shading a crown casts.
  static double leaf_above_deep(const CanopyShape& c, double z_over_height) {
    return c.Q(z_over_height);
  }

  // FlatTopBox: all leaf area collapsed into the thin crown-centre layer, so the
  // crown fully shades everything below z = H*eta_c and nothing above. A step.
  static double leaf_above_box(const CanopyShape& c, double z_over_height) {
    return z_over_height < c.eta_c_ ? 1.0 : 0.0;
  }

  // FlatTopSoftBox: the hard step softened into a monotone C1 drop, full shade up
  // to lo = max(0, 2*eta_c - 1) then a cubic-smoothstep fall to zero at the crown
  // top (so the transition is centred on the crown centre eta_c and the profile
  // is continuous -- buildable -- but still box-like, not the Yokozawa taper).
  static double leaf_above_softbox(const CanopyShape& c, double z_over_height) {
    const double lo = c.eta_c_ > 0.5 ? 2.0 * c.eta_c_ - 1.0 : 0.0;
    if (z_over_height <= lo) {
      return 1.0;
    }
    if (z_over_height >= 1.0) {
      return 0.0;
    }
    const double t = (z_over_height - lo) / (1.0 - lo);
    return 1.0 - t * t * (3.0 - 2.0 * t);
  }

  // Divided by the height rather than by z, which is what makes the crown base a
  // value where q(u, z) is 0 / 0. Same limit q_from_height takes.
  static double leaf_above_q_deep(const CanopyShape& c, double z_over_height,
                                  double height_inverse) {
    if (z_over_height > 1.0) {
      return 0.0;
    }
    if (z_over_height <= 0.0) {
      return c.eta_ == 1.0 ? 2.0 * height_inverse : 0.0;
    }
    const double u_eta = c.pow_eta_(z_over_height, c.eta_);
    return 2.0 * c.eta_ * (1.0 - u_eta) * u_eta * height_inverse / z_over_height;
  }

  // A step is flat either side of the crown centre, and the jump there is not a
  // number a field carries.
  static double leaf_above_q_box(const CanopyShape&, double, double) {
    return 0.0;
  }

  static double leaf_above_q_softbox(const CanopyShape& c, double z_over_height,
                                     double height_inverse) {
    const double lo = c.eta_c_ > 0.5 ? 2.0 * c.eta_c_ - 1.0 : 0.0;
    if (z_over_height <= lo || z_over_height >= 1.0) {
      return 0.0;
    }
    const double t = (z_over_height - lo) / (1.0 - lo);
    return 6.0 * t * (1.0 - t) * height_inverse / (1.0 - lo);
  }

  static pow_eta_fn select_pow_eta(double eta) {
    if (eta == 1.0) {
      return &pow_eta_1;
    } else if (eta == 2.0) {
      return &pow_eta_2;
    } else if (eta == 4.0) {
      return &pow_eta_4;
    } else if (eta == 8.0) {
      return &pow_eta_8;
    } else if (eta == 10.0) {
      return &pow_eta_10;
    } else if (eta == 12.0) {
      return &pow_eta_12;
    } else {
      return &pow_eta_general;
    }
  }

  static double pow_eta_general(double u, double eta) {
    return std::pow(u, eta);
  }

  static double pow_eta_1(double u, double) {
    return u;
  }

  static double pow_eta_2(double u, double) {
    return u * u;
  }

  static double pow_eta_4(double u, double) {
    const double u2 = u * u;
    return u2 * u2;
  }

  static double pow_eta_8(double u, double) {
    const double u2 = u * u;
    const double u4 = u2 * u2;
    return u4 * u4;
  }

  static double pow_eta_10(double u, double) {
    const double u2 = u * u;
    const double u4 = u2 * u2;
    const double u8 = u4 * u4;
    return u8 * u2;
  }

  static double pow_eta_12(double u, double) {
    const double u2 = u * u;
    const double u4 = u2 * u2;
    const double u8 = u4 * u4;
    return u8 * u4;
  }

  double eta_;
  double eta_inverse_;
  double eta_c_;
  pow_eta_fn pow_eta_;
  leaf_above_fn leaf_above_;
  leaf_above_q_fn leaf_above_q_;
};

}

#endif
