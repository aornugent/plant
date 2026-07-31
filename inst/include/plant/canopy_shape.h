// -*-c++-*-
#ifndef PLANT_PLANT_CANOPY_SHAPE_H_
#define PLANT_PLANT_CANOPY_SHAPE_H_

#include <cmath>
#include <string>
#include <stdexcept>

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
//   q(z, H) = 2 eta (1 - u^eta) u^(eta - 1) / H
//   Q(z, H) = (1 - u^eta)^2
//   Qp(x, H) = (1 - sqrt(x))^(1 / eta) H
//
// In the solvers these functions are called very frequently from two hot
// paths: crown assimilation quadrature and competition/environment rebuilds.
// The primary q() and Q() methods therefore take the height-normalised ratio
// u = z / H directly, so callers can hoist the z / H division out of inner
// loops; q() takes 1 / H alongside it for the same reason.
//
// q is written over u^(eta - 1) / H rather than the equivalent u^eta / z, which
// divides by zero at the crown base. The u -> 0 limit is 0 for eta > 1 and 2 / H
// at eta = 1, and u^(eta - 1) supplies both with no branch.
//
// The *_from_height() helpers keep the full-height form available for
// less performance-sensitive code and for reading the original equations.
// initialise() picks a multiplication chain for u^(eta - 1) once, and caches
// 1 / eta for Qp(). Q multiplies that chain by u to get u^eta, so one chain
// serves both and there is a single place where the exponent is applied. When a
// model can choose eta without changing its intended biology, prefer one of the
// values with a chain (1, 2, 4, 8, 10, 12); any other eta works and calls
// std::pow.

class CanopyShape {
public:
  CanopyShape()
    : eta_(12.0), eta_inverse_(1.0 / 12.0), eta_minus_1_(11.0),
      eta_c_(eta_c(12.0)), chain_(&u_pow_11), leaf_above_(&leaf_above_deep) {
  }

  explicit CanopyShape(double eta) {
    initialise(eta);
  }

  void initialise(double eta, ShadingModel shading_model = ShadingModel::DeepCrown) {
    eta_ = eta;
    eta_inverse_ = 1.0 / eta;
    eta_minus_1_ = eta - 1.0;
    eta_c_ = eta_c(eta);
    chain_ = select_chain(eta);
    // Most models cast shade via the smooth Yokozawa Q (leaf_area_above == Q).
    // FlatTopBox collapses it to a hard step; FlatTopSoftBox to a smoothed step.
    switch (shading_model) {
    case ShadingModel::FlatTopBox:     leaf_above_ = &leaf_above_box;     break;
    case ShadingModel::FlatTopSoftBox: leaf_above_ = &leaf_above_softbox; break;
    default:                           leaf_above_ = &leaf_above_deep;    break;
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

  // u^(eta - 1) and u^eta. One chain serves both, because the second is the
  // first times u -- and q needs the first, which is what keeps it finite at
  // u = 0.
  double u_pow_eta_minus_1(double u) const { return chain_(u, eta_minus_1_); }
  double u_pow_eta(double u) const { return u_pow_eta_minus_1(u) * u; }

  // [eqn 9] Leaf area density at u = z / H.
  double q(double z_over_height, double height_inverse) const {
    const double u_eta_minus_1 = u_pow_eta_minus_1(z_over_height);
    return 2.0 * eta_ * (1.0 - u_eta_minus_1 * z_over_height) * u_eta_minus_1 *
           height_inverse;
  }

  double q_from_height(double z, double height) const {
    const double height_inverse = 1.0 / height;
    return q(z * height_inverse, height_inverse);
  }

  double Q(double z_over_height) const {
    if (z_over_height > 1.0) {
      return 0.0;
    }
    const double tmp = 1.0 - u_pow_eta(z_over_height);
    return tmp * tmp;
  }

  double Q_from_height(double z, double height) const {
    if (z > height) {
      return 0.0;
    }
    return Q(z / height);
  }

  // [eqn 12] Crown-centre coordinate u = z / H. Static because the strategies
  // need the same number for their sapwood and conductance terms, and one
  // formula is better than three.
  static double eta_c(double eta) {
    return 1.0 - 2.0 / (1.0 + eta) + 1.0 / (1.0 + 2.0 * eta);
  }

  double Qp(double x, double height) const {
    return std::pow(1.0 - std::sqrt(x), eta_inverse_) * height;
  }

private:
  typedef double (*chain_fn)(double u, double exponent);
  typedef double (*leaf_above_fn)(const CanopyShape&, double);

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

  // u^(eta - 1), by a multiplication chain for the etas a model is likely to
  // choose. Each function is named for the exponent it computes, so eta = 12
  // selects u_pow_11. Everything else falls back to std::pow.
  static chain_fn select_chain(double eta) {
    if (eta == 1.0)  return &u_pow_0;
    if (eta == 2.0)  return &u_pow_1;
    if (eta == 4.0)  return &u_pow_3;
    if (eta == 8.0)  return &u_pow_7;
    if (eta == 10.0) return &u_pow_9;
    if (eta == 12.0) return &u_pow_11;
    return &u_pow_general;
  }

  static double u_pow_general(double u, double exponent) { return std::pow(u, exponent); }
  static double u_pow_0(double, double)  { return 1.0; }
  static double u_pow_1(double u, double) { return u; }
  static double u_pow_3(double u, double) { return u * u * u; }

  static double u_pow_7(double u, double) {
    const double u2 = u * u;
    return u2 * u2 * u2 * u;
  }

  static double u_pow_9(double u, double) {
    const double u2 = u * u;
    const double u4 = u2 * u2;
    return u4 * u4 * u;
  }

  static double u_pow_11(double u, double) {
    const double u2 = u * u;
    const double u4 = u2 * u2;
    return u4 * u4 * u2 * u;
  }

  double eta_;
  double eta_inverse_;
  double eta_minus_1_;
  double eta_c_;
  chain_fn chain_;
  leaf_above_fn leaf_above_;
};

}

#endif
