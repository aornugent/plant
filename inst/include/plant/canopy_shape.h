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
    : eta_(12.0), eta_inverse_(1.0 / 12.0), eta_c_(compute_eta_c(12.0)),
      pow_eta_kind_(PowEtaKind::P12), leaf_above_kind_(LeafAboveKind::Deep) {
  }

  explicit CanopyShape(double eta) {
    initialise(eta);
  }

  void initialise(double eta, ShadingModel shading_model = ShadingModel::DeepCrown) {
    eta_ = eta;
    eta_inverse_ = 1.0 / eta;
    eta_c_ = compute_eta_c(eta);
    pow_eta_kind_ = select_pow_eta_kind(eta);
    // Most models cast shade via the smooth Yokozawa Q (leaf_area_above == Q).
    // FlatTopBox collapses it to a hard step; FlatTopSoftBox to a smoothed step.
    switch (shading_model) {
    case ShadingModel::FlatTopBox:     leaf_above_kind_ = LeafAboveKind::Box;     break;
    case ShadingModel::FlatTopSoftBox: leaf_above_kind_ = LeafAboveKind::SoftBox; break;
    default:                           leaf_above_kind_ = LeafAboveKind::Deep;    break;
    }
  }

  // The profile methods below are templated on the position scalar Z: Z = double
  // is the production path (bit-identical -- the eta-specialised multiply chains
  // and the profile expressions are unchanged, only the double(*)(...) dispatch
  // is replaced by a switch on a stored kind), while an active Z tapes the crown
  // integral's moving quadrature nodes (plan §4.4). The shape coefficients
  // (eta_, eta_c_) stay double: the profile is not differentiated w.r.t. eta.

  // [eqn 11] Fraction of projected leaf area above the height-normalised
  // coordinate u = z / H -- the shading a plant casts at u. Smooth Yokozawa Q
  // for every model except FlatTopBox (a hard step) and FlatTopSoftBox (smoothed).
  template <class Z>
  Z leaf_area_above(Z z_over_height) const {
    switch (leaf_above_kind_) {
    case LeafAboveKind::Box:     return z_over_height < eta_c_ ? Z(1.0) : Z(0.0);
    case LeafAboveKind::SoftBox: return leaf_above_softbox(z_over_height);
    default:                     return Q(z_over_height); // Deep
    }
  }

  template <class Z>
  Z q(Z z_over_height, Z z) const {
    const Z u_eta = pow_eta(z_over_height);
    return 2.0 * eta_ * (1.0 - u_eta) * u_eta / z;
  }

  template <class Z>
  Z q_from_height(Z z, Z height) const {
    return q(z / height, z);
  }

  template <class Z>
  Z Q(Z z_over_height) const {
    if (z_over_height > 1.0) {
      return Z(0.0);
    }
    const Z tmp = 1.0 - pow_eta(z_over_height);
    return tmp * tmp;
  }

  template <class Z>
  Z Q_from_height(Z z, Z height) const {
    if (z > height) {
      return Z(0.0);
    }
    return Q(z / height);
  }

  // Qp is an R-facing/diagnostic inverse, off the differentiated rate path, so
  // it stays double.
  double Qp(double x, double height) const {
    return std::pow(1.0 - std::sqrt(x), eta_inverse_) * height;
  }

private:
  enum class PowEtaKind { P1, P2, P4, P8, P10, P12, General };
  enum class LeafAboveKind { Deep, Box, SoftBox };

  static double compute_eta_c(double eta) {
    return 1.0 - 2.0 / (1.0 + eta) + 1.0 / (1.0 + 2.0 * eta);
  }

  // FlatTopSoftBox: the hard step softened into a monotone C1 drop, full shade up
  // to lo = max(0, 2*eta_c - 1) then a cubic-smoothstep fall to zero at the crown
  // top (so the transition is centred on the crown centre eta_c and the profile
  // is continuous -- buildable -- but still box-like, not the Yokozawa taper).
  template <class Z>
  Z leaf_above_softbox(Z z_over_height) const {
    const double lo = eta_c_ > 0.5 ? 2.0 * eta_c_ - 1.0 : 0.0;
    if (z_over_height <= lo) {
      return Z(1.0);
    }
    if (z_over_height >= 1.0) {
      return Z(0.0);
    }
    const Z t = (z_over_height - lo) / (1.0 - lo);
    return 1.0 - t * t * (3.0 - 2.0 * t);
  }

  static PowEtaKind select_pow_eta_kind(double eta) {
    if (eta == 1.0) {
      return PowEtaKind::P1;
    } else if (eta == 2.0) {
      return PowEtaKind::P2;
    } else if (eta == 4.0) {
      return PowEtaKind::P4;
    } else if (eta == 8.0) {
      return PowEtaKind::P8;
    } else if (eta == 10.0) {
      return PowEtaKind::P10;
    } else if (eta == 12.0) {
      return PowEtaKind::P12;
    } else {
      return PowEtaKind::General;
    }
  }

  // u^eta via the eta-specialised multiply chains (bit-identical to the former
  // pow_eta_* function pointers); the general case is std::pow, made ADL-friendly
  // so an active Z resolves to xad::pow.
  template <class Z>
  Z pow_eta(Z u) const {
    switch (pow_eta_kind_) {
    case PowEtaKind::P1:  return u;
    case PowEtaKind::P2:  return u * u;
    case PowEtaKind::P4:  { const Z u2 = u * u; return u2 * u2; }
    case PowEtaKind::P8:  { const Z u2 = u * u; const Z u4 = u2 * u2; return u4 * u4; }
    case PowEtaKind::P10: { const Z u2 = u * u; const Z u4 = u2 * u2; const Z u8 = u4 * u4; return u8 * u2; }
    case PowEtaKind::P12: { const Z u2 = u * u; const Z u4 = u2 * u2; const Z u8 = u4 * u4; return u8 * u4; }
    default:             { using std::pow; return pow(u, eta_); }
    }
  }

  double eta_;
  double eta_inverse_;
  double eta_c_;
  PowEtaKind pow_eta_kind_;
  LeafAboveKind leaf_above_kind_;
};

}

#endif
