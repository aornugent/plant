// -*-c++-*-
#ifndef PLANT_PLANT_FF16_PRODUCTION_KERNEL_H_
#define PLANT_PLANT_FF16_PRODUCTION_KERNEL_H_

#include <vector>
#include <cstddef>
#include <cmath>   // std::pow; XAD provides pow for active types via ADL

// Scalar-templated core of the FF16 single-plant net-mass-production chain
// (#472 scope B / traitecoevo/plant#537, Milestone A). The pieces below [eqn
// 12-15] are elementary arithmetic, so templating them on the scalar S makes net
// production differentiable w.r.t. traits by reverse-mode AD with no special
// handling. They are the SINGLE SOURCE OF TRUTH: FF16_Strategy's double methods
// delegate to them (so the existing FF16 test suite validates faithfulness), and
// the AD calibration path instantiates them with an active scalar.

namespace plant {

// [eqn 2] Leaf area as a function of height (inverse of [eqn 3]). Carries the
// allometric traits a_l1, a_l2, so it is the entry point through which a trait
// reaches both the mass cascade and the resident competition / light spline.
template <typename S>
S ff16_area_leaf(S a_l1, S a_l2, S height) {
  using std::pow;
  return pow(height / a_l1, 1.0 / a_l2);
}

// [eqn 12] Photosynthetic rate per leaf area; x is openness in [0,1].
template <typename S>
S ff16_assimilation_leaf(S a_p1, S a_p2, S x) {
  return a_p1 * x / (x + a_p2);
}

// [eqn 13] Total maintenance respiration (linear in the mass cascade).
template <typename S>
S ff16_respiration(S mass_leaf, S mass_sapwood, S mass_bark, S mass_root,
                   S r_l, S r_s, S r_b, S r_r) {
  return r_l * mass_leaf + r_b * mass_bark + r_s * mass_sapwood + r_r * mass_root;
}

// [eqn 14] Total turnover.
template <typename S>
S ff16_turnover(S mass_leaf, S mass_bark, S mass_sapwood, S mass_root,
                S k_l, S k_b, S k_s, S k_r) {
  return k_l * mass_leaf + k_b * mass_bark + k_s * mass_sapwood + k_r * mass_root;
}

// [eqn 15] Net production from assimilation/respiration/turnover.
template <typename S>
S ff16_net_production_A(S a_bio, S a_y, S assimilation, S respiration, S turnover) {
  return a_bio * a_y * (assimilation - respiration) - turnover;
}

// The FF16 parameters the crown-top net-production chain reads, plus the
// prepare_strategy()-derived eta_c. Used by the all-in-one convenience below
// (the AD calibration entry point); the hot double path uses the per-piece
// functions directly with pars.* members.
template <typename S>
struct FF16ProdPars {
  S lma, rho, theta, a_b1, a_r1, eta_c;
  S a_p1, a_p2;
  S r_l, r_s, r_b, r_r;
  S k_l, k_b, k_s, k_r;
  S a_bio, a_y;
  // Allometry + allocation parameters for the height-growth rate (Milestone C).
  S a_l1, a_l2;       // height <-> leaf-area allometry [eqn 2/3]
  S a_f1, a_f2, hmat; // reproduction-allocation logistic [eqn 16]
};

// Whole single-plant net production under the CROWN-TOP assimilation variant (a
// single light evaluation light_E at the canopy top; no quadrature -- the
// adaptive crown integral is a separate frozen-replay concern, Milestone B).
// Mirrors FF16_Strategy's mass cascade + assimilation_crown_top + the pieces
// above exactly. light_E is the (fixed, double on the resident path) light
// fraction at the crown.
// Mass cascade -> respiration/turnover -> net production, GIVEN the assimilation
// rate. Shared by every assimilation variant (crown-top, deep-crown, ...), so
// the variants differ only in how they compute `assimilation`.
template <typename S>
S ff16_net_from_components(const FF16ProdPars<S>& p, S height, S area_leaf,
                           S assimilation) {
  const S mass_leaf    = area_leaf * p.lma;
  const S area_sapwood = area_leaf * p.theta;
  const S mass_sapwood = area_sapwood * height * p.eta_c * p.rho;
  const S area_bark    = p.a_b1 * area_leaf * p.theta;
  const S mass_bark    = area_bark * height * p.eta_c * p.rho;
  const S mass_root    = p.a_r1 * area_leaf;

  const S respiration = ff16_respiration(mass_leaf, mass_sapwood, mass_bark, mass_root,
                                         p.r_l, p.r_s, p.r_b, p.r_r);
  const S turnover    = ff16_turnover(mass_leaf, mass_bark, mass_sapwood, mass_root,
                                      p.k_l, p.k_b, p.k_s, p.k_r);
  return ff16_net_production_A(p.a_bio, p.a_y, assimilation, respiration, turnover);
}

template <typename S>
S ff16_net_mass_production_crown_top(const FF16ProdPars<S>& p,
                                     S height, S area_leaf, S light_E) {
  const S assimilation = area_leaf * ff16_assimilation_leaf(p.a_p1, p.a_p2, light_E);
  return ff16_net_from_components(p, height, area_leaf, assimilation);
}

// Deep-crown assimilation as a FROZEN-REPLAY weighted sum (#472 scope B,
// Milestone B). The production model's default assimilation integrates
// assimilation_leaf(light(z)) * q(z/height, z) over crown depth with adaptive
// Gauss-Kronrod. Per the two-pass plan we do NOT differentiate the adaptive
// controller: pass 1 (double) discovers the nodes z_j and the COMBINED weights
// wq_j = w_j * q(z_j/height, z_j) (the leaf-area density q is constant in the
// physiology traits, so it folds into the frozen weight); pass 2 replays
//   A = area_leaf * sum_j wq_j * assimilation_leaf(a_p1, a_p2, light(z_j)).
// `light` is any callable z -> S (e.g. an AD-capable resident light spline,
// odelia::interpolator::basic_interpolator<S>), so A is differentiable w.r.t.
// a_p1/a_p2 and w.r.t. light's knot values (the resident self-shading coupling).
template <typename S, typename LightFn>
S ff16_assimilation_deep_crown_replay(S a_p1, S a_p2, S area_leaf,
                                      const std::vector<double>& z,
                                      const std::vector<double>& wq,
                                      LightFn&& light) {
  S A = S(0.0);
  for (std::size_t j = 0; j < z.size(); ++j) {
    A += wq[j] * ff16_assimilation_leaf(a_p1, a_p2, S(light(z[j])));
  }
  return area_leaf * A;
}

// ---------------------------------------------------------------------------
// Height-growth rate pieces (#472 scope B, Milestone C). Mirror
// FF16_Strategy::{fraction_allocation_growth, dheight_darea_leaf,
// dmass_*_darea_leaf, darea_leaf_dmass_live} and the dheight/dt assembly in
// compute_rates. Elementary, so dheight/dt is differentiable w.r.t. height
// (A1 -- the exact growth-rate gradient now done by finite difference in
// Node::growth_rate_gradient) and w.r.t. traits.
// ---------------------------------------------------------------------------

// [eqn 16] Fraction of production allocated to growth = 1 - reproduction.
template <typename S>
S ff16_fraction_allocation_growth(S a_f1, S a_f2, S hmat, S height) {
  using std::exp;
  return 1.0 - a_f1 / (1.0 + exp(a_f2 * (1.0 - height / hmat)));
}

// d(height)/d(area_leaf): derivative of the [eqn 2] allometry.
template <typename S>
S ff16_dheight_darea_leaf(S a_l1, S a_l2, S area_leaf) {
  using std::pow;
  return a_l1 * a_l2 * pow(area_leaf, a_l2 - 1.0);
}

// d(area_leaf)/d(mass_live): reciprocal of the summed per-component mass
// derivatives (leaf + sapwood + bark + root) w.r.t. area_leaf.
template <typename S>
S ff16_darea_leaf_dmass_live(const FF16ProdPars<S>& p, S area_leaf) {
  using std::pow;
  const S dmass_leaf    = p.lma;                                   // d(area_leaf*lma)
  const S dmass_sapwood = p.rho * p.eta_c * p.a_l1 * p.theta *
                          (p.a_l2 + 1.0) * pow(area_leaf, p.a_l2);
  const S dmass_bark    = p.a_b1 * dmass_sapwood;
  const S dmass_root    = p.a_r1;
  return 1.0 / (dmass_leaf + dmass_sapwood + dmass_bark + dmass_root);
}

// dheight/dt for a plant of the given height under the CROWN-TOP assimilation
// variant, in light light_E. Returns 0 when net production is non-positive
// (the compute_rates growth clamp). area_leaf is derived from height so the
// gradient w.r.t. height flows through the whole chain.
template <typename S>
S ff16_height_dt_crown_top(const FF16ProdPars<S>& p, S height, S light_E) {
  const S area_leaf = ff16_area_leaf(p.a_l1, p.a_l2, height);
  const S net = ff16_net_mass_production_crown_top(p, height, area_leaf, light_E);
  if (net <= 0.0) return S(0.0);
  const S frac_growth = ff16_fraction_allocation_growth(p.a_f1, p.a_f2, p.hmat, height);
  const S darea_dmass = ff16_darea_leaf_dmass_live(p, area_leaf);
  const S area_leaf_dt = net * frac_growth * darea_dmass;
  return ff16_dheight_darea_leaf(p.a_l1, p.a_l2, area_leaf) * area_leaf_dt;
}

}  // namespace plant

#endif
