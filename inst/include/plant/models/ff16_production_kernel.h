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
  // Demographic rate parameters for the full compute_rates fill (Milestone C):
  // fecundity [eqn 17] and mortality [eqn 21].
  S omega, a_f3;       // seed mass + accessory reproduction cost (fecundity_dt)
  S d_I, a_dG1, a_dG2; // mortality: growth-independent + growth-dependent
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

// [eqn 16] Fraction of production allocated to reproduction (logistic in height).
template <typename S>
S ff16_fraction_allocation_reproduction(S a_f1, S a_f2, S hmat, S height) {
  using std::exp;
  return a_f1 / (1.0 + exp(a_f2 * (1.0 - height / hmat)));
}

// [eqn 16] Fraction of production allocated to growth = 1 - reproduction.
template <typename S>
S ff16_fraction_allocation_growth(S a_f1, S a_f2, S hmat, S height) {
  return 1.0 - ff16_fraction_allocation_reproduction(a_f1, a_f2, hmat, height);
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

// dheight/dt given net production (the compute_rates growth assembly): returns
// 0 when net is non-positive (the growth clamp), else dheight_darea_leaf *
// area_leaf_dt. Shared by every assimilation variant.
template <typename S>
S ff16_height_dt_from_net(const FF16ProdPars<S>& p, S height, S area_leaf, S net) {
  if (net <= 0.0) return S(0.0);
  const S frac_growth = ff16_fraction_allocation_growth(p.a_f1, p.a_f2, p.hmat, height);
  const S darea_dmass = ff16_darea_leaf_dmass_live(p, area_leaf);
  const S area_leaf_dt = net * frac_growth * darea_dmass;
  return ff16_dheight_darea_leaf(p.a_l1, p.a_l2, area_leaf) * area_leaf_dt;
}

// dheight/dt for a plant of the given height under the CROWN-TOP assimilation
// variant, in light light_E. area_leaf is derived from height so the gradient
// w.r.t. height flows through the whole chain.
template <typename S>
S ff16_height_dt_crown_top(const FF16ProdPars<S>& p, S height, S light_E) {
  const S area_leaf = ff16_area_leaf(p.a_l1, p.a_l2, height);
  const S net = ff16_net_mass_production_crown_top(p, height, area_leaf, light_E);
  return ff16_height_dt_from_net(p, height, area_leaf, net);
}

// The five ODE state rates FF16_Strategy::compute_rates writes, plus the net
// production aux. Scalar-templated (#472 scope B, Milestone C) so the whole
// demographic rate fill differentiates w.r.t. a trait by reverse-mode AD.
template <typename S>
struct FF16Rates {
  S net_mass_production_dt;
  S height_dt;
  S fecundity_dt;
  S area_heartwood_dt;
  S mass_heartwood_dt;
  S mortality_dt;
};

// Full FF16 compute_rates fill for the CROWN-TOP assimilation variant (single
// light evaluation light_E). Mirrors FF16_Strategy::compute_rates EXACTLY: the
// net>0 growth clamp gates the growth/fecundity/heartwood rates, and mortality
// is the [eqn 21] growth-independent + growth-dependent sum (productivity =
// net/area_leaf). `mortality_finite` is the frozen util::is_finite(cumulative
// mortality) branch -- a pass-1 (double) control-flow decision, passed in so the
// taped replay is branch-free (it never differentiates the is_finite test).
// Deep-crown differs only in how `net` is formed (the frozen-replay crown
// integral, ff16_assimilation_deep_crown_replay -> ff16_net_from_components),
// so a deep-crown fill reuses everything below by substituting that `net`.
template <typename S>
FF16Rates<S> ff16_compute_rates_crown_top(const FF16ProdPars<S>& p, S height,
                                          S light_E, bool mortality_finite) {
  const S area_leaf = ff16_area_leaf(p.a_l1, p.a_l2, height);
  const S net = ff16_net_mass_production_crown_top(p, height, area_leaf, light_E);
  FF16Rates<S> r;
  r.net_mass_production_dt = net;
  if (net > 0.0) {
    const S frac_repro = ff16_fraction_allocation_reproduction(p.a_f1, p.a_f2,
                                                               p.hmat, height);
    r.height_dt          = ff16_height_dt_from_net(p, height, area_leaf, net);
    r.fecundity_dt       = net * frac_repro / (p.omega + p.a_f3);
    const S area_sapwood = area_leaf * p.theta;            // [eqn 4]
    r.area_heartwood_dt  = p.k_s * area_sapwood;           // turnover of sapwood area
    const S mass_sapwood = area_sapwood * height * p.eta_c * p.rho;
    r.mass_heartwood_dt  = p.k_s * mass_sapwood;           // turnover_sapwood(mass)
  } else {
    r.height_dt = S(0.0); r.fecundity_dt = S(0.0);
    r.area_heartwood_dt = S(0.0); r.mass_heartwood_dt = S(0.0);
  }
  // [eqn 21] instantaneous mortality rate; productivity_area = net / area_leaf.
  using std::exp;
  if (mortality_finite) {
    const S productivity_area = net / area_leaf;
    r.mortality_dt = p.d_I + p.a_dG1 * exp(-p.a_dG2 * productivity_area);
  } else {
    r.mortality_dt = S(0.0);
  }
  return r;
}

// [eqn] Yokozawa leaf-area density q(z,H) = 2 eta (1 - u^eta) u^eta / z,
// u = z/H. Mirrors CanopyShape::q exactly (eta is a fixed double; the gradient
// flows through the active u and z). The deep-crown crown integral weights the
// per-depth assimilation by this.
template <typename S>
S ff16_canopy_q(double eta, S u, S z) {
  using std::pow;
  const S u_eta = pow(u, eta);
  return 2.0 * eta * (1.0 - u_eta) * u_eta / z;
}

// Single-plant height TRAJECTORY: integrate dheight/dt from h0 over n fixed RK4
// steps to age t_end, in a fixed crown-top light light_E (#472 scope B). This is
// the bridge from instantaneous-rate gradients to time-integrated emergent
// outputs: the whole trajectory is scalar-templated, so reverse/forward AD gives
// d(height at age t_end)/d(trait) -- a calibration gradient through the growth
// ODE. A fixed step schedule is exactly the frozen-schedule formulation
// end-to-end AD needs (the adaptive stepper is a pass-1 discovery, replayed).
template <typename S>
S ff16_grow_height(const FF16ProdPars<S>& p, S h0, S light_E,
                   double t_end, int n_steps) {
  S h = h0;
  const double dt = t_end / n_steps;
  for (int i = 0; i < n_steps; ++i) {
    // Materialise each stage as S: with XAD expression templates `h + c*k`
    // is an expression type, not S, which would break template deduction of
    // ff16_height_dt_crown_top's scalar.
    const S k1 = ff16_height_dt_crown_top(p, h, light_E);
    const S h2 = h + S(0.5 * dt) * k1;
    const S k2 = ff16_height_dt_crown_top(p, h2, light_E);
    const S h3 = h + S(0.5 * dt) * k2;
    const S k3 = ff16_height_dt_crown_top(p, h3, light_E);
    const S h4 = h + S(dt) * k3;
    const S k4 = ff16_height_dt_crown_top(p, h4, light_E);
    h = h + S(dt / 6.0) * (k1 + S(2.0) * k2 + S(2.0) * k3 + k4);
  }
  return h;
}

}  // namespace plant

#endif
