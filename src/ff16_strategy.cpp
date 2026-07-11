#include <plant/models/ff16_strategy.h>
#include <plant/models/ff16_production_kernel.h>
#include <plant/individual.h>
#include <plant/patch.h>
#include <odelia/ode_solver.hpp>
#include <XAD/XAD.hpp>
#include <algorithm>
#include <cmath>

namespace plant {

template <class S>
FF16_Strategy_<S>::FF16_Strategy_() {
  collect_all_auxiliary = false;
  // build the string state/aux name to index map
  refresh_indices();
  name = "FF16";
}

template <class S>
void FF16_Strategy_<S>::refresh_indices () {
    // Create and fill the name to state index maps
  state_index = std::map<std::string,int>();
  aux_index   = std::map<std::string,int>();
  std::vector<std::string> aux_names_vec = aux_names();
  std::vector<std::string> state_names_vec = state_names();
  for (size_t i = 0; i < state_names_vec.size(); i++) {
    state_index[state_names_vec[i]] = i;
  }
  for (size_t i = 0; i < aux_names_vec.size(); i++) {
    aux_index[aux_names_vec[i]] = i;
  }
}

// area_leaf() is defined inline in ff16_strategy.h (hot path).

// [eqn 1] mass_leaf (inverse of [eqn 2])
template <class S>
S FF16_Strategy_<S>::mass_leaf(S area_leaf) const {
  return area_leaf * pars.lma;
}

// [eqn 4] area and mass of sapwood
template <class S>
S FF16_Strategy_<S>::area_sapwood(S area_leaf) const {
  return area_leaf * pars.theta;
}

template <class S>
S FF16_Strategy_<S>::mass_sapwood(S area_sapwood, S height) const {
  return area_sapwood * height * eta_c * pars.rho;
}

// [eqn 5] area and mass of bark
template <class S>
S FF16_Strategy_<S>::area_bark(S area_leaf) const {
  return pars.a_b1 * area_leaf * pars.theta;
}

template <class S>
S FF16_Strategy_<S>::mass_bark(S area_bark, S height) const {
  return area_bark * height * eta_c * pars.rho;
}

template <class S>
S FF16_Strategy_<S>::area_stem(S area_bark, S area_sapwood,
                            S area_heartwood) const {
  return area_bark + area_sapwood + area_heartwood;
}

template <class S>
S FF16_Strategy_<S>::diameter_stem(S area_stem) const {
  using std::sqrt;  // ADL picks up XAD's sqrt for the active scalar
  return sqrt(4 * area_stem / M_PI);
}

// [eqn 7] Mass of (fine) roots
template <class S>
S FF16_Strategy_<S>::mass_root(S area_leaf) const {
  return pars.a_r1 * area_leaf;
}

// [eqn 8] Total mass
template <class S>
S FF16_Strategy_<S>::mass_live(S mass_leaf, S mass_bark,
                           S mass_sapwood, S mass_root) const {
  return mass_leaf + mass_sapwood + mass_bark + mass_root;
}

template <class S>
S FF16_Strategy_<S>::mass_total(S mass_leaf, S mass_bark,
                            S mass_sapwood, S mass_heartwood,
                            S mass_root) const {
  return mass_leaf + mass_bark + mass_sapwood +  mass_heartwood + mass_root;
}

template <class S>
S FF16_Strategy_<S>::mass_above_ground(S mass_leaf, S mass_bark,
                            S mass_sapwood, S mass_heartwood) const {
  return mass_leaf + mass_bark + mass_sapwood + mass_heartwood;
}

// update_dependent_aux() is defined inline in ff16_strategy.h (hot path).

// one-shot update of the scm variables
// i.e. setting rates of ode vars from the state and updating aux vars
template <class S>
void FF16_Strategy_<S>::compute_rates(const FF16_Environment& environment,  Internals_<S>& vars) {

  S height = vars.state(HEIGHT_INDEX);
  S area_leaf_ = vars.aux(COMPETITION_EFFECT_AUX_INDEX);
  S height_inverse = vars.aux(HEIGHT_INVERSE_AUX_INDEX);

  // Reuse the sapwood intermediates the worker already computes (for
  // respiration/turnover) rather than recomputing them below; bit-identical.
  S area_sapwood_, mass_sapwood_;
  const S net_mass_production_dt_ =
    net_mass_production_dt(environment, height, area_leaf_, height_inverse,
                           area_sapwood_, mass_sapwood_);

  // store the aux sate
  vars.set_aux(NET_MASS_PRODUCTION_DT_AUX_INDEX, net_mass_production_dt_);

  if (net_mass_production_dt_ > 0) {

    const S fraction_allocation_reproduction_ = fraction_allocation_reproduction(height);
    // dheight_darea_leaf and the sapwood/bark terms in darea_leaf_dmass_live all
    // share pow(area_leaf, pars.a_l2); evaluate this libm pow once and reuse it for
    // both rates rather than paying it twice per node per step (issue #361).
    const S area_leaf_pow_a_l2 = pow(area_leaf_, pars.a_l2);
    const S darea_leaf_dmass_live_ = darea_leaf_dmass_live(area_leaf_, area_leaf_pow_a_l2);
    const S fraction_allocation_growth_ = fraction_allocation_growth(height);
    const S area_leaf_dt = net_mass_production_dt_ * fraction_allocation_growth_ * darea_leaf_dmass_live_;

    vars.set_rate(HEIGHT_INDEX, dheight_darea_leaf(area_leaf_, area_leaf_pow_a_l2) * area_leaf_dt);
    vars.set_rate(FECUNDITY_INDEX,
      fecundity_dt(net_mass_production_dt_, fraction_allocation_reproduction_));

    vars.set_rate(AREA_HEARTWOOD_INDEX, area_heartwood_dt(area_leaf_));
    vars.set_rate(MASS_HEARTWOOD_INDEX, mass_heartwood_dt(mass_sapwood_));

    if (collect_all_auxiliary) {
      vars.set_aux(AREA_SAPWOOD_AUX_INDEX, area_sapwood_);
    }
  } else {
    vars.set_rate(HEIGHT_INDEX, 0.0);
    vars.set_rate(FECUNDITY_INDEX, 0.0);
    vars.set_rate(AREA_HEARTWOOD_INDEX, 0.0);
    vars.set_rate(MASS_HEARTWOOD_INDEX, 0.0);
  }
  // [eqn 21] - Instantaneous mortality rate
  vars.set_rate(MORTALITY_INDEX,
      mortality_dt(net_mass_production_dt_ / area_leaf_, vars.state(MORTALITY_INDEX)));
}

// [eqn 12] Gross annual CO2 assimilation -- deep-crown model.
// Integrate photosynthesis over crown depth: for a given height in the crown,
// take photosynthesis at that depth multiplied by the amount of leaf there.
template <class S>
S FF16_Strategy_<S>::assimilation_deep_crown(const FF16_Environment& environment,
                                              S height,
                                              S area_leaf,
                                              S height_inverse) {

  S A = 0.0;

  // Define an anonymous function to integrate.
  // Keep the lambda's own closure type (do not wrap in std::function) so the
  // templated QK::integrate inlines the integrand at each quadrature point
  // instead of making a type-erased indirect call.
  // Hoist the light-spline upper bound (canopy top) out of the integrand: it
  // is invariant across the quadrature, so fetch it once and pass it into the
  // capped get_environment_at_height() overload rather than re-reading
  // spline.max() at every quadrature point.
  // The crown integral runs in double: the light profile is frozen (invasion)
  // and canopy_shape/QK are the precomputed double machinery. area_leaf carries
  // the active derivative; the integrand's own trait/height sensitivity is
  // recovered once the crown quadrature runs at the active scalar.
  const double canopy_top = environment.max_environment_height();
  auto f = [&](double z) -> double {
    return ad_value(assimilation_leaf(environment.get_environment_at_height(z, canopy_top))) *
      canopy_shape.q(ad_value(z * height_inverse), z);
  };

  // Integrate over crown depth using using Gauss-Kronrod quadrature.
  // The number of points used in the integration is determined by the control parameter
  // function_integration_rule. Rules defined in qk_rules.cpp
  A = function_integrator.integrate(f, 0.0, ad_value(height));

  return area_leaf * A;
}

// [eqn 12] Gross annual CO2 assimilation -- mean-light model.
// Integrate the *light* over crown depth, weighted by the leaf-area density q
// (which integrates to one over the crown), to get the leaf-area-weighted mean
// light the crown experiences, then take a single photosynthesis evaluation of
// that mean. This sits between deep-crown (which integrates the concave
// photosynthetic rate itself) and crown-centre (a single point evaluation): it
// captures the mean light exactly but ignores the curvature of photosynthesis
// across the within-crown light distribution.
template <class S>
S FF16_Strategy_<S>::assimilation_average_light(const FF16_Environment& environment,
                                                 S height,
                                                 S area_leaf,
                                                 S height_inverse) {
  const double canopy_top = environment.max_environment_height();
  auto f = [&](double z) -> double {
    return environment.get_environment_at_height(z, canopy_top) *
      canopy_shape.q(ad_value(z * height_inverse), z);
  };
  const double mean_light = function_integrator.integrate(f, 0.0, ad_value(height));
  return area_leaf * assimilation_leaf(mean_light);
}

// [eqn 12] Gross annual CO2 assimilation -- crown-top model (crown-centre and PPA).
// Leaf area is treated as a thin layer at the crown centre, so a single
// evaluation of the light environment there replaces the crown-depth integral.
// The crown-centre and PPA models share this code: under crown-centre the environment
// returns the smooth light profile, under PPA it returns the stepped (layered)
// profile, so the only difference between them lives in the environment build.
// (height_inverse is unused but kept for a common dispatch signature.)
template <class S>
S FF16_Strategy_<S>::assimilation_crown_top(const FF16_Environment& environment,
                                             S height,
                                             S area_leaf,
                                             S /* height_inverse */) {
  const double E = environment.get_environment_at_height(ad_value(height * eta_c));
  return area_leaf * assimilation_leaf(E);
}

// Photosynthetic rate per leaf area
// `x` is openness, ranging from 0 to 1.
template <class S>
S FF16_Strategy_<S>::assimilation_leaf(S x) const {
  // Single source: scalar-templated kernel (#472 scope B, Milestone A).
  return ff16_assimilation_leaf(pars.a_p1, pars.a_p2, x);
}

// [eqn 13] Total maintenance respiration
// NOTE: In contrast with Falster ref model, we do not normalise by pars.a_y*pars.a_bio.
template <class S>
S FF16_Strategy_<S>::respiration(S mass_leaf, S mass_sapwood,
                             S mass_bark, S mass_root) const {
  // Single source: scalar-templated kernel (#472 scope B, Milestone A).
  return ff16_respiration(mass_leaf, mass_sapwood, mass_bark, mass_root,
                          pars.r_l, pars.r_s, pars.r_b, pars.r_r);
}

template <class S>
S FF16_Strategy_<S>::respiration_leaf(S mass) const {
  return pars.r_l * mass;
}

template <class S>
S FF16_Strategy_<S>::respiration_bark(S mass) const {
  return pars.r_b * mass;
}

template <class S>
S FF16_Strategy_<S>::respiration_sapwood(S mass) const {
  return pars.r_s * mass;
}

template <class S>
S FF16_Strategy_<S>::respiration_root(S mass) const {
  return pars.r_r * mass;
}

// [eqn 14] Total turnover
template <class S>
S FF16_Strategy_<S>::turnover(S mass_leaf, S mass_bark,
                          S mass_sapwood, S mass_root) const {
   // Single source: scalar-templated kernel (#472 scope B, Milestone A).
   return ff16_turnover(mass_leaf, mass_bark, mass_sapwood, mass_root,
                        pars.k_l, pars.k_b, pars.k_s, pars.k_r);
}

template <class S>
S FF16_Strategy_<S>::turnover_leaf(S mass) const {
  return pars.k_l * mass;
}

template <class S>
S FF16_Strategy_<S>::turnover_bark(S mass) const {
  return pars.k_b * mass;
}

template <class S>
S FF16_Strategy_<S>::turnover_sapwood(S mass) const {
  return pars.k_s * mass;
}

template <class S>
S FF16_Strategy_<S>::turnover_root(S mass) const {
  return pars.k_r * mass;
}

// [eqn 15] Net production
//
// NOTE: Translation of variable names from the Falster 2011.  Everything
// before the minus sign is SCM's N, our `net_mass_production_dt` is SCM's P.
template <class S>
S FF16_Strategy_<S>::net_mass_production_dt_A(S assimilation, S respiration,
                                S turnover) const {
  // Single source: scalar-templated kernel (#472 scope B, Milestone A).
  return ff16_net_production_A(pars.a_bio, pars.a_y, assimilation, respiration, turnover);
}

// One shot calculation of net_mass_production_dt
// Used by establishment_probability() and compute_rates().
template <class S>
S FF16_Strategy_<S>::net_mass_production_dt(const FF16_Environment& environment,
                                S height, S area_leaf_) {
  return net_mass_production_dt(environment, height, area_leaf_,
                                1.0 / height);
}

template <class S>
S FF16_Strategy_<S>::net_mass_production_dt(const FF16_Environment& environment,
                                S height, S area_leaf_,
                                S height_inverse) {
  S area_sapwood_, mass_sapwood_;
  return net_mass_production_dt(environment, height, area_leaf_, height_inverse,
                                area_sapwood_, mass_sapwood_);
}

template <class S>
S FF16_Strategy_<S>::net_mass_production_dt(const FF16_Environment& environment,
                                S height, S area_leaf_,
                                S height_inverse,
                                S& area_sapwood_, S& mass_sapwood_) {
  const S mass_leaf_    = mass_leaf(area_leaf_);
  area_sapwood_ = area_sapwood(area_leaf_);
  mass_sapwood_ = mass_sapwood(area_sapwood_, height);
  const S area_bark_    = area_bark(area_leaf_);
  const S mass_bark_    = mass_bark(area_bark_, height);
  const S mass_root_    = mass_root(area_leaf_);
  const S assimilation_ =
    assimilation(environment, height, area_leaf_, height_inverse);
  const S respiration_ =
    respiration(mass_leaf_, mass_sapwood_, mass_bark_, mass_root_);
  const S turnover_ =
    turnover(mass_leaf_, mass_bark_, mass_sapwood_, mass_root_);
  return net_mass_production_dt_A(assimilation_, respiration_, turnover_);
}

// [eqn 16] Fraction of production allocated to reproduction
template <class S>
S FF16_Strategy_<S>::fraction_allocation_reproduction(S height) const {
  return pars.a_f1 / (1.0 + exp(pars.a_f2 * (1.0 - height / pars.hmat)));
}

// Fraction of production allocated to growth
template <class S>
S FF16_Strategy_<S>::fraction_allocation_growth(S height) const {
  return 1.0 - fraction_allocation_reproduction(height);
}

// [eqn 17] Rate of offspring production
template <class S>
S FF16_Strategy_<S>::fecundity_dt(S net_mass_production_dt,
                               S fraction_allocation_reproduction) const {
  return net_mass_production_dt * fraction_allocation_reproduction /
    (pars.omega + pars.a_f3);
}

template <class S>
S FF16_Strategy_<S>::darea_leaf_dmass_live(S area_leaf) const {
  return darea_leaf_dmass_live(area_leaf, pow(area_leaf, pars.a_l2));
}

template <class S>
S FF16_Strategy_<S>::darea_leaf_dmass_live(S area_leaf,
                                            S area_leaf_pow_a_l2) const {
  // dmass_bark_darea_leaf(area_leaf) == pars.a_b1 * dmass_sapwood_darea_leaf(area_leaf),
  // so compute the shared pow(area_leaf, pars.a_l2) term once rather than twice.
  const S dmass_sapwood_darea_leaf_ =
    dmass_sapwood_darea_leaf(area_leaf, area_leaf_pow_a_l2);
  return 1.0/(  dmass_leaf_darea_leaf(area_leaf)
              + dmass_sapwood_darea_leaf_
              + pars.a_b1 * dmass_sapwood_darea_leaf_
              + dmass_root_darea_leaf(area_leaf));
}

template <class S>
S FF16_Strategy_<S>::dheight_darea_leaf(S area_leaf) const {
  return pars.a_l1 * pars.a_l2 * pow(area_leaf, pars.a_l2 - 1);
}

template <class S>
S FF16_Strategy_<S>::dheight_darea_leaf(S area_leaf,
                                         S area_leaf_pow_a_l2) const {
  // pow(area_leaf, pars.a_l2 - 1) == pow(area_leaf, pars.a_l2) / area_leaf, so reuse the
  // shared power instead of a second libm pow (a reciprocal-multiply reorder,
  // so not bit-identical to the single-argument form).
  return pars.a_l1 * pars.a_l2 * area_leaf_pow_a_l2 / area_leaf;
}

// Mass of leaf needed for new unit area leaf, d m_s / d a_l
template <class S>
S FF16_Strategy_<S>::dmass_leaf_darea_leaf(S /* area_leaf */) const {
  return pars.lma;
}

// Mass of stem needed for new unit area leaf, d m_s / d a_l
template <class S>
S FF16_Strategy_<S>::dmass_sapwood_darea_leaf(S area_leaf) const {
  return dmass_sapwood_darea_leaf(area_leaf, pow(area_leaf, pars.a_l2));
}

template <class S>
S FF16_Strategy_<S>::dmass_sapwood_darea_leaf(S /* area_leaf */,
                                               S area_leaf_pow_a_l2) const {
  return pars.rho * eta_c * pars.a_l1 * pars.theta * (pars.a_l2 + 1.0) * area_leaf_pow_a_l2;
}

// Mass of bark needed for new unit area leaf, d m_b / d a_l
template <class S>
S FF16_Strategy_<S>::dmass_bark_darea_leaf(S area_leaf) const {
  return pars.a_b1 * dmass_sapwood_darea_leaf(area_leaf);
}

// Mass of root needed for new unit area leaf, d m_r / d a_l
template <class S>
S FF16_Strategy_<S>::dmass_root_darea_leaf(S /* area_leaf */) const {
  return pars.a_r1;
}

// Growth rate of basal diameter_stem per unit time
template <class S>
S FF16_Strategy_<S>::ddiameter_stem_darea_stem(S area_stem) const {
  return pow(M_PI * area_stem, -0.5);
}

// Growth rate of sapwood area at base per unit time
template <class S>
S FF16_Strategy_<S>::area_sapwood_dt(S area_leaf_dt) const {
  return area_leaf_dt * pars.theta;
}

// Note, unlike others, heartwood growth does not depend on leaf area growth, but
// rather existing sapwood
template <class S>
S FF16_Strategy_<S>::area_heartwood_dt(S area_leaf) const {
  return pars.k_s * area_sapwood(area_leaf);
}

// Growth rate of bark area at base per unit time
template <class S>
S FF16_Strategy_<S>::area_bark_dt(S area_leaf_dt) const {
  return pars.a_b1 * area_leaf_dt * pars.theta;
}

// Growth rate of stem basal area per unit time
template <class S>
S FF16_Strategy_<S>::area_stem_dt(S area_leaf,
                               S area_leaf_dt) const {
  return area_sapwood_dt(area_leaf_dt) +
    area_bark_dt(area_leaf_dt) +
    area_heartwood_dt(area_leaf);
}

// Growth rate of basal diameter_stem per unit time
template <class S>
S FF16_Strategy_<S>::diameter_stem_dt(S area_stem, S area_stem_dt) const {
  return ddiameter_stem_darea_stem(area_stem) * area_stem_dt;
}

// Growth rate of root mass per unit time
template <class S>
S FF16_Strategy_<S>::mass_root_dt(S area_leaf,
                               S area_leaf_dt) const {
  return area_leaf_dt * dmass_root_darea_leaf(area_leaf);
}

template <class S>
S FF16_Strategy_<S>::mass_live_dt(S fraction_allocation_reproduction,
                               S net_mass_production_dt) const {
  return (1 - fraction_allocation_reproduction) * net_mass_production_dt;
}

template <class S>
S FF16_Strategy_<S>::mass_total_dt(S fraction_allocation_reproduction,
                                     S net_mass_production_dt,
                                     S mass_heartwood_dt) const {
  return mass_live_dt(fraction_allocation_reproduction, net_mass_production_dt) +
    mass_heartwood_dt;
}

// TODO(#480): Do we not track root mass change?
template <class S>
S FF16_Strategy_<S>::mass_above_ground_dt(S area_leaf,
                                       S fraction_allocation_reproduction,
                                       S net_mass_production_dt,
                                       S mass_heartwood_dt,
                                       S area_leaf_dt) const {
  const S mass_root_dt =
    area_leaf_dt * dmass_root_darea_leaf(area_leaf);
  return mass_total_dt(fraction_allocation_reproduction, net_mass_production_dt,
                        mass_heartwood_dt) - mass_root_dt;
}

template <class S>
S FF16_Strategy_<S>::mass_heartwood_dt(S mass_sapwood) const {
  return turnover_sapwood(mass_sapwood);
}


template <class S>
S FF16_Strategy_<S>::mass_live_given_height(S height) const {
  S area_leaf_ = area_leaf(height);
  return mass_leaf(area_leaf_) +
         mass_bark(area_bark(area_leaf_), height) +
         mass_sapwood(area_sapwood(area_leaf_), height) +
         mass_root(area_leaf_);
}

template <class S>
S FF16_Strategy_<S>::height_given_mass_leaf(S mass_leaf) const {
  return pars.a_l1 * pow(mass_leaf / pars.lma, pars.a_l2);
}

template <class S>
S FF16_Strategy_<S>::mortality_dt(S productivity_area,
                              S cumulative_mortality) const {

  // NOTE: When plants are extremely inviable, the rate of change in
  // mortality can be Inf, because net production is negative, leaf
  // area is small and so we get exp(big number).  However, most of
  // the time that happens we should get infinite mortality variable
  // levels and the rate of change won't matter.  It is possible that
  // we will need to trim this to some large finite value, but for
  // now, just checking that the actual mortality rate is finite.
  if (util::is_finite(ad_value(cumulative_mortality))) {
    return
      mortality_growth_independent_dt() +
      mortality_growth_dependent_dt(productivity_area);
 } else {
    // If mortality probability is 1 (latency = Inf) then the rate
    // calculations break.  Setting them to zero gives the correct
    // behaviour.
    return 0.0;
  }
}

template <class S>
S FF16_Strategy_<S>::mortality_growth_independent_dt() const {
  return pars.d_I;
}

template <class S>
S FF16_Strategy_<S>::mortality_growth_dependent_dt(S productivity_area) const {
  return pars.a_dG1 * exp(-pars.a_dG2 * productivity_area);
}

// [eqn 20] Survival of seedlings during establishment
template <class S>
S FF16_Strategy_<S>::establishment_probability(const FF16_Environment& environment) {

  S decay_over_time = exp(-pars.recruitment_decay * environment.time);

  const S net_mass_production_dt_ =
    net_mass_production_dt(environment, height_0, area_leaf_0,
                           height_0_inverse);
  if (net_mass_production_dt_ > 0) {
    const S tmp = pars.a_d0 * area_leaf_0 / net_mass_production_dt_;
    return 1.0 / (tmp * tmp + 1.0) * decay_over_time;
  } else {
    return 0.0;
  }
}

// compute_competition() overloads and compute_competition_by_ratio() are
// defined inline in ff16_strategy.h (per-node hot path).

// (inverse of [eqn 10]; return the height above which fraction 'x' of
// the leaf mass would be found).
template <class S>
S FF16_Strategy_<S>::Qp(S x, S height) const { // x in [0,1], unchecked.
  return canopy_shape.Qp(ad_value(x), ad_value(height));
}

// The aim is to find a plant height that gives the correct seed mass.
template <class S>
S FF16_Strategy_<S>::height_seed(void) const {

  // Note, these are not entirely correct bounds. Ideally we would use height
  // given *total* mass, not leaf mass, but that is difficult to calculate.
  // Using "height given leaf mass" will expand upper bound, but that's ok
  // most of time. Only issue is that could break with obscure parameter
  // values for LMA or height-leaf area scaling. Could instead use some
  // absolute maximum height for new seedling, e.g. 1m?
  const double
    h0 = ad_value(height_given_mass_leaf(std::numeric_limits<double>::min())),
    h1 = ad_value(height_given_mass_leaf(pars.omega));

  const double tol = control.offspring_production_tol;
  const size_t max_iterations = control.offspring_production_iterations;

  // Converge the seed height in double -- the bracketing and convergence test
  // are genuine double control-flow. mass_live_given_height stays at the active
  // scalar so its value is exact; ad_value takes the passive residual the
  // bisection compares.
  auto target = [&] (double x) mutable -> double {
    return ad_value(mass_live_given_height(x) - pars.omega);
  };

  const double h_root = util::uniroot(target, h0, h1, tol, max_iterations);

  if constexpr (std::is_same_v<S, double>) {
    return h_root;
  } else {
    // Reattach the trait derivative the double solve discards, by the implicit
    // function theorem. The root h*(theta) solves g(h*, theta) =
    // mass_live_given_height(h*) - omega = 0, so dh*/dtheta = -(dg/dtheta) /
    // (dg/dh). Evaluated at the converged root g's value is ~0, so the corrected
    // scalar keeps h_root's value while carrying the seed's derivative. dg/dh is
    // a passive slope from a central difference in double.
    const double dh = std::max(std::abs(h_root), 1.0) * 1e-6;
    const double dg_dh =
      (ad_value(mass_live_given_height(h_root + dh)) -
       ad_value(mass_live_given_height(h_root - dh))) / (2.0 * dh);
    const S g = mass_live_given_height(h_root) - pars.omega;
    const S g_deriv = g - ad_value(g);
    return S(h_root) - g_deriv / dg_dh;
  }
}

template <class S>
void FF16_Strategy_<S>::prepare_strategy() {

  // Set up the function_integrator
  function_integrator = quadrature::QK(
      // Gauss-Kronrod quadrature integeration rule (see qkrules)
      control.function_integration_rule);

  // Resolve the crown shading model once (string -> enum), then bind both hot
  // paths to it: canopy_shape handles competition (leaf_area_above), and
  // assimilation_fn selects the matching assimilation implementation. After
  // this, neither path compares the model string per call.
  const ShadingModel shading_model =
    shading_model_from_string(control.shading_model, ShadingModel::DeepCrown);
  // canopy_shape also selects the competition contribution: smooth Q for every
  // model except flat-top-box, which casts a step (see CanopyShape). It is the
  // precomputed double shape; its eta derivative is the resident self-shading
  // term added when the environment is templated on the scalar. eta_c above
  // keeps the active eta for the mass cascade.
  canopy_shape.initialise(ad_value(pars.eta), shading_model);
  switch (shading_model) {
  case ShadingModel::DeepCrown:
    assimilation_fn = &FF16_Strategy_<S>::assimilation_deep_crown;
    break;
  case ShadingModel::MeanLight:
    assimilation_fn = &FF16_Strategy_<S>::assimilation_average_light;
    break;
  case ShadingModel::CrownCentre:
  case ShadingModel::FlatTopBox:
  case ShadingModel::FlatTopSoftBox:
  case ShadingModel::PPA:
    // All evaluate assimilation at the crown centre. They differ in the light
    // profile they read: crown-centre from the smooth profile, flat-top-box / -soft-
    // box from a profile built with (hard / smoothed) box competition, PPA from
    // a stepped profile.
    assimilation_fn = &FF16_Strategy_<S>::assimilation_crown_top;
    break;
  }

  // NOTE: this pre-computes something to save a very small amount of time
  eta_c = 1 - 2/(1 + pars.eta) + 1/(1 + 2*pars.eta);
  // NOTE: Also pre-computing, though less trivial
  height_0 = height_seed();
  height_0_inverse = 1.0 / height_0;
  area_leaf_0 = area_leaf(height_0);

  if (is_variable_birth_rate) {
    extrinsic_drivers.set_variable("birth_rate", birth_rate_x, birth_rate_y);
  } else {
    extrinsic_drivers.set_constant("birth_rate", birth_rate_y[0]);
  }
}

FF16_Strategy::ptr make_strategy_ptr(FF16_Strategy s) {
  s.prepare_strategy();
  return std::make_shared<FF16_Strategy>(s);
}

// The resident numerics cross the R boundary at double. The forward-mode active
// scalar is instantiated here too, so plant.so carries the differentiable FF16
// physiology: prepare_strategy -> height_seed and the mass cascade differentiate
// exactly w.r.t. traits; the crown integral and canopy shape carry the physiology
// derivative with the double-frozen light, their own integral/shape derivatives
// completed once the crown quadrature and environment run at the active scalar.
template class FF16_Strategy_<double>;
template class FF16_Strategy_<xad::fwd<double>::active_type>;

// Reverse-mode active scalar via odelia's Solver alias (plant names no XAD
// reverse primitive); emits FF16's out-of-line virtuals for the active vtable.
using ad_reverse =
    odelia::ode::Solver<Patch<FF16_Strategy_<double>, FF16_Environment>>::active_scalar;
template class FF16_Strategy_<ad_reverse>;
}
