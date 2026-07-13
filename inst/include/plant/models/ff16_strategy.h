// -*-c++-*-
#ifndef PLANT_PLANT_FF16_STRATEGY_H_
#define PLANT_PLANT_FF16_STRATEGY_H_

#include <plant/strategy.h>
#include <plant/models/ff16_environment.h>
#include <plant/qag.h>
#include <plant/canopy_shape.h>
#include <cmath>
#include <limits> // std::numeric_limits (height_seed bounds)

namespace plant {

// Biological (user-settable) parameters for the FF16 strategy. Held as a value
// member `pars` on FF16_Strategy and exposed to R as a nested RcppR6 list class
// (so R access is `s$pars$lma`). Templated on the scalar S so a trait derivative
// flows when a field is seeded; S = double is the production path (the `FF16_Pars`
// alias below). Derived/precomputed quantities (eta_c, height_0, canopy_shape,
// ...) are NOT here -- they are outputs of prepare_strategy() and stay as plain
// members on the strategy.
template <class S = double>
struct FF16_Pars_ {
  // * Core traits
  S lma       = 0.1978791;  // Leaf mass per area [kg / m2]
  S rho       = 608.0;      // Wood density [kg/m3]
  S hmat      = 16.5958691; // Height at maturation [m]
  S omega     = 3.8e-5;     // Seed mass [kg]
  // * Individual allometry
  // Canopy shape parameter
  S eta       = 12.0; // [dimensionless]
  // Sapwood area per leaf area
  // Ratio sapwood area area to leaf area
  S theta     = 1.0/4669; // [dimensionless]
  // Height - leaf mass scaling
  S a_l1        = 5.44; // height with 1m2 leaf [m]
  S a_l2        = 0.306; // dimensionless scaling of height with leaf area
  // Root mass per leaf area
  S a_r1        = 0.07;  //[kg / m]
  // Ratio of bark area : sapwood area
  S a_b1         = 0.17; // [dimensionless]

  // * Production
  // Ratio of leaf dark respiration to leaf mass [mol CO2 / yr  / kg]
  // =  [mol CO2 / m2 / yr]  |  (39.27 = 2100 * 0.00187)  | narea * photosynthesis_per_nitrogen
  //    / [kg(leaf) / m2 ]   |    / (0.1978791)           | lma
  // Hard coded in value of lma here so that this value doesn't change
  // if that trait changes above.
  S r_l    = 39.27 / 0.1978791;
  // Root respiration per mass [mol CO2 / yr / kg]
  S r_r    = 217.0;
  // Sapwood respiration per stem mass  [mol CO2 / yr / kg]
  // = respiration per volume [mol CO2 / m3 / yr]
  // /  wood density [kg/m3]
  S r_s    = 4012.0 / 608.0;
  // Bark respiration per stem mass
  // assumed to be twice rate of sapwood
  // (NOTE that there is a re-parametrisation here relative to the paper
  // -- r_b is defined (new) as 2*r_s, whereas the paper assumes a
  // fixed multiplication by 2)
  S r_b    = 2.0 * r_s;
  // Carbon conversion parameter
  S a_y    = 0.7;
  // Constant converting assimilated CO2 to dry mass [kg / mol]
  // (12E-3 / 0.49)
  S a_bio  = 2.45e-2;
  // Leaf turnover [/yr]
  S k_l    =  0.4565855;
  // Bark turnover [/yr]
  S k_b    = 0.2;
  // Sapwood turnover [/yr]
  S k_s           = 0.2;
  // Root turnover [/yr]
  S k_r    = 1.0;
  // Parameters of the hyperbola for annual LRC
  S a_p1   = 151.177775377968; // [mol CO2 / yr / m2]
  S a_p2   = 0.204716166503633; // [dimensionless]

  // * Seed production
  // Accessory cost of reproduction
  S a_f3  = 3.0 *  3.8e-5; // [kg per seed]
  // Maximum allocation to reproduction
  S a_f1   = 1.0; //[dimensionless]
  // Size range across which individuals mature
  S a_f2   = 50; // [dimensionless]

  // * Mortality parameters
  // Probability of survival during dispersal
  S S_D   = 0.25; // [dimensionless]
  // Parameter for seedling survival
  S a_d0    = 0.1; //[kg / yr / m2]
  // Baseline for intrinsic mortality
  S d_I    = 0.01; // [ / yr]
  // Baseline rate for growth-related mortality
  S a_dG1    = 5.5; // [ / yr]
  // Risk coefficient for dry mass production (per area)
  S a_dG2    = 20.0;// [yr m2 / kg ]

  // Germination
  S recruitment_decay = 0.0;

  // * Light capture parameters
  S k_I = 0.5;
};

using FF16_Pars = FF16_Pars_<double>;

// The single list of low-level FF16 parameters the gradient is taken with
// respect to (§8.1): field_ptrs() and field_names() on FF16_Strategy_ both
// expand it, so a parameter can never appear in one and not the other, and no
// field count is hand-typed.
#define FF16_AD_FIELDS(X)                                              \
  X(lma) X(rho) X(hmat) X(omega) X(eta) X(theta) X(a_l1) X(a_l2)       \
  X(a_r1) X(a_b1) X(r_l) X(r_r) X(r_s) X(r_b) X(a_y) X(a_bio) X(k_l)   \
  X(k_b) X(k_s) X(k_r) X(a_p1) X(a_p2) X(a_f3) X(a_f1) X(a_f2) X(S_D)  \
  X(a_d0) X(d_I) X(a_dG1) X(a_dG2) X(recruitment_decay) X(k_I)

// The FF16 (Falster et al 2016) strategy. Templated on the scalar S carried by
// its physiology; S = double is the production path (the `FF16_Strategy` alias
// below). All method bodies are defined here (inline template members) so the
// active (AD) instantiation can see them. Base members of Strategy<E> are
// reached through this-> because the base is dependent on S.
template <class S = double>
class FF16_Strategy_ : public Strategy<FF16_Environment_<S>> {
public:
  using environment_type = FF16_Environment_<S>;
  using value_type = S;
  typedef std::shared_ptr<FF16_Strategy_<S>> ptr;

  FF16_Strategy_() {
    this->collect_all_auxiliary = false;
    // build the string state/aux name to index map
    refresh_indices();
    this->name = "FF16";
  }

  // Fixed integer slots for the hot ODE rate path, used instead of
  // state_index.at("...") / aux_index.at("...") string-map lookups (those map
  // lookups showed up in profiling, see #466). These MUST stay in sync with
  // the order of state_names() and aux_names() below.
  static constexpr int AREA_HEARTWOOD_INDEX = 3;
  static constexpr int MASS_HEARTWOOD_INDEX = 4;
  static constexpr int COMPETITION_EFFECT_AUX_INDEX = 0;
  static constexpr int HEIGHT_INVERSE_AUX_INDEX = 1;
  static constexpr int NET_MASS_PRODUCTION_DT_AUX_INDEX = 2;
  static constexpr int AREA_SAPWOOD_AUX_INDEX = 3;

  // Overrides ----------------------------------------------

  // update this when the length of state_names changes
  static size_t state_size () { return 5; }
  // update this when the length of aux_names changes
  size_t aux_size () { return aux_names().size(); }

  static std::vector<std::string> state_names() {
    return  std::vector<std::string>({
      "height",
      "mortality",
      "fecundity",
      "area_heartwood",
      "mass_heartwood"
      });
  }

  std::vector<std::string> aux_names() {
    std::vector<std::string> ret({
      "competition_effect",
      "height_inverse",
      "net_mass_production_dt"
    });
    // add the associated computation to compute_rates and compute there
    if (this->collect_all_auxiliary) {
      ret.push_back("area_sapwood");
    }
    return ret;
  }

  // Translate generic methods to FF16 strategy leaf area methods

  S competition_effect(S height) const {
    return area_leaf(height);
  }

  void refresh_indices () {
    // Create and fill the name to state index maps
    this->state_index = std::map<std::string,int>();
    this->aux_index   = std::map<std::string,int>();
    std::vector<std::string> aux_names_vec = aux_names();
    std::vector<std::string> state_names_vec = state_names();
    for (size_t i = 0; i < state_names_vec.size(); i++) {
      this->state_index[state_names_vec[i]] = i;
    }
    for (size_t i = 0; i < aux_names_vec.size(); i++) {
      this->aux_index[aux_names_vec[i]] = i;
    }
  }


  // FF16 Methods  ----------------------------------------------

  // [eqn 2] area_leaf (inverse of [eqn 3])
  S area_leaf(S height) const {
    return pow(height / pars.a_l1, 1.0 / pars.a_l2);
  }

  // [eqn 1] mass_leaf (inverse of [eqn 2])
  S mass_leaf(S area_leaf) const {
    return area_leaf * pars.lma;
  }

  // [eqn 4] area and mass of sapwood
  S area_sapwood(S area_leaf) const {
    return area_leaf * pars.theta;
  }

  S mass_sapwood(S area_sapwood, S height) const {
    return area_sapwood * height * eta_c * pars.rho;
  }

  // [eqn 5] area and mass of bark
  S area_bark(S area_leaf) const {
    return pars.a_b1 * area_leaf * pars.theta;
  }

  S mass_bark (S area_bark, S height) const {
    return area_bark * height * eta_c * pars.rho;
  }

  S area_stem(S area_bark, S area_sapwood,
              S area_heartwood) const {
    return area_bark + area_sapwood + area_heartwood;
  }

  S diameter_stem(S area_stem) const {
    using std::sqrt;
    return sqrt(4 * area_stem / M_PI);
  }

  // [eqn 7] Mass of (fine) roots
  S mass_root(S area_leaf) const {
    return pars.a_r1 * area_leaf;
  }

  // [eqn 8] Total Mass
  S mass_live(S mass_leaf, S mass_bark,
              S mass_sapwood, S mass_root) const {
    return mass_leaf + mass_sapwood + mass_bark + mass_root;
  }

  S mass_total(S mass_leaf, S mass_bark, S mass_sapwood,
               S mass_heartwood, S mass_root) const {
    return mass_leaf + mass_bark + mass_sapwood +  mass_heartwood + mass_root;
  }

  // Above-ground mass = leaf + all stem components (bark + sapwood +
  // heartwood); excludes roots.
  S mass_above_ground(S mass_leaf, S mass_bark,
                      S mass_sapwood, S mass_heartwood) const {
    return mass_leaf + mass_bark + mass_sapwood + mass_heartwood;
  }

  // Inline (header): called per state-set / ODE-state update from templated
  // Individual<FF16> code.
  void update_dependent_aux(const int index, Internals_<S>& vars) {
    if (index == HEIGHT_INDEX) {
      S height = vars.state(HEIGHT_INDEX);
      vars.set_aux(COMPETITION_EFFECT_AUX_INDEX, area_leaf(height));
      vars.set_aux(HEIGHT_INVERSE_AUX_INDEX, 1.0 / height);
    }
  }

  // one-shot update of the scm variables
  // i.e. setting rates of ode vars from the state and updating aux vars
  void compute_rates(const environment_type& environment, Internals_<S>& vars) {

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

      if (this->collect_all_auxiliary) {
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

  // * Mass production
  // [eqn 12] Gross annual CO2 assimilation. Thin dispatcher: forwards to the
  // shading-model implementation bound once in prepare_strategy(), so the
  // choice of crown model costs a single predicted indirect call per
  // derivative evaluation and never a string comparison.
  S assimilation(const environment_type& environment, S height,
                 S area_leaf, S height_inverse) {
    return (this->*assimilation_fn)(environment, height, area_leaf,
                                    height_inverse);
  }

  // [eqn 12] Gross annual CO2 assimilation -- deep-crown model.
  // Integrate photosynthesis over crown depth: for a given height in the crown,
  // take photosynthesis at that depth multiplied by the amount of leaf there.
  S assimilation_deep_crown(const environment_type& environment,
                            S height, S area_leaf, S height_inverse) {

    S A = 0.0;

    // Keep the lambda's own closure type (do not wrap in std::function) so the
    // templated QK::integrate inlines the integrand at each quadrature point.
    // Hoist the light-spline upper bound (canopy top) out of the integrand.
    const double canopy_top = environment.max_environment_height();
    auto f = [&](S z) -> S {
      // Collapse the XAD expression z * height_inverse to a concrete S so both
      // arguments of the templated q(Z, Z) deduce the same Z at an active scalar.
      const S z_over_height = z * height_inverse;
      return assimilation_leaf(environment.get_environment_at_height(z, static_cast<S>(canopy_top))) *
        canopy_shape.q(z_over_height, z);
    };

    // Integrate over crown depth using Gauss-Kronrod quadrature. A fixed rule,
    // so the active bound tapes exactly (Cluster 4).
    A = function_integrator.integrate(f, static_cast<S>(0.0), height);

    return area_leaf * A;
  }

  // [eqn 12] Gross annual CO2 assimilation -- mean-light model (TF24's default).
  // Integrate the *light* over crown depth, weighted by the leaf-area density q,
  // to get the leaf-area-weighted mean light the crown experiences, then take a
  // single photosynthesis evaluation of that mean. Re-bodied to S so the crown
  // derivative is not silently dropped when MeanLight is selected (plan §4.2 M3).
  S assimilation_average_light(const environment_type& environment,
                               S height, S area_leaf, S height_inverse) {
    const double canopy_top = environment.max_environment_height();
    auto f = [&](S z) -> S {
      const S z_over_height = z * height_inverse;
      return environment.get_environment_at_height(z, static_cast<S>(canopy_top)) *
        canopy_shape.q(z_over_height, z);
    };
    const S mean_light = function_integrator.integrate(f, static_cast<S>(0.0), height);
    return area_leaf * assimilation_leaf(mean_light);
  }

  // [eqn 12] Gross annual CO2 assimilation -- crown-top model (crown-centre and PPA).
  // (height_inverse is unused but kept for a common dispatch signature.)
  S assimilation_crown_top(const environment_type& environment,
                           S height, S area_leaf, S /* height_inverse */) {
    const S E = environment.get_environment_at_height(height * eta_c);
    return area_leaf * assimilation_leaf(E);
  }

  typedef S (FF16_Strategy_::*assimilation_fn_t)(const environment_type&,
                                                 S, S, S);
  // Bound once in prepare_strategy(); defaults to the deep-crown integral.
  assimilation_fn_t assimilation_fn = &FF16_Strategy_::assimilation_deep_crown;

  // [Appendix S6] Per-leaf photosynthetic rate. `x` is openness, 0 to 1.
  S assimilation_leaf(S x) const {
    return pars.a_p1 * x / (x + pars.a_p2);
  }

  // [eqn 13] Total maintenance respiration
  // NOTE: In contrast with Falster ref model, we do not normalise by pars.a_y*pars.a_bio.
  S respiration(S mass_leaf, S mass_sapwood,
                S mass_bark, S mass_root) const {
    return respiration_leaf(mass_leaf) +
           respiration_bark(mass_bark) +
           respiration_sapwood(mass_sapwood) +
           respiration_root(mass_root);
  }

  S respiration_leaf(S mass) const { return pars.r_l * mass; }
  S respiration_bark(S mass) const { return pars.r_b * mass; }
  S respiration_sapwood(S mass) const { return pars.r_s * mass; }
  S respiration_root(S mass) const { return pars.r_r * mass; }

  // [eqn 14] Total turnover
  S turnover(S mass_leaf, S mass_bark,
             S mass_sapwood, S mass_root) const {
     return turnover_leaf(mass_leaf) +
            turnover_bark(mass_bark) +
            turnover_sapwood(mass_sapwood) +
            turnover_root(mass_root);
  }
  S turnover_leaf(S mass) const { return pars.k_l * mass; }
  S turnover_bark(S mass) const { return pars.k_b * mass; }
  S turnover_sapwood(S mass) const { return pars.k_s * mass; }
  S turnover_root(S mass) const { return pars.k_r * mass; }

  // [eqn 15] Net production
  //
  // NOTE: Translation of variable names from the Falster 2011.  Everything
  // before the minus sign is SCM's N, our `net_mass_production_dt` is SCM's P.
  S net_mass_production_dt_A(S assimilation, S respiration,
                             S turnover) const {
    return pars.a_bio * pars.a_y * (assimilation - respiration) - turnover;
  }

  // One shot calculation of net_mass_production_dt.
  // Used by establishment_probability() and compute_rates().
  virtual S net_mass_production_dt(const environment_type& environment,
                                   S height, S area_leaf_) {
    return net_mass_production_dt(environment, height, area_leaf_,
                                  1.0 / height);
  }
  S net_mass_production_dt(const environment_type& environment,
                           S height, S area_leaf_,
                           S height_inverse) {
    S area_sapwood_, mass_sapwood_;
    return net_mass_production_dt(environment, height, area_leaf_, height_inverse,
                                  area_sapwood_, mass_sapwood_);
  }
  // Worker overload that also reports the sapwood intermediates so callers
  // (compute_rates) can reuse them instead of recomputing. Bit-identical: the
  // out-refs receive exactly the values the body already computed.
  S net_mass_production_dt(const environment_type& environment,
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
  // Strategy-agnostic entry point used by Individual<FF16> (#266): reads the
  // height state and the cached aux slots itself.
  S net_mass_production_dt(const environment_type& environment,
                           const Internals_<S>& vars) {
    return net_mass_production_dt(environment, vars.state(HEIGHT_INDEX),
                                  vars.aux(COMPETITION_EFFECT_AUX_INDEX),
                                  vars.aux(HEIGHT_INVERSE_AUX_INDEX));
  }

  // [eqn 16] Fraction of production allocated to reproduction
  virtual S fraction_allocation_reproduction(S height) const {
    return pars.a_f1 / (1.0 + exp(pars.a_f2 * (1.0 - height / pars.hmat)));
  }
  // Fraction of production allocated to growth
  S fraction_allocation_growth(S height) const {
    return 1.0 - fraction_allocation_reproduction(height);
  }
  // [eqn 17] Rate of offspring production
  S fecundity_dt(S net_mass_production_dt,
                 S fraction_allocation_reproduction) const {
    return net_mass_production_dt * fraction_allocation_reproduction /
      (pars.omega + pars.a_f3);
  }

  // [eqn 18] Fraction of mass growth that is leaves
  S darea_leaf_dmass_live(S area_leaf) const {
    return darea_leaf_dmass_live(area_leaf, pow(area_leaf, pars.a_l2));
  }
  S darea_leaf_dmass_live(S area_leaf,
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

  // change in height per change in leaf area
  S dheight_darea_leaf(S area_leaf) const {
    return pars.a_l1 * pars.a_l2 * pow(area_leaf, pars.a_l2 - 1);
  }
  S dheight_darea_leaf(S area_leaf,
                       S area_leaf_pow_a_l2) const {
    // pow(area_leaf, pars.a_l2 - 1) == pow(area_leaf, pars.a_l2) / area_leaf, so reuse the
    // shared power instead of a second libm pow (a reciprocal-multiply reorder,
    // so not bit-identical to the single-argument form).
    return pars.a_l1 * pars.a_l2 * area_leaf_pow_a_l2 / area_leaf;
  }

  // Mass of leaf needed for new unit area leaf, d m_s / d a_l
  S dmass_leaf_darea_leaf(S /* area_leaf */) const {
    return pars.lma;
  }

  // Mass of stem needed for new unit area leaf, d m_s / d a_l
  S dmass_sapwood_darea_leaf(S area_leaf) const {
    return dmass_sapwood_darea_leaf(area_leaf, pow(area_leaf, pars.a_l2));
  }
  S dmass_sapwood_darea_leaf(S /* area_leaf */,
                             S area_leaf_pow_a_l2) const {
    return pars.rho * eta_c * pars.a_l1 * pars.theta * (pars.a_l2 + 1.0) * area_leaf_pow_a_l2;
  }

  // Mass of bark needed for new unit area leaf, d m_b / d a_l
  S dmass_bark_darea_leaf(S area_leaf) const {
    return pars.a_b1 * dmass_sapwood_darea_leaf(area_leaf);
  }

  // Mass of root needed for new unit area leaf, d m_r / d a_l
  S dmass_root_darea_leaf(S /* area_leaf */) const {
    return pars.a_r1;
  }

  // Growth rate of basal diameter_stem per unit stem area
  S ddiameter_stem_darea_stem(S area_stem) const {
    return pow(M_PI * area_stem, -0.5);
  }

  // Growth rate of sapwood area at base per unit time
  S area_sapwood_dt(S area_leaf_dt) const {
    return area_leaf_dt * pars.theta;
  }

  // Note, unlike others, heartwood growth does not depend on leaf area growth, but
  // rather existing sapwood
  S area_heartwood_dt(S area_leaf) const {
    return pars.k_s * area_sapwood(area_leaf);
  }

  // Growth rate of bark area at base per unit time
  S area_bark_dt(S area_leaf_dt) const {
    return pars.a_b1 * area_leaf_dt * pars.theta;
  }

  // Growth rate of stem basal area per unit time
  S area_stem_dt(S area_leaf,
                 S area_leaf_dt) const {
    return area_sapwood_dt(area_leaf_dt) +
      area_bark_dt(area_leaf_dt) +
      area_heartwood_dt(area_leaf);
  }

  // Growth rate of basal diameter_stem per unit time
  S diameter_stem_dt(S area_stem, S area_stem_dt) const {
    return ddiameter_stem_darea_stem(area_stem) * area_stem_dt;
  }

  // Growth rate of root mass per unit time
  S mass_root_dt(S area_leaf,
                 S area_leaf_dt) const {
    return area_leaf_dt * dmass_root_darea_leaf(area_leaf);
  }

  S mass_live_dt(S fraction_allocation_reproduction,
                 S net_mass_production_dt) const {
    return (1 - fraction_allocation_reproduction) * net_mass_production_dt;
  }

  S mass_total_dt(S fraction_allocation_reproduction,
                  S net_mass_production_dt,
                  S mass_heartwood_dt) const {
    return mass_live_dt(fraction_allocation_reproduction, net_mass_production_dt) +
      mass_heartwood_dt;
  }

  // TODO(#480): Do we not track root mass change?
  S mass_above_ground_dt(S area_leaf,
                         S fraction_allocation_reproduction,
                         S net_mass_production_dt,
                         S mass_heartwood_dt,
                         S area_leaf_dt) const {
    const S mass_root_dt =
      area_leaf_dt * dmass_root_darea_leaf(area_leaf);
    return mass_total_dt(fraction_allocation_reproduction, net_mass_production_dt,
                          mass_heartwood_dt) - mass_root_dt;
  }

  S mass_heartwood_dt(S mass_sapwood) const {
    return turnover_sapwood(mass_sapwood);
  }

  S mass_live_given_height(S height) const {
    S area_leaf_ = area_leaf(height);
    return mass_leaf(area_leaf_) +
           mass_bark(area_bark(area_leaf_), height) +
           mass_sapwood(area_sapwood(area_leaf_), height) +
           mass_root(area_leaf_);
  }

  S height_given_mass_leaf(S mass_leaf) const {
    return pars.a_l1 * pow(mass_leaf / pars.lma, pars.a_l2);
  }

  S mortality_dt(S productivity_area, S cumulative_mortality) const {

    // NOTE: When plants are extremely inviable, the rate of change in
    // mortality can be Inf, because net production is negative, leaf
    // area is small and so we get exp(big number).  However, most of
    // the time that happens we should get infinite mortality variable
    // levels and the rate of change won't matter.  It is possible that
    // we will need to trim this to some large finite value, but for
    // now, just checking that the actual mortality rate is finite.
    if (util::is_finite(cumulative_mortality)) {
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

  S mortality_growth_independent_dt() const {
    return pars.d_I;
  }

  S mortality_growth_dependent_dt(S productivity_area) const {
    return pars.a_dG1 * exp(-pars.a_dG2 * productivity_area);
  }

  // [eqn 20] Survival of seedlings during establishment
  S establishment_probability(const environment_type& environment) {

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

  // * Competitive environment
  // [eqn 11] total projected leaf area above height above height `z` for given plant
  S compute_competition(double z, S height) const {
    return compute_competition(z, area_leaf(height), 1.0 / height);
  }
  S compute_competition(double z, S area_leaf_,
                        S height_inverse) const {
    return compute_competition_by_ratio(z * height_inverse, area_leaf_);
  }
  // Templated on the ratio scalar R: only the resident field-build hot path
  // (deferred) drives an active ratio into canopy_shape.leaf_area_above.
  template <typename R>
  S compute_competition_by_ratio(R z_over_height,
                                 S area_leaf_) const {
    return pars.k_I * area_leaf_ * canopy_shape.leaf_area_above(z_over_height);
  }
  // Strategy-agnostic entry point used by Individual<FF16> (#266).
  S compute_competition(double z, const Internals_<S>& vars) const {
    return compute_competition(z, vars.aux(COMPETITION_EFFECT_AUX_INDEX),
                               vars.aux(HEIGHT_INVERSE_AUX_INDEX));
  }

  // [      ] Inverse of Q: height above which fraction 'x' of leaf found.
  // R-facing/diagnostic, off the differentiated rate path -> double.
  double Qp(double x, double height) const {
    return canopy_shape.Qp(x, height);
  }

  // The aim is to find a plant height that gives the correct seed mass. Kept as a
  // double root-find (the supplied_derivative IFT seam is deferred, plan §7.2);
  // xad::value narrowings are no-ops at S = double.
  double height_seed(void) const {

    // Note, these are not entirely correct bounds. Ideally we would use height
    // given *total* mass, not leaf mass, but that is difficult to calculate.
    const double
      h0 = xad::value(height_given_mass_leaf(std::numeric_limits<double>::min())),
      h1 = xad::value(height_given_mass_leaf(pars.omega));

    const double tol = this->control.offspring_production_tol;
    const size_t max_iterations = this->control.offspring_production_iterations;

    auto target = [&] (double x) mutable -> double {
      return xad::value(mass_live_given_height(x)) - xad::value(pars.omega);
    };

    return util::uniroot(target, h0, h1, tol, max_iterations);
  }

  // Set constants within FF16_Strategy
  void prepare_strategy() {

    // Set up the function_integrator
    function_integrator = quadrature::QK(
        // Gauss-Kronrod quadrature integeration rule (see qkrules)
        this->control.function_integration_rule);

    // Resolve the crown shading model once (string -> enum), then bind both hot
    // paths to it: canopy_shape handles competition (leaf_area_above), and
    // assimilation_fn selects the matching assimilation implementation.
    const ShadingModel shading_model =
      shading_model_from_string(this->control.shading_model, ShadingModel::DeepCrown);
    // canopy_shape also selects the competition contribution: smooth Q for every
    // model except flat-top-box, which casts a step (see CanopyShape). eta is a
    // shape coefficient, not differentiated -> narrow to double.
    canopy_shape.initialise(xad::value(pars.eta), shading_model);
    switch (shading_model) {
    case ShadingModel::DeepCrown:
      assimilation_fn = &FF16_Strategy_::assimilation_deep_crown;
      break;
    case ShadingModel::MeanLight:
      assimilation_fn = &FF16_Strategy_::assimilation_average_light;
      break;
    case ShadingModel::CrownCentre:
    case ShadingModel::FlatTopBox:
    case ShadingModel::FlatTopSoftBox:
    case ShadingModel::PPA:
      // All evaluate assimilation at the crown centre. They differ in the light
      // profile they read (see FF16_Environment::compute_environment).
      assimilation_fn = &FF16_Strategy_::assimilation_crown_top;
      break;
    }

    // NOTE: this pre-computes something to save a very small amount of time
    eta_c = 1 - 2/(1 + pars.eta) + 1/(1 + 2*pars.eta);
    // NOTE: Also pre-computing, though less trivial
    height_0 = height_seed();
    height_0_inverse = 1.0 / height_0;
    area_leaf_0 = area_leaf(height_0);
    initial_height_ = lift_birth_height(height_0);

    if (this->is_variable_birth_rate) {
      this->extrinsic_drivers.set_variable("birth_rate", this->birth_rate_x, this->birth_rate_y);
    } else {
      this->extrinsic_drivers.set_constant("birth_rate", this->birth_rate_y[0]);
    }
  }

  // Birth height of a (germinated) seed. Strategy-agnostic accessor used by
  // the templated Individual; derived in prepare_strategy() and lifted to carry
  // its parameter derivative (see lift_birth_height).
  S initial_height() const { return initial_height_; }

  // Lift the double birth height to S carrying its parameter derivative. The
  // birth height solves g(h,theta) = mass_live_given_height(h) - omega = 0, whose
  // root is found in double (height_seed). mass_live_given_height is retapeable,
  // so one Newton step from that root evaluated at the active parameters yields
  // both the exact value (the double path is bit-identical) and the implicit-
  // function-theorem derivative dh*/dtheta = -(dg/dtheta)/(dg/dh) for every
  // parameter at once. dg/dh is a double central difference at the root.
  S lift_birth_height(double h_star) const {
    if constexpr (std::is_same_v<S, double>) {
      return h_star;
    } else {
      // dg/dh at the root, as a double constant (central difference; g is smooth
      // and near-linear in h here, so a 2-point rule is not the accuracy limit).
      const double eps = 1e-6 * (h_star + 1.0);
      auto gv = [&](double h) { return xad::value(mass_live_given_height(S(h))); };
      const double dgdh = (gv(h_star + eps) - gv(h_star - eps)) / (2.0 * eps);
      // Newton correction from the double root. Subtract its own value so the
      // birth height's *value* stays exactly h_star -- bit-identical, and no value
      // shift to perturb the gradient of parameters that height does not depend on
      // -- while its derivative is the IFT term dh*/dtheta = -(dg/dtheta)/(dg/dh).
      const S corr = (mass_live_given_height(S(h_star)) - pars.omega) / dgdh;
      return S(h_star) - corr + xad::value(corr);
    }
  }

  // Biological (user-settable) parameters; see FF16_Pars above.
  FF16_Pars_<S> pars;

  // Pointers to the low-level parameters the gradient is taken with respect to
  // (§8.1), in FF16_AD_FIELDS order; field_names() gives the matching column
  // labels. Both are generated from the one FF16_AD_FIELDS list, so they cannot
  // disagree in membership or order.
#define PLANT_AD_PTR(f) &pars.f,
#define PLANT_AD_NAME(f) #f,
  std::vector<S*> field_ptrs() { return { FF16_AD_FIELDS(PLANT_AD_PTR) }; }
  static std::vector<std::string> field_names() {
    return { FF16_AD_FIELDS(PLANT_AD_NAME) };
  }
#undef PLANT_AD_PTR
#undef PLANT_AD_NAME

  // Derived / precomputed in prepare_strategy() (NOT user-set) -------------
  // Crown shape factor, precomputed from pars.eta
  S eta_c     = NA_REAL; // [dimensionless]
  CanopyShape canopy_shape;
  // Height and leaf area of a (germinated) seed. height_0 stays double (the
  // root-find is double until the height_seed IFT seam lands).
  double height_0  = NA_REAL;
  double height_0_inverse = NA_REAL;
  S area_leaf_0;
  // Birth height lifted to S (value == height_0; carries the parameter
  // derivative on a gradient pass). Returned by initial_height().
  S initial_height_ = NA_REAL;

  // For integrating functions with using Gauss-Kronrod quadrature
  quadrature::QK function_integrator;
};

using FF16_Strategy = FF16_Strategy_<double>;

}

#endif
