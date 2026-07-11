// -*-c++-*-
#ifndef PLANT_PLANT_FF16_STRATEGY_H_
#define PLANT_PLANT_FF16_STRATEGY_H_

#include <plant/strategy.h>
#include <plant/models/ff16_environment.h>
#include <plant/qag.h>
#include <plant/canopy_shape.h>
#include <plant/models/ff16_production_kernel.h>

namespace plant {

// Biological (user-settable) parameters for the FF16 strategy. Held as a value
// member `pars` on FF16_Strategy and exposed to R as a nested RcppR6 list class
// (so R access is `s$pars$lma`). Derived/precomputed quantities (eta_c,
// height_0, canopy_shape, ...) are NOT here -- they are outputs of
// prepare_strategy() and stay as plain members on the strategy.
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

// The double parameters cross the R boundary; the template default keeps this
// name a concrete type for RcppR6 and every existing caller.
using FF16_Pars = FF16_Pars_<double>;

template <class S = double>
class FF16_Strategy_: public Strategy<FF16_Environment> {
public:
  using value_type = S;
  typedef std::shared_ptr<FF16_Strategy_> ptr;
  FF16_Strategy_();

  // Fixed integer slots for the hot ODE rate path, used instead of
  // state_index.at("...") / aux_index.at("...") string-map lookups (those map
  // lookups showed up in profiling, see #466). These MUST stay in sync with
  // the order of state_names() and aux_names() below: *_INDEX is the position
  // of that name in state_names(), *_AUX_INDEX the position in aux_names(). If
  // you add/reorder a name there, update these constants (refresh_indices()
  // still validates the named maps used by the R-facing paths).
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
    if (collect_all_auxiliary) {
      ret.push_back("area_sapwood");
    }
    return ret;
  }

  // Translate generic methods to FF16 strategy leaf area methods

  S competition_effect(S height) const {
    return area_leaf(height);
  }

  void refresh_indices();


  // FF16 Methods  ----------------------------------------------

  // [eqn 2] area_leaf (inverse of [eqn 3])
  // Inline (header) so it can inline into the hot competition/assimilation
  // paths that reach it from templated Individual<FF16> code (no LTO build).
  S area_leaf(S height) const {
    // Single source: scalar-templated kernel (#472 scope B). Bit-identical to
    // the previous std::pow(height/a_l1, 1/a_l2); validated by the FF16
    // reference-comparison test.
    return ff16_area_leaf(pars.a_l1, pars.a_l2, height);
  }

  // [eqn 1] mass_leaf (inverse of [eqn 2])
  S mass_leaf(S area_leaf) const;

  // [eqn 4] area and mass of sapwood
  S area_sapwood(S area_leaf) const;
  S mass_sapwood(S area_sapwood, S height) const;

  // [eqn 5] area and mass of bark
  S area_bark(S area_leaf) const;
  S mass_bark (S area_bark, S height) const;

  S area_stem(S area_bark, S area_sapwood,
                            S area_heartwood) const;
  S diameter_stem(S area_stem) const;

  // [eqn 7] Mass of (fine) roots
  S mass_root(S area_leaf) const;

  // [eqn 8] Total Mass
  S mass_live(S mass_leaf, S mass_bark,
                   S mass_sapwood, S mass_root) const;

  S mass_total(S mass_leaf, S mass_bark, S mass_sapwood,
                    S mass_heartwood, S mass_root) const;

  // Above-ground mass = leaf + all stem components (bark + sapwood +
  // heartwood); excludes roots.
  S mass_above_ground(S mass_leaf, S mass_bark,
                           S mass_sapwood, S mass_heartwood) const;

  void compute_rates(const FF16_Environment& environment,
                Internals_<S>& vars);

  // Inline (header): called per state-set / ODE-state update from templated
  // Individual<FF16> code, so inlining avoids a cross-TU call (no LTO build)
  // and lets the now-inline area_leaf fold in.
  void update_dependent_aux(const int index, Internals_<S>& vars) {
    if (index == HEIGHT_INDEX) {
      S height = vars.state(HEIGHT_INDEX);
      vars.set_aux(COMPETITION_EFFECT_AUX_INDEX, area_leaf(height));
      vars.set_aux(HEIGHT_INVERSE_AUX_INDEX, 1.0 / height);
    }
  }

  // * Mass production
  // [eqn 12] Gross annual CO2 assimilation. Thin dispatcher: forwards to the
  // shading-model implementation bound once in prepare_strategy(), so the
  // choice of crown model costs a single predicted indirect call per
  // derivative evaluation and never a string comparison.
  S assimilation(const FF16_Environment& environment, S height,
                      S area_leaf, S height_inverse) {
    return (this->*assimilation_fn)(environment, height, area_leaf,
                                    height_inverse);
  }
  // Shading-model implementations of assimilation (selected in
  // prepare_strategy via assimilation_fn).
  //  - deep crown: integrate photosynthesis over crown depth (Yokozawa q).
  //  - crown top:  single evaluation of the light at the crown centre. Used by
  //                both crown-centre and PPA; they differ only in how the patch
  //                light profile is built (smooth vs stepped), which this read
  //                picks up transparently through the environment.
  S assimilation_deep_crown(const FF16_Environment& environment,
                                 S height, S area_leaf,
                                 S height_inverse);
  //  - average light: integrate the light over crown depth to a leaf-area-
  //                weighted mean, then a single photosynthesis evaluation.
  S assimilation_average_light(const FF16_Environment& environment,
                                    S height, S area_leaf,
                                    S height_inverse);
  S assimilation_crown_top(const FF16_Environment& environment,
                                S height, S area_leaf,
                                S height_inverse);

  typedef S (FF16_Strategy_::*assimilation_fn_t)(const FF16_Environment&,
                                                     S, S, S);
  // Bound once in prepare_strategy(); defaults to the deep-crown integral.
  assimilation_fn_t assimilation_fn = &FF16_Strategy_::assimilation_deep_crown;

  // [Appendix S6] Per-leaf photosynthetic rate.
  S assimilation_leaf(S x) const;

  // [eqn 13] Total maintenance respiration
  S respiration(S mass_leaf, S mass_sapwood,
                     S mass_bark, S mass_root) const;

  S respiration_leaf(S mass) const;
  S respiration_bark(S mass) const;
  S respiration_sapwood(S mass) const;
  S respiration_root(S mass) const;

  // [eqn 14] Total turnover
  S turnover(S mass_leaf, S mass_bark,
                  S mass_sapwood, S mass_root) const;
  S turnover_leaf(S mass) const;
  S turnover_bark(S mass) const;
  S turnover_sapwood(S mass) const;
  S turnover_root(S mass) const;

  // [eqn 15] Net production
  S net_mass_production_dt_A(S assimilation, S respiration,
                                  S turnover) const;

  virtual S net_mass_production_dt(const FF16_Environment& environment,
                                S height, S area_leaf_);
  S net_mass_production_dt(const FF16_Environment& environment,
                                S height, S area_leaf_,
                                S height_inverse);
  // Worker overload that also reports the sapwood intermediates so callers
  // (compute_rates) can reuse them instead of recomputing. Bit-identical: the
  // out-refs receive exactly the values the body already computed.
  S net_mass_production_dt(const FF16_Environment& environment,
                                S height, S area_leaf_,
                                S height_inverse,
                                S& area_sapwood_, S& mass_sapwood_);
  // Strategy-agnostic entry point used by Individual<FF16> (#266): reads the
  // height state and the cached aux slots itself, so the generic Individual
  // does not need to know FF16's state/aux layout.
  S net_mass_production_dt(const FF16_Environment& environment,
                                const Internals_<S>& vars) {
    return net_mass_production_dt(environment, vars.state(HEIGHT_INDEX),
                                  vars.aux(COMPETITION_EFFECT_AUX_INDEX),
                                  vars.aux(HEIGHT_INVERSE_AUX_INDEX));
  }

  // [eqn 16] Fraction of whole plan growth that is leaf
  virtual S fraction_allocation_reproduction(S height) const;
  S fraction_allocation_growth(S height) const;
  // [eqn 17] Rate of offspring production
  S fecundity_dt(S net_mass_production_dt,
                      S fraction_allocation_reproduction) const;

  // [eqn 18] Fraction of mass growth that is leaves
  S darea_leaf_dmass_live(S area_leaf) const;

  // change in height per change in leaf area
  S dheight_darea_leaf(S area_leaf) const;
  // Mass of leaf needed for new unit area leaf, d m_s / d a_l
  S dmass_leaf_darea_leaf(S area_leaf) const;
  // Mass of stem needed for new unit area leaf, d m_s / d a_l
  S dmass_sapwood_darea_leaf(S area_leaf) const;
  // Overloads taking a precomputed pow(area_leaf, a_l2): compute_rates needs
  // this term for both the height rate and the live-mass partition, so it
  // evaluates the (non-integer, libm) pow once and shares it -- see issue #361.
  S darea_leaf_dmass_live(S area_leaf,
                               S area_leaf_pow_a_l2) const;
  S dheight_darea_leaf(S area_leaf, S area_leaf_pow_a_l2) const;
  S dmass_sapwood_darea_leaf(S area_leaf,
                                  S area_leaf_pow_a_l2) const;
  // Mass of bark needed for new unit area leaf, d m_b / d a_l
  S dmass_bark_darea_leaf(S area_leaf) const;
  // Mass of root needed for new unit area leaf, d m_r / d a_l
  S dmass_root_darea_leaf(S area_leaf) const;
  // Growth rate of basal diameter_stem per unit stem area
  S ddiameter_stem_darea_stem(S area_stem) const;
  // Growth rate of components per unit time:
  S area_leaf_dt(S area_leaf_dt) const;
  S area_sapwood_dt(S area_leaf_dt) const;
  S area_heartwood_dt(S area_leaf) const;
  S area_bark_dt(S area_leaf_dt) const;
  S area_stem_dt(S area_leaf, S area_leaf_dt) const;
  S diameter_stem_dt(S area_stem, S area_stem_dt) const;
  S mass_root_dt(S area_leaf,
                       S area_leaf_dt) const;
  S mass_live_dt(S fraction_allocation_reproduction,
                       S net_mass_production_dt) const;
  S mass_total_dt(S fraction_allocation_reproduction,
                        S net_mass_production_dt,
                        S mass_heartwood_dt) const;
  S mass_above_ground_dt(S area_leaf,
                               S fraction_allocation_reproduction,
                               S net_mass_production_dt,
                               S mass_heartwood_dt,
                               S area_leaf_dt) const;

  S mass_heartwood_dt(S mass_sapwood) const;

  S mass_live_given_height(S height) const;
  S height_given_mass_leaf(S mass_leaf_) const;


  S mortality_dt(S productivity_area, S cumulative_mortality) const;
  S mortality_growth_independent_dt()const ;
  S mortality_growth_dependent_dt(S productivity_area) const;
  // [eqn 20] Survival of seedlings during establishment
  S establishment_probability(const FF16_Environment& environment);

  // * Competitive environment
  // [eqn 11] total projected leaf area above height above height `z` for given plant
  // Inline (header) so the per-node hot competition path called from
  // Individual<FF16>::compute_competition can inline these tiny helpers
  // instead of paying a cross-TU call each iteration (no LTO build).
  S compute_competition(double z, S height) const {
    return compute_competition(z, area_leaf(height), 1.0 / height);
  }
  S compute_competition(double z, S area_leaf_,
                             S height_inverse) const {
    return compute_competition_by_ratio(z * height_inverse, area_leaf_);
  }
  S compute_competition_by_ratio(S z_over_height,
                                      S area_leaf_) const {
    return pars.k_I * area_leaf_ * canopy_shape.leaf_area_above(z_over_height);
  }
  // Strategy-agnostic entry point used by Individual<FF16> (#266): reads the
  // cached competition_effect (= area_leaf) and height_inverse aux slots
  // itself. Inline (header) to keep the per-node hot competition path free of
  // a cross-TU call (no LTO build).
  S compute_competition(double z, const Internals_<S>& vars) const {
    return compute_competition(z, vars.aux(COMPETITION_EFFECT_AUX_INDEX),
                               vars.aux(HEIGHT_INVERSE_AUX_INDEX));
  }

  // [      ] Inverse of Q: height above which fraction 'x' of leaf found
  S Qp(S x, S height) const;

  // The aim is to find a plant height that gives the correct seed mass.
  S height_seed(void) const;

  // Set constants within FF16_Strategy
  void prepare_strategy();

  // Birth height of a (germinated) seed. Strategy-agnostic accessor used by
  // the templated Individual; here height_0 is derived in prepare_strategy().
  S initial_height() const { return height_0; }

  // Biological (user-settable) parameters; see FF16_Pars above.
  FF16_Pars_<S> pars;

  // Derived / precomputed in prepare_strategy() (NOT user-set) -------------
  // Crown shape factor, precomputed from pars.eta
  S eta_c     = NA_REAL; // [dimensionless]
  CanopyShape canopy_shape;
  // Height and leaf area of a (germinated) seed
  S height_0  = NA_REAL;
  S height_0_inverse = NA_REAL;
  S area_leaf_0;

  // For integrating functions with using Gauss-Kronrod quadrature
  quadrature::QK function_integrator;
};

// S=double is the resident numerics that cross the R boundary; the alias keeps
// this name a concrete type for RcppR6 and every existing caller.
using FF16_Strategy = FF16_Strategy_<double>;

FF16_Strategy::ptr make_strategy_ptr(FF16_Strategy s);

}

#endif
