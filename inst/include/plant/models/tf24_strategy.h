// Built from  inst/include/plant/models/ff16_strategy.h on Mon Feb 12 09:52:27 2024 using the scaffolder, from the strategy:  FF16
// -*-c++-*-
#ifndef PLANT_PLANT_TF24_STRATEGY_H_
#define PLANT_PLANT_TF24_STRATEGY_H_

#include <plant/strategy.h>
#include <plant/models/tf24_environment.h>
#include <plant/qag.h>
#include <plant/leaf_model.h>
#include <plant/canopy_shape.h> // ShadingModel
#include <odelia/ode_util.hpp>  // odelia::util::to_passive (Leaf boundary is double)

namespace plant {

// Biological (user-settable) parameters for the TF24 strategy. Held as a value
// member `pars` on TF24_Strategy and exposed to R as a nested RcppR6 list class
// (access as `s$pars$lma`). Templated on the scalar S so a trait derivative
// flows when a field is seeded; S = double is the production path (the `TF24_Pars`
// alias below). Only parameters that were previously exposed to R live here;
// derived quantities (eta_c, height_0, ...), the embedded Leaf model, solver
// tolerances and hard-coded hydraulic-root constants stay as plain members on
// the strategy.
template <class S = double>
struct TF24_Pars_ {
  // * Core traits
  S lma       = 0.1978791;  // Leaf mass per area [kg / m2]
  S rho       = 608.0;      // Wood density [kg/m3]
  S hmat      = 16.5958691; // Height at maturation [m]
  S omega     = 3.8e-5;     // Seed mass [kg]
  // * Individual allometry
  S eta       = 12.0;       // Canopy shape parameter [dimensionless]
  S theta     = 1.0/4669;   // Sapwood area per leaf area [dimensionless]
  S a_l1      = 5.44;       // height with 1m2 leaf [m]
  S a_l2      = 0.306;      // scaling of height with leaf area
  S a_r1      = 0.07;       // Root mass per leaf area [kg / m]
  S a_b1      = 0.17;       // Ratio of bark area : sapwood area
  // * Production
  S r_s    = 4012.0 / 608.0; // Sapwood respiration per stem mass
  S r_b    = 2.0 * (4012.0 / 608.0); // Bark respiration (assumed 2 x sapwood)
  S r_r    = 217.0;          // Root respiration per mass
  S r_l    = 39.27 / 0.1978791; // Leaf dark respiration per leaf mass
  S a_y    = 0.7;            // Carbon conversion parameter
  S a_bio  = 2.45e-2;        // CO2 -> dry mass [kg / mol]
  S k_l    = 0.4565855;      // Leaf turnover [/yr]
  S k_b    = 0.2;            // Bark turnover [/yr]
  S k_s    = 0.2;            // Sapwood turnover [/yr]
  S k_r    = 1.0;            // Root turnover [/yr]
  S a_p1   = 151.177775377968;   // LRC hyperbola [mol CO2 / yr / m2]
  S a_p2   = 0.204716166503633;  // LRC hyperbola shape
  // * Seed production
  S a_f3   = 3.0 *  3.8e-5;  // Accessory cost of reproduction [kg/seed]
  S a_f1   = 1.0;            // Maximum allocation to reproduction
  S a_f2   = 50;             // Size range across which individuals mature
  // * Mortality parameters
  S S_D    = 0.25;           // Probability of survival during dispersal
  S a_d0   = 0.1;            // Parameter for seedling survival
  S d_I    = 0.01;           // Baseline intrinsic mortality [/yr]
  S a_dG1  = 5.5;            // Baseline growth-related mortality [/yr]
  S a_dG2  = 20.0;           // Risk coefficient for dry mass production
  // * Light capture
  S k_I = 0.5;
  // * Leaf hydraulic / photosynthesis traits (default Eucalyptus saligna)
  S vcmax_25 = 96;
  S p_50 = 1.85;
  S K_s = 1;
  S c = log(log(1-0.5)/log(1-0.88))/(log(p_50) - log(5.16));
  S b = p_50 / pow(-log(1 - 50.0 / 100.0), 1 / c);
  S psi_crit = b*pow(log(1/0.05),1/c); // derived from b and c
  S beta1 = 20000;
  S beta2 = 1.5;
  S g1_TF24 = 7.5;
  S jmax_25 = vcmax_25*1.64;
  S a = 0.30; // effective quantum yield of electron transport
  S curv_fact_elec_trans = 0.7;
  S curv_fact_colim = 0.99;
  S var_sapwood_volume_cost = 1;
  // nitrogen allocation traits (parameterised from Austraits 4.1.0)
  S nmass_l = 13e-3; // kg N kg^-1 mass
  S nmass_s = 1.98e-3; // kg N kg^-1 mass
  S nmass_b = 3.40e-3; // kg N kg^-1 mass
  S nmass_r = 3.35e-3; // kg N kg^-1 mass
  S dmass_dN = 0; // change in mass per change in kg kg^-1 N
  // shape exponent for the Q() root-fraction-with-depth profile
  S root_depth_shape_eta = 0.2;
  // Germination
  S recruitment_decay = 0.0;
};

using TF24_Pars = TF24_Pars_<double>;

// The TF24 strategy. Templated on the scalar S carried by its physiology; S =
// double is the production path (the `TF24_Strategy` alias below). Method bodies
// live in tf24_strategy.cpp with explicit instantiations for double and the
// active reverse scalar: TF24's density-transport dg/dh is delivered by
// supplied_derivative (the leaf optimiser, envelope theorem), not the nested
// forward-over-reverse path, so its instantiation set is closed (no rebind).
// The embedded Leaf hydraulic sub-model stays double throughout; the S mass
// cascade crosses to it via to_passive at the boundary. Base members of
// Strategy<E> are reached through this-> because the base is dependent on S.
template <class S = double>
class TF24_Strategy_: public Strategy<TF24_Environment_<S>> {
public:
  using environment_type = TF24_Environment_<S>;
  using value_type = S;
  typedef std::shared_ptr<TF24_Strategy_<S>> ptr;
  TF24_Strategy_();

  double compute_average_light_environment(double z, double height,
                                           const environment_type &environment);

  // calculate the amount of water transpired relativised by leaf area index.

  S evapotranspiration_dt(S area_leaf_, int soil_layer);


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
      "net_mass_production_dt",
      "root_mass",
      "opt_psi_stem",
      "opt_root_psi",
      "transpiration",
      "E_up_",
      "profit",
      "stom_cond_CO2",
      "assimilation"
    });
    // add the associated computation to compute_rates and compute there
    if (this->collect_all_auxiliary) {
      ret.push_back("area_sapwood");
    }
    return ret;
  }

  // Translate generic methods to TF24 strategy leaf area methods

  S competition_effect(S height) const {
    return area_leaf(height);
  }

  void refresh_indices();


  // TF24 Methods  ----------------------------------------------

  // [eqn 2] area_leaf (inverse of [eqn 3])
  S area_leaf(S height) const;

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

  void compute_rates(const environment_type& environment,
                Internals_<S>& vars);

  void update_dependent_aux(const int index, Internals_<S>& vars);

  // * Mass production
  // [eqn 12] Gross annual CO2 assimilation
  S assimilation(const environment_type& environment, S height,
                 S area_leaf);
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

  virtual S net_mass_production_dt(const environment_type& environment,
                                S height, S area_leaf_,
                                S height_inverse);

  // Resolve the leaf operating point on the already-set-up `leaf` (i.e. after
  // leaf.set_physiology(...)). Base TF24 optimises the root-collar psi via
  // golden-section search; the TF24f variant overrides this to make the optimum
  // chase a tracked ODE state (#525). Called per crown light point from
  // net_mass_production_dt, so it must be virtual to dispatch to the override
  // when net_mass_production_dt is reused unchanged by the subclass.
  virtual void solve_leaf();
  // Strategy-agnostic entry point used by Individual<TF24> (#266): reads the
  // height state and the cached aux slots itself, so the generic Individual
  // does not need to know TF24's state/aux layout.
  S net_mass_production_dt(const environment_type& environment,
                           const Internals_<S>& vars) {
    return net_mass_production_dt(environment, vars.state(HEIGHT_INDEX),
                                  vars.aux(aux_idx_competition_effect),
                                  vars.aux(aux_idx_height_inverse));
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
  S establishment_probability(const environment_type& environment);

  // * Competitive environment
  // [eqn 11] total projected leaf area above height above height `z` for given plant
  S compute_competition(double z, S height) const;
  // Optimised overload called from Individual<TF24>::compute_competition with the
  // cached competition_effect (= area_leaf(height)) and height_inverse (= 1/height)
  // aux values, matching the shared individual.h interface (no recompute per call).
  S compute_competition(double z, S area_leaf_,
                        S height_inverse) const;
  // Strategy-agnostic entry point used by Individual<TF24> (#266): reads the
  // cached competition_effect and height_inverse aux slots itself.
  S compute_competition(double z, const Internals_<S>& vars) const {
    return compute_competition(z, vars.aux(aux_idx_competition_effect),
                               vars.aux(aux_idx_height_inverse));
  }

  // [eqn  9] Probability density of leaf area at height `z`. Off the S rate
  // path: only used to weight the double Leaf's crown aggregation, so double.
  double q(double z, double height) const;
  // [eqn 10] Fraction of leaf area above height `z`. Double: used for the root
  // depth distribution feeding the double Leaf (compute_competition inlines its
  // own S form for the light field).
  double Q(double z, double height, double eta_x) const;
  // [      ] Inverse of Q: height above which fraction 'x' of leaf found
  double Qp(double x, double height) const;

  // The aim is to find a plant height that gives the correct seed mass.
  double height_seed(void) const;

  // Set constants within TF24_Strategy
  void prepare_strategy();

  // Birth height of a (germinated) seed. Strategy-agnostic accessor used by
  // the templated Individual; here height_0 is derived in prepare_strategy().
  S initial_height() const { return height_0; }

  // Crown shading model, resolved once from control.shading_model in
  // prepare_strategy(). TF24 supports deep-crown, mean-light (its default)
  // and crown-centre; PPA is not available for TF24.
  ShadingModel shading_model_ = ShadingModel::MeanLight;

  // Biological (user-settable) parameters; see TF24_Pars above.
  TF24_Pars_<S> pars;

  // Derived / precomputed in prepare_strategy() (NOT user-set) -------------
  S eta_c     = NA_REAL; // crown shape factor, precomputed from pars.eta
  // Height and leaf area of a (germinated) seed
  S height_0  = NA_REAL;
  S area_leaf_0;

  // Embedded leaf hydraulic/photosynthesis sub-model, built in prepare_strategy().
  // Always double: the Leaf is never templated on S (design 4.3); its parameter
  // sensitivity reaches the tape via supplied_derivative, not by templating.
  Leaf leaf;

  // Hydraulic root parameters (not currently exposed to R; see review #9)
  double root_c = 2.680147;
  double root_b = 3.898245;
  double root_psi_crit = root_b*std::pow(log(1/0.05),1/root_c); // derived from root_b and root_c

  // Solver tolerances and other constants not currently exposed to R
  double newton_tol_abs = 0.001;
  double GSS_tol_abs = 1e-3;
  double vulnerability_curve_ncontrol = 100;
  double ci_abs_tol = 1e-6;
  double ci_niter = 1000;
  double beta_R_H = 3.4e2;
  double beta_R_V = 9.4e3;

  // Cached aux/state indices, resolved once in refresh_indices(), so the hot
  // compute_rates path does not do a std::map<string,int>::at (string compare)
  // lookup per ODE derivs evaluation per individual (profile hot spot).
  int aux_idx_competition_effect = -1;
  int aux_idx_height_inverse = -1;
  int aux_idx_net_mass_production_dt = -1;
  int aux_idx_root_mass = -1;
  int aux_idx_opt_psi_stem = -1;
  int aux_idx_opt_root_psi = -1;
  int aux_idx_transpiration = -1;
  int aux_idx_E_up = -1;
  int aux_idx_profit = -1;
  int aux_idx_stom_cond_CO2 = -1;
  int aux_idx_area_sapwood = -1;       // only present when collect_all_auxiliary
  int state_idx_area_heartwood = -1;
  int state_idx_mass_heartwood = -1;

  // For integrating functions with using Gauss-Kronrod quadrature
  quadrature::QK function_integrator;

  // Reusable per-layer root-mass buffer, refilled (not reallocated) each
  // net_mass_production_dt call to avoid a heap allocation per derivs eval.
  std::vector<double> mass_root_prop_;
};

using TF24_Strategy = TF24_Strategy_<double>;

}

#endif
