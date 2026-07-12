// Built from  inst/include/plant/models/ff16_strategy.h on Mon Feb 12 09:52:27 2024 using the scaffolder, from the strategy:  FF16
// -*-c++-*-
#ifndef PLANT_PLANT_TF24_STRATEGY_H_
#define PLANT_PLANT_TF24_STRATEGY_H_

#include <plant/strategy.h>
#include <plant/models/tf24_environment.h>
#include <plant/qag.h>
#include <plant/leaf_model.h>
#include <plant/canopy_shape.h> // ShadingModel
#include <plant/ad_value.h>
#include <cmath>

namespace plant {

// Biological (user-settable) parameters for the TF24 strategy. Templated on the
// scalar S (double for the resident numerics, active for a trait gradient); the
// RcppR6 alias TF24_Pars = TF24_Pars_<double> keeps the R boundary at double.
// Derived quantities (eta_c, height_0, ...), the embedded Leaf model, solver
// tolerances and hard-coded hydraulic-root constants stay off `pars`.
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
  S r_b    = 2.0 * r_s;      // Bark respiration (assumed 2 x sapwood)
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

  // Differentiable-parameter handles in a FIXED, documented order -- this
  // declaration order IS the AD column-order contract, and rebind() carries every
  // field so the active twin's Leaf is built from the resident's values.
  std::vector<S*> field_ptrs() {
    return {&lma,   &rho,   &hmat,  &omega, &eta,   &theta, &a_l1,  &a_l2,
            &a_r1,  &a_b1,  &r_s,   &r_b,   &r_r,   &r_l,   &a_y,   &a_bio,
            &k_l,   &k_b,   &k_s,   &k_r,   &a_p1,  &a_p2,  &a_f3,  &a_f1,
            &a_f2,  &S_D,   &a_d0,  &d_I,   &a_dG1, &a_dG2, &k_I,
            &vcmax_25, &p_50, &K_s, &c, &b, &psi_crit, &beta1, &beta2,
            &g1_TF24, &jmax_25, &a, &curv_fact_elec_trans, &curv_fact_colim,
            &var_sapwood_volume_cost, &nmass_l, &nmass_s, &nmass_b, &nmass_r,
            &dmass_dN, &root_depth_shape_eta, &recruitment_decay};
  }
  static std::vector<std::string> field_names() {
    return {"lma",  "rho",  "hmat", "omega", "eta",  "theta", "a_l1", "a_l2",
            "a_r1", "a_b1", "r_s",  "r_b",   "r_r",  "r_l",   "a_y",  "a_bio",
            "k_l",  "k_b",  "k_s",  "k_r",   "a_p1", "a_p2",  "a_f3", "a_f1",
            "a_f2", "S_D",  "a_d0", "d_I",   "a_dG1","a_dG2", "k_I",
            "vcmax_25", "p_50", "K_s", "c", "b", "psi_crit", "beta1", "beta2",
            "g1_TF24", "jmax_25", "a", "curv_fact_elec_trans", "curv_fact_colim",
            "var_sapwood_volume_cost", "nmass_l", "nmass_s", "nmass_b", "nmass_r",
            "dmass_dN", "root_depth_shape_eta", "recruitment_decay"};
  }

  // Config-only copy onto another scalar (values via ad_value, no tape identity).
  template <class S2>
  TF24_Pars_<S2> rebind() const {
    TF24_Pars_<S2> out;
    std::vector<S*> src = const_cast<TF24_Pars_*>(this)->field_ptrs();
    std::vector<S2*> dst = out.field_ptrs();
    for (size_t i = 0; i < src.size(); ++i) {
      *dst[i] = S2(ad_value(*src[i]));
    }
    return out;
  }
};

// The double parameters cross the R boundary; the alias keeps this name a concrete
// type for RcppR6 and every existing caller.
using TF24_Pars = TF24_Pars_<double>;

template <class S = double>
class TF24_Strategy_: public Strategy<TF24_Environment> {
public:
  using value_type = S;
  typedef std::shared_ptr<TF24_Strategy_> ptr;
  TF24_Strategy_();

  // Leaf-area-weighted mean canopy openness at crown depth z (double; used only in
  // the off-tape leaf block). Reads the frozen canopy via ad_value.
  double compute_average_light_environment(double z, double height,
                                           const TF24_Environment_<S> &environment);

  // calculate the amount of water transpired relativised by leaf area index.
  double evapotranspiration_dt(double area_leaf_, int soil_layer);

  // Double replicas of the crown-shape kernels for the off-tape leaf solve (they
  // read ad_value(pars.eta) so they are bit-identical to q/Q/area_leaf at double).
  double q_double(double z, double height) const;
  double Q_double(double z, double height, double eta_x) const;
  double area_leaf_double(double height) const;


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
    if (collect_all_auxiliary) {
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

  S area_stem(S area_bark, S area_sapwood, S area_heartwood) const;
  S diameter_stem(S area_stem) const;

  // [eqn 7] Mass of (fine) roots
  S mass_root(S area_leaf) const;

  // [eqn 8] Total Mass
  S mass_live(S mass_leaf, S mass_bark, S mass_sapwood, S mass_root) const;

  S mass_total(S mass_leaf, S mass_bark, S mass_sapwood,
               S mass_heartwood, S mass_root) const;

  // Above-ground mass = leaf + all stem components (bark + sapwood +
  // heartwood); excludes roots.
  S mass_above_ground(S mass_leaf, S mass_bark,
                      S mass_sapwood, S mass_heartwood) const;

  void compute_rates(const TF24_Environment_<S>& environment, Internals_<S>& vars);

  void update_dependent_aux(const int index, Internals_<S>& vars) {
    if (index == HEIGHT_INDEX) {
      S height = vars.state(HEIGHT_INDEX);
      vars.set_aux(aux_idx_competition_effect, area_leaf(height));
      vars.set_aux(aux_idx_height_inverse, 1.0 / height);
    }
  }

  // TF24 carries no acclimating/tracked initial cohort state; this scalar-templated
  // no-op overrides the (double-only) base. (TF24f, which carries one, overrides it.)
  void set_initial_states(const TF24_Environment_<S>&, Internals_<S>&) {}

  // * Mass production
  // [eqn 12] Gross annual CO2 assimilation (double diagnostic; not on the rate path)
  double assimilation(const TF24_Environment_<S>& environment, double height,
                      double area_leaf);
  // [Appendix S6] Per-leaf photosynthetic rate.
  double assimilation_leaf(double x) const;

  // [eqn 13] Total maintenance respiration
  S respiration(S mass_leaf, S mass_sapwood, S mass_bark, S mass_root) const;

  S respiration_leaf(S mass) const;
  S respiration_bark(S mass) const;
  S respiration_sapwood(S mass) const;
  S respiration_root(S mass) const;

  // [eqn 14] Total turnover
  S turnover(S mass_leaf, S mass_bark, S mass_sapwood, S mass_root) const;
  S turnover_leaf(S mass) const;
  S turnover_bark(S mass) const;
  S turnover_sapwood(S mass) const;
  S turnover_root(S mass) const;

  // [eqn 15] Net production
  S net_mass_production_dt_A(S assimilation, S respiration, S turnover) const;

  virtual S net_mass_production_dt(const TF24_Environment_<S>& environment,
                                   S height, S area_leaf_, S height_inverse);

  // Solve the embedded (double) leaf submodel at the given size and the frozen
  // canopy the environment holds, leaving the leaf.* outputs at the operating
  // point (read by compute_rates' aux) and returning the per-leaf-area carbon
  // profit as a double. The leaf's own optimum is a plant-local, off-tape solve
  // (design decision 6); the trait/size sensitivity re-enters the reverse tape
  // linearised about the operating point in net_mass_production_dt. Virtual so
  // TF24f's tracked solve overrides the inner optimum.
  double solve_leaf_at_size(const TF24_Environment_<S>& environment,
                            double height, double area_leaf_);
  // The inner leaf-operating-point resolve (base TF24: golden-section on the
  // root-collar psi). Kept virtual so TF24f can override it while reusing
  // solve_leaf_at_size unchanged.
  virtual void solve_leaf();

  // Strategy-agnostic entry point used by Individual<TF24> (#266): reads the
  // height state and the cached aux slots itself.
  S net_mass_production_dt(const TF24_Environment_<S>& environment,
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
  S dmass_leaf_darea_leaf(S area_leaf) const;
  S dmass_sapwood_darea_leaf(S area_leaf) const;
  S dmass_bark_darea_leaf(S area_leaf) const;
  S dmass_root_darea_leaf(S area_leaf) const;
  S ddiameter_stem_darea_stem(S area_stem) const;
  S area_leaf_dt(S area_leaf_dt) const;
  S area_sapwood_dt(S area_leaf_dt) const;
  S area_heartwood_dt(S area_leaf) const;
  S area_bark_dt(S area_leaf_dt) const;
  S area_stem_dt(S area_leaf, S area_leaf_dt) const;
  S diameter_stem_dt(S area_stem, S area_stem_dt) const;
  S mass_root_dt(S area_leaf, S area_leaf_dt) const;
  S mass_live_dt(S fraction_allocation_reproduction,
                 S net_mass_production_dt) const;
  S mass_total_dt(S fraction_allocation_reproduction,
                  S net_mass_production_dt, S mass_heartwood_dt) const;
  S mass_above_ground_dt(S area_leaf, S fraction_allocation_reproduction,
                         S net_mass_production_dt, S mass_heartwood_dt,
                         S area_leaf_dt) const;

  S mass_heartwood_dt(S mass_sapwood) const;

  S mass_live_given_height(S height) const;
  S height_given_mass_leaf(S mass_leaf_) const;

  S mortality_dt(S productivity_area, S cumulative_mortality) const;
  S mortality_growth_independent_dt() const;
  S mortality_growth_dependent_dt(S productivity_area) const;
  // [eqn 20] Survival of seedlings during establishment
  S establishment_probability(const TF24_Environment_<S>& environment);

  // * Competitive environment
  S compute_competition(double z, S height) const;
  S compute_competition(double z, S area_leaf_, S height_inverse) const;
  // Strategy-agnostic entry point used by Individual<TF24> (#266).
  S compute_competition(double z, const Internals_<S>& vars) const {
    return compute_competition(z, vars.aux(aux_idx_competition_effect),
                               vars.aux(aux_idx_height_inverse));
  }

  // Per-plant census quantities the emergent kernels reduce over the cohorts,
  // reusing the allocation model: live+heartwood mass and stem basal area.
  S census_biomass(const Internals_<S>& vars) const {
    return mass_live_given_height(vars.state(HEIGHT_INDEX)) +
           vars.state(state_idx_mass_heartwood);
  }
  S census_basal_area(const Internals_<S>& vars) const {
    const S area_leaf_ = area_leaf(vars.state(HEIGHT_INDEX));
    return area_stem(area_bark(area_leaf_), area_sapwood(area_leaf_),
                     vars.state(state_idx_area_heartwood));
  }

  // [eqn  9] Probability density of leaf area at height `z`
  S q(S z, S height) const;
  // [eqn 10] Fraction of leaf area above height `z`
  S Q(S z, S height, S eta_x) const;
  // [      ] Inverse of Q: height above which fraction 'x' of leaf found
  S Qp(S x, S height) const;

  // The aim is to find a plant height that gives the correct seed mass.
  S height_seed(void) const;

  // Set constants within TF24_Strategy
  void prepare_strategy();

  // Birth height of a (germinated) seed.
  S initial_height() const { return height_0; }

  // Crown shading model, resolved once from control.shading_model.
  ShadingModel shading_model_ = ShadingModel::MeanLight;

  // Whether the reverse-mode leaf edge also linearises the per-area profit in the
  // (active) plant height (the size-feedback channel). On by default: the leaf's
  // per-area profit depends strongly on plant size through the hydraulic supply
  // (leaf-specific conductance ~ 1/height, sapwood volume ~ height) and the crown
  // radiation, so this channel dominates a trait's effect on growth.
  bool add_size_feedback_ = true;

  // odelia differentiable-System handles ------------------------------------
  template <class S2> using rebind = TF24_Strategy_<S2>;

  template <class S2>
  TF24_Strategy_<S2> rebind_from() const {
    TF24_Strategy_<S2> out;
    out.pars = pars.template rebind<S2>();
    out.control = control;
    out.name = name;
    out.birth_rate_x = birth_rate_x;
    out.birth_rate_y = birth_rate_y;
    out.is_variable_birth_rate = is_variable_birth_rate;
    out.collect_all_auxiliary = collect_all_auxiliary;
    out.refresh_indices();
    return out;
  }

  // Handles to the differentiable trait fields of this one shared TF24_Pars.
  std::vector<S*> ad_parameters() { return pars.field_ptrs(); }

  // Biological (user-settable) parameters; see TF24_Pars above.
  TF24_Pars_<S> pars;

  // Derived / precomputed in prepare_strategy() (NOT user-set) -------------
  S eta_c     = NA_REAL; // crown shape factor, precomputed from pars.eta
  S height_0  = NA_REAL;
  S area_leaf_0;

  // Embedded leaf hydraulic/photosynthesis sub-model, built in prepare_strategy().
  // Always double: its own optimum is a plant-local forward-mode solve and never
  // goes reverse-active (design decision 6).
  Leaf leaf;

  // Hydraulic root parameters (not exposed to R; always double)
  double root_c = 2.680147;
  double root_b = 3.898245;
  double root_psi_crit = root_b*std::pow(log(1/0.05),1/root_c);

  // Solver tolerances and other constants not exposed to R (always double)
  double newton_tol_abs = 0.001;
  double GSS_tol_abs = 1e-3;
  double vulnerability_curve_ncontrol = 100;
  double ci_abs_tol = 1e-6;
  double ci_niter = 1000;
  double beta_R_H = 3.4e2;
  double beta_R_V = 9.4e3;

  // Cached aux/state indices, resolved once in refresh_indices().
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

  // Reusable per-layer root-mass buffer, refilled each solve_leaf_at_size call.
  std::vector<double> mass_root_prop_;
};

// S=double is the resident numerics that cross the R boundary; the alias keeps
// this name a concrete type for RcppR6 and every existing caller.
using TF24_Strategy = TF24_Strategy_<double>;

TF24_Strategy::ptr make_strategy_ptr(TF24_Strategy s);

// Active strategies are instantiated only where a trait gradient is taken.
template <class S>
typename TF24_Strategy_<S>::ptr make_strategy_ptr(TF24_Strategy_<S> s) {
  s.prepare_strategy();
  return std::make_shared<TF24_Strategy_<S>>(s);
}

}

#endif
