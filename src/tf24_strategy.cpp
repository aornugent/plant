// Built from  src/ff16_strategy.cpp on Mon Feb 12 09:52:27 2024 using the scaffolder, from the strategy:  FF16
#include <plant/models/tf24_strategy.h>
#include <odelia/supplied_derivative.hpp> // envelope-theorem leaf seam (§7.3)
#include <limits> // std::numeric_limits (height_seed bounds)
#include <cmath>  // std::abs (FD step)

namespace plant {

// Full AD strip to double at the Leaf boundary (the Leaf is never templated on
// S); a no-op on the double path.
using odelia::util::to_passive;

// --- Hard-coded root-distribution constants (review #9) ---------------------
// Named here for clarity; promotion to user-tunable traits (RcppR6) is a
// deliberate follow-up (see vignettes/models/code_review_leaf_tf24.qmd #9).
// rescales total fine-root mass into the per-layer carbon units expected by the
// root hydraulic network in Leaf::set_physiology.
static const double root_mass_carbon_scale = 83.26 * 0.5;
// rooting depth cap (m), i.e. the depth of the soil column.
static const double rooting_depth_max = 1.5;

// NOTE (review #9): the per-second -> annual factor 60*60*12*365 (seconds of
// daylight per year, 12 h day x 365 d) recurs in compute_rates and
// net_mass_production_dt below. It is deliberately left inline rather than
// hoisted to a constant: collapsing the 4-step integer product into one double
// changes the floating-point rounding, and the adaptive ODE amplifies it
// (offspring_production shifts ~0.2%). Kept inline to preserve bit-identical
// results.

// TODO: Document consistent argument order: l, b, s, h, r
// TODO: Document ordering of different types of variables (size
// before physiology, before compound things?)
// TODO: Consider moving to activating as an initialisation list?
template <class S>
TF24_Strategy_<S>::TF24_Strategy_() {
  this->collect_all_auxiliary = false;
  // build the string state/aux name to index map
  refresh_indices();
  this->name = "TF24";
}

// not sure 'average' is the right term here..
// Off the S rate path: this weights the double Leaf's crown light aggregation,
// so it strips the (possibly active) light read to double.
template <class S>
double TF24_Strategy_<S>::compute_average_light_environment(
    double z, double height, const environment_type &environment) {
// NOTE: the light environment is clamped to a small positive floor (1e-4)
// rather than allowed to reach 0 (original rationale was never recorded;
// preserved as-is).

     return std::max(to_passive(environment.get_environment_at_height(z)), 0.0001) * q(z, height);
}

// assumes optimise_psi_stem_TF has been run for optimal psi_stem
template <class S>
S TF24_Strategy_<S>::evapotranspiration_dt(S area_leaf_, int soil_layer) {
  return leaf.soil_consumption_[soil_layer] * area_leaf_;
}

template <class S>
void TF24_Strategy_<S>::refresh_indices () {
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

  // Cache integer indices for the keys used in the hot compute_rates path, so
  // it no longer does a std::map<string,int> lookup per derivs evaluation.
  aux_idx_competition_effect    = this->aux_index.at("competition_effect");
  aux_idx_height_inverse        = this->aux_index.at("height_inverse");
  aux_idx_net_mass_production_dt = this->aux_index.at("net_mass_production_dt");
  aux_idx_root_mass             = this->aux_index.at("root_mass");
  aux_idx_opt_psi_stem          = this->aux_index.at("opt_psi_stem");
  aux_idx_opt_root_psi          = this->aux_index.at("opt_root_psi");
  aux_idx_transpiration         = this->aux_index.at("transpiration");
  aux_idx_E_up                  = this->aux_index.at("E_up_");
  aux_idx_profit                = this->aux_index.at("profit");
  aux_idx_stom_cond_CO2         = this->aux_index.at("stom_cond_CO2");
  // area_sapwood is only registered when collect_all_auxiliary is set.
  aux_idx_area_sapwood = this->aux_index.count("area_sapwood") ? this->aux_index.at("area_sapwood") : -1;
  state_idx_area_heartwood      = this->state_index.at("area_heartwood");
  state_idx_mass_heartwood      = this->state_index.at("mass_heartwood");
}

// [eqn 2] area_leaf (inverse of [eqn 3])
template <class S>
S TF24_Strategy_<S>::area_leaf(S height) const {
  return pow(height / pars.a_l1, 1.0 / pars.a_l2);
}

// [eqn 1] mass_leaf (inverse of [eqn 2])
template <class S>
S TF24_Strategy_<S>::mass_leaf(S area_leaf) const {
  return area_leaf * pars.lma;
}

// [eqn 4] area and mass of sapwood
template <class S>
S TF24_Strategy_<S>::area_sapwood(S area_leaf) const {
  return area_leaf * pars.theta;
}

template <class S>
S TF24_Strategy_<S>::mass_sapwood(S area_sapwood, S height) const {
  return area_sapwood * height * eta_c * pars.rho;
}

// [eqn 5] area and mass of bark
template <class S>
S TF24_Strategy_<S>::area_bark(S area_leaf) const {
  return pars.a_b1 * area_leaf * pars.theta;
}

template <class S>
S TF24_Strategy_<S>::mass_bark(S area_bark, S height) const {
  return area_bark * height * eta_c * pars.rho;
}

template <class S>
S TF24_Strategy_<S>::area_stem(S area_bark, S area_sapwood,
                            S area_heartwood) const {
  return area_bark + area_sapwood + area_heartwood;
}

template <class S>
S TF24_Strategy_<S>::diameter_stem(S area_stem) const {
  using std::sqrt;
  return sqrt(4 * area_stem / M_PI);
}

// [eqn 7] Mass of (fine) roots
template <class S>
S TF24_Strategy_<S>::mass_root(S area_leaf) const {
  return pars.a_r1 * area_leaf;
}

// [eqn 8] Total mass
template <class S>
S TF24_Strategy_<S>::mass_live(S mass_leaf, S mass_bark,
                           S mass_sapwood, S mass_root) const {
  return mass_leaf + mass_sapwood + mass_bark + mass_root;
}

template <class S>
S TF24_Strategy_<S>::mass_total(S mass_leaf, S mass_bark,
                            S mass_sapwood, S mass_heartwood,
                            S mass_root) const {
  return mass_leaf + mass_bark + mass_sapwood +  mass_heartwood + mass_root;
}

template <class S>
S TF24_Strategy_<S>::mass_above_ground(S mass_leaf, S mass_bark,
                            S mass_sapwood, S mass_heartwood) const {
  return mass_leaf + mass_bark + mass_sapwood + mass_heartwood;
}

// for updating auxiliary state
template <class S>
void TF24_Strategy_<S>::update_dependent_aux(const int index, Internals_<S>& vars) {
  if (index == HEIGHT_INDEX) {
    S height = vars.state(HEIGHT_INDEX);
    vars.set_aux(aux_idx_competition_effect, area_leaf(height));
    vars.set_aux(aux_idx_height_inverse, 1.0 / height);
  }
}


// one-shot update of the scm variables
// i.e. setting rates of ode vars from the state and updating aux vars
template <class S>
void TF24_Strategy_<S>::compute_rates(const environment_type& environment,  Internals_<S>& vars) {
  S height = vars.state(HEIGHT_INDEX);
  S area_leaf_ = vars.aux(aux_idx_competition_effect);

  const S net_mass_production_dt_ =
    net_mass_production_dt(environment, height, area_leaf_,
                           vars.aux(aux_idx_height_inverse));

  // store the aux sate
  vars.set_aux(aux_idx_net_mass_production_dt, net_mass_production_dt_);
  vars.set_aux(aux_idx_root_mass, mass_root(area_leaf_));
  // Leaf diagnostics (double); recorded as aux for reporting only.
  vars.set_aux(aux_idx_opt_psi_stem, leaf.opt_psi_stem_);
  vars.set_aux(aux_idx_opt_root_psi, leaf.root_collar_psi_);
  vars.set_aux(aux_idx_transpiration, leaf.transpiration_);
  vars.set_aux(aux_idx_E_up, leaf.E_up_);
  vars.set_aux(aux_idx_profit, leaf.profit_);
  vars.set_aux(aux_idx_stom_cond_CO2, leaf.stom_cond_CO2_);




  // consumption rates should be emerging from net_mass_produciton_dt
  // convert evapotranspiration per leaf area per soil layer (mol H20 m^-2 s^-1) to canopy-level total
  // yearly evapotranspiration per soil layer (m yr^-1)
  // stubbing out E_p for integration
  int soil_number_of_depths_ = environment.get_soil_number_of_depths();


  for (int i = 0; i < soil_number_of_depths_; i++) {

    // evapotranspiration (mol H20 m^-2 s^-1 layer^-1)
    // consumption rate (m yr^-1 layer ^-1)
    vars.set_consumption_rate(i, evapotranspiration_dt(area_leaf_, i)*60*60*12*365/1000*kg_per_mol_h2o);
  }

  if (net_mass_production_dt_ > 0) {

    const S fraction_allocation_reproduction_ = fraction_allocation_reproduction(height);
    const S darea_leaf_dmass_live_ = darea_leaf_dmass_live(area_leaf_);
    const S fraction_allocation_growth_ = fraction_allocation_growth(height);
    const S area_leaf_dt = net_mass_production_dt_ * fraction_allocation_growth_ * darea_leaf_dmass_live_;

    vars.set_rate(HEIGHT_INDEX, dheight_darea_leaf(area_leaf_) * area_leaf_dt);
    vars.set_rate(FECUNDITY_INDEX,
      fecundity_dt(net_mass_production_dt_, fraction_allocation_reproduction_));

    vars.set_rate(state_idx_area_heartwood, area_heartwood_dt(area_leaf_));
    const S area_sapwood_ = area_sapwood(area_leaf_);
    const S mass_sapwood_ = mass_sapwood(area_sapwood_, height);
    vars.set_rate(state_idx_mass_heartwood, mass_heartwood_dt(mass_sapwood_));

    if (this->collect_all_auxiliary) {
      vars.set_aux(aux_idx_area_sapwood, area_sapwood_);
    }
  } else {
    vars.set_rate(HEIGHT_INDEX, 0.0);
    vars.set_rate(FECUNDITY_INDEX, 0.0);
    vars.set_rate(state_idx_area_heartwood, 0.0);
    vars.set_rate(state_idx_mass_heartwood, 0.0);
  }
  // [eqn 21] - Instantaneous mortality rate
  vars.set_rate(MORTALITY_INDEX,
      mortality_dt(net_mass_production_dt_ / area_leaf_, vars.state(MORTALITY_INDEX)));

}

// [eqn 12] Gross annual CO2 assimilation (!!not in use for TF24 model!!). Kept
// double: it is off the rate path (net production uses the Leaf submodel), and
// the light read is stripped so the unused overload compiles at an active S.
template <class S>
S TF24_Strategy_<S>::assimilation(const environment_type& environment,
                                    S height,
                                    S area_leaf) {


  double A = 0.0;
  const double height_d = to_passive(height);

  // Define an anonymous function to integrate
  // For given height in crown, take photosynthesis at depth multipled by
  //   amount of leaf at that depth
  std::function<double(double)> f = [&](double z) -> double {
    return to_passive(assimilation_leaf(S(to_passive(environment.get_environment_at_height(S(z)))))) * q(z, height_d);
  };

  // Integrate over crown depth using using Gauss-Kronrod quadrature.
  // The number of points used in the integration is determined by the control parameter
  // function_integration_rule. Rules defined in qk_rules.cpp
  A = function_integrator.integrate(f, 0.0, height_d);

  return area_leaf * A;
}

// Photosynthetic rate per leaf area
// `x` is openness, ranging from 0 to 1.
template <class S>
S TF24_Strategy_<S>::assimilation_leaf(S x) const {
  return pars.a_p1 * x / (x + pars.a_p2);
}

// [eqn 13] Total maintenance respiration
// NOTE: In contrast with Falster ref model, we do not normalise by pars.a_y*pars.a_bio.
template <class S>
S TF24_Strategy_<S>::respiration(S mass_leaf, S mass_sapwood,
                             S mass_bark, S mass_root) const {
  return respiration_leaf(mass_leaf) +
         respiration_bark(mass_bark) +
         respiration_sapwood(mass_sapwood) +
         respiration_root(mass_root);
}

template <class S>
S TF24_Strategy_<S>::respiration_leaf(S mass) const {
  return pars.r_l * mass;
}

template <class S>
S TF24_Strategy_<S>::respiration_bark(S mass) const {
  return pars.r_b * mass;
}

template <class S>
S TF24_Strategy_<S>::respiration_sapwood(S mass) const {
  return pars.r_s * mass;
}

template <class S>
S TF24_Strategy_<S>::respiration_root(S mass) const {
  return pars.r_r * mass;
}

// [eqn 14] Total turnover
template <class S>
S TF24_Strategy_<S>::turnover(S mass_leaf, S mass_bark,
                          S mass_sapwood, S mass_root) const {
   return turnover_leaf(mass_leaf) +
          turnover_bark(mass_bark) +
          turnover_sapwood(mass_sapwood) +
          turnover_root(mass_root);
}

template <class S>
S TF24_Strategy_<S>::turnover_leaf(S mass) const {
  return pars.k_l * mass;
}

template <class S>
S TF24_Strategy_<S>::turnover_bark(S mass) const {
  return pars.k_b * mass;
}

template <class S>
S TF24_Strategy_<S>::turnover_sapwood(S mass) const {
  return pars.k_s * mass;
}

template <class S>
S TF24_Strategy_<S>::turnover_root(S mass) const {
  return pars.k_r * mass;
}

// [eqn 15] Net production
//
// NOTE: Translation of variable names from the Falster 2011.  Everything
// before the minus sign is SCM's N, our `net_mass_production_dt` is SCM's P.
template <class S>
S TF24_Strategy_<S>::net_mass_production_dt_A(S assimilation, S respiration,
                                S turnover) const {
  return pars.a_bio * pars.a_y * (assimilation - respiration) - turnover;
}

// One shot calculation of net_mass_production_dt
// Used by establishment_probability() and compute_rates().
template <class S>
S TF24_Strategy_<S>::net_mass_production_dt(const environment_type& environment,
                                S height, S area_leaf_,
                                S height_inverse) {
  // height_inverse (= 1/height) is supplied by the shared individual.h interface
  // (cached aux); unused here as the TF24 root-water path works in height directly.
  (void)height_inverse;
  const S mass_leaf_    = mass_leaf(area_leaf_);
  const S area_sapwood_ = area_sapwood(area_leaf_);
  const S mass_sapwood_ = mass_sapwood(area_sapwood_, height);
  const S area_bark_    = area_bark(area_leaf_);
  const S mass_bark_    = mass_bark(area_bark_, height);
  const S mass_root_    = mass_root(area_leaf_);

  int soil_number_of_depths_ = environment.get_soil_number_of_depths();
  const std::vector<double>& soil_depths_ = environment.z;

  // --------------------------------------------------------------------------
  // Leaf hydraulic sub-model. Everything from here to the assimilation_ line is
  // DOUBLE: the Leaf is never templated on S (design 4.3). The S mass cascade
  // re-enters only via the area_leaf_ factor in assimilation_ below; every
  // seeded-parameter sensitivity of the Leaf outputs is injected onto the tape
  // by supplied_derivative in a later stage. On the double path this is
  // bit-identical to the pre-template code (to_passive is the identity).
  // --------------------------------------------------------------------------
  const double height_d    = to_passive(height);
  const double area_leaf_d = to_passive(area_leaf_);
  const double mass_root_d = to_passive(mass_root_);

  // psi_soil (-MPa): strip AD off the (possibly active) soil state for the leaf.
  const std::vector<S>& psi_soil_S = environment.get_soil_water_potential_state();
  std::vector<double> psi_soil(psi_soil_S.size());
  for (size_t i = 0; i < psi_soil_S.size(); ++i) {
    psi_soil[i] = to_passive(psi_soil_S[i]);
  }

  // find leaf specific max hydraulic conductance (kg m^-2 LA s^-1 MPa ^-1)
  const double leaf_specific_conductance_max =
      to_passive(pars.K_s * pars.theta / (height * eta_c));

  // find sapwood volume per leaf area
  const double sapwood_volume_per_leaf_area =
      to_passive(pars.theta * (height * eta_c));

  // ----------------------------------------------------------------------
  // ROOT MASS DISTRIBUTION ACROSS SOIL LAYERS (double; feeds the double leaf)
  // ----------------------------------------------------------------------
  // Total fine-root mass (mass_root_) is distributed over depth using the same
  // cumulative shape function Q() used for the leaf canopy, but parameterised
  // over soil depth. See the original derivation for the units rescaling.
  mass_root_prop_.assign(soil_number_of_depths_, 0.0);

  double rooting_depth = std::min(height_d, rooting_depth_max);
  const double root_mass_scale = root_mass_carbon_scale * mass_root_d;
  const double root_depth_shape_eta_d = to_passive(pars.root_depth_shape_eta);

  double prev_q = 1.0;
  for (int a = 0; a < soil_number_of_depths_; ++a) {
    if(prev_q == 0){
      break;
    }
    const double qd = Q(soil_depths_[a], rooting_depth, root_depth_shape_eta_d);

    mass_root_prop_[a] = root_mass_scale * (prev_q - qd);
    prev_q = qd;
  }

  // Reuse geometry precomputed by environment; avoids rebuilding z midpoints each call.
  leaf.z_soil_mid_ = environment.get_soil_mid_depths();
  leaf.use_precomputed_z_soil_mid_ = true;

  const double rho_d   = to_passive(pars.rho);
  const double a_bio_d = to_passive(pars.a_bio);
  const double k_I_d   = to_passive(pars.k_I);

  // Optimise the leaf at a given absorbed radiation (see original notes).
  auto optimise_at = [&](double radiation) {
    leaf.set_physiology(area_leaf_d, mass_root_prop_, rho_d, a_bio_d, radiation, psi_soil, soil_depths_, leaf_specific_conductance_max, environment.get_atm_vpd(), environment.get_ca(), sapwood_volume_per_leaf_area, environment.get_leaf_temp(), environment.get_atm_o2_kpa(), environment.get_atm_kpa());
    solve_leaf();
  };

  // Convert canopy openness (0-1) into absorbed radiation.
  const double PPFD = environment.get_PPFD();
  auto radiation_at = [&](double light) {
    return k_I_d * std::max(light, 0.0001) * PPFD;
  };

  // Aggregate the leaf submodel over the crown according to the shading model.
  // single_solve models run one leaf optimisation at a scalar radiation; the
  // active leaf seam below re-evaluates the leaf at that same radiation.
  double radiation_used = NA_REAL;
  bool single_solve = true;
  if (shading_model_ == ShadingModel::CrownCentre) {
    radiation_used = radiation_at(to_passive(environment.get_environment_at_height(height * eta_c)));
    optimise_at(radiation_used);
  } else if (shading_model_ == ShadingModel::MeanLight) {
    // Leaf-area-weighted mean canopy openness = integral of (light * q) over the
    // crown (q integrates to one).
    auto f = [&](double x) -> double {
      return compute_average_light_environment(x, height_d, environment);
    };
    radiation_used = radiation_at(function_integrator.integrate(f, 0.0, height_d));
    optimise_at(radiation_used);
  } else { // DeepCrown
    single_solve = false;
    const std::vector<double> nodes =
      function_integrator.integrate_vector_x(0.0, height_d);
    const size_t nn = nodes.size();
    std::vector<double> profit_y(nn), trans_y(nn), eup_y(nn), psi_y(nn),
      root_psi_y(nn), gco2_y(nn);
    std::vector<std::vector<double>> soil_y(
      soil_number_of_depths_, std::vector<double>(nn));
    for (size_t i = 0; i < nn; ++i) {
      const double qi = q(nodes[i], height_d);
      optimise_at(radiation_at(to_passive(environment.get_environment_at_height(S(nodes[i])))));
      profit_y[i]   = leaf.profit_ * qi;
      trans_y[i]    = leaf.transpiration_ * qi;
      eup_y[i]      = leaf.E_up_ * qi;
      psi_y[i]      = leaf.opt_psi_stem_ * qi;
      root_psi_y[i] = leaf.root_collar_psi_ * qi;
      gco2_y[i]     = leaf.stom_cond_CO2_ * qi;
      for (int a = 0; a < soil_number_of_depths_; ++a) {
        soil_y[a][i] = leaf.soil_consumption_[a] * qi;
      }
    }
    // Integrate each leaf output to its leaf-area-weighted crown mean.
    leaf.profit_          = function_integrator.integrate_vector(profit_y, 0.0, height_d);
    leaf.transpiration_   = function_integrator.integrate_vector(trans_y, 0.0, height_d);
    leaf.E_up_            = function_integrator.integrate_vector(eup_y, 0.0, height_d);
    leaf.opt_psi_stem_    = function_integrator.integrate_vector(psi_y, 0.0, height_d);
    leaf.root_collar_psi_ = function_integrator.integrate_vector(root_psi_y, 0.0, height_d);
    leaf.stom_cond_CO2_   = function_integrator.integrate_vector(gco2_y, 0.0, height_d);
    for (int a = 0; a < soil_number_of_depths_; ++a) {
      leaf.soil_consumption_[a] =
        function_integrator.integrate_vector(soil_y[a], 0.0, height_d);
    }
  }


  //TODO: one point constant ratio and integral width for daylength
  // convert assimilation per leaf area per second (umol m^-2 s^-1) to canopy-level total yearly assimilation (mol yr^-1)
  // converts to canopy area, then years, then mols. leaf.profit_ is double.
  //
  // The active leaf supplied_derivative seam (§7.3): the Leaf is double, so its
  // profit carries no parameter derivative on its own. On the active path we wrap
  // profit_ in a supplied_derivative carrying d(profit)/d(input) for every seeded
  // input that reaches the leaf, obtained by central FD of leaf_profit_frozen at
  // the FROZEN optimum collar psi (envelope theorem: d(profit)/d(psi*) = 0, so no
  // re-optimise). The §7.4 pair-filter (shouldRecord()) drops frozen/mutant inputs
  // whose slot is INVALID. The S mass cascade re-enters through the area_leaf_
  // factor below, independent of this seam. On the double path this is exactly
  // leaf.profit_ (to_passive/if constexpr collapse to the identity).
  S profit_s = S(leaf.profit_);
  if constexpr (!std::is_same_v<S, double>) {
    auto* tape = xad::Tape<double>::getActive();
    if (tape != nullptr && tape->isActive()) {
      if (!single_solve) {
        util::stop("TF24 active leaf gradient: the deep-crown seam is not yet "
                   "implemented; use shading_model 'mean-light' or 'crown-centre'.");
      }
      const double frozen_collar_psi = -leaf.root_collar_psi_;
      const double profit0 = leaf.profit_;
      const double atm_vpd_d = environment.get_atm_vpd();
      const double ca_d = environment.get_ca();
      const double leaf_temp_d = environment.get_leaf_temp();
      const double atm_o2_d = environment.get_atm_o2_kpa();
      const double atm_kpa_d = environment.get_atm_kpa();
      const std::vector<double> z_soil_mid_ = environment.get_soil_mid_depths();

      // Baseline double parameter set + double field pointers, in the SAME
      // TF24_AD_FIELDS order as field_ptrs(), so input i pairs with partial i.
      TF24_Pars_<double> pd;
#define PLANT_AD_PDPTR(f) &pd.f,
      std::vector<double*> pdp = { TF24_AD_FIELDS(PLANT_AD_PDPTR) };
#undef PLANT_AD_PDPTR
      std::vector<S*> ap = field_ptrs();
      for (std::size_t i = 0; i < ap.size(); ++i) *pdp[i] = to_passive(*ap[i]);

      auto profit_frozen = [&](void) {
        return leaf_profit_frozen(pd, height_d, radiation_used, psi_soil,
                                  soil_depths_, z_soil_mid_, atm_vpd_d, ca_d,
                                  leaf_temp_d, atm_o2_d, atm_kpa_d,
                                  frozen_collar_psi);
      };

      std::vector<S*> inputs;
      std::vector<double> partials;
      inputs.reserve(ap.size() + 1);
      partials.reserve(ap.size() + 1);

      // Per-parameter central FD at the frozen optimum. Only seeded inputs (valid
      // tape slot) are wired; the rest are the mutant/background frozen inputs the
      // §7.4 pair-filter drops. Every field gets a partial computed, so a seeded
      // parameter is never an implicit zero (M1).
      for (std::size_t i = 0; i < ap.size(); ++i) {
        const double base = *pdp[i];
        const double hstep = 1e-6 * (std::abs(base) + 1e-6);
        *pdp[i] = base + hstep;
        const double pplus = profit_frozen();
        *pdp[i] = base - hstep;
        const double pminus = profit_frozen();
        *pdp[i] = base;
        if (ap[i]->shouldRecord()) {
          inputs.push_back(ap[i]);
          partials.push_back((pplus - pminus) / (2.0 * hstep));
        }
      }

      // Height channel: the leaf profit depends on the current height state
      // (conductance, sapwood volume, area_leaf), so the state coupling through
      // the leaf must be on the tape. Perturb the frozen height directly.
      if (height.shouldRecord()) {
        const double hbase = height_d;
        const double hstep = 1e-6 * (std::abs(hbase) + 1e-6);
        const double pplus = leaf_profit_frozen(pd, hbase + hstep, radiation_used,
            psi_soil, soil_depths_, z_soil_mid_, atm_vpd_d, ca_d, leaf_temp_d,
            atm_o2_d, atm_kpa_d, frozen_collar_psi);
        const double pminus = leaf_profit_frozen(pd, hbase - hstep, radiation_used,
            psi_soil, soil_depths_, z_soil_mid_, atm_vpd_d, ca_d, leaf_temp_d,
            atm_o2_d, atm_kpa_d, frozen_collar_psi);
        inputs.push_back(&height);
        partials.push_back((pplus - pminus) / (2.0 * hstep));
      }

      if (!inputs.empty()) {
        profit_s = odelia::ode::supplied_derivative(*tape, profit0, inputs, partials);
      }
    }
  }
  const S assimilation_ = profit_s * area_leaf_* 60*60*12*365/1e6;
  // const double assimilation_ = assimilation(environment, height, area_leaf_);
  const S respiration_ =
    respiration(mass_leaf_, mass_sapwood_, mass_bark_, mass_root_);
  const S turnover_ =
    turnover(mass_leaf_, mass_bark_, mass_sapwood_, mass_root_);
  return net_mass_production_dt_A(assimilation_, respiration_, turnover_);
}

// Base TF24: optimise the root-collar water potential from scratch each call.
template <class S>
void TF24_Strategy_<S>::solve_leaf() {
  leaf.find_root_collar_psi();
}

// Lift the double birth height to S carrying its parameter derivative via the
// implicit function theorem (§7.2). Mirrors FF16_Strategy_::lift_birth_height.
template <class S>
S TF24_Strategy_<S>::lift_birth_height(double h_star) const {
  if constexpr (std::is_same_v<S, double>) {
    return h_star;
  } else {
    // dg/dh at the root, as a double constant (central difference; g is smooth
    // and near-linear in h here, so a 2-point rule is not the accuracy limit).
    const double eps = 1e-6 * (h_star + 1.0);
    auto gv = [&](double h) { return to_passive(mass_live_given_height(S(h))); };
    const double dgdh = (gv(h_star + eps) - gv(h_star - eps)) / (2.0 * eps);
    // Newton correction from the double root. Subtract its own value so the
    // birth height's value stays exactly h_star (bit-identical), while its
    // derivative is the IFT term dh*/dtheta = -(dg/dtheta)/(dg/dh).
    const S corr = (mass_live_given_height(S(h_star)) - pars.omega) / dgdh;
    return S(h_star) - corr + to_passive(corr);
  }
}

// FD engine for the leaf supplied_derivative seam (§7.3). Reconstruct a
// throwaway double Leaf from `pd` and the frozen crown context, evaluate profit
// at the FROZEN optimum collar psi, and return it. Mirrors the single-solve
// (mean-light / crown-centre) leaf setup in net_mass_production_dt exactly, but
// with the parameters taken from `pd` (so a perturbed field flows through every
// channel it reaches) and the collar psi held fixed at the forward optimum.
template <class S>
double TF24_Strategy_<S>::leaf_profit_frozen(
    const TF24_Pars_<double>& pd, double height_d, double radiation,
    const std::vector<double>& psi_soil_d,
    const std::vector<double>& soil_depths_,
    const std::vector<double>& z_soil_mid_, double atm_vpd, double ca,
    double leaf_temp, double atm_o2_kpa, double atm_kpa,
    double frozen_collar_psi) {

  const int soil_number_of_depths_ = static_cast<int>(soil_depths_.size());

  // Derived crown/allometry quantities recomputed from pd so a perturbed a_l1 /
  // a_l2 / eta / theta / K_s / a_r1 / root_depth_shape_eta flows into the leaf.
  const double eta_c_d = 1 - 2 / (1 + pd.eta) + 1 / (1 + 2 * pd.eta);
  const double area_leaf_d = std::pow(height_d / pd.a_l1, 1.0 / pd.a_l2);
  const double mass_root_d = pd.a_r1 * area_leaf_d;
  const double leaf_specific_conductance_max =
      pd.K_s * pd.theta / (height_d * eta_c_d);
  const double sapwood_volume_per_leaf_area = pd.theta * (height_d * eta_c_d);

  // Root-mass distribution across soil layers (mirrors net_mass_production_dt).
  std::vector<double> mass_root_prop(soil_number_of_depths_, 0.0);
  const double rooting_depth = std::min(height_d, rooting_depth_max);
  const double root_mass_scale = root_mass_carbon_scale * mass_root_d;
  double prev_q = 1.0;
  for (int a = 0; a < soil_number_of_depths_; ++a) {
    if (prev_q == 0) break;
    const double qd = Q(soil_depths_[a], rooting_depth, pd.root_depth_shape_eta);
    mass_root_prop[a] = root_mass_scale * (prev_q - qd);
    prev_q = qd;
  }

  Leaf L(pd.vcmax_25, pd.c, pd.b, pd.psi_crit, root_c, root_b, root_psi_crit,
         pd.beta2, pd.jmax_25, pd.a, pd.curv_fact_elec_trans,
         pd.curv_fact_colim, GSS_tol_abs, vulnerability_curve_ncontrol,
         ci_abs_tol, ci_niter, pd.g1_TF24, beta_R_H, beta_R_V);
  L.z_soil_mid_ = z_soil_mid_;
  L.use_precomputed_z_soil_mid_ = true;
  L.set_physiology(area_leaf_d, mass_root_prop, pd.rho, pd.a_bio, radiation,
                   psi_soil_d, soil_depths_, leaf_specific_conductance_max,
                   atm_vpd, ca, sapwood_volume_per_leaf_area, leaf_temp,
                   atm_o2_kpa, atm_kpa);
  L.evaluate_root_collar_psi(frozen_collar_psi);
  return L.profit_;
}

// [eqn 16] Fraction of production allocated to reproduction
template <class S>
S TF24_Strategy_<S>::fraction_allocation_reproduction(S height) const {
  return pars.a_f1 / (1.0 + exp(pars.a_f2 * (1.0 - height / pars.hmat)));
}

// Fraction of production allocated to growth
template <class S>
S TF24_Strategy_<S>::fraction_allocation_growth(S height) const {
  return 1.0 - fraction_allocation_reproduction(height);
}

// [eqn 17] Rate of offspring production
template <class S>
S TF24_Strategy_<S>::fecundity_dt(S net_mass_production_dt,
                               S fraction_allocation_reproduction) const {
  return net_mass_production_dt * fraction_allocation_reproduction /
    (pars.omega + pars.a_f3);
}

template <class S>
S TF24_Strategy_<S>::darea_leaf_dmass_live(S area_leaf) const {
  return 1.0/(  dmass_leaf_darea_leaf(area_leaf)
              + dmass_sapwood_darea_leaf(area_leaf)
              + dmass_bark_darea_leaf(area_leaf)
              + dmass_root_darea_leaf(area_leaf));
}

template <class S>
S TF24_Strategy_<S>::dheight_darea_leaf(S area_leaf) const {
  return pars.a_l1 * pars.a_l2 * pow(area_leaf, pars.a_l2 - 1);
}

// Mass of leaf needed for new unit area leaf, d m_s / d a_l
template <class S>
S TF24_Strategy_<S>::dmass_leaf_darea_leaf(S /* area_leaf */) const {
  return pars.lma;
}

// Mass of stem needed for new unit area leaf, d m_s / d a_l
template <class S>
S TF24_Strategy_<S>::dmass_sapwood_darea_leaf(S area_leaf) const {
  return pars.rho * eta_c * pars.a_l1 * pars.theta * (pars.a_l2 + 1.0) * pow(area_leaf, pars.a_l2);
}

// Mass of bark needed for new unit area leaf, d m_b / d a_l
template <class S>
S TF24_Strategy_<S>::dmass_bark_darea_leaf(S area_leaf) const {
  return pars.a_b1 * dmass_sapwood_darea_leaf(area_leaf);
}

// Mass of root needed for new unit area leaf, d m_r / d a_l
template <class S>
S TF24_Strategy_<S>::dmass_root_darea_leaf(S /* area_leaf */) const {
  return pars.a_r1;
}

// Growth rate of basal diameter_stem per unit time
template <class S>
S TF24_Strategy_<S>::ddiameter_stem_darea_stem(S area_stem) const {
  return pow(M_PI * area_stem, -0.5);
}

// Growth rate of sapwood area at base per unit time
template <class S>
S TF24_Strategy_<S>::area_sapwood_dt(S area_leaf_dt) const {
  return area_leaf_dt * pars.theta;
}

// Note, unlike others, heartwood growth does not depend on leaf area growth, but
// rather existing sapwood
template <class S>
S TF24_Strategy_<S>::area_heartwood_dt(S area_leaf) const {
  return pars.k_s * area_sapwood(area_leaf);
}

// Growth rate of bark area at base per unit time
template <class S>
S TF24_Strategy_<S>::area_bark_dt(S area_leaf_dt) const {
  return pars.a_b1 * area_leaf_dt * pars.theta;
}

// Growth rate of stem basal area per unit time
template <class S>
S TF24_Strategy_<S>::area_stem_dt(S area_leaf,
                               S area_leaf_dt) const {
  return area_sapwood_dt(area_leaf_dt) +
    area_bark_dt(area_leaf_dt) +
    area_heartwood_dt(area_leaf);
}

// Growth rate of basal diameter_stem per unit time
template <class S>
S TF24_Strategy_<S>::diameter_stem_dt(S area_stem, S area_stem_dt) const {
  return ddiameter_stem_darea_stem(area_stem) * area_stem_dt;
}

// Growth rate of root mass per unit time
template <class S>
S TF24_Strategy_<S>::mass_root_dt(S area_leaf,
                               S area_leaf_dt) const {
  return area_leaf_dt * dmass_root_darea_leaf(area_leaf);
}

template <class S>
S TF24_Strategy_<S>::mass_live_dt(S fraction_allocation_reproduction,
                               S net_mass_production_dt) const {
  return (1 - fraction_allocation_reproduction) * net_mass_production_dt;
}

template <class S>
S TF24_Strategy_<S>::mass_total_dt(S fraction_allocation_reproduction,
                                     S net_mass_production_dt,
                                     S mass_heartwood_dt) const {
  return mass_live_dt(fraction_allocation_reproduction, net_mass_production_dt) +
    mass_heartwood_dt;
}

// TODO: Do we not track root mass change?
template <class S>
S TF24_Strategy_<S>::mass_above_ground_dt(S area_leaf,
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
S TF24_Strategy_<S>::mass_heartwood_dt(S mass_sapwood) const {
  return turnover_sapwood(mass_sapwood);
}


template <class S>
S TF24_Strategy_<S>::mass_live_given_height(S height) const {
  S area_leaf_ = area_leaf(height);
  return mass_leaf(area_leaf_) +
         mass_bark(area_bark(area_leaf_), height) +
         mass_sapwood(area_sapwood(area_leaf_), height) +
         mass_root(area_leaf_);
}

template <class S>
S TF24_Strategy_<S>::height_given_mass_leaf(S mass_leaf) const {
  return pars.a_l1 * pow(mass_leaf / pars.lma, pars.a_l2);
}

template <class S>
S TF24_Strategy_<S>::mortality_dt(S productivity_area,
                              S cumulative_mortality) const {

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

template <class S>
S TF24_Strategy_<S>::mortality_growth_independent_dt() const {
  return pars.d_I;
}

template <class S>
S TF24_Strategy_<S>::mortality_growth_dependent_dt(S productivity_area) const {
  return pars.a_dG1 * exp(-pars.a_dG2 * productivity_area);
}

// [eqn 20] Survival of seedlings during establishment
template <class S>
S TF24_Strategy_<S>::establishment_probability(const environment_type& environment) {

  S decay_over_time = exp(-pars.recruitment_decay * environment.time);

  const S net_mass_production_dt_ =
    net_mass_production_dt(environment, height_0, area_leaf_0, 1.0 / height_0);
  if (net_mass_production_dt_ > 0) {
    const S tmp = pars.a_d0 * area_leaf_0 / net_mass_production_dt_;
    return 1.0 / (tmp * tmp + 1.0) * decay_over_time;
  } else {
    return 0.0;
  }
}

template <class S>
S TF24_Strategy_<S>::compute_competition(double z, S height) const {
  return compute_competition(z, area_leaf(height), 1.0 / height);
}

// Ratio-first hot-path overload (see header): receives the cached
// competition_effect (= area_leaf(height)) and height_inverse (= 1/height), so the
// per-call area_leaf() evaluation and z/height division are hoisted out of the
// inner competition loop. Reproduces pars.k_I * area_leaf(height) * Q(z, height, pars.eta).
template <class S>
S TF24_Strategy_<S>::compute_competition(double z, S area_leaf_,
                                          S height_inverse) const {
  const S u = z * height_inverse;  // z / height
  if (u > 1.0) {
    return 0.0;
  }
  const S tmp = 1.0 - pow(u, pars.eta);
  return pars.k_I * area_leaf_ * tmp * tmp;
}

// [eqn  9] Probability density of leaf area at height `z` (double: leaf weighting)
template <class S>
double TF24_Strategy_<S>::q(double z, double height) const {
  const double tmp = pow(z / height, to_passive(pars.eta));
  return 2 * to_passive(pars.eta) * (1 - tmp) * tmp / z;
}

// [eqn 10] ... Fraction of leaf area above height 'z' for an
//              individual of height 'height' (double: root-depth distribution)
template <class S>
double TF24_Strategy_<S>::Q(double z, double height, double eta_x) const {
  if (z > height) {
    return 0.0;
  }
  const double tmp = 1.0-pow(z / height, eta_x);
  return tmp * tmp;
}

// (inverse of [eqn 10]; return the height above which fraction 'x' of
// the leaf mass would be found).
template <class S>
double TF24_Strategy_<S>::Qp(double x, double height) const { // x in [0,1], unchecked.
  return pow(1 - sqrt(x), (1/to_passive(pars.eta))) * height;
}

// The aim is to find a plant height that gives the correct seed mass. Kept a
// double root-find (the birth-height IFT seam is deferred); the narrowings
// strip every AD layer, so they are no-ops at S = double.
template <class S>
double TF24_Strategy_<S>::height_seed(void) const {

  // Note, these are not entirely correct bounds. Ideally we would use height
  // given *total* mass, not leaf mass, but that is difficult to calculate.
  const double
    h0 = to_passive(height_given_mass_leaf(std::numeric_limits<double>::min())),
    h1 = to_passive(height_given_mass_leaf(pars.omega));

  const double tol = this->control.offspring_production_tol;
  const size_t max_iterations = this->control.offspring_production_iterations;

  auto target = [&] (double x) mutable -> double {
    return to_passive(mass_live_given_height(S(x))) - to_passive(pars.omega);
  };

  return util::uniroot(target, h0, h1, tol, max_iterations);
}

template <class S>
void TF24_Strategy_<S>::prepare_strategy() {

  // Set up the function_integrator
  function_integrator = quadrature::QK(
      // Gauss-Kronrod quadrature integeration rule (see qkrules)
      this->control.function_integration_rule);

  // Resolve the crown shading model once. The empty Control default maps to
  // TF24's own default (mean-light, its established behaviour); PPA is an
  // FF16-only stepped-light model and is rejected here.
  shading_model_ =
    shading_model_from_string(this->control.shading_model, ShadingModel::MeanLight);
  // PPA and the flat-top-box variants are FF16-only (they reshape the FF16
  // competition / light profile, which TF24 does not use).
  if (shading_model_ == ShadingModel::PPA ||
      shading_model_ == ShadingModel::FlatTopBox ||
      shading_model_ == ShadingModel::FlatTopSoftBox) {
    throw std::invalid_argument(
      "shading_model '" + this->control.shading_model +
      "' is not supported for the TF24 strategy");
  }

  // NOTE: this pre-computes something to save a very small amount of time
  eta_c = 1 - 2/(1 + pars.eta) + 1/(1 + 2*pars.eta);
  // NOTE: Also pre-computing, though less trivial. Birth height is a double
  // root-find lifted to S; height_0 keeps the value, initial_height_ additionally
  // carries the birth-height IFT derivative (§7.2, mirrors FF16).
  height_0 = height_seed();
  area_leaf_0 = area_leaf(height_0);
  initial_height_ = lift_birth_height(to_passive(height_0));

  if (this->is_variable_birth_rate) {
    this->extrinsic_drivers.set_variable("birth_rate", this->birth_rate_x, this->birth_rate_y);
  } else {
    this->extrinsic_drivers.set_constant("birth_rate", this->birth_rate_y[0]);
  }
  leaf = Leaf(to_passive(pars.vcmax_25), to_passive(pars.c), to_passive(pars.b),
              to_passive(pars.psi_crit), root_c, root_b,
              root_psi_crit, to_passive(pars.beta2), to_passive(pars.jmax_25),
              to_passive(pars.a), to_passive(pars.curv_fact_elec_trans),
              to_passive(pars.curv_fact_colim),
              this->control.GSS_tol_abs, this->control.vulnerability_curve_ncontrol,
              this->control.ci_abs_tol, this->control.ci_niter,
              to_passive(pars.g1_TF24), beta_R_H, beta_R_V);
}

// Explicit instantiations. TF24's dg/dh is delivered by supplied_derivative
// (not the nested forward-over-reverse path), so the instantiation set is
// closed: the double production path, and the active reverse scalar for AD.
template class TF24_Strategy_<double>;
template class TF24_Strategy_<xad::adj<double>::active_type>;

}
