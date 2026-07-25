// Built from  src/ff16_strategy.cpp on Mon Feb 12 09:52:27 2024 using the scaffolder, from the strategy:  FF16
#include <plant/models/tf24_strategy.h>
#include <odelia/implicit_node.hpp> // odelia::implicit_value (leaf roots), lift_root (birth-height IFT lift)
#include <odelia/ode_util.hpp>      // odelia::util::graft_value (leaf output value-graft)
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
  // Tier-B: prefer the active per-layer uptake (carries theta/soil-psi
  // sensitivity assembled by assemble_leaf_from); fall back to the leaf value
  // when it has not been populated (e.g. the double production path before
  // net_mass_production_dt, or shading models without the single-solve assembly).
  const S uptake =
      (soil_layer < static_cast<int>(soil_consumption_active_.size()))
          ? soil_consumption_active_[soil_layer]
          : S(leaf.soil_consumption_[soil_layer]);
  return uptake * area_leaf_;
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

// Census weights (mirror FF16::census_*), reading TF24's runtime aux/state
// slots. area_leaf is the competition-effect aux; heartwood mass/area are ODE
// states.
template <class S>
S TF24_Strategy_<S>::census_leaf_area(const Internals_<S>& vars) const {
  return vars.aux(aux_idx_competition_effect);  // = area_leaf(height)
}
template <class S>
S TF24_Strategy_<S>::census_mass(const Internals_<S>& vars) const {
  const S height        = vars.state(HEIGHT_INDEX);
  const S area_leaf_     = vars.aux(aux_idx_competition_effect);
  const S mass_leaf_     = mass_leaf(area_leaf_);
  const S mass_sapwood_  = mass_sapwood(area_sapwood(area_leaf_), height);
  const S mass_bark_     = mass_bark(area_bark(area_leaf_), height);
  const S mass_heartwood_ = vars.state(state_idx_mass_heartwood);
  return mass_above_ground(mass_leaf_, mass_bark_, mass_sapwood_, mass_heartwood_);
}
template <class S>
S TF24_Strategy_<S>::census_basal_area(const Internals_<S>& vars) const {
  const S area_leaf_      = vars.aux(aux_idx_competition_effect);
  const S area_heartwood_ = vars.state(state_idx_area_heartwood);
  return area_stem(area_bark(area_leaf_), area_sapwood(area_leaf_), area_heartwood_);
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
  // Leaf hydraulic sub-model. The Leaf is solved in DOUBLE (its golden-section /
  // bracketing solvers have double-only signatures). The S mass cascade re-enters
  // via the area_leaf_ factor in assimilation_ below; every seeded-parameter /
  // soil-state / self-shading sensitivity of the Leaf outputs re-enters through
  // assemble_leaf_from (the closed-form output map + implicit_value leaf roots),
  // just below. On the double path this is bit-identical to the pre-template code
  // (to_passive is the identity).
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
  // single_solve models run one leaf optimisation at a scalar openness; the
  // active leaf seam below re-evaluates the leaf at that same openness and
  // differentiates through it (the resident self-shading light channel).
  //   light_openness_double : the openness fed to radiation_at (double, forward)
  //   light_active          : the same openness as an active scalar whose slot
  //                           connects to the resident light field (self-shading)
  double radiation_used = NA_REAL;
  double light_openness_double = NA_REAL;
  S light_active = S(0.0);
  bool single_solve = true;
  if (shading_model_ == ShadingModel::CrownCentre) {
    light_active = environment.get_environment_at_height(height * eta_c);
    light_openness_double = to_passive(light_active);
    radiation_used = radiation_at(light_openness_double);
    optimise_at(radiation_used);
  } else if (shading_model_ == ShadingModel::MeanLight) {
    // Leaf-area-weighted mean canopy openness = integral of (light * q) over the
    // crown (q integrates to one).
    auto f = [&](double x) -> double {
      return compute_average_light_environment(x, height_d, environment);
    };
    light_openness_double = function_integrator.integrate(f, 0.0, height_d);
    // The same mean openness as an active scalar for the seam's light channel.
    // Fixed integration bounds (double height_d cast to S), so no plant-height
    // derivative is entangled; the active light reads carry the self-shading
    // derivative and the value matches light_openness_double (same quadrature).
    // Active path only -- the double forward never reads light_active.
    if constexpr (!std::is_same_v<S, double>) {
      auto f_active = [&](S z) -> S {
        const S er = environment.get_environment_at_height(z);
        const S lit = (er > 0.0001) ? er : S(0.0001);
        return lit * q(to_passive(z), height_d);
      };
      light_active = function_integrator.integrate(f_active, S(0.0), S(height_d));
    }
    radiation_used = radiation_at(light_openness_double);
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
  // converts to canopy area, then years, then mols.
  //
  // Reduce the leaf submodel to net photosynthetic profit per leaf area. The Leaf
  // is solved in double (its golden-section/bracketing solvers have double-only
  // signatures -- an active scalar cannot enter the iteration). assemble_leaf_from
  // re-expresses the converged operating point as active scalars so every
  // seeded-parameter / soil-state / self-shading sensitivity of profit + per-layer
  // soil uptake reaches the run tape through the closed-form output map and the
  // implicit_value nodes at the leaf's roots (ci, psi_stem, p*). On the double path
  // it collapses to leaf.profit_ (the whole branch is compiled out).
  S profit_s = S(leaf.profit_);
  if constexpr (!std::is_same_v<S, double>) {
    if (!single_solve) {
      util::stop("TF24 active leaf gradient: the deep-crown assembly is not yet "
                 "implemented; use shading_model 'mean-light' or 'crown-centre'.");
    }
    const int nsoil = environment.get_soil_number_of_depths();
    // Leaf signed-potential convention: soil_uptake reads the negated soil state.
    std::vector<S> psi_soil_signed(psi_soil_S.size());
    for (std::size_t i = 0; i < psi_soil_S.size(); ++i)
      psi_soil_signed[i] = -psi_soil_S[i];
    std::vector<S> cons;
    profit_s = assemble_leaf_from(pars, height, light_active,
                                  light_openness_double, psi_soil_signed,
                                  environment, cons);
    const int mlayer = static_cast<int>(cons.size());
    soil_consumption_active_.assign(nsoil, S(0.0));
    for (int L = 0; L < nsoil; ++L)
      soil_consumption_active_[L] =
          (L < mlayer) ? cons[L] : S(to_passive(leaf.soil_consumption_[L]));
  }
  const S assimilation_ = profit_s * area_leaf_* 60*60*12*365/1e6;
  // const double assimilation_ = assimilation(environment, height, area_leaf_);
  const S respiration_ =
    respiration(mass_leaf_, mass_sapwood_, mass_bark_, mass_root_);
  const S turnover_ =
    turnover(mass_leaf_, mass_bark_, mass_sapwood_, mass_root_);
  return net_mass_production_dt_A(assimilation_, respiration_, turnover_);
}

template <class S>
S TF24_Strategy_<S>::assemble_leaf_from(const TF24_Pars_<S>& p, S height,
                                        S light_active,
                                        double light_openness_double,
                                        const std::vector<S>& psi_soil_signed,
                                        const environment_type& environment,
                                        std::vector<S>& cons_out) {
  using std::pow;
  // Value from the converged double leaf, derivative from the assembled active
  // expression. odelia::graft_value owns this idiom precisely because writing it
  // inline is a memory-safety trap (a deduced-return-type lambda returns a
  // dangling expression template); see the note on the primitive.
  auto anchor = [](double v, const S& x) -> S { return odelia::util::graft_value<S>(v, x); };

  // Geometry active in height + the seeded allometry/root traits.
  const S eta_c_local = 1.0 - 2.0 / (1.0 + p.eta) + 1.0 / (1.0 + 2.0 * p.eta);
  const S area_leaf_local = pow(height / p.a_l1, 1.0 / p.a_l2);
  const S mass_root_local = p.a_r1 * area_leaf_local;

  // Root-mass distribution across soil layers (mirrors net_mass_production_dt's
  // double loop); the resistances soil_uptake divides by derive from it, so this
  // is the height/a_r1/a_l1/a_l2/root_depth_shape_eta channel into the uptake.
  const std::vector<double>& soil_depths_ = environment.z;
  const int nsoil = environment.get_soil_number_of_depths();
  const S rooting_depth =
      (to_passive(height) < rooting_depth_max) ? height : S(rooting_depth_max);
  const S root_mass_scale = root_mass_carbon_scale * mass_root_local;
  std::vector<S> mrp(nsoil, S(0.0));
  S prev_q(1.0);
  for (int a = 0; a < nsoil; ++a) {
    if (to_passive(prev_q) == 0.0) break;
    S qd(0.0);
    if (soil_depths_[a] <= to_passive(rooting_depth)) {
      const S tmp = S(1.0) - pow(S(soil_depths_[a]) / rooting_depth,
                                 p.root_depth_shape_eta);
      qd = tmp * tmp;
    }
    mrp[a] = root_mass_scale * (prev_q - qd);
    prev_q = qd;
  }

  const int mlayer = static_cast<int>(leaf.r_R_H_min.size());
  const double dz_sq = leaf.dz_ * leaf.dz_;
  std::vector<S> rH(mlayer), rVsum(mlayer);
  S vsum(0.0);
  for (int i = 0; i < mlayer; ++i) {
    const S rm = mrp[i];
    if (to_passive(rm) == 0.0) {
      rH[i] = S(0.0);
      rVsum[i] = vsum;
      continue;
    }
    rH[i] = leaf.beta_R_H / (rm * 2.0 / 3.0);
    vsum += leaf.beta_R_V * dz_sq / (rm / 3.0);
    rVsum[i] = vsum;
  }

  // Physiology from the active parameters, mirroring set_physiology. gstar_Pa,
  // km, ca, atm_* are passive (temperature/O2/environment only).
  const double leaf_temp = environment.get_leaf_temp();
  const double PPFD_env = environment.get_PPFD();
  const S vcmax = Leaf::peak_arrh_curve<S>(vcmax_ha, p.vcmax_25, leaf_temp,
                                                  vcmax_H_d, vcmax_d_S);
  const S jmax = Leaf::peak_arrh_curve<S>(jmax_ha, p.jmax_25, leaf_temp,
                                                 jmax_H_d, jmax_d_S);
  const S lit = (light_openness_double > 0.0001) ? light_active : S(0.0001);
  const S radiation = p.k_I * lit * PPFD_env;
  const S et = Leaf::electron_transport<S>(p.a, radiation, jmax,
                                                  p.curv_fact_elec_trans);
  const S R_d = vcmax * 0.015;
  const double gstar_Pa = leaf.gamma_ * umol_per_mol_to_Pa;
  const double km = leaf.km_;
  const double ca = environment.get_ca();
  const double atm_kpa = environment.get_atm_kpa();
  const double atm_vpd = environment.get_atm_vpd();
  const S kmax = p.K_s * p.theta / (height * eta_c_local);
  const S b = p.b, c = p.c, psi_crit = p.psi_crit;
  const S g1 = p.g1_TF24, beta2 = p.beta2, curv_colim = p.curv_fact_colim;

  // Converged double anchors.
  const double p_star_d = -leaf.root_collar_psi_;
  const double psi_stem_star_d = leaf.opt_psi_stem_;
  const double ci_star_d = leaf.ci_;
  const double root_crit_d = leaf.root_collar_psi_;  // signed; == -p_star_d

  // Assemble the outputs at a given active collar magnitude p_collar, anchoring
  // psi_stem/ci at the double optimum. Fills the per-layer uptake if asked.
  auto assemble = [&](S p_collar, std::vector<S>* out) -> S {
    std::vector<S> cons(mlayer);
    S eup = Leaf::soil_uptake<S>(psi_soil_signed, -p_collar, area_leaf_local,
                                        rH, rVsum, leaf.grav_head_z_, leaf.root_b,
                                        leaf.root_c, cons);
    S psi_stem =
        Leaf::psistem_node<S>(psi_stem_star_d, p_collar, eup, kmax, b, c);
    S transp = Leaf::transpiration<S>(psi_stem, p_collar, kmax, b, c);
    S gc = Leaf::stom_cond_CO2<S>(transp, atm_kpa, atm_vpd);
    S ci = Leaf::ci_node<S>(ci_star_d, vcmax, et, gstar_Pa, km, R_d,
                                   curv_colim, gc, ca, atm_kpa);
    S assim = Leaf::assim_colimited<S>(ci, vcmax, et, gstar_Pa, km, R_d,
                                              curv_colim);
    S cost = Leaf::hydraulic_cost_TF<S>(psi_stem, g1, beta2, b, c);
    if (out) *out = cons;
    return assim - cost;
  };

  // --------------------------------------------------------------------------
  // Stateless double anchors for the interior dp* central difference.
  //
  // The Leaf's own solvers (find_psi_stem_from_psi_root / psi_stem_to_ci) mutate
  // leaf scratch and run bracketing solvers behind std::function. Re-solving the
  // anchors here from the SAME closed forms at T=double keeps the difference a
  // pure function: no shared state to save and restore, and nothing between the
  // active expressions and the tape. Both inner problems are monotone, so a
  // bisection on the closed form is sufficient and needs no bracket heuristics.
  // --------------------------------------------------------------------------
  const double kmax_d  = to_passive(kmax);
  const double b_d     = to_passive(b);
  const double c_d     = to_passive(c);
  const double curv_d  = to_passive(curv_colim);
  const double vcmax_d = to_passive(vcmax);
  const double et_d    = to_passive(et);
  const double R_d_d   = to_passive(R_d);
  const double area_d  = to_passive(area_leaf_local);
  const double psi_crit_d = to_passive(psi_crit);
  std::vector<double> psi_soil_d(psi_soil_signed.size());
  for (std::size_t i = 0; i < psi_soil_signed.size(); ++i)
    psi_soil_d[i] = to_passive(psi_soil_signed[i]);
  std::vector<double> rH_d(mlayer), rV_d(mlayer);
  for (int i = 0; i < mlayer; ++i) {
    rH_d[i] = to_passive(rH[i]);
    rV_d[i] = to_passive(rVsum[i]);
  }

  // psi_stem from the supply balance transpiration(psi_stem, p) = E_up, and ci
  // from A(ci) = gc (ca - ci); both strictly monotone in their unknown.
  auto anchors_at = [&](double p_pass, double& psi_stem_out, double& ci_out) {
    std::vector<double> cons_d(mlayer);
    const double eup = Leaf::soil_uptake<double>(psi_soil_d, -p_pass, area_d, rH_d,
                                                 rV_d, leaf.grav_head_z_,
                                                 leaf.root_b, leaf.root_c, cons_d);
    double lo = p_pass, hi = std::max(psi_crit_d, p_pass + 1e-6);
    for (int k = 0; k < 100; ++k) {
      const double m = 0.5 * (lo + hi);
      if (Leaf::transpiration<double>(m, p_pass, kmax_d, b_d, c_d) < eup) lo = m;
      else hi = m;
    }
    psi_stem_out = 0.5 * (lo + hi);

    const double transp =
        Leaf::transpiration<double>(psi_stem_out, p_pass, kmax_d, b_d, c_d);
    const double gc_d = Leaf::stom_cond_CO2<double>(transp, atm_kpa, atm_vpd);
    const double inv_atm = 1.0 / (atm_kpa * kPa_to_Pa);
    double clo = gstar_Pa, chi = ca;
    for (int k = 0; k < 100; ++k) {
      const double m = 0.5 * (clo + chi);
      const double f = Leaf::assim_colimited<double>(m, vcmax_d, et_d, gstar_Pa, km,
                                                     R_d_d, curv_d) * umol_to_mol -
                       gc_d * (ca - m) * inv_atm;
      if (f < 0.0) clo = m; else chi = m;
    }
    ci_out = 0.5 * (clo + chi);
  };

  // The regime probe (E_column) still re-solves the double leaf off its optimum
  // and mutates leaf scratch (E_up_, root_collar_psi_, ci_, soil_consumption_).
  // Save the converged operating point and restore it before the output anchors
  // read it -- so the returned value stays the converged double leaf and no
  // scratch leaks into the post-call aux reads. This is the localized replacement
  // for the old whole-leaf snapshot in net_mass_production_dt.
  const double saved_E_up = leaf.E_up_;
  const double saved_rcp  = leaf.root_collar_psi_;
  const double saved_ci   = leaf.ci_;
  const std::vector<double> saved_cons = leaf.soil_consumption_;

  // The active pivot p*.
  S p_star;
  if (S* collar = seam_collar_psi_input()) {
    // Tracked-collar variant: operating collar is the double point (clamped), its
    // derivative rides the tracked ODE state.
    p_star = anchor(p_star_d, *collar);
  } else {
    // Regime from the converged point: p* sits at the stem-critical bound iff the
    // continuity residual there is zero.
    const double e_col =
        leaf.E_column(root_crit_d, leaf.psi_soil_inverted_, leaf.psi_crit);
    if (std::abs(e_col) < 1e-6) {
      // Bound regime: p* = -root_crit; root_crit solves the continuity residual
      // at psi_crit (no re-solve -- it is the converged signed collar).
      p_star = -odelia::implicit_value<S>(root_crit_d, [&](S x) -> S {
        std::vector<S> cons(mlayer);
        S eup = Leaf::soil_uptake<S>(psi_soil_signed, x, area_leaf_local, rH,
                                            rVsum, leaf.grav_head_z_, leaf.root_b,
                                            leaf.root_c, cons);
        S transp = Leaf::transpiration<S>(psi_crit, -x, kmax, b, c);
        return eup - transp;
      });
    } else {
      // Interior optimum: dp*/dstate = -P_ps/P_pp via a central difference of the
      // reduced profit (anchors re-solved off-tape per p), XAD supplying P_ps.
      const double eps = 1e-2 * (std::abs(p_star_d) + 1.0);
      auto profit_reduced = [&](S p_collar) -> S {
        const double p_pass = to_passive(p_collar);
        double psi_stem_star = 0.0, ci_star = 0.0;
        anchors_at(p_pass, psi_stem_star, ci_star);  // stateless re-solve
        std::vector<S> cons(mlayer);
        S eup = Leaf::soil_uptake<S>(psi_soil_signed, -p_collar,
                                            area_leaf_local, rH, rVsum,
                                            leaf.grav_head_z_, leaf.root_b,
                                            leaf.root_c, cons);
        S psi_stem =
            Leaf::psistem_node<S>(psi_stem_star, p_collar, eup, kmax, b, c);
        S transp = Leaf::transpiration<S>(psi_stem, p_collar, kmax, b, c);
        S gc = Leaf::stom_cond_CO2<S>(transp, atm_kpa, atm_vpd);
        S ci = Leaf::ci_node<S>(ci_star, vcmax, et, gstar_Pa, km, R_d,
                                       curv_colim, gc, ca, atm_kpa);
        S assim = Leaf::assim_colimited<S>(ci, vcmax, et, gstar_Pa, km, R_d,
                                                  curv_colim);
        S cost = Leaf::hydraulic_cost_TF<S>(psi_stem, g1, beta2, b, c);
        return assim - cost;
      };
      p_star = odelia::implicit_value<S>(p_star_d, [&](S p_collar) -> S {
        return (profit_reduced(p_collar + eps) - profit_reduced(p_collar - eps)) /
               (2.0 * eps);
      });
    }
  }

  // Restore the converged operating point clobbered by the regime re-solves above.
  leaf.E_up_ = saved_E_up;
  leaf.root_collar_psi_ = saved_rcp;
  leaf.ci_ = saved_ci;
  leaf.soil_consumption_ = saved_cons;

  std::vector<S> cons;
  const S profit_assembled = assemble(p_star, &cons);
  cons_out.assign(mlayer, S(0.0));
  for (int i = 0; i < mlayer; ++i)
    cons_out[i] = anchor(leaf.soil_consumption_[i], cons[i]);
  return anchor(leaf.profit_, profit_assembled);
}

// Base TF24: optimise the root-collar water potential from scratch each call.
template <class S>
void TF24_Strategy_<S>::solve_leaf() {
  leaf.find_root_collar_psi();
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
  // At the ground query z==0 (u==0), pow(u, eta) is identically 0 for every
  // eta>0, so tmp==1; but pow(active 0, active eta) yields a NaN *tangent* via
  // log(base)=log(0) in XAD's derivative formula (0*(-Inf)). Short-circuit to the
  // exact value (bit-identical on the double path) with the correct zero tangent.
  const S tmp = (z == 0.0) ? S(1.0) : S(1.0 - pow(u, pars.eta));
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
  // carries the birth-height IFT derivative (mirrors FF16).
  height_0 = height_seed();
  area_leaf_0 = area_leaf(height_0);
  initial_height_ = odelia::implicit_value<S>(
      to_passive(height_0), [this](S h) -> S { return mass_live_given_height(h) - pars.omega; });

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

// Explicit instantiations. The instantiation set is closed: the double
// production path and the active reverse scalar for AD (the leaf gradient
// re-enters through assemble_leaf_from's closed forms + implicit_value nodes).
template class TF24_Strategy_<double>;
template class TF24_Strategy_<xad::adj<double>::active_type>;

}
