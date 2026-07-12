// Built from  src/ff16_strategy.cpp on Mon Feb 12 09:52:27 2024 using the scaffolder, from the strategy:  FF16
#include <plant/models/tf24_strategy.h>
#include <plant/individual.h>
#include <plant/patch.h>
#include <odelia/ode_solver.hpp>
#include <XAD/XAD.hpp>
#include <algorithm>
#include <cmath>

namespace plant {

// --- Hard-coded root-distribution constants (review #9) ---------------------
static const double root_mass_carbon_scale = 83.26 * 0.5;
// rooting depth cap (m), i.e. the depth of the soil column.
static const double rooting_depth_max = 1.5;

// NOTE (review #9): the per-second -> annual factor 60*60*12*365 recurs below and
// is deliberately left inline to preserve bit-identical results (collapsing the
// integer product changes rounding, which the adaptive ODE amplifies).

template <class S>
TF24_Strategy_<S>::TF24_Strategy_() {
  collect_all_auxiliary = false;
  refresh_indices();
  name = "TF24";
}

// not sure 'average' is the right term here..
template <class S>
double TF24_Strategy_<S>::compute_average_light_environment(
    double z, double height, const TF24_Environment_<S> &environment) {
  // NOTE: the light environment is clamped to a small positive floor (1e-4).
  return std::max(ad_value(environment.get_environment_at_height(z)), 0.0001) *
         q_double(z, height);
}

// assumes the leaf has been optimised for the operating point
template <class S>
double TF24_Strategy_<S>::evapotranspiration_dt(double area_leaf_, int soil_layer) {
  return leaf.soil_consumption_[soil_layer] * area_leaf_;
}

// Double replicas of the crown-shape kernels used by the off-tape leaf block.
template <class S>
double TF24_Strategy_<S>::q_double(double z, double height) const {
  const double eta = ad_value(pars.eta);
  const double tmp = pow(z / height, eta);
  return 2 * eta * (1 - tmp) * tmp / z;
}

template <class S>
double TF24_Strategy_<S>::Q_double(double z, double height, double eta_x) const {
  if (z > height) {
    return 0.0;
  }
  const double tmp = 1.0 - pow(z / height, eta_x);
  return tmp * tmp;
}

template <class S>
double TF24_Strategy_<S>::area_leaf_double(double height) const {
  return pow(height / ad_value(pars.a_l1), 1.0 / ad_value(pars.a_l2));
}

template <class S>
void TF24_Strategy_<S>::refresh_indices () {
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

  aux_idx_competition_effect    = aux_index.at("competition_effect");
  aux_idx_height_inverse        = aux_index.at("height_inverse");
  aux_idx_net_mass_production_dt = aux_index.at("net_mass_production_dt");
  aux_idx_root_mass             = aux_index.at("root_mass");
  aux_idx_opt_psi_stem          = aux_index.at("opt_psi_stem");
  aux_idx_opt_root_psi          = aux_index.at("opt_root_psi");
  aux_idx_transpiration         = aux_index.at("transpiration");
  aux_idx_E_up                  = aux_index.at("E_up_");
  aux_idx_profit                = aux_index.at("profit");
  aux_idx_stom_cond_CO2         = aux_index.at("stom_cond_CO2");
  aux_idx_area_sapwood = aux_index.count("area_sapwood") ? aux_index.at("area_sapwood") : -1;
  state_idx_area_heartwood      = state_index.at("area_heartwood");
  state_idx_mass_heartwood      = state_index.at("mass_heartwood");
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

// one-shot update of the scm variables
template <class S>
void TF24_Strategy_<S>::compute_rates(const TF24_Environment_<S>& environment,
                                      Internals_<S>& vars) {
  S height = vars.state(HEIGHT_INDEX);
  S area_leaf_ = vars.aux(aux_idx_competition_effect);

  const S net_mass_production_dt_ =
    net_mass_production_dt(environment, height, area_leaf_,
                           vars.aux(aux_idx_height_inverse));

  vars.set_aux(aux_idx_net_mass_production_dt, net_mass_production_dt_);
  vars.set_aux(aux_idx_root_mass, mass_root(area_leaf_));
  vars.set_aux(aux_idx_opt_psi_stem, S(leaf.opt_psi_stem_));
  vars.set_aux(aux_idx_opt_root_psi, S(leaf.root_collar_psi_));
  vars.set_aux(aux_idx_transpiration, S(leaf.transpiration_));
  vars.set_aux(aux_idx_E_up, S(leaf.E_up_));
  vars.set_aux(aux_idx_profit, S(leaf.profit_));
  vars.set_aux(aux_idx_stom_cond_CO2, S(leaf.stom_cond_CO2_));

  // consumption rates: canopy-level yearly evapotranspiration per soil layer.
  const double area_leaf_d = ad_value(area_leaf_);
  int soil_number_of_depths_ = environment.get_soil_number_of_depths();
  for (int i = 0; i < soil_number_of_depths_; i++) {
    vars.set_consumption_rate(
        i, S(evapotranspiration_dt(area_leaf_d, i) * 60 * 60 * 12 * 365 / 1000 *
             kg_per_mol_h2o));
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

    if (collect_all_auxiliary) {
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

// [eqn 12] Gross annual CO2 assimilation (double diagnostic; not on the rate path)
template <class S>
double TF24_Strategy_<S>::assimilation(const TF24_Environment_<S>& environment,
                                       double height, double area_leaf) {
  std::function<double(double)> f = [&](double z) -> double {
    return assimilation_leaf(ad_value(environment.get_environment_at_height(z))) *
           q_double(z, height);
  };
  const double A = function_integrator.integrate(f, 0.0, height);
  return area_leaf * A;
}

// Photosynthetic rate per leaf area; `x` is openness in [0, 1].
template <class S>
double TF24_Strategy_<S>::assimilation_leaf(double x) const {
  return ad_value(pars.a_p1) * x / (x + ad_value(pars.a_p2));
}

// [eqn 13] Total maintenance respiration
template <class S>
S TF24_Strategy_<S>::respiration(S mass_leaf, S mass_sapwood,
                                 S mass_bark, S mass_root) const {
  return respiration_leaf(mass_leaf) +
         respiration_bark(mass_bark) +
         respiration_sapwood(mass_sapwood) +
         respiration_root(mass_root);
}

template <class S>
S TF24_Strategy_<S>::respiration_leaf(S mass) const { return pars.r_l * mass; }
template <class S>
S TF24_Strategy_<S>::respiration_bark(S mass) const { return pars.r_b * mass; }
template <class S>
S TF24_Strategy_<S>::respiration_sapwood(S mass) const { return pars.r_s * mass; }
template <class S>
S TF24_Strategy_<S>::respiration_root(S mass) const { return pars.r_r * mass; }

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
S TF24_Strategy_<S>::turnover_leaf(S mass) const { return pars.k_l * mass; }
template <class S>
S TF24_Strategy_<S>::turnover_bark(S mass) const { return pars.k_b * mass; }
template <class S>
S TF24_Strategy_<S>::turnover_sapwood(S mass) const { return pars.k_s * mass; }
template <class S>
S TF24_Strategy_<S>::turnover_root(S mass) const { return pars.k_r * mass; }

// [eqn 15] Net production
template <class S>
S TF24_Strategy_<S>::net_mass_production_dt_A(S assimilation, S respiration,
                                S turnover) const {
  return pars.a_bio * pars.a_y * (assimilation - respiration) - turnover;
}

// One shot calculation of net_mass_production_dt.
// The size cascade (respiration/turnover) runs at the active scalar and flows
// through the reverse tape natively. The leaf optimum is a plant-local, off-tape
// double solve; its per-area profit re-enters the tape linearised about the operating
// point, carrying d(profit)/d(vcmax_25) (envelope + IFT) and d(profit)/d(height) (the
// size feedback, central difference of the off-tape solve), so the reverse sweep
// traverses density->optimum->trait without recording the leaf solve.
template <class S>
S TF24_Strategy_<S>::net_mass_production_dt(const TF24_Environment_<S>& environment,
                                S height, S area_leaf_, S height_inverse) {
  (void)height_inverse;  // TF24 works in height directly
  const S mass_leaf_    = mass_leaf(area_leaf_);
  const S area_sapwood_ = area_sapwood(area_leaf_);
  const S mass_sapwood_ = mass_sapwood(area_sapwood_, height);
  const S area_bark_    = area_bark(area_leaf_);
  const S mass_bark_    = mass_bark(area_bark_, height);
  const S mass_root_    = mass_root(area_leaf_);

  const double height_d = ad_value(height);
  const double area_leaf_d = ad_value(area_leaf_);
  const double profit0 = solve_leaf_at_size(environment, height_d, area_leaf_d);

  S profit_active;
  if constexpr (std::is_same_v<S, double>) {
    profit_active = profit0;  // resident numerics, bit-identical
  } else {
    // Active scalar (the SCM trait gradient runs the reverse twin): re-attach the
    // off-tape leaf's per-area profit to the tape linearised about the operating
    // point. profit0 is the frozen value; each partial * (active_input -
    // value(active_input)) is value-zero but carries that input's derivative, so the
    // reverse sweep traverses density->optimum->trait on the one tape using only the
    // model's own active arithmetic -- no off-tape leaf is registered as a tape input.
    const double opt = -leaf.root_collar_psi_;  // optimised collar potential (+ve magnitude)
    const double dP_dvcmax25 = leaf.dprofit_dvcmax25(opt);
    profit_active = S(profit0) + dP_dvcmax25 * (pars.vcmax_25 - ad_value(pars.vcmax_25));
    // Size-feedback channel: the per-area profit also responds to the (active) plant
    // height through the hydraulic supply (~1/height), sapwood volume (~height), crown
    // radiation and leaf area. Add its central-difference partial (a further off-tape
    // solve pair), linearised in the active height.
    if (add_size_feedback_ && height_d != 0.0) {
      // Balanced step: large enough that the leaf-solver residual noise (~ci_abs_tol)
      // does not swamp the slope, small enough that the O(h^2) truncation stays small.
      const double dh = std::max(std::abs(height_d), 0.5) * 3e-4;
      const double p_plus =
          solve_leaf_at_size(environment, height_d + dh, area_leaf_double(height_d + dh));
      const double p_minus =
          solve_leaf_at_size(environment, height_d - dh, area_leaf_double(height_d - dh));
      const double dP_dheight = (p_plus - p_minus) / (2.0 * dh);
      solve_leaf_at_size(environment, height_d, area_leaf_d);  // restore leaf.* outputs
      profit_active += dP_dheight * (height - ad_value(height));
    }
  }

  const S assimilation_ = profit_active * area_leaf_ * 60 * 60 * 12 * 365 / 1e6;
  const S respiration_ =
    respiration(mass_leaf_, mass_sapwood_, mass_bark_, mass_root_);
  const S turnover_ =
    turnover(mass_leaf_, mass_bark_, mass_sapwood_, mass_root_);
  return net_mass_production_dt_A(assimilation_, respiration_, turnover_);
}

// Off-tape (double) leaf optimum at a given size and the frozen canopy the
// environment holds. Mirrors develop's net_mass_production_dt body up to the leaf
// solve; leaves leaf.* at the operating point and returns per-leaf-area profit.
template <class S>
double TF24_Strategy_<S>::solve_leaf_at_size(const TF24_Environment_<S>& environment,
                                             double height, double area_leaf_) {
  const double mass_root_ = ad_value(pars.a_r1) * area_leaf_;  // = mass_root(area_leaf_)

  const int soil_number_of_depths_ = environment.get_soil_number_of_depths();
  const std::vector<double>& soil_depths_ = environment.z;
  const std::vector<double>& psi_soil = environment.get_soil_water_potential_state();

  const double eta_c_d = ad_value(eta_c);
  // leaf specific max hydraulic conductance (kg m^-2 LA s^-1 MPa^-1)
  const double leaf_specific_conductance_max =
      ad_value(pars.K_s) * ad_value(pars.theta) / (height * eta_c_d);
  // sapwood volume per leaf area
  const double sapwood_volume_per_leaf_area = ad_value(pars.theta) * (height * eta_c_d);

  // ROOT MASS DISTRIBUTION ACROSS SOIL LAYERS (Q with the root shape exponent).
  mass_root_prop_.assign(soil_number_of_depths_, 0.0);
  const double rooting_depth = std::min(height, rooting_depth_max);
  const double root_mass_scale = root_mass_carbon_scale * mass_root_;
  const double root_eta = ad_value(pars.root_depth_shape_eta);
  double prev_q = 1.0;
  for (int a = 0; a < soil_number_of_depths_; ++a) {
    if (prev_q == 0) {
      break;
    }
    const double q_ = Q_double(soil_depths_[a], rooting_depth, root_eta);
    mass_root_prop_[a] = root_mass_scale * (prev_q - q_);
    prev_q = q_;
  }

  leaf.z_soil_mid_ = environment.get_soil_mid_depths();
  leaf.use_precomputed_z_soil_mid_ = true;

  auto optimise_at = [&](double radiation) {
    leaf.set_physiology(area_leaf_, mass_root_prop_, ad_value(pars.rho),
                        ad_value(pars.a_bio), radiation, psi_soil, soil_depths_,
                        leaf_specific_conductance_max, environment.get_atm_vpd(),
                        environment.get_ca(), sapwood_volume_per_leaf_area,
                        environment.get_leaf_temp(), environment.get_atm_o2_kpa(),
                        environment.get_atm_kpa());
    solve_leaf();
  };

  const double PPFD = environment.get_PPFD();
  const double k_I_d = ad_value(pars.k_I);
  auto radiation_at = [&](double light) {
    return k_I_d * std::max(light, 0.0001) * PPFD;
  };

  if (shading_model_ == ShadingModel::CrownCentre) {
    optimise_at(radiation_at(
        ad_value(environment.get_environment_at_height(height * eta_c_d))));
  } else if (shading_model_ == ShadingModel::MeanLight) {
    auto f = [&](double x) -> double {
      return compute_average_light_environment(x, height, environment);
    };
    optimise_at(radiation_at(function_integrator.integrate(f, 0.0, height)));
  } else { // DeepCrown
    const std::vector<double> nodes =
      function_integrator.integrate_vector_x(0.0, height);
    const size_t nn = nodes.size();
    std::vector<double> profit_y(nn), trans_y(nn), eup_y(nn), psi_y(nn),
      root_psi_y(nn), gco2_y(nn);
    std::vector<std::vector<double>> soil_y(
      soil_number_of_depths_, std::vector<double>(nn));
    for (size_t i = 0; i < nn; ++i) {
      const double qi = q_double(nodes[i], height);
      optimise_at(radiation_at(ad_value(environment.get_environment_at_height(nodes[i]))));
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
    leaf.profit_          = function_integrator.integrate_vector(profit_y, 0.0, height);
    leaf.transpiration_   = function_integrator.integrate_vector(trans_y, 0.0, height);
    leaf.E_up_            = function_integrator.integrate_vector(eup_y, 0.0, height);
    leaf.opt_psi_stem_    = function_integrator.integrate_vector(psi_y, 0.0, height);
    leaf.root_collar_psi_ = function_integrator.integrate_vector(root_psi_y, 0.0, height);
    leaf.stom_cond_CO2_   = function_integrator.integrate_vector(gco2_y, 0.0, height);
    for (int a = 0; a < soil_number_of_depths_; ++a) {
      leaf.soil_consumption_[a] =
        function_integrator.integrate_vector(soil_y[a], 0.0, height);
    }
  }
  return leaf.profit_;
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

template <class S>
S TF24_Strategy_<S>::dmass_leaf_darea_leaf(S /* area_leaf */) const {
  return pars.lma;
}

template <class S>
S TF24_Strategy_<S>::dmass_sapwood_darea_leaf(S area_leaf) const {
  return pars.rho * eta_c * pars.a_l1 * pars.theta * (pars.a_l2 + 1.0) * pow(area_leaf, pars.a_l2);
}

template <class S>
S TF24_Strategy_<S>::dmass_bark_darea_leaf(S area_leaf) const {
  return pars.a_b1 * dmass_sapwood_darea_leaf(area_leaf);
}

template <class S>
S TF24_Strategy_<S>::dmass_root_darea_leaf(S /* area_leaf */) const {
  return pars.a_r1;
}

template <class S>
S TF24_Strategy_<S>::ddiameter_stem_darea_stem(S area_stem) const {
  return pow(M_PI * area_stem, -0.5);
}

template <class S>
S TF24_Strategy_<S>::area_sapwood_dt(S area_leaf_dt) const {
  return area_leaf_dt * pars.theta;
}

template <class S>
S TF24_Strategy_<S>::area_heartwood_dt(S area_leaf) const {
  return pars.k_s * area_sapwood(area_leaf);
}

template <class S>
S TF24_Strategy_<S>::area_bark_dt(S area_leaf_dt) const {
  return pars.a_b1 * area_leaf_dt * pars.theta;
}

template <class S>
S TF24_Strategy_<S>::area_stem_dt(S area_leaf, S area_leaf_dt) const {
  return area_sapwood_dt(area_leaf_dt) +
    area_bark_dt(area_leaf_dt) +
    area_heartwood_dt(area_leaf);
}

template <class S>
S TF24_Strategy_<S>::diameter_stem_dt(S area_stem, S area_stem_dt) const {
  return ddiameter_stem_darea_stem(area_stem) * area_stem_dt;
}

template <class S>
S TF24_Strategy_<S>::mass_root_dt(S area_leaf, S area_leaf_dt) const {
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
  if (util::is_finite(ad_value(cumulative_mortality))) {
    return
      mortality_growth_independent_dt() +
      mortality_growth_dependent_dt(productivity_area);
 } else {
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
S TF24_Strategy_<S>::establishment_probability(const TF24_Environment_<S>& environment) {
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

// Ratio-first hot-path overload; reproduces pars.k_I * area_leaf(height) * Q(z, height, pars.eta).
template <class S>
S TF24_Strategy_<S>::compute_competition(double z, S area_leaf_,
                                          S height_inverse) const {
  const S u = z * height_inverse;  // z / height
  if (ad_value(u) > 1.0) {
    return S(0.0);
  }
  const S tmp = 1.0 - pow(u, pars.eta);
  return pars.k_I * area_leaf_ * tmp * tmp;
}

// [eqn  9] Probability density of leaf area at height `z`
template <class S>
S TF24_Strategy_<S>::q(S z, S height) const {
  const S tmp = pow(z / height, pars.eta);
  return 2 * pars.eta * (1 - tmp) * tmp / z;
}

// [eqn 10] Fraction of leaf area above height 'z'
template <class S>
S TF24_Strategy_<S>::Q(S z, S height, S eta_x) const {
  if (ad_value(z) > ad_value(height)) {
    return S(0.0);
  }
  const S tmp = 1.0-pow(z / height, eta_x);
  return tmp * tmp;
}

// (inverse of [eqn 10])
template <class S>
S TF24_Strategy_<S>::Qp(S x, S height) const { // x in [0,1], unchecked.
  using std::sqrt;
  return pow(1 - sqrt(x), (1/pars.eta)) * height;
}

// The aim is to find a plant height that gives the correct seed mass.
template <class S>
S TF24_Strategy_<S>::height_seed(void) const {
  const double
    h0 = ad_value(height_given_mass_leaf(std::numeric_limits<double>::min())),
    h1 = ad_value(height_given_mass_leaf(pars.omega));

  const double tol = control.offspring_production_tol;
  const size_t max_iterations = control.offspring_production_iterations;

  auto target = [&] (double x) mutable -> double {
    return ad_value(mass_live_given_height(x) - pars.omega);
  };

  const double h_root = util::uniroot(target, h0, h1, tol, max_iterations);

  if constexpr (std::is_same_v<S, double>) {
    return h_root;
  } else {
    // Reattach the trait derivative the double solve discards, by the IFT: the
    // root h*(theta) solves g(h*, theta) = mass_live_given_height(h*) - omega = 0.
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
void TF24_Strategy_<S>::prepare_strategy() {
  function_integrator = quadrature::QK(control.function_integration_rule);

  // Resolve the crown shading model once. The empty Control default maps to
  // TF24's own default (mean-light); PPA / flat-top-box are FF16-only and rejected.
  shading_model_ =
    shading_model_from_string(control.shading_model, ShadingModel::MeanLight);
  if (shading_model_ == ShadingModel::PPA ||
      shading_model_ == ShadingModel::FlatTopBox ||
      shading_model_ == ShadingModel::FlatTopSoftBox) {
    throw std::invalid_argument(
      "shading_model '" + control.shading_model +
      "' is not supported for the TF24 strategy");
  }

  eta_c = 1 - 2/(1 + pars.eta) + 1/(1 + 2*pars.eta);
  height_0 = height_seed();
  area_leaf_0 = area_leaf(height_0);

  if (is_variable_birth_rate) {
    extrinsic_drivers.set_variable("birth_rate", birth_rate_x, birth_rate_y);
  } else {
    extrinsic_drivers.set_constant("birth_rate", birth_rate_y[0]);
  }
  // The embedded leaf is plant-local double: seed it from the passive trait values.
  leaf = Leaf(ad_value(pars.vcmax_25), ad_value(pars.c), ad_value(pars.b),
              ad_value(pars.psi_crit), root_c, root_b, root_psi_crit,
              ad_value(pars.beta2), ad_value(pars.jmax_25), ad_value(pars.a),
              ad_value(pars.curv_fact_elec_trans), ad_value(pars.curv_fact_colim),
              control.GSS_tol_abs, control.vulnerability_curve_ncontrol,
              control.ci_abs_tol, control.ci_niter, ad_value(pars.g1_TF24),
              beta_R_H, beta_R_V);
}

TF24_Strategy::ptr make_strategy_ptr(TF24_Strategy s) {
  s.prepare_strategy();
  return std::make_shared<TF24_Strategy>(s);
}

// The resident numerics cross the R boundary at double; the reverse-mode active
// scalar is instantiated here too so plant.so carries the differentiable TF24
// census path (the leaf edge injects the off-tape optimum's sensitivity).
template class TF24_Strategy_<double>;

using ad_reverse =
    odelia::ode::Solver<Patch<TF24_Strategy_<double>, TF24_Environment>>::active_scalar;
template class TF24_Strategy_<ad_reverse>;
}
