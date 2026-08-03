#include <plant/leaf_model.h>
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <limits>
#include <string>
#include <exception>
#include <boost/math/special_functions/gamma.hpp>
#include <plant/models/tf24_environment.h>
#include <plant/collar_census.h>
#include <odelia/gradient.hpp>

namespace plant {

namespace {
// Templated replicas of the analytic profit algebra, so forward-mode AD gives
// their exact derivatives (used by Leaf::dprofit_droot_collar_psi). They mirror
// Leaf::assim_colimited and Leaf::hydraulic_cost_TF exactly.
// Every argument carries the scalar so a parameter is differentiable as well as
// ci or psi_stem; call sites name T explicitly and the passive arguments
// convert.
template <typename T>
T assim_colimited_ad(T ci, T vcmax, T et, T gstar_Pa, T km, T R_d, T curv) {
  T ar = vcmax * (ci - gstar_Pa) / (ci + km);
  T ae = et / 4.0 * (ci - gstar_Pa) / (ci + 2.0 * gstar_Pa);
  T s = ar + ae;
  return (s - sqrt(s * s - 4.0 * curv * ar * ae)) / (2.0 * curv) - R_d;
}
template <typename T>
T hydraulic_cost_ad(T psi_stem, T b, T c, T g1, T beta2) {
  return g1 * pow(1.0 - exp(-pow(psi_stem / b, c)), beta2);
}
// Leaf::electron_transport's algebra, templated for the same reason.
template <typename T>
T electron_transport_ad(T a, T PPFD, T jmax, T curv) {
  T s = a * PPFD + jmax;
  return (s - sqrt(s * s - 4.0 * curv * a * PPFD * jmax)) / (2.0 * curv);
}
}  // namespace
Leaf::Leaf()
    :
    vcmax_25(96), // umol m^-2 s^-1 
    c(2.680147), //unitless
    b(3.898245), //-MPa
    psi_crit(5.870283), //-MPa 
    root_c(2.680147), //unitless
    root_b(3.898245), //-MPa
    root_psi_crit(5.870283), //-MPa 
    beta2(1.5), //exponent for effect of hydraulic risk (unitless)
    jmax_25(157.44), // maximum electron transport rate umol m^-2 s^-1
    a(0.30), //quantum yield of photosynthetic electron transport (mol mol^-1)
    curv_fact_elec_trans(0.7), //curvature factor for the light response curve (unitless)
    curv_fact_colim(0.99), //curvature factor for the colimited photosythnthesis equatiom
    GSS_tol_abs(1e-3),
    vulnerability_curve_ncontrol(100),
    ci_abs_tol(1e-3),
    ci_niter(1000),
    g1_TF24(7.5), //cost parameter for TF24 profit model umol m^-2 s^-1
    beta_R_H(3.4e2), //proportionality constant between minimum horizontal (intraleyer) root hydraulic resistance and C_r^-1 in [MPa * s * (mol C) / (mol H2O)]
    beta_R_V(9.4e3) //proportionality constant between minimum vertical (interlayer) root hydraulic resistance and dz^2/C_r in [MPa * (mol C) * s / (mol H2O) / m^2]
   {
      setup_transpiration(100); // arg: num control points for integration
      setup_root_vulnerability(100);
      setup_clean_leaf();
}

Leaf::Leaf(double vcmax_25, double c, double b,
           double psi_crit, // derived from b and c,
           double root_c,
           double root_b,
           double root_psi_crit,
           double beta2, double jmax_25,
           double a, double curv_fact_elec_trans, double curv_fact_colim, 
           double GSS_tol_abs,
           double vulnerability_curve_ncontrol,
           double ci_abs_tol,
           double ci_niter,
           double g1_TF24,
           double beta_R_H,
           double beta_R_V)
    : vcmax_25(vcmax_25), // umol m^-2 s^-1 
    c(c), //unitless
    b(b), //-MPa
    psi_crit(psi_crit), //-MPa 
    root_c(root_c), //unitless
    root_b(root_b), //-MPa
    root_psi_crit(root_psi_crit), //-MPa 
    beta2(beta2), //exponent for effect of hydraulic risk (unitless)
    jmax_25(jmax_25), // maximum electron transport rate umol m^-2 s^-1
    a(a), //quantum yield of photosynthetic electron transport (mol mol^-1)
    curv_fact_elec_trans(curv_fact_elec_trans), //curvature factor for the light response curve (unitless)
    curv_fact_colim(curv_fact_colim), //curvature factor for the colimited photosythnthesis equation
    GSS_tol_abs(GSS_tol_abs),
    vulnerability_curve_ncontrol(vulnerability_curve_ncontrol),
    ci_abs_tol(ci_abs_tol),
    ci_niter(ci_niter),
    g1_TF24(g1_TF24), //cost parameter for TF24 profit model umol m^-2 s^-1
    beta_R_H(beta_R_H),
    beta_R_V(beta_R_V)
   {
      setup_transpiration(vulnerability_curve_ncontrol); // arg: num control points for integration
      setup_root_vulnerability(vulnerability_curve_ncontrol);
      setup_clean_leaf();
}

// set various states and physiology parameters obtained from TF24 to NA to clean leaf object
void Leaf::setup_clean_leaf() {
  ci_ = NA_REAL; // Pa
  stom_cond_CO2_= NA_REAL; //mol Co2 m^-2 s^-1 
  assim_colimited_= NA_REAL; // umol C m^-2 s^-1 
  transpiration_= NA_REAL; // kg m^-2 s^-1 
  profit_= NA_REAL; // umol C m^-2 s^-1 
  lambda_= NA_REAL; // umol C m^-2 s^-1 kg^-1 m^2 s^1
  lambda_analytical_= NA_REAL; // umol C m^-2 s^-1 kg^-1 m^2 s^1
  hydraulic_cost_= NA_REAL; // umol C m^-2 s^-1 
  electron_transport_= NA_REAL; //electron transport rate umol m^-2 s^-1
  gamma_= NA_REAL;
  ko_= NA_REAL;
  kc_= NA_REAL;
  km_= NA_REAL;
  R_d_= NA_REAL;
  leaf_specific_conductance_max_= NA_REAL; //kg m^-2 s^-1 MPa^-1 
  sapwood_volume_per_leaf_area_ = NA_REAL; //m^3 SA m^-2 LA
  area_leaf_ = NA_REAL;
  rho_= NA_REAL; //kg m^-3
  vcmax_= NA_REAL; //kg m^-3
  jmax_= NA_REAL; //kg m^-3
  a_bio_= NA_REAL; //kg mol^-1
  root_collar_psi_ = NA_REAL; //-MPa
  leaf_temp_= NA_REAL; // deg C
  PPFD_= NA_REAL; //umol m^-2 s^-1
  atm_vpd_= NA_REAL; //kPa 
  atm_o2_kpa_= NA_REAL; // kPa
  atm_kpa_= NA_REAL; // kPa
  ca_= NA_REAL; //Pa
  opt_psi_stem_= NA_REAL; //-MPa 
  opt_ci_= NA_REAL; //Pa
  E_up_ = NA_REAL;
  collar_pinned_ = false;
  medlyn_model_gs_ = NA_REAL; // mol CO2 m^-2 s^-1 (Medlyn model, develop #450)
  theta_w_ = NA_REAL;
  theta_fc_ = NA_REAL;
  theta_ = NA_REAL;
  psi_soil_.clear();
  soil_depth_.clear();
  z_soil_mid_.clear();  // ADD THIS LINE
  grav_head_z_.clear();
  use_precomputed_z_soil_mid_ = false;
  c_r_V_.clear(); // carbon per layer dedicated to vertical transport (kg m^-2);
  c_r_H_.clear(); // carbon per layer dedicated to horizantal transport (kg m^-2);
  r_R_H_min.clear(); //minimum horizontal portion of root resistance in each soil-layer in [MPa * s * (mol H2O)^-1 m^-2];
  r_R_V.clear(); // vertical root resitance [MPa * s * (mol H2O)^-1 m^-2];
  r_R_V_sum.clear(); // summed vertical root resistance as depth increase;
  soil_consumption_.clear(); // soil consumption mol  m^-2 s^-1;

  soil_number_of_depths_ = NA_INTEGER;
  max_soil_layer = NA_INTEGER; // number of soil layers with root mass greater than 0;

  transpiration_cached_ = false; // invalidate transpiration() memo
  photo_temp_cached_ = false;    // members above set to NA; force recompute
}

// Set the per-individual, per-timestep physiology that stays constant during
// the subsequent root-collar/stem optimisation.
//
// Two groups of quantities are computed here:
//   1. Temperature-dependent photosynthetic parameters (vcmax_, jmax_,
//      electron_transport_, gamma_, ko_, kc_, km_, R_d_) via Arrhenius/peaked
//      Arrhenius functions, plus assim_max_ (assimilation at ci = ca).
//   2. The root hydraulic-resistance network across soil layers. Each layer's
//      root carbon (mass_root_prop[i], kg) is split 1/3 vertical : 2/3
//      horizontal (c_r_V_, c_r_H_). From these:
//        r_R_H_min[i] = beta_R_H / c_r_h      (min horizontal resistance,
//                                              i.e. reciprocal of max conductance)
//        r_R_V[i]     = beta_R_V * dz^2 / c_r_v (vertical resistance; dz^2 because
//                                              vertical conductivity scales with
//                                              root cross-sectional area)
//        r_R_V_sum[i] = cumulative vertical resistance from surface to layer i.
//      max_soil_layer is the deepest layer with non-zero root mass; all
//      resistance vectors are sized to it so the hot E_from_Soil loop only
//      iterates over layers that actually contain roots.
//
// NOTE: the temperature-dependent block (group 1) depends only on leaf_temp_
// (constant across the run in the current driver setup) yet is recomputed on
// every call; see the optimisation notes / caching opportunity.
//
//sets various parameters which are constant for a given node at a given time
void Leaf::set_physiology(double area_leaf, const std::vector<double>& mass_root_prop, double rho, double a_bio, double PPFD, const std::vector<double>& psi_soil, const std::vector<double>& soil_depth, double leaf_specific_conductance_max, double atm_vpd, double ca, double sapwood_volume_per_leaf_area, double leaf_temp, double atm_o2_kpa, double atm_kpa) {
    if (psi_soil.size() != soil_depth.size() || mass_root_prop.size() != soil_depth.size()) {
    util::stop("soil_depth, psi_soil and mass_root_prop must have the same number of elements");
  }
  if (!std::isfinite(area_leaf) || !std::isfinite(rho) || !std::isfinite(a_bio) ||
      !std::isfinite(PPFD) || !std::isfinite(leaf_specific_conductance_max) ||
      !std::isfinite(atm_vpd) || !std::isfinite(ca) ||
      !std::isfinite(sapwood_volume_per_leaf_area) || !std::isfinite(leaf_temp) ||
      !std::isfinite(atm_o2_kpa) || !std::isfinite(atm_kpa)) {
    util::stop("set_physiology received non-finite scalar input");
  }
  for (size_t i = 0; i < psi_soil.size(); ++i) {
    if (!std::isfinite(psi_soil[i])) {
      util::stop("set_physiology received non-finite psi_soil at layer=" + std::to_string(i) +
                 "; psi_soil=" + util::to_string(psi_soil[i]));
    }
    if (!std::isfinite(soil_depth[i])) {
      util::stop("set_physiology received non-finite soil_depth at layer=" + std::to_string(i) +
                 "; soil_depth=" + util::to_string(soil_depth[i]));
    }
    if (!std::isfinite(mass_root_prop[i])) {
      util::stop("set_physiology received non-finite mass_root_prop at layer=" + std::to_string(i) +
                 "; mass_root_prop=" + util::to_string(mass_root_prop[i]));
    }
  }
  area_leaf_ = area_leaf;
  rho_ = rho;
   a_bio_ = a_bio;
   atm_vpd_ = atm_vpd;
   leaf_temp_ = leaf_temp;
   atm_kpa_ = atm_kpa;
   atm_o2_kpa_ = atm_o2_kpa;
   PPFD_ = PPFD;
   psi_soil_ = psi_soil;
   soil_depth_ = soil_depth;
   soil_number_of_depths_ = soil_depth_.size();
   
   if (!(use_precomputed_z_soil_mid_ &&
         z_soil_mid_.size() == static_cast<size_t>(soil_number_of_depths_))) {
     // Fallback for paths that do not provide environment-precomputed midpoints.
     z_soil_mid_.resize(soil_number_of_depths_);
     for (size_t i = 0; i < soil_number_of_depths_; ++i) {
       if (i == 0) {
         z_soil_mid_[i] = (soil_depth_[i] / 2.0);
       } else {
         z_soil_mid_[i] = ((soil_depth_[i - 1] + soil_depth_[i]) / 2.0);
       }
     }
   }

   use_precomputed_z_soil_mid_ = false;

   // Precompute the per-layer gravitational head gravity_head * z_soil_mid_[i];
   // z_soil_mid_ is fixed for the whole solve, so E_from_Soil_to_Root_Collar can
   // read it instead of recomputing the product on every (re)evaluation.
   grav_head_z_.resize(soil_number_of_depths_);
   for (size_t i = 0; i < soil_number_of_depths_; ++i) {
     grav_head_z_[i] = gravity_head * z_soil_mid_[i];
   }

   leaf_specific_conductance_max_ = leaf_specific_conductance_max;
   // conductance changed -> invalidate the transpiration() memo
   transpiration_cached_ = false;
   sapwood_volume_per_leaf_area_ = sapwood_volume_per_leaf_area;
   ca_ = ca;
   // Temperature/O2-dependent block: recomputed only when (leaf_temp_,
   // atm_o2_kpa_) changes from the previous call (see photo_temp_cache_ in the
   // header). Same inputs -> bit-identical outputs, so reusing is exact.
   if (!(photo_temp_cached_ &&
         leaf_temp_ == photo_temp_cache_leaf_temp_ &&
         atm_o2_kpa_ == photo_temp_cache_atm_o2_kpa_)) {
     vcmax_ = peak_arrh_curve(vcmax_ha, vcmax_25, leaf_temp_, vcmax_H_d, vcmax_d_S);
     jmax_ = peak_arrh_curve(jmax_ha, jmax_25, leaf_temp_, jmax_H_d, jmax_d_S);
     gamma_ = arrh_curve(gamma_ha, gamma_25, leaf_temp_);
     ko_ = arrh_curve(ko_ha, ko_25, leaf_temp_);
     kc_ = arrh_curve(kc_ha, kc_25, leaf_temp_);
     R_d_ = vcmax_*0.015;
     km_ = (kc_*umol_per_mol_to_Pa)*(1 + (atm_o2_kpa_*kPa_to_Pa)/(ko_*umol_per_mol_to_Pa));
     photo_temp_cache_leaf_temp_ = leaf_temp_;
     photo_temp_cache_atm_o2_kpa_ = atm_o2_kpa_;
     photo_temp_cached_ = true;
   }
   // depends on the per-call PPFD_ (and cached jmax_), so always recomputed
   electron_transport_ = electron_transport();

   dz_ = soil_depth_.back()/soil_number_of_depths_;


  // find max soil layer as last iteration with mass_root_prop greater than 0
  max_soil_layer = 0;
  for (size_t i = 0; i < soil_number_of_depths_; ++i) {
    if (mass_root_prop[i] != 0) {
      max_soil_layer = i + 1;
    }
  }
  c_r_V_.assign(max_soil_layer, 0.0);
  c_r_H_.assign(max_soil_layer, 0.0);
  r_R_H_min.resize(max_soil_layer);
  r_R_V.resize(max_soil_layer);
  r_R_V_sum.resize(max_soil_layer);

  const double dz_sq = dz_ * dz_;
  double vertical_resistance_sum = 0.0;
  for (size_t i = 0; i < max_soil_layer; ++i) {
    if(mass_root_prop[i] < 0){
            util::stop("Root mass lower than 0");
    }
    const double root_mass = mass_root_prop[i];
    if (root_mass == 0.0) {
      r_R_H_min[i] = 0.0;
      r_R_V[i] = 0.0;
      r_R_V_sum[i] = vertical_resistance_sum;
      continue;
    }

    const double c_r_v = root_mass / 3.0;
    const double c_r_h = root_mass * 2.0 / 3.0;
    c_r_V_[i] = c_r_v;
    c_r_H_[i] = c_r_h;

    // Set horizantal minimum resistance per soil layer (i.e. reciprocal of maximum conductance).
    r_R_H_min[i] = beta_R_H / c_r_h;
    // The vertical conductivity is likely linearly proportional to the root area projected onto the horizontal plane, hence dz^2.
    r_R_V[i] = beta_R_V * dz_sq / c_r_v;
    vertical_resistance_sum += r_R_V[i];
    r_R_V_sum[i] = vertical_resistance_sum;
  }

  // Set up vector of root water uptake from layer
  // .assign, not .resize: the solve writes only up to max_soil_layer and resize's
  // fill reaches only new elements, so deeper layers would hold the previous solve's.
  soil_consumption_.assign(soil_number_of_depths_, 0.0);

  // Soil-moisture state for the Medlyn beta_ stress factor (develop #450). The
  // root-water compute path does not use these; they make the standalone,
  // R-callable Medlyn methods well-defined with the default soil-moisture values.
  theta_w_ = theta_w;
  theta_fc_ = theta_fc;
  theta_ = theta;

  // Find maximum assimilation assuming ci = ca
  assim_max_ = assim_colimited(ca_);
}

// ===========================================================================
// SIGN CONVENTIONS FOR WATER POTENTIAL (psi)  [review #7]
// ---------------------------------------------------------------------------
// This file deliberately uses TWO psi conventions, each natural to its domain.
// They meet at a few clearly-marked "bridge" points that flip with a leading
// minus sign; read those flips with this map in hand:
//
//   * SIGNED (negative) potentials -- the soil -> root-collar transport.
//     psi_soil arrives as positive magnitudes and is flipped once into
//     psi_soil_inverted_ (<= 0). From there P_x_r, the find_root_psi / E_column
//     root variable `x`, find_psi_stem_from_psi_root's psi_root, and
//     transpiration_to_psi_stem's psi_upstream are all SIGNED (<= 0). The
//     physics here uses real signed gradients (psi_soil - P_x_r - gravity*z).
//     The vulnerability splines take a magnitude, so these sites flip back with
//     a leading `-` (e.g. root_vuln_from_psi.eval(-P_src_min)).
//
//   * POSITIVE magnitudes -- the root-collar -> leaf supply. transpiration(),
//     proportion_of_conductivity, hydraulic_cost_TF, psi_stem_to_ci,
//     profit_psi_stem_TF, opt_psi_stem_, psi_crit and the four splines all take
//     a positive magnitude. NB: transpiration() reads eval(psi_upstream)
//     directly while its inverse transpiration_to_psi_stem() reads
//     eval(-psi_upstream): NOT a bug -- they are called with psi_upstream of
//     OPPOSITE sign (positive vs signed), so each is internally consistent.
//
//   * root_collar_psi_ (exported as the opt_root_psi aux) is stored as a SIGNED
//     (negative) potential in ALL branches of find_root_collar_psi (#7 made the
//     Brent / collapsed / root_psi_crit exits agree with the shut-down exits).
//
//   * opt_psi_stem_ (exported as the opt_psi_stem aux) is a POSITIVE magnitude in
//     ALL branches. The assim_max_ < 0 early-exit previously stored the signed
//     root_zero_E here (the lone exception, out of #7 scope); it now stores
//     -root_zero_E so the aux never flips sign by code path.
// ===========================================================================
//
// ---------------------------------------------------------------------------
// SOIL -> ROOT-COLLAR WATER TRANSPORT
// ---------------------------------------------------------------------------
// Scientific model (after Potkay et al. 2021; prototyped in
// vignettes/models/root_water_uptake.Rmd as E_from_Soil_to_Root_Collar):
//
// The root system is represented as a set of parallel soil layers, each
// connected to a single root collar (the point where roots join the stem).
// Within each layer i, water flows from soil to collar driven by the water
// potential gradient (psi_soil[i] - P_x_r), corrected for the gravitational
// head needed to lift water to the layer midpoint (gravity_head * z_soil_mid).
//
// The hydraulic resistance of each layer is the sum of two terms:
//   * r_R_H : horizontal (intra-layer, soil->root) resistance. Set during
//             set_physiology as r_R_H_min[i] / f_r, where r_R_H_min scales
//             with the carbon invested in horizontal roots and f_r is the
//             fractional loss of conductivity from the root vulnerability
//             curve at the operating potential.
//   * r_R_V : vertical (inter-layer, along the root axis to the collar)
//             resistance, accumulated from the surface down to layer i
//             (r_R_V_sum). It scales with dz^2 / carbon-in-vertical-roots.
//
// Because the root vulnerability curve f_r is non-linear in psi, the
// horizontal resistance is evaluated using the *average* fractional
// conductivity over the potential interval spanned between the soil and the
// collar (P_src_min..P_src_max). This mean is obtained as
// (1/(b-a)) * integral_a^b f_r dpsi from a pre-integrated curve
// (root_vuln_integral_from_psi) with two spline evals, the same technique used
// for stem transpiration in setup_transpiration.
//
// Output: E_up_ = total water drawn from all layers to the collar
//                 (converted to kg H2O m^-2 leaf s^-1), and soil_consumption_[i]
//                 = per-layer uptake (mol H2O m^-2 leaf s^-1). Negative E_i in a
//                 layer means that layer is *gaining* water (hydraulic redistribution).
//
// Implementation decisions:
//   * f_r and its running integral are read from pre-computed splines
//     (root_vuln_from_psi, root_vuln_integral_from_psi) instead of repeatedly
//     evaluating exp(-(psi/b)^c); see setup_root_vulnerability.
//   * Two special cases are handled exactly to avoid division/round-off issues:
//     (a) collar potential equals layer potential, and (b) the gradient
//     exactly balances gravity (E_i = 0).
//   * Extensive isfinite() guards are present because this function is called
//     from within nested root-finders where bad brackets can produce NaNs;
//     they fail fast with diagnostic context rather than propagating NaN.
//
// This function calculates the total transpiration from the soil based on the
// root collar pressure and the respective soil layer pressures
void Leaf::E_from_Soil_to_Root_Collar(double P_x_r, const std::vector<double>& psi_soil){

    if (!std::isfinite(P_x_r) || !std::isfinite(area_leaf_)) {
      util::stop("E_from_Soil_to_Root_Collar invalid input; P_x_r=" + util::to_string(P_x_r) +
                 "; area_leaf_=" + util::to_string(area_leaf_));
    }

    E_up_ = 0;

    // area_leaf_ is constant across the whole solve; fold its reciprocal into a
    // per-layer multiply instead of a per-layer division (1 fdiv/call vs 15).
    const double inv_area_leaf = 1.0 / area_leaf_;

    // Cumulative-integral spline caching (bit-identical fast path). The only two
    // arguments ever passed to root_vuln_integral_from_psi in the loop below are
    // -P_src_min and -hi_neg, each of which resolves to exactly one of
    // {-psi_soil[i], -P_x_r, 0}. -P_x_r is constant across all layers (compute
    // once), and -psi_soil[i] is constant across the whole solve (precomputed in
    // find_root_collar_psi). We only take this path when psi_soil is the cached
    // psi_soil_inverted_ vector; any other caller falls back to direct evals.
    const bool use_integral_cache =
        (&psi_soil == &psi_soil_inverted_) &&
        root_vuln_integral_soil_.size() == static_cast<size_t>(max_soil_layer);
    const double neg_P_x_r = -P_x_r;
    const double G_at_P_x_r =
        use_integral_cache ? root_vuln_integral_from_psi.eval(neg_P_x_r) : 0.0;

    // GUARD POLICY (the per-layer isfinite/stop guards here were added while
    // debugging the #485 drought-NaN, now fixed at source by the soil residual-
    // moisture floor). Most were defensive and redundant, so they have been
    // removed from this hot loop; the remaining two are load-bearing:
    //   * the equal-potentials f_ri <= 0 check below: root_vuln_from_psi
    //     LINEARLY extrapolates NEGATIVE beyond its domain, so a deep-drought
    //     layer can produce negative conductivity -> negative-but-FINITE r_R ->
    //     wrong-sign E_i that the post-loop isfinite(E_up_) net would NOT catch.
    //   * the post-loop isfinite(E_up_) check: any non-finite produced anywhere
    //     in the loop propagates into the sum and is caught there once per call.
    // Everything else is provably safe to drop on the valid path: psi_soil is
    // validated in set_physiology; P_src_min<=P_src_max by construction; the
    // general-branch integral comes from a monotone-increasing spline so it is
    // strictly > 0 (span>0), giving r_R>0 and finite E_i; and any stray NaN/Inf
    // still reaches the post-loop net.
    for(size_t i = 0; i < max_soil_layer; i++){

    // Find the most negative soil potential out of the given soil layer and the root collar
    double P_src_min = std::min(psi_soil[i], P_x_r);

    // Find the least negative soil potential out of the given soil layer and the root collar
    double P_src_max = std::max(psi_soil[i], P_x_r);

     // If root collar soil water potential equals the soil water potential in a given layer
    if(std::abs(P_x_r - psi_soil[i]) < 1e-8){

      // Fraction of conductance in roots in a given layer at most negative soil water potential (but actually is equal to root collar)
      // root_vuln_from_psi is a pre-built spline of exp(-(|psi|/b_root)^c_root)
      double f_ri = root_vuln_from_psi.eval(-P_src_min);
      if (!std::isfinite(f_ri) || f_ri <= 0.0) {
        util::stop("E_from_Soil_to_Root_Collar invalid f_ri; layer=" + std::to_string(i) +
                   "; f_ri=" + util::to_string(f_ri) +
                   "; P_src_min=" + util::to_string(P_src_min) +
                   "; P_x_r=" + util::to_string(P_x_r));
      }

      // Fraction of conductance in roots in a given layer at most negative soil water potential
      double r_R_H = r_R_H_min[i] / f_ri; // [MPa * s * (mol H2O)^-1]

      // Total root resistance (horizantal plus vertical)
      double r_R = r_R_H + r_R_V_sum[i];

      // Transpiration is equivalent to gravitational water loss (i.e. layer gains water)
      double E_i = -grav_head_z_[i] * inv_area_leaf / r_R ;

      soil_consumption_[i] = E_i;
      E_up_ += E_i;

    }
    else if(std::abs((psi_soil[i] - P_x_r) - grav_head_z_[i]) < 1e-8){
      // If pressure difference perfectly balances gravity transpiration is equal to zero
      double E_i = 0.0; // [mol H2O / m^2 / s]
      
      soil_consumption_[i] = E_i;

      E_up_ += E_i;

    } else{

      // Mean fractional root conductivity over the potential interval
      // [P_src_min, P_src_max], i.e. (1/(b-a)) * integral_a^b f_r dpsi.
      // Computed from the pre-integrated curve G(m) = integral_0^m f_r(s) ds
      // (root_vuln_integral_from_psi, indexed by magnitude m = -psi) with 2
      // evals instead of the old (n+1)-point sample mean. The interval is split
      // at psi = 0: for psi > 0 (above-atmospheric) vulnerability is 1.
      double hi_neg = std::min(P_src_max, 0.0); // boundary of the psi<=0 part
      double lo_pos = std::max(P_src_min, 0.0); // boundary of the psi>0 part

      // Memoised cumulative-integral lookup. Returns the exact same double the
      // spline would (same input -> same output); the comparisons select the
      // precomputed value because -P_src_min / -hi_neg are bit-for-bit equal to
      // one of the cached arguments in the common (psi<=0) case.
      const double neg_psi_soil_i = -psi_soil[i];
      auto G_integral = [&](double arg) -> double {
        if (use_integral_cache) {
          if (arg == neg_P_x_r) return G_at_P_x_r;
          if (arg == neg_psi_soil_i) return root_vuln_integral_soil_[i];
        }
        return root_vuln_integral_from_psi.eval(arg);
      };

      double integral = 0.0;
      if (hi_neg > P_src_min) {
        // psi<=0 part: magnitude m runs from -hi_neg up to -P_src_min
        integral += G_integral(-P_src_min) - G_integral(-hi_neg);
      }
      if (P_src_max > lo_pos) {
        // psi>0 part: f_r == 1 over its length
        integral += (P_src_max - lo_pos);
      }

    // span = P_src_max - P_src_min > 0 here (the equal-potentials case is
    // handled in the branch above). integral comes from the monotone-increasing
    // cumulative-vulnerability spline so it is strictly > 0 over a span>0
    // interval; forming r_R_H as r_R_H_min * span / integral is one division
    // (vs the old f_r_average = integral/span then r_R_H_min/f_r_average two),
    // and needs no per-layer finiteness guard (any stray NaN/Inf propagates to
    // the post-loop isfinite(E_up_) net).
    const double span = P_src_max - P_src_min;

    // Find the horizantal resistance in a given layer by dividing the minimum resistance (i.e. maximum conductivity) by the fractional loss of conductivity
    double r_R_H = r_R_H_min[i] * span / integral; // [MPa * s * (mol H2O)^-1]

    // Find the total resistance in a given layer by adding the vertical resistance in that layer
    double r_R = r_R_H + r_R_V_sum[i]; // [MPa * s * (mol H2O)^-1]

    // Transpiration is equal to the potentail gradient between the root collar and the soil, accounting for gravitational potential
    double E_i = (psi_soil[i] - P_x_r - grav_head_z_[i]) * inv_area_leaf / r_R; // [mol H2O / m^2 / s]

    soil_consumption_[i] = E_i;
    E_up_ += E_i;

    }
  }
  // Convert the summed uptake to kg H2O m^-2 s^-1, consistent with the rest of
  // the leaf model and environment. NOTE (review #10): only the aggregate E_up_
  // is converted to kg here; the per-layer soil_consumption_[i] above is left in
  // mol H2O m^-2 s^-1 and converted downstream in TF24_Strategy::compute_rates.
  // The two siblings therefore carry different units by design.
  E_up_ = E_up_ * kg_per_mol_h2o;
  if (!std::isfinite(E_up_)) {
    util::stop("E_from_Soil_to_Root_Collar non-finite E_up_; P_x_r=" + util::to_string(P_x_r) +
               "; max_soil_layer=" + std::to_string(max_soil_layer) +
               "; area_leaf_=" + util::to_string(area_leaf_));
  }
}



// This function is used to find root collar pressure which equilibrates the soil-root-stem water continuuum
double Leaf::E_column(double x, const std::vector<double>& psi_soil, double psi_leaf) {

  E_from_Soil_to_Root_Collar(x, psi_soil);
  root_collar_psi_ = -x;
  double E_root_to_leaf = transpiration(psi_leaf, root_collar_psi_);
  return E_up_ - E_root_to_leaf;
}

// This function is used to find root collar pressure where water form soil is equal to zero
double Leaf::E_column_zero(double x, const std::vector<double>& psi_soil) {

  E_from_Soil_to_Root_Collar(x, psi_soil);

  return E_up_;
}

// find root psi based on required condition, i.e. equilibrated continuum, zero water from soil
//
// #486: both targets (E_column / E_column_zero, the soil->collar continuity
// residual over the collar potential x in [-psi_crit, wettest_soil_layer]) are
// smooth and strictly monotone in their *normal operating regime* -- a clean
// single sign-change with derivatives continuous across every x == psi_soil[i]
// layer crossing (the per-layer branch switches in E_from_Soil_to_Root_Collar
// are bit-level kinks, relative slope jump ~1e-6, not real corners). So a
// superlinear bracketing solver is safe and faster here: TOMS748 reaches the
// same root in ~6-8 E_from_Soil evals vs bisection's ~15-16 at the same 1e-4
// tol (see the test-leaf.r "find_root_psi soil->collar continuity solve"
// contract block). This directly attacks the dominant E_from_Soil per-layer
// arithmetic hot-spot, evaluated ~270x per collar solve through this finder.
//
// SCOPE/CAVEAT: the brackets here are guaranteed valid (opposite-sign, finite
// endpoints) by find_root_collar_psi's preceding early-exits -- the crit=1
// lower endpoint is exactly the E_column(-psi_crit) < 0 shutdown test, and
// crit=0 is only reached for soil wetter than psi_crit. The genuinely
// non-smooth failure mode in the earlier blanket-swap rejection (the root
// vulnerability spline extrapolating negative beyond its ~root_psi_crit domain)
// lives in E_from_Soil_to_Root_Collar itself and bites both solvers identically;
// it does not arise on the brackets this finder is actually handed. Like
// psi_stem_to_ci (Phase 6) this is a same-tolerance method swap, NOT a tolerance
// loosening: same root, fewer evals.
double Leaf::find_root_psi(double wettest_soil_layer, const std::vector<double>& psi_soil, int find_root_crit) {
  // tol and iterations copied from control defaults (for now) - changed recently to 1e-6
  if (find_root_crit == 1) {
    auto target = [&](double x) -> double {
      return E_column(x, psi_soil, psi_crit);
    };
    try {
      return util::uniroot_smooth(target, -psi_crit, wettest_soil_layer, 1e-4, ci_niter);
    } catch (const std::exception& e) {
      util::stop("find_root_psi(find_root_crit=1) failed: " + std::string(e.what()) +
                 "; min=" + util::to_string(-psi_crit) +
                 "; max=" + util::to_string(wettest_soil_layer));
    }
  }

  auto target = [&](double x) -> double {
    return E_column_zero(x, psi_soil);
  };
  try {
    return util::uniroot_smooth(target, -psi_crit, wettest_soil_layer, 1e-4, ci_niter);
  } catch (const std::exception& e) {
    util::stop("find_root_psi(find_root_crit=0) failed: " + std::string(e.what()) +
               "; min=" + util::to_string(-psi_crit) +
               "; max=" + util::to_string(wettest_soil_layer));
  }

}

// When root pressure is known, find E from soil, then use E from soil to find psi stem
double Leaf::find_psi_stem_from_psi_root(double psi_root, const std::vector<double>& psi_soil){
  E_from_Soil_to_Root_Collar(psi_root, psi_soil);

  double psi_stem = transpiration_to_psi_stem(E_up_, psi_root);
  return psi_stem;
}

// ---------------------------------------------------------------------------
// MASTER SOLVER: optimal root-collar (and stem) water potential
// ---------------------------------------------------------------------------
// This is the entry point called once per individual per environment update
// (from TF24_Strategy::net_mass_production_dt). It solves the whole
// soil -> root -> stem -> leaf hydraulic continuum and stores the optimal
// operating point in opt_psi_stem_, root_collar_psi_ and profit_.
//
// The solve has two nested levels:
//
//   1. CONTINUITY (find_root_psi / E_column): for any candidate root-collar
//      potential, water supplied from the soil (E_from_Soil_to_Root_Collar)
//      must equal water transpired through the stem (transpiration()). This is
//      a 1-D root-find on the collar potential.
//
//   2. OPTIMISATION (Golden-Section Search): among feasible collar potentials,
//      choose the one that maximises carbon profit = assimilation - hydraulic
//      cost (profit_psi_stem_TF). The collar potential is bracketed between
//      `root_zero_E` (collar where soil uptake is zero, the wettest feasible
//      point) and `root_crit` (collar at which the stem reaches psi_crit, the
//      driest feasible point), clamped to root_psi_crit.
//
// Several early-exit short-circuits avoid the (expensive) GSS loop when no
// meaningful optimisation is possible. In each case the plant is effectively
// shut down (operating at psi_crit, paying only respiration + hydraulic cost):
//   * wettest soil layer is already drier than psi_crit -> no transpiration;
//   * even at psi_crit the soil cannot supply the demanded flux (E_column<0);
//   * the continuity root would require the collar drier than psi_crit;
//   * maximum possible assimilation (at ci = ca) is negative.
//
// Implementation note: psi_soil arrives as positive magnitudes and is used
// here as negative potentials, hence psi_soil_inverted_. The GSS reuses one
// profit evaluation per iteration (golden ratio) to halve function calls, and
// a collapsed-interval branch handles the degenerate single-feasible-point case.
// Shut-down operating point shared by find_root_collar_psi's early-exits: the
// stem is held at psi_crit (transpiration not possible), so the plant pays only
// respiration (R_d_) plus the hydraulic cost at psi_crit. Only the recorded
// root-collar potential differs between the calling cases.
void Leaf::set_shutdown_state(double root_collar) {
  root_collar_psi_ = root_collar;
  opt_psi_stem_ = psi_crit;
  profit_ = -R_d_ - hydraulic_cost_TF(psi_crit);
  // Stomata are shut, so gross assimilation is zero and the reported net rate
  // is -R_d_. Set explicitly: this branch does not go through
  // profit_psi_stem_TF, so assim_colimited_ would otherwise be left at
  // whatever the last probe wrote, and it is now reported as an aux variable.
  // Keeps profit_ == assim_colimited_ - hydraulic_cost_TF() in every branch.
  assim_colimited_ = -R_d_;
  // Shut down means no water movement, so zero the whole transport chain.
  // This matters beyond diagnostics: the first caller below returns before any
  // E_from_Soil_to_Root_Collar call in this solve, and `Leaf` is a value member
  // reused across every compute_rates call for an individual. Left alone,
  // soil_consumption_ therefore keeps the *previous* step's values, and
  // TF24_Strategy::evapotranspiration_dt feeds those straight into the patch
  // water balance -- a plant that has closed its stomata carries on drawing its
  // last wet-step uptake out of the soil. Note the water budget still *closes*
  // in that state (what is recorded as depleted is what is removed), so the
  // conservation tests cannot catch it; only the physics is wrong.
  transpiration_ = 0.0;
  stom_cond_CO2_ = 0.0;
  E_up_ = 0.0;
  std::fill(soil_consumption_.begin(), soil_consumption_.end(), 0.0);
}

// psi_soil_ arrives as positive magnitudes; flip once to the signed convention the
// soil->collar transport uses, and precompute each layer's cumulative-integral
// lookup at that same argument. Returns the wettest layer's potential.
double Leaf::refresh_soil_potentials() {
  psi_soil_inverted_.resize(max_soil_layer);
  root_vuln_integral_soil_.resize(max_soil_layer);
  double wettest_soil_layer = -std::numeric_limits<double>::infinity();
  for (size_t i = 0; i < max_soil_layer; ++i) {
    const double psi_inverted = -psi_soil_[i];
    psi_soil_inverted_[i] = psi_inverted;
    root_vuln_integral_soil_[i] =
        root_vuln_integral_from_psi.eval(-psi_inverted);
    wettest_soil_layer = std::max(wettest_soil_layer, psi_inverted);
  }
  return wettest_soil_layer;
}

// Shared setup + feasibility handling for the root-collar solve. Extracted
// verbatim from find_root_collar_psi (no reordering of floating-point ops, so
// the TF24 optimisation path stays bit-identical) and reused by
// evaluate_root_collar_psi. Returns false when the final operating point is
// already determined here (shutdown / assim<0 / collapsed interval) and the
// caller should stop; returns true with [bound_a, bound_b] set to the feasible
// collar-potential interval (positive magnitudes) otherwise.
bool Leaf::prepare_collar_solve(double& bound_a, double& bound_b){

  // Only polish_root_collar_psi sets this, so an exit below would otherwise
  // leave the previous step's pinned flag in place.
  collar_pinned_ = false;

  const double wettest_soil_layer = refresh_soil_potentials();

  // Avoid loop if the wettest psi layer is drier than psi_crit in stem, transpiration not possible and so all variables set to
  // shut down

  if (-wettest_soil_layer >= psi_crit){
    set_shutdown_state(-psi_crit);
    return false;
  }

if(E_column(-psi_crit, psi_soil_inverted_, psi_crit) < 0){
      // root_collar_psi_ is reported as a signed (negative) potential, so store
      // -root_psi_crit rather than the positive magnitude root_psi_crit.
      set_shutdown_state(-root_psi_crit);
      return false;
}

  // Avoid loop if the wettest psi layer is drier than psi_crit in stem, transpiration not possible and so all variables set to
  // shut down
double root_crit = find_root_psi(wettest_soil_layer, psi_soil_inverted_, 1);

// If root crit would have to be larger than psi crit, also avoid loop as above

    if (-root_crit >= psi_crit){
    set_shutdown_state(root_crit);
    return false;
  }

// Find root collar where transpiration from soil is 0
double root_zero_E = find_root_psi(wettest_soil_layer, psi_soil_inverted_, 0);

// If assimilation would be less than 0 even at Ca, also end loop
if(assim_max_ < 0){
    // At zero transpiration the stem equilibrates with the collar (no flux, no
    // gradient), so the operating point is root_zero_E for both. root_collar_psi_
    // is the signed (negative) potential (#7); opt_psi_stem_ is the matching
    // positive magnitude (-root_zero_E), keeping it sign-consistent with every
    // other branch of this solver.
    opt_psi_stem_ = -root_zero_E;
    root_collar_psi_ = root_zero_E;
    E_from_Soil_to_Root_Collar(root_collar_psi_, psi_soil_inverted_);

    profit_ = - R_d_ - hydraulic_cost_TF(-root_collar_psi_);
    // As in set_shutdown_state: transpiration is zero here, so gross
    // assimilation is zero and the reported net rate is -R_d_.
    assim_colimited_ = -R_d_;
    // E_up_ and soil_consumption_ are already correct: the
    // E_from_Soil_to_Root_Collar call above evaluates them at root_zero_E, the
    // collar potential at which uptake is zero. The leaf-side pair is not set
    // anywhere on this path, though, so zero it here rather than leave the
    // previous step's values (see set_shutdown_state for why that matters).
    transpiration_ = 0.0;
    stom_cond_CO2_ = 0.0;

        if(std::isnan(profit_)){
          util::stop("Error: profit nan");
    }

    return false;
}
// opt_psi_stem_ = psi_soil_;


  // optimise for stem water potential
    bound_a = -root_zero_E;
    bound_b = std::max(-root_crit,-root_psi_crit);

    // If no interval exists (single feasible root-collar value), use that
    // point directly as the alternative solution instead of running GSS.
    if (std::abs(bound_b - bound_a) <= GSS_tol_abs) {
      const double opt_root_psi = 0.5 * (bound_a + bound_b);
      const double psi_stem_single = find_psi_stem_from_psi_root(-opt_root_psi, psi_soil_inverted_);

      if (!std::isfinite(psi_stem_single)) {
        util::stop("Error: non-finite psi_stem_single in collapsed-root interval; "
                   "opt_root_psi=" + util::to_string(opt_root_psi) +
                   "; bound_a=" + util::to_string(bound_a) +
                   "; bound_b=" + util::to_string(bound_b) +
                   "; root_crit=" + util::to_string(root_crit) +
                   "; root_zero_E=" + util::to_string(root_zero_E) +
                   "; E_up_=" + util::to_string(E_up_));
      }

      opt_psi_stem_ = psi_stem_single;
      // profit_psi_stem_TF takes psi_upstream as a positive magnitude, so feed
      // it opt_root_psi; root_collar_psi_ is stored as the signed (negative)
      // potential for a sign-consistent aux output.
      profit_ = profit_psi_stem_TF(opt_psi_stem_, opt_root_psi);
      root_collar_psi_ = -opt_root_psi;

      if (!std::isfinite(profit_)) {
        util::stop("Error: non-finite profit in collapsed-root interval; "
                   "opt_psi_stem_=" + util::to_string(opt_psi_stem_) +
                   "; root_collar_psi_=" + util::to_string(root_collar_psi_) +
                   "; bound_a=" + util::to_string(bound_a) +
                   "; bound_b=" + util::to_string(bound_b) +
                   "; root_crit=" + util::to_string(root_crit) +
                   "; root_zero_E=" + util::to_string(root_zero_E) +
                   "; E_up_=" + util::to_string(E_up_) +
                   "; assim_colimited_=" + util::to_string(assim_colimited_) +
                   "; hydraulic_cost_=" + util::to_string(hydraulic_cost_));
      }
      return false;
    }

    return true;
}

void Leaf::find_root_collar_psi(){
    double bound_a, bound_b;
    if (!prepare_collar_solve(bound_a, bound_b)) {
      if (internals::collar_census_active()) {
        internals::the_collar_census().add_early_exit();
      }
      return;
    }
    if (internals::collar_census_active()) {
      internals::the_collar_census().add_bracket(std::abs(bound_b - bound_a));
    }
    // root_crit / root_zero_E were consumed inside prepare_collar_solve; recover
    // them for the diagnostic message only if the profit check below fails.

    // Maximise carbon profit over the feasible collar-potential interval via
    // golden-section search (util::golden_section_max). Unlike Brent, its argmax
    // is a smooth (fixed-iteration) function of the inputs, so the operating
    // point varies smoothly with plant height -- the demographic growth-rate
    // gradient relies on this. The objective maps a candidate collar potential
    // `bound` to its profit (find the stem psi it implies, then evaluate profit).
    const double search_root_psi = util::golden_section_max(
        [&](double bound) {
          const double psi_stem =
              find_psi_stem_from_psi_root(-bound, psi_soil_inverted_);
          return profit_psi_stem_TF(psi_stem, bound);
        },
        bound_a, bound_b, GSS_tol_abs);

    // The search stops at its bracket tolerance, where d(profit)/d(collar) is
    // still of the bracket's size; Newton takes it to a stationary point.
    const double opt_root_psi =
        polish_root_collar_psi(search_root_psi, bound_a, bound_b);

    // Every operating-point output is set from opt_root_psi below, so the probe
    // points polish_root_collar_psi left behind are overwritten here.
    opt_psi_stem_ = find_psi_stem_from_psi_root(-opt_root_psi, psi_soil_inverted_);

    // store as the signed (negative) potential for a sign-consistent aux output;
    // profit_psi_stem_TF takes psi_upstream as a positive magnitude.
    root_collar_psi_ = -opt_root_psi;
    profit_ = profit_psi_stem_TF(opt_psi_stem_, opt_root_psi);

    if(!std::isfinite(profit_)){
        util::stop("Error: non-finite profit; opt_psi_stem_=" + util::to_string(opt_psi_stem_) +
             "; root_collar_psi_=" + util::to_string(root_collar_psi_) +
             "; bound_a=" + util::to_string(bound_a) +
             "; bound_b=" + util::to_string(bound_b) +
             "; E_up_=" + util::to_string(E_up_) +
             "; assim_colimited_=" + util::to_string(assim_colimited_) +
             "; hydraulic_cost_=" + util::to_string(hydraulic_cost_));
    }
}

// Newton on R = dprofit_droot_collar_psi (see header). The envelope relation the
// reverse pass uses holds only where R is zero, and the bracket search returns a
// point whose R is of the bracket's size.
//
// dR/d(collar) is a central difference of the analytic R. It is all double and it
// enters only as the divisor of a step length, so its error shortens or lengthens
// a step and does not move the point the steps converge to.
//
// Leaves the operating-point outputs at the last probe point of
// dprofit_droot_collar_psi; the caller must set them from the returned potential.
double Leaf::polish_root_collar_psi(double opt_root_psi, double bound_a,
                                   double bound_b) {
  const double h = 1e-6;          // step of the dR/d(collar) difference
  const double R_tol = 1e-11;
  // Exhausting max_iter returns the last iterate, whose R is of the bracket's
  // size rather than R_tol, so the envelope relation does not hold there.
  const int max_iter = 20;

  double psi = opt_root_psi;
  double dR_dcollar = 0.0;
  bool have_dR = false;   // dR_dcollar holds a value, taken at psi or before it
  int n_eval = 0;
  int cls = internals::COLLAR_EXHAUSTED;
  collar_pinned_ = false;
  for (int iter = 0; iter < max_iter; ++iter) {
    const double R = dprofit_droot_collar_psi(psi);
    ++n_eval;
    if (!std::isfinite(R) || std::abs(R) <= R_tol) {
      cls = std::isfinite(R) ? internals::COLLAR_INTERIOR
                             : internals::COLLAR_R_NONFINITE;
      break;
    }
    // Every remaining step needs room inside the bracket. Running out of it is
    // the pinned case: the maximum is the bound, profit is not stationary there,
    // and the point the search returned is the answer.
    if (psi - h <= bound_a || psi + h >= bound_b) {
      cls = (psi - h <= bound_a) ? internals::COLLAR_BOUND_A
                                 : internals::COLLAR_BOUND_B;
      break;
    }
    // A step after the first reuses the derivative already held: it divides a
    // step length, so a stale one changes the step and not the point reached.
    bool stepped = false;
    bool dR_at_psi = false;
    while (!stepped) {
      if (!have_dR) {
        dR_dcollar = (dprofit_droot_collar_psi(psi + h) -
                      dprofit_droot_collar_psi(psi - h)) / (2.0 * h);
        n_eval += 2;
        have_dR = true;
        dR_at_psi = true;
      }
      // Profit is concave at an interior maximum, so a non-negative curvature
      // says this is not one and a Newton step would leave it.
      if (std::isfinite(dR_dcollar) && dR_dcollar < 0.0) {
        const double psi_next = psi - R / dR_dcollar;
        if (std::isfinite(psi_next) && psi_next > bound_a &&
            psi_next < bound_b) {
          psi = psi_next;
          stepped = true;
        }
      }
      // A rejected step is retried once against a derivative taken at psi; a
      // fresh derivative that is still rejected is the bound case.
      if (!stepped) {
        if (dR_at_psi) {
          cls = (std::isfinite(dR_dcollar) && dR_dcollar < 0.0)
                    ? internals::COLLAR_BOUND_STEP
                    : internals::COLLAR_BOUND_CURVATURE;
          break;
        }
        have_dR = false;
      }
    }
    if (!stepped) {
      break;
    }
  }

  // Recorded whether or not a gradient is taken. The four exits that ran out of
  // bracket leave the collar potential at a bound of the feasible interval.
  collar_pinned_ = (cls == internals::COLLAR_BOUND_A ||
                    cls == internals::COLLAR_BOUND_B ||
                    cls == internals::COLLAR_BOUND_STEP ||
                    cls == internals::COLLAR_BOUND_CURVATURE);

  if (std::getenv("PLANT_POLISH_TRACE") != nullptr) {
    std::fprintf(stderr, "polish_root_collar_psi: n_eval=%d displacement=%.3e\n",
                 n_eval, psi - opt_root_psi);
  }

  if (internals::collar_census_active()) {
    // R at the returned point costs one more evaluation, and it leaves the
    // operating-point outputs at psi, which is where the caller resets them from.
    double psi_wet = psi_soil_.empty() ? 0.0 : psi_soil_[0];
    for (size_t i = 1; i < psi_soil_.size(); ++i) {
      if (psi_soil_[i] < psi_wet) psi_wet = psi_soil_[i];
    }
    internals::the_collar_census().add(
        cls, std::abs(dprofit_droot_collar_psi(psi)), psi_wet, PPFD_,
        area_leaf_);
  }
  return psi;
}

// Evaluate the operating point at a given collar potential rather than
// optimising it (see header). Clamps the target into the feasible interval so a
// tracked state that has drifted outside it still yields a finite operating
// point; the gradient computed by the caller then pulls it back inside.
double Leaf::evaluate_root_collar_psi(double target_opt_root_psi){
    double bound_a, bound_b;
    if (!prepare_collar_solve(bound_a, bound_b)) {
      // Operating point fully determined by feasibility handling (shutdown /
      // assim<0 / collapsed interval); profit_ is already set.
      return profit_;
    }

    return profit_at_collar_psi(target_opt_root_psi, bound_a, bound_b);
}

// Post-prepare body of evaluate_root_collar_psi (see header). Kept as a separate
// entry point so callers that evaluate several collar potentials within one step
// (the centred finite difference, #530) can run prepare_collar_solve once and
// reuse the soil-side caches across every profit eval. The clamp into
// [bound_a, bound_b] is identical to evaluate_root_collar_psi's, so near a
// boundary a perturbed potential collapses onto the boundary -- which is exactly
// how the FD path degrades gracefully to a one-sided difference.
double Leaf::profit_at_collar_psi(double target_opt_root_psi,
                                  double bound_a, double bound_b){
    const double opt_root_psi =
        std::min(std::max(target_opt_root_psi, bound_a), bound_b);

    opt_psi_stem_ = find_psi_stem_from_psi_root(-opt_root_psi, psi_soil_inverted_);
    root_collar_psi_ = -opt_root_psi;
    profit_ = profit_psi_stem_TF(opt_psi_stem_, opt_root_psi);

    if(!std::isfinite(profit_)){
        util::stop("Error: non-finite profit in evaluate_root_collar_psi; "
             "target=" + util::to_string(target_opt_root_psi) +
             "; opt_root_psi=" + util::to_string(opt_root_psi) +
             "; bound_a=" + util::to_string(bound_a) +
             "; bound_b=" + util::to_string(bound_b));
    }
    return profit_;
}

// Exact d(profit)/d(opt_root_psi). profit(psi) = assim_colimited(ci) -
// hydraulic_cost_TF(psi_stem), with psi_stem = find_psi_stem_from_psi_root(-psi)
// (smooth spline transport) and ci = psi_stem_to_ci(psi_stem, psi) (root-find).
// Chain rule:
//   dprofit/dpsi = A'(ci) dci/dpsi - C'(psi_stem) dpsi_stem/dpsi
// where dci/dpsi = (dci/dpsi_stem) dpsi_stem/dpsi + (dci/dpsi)|_explicit, and the
// dci/d* terms come from the implicit-function theorem on the residual
//   g(ci; psi_stem, psi) = A(ci) umol_to_mol - gc(psi_stem,psi) (ca-ci)/(atm kPa)
// with gc = const * transpiration(psi_stem,psi). A'/C' are obtained by forward
// AD; the gc partials use the analytic spline derivative (transpiration_from_psi
// .deriv); dpsi_stem/dpsi by a tight central difference on the smooth transport.
double Leaf::dprofit_droot_collar_psi(double opt_root_psi) {
  const double psi = opt_root_psi;
  const double gstar_Pa = gamma_ * umol_per_mol_to_Pa;

  // Every transport evaluation below reads psi_soil_inverted_, so seat it on the
  // current psi_soil_ here rather than depending on the caller's last solve.
  refresh_soil_potentials();

  // Operating point in double.
  const double psi_stem = find_psi_stem_from_psi_root(-psi, psi_soil_inverted_);
  // Shut down before the ci solve, not after. psi and psi_stem are positive
  // magnitudes here, so psi >= psi_stem is the no-flow / reversed-gradient case
  // -- the same condition set_leaf_states_rates_from_psi_stem() treats as zero
  // transpiration. It has to be caught *here* because psi_stem_to_ci() does not
  // return non-finite in that state, it throws: gc = const * transpiration goes
  // negative, which flips the sign of the supply term so the residual no longer
  // crosses zero over (gamma*, ca] and TOMS748 reports "a and b do not bracket
  // the root". The isfinite check below was written to cover shut-down but
  // cannot see a thrown exception, so a dry patch killed the whole run:
  // reproduced on TF24f at 5 layers, theta = 0.005-0.03 with 1 m/yr rainfall,
  // at psi_stem = 1.23 against psi_upstream = 5.92 MPa.
  if (!std::isfinite(psi_stem) || psi >= psi_stem) {
    return 0.0;  // shut-down / infeasible: no informative gradient
  }
  const double ci = psi_stem_to_ci(psi_stem, psi);
  if (!std::isfinite(ci)) {
    return 0.0;
  }

  // A'(ci) and C'(psi_stem) via forward-mode AD of the analytic algebra.
  const double A_prime =
      odelia::ode::forward_derivative(ci, [&](auto x) -> decltype(x) {
        return assim_colimited_ad<decltype(x)>(x, vcmax_, electron_transport_,
                                               gstar_Pa, km_, R_d_,
                                               curv_fact_colim);
      });
  const double C_prime =
      odelia::ode::forward_derivative(psi_stem, [&](auto x) -> decltype(x) {
        return hydraulic_cost_ad<decltype(x)>(x, b, c, g1_TF24, beta2);
      });

  // Stomatal-conductance supply coefficient gc and its partials. gc =
  // gc_const * transpiration(psi_stem, psi); transpiration is conductance_max *
  // (transp_from_psi(psi_stem) - transp_from_psi(psi)), so the partials use the
  // analytic spline derivative.
  const double gc_const =
      atm_kpa_ * kg_to_mol_h2o / atm_vpd_ / H2O_CO2_stom_diff_ratio;
  const double gc = gc_const * transpiration(psi_stem, psi);
  const double dgc_dpsistem =
      gc_const * leaf_specific_conductance_max_ * transpiration_from_psi.deriv(psi_stem);
  const double dgc_dpsi =
      gc_const * leaf_specific_conductance_max_ * (-transpiration_from_psi.deriv(psi));

  // IFT on g(ci; psi_stem, psi): dci/dp = -(dg/dp)/(dg/dci).
  const double inv_atm = 1.0 / (atm_kpa_ * kPa_to_Pa);
  const double g_ci = A_prime * umol_to_mol + gc * inv_atm;      // dg/dci
  const double dci_dpsistem = -(-dgc_dpsistem * (ca_ - ci) * inv_atm) / g_ci;
  const double dci_dpsi_expl = -(-dgc_dpsi * (ca_ - ci) * inv_atm) / g_ci;

  // dpsi_stem/dpsi: psi_stem = P(E_psi_stem) with
  //   E_psi_stem = E_up_(r)/k_max + S(psi),   r = -psi,
  // S = transpiration_from_psi, P = psi_from_transpiration (both C2 splines), and
  // E_up_(r) the soil->collar uptake. Chain rule, with dr/dpsi = -1:
  //   dE_psi_stem/dpsi = -E_up_'(r)/k_max + S'(psi)
  //   dpsi_stem/dpsi   = P'(E_psi_stem) * dE_psi_stem/dpsi.
  // E_up_'(r) is analytic (dE_from_soil_dpsi_collar); near a branch kink it
  // returns NaN and we fall back to the central difference on the transport.
  const double r = -psi;
  const double dEup_dr = dE_from_soil_dpsi_collar(r, psi_soil_inverted_);
  double dpsistem_dpsi;
  if (std::isfinite(dEup_dr)) {
    E_from_Soil_to_Root_Collar(r, psi_soil_inverted_);  // refresh E_up_ at r
    const double E_psi_stem =
        E_up_ / leaf_specific_conductance_max_ + transpiration_from_psi.eval(psi);
    const double dEpsistem_dpsi =
        -dEup_dr / leaf_specific_conductance_max_ + transpiration_from_psi.deriv(psi);
    dpsistem_dpsi = psi_from_transpiration.deriv(E_psi_stem) * dEpsistem_dpsi;
  } else {
    const double h = 1e-6;
    dpsistem_dpsi =
        (find_psi_stem_from_psi_root(-(psi + h), psi_soil_inverted_) -
         find_psi_stem_from_psi_root(-(psi - h), psi_soil_inverted_)) / (2.0 * h);
  }

  const double dci_dpsi = dci_dpsistem * dpsistem_dpsi + dci_dpsi_expl;
  return A_prime * dci_dpsi - C_prime * dpsistem_dpsi;
}

// Analytic d(E_up_)/d(P_x_r): the signed-collar-potential derivative of the
// soil->root-collar uptake, mirroring the general branch of
// E_from_Soil_to_Root_Collar. Per layer, with span = |psi_soil[i] - P_x_r| and
// integral = \int f_r over [P_src_min, P_src_max] (the cumulative-vulnerability
// curve root_vuln_integral_from_psi, whose integrand is root_vuln_from_psi):
//   E_i        = (psi_soil[i] - P_x_r - grav) / area_leaf / r_R,
//   r_R        = r_R_H_min[i] * span / integral + r_R_V_sum[i],
//   dspan/dP   = sign_var   (+1 if P_x_r is the upper bound, else -1),
//   dinteg/dP  = sign_var * f_r(-P_x_r)  for P_x_r<0  (else sign_var, f_r==1),
// and dE_i/dP follows by the quotient rule. Returns NaN on any branch kink so
// the caller falls back to finite differences.
double Leaf::dE_from_soil_dpsi_collar(double P_x_r, const std::vector<double>& psi_soil) {
  const double inv_area_leaf = 1.0 / area_leaf_;
  const double kink_tol = 1e-8;
  double dEup_dr_mol = 0.0;

  for (int i = 0; i < max_soil_layer; i++) {
    // Branch kinks: equal potentials, gravity-balance, and the psi==0 split of
    // the vulnerability integral. The analytic general-branch derivative is not
    // valid across these, so signal a fallback.
    if (std::abs(P_x_r - psi_soil[i]) < kink_tol ||
        std::abs((psi_soil[i] - P_x_r) - grav_head_z_[i]) < kink_tol ||
        std::abs(P_x_r) < kink_tol) {
      return std::numeric_limits<double>::quiet_NaN();
    }

    const double P_src_min = std::min(psi_soil[i], P_x_r);
    const double P_src_max = std::max(psi_soil[i], P_x_r);
    const double span = P_src_max - P_src_min;
    const double sign_var = (P_x_r > psi_soil[i]) ? 1.0 : -1.0;  // = dspan/dP_x_r

    // integral, replicated bit-for-bit from E_from_Soil_to_Root_Collar.
    const double hi_neg = std::min(P_src_max, 0.0);
    const double lo_pos = std::max(P_src_min, 0.0);
    double integral = 0.0;
    if (hi_neg > P_src_min) {
      integral += root_vuln_integral_from_psi.eval(-P_src_min) -
                  root_vuln_integral_from_psi.eval(-hi_neg);
    }
    if (P_src_max > lo_pos) {
      integral += (P_src_max - lo_pos);
    }

    // d(integral)/d(P_x_r): for P_x_r<0 the moving bound is in the vulnerable
    // region. The integrand is the derivative of the *same* cumulative spline
    // that produced `integral` (root_vuln_integral_from_psi.deriv), NOT the
    // separate root_vuln_from_psi spline: the two agree on the knot domain but
    // extrapolate independently (both clamp-to-last-value, #527), so beyond the
    // domain only the integral spline's own derivative stays consistent with its
    // value. For P_x_r>0 the moving bound is in the above-atmospheric part
    // (f_r==1), contributed linearly, so the slope is 1.
    const double fr_at =
        (P_x_r < 0.0) ? root_vuln_integral_from_psi.deriv(-P_x_r) : 1.0;
    const double dinteg_dr = sign_var * fr_at;

    const double r_R_H = r_R_H_min[i] * span / integral;
    const double r_R = r_R_H + r_R_V_sum[i];
    const double dr_R_H_dr =
        r_R_H_min[i] * (sign_var * integral - span * dinteg_dr) / (integral * integral);
    const double dr_R_dr = dr_R_H_dr;

    const double num = (psi_soil[i] - P_x_r - grav_head_z_[i]) * inv_area_leaf;
    const double dnum_dr = -inv_area_leaf;
    // E_i = num / r_R  ->  quotient rule.
    dEup_dr_mol += (dnum_dr * r_R - num * dr_R_dr) / (r_R * r_R);
  }

  return dEup_dr_mol * kg_per_mol_h2o;  // match E_up_'s kg units
}

// The general branch of E_from_Soil_to_Root_Collar carried one derivative
// further than dE_from_soil_dpsi_collar takes it, on shared intermediates.
void Leaf::layer_flux_partials(double P_x_r, const std::vector<double>& psi_soil,
                               LayerFlux& out) {
  const double nan = std::numeric_limits<double>::quiet_NaN();
  out.E.assign(max_soil_layer, nan);
  out.dE_dr.assign(max_soil_layer, nan);
  out.dE_dpsi.assign(max_soil_layer, nan);
  out.d2E_dr_dpsi.assign(max_soil_layer, nan);
  out.r_R.assign(max_soil_layer, nan);
  out.dr_R_dr.assign(max_soil_layer, nan);
  out.num.assign(max_soil_layer, nan);
  out.integral.assign(max_soil_layer, nan);

  const double inv_area_leaf = 1.0 / area_leaf_;
  const double kink_tol = 1e-8;

  for (int i = 0; i < max_soil_layer; i++) {
    const double x = psi_soil[i];
    if (std::abs(P_x_r - x) < kink_tol ||
        std::abs((x - P_x_r) - grav_head_z_[i]) < kink_tol ||
        std::abs(P_x_r) < kink_tol) {
      return;  // a branch kink: every entry stays NaN
    }

    const double P_src_min = std::min(x, P_x_r);
    const double P_src_max = std::max(x, P_x_r);
    const double span = P_src_max - P_src_min;
    const double sv = (P_x_r > x) ? 1.0 : -1.0;  // = dspan/dP_x_r

    const double hi_neg = std::min(P_src_max, 0.0);
    const double lo_pos = std::max(P_src_min, 0.0);
    double integral = 0.0;
    if (hi_neg > P_src_min) {
      integral += root_vuln_integral_from_psi.eval(-P_src_min) -
                  root_vuln_integral_from_psi.eval(-hi_neg);
    }
    if (P_src_max > lo_pos) {
      integral += (P_src_max - lo_pos);
    }

    // Slopes of the same cumulative spline at the two endpoints; the collar's
    // endpoint moves with P_x_r and the layer's with x.
    const double fr_r = (P_x_r < 0.0)
                            ? root_vuln_integral_from_psi.deriv(-P_x_r) : 1.0;
    const double fr_x = (x < 0.0)
                            ? root_vuln_integral_from_psi.deriv(-x) : 1.0;
    const double dI_dr = sv * fr_r;
    const double dI_dx = -sv * fr_x;

    const double h = r_R_H_min[i];
    const double integral_sq = integral * integral;
    const double r_R = h * span / integral + r_R_V_sum[i];
    const double U = sv * integral - span * dI_dr;
    const double dr_R_dr = h * U / integral_sq;
    const double dr_R_dx = h * (-sv * integral - span * dI_dx) / integral_sq;
    const double dU_dx = sv * dI_dx + fr_r;
    const double d_dr_R_dr_dx =
        h * (dU_dx * integral - 2.0 * U * dI_dx) / (integral_sq * integral);

    const double num = (x - P_x_r - grav_head_z_[i]) * inv_area_leaf;
    const double dnum_dr = -inv_area_leaf;
    const double dnum_dx = inv_area_leaf;
    const double N = dnum_dr * r_R - num * dr_R_dr;  // numerator of dE/dr
    const double dN_dx =
        dnum_dr * dr_R_dx - dnum_dx * dr_R_dr - num * d_dr_R_dr_dx;

    out.E[i] = num / r_R;
    out.dE_dr[i] = N / (r_R * r_R);
    out.dE_dpsi[i] = (dnum_dx * r_R - num * dr_R_dx) / (r_R * r_R);
    out.d2E_dr_dpsi[i] = (dN_dx * r_R - 2.0 * N * dr_R_dx) / (r_R * r_R * r_R);
    out.r_R[i] = r_R;
    out.dr_R_dr[i] = dr_R_dr;
    out.num[i] = num;
    out.integral[i] = integral;
  }
}

double Leaf::dR_dcollar_at(double opt_root_psi, double h) {
  const double R_plus = dprofit_droot_collar_psi(opt_root_psi + h);
  const double R_minus = dprofit_droot_collar_psi(opt_root_psi - h);
  return (R_plus - R_minus) / (2.0 * h);
}

// Under a uniform drying every potential difference and the gravitational head
// are unchanged, leaving the vulnerability integral over a sliding interval.
void Leaf::translation_partials(std::vector<double>& dE_dd, double& dpsistem_dd) {
  const double p = -root_collar_psi_;
  const double r = root_collar_psi_;
  refresh_soil_potentials();

  LayerFlux fx;
  layer_flux_partials(r, psi_soil_inverted_, fx);

  dE_dd.assign(max_soil_layer, std::numeric_limits<double>::quiet_NaN());
  double dEup_dd_mol = 0.0;
  for (int i = 0; i < max_soil_layer; ++i) {
    const double x = psi_soil_inverted_[i];
    const double P_src_min = std::min(x, r);
    const double P_src_max = std::max(x, r);
    const double fr_min = (P_src_min < 0.0)
        ? root_vuln_integral_from_psi.deriv(-P_src_min) : 1.0;
    const double fr_max = (P_src_max < 0.0)
        ? root_vuln_integral_from_psi.deriv(-P_src_max) : 1.0;
    const double r_R_H = fx.r_R[i] - r_R_V_sum[i];
    dE_dd[i] = fx.E[i] * (r_R_H / fx.r_R[i]) * (fr_min - fr_max) / fx.integral[i];
    dEup_dd_mol += dE_dd[i];
  }
  const double dEup_dd = dEup_dd_mol * kg_per_mol_h2o;

  E_from_Soil_to_Root_Collar(r, psi_soil_inverted_);
  const double S_p = transpiration_from_psi.eval(p);
  const double E_psi_stem = E_up_ / leaf_specific_conductance_max_ + S_p;
  const double P_prime = psi_from_transpiration.deriv(E_psi_stem);
  const double S_prime = transpiration_from_psi.deriv(p);
  // psi_from_transpiration inverts transpiration_from_psi, so the leading 1 is
  // exact and the bracket is a difference of two slopes of the same pair.
  dpsistem_dd = 1.0 + P_prime * dEup_dd / leaf_specific_conductance_max_ +
                S_prime * (P_prime - psi_from_transpiration.deriv(S_p));
}

double Leaf::dR_dflux_from_layer(int layer, double delta,
                                 double dR_dflux_slope) {
  const double p = -root_collar_psi_;
  const double r = root_collar_psi_;
  const double psi_soil_layer = psi_soil_[layer];

  double R[2], E_up[2], slope[2];
  for (int k = 0; k < 2; ++k) {
    psi_soil_[layer] = psi_soil_layer + (k == 0 ? delta : -delta);
    transpiration_cached_ = false;
    R[k] = dprofit_droot_collar_psi(p);
    refresh_soil_potentials();
    E_from_Soil_to_Root_Collar(r, psi_soil_inverted_);
    E_up[k] = E_up_;
    slope[k] = dE_from_soil_dpsi_collar(r, psi_soil_inverted_);
  }
  psi_soil_[layer] = psi_soil_layer;
  transpiration_cached_ = false;
  refresh_soil_potentials();

  const double dR = (R[0] - R[1]) / (2.0 * delta);
  const double dE_up = (E_up[0] - E_up[1]) / (2.0 * delta);
  const double dslope = (slope[0] - slope[1]) / (2.0 * delta);
  return (dR - dR_dflux_slope * dslope) / dE_up;
}

// The leaf's parameter inputs, after the state-dependent block, each beside the
// member input_adjoints writes while it differentiates that row. beta_R_H and
// beta_R_V are the only multiplicative scale on the root resistance network and
// have no row here, so a strategy varying either reads exactly zero.
static const char* const leaf_parameter_inputs[] = {
    "vcmax_25", "jmax_25", "a",      "curv_fact_elec_trans",
    "curv_fact_colim", "b", "c",     "psi_crit",
    "beta2",    "g1_TF24", "rho",    "a_bio",
    "root_b",   "root_c",  "root_psi_crit"};
static double Leaf::*const leaf_parameter_slots[] = {
    &Leaf::vcmax_25, &Leaf::jmax_25, &Leaf::a,       &Leaf::curv_fact_elec_trans,
    &Leaf::curv_fact_colim, &Leaf::b, &Leaf::c,      &Leaf::psi_crit,
    &Leaf::beta2,    &Leaf::g1_TF24, &Leaf::rho_,    &Leaf::a_bio_,
    &Leaf::root_b,   &Leaf::root_c,  &Leaf::root_psi_crit};
static const int n_leaf_parameter_inputs = 15;
// Index into both tables above.
enum LeafParameter {
  PAR_VCMAX_25 = 0, PAR_JMAX_25, PAR_A, PAR_CURV_ELEC, PAR_CURV_COLIM,
  PAR_B, PAR_C, PAR_PSI_CRIT, PAR_BETA2, PAR_G1, PAR_RHO, PAR_A_BIO,
  PAR_ROOT_B, PAR_ROOT_C, PAR_ROOT_PSI_CRIT};

// Where the operating point is a bound of the feasible interval, dp*/du is that
// bound's derivative. Each bound is the root of a residual in the collar
// potential, so each row is the implicit function theorem on it, at the bound.
void Leaf::bound_partials(std::vector<double>& out) {
  const int n = max_soil_layer;
  out.assign(inputs().size(), 0.0);
  const int i_psi0 = 1, i_area = 1 + n, i_mass0 = 2 + n, i_kappa = 2 + 2 * n,
            i_par0 = 3 + 2 * n;

  double bound_a = 0.0, bound_b = 0.0;
  if (!prepare_collar_solve(bound_a, bound_b)) {
    util::stop("bound_partials: the collar potential is pinned to an interval "
               "that prepare_collar_solve no longer produces");
  }
  const double p = -root_collar_psi_;
  const bool at_bound_a = (p - bound_a) < (bound_b - p);
  const double p_bound = at_bound_a ? bound_a : bound_b;

  // bound_b takes the drier of the stem's and the root's ceilings. Where the
  // root's wins it is an input in its own right and nothing else moves it.
  if (!at_bound_a && p_bound == -root_psi_crit) {
    out[i_par0 + PAR_ROOT_PSI_CRIT] = -1.0;
    return;
  }

  const double x = -p_bound;  // the same potential, signed
  const double kappa = leaf_specific_conductance_max_;
  refresh_soil_potentials();
  LayerFlux fx;
  layer_flux_partials(x, psi_soil_inverted_, fx);

  double dEup_dx = 0.0, E_up_here = 0.0;
  std::vector<double> dEup_dm(n, 0.0);
  for (int i = 0; i < n; ++i) {
    dEup_dx += fx.dE_dr[i];
    E_up_here += fx.E[i];
  }
  dEup_dx *= kg_per_mol_h2o;
  E_up_here *= kg_per_mol_h2o;
  for (int j = 0; j < n; ++j) {
    const double m_j = 3.0 * c_r_V_[j];
    for (int i = j; i < n; ++i) {
      const double m_i = 3.0 * c_r_V_[i];
      const double r_R = fx.r_R[i];
      const double q = -r_R_V[j] / m_j +
                       (i == j ? -(r_R - r_R_V_sum[i]) / m_i : 0.0);
      dEup_dm[j] += -fx.num[i] * q / (r_R * r_R);
    }
    dEup_dm[j] *= kg_per_mol_h2o;
  }

  // bound_a is where uptake is zero; bound_b is where the stem reaches psi_crit,
  // which adds the supply-side term. Both residuals are in kg.
  const double S_bound = transpiration_from_psi.eval(p_bound);
  const double S_crit = transpiration_from_psi.eval(psi_crit);
  double dF_dx = dEup_dx;
  if (!at_bound_a) {
    dF_dx -= kappa * transpiration_from_psi.deriv(p_bound);
    out[i_kappa] = -(S_crit - S_bound);
    out[i_par0 + PAR_PSI_CRIT] = -kappa * transpiration_from_psi.deriv(psi_crit);
  }
  for (int j = 0; j < n; ++j) {
    out[i_psi0 + j] = -kg_per_mol_h2o * fx.dE_dpsi[j];
    out[i_mass0 + j] = dEup_dm[j];
  }
  out[i_area] = -E_up_here / area_leaf_;

  // The interpolant parameters, on knot grids held still while the parameter
  // moves (set_transpiration_at, set_root_vulnerability_at).
  std::vector<double> knots_stem, knots_root, knot_values;
  build_cumulative_vulnerability_integral(b, c, vulnerability_curve_ncontrol,
                                          knots_stem, knot_values);
  build_cumulative_vulnerability_integral(root_b, root_c,
                                          vulnerability_curve_ncontrol,
                                          knots_root, knot_values);
  const int transport_pars[4] = {PAR_B, PAR_C, PAR_ROOT_B, PAR_ROOT_C};
  for (int t = 0; t < 4; ++t) {
    const int k = transport_pars[t];
    const bool root = (k == PAR_ROOT_B || k == PAR_ROOT_C);
    if (!par_wanted(k)) {
      continue;  // SPIKE: row not requested; input_adjoints poisons it to NaN
    }
    if (at_bound_a && !root) {
      continue;  // the supply-side term the stem parameters reach is bound_b's
    }
    const double keep = this->*leaf_parameter_slots[k];
    const double h = std::abs(keep) * 1e-6;
    if (!(h > 0.0)) {
      continue;
    }
    double F[2];
    for (int side = 0; side < 2; ++side) {
      this->*leaf_parameter_slots[k] = keep + (side == 0 ? h : -h);
      if (root) {
        set_root_vulnerability_at(root_b, root_c, knots_root);
      } else {
        set_transpiration_at(b, c, knots_stem);
      }
      refresh_soil_potentials();
      E_from_Soil_to_Root_Collar(x, psi_soil_inverted_);
      F[side] = E_up_;
      if (!at_bound_a) {
        F[side] -= kappa * (transpiration_from_psi.eval(psi_crit) -
                            transpiration_from_psi.eval(p_bound));
      }
    }
    this->*leaf_parameter_slots[k] = keep;
    if (root) {
      set_root_vulnerability_at(root_b, root_c, knots_root);
    } else {
      set_transpiration_at(b, c, knots_stem);
    }
    out[i_par0 + k] = (F[0] - F[1]) / (2.0 * h);
  }
  refresh_soil_potentials();

  // The theorem gives dx/du = -(dF/du)/(dF/dx) and the bound is -x, so the two
  // sign flips cancel.
  for (size_t k = 0; k < out.size(); ++k) {
    out[k] /= dF_dx;
  }
}

std::vector<std::string> Leaf::inputs() const {
  std::vector<std::string> out;
  out.reserve(2 * max_soil_layer + 3 + n_leaf_parameter_inputs);
  out.push_back("PPFD");
  for (int i = 0; i < max_soil_layer; ++i) {
    out.push_back("psi_soil[" + std::to_string(i) + "]");
  }
  out.push_back("area_leaf");
  for (int i = 0; i < max_soil_layer; ++i) {
    out.push_back("mass_root[" + std::to_string(i) + "]");
  }
  out.push_back("leaf_specific_conductance_max");
  for (int i = 0; i < n_leaf_parameter_inputs; ++i) {
    out.push_back(leaf_parameter_inputs[i]);
  }
  return out;
}

void Leaf::input_adjoints(double lambda_profit,
                          const std::vector<double>& lambda_uptake,
                          std::vector<double>& input_adjoints) {
  const int n = max_soil_layer;
  if (static_cast<int>(lambda_uptake.size()) != n) {
    util::stop("input_adjoints: one uptake adjoint per rooted layer; got " +
               std::to_string(lambda_uptake.size()) + " for max_soil_layer=" +
               std::to_string(n));
  }
  const std::vector<std::string> names = inputs();
  input_adjoints.assign(names.size(), 0.0);

  // Every evaluator below writes the operating-point outputs, so hold them and
  // put them back: a caller's forward values must survive this call.
  const double keep_ci = ci_, keep_stom = stom_cond_CO2_,
               keep_assim = assim_colimited_, keep_transpiration = transpiration_,
               keep_profit = profit_, keep_cost = hydraulic_cost_,
               keep_E_up = E_up_, keep_stem = opt_psi_stem_,
               keep_collar = root_collar_psi_;
  const std::vector<double> keep_consumption = soil_consumption_;
  auto restore_operating_point = [&]() {
    transpiration_cached_ = false;
    refresh_soil_potentials();
    ci_ = keep_ci;
    stom_cond_CO2_ = keep_stom;
    assim_colimited_ = keep_assim;
    transpiration_ = keep_transpiration;
    profit_ = keep_profit;
    hydraulic_cost_ = keep_cost;
    E_up_ = keep_E_up;
    opt_psi_stem_ = keep_stem;
    root_collar_psi_ = keep_collar;
    soil_consumption_ = keep_consumption;
  };

  const int i_ppfd = 0, i_psi0 = 1, i_area = 1 + n, i_mass0 = 2 + n,
            i_kappa = 2 + 2 * n, i_par0 = 3 + 2 * n;

  // SPIKE (p3/trait-mask): a leaf parameter row the caller did not ask for is
  // left ABSENT, not zero. Whatever the branches below happened to write into
  // it is replaced by a quiet NaN on every exit path, so a masked row can never
  // be mistaken for a computed derivative that came out at zero.
  auto poison_unwanted = [&]() -> void {
    for (int k = 0; k < n_leaf_parameter_inputs; ++k) {
      if (!par_wanted(k)) {
        input_adjoints[i_par0 + k] = std::numeric_limits<double>::quiet_NaN();
      }
    }
  };

  const double p = -keep_collar;   // collar potential, positive magnitude
  const double r = keep_collar;    // the same potential, signed
  const double psi_stem = keep_stem;
  const double kappa = leaf_specific_conductance_max_;
  const double gstar_Pa = gamma_ * umol_per_mol_to_Pa;

  refresh_soil_potentials();
  const double ci = psi_stem_to_ci(psi_stem, p);

  const double A_prime =
      odelia::ode::forward_derivative(ci, [&](auto x) -> decltype(x) {
        return assim_colimited_ad<decltype(x)>(x, vcmax_, electron_transport_,
                                               gstar_Pa, km_, R_d_,
                                               curv_fact_colim);
      });
  const double C_prime =
      odelia::ode::forward_derivative(psi_stem, [&](auto x) -> decltype(x) {
        return hydraulic_cost_ad<decltype(x)>(x, b, c, g1_TF24, beta2);
      });

  const double gc_const =
      atm_kpa_ * kg_to_mol_h2o / atm_vpd_ / H2O_CO2_stom_diff_ratio;
  const double gc = gc_const * transpiration(psi_stem, p);
  const double dgc_dpsistem =
      gc_const * kappa * transpiration_from_psi.deriv(psi_stem);
  const double inv_atm = 1.0 / (atm_kpa_ * kPa_to_Pa);
  const double g_ci = A_prime * umol_to_mol + gc * inv_atm;
  const double dci_dpsistem = dgc_dpsistem * (ca_ - ci) * inv_atm / g_ci;

  E_from_Soil_to_Root_Collar(r, psi_soil_inverted_);
  const double E_up = E_up_;  // kg H2O m^-2 leaf s^-1
  const double E_psi_stem = E_up / kappa + transpiration_from_psi.eval(p);
  const double P_prime = psi_from_transpiration.deriv(E_psi_stem);
  const double dprofit_dpsistem = A_prime * dci_dpsistem - C_prime;

  LayerFlux fx;
  layer_flux_partials(r, psi_soil_inverted_, fx);

  // Root mass enters layer i's resistance through its own horizontal term and
  // every shallower vertical term, so the columns sum over i >= j.
  std::vector<double> dEup_dm(n, 0.0), dslope_dm(n, 0.0), uptake_dm(n, 0.0);
  for (int j = 0; j < n; ++j) {
    const double m_j = 3.0 * c_r_V_[j];
    for (int i = j; i < n; ++i) {
      const double m_i = 3.0 * c_r_V_[i];
      const double r_R = fx.r_R[i];
      const double q = -r_R_V[j] / m_j +
                       (i == j ? -(r_R - r_R_V_sum[i]) / m_i : 0.0);
      const double d_dr_R_dr_dm = (i == j) ? -fx.dr_R_dr[i] / m_i : 0.0;
      const double N = fx.dE_dr[i] * r_R * r_R;
      const double dN_dm = (-1.0 / area_leaf_) * q - fx.num[i] * d_dr_R_dr_dm;
      const double dE_dm = -fx.num[i] * q / (r_R * r_R);
      const double d2E_dm = (dN_dm * r_R - 2.0 * N * q) / (r_R * r_R * r_R);
      dEup_dm[j] += dE_dm;
      dslope_dm[j] += d2E_dm;
      uptake_dm[j] += lambda_uptake[i] * dE_dm;
    }
    dEup_dm[j] *= kg_per_mol_h2o;
    dslope_dm[j] *= kg_per_mol_h2o;
  }

  // profit_ sits at its own maximiser, so its sensitivity to every input is the
  // direct partial with the collar potential held still.
  const double phi = (ci - gstar_Pa) / (ci + 2.0 * gstar_Pa);
  const double ar = vcmax_ * (ci - gstar_Pa) / (ci + km_);
  const double ae = electron_transport_ / 4.0 * phi;
  const double sum_a = ar + ae;
  const double disc = sum_a * sum_a - 4.0 * curv_fact_colim * ar * ae;
  const double dA_dae =
      (1.0 - (2.0 * sum_a - 4.0 * curv_fact_colim * ar) / (2.0 * sqrt(disc))) /
      (2.0 * curv_fact_colim);
  const double aP = a * PPFD_;
  const double disc_J =
      (aP + jmax_) * (aP + jmax_) - 4.0 * curv_fact_elec_trans * aP * jmax_;
  const double dJ_dPPFD =
      (a - (a * (aP + jmax_) - 2.0 * curv_fact_elec_trans * a * jmax_) /
               sqrt(disc_J)) /
      (2.0 * curv_fact_elec_trans);
  // Radiation moves ci as well as assimilation; the two combine to the direct
  // partial times the supply share of the ci residual's slope.
  const double dprofit_dPPFD =
      dA_dae * phi / 4.0 * dJ_dPPFD * gc * inv_atm / g_ci;
  const double dci_dkappa = (gc / kappa) * (ca_ - ci) * inv_atm / g_ci;

  input_adjoints[i_ppfd] += lambda_profit * dprofit_dPPFD;
  for (int j = 0; j < n; ++j) {
    input_adjoints[i_psi0 + j] += lambda_profit * dprofit_dpsistem * P_prime *
                                  (-kg_per_mol_h2o * fx.dE_dpsi[j]) / kappa;
    input_adjoints[i_mass0 + j] +=
        lambda_profit * dprofit_dpsistem * P_prime * dEup_dm[j] / kappa;
  }
  input_adjoints[i_area] +=
      lambda_profit * dprofit_dpsistem * P_prime * (-E_up / area_leaf_) / kappa;
  input_adjoints[i_kappa] +=
      lambda_profit * (dprofit_dpsistem * P_prime * (-E_up / (kappa * kappa)) +
                       A_prime * dci_dkappa);

  // Uptake reads its own layer's potential and the root masses at or above it;
  // 1/area_leaf is a single factor and the resistances carry no leaf area.
  for (int j = 0; j < n; ++j) {
    input_adjoints[i_psi0 + j] += lambda_uptake[j] * (-fx.dE_dpsi[j]);
    input_adjoints[i_mass0 + j] += uptake_dm[j];
    input_adjoints[i_area] += lambda_uptake[j] * (-fx.E[j] / area_leaf_);
  }

  // Uptake consumes the collar potential, so the flux adjoints collapse onto one
  // scalar; profit adds nothing to it at an interior maximum, where dprofit/dp
  // is zero, and does at a bound, where it is not.
  double s_adjoint = 0.0;
  for (int j = 0; j < n; ++j) {
    s_adjoint += lambda_uptake[j] * (-fx.dE_dr[j]);
  }

  // The parameter rows. Photosynthesis and cost reach profit only through the
  // two templated functions, so those partials are exact; the four that rebuild
  // an interpolant, and every parameter's partial of R, are central differences
  // of the leaf's own evaluators at the collar potential the solve left.
  //
  // R is itself a derivative, so its parameter partial is a second derivative
  // the templated pair does not carry -- the same substitution radiation's and
  // the conductance's rows above already make.
  //
  // The knot grids are captured here and held still across every perturbation.
  // The builder's grid ends at b*log(100)^(1/c) and its knot COUNT steps as b or
  // c moves; a difference across that step reads 47x the derivative of
  // d(R)/d(root_b) at a 1e-6 step and nothing announces it.
  std::vector<double> knots_stem, knots_root, knot_values;
  build_cumulative_vulnerability_integral(b, c, vulnerability_curve_ncontrol,
                                          knots_stem, knot_values);
  build_cumulative_vulnerability_integral(root_b, root_c,
                                          vulnerability_curve_ncontrol,
                                          knots_root, knot_values);
  // Writes parameter k and everything the leaf derives from it, so every
  // evaluator below reads a consistent leaf.
  auto set_parameter = [&](int k, double v) {
    this->*leaf_parameter_slots[k] = v;
    switch (k) {
    case PAR_VCMAX_25:
      vcmax_ = peak_arrh_curve(vcmax_ha, vcmax_25, leaf_temp_, vcmax_H_d,
                               vcmax_d_S);
      R_d_ = vcmax_ * 0.015;
      break;
    case PAR_JMAX_25:
      jmax_ = peak_arrh_curve(jmax_ha, jmax_25, leaf_temp_, jmax_H_d, jmax_d_S);
      electron_transport_ = electron_transport();
      break;
    case PAR_A:
    case PAR_CURV_ELEC:
      electron_transport_ = electron_transport();
      break;
    case PAR_B:
    case PAR_C:
      set_transpiration_at(b, c, knots_stem);
      break;
    case PAR_ROOT_B:
    case PAR_ROOT_C:
      set_root_vulnerability_at(root_b, root_c, knots_root);
      break;
    default:
      break;
    }
    transpiration_cached_ = false;
  };
  // psi_crit and root_psi_crit set the feasible interval and appear in no
  // evaluator inside it; rho and a_bio are stored by set_physiology and read
  // nowhere, as sapwood_volume_per_leaf_area is.
  auto reaches_operating_point = [](int k) -> bool {
    return k != PAR_PSI_CRIT && k != PAR_ROOT_PSI_CRIT && k != PAR_RHO &&
           k != PAR_A_BIO;
  };
  auto rebuilds_transport = [](int k) -> bool {
    return k == PAR_B || k == PAR_C || k == PAR_ROOT_B || k == PAR_ROOT_C;
  };

  const double supply_share = gc * inv_atm / g_ci;
  std::vector<double> dprofit_dpar(n_leaf_parameter_inputs, 0.0);
  std::vector<double> dR_dpar(n_leaf_parameter_inputs, 0.0);
  std::vector<std::vector<double>> dE_dpar(n_leaf_parameter_inputs);

  // A parameter of the assimilation algebra also moves ci, and the two combine
  // to the direct partial times the supply share of the ci residual's slope.
  dprofit_dpar[PAR_VCMAX_25] =
      supply_share *
      odelia::ode::forward_derivative(vcmax_25, [&](auto x) -> decltype(x) {
        const decltype(x) vc = x * (vcmax_ / vcmax_25);
        return assim_colimited_ad<decltype(x)>(ci, vc, electron_transport_,
                                               gstar_Pa, km_, vc * 0.015,
                                               curv_fact_colim);
      });
  dprofit_dpar[PAR_JMAX_25] =
      supply_share *
      odelia::ode::forward_derivative(jmax_25, [&](auto x) -> decltype(x) {
        const decltype(x) jm = x * (jmax_ / jmax_25);
        const decltype(x) et = electron_transport_ad<decltype(x)>(
            a, PPFD_, jm, curv_fact_elec_trans);
        return assim_colimited_ad<decltype(x)>(ci, vcmax_, et, gstar_Pa, km_,
                                               R_d_, curv_fact_colim);
      });
  dprofit_dpar[PAR_A] =
      supply_share *
      odelia::ode::forward_derivative(a, [&](auto x) -> decltype(x) {
        const decltype(x) et = electron_transport_ad<decltype(x)>(
            x, PPFD_, jmax_, curv_fact_elec_trans);
        return assim_colimited_ad<decltype(x)>(ci, vcmax_, et, gstar_Pa, km_,
                                               R_d_, curv_fact_colim);
      });
  dprofit_dpar[PAR_CURV_ELEC] =
      supply_share *
      odelia::ode::forward_derivative(
          curv_fact_elec_trans, [&](auto x) -> decltype(x) {
            const decltype(x) et =
                electron_transport_ad<decltype(x)>(a, PPFD_, jmax_, x);
            return assim_colimited_ad<decltype(x)>(ci, vcmax_, et, gstar_Pa,
                                                   km_, R_d_, curv_fact_colim);
          });
  dprofit_dpar[PAR_CURV_COLIM] =
      supply_share *
      odelia::ode::forward_derivative(
          curv_fact_colim, [&](auto x) -> decltype(x) {
            return assim_colimited_ad<decltype(x)>(
                ci, vcmax_, electron_transport_, gstar_Pa, km_, R_d_, x);
          });
  // The cost is outside the ci residual, so it reaches profit and nothing else.
  dprofit_dpar[PAR_G1] =
      -odelia::ode::forward_derivative(g1_TF24, [&](auto x) -> decltype(x) {
        return hydraulic_cost_ad<decltype(x)>(psi_stem, b, c, x, beta2);
      });
  dprofit_dpar[PAR_BETA2] =
      -odelia::ode::forward_derivative(beta2, [&](auto x) -> decltype(x) {
        return hydraulic_cost_ad<decltype(x)>(psi_stem, b, c, g1_TF24, x);
      });

  for (int k = 0; k < n_leaf_parameter_inputs; ++k) {
    if (!reaches_operating_point(k)) {
      continue;
    }
    if (!par_wanted(k)) {
      continue;  // SPIKE: row not requested; poisoned to NaN below
    }
    const double keep = this->*leaf_parameter_slots[k];
    const double h = std::abs(keep) * 1e-6;
    if (!(h > 0.0)) {
      continue;  // a parameter at zero has no scale to step on
    }
    double R_pm[2], profit_pm[2] = {0.0, 0.0};
    std::vector<double> uptake_pm[2];
    for (int side = 0; side < 2; ++side) {
      set_parameter(k, keep + (side == 0 ? h : -h));
      R_pm[side] = dprofit_droot_collar_psi(p);
      if (rebuilds_transport(k)) {
        refresh_soil_potentials();
        profit_pm[side] = profit_psi_stem_TF(
            find_psi_stem_from_psi_root(r, psi_soil_inverted_), p);
        E_from_Soil_to_Root_Collar(r, psi_soil_inverted_);
        uptake_pm[side] = soil_consumption_;
      }
    }
    set_parameter(k, keep);
    dR_dpar[k] = (R_pm[0] - R_pm[1]) / (2.0 * h);
    if (rebuilds_transport(k)) {
      dprofit_dpar[k] = (profit_pm[0] - profit_pm[1]) / (2.0 * h);
      dE_dpar[k].assign(n, 0.0);
      for (int j = 0; j < n; ++j) {
        dE_dpar[k][j] = (uptake_pm[0][j] - uptake_pm[1][j]) / (2.0 * h);
      }
    }
  }
  refresh_soil_potentials();

  for (int k = 0; k < n_leaf_parameter_inputs; ++k) {
    input_adjoints[i_par0 + k] = lambda_profit * dprofit_dpar[k];
    for (int j = 0; j < n && !dE_dpar[k].empty(); ++j) {
      input_adjoints[i_par0 + k] += lambda_uptake[j] * dE_dpar[k][j];
    }
  }

  if (collar_pinned_) {
    // p* is the bound, so dp*/du is the bound's derivative and none of the
    // interior argmax machinery applies. profit is not stationary at a bound,
    // so it carries the operating point's motion here where the interior
    // branch's envelope argument leaves it out.
    const double w = lambda_profit * dprofit_droot_collar_psi(p) + s_adjoint;
    std::vector<double> dbound;
    bound_partials(dbound);
    collar_pinned_ = true;  // prepare_collar_solve inside bound_partials clears it
    for (size_t k = 0; k < dbound.size(); ++k) {
      input_adjoints[k] += w * dbound[k];
    }
    poison_unwanted();
    restore_operating_point();
    return;
  }

  // Taken here rather than kept from the solve: the Newton loop's divisor is
  // allowed to lag an iterate, which is a factor of about 1.02 on every row.
  const double Pi_pp = dR_dcollar_at(p, 1e-6);
  const double mu = -s_adjoint / Pi_pp;
  for (int k = 0; k < n_leaf_parameter_inputs; ++k) {
    if (!par_wanted(k)) {
      continue;
    }
    input_adjoints[i_par0 + k] += mu * dR_dpar[k];
  }

  // R reads the potentials, root masses and leaf area only through the flux and
  // its collar slope, so those 2n+1 directions cost these two scalars.
  const double dR_dflux_slope = -dprofit_dpsistem * P_prime / kappa;
  const double dR_dflux = dR_dflux_from_layer(0, 1e-6, dR_dflux_slope);
  double dEup_dr = 0.0;
  for (int j = 0; j < n; ++j) {
    dEup_dr += fx.dE_dr[j];
  }
  dEup_dr *= kg_per_mol_h2o;

  for (int j = 0; j < n; ++j) {
    input_adjoints[i_psi0 + j] +=
        mu * (dR_dflux * (-kg_per_mol_h2o * fx.dE_dpsi[j]) +
              dR_dflux_slope * (-kg_per_mol_h2o * fx.d2E_dr_dpsi[j]));
    input_adjoints[i_mass0 + j] +=
        mu * (dR_dflux * dEup_dm[j] + dR_dflux_slope * dslope_dm[j]);
  }
  input_adjoints[i_area] +=
      mu * (dR_dflux * (-E_up / area_leaf_) +
            dR_dflux_slope * (-dEup_dr / area_leaf_));

  // Radiation and the conductance reach R outside the flux pair, so each takes
  // its own residual pair.
  {
    const double h = PPFD_ * 1e-6;
    double R_pm[2];
    for (int k = 0; k < 2; ++k) {
      PPFD_ += (k == 0 ? h : -2.0 * h);
      electron_transport_ = electron_transport();
      transpiration_cached_ = false;
      R_pm[k] = dprofit_droot_collar_psi(p);
    }
    PPFD_ += h;
    electron_transport_ = electron_transport();
    input_adjoints[i_ppfd] += mu * (R_pm[0] - R_pm[1]) / (2.0 * h);
  }
  {
    const double h = kappa * 1e-6;
    double R_pm[2];
    for (int k = 0; k < 2; ++k) {
      leaf_specific_conductance_max_ += (k == 0 ? h : -2.0 * h);
      transpiration_cached_ = false;
      R_pm[k] = dprofit_droot_collar_psi(p);
    }
    leaf_specific_conductance_max_ = kappa;
    input_adjoints[i_kappa] += mu * (R_pm[0] - R_pm[1]) / (2.0 * h);
  }

  poison_unwanted();
  restore_operating_point();
}

double Leaf::arrh_curve(double Ea, double ref_value, double leaf_temp) const {


  return ref_value*exp(Ea*((leaf_temp+C_to_K) - (25 + C_to_K))/((25 + C_to_K)*R*(leaf_temp+C_to_K)));
}

double Leaf::peak_arrh_curve(double Ea, double ref_value, double leaf_temp, double H_d, double d_S) const {
  double arrh = arrh_curve(Ea, ref_value, leaf_temp);
  double arg2 = 1 + exp((d_S*(25 + C_to_K) - H_d)/(R*(25 + C_to_K)));
  double arg3 = 1 + exp((d_S*(leaf_temp + C_to_K) - H_d)/(R*(leaf_temp + C_to_K)));

  return arrh * arg2/arg3;
}


// transpiration supply functions

// returns proportion of conductance taken from hydraulic vulnerability curve (unitless)
double Leaf::proportion_of_conductivity(double psi) const {

  return exp(-pow((psi / b), c));
}

// Build the knot grid {0, step, 2*step, .., <= psi_max} (psi_max = the potential
// magnitude at which conductivity drops to 1%, step = psi_max/resolution) and
// the cumulative vulnerability integral
//   G(m) = int_0^m exp(-(s/b)^c) ds = (b/c) * gamma_lower(1/c, (m/b)^c)
// (lower incomplete gamma) seeded from this closed form. Seeding knots with the
// closed form instead of a running trapezoid sum removes the dominant quadrature
// bias at no hot-path cost -- same knots, same tk::spline, same O(1) eval. See
// issue #468 and scripts/validate_gamma_transform.R.
//
// Shared by setup_transpiration (xylem) and setup_root_vulnerability (roots);
// each caller wires the resulting knots into its own interpolator(s).
void Leaf::build_cumulative_vulnerability_integral(double b, double c,
                                                   double resolution,
                                                   std::vector<double>& x,
                                                   std::vector<double>& y_integral) {
  x = std::vector<double>{0.0};
  y_integral = std::vector<double>{0.0}; // G(0) = 0
  double psi_max = b * pow(log(1.0 / 0.01), 1.0 / c);
  double step = psi_max / resolution;
  for (double psi = step; psi <= psi_max; psi += step) {
    x.push_back(psi);
    y_integral.push_back((b / c) *
                         boost::math::tgamma_lower(1.0 / c, pow(psi / b, c)));
  }
}

// pre-compute root vulnerability curve f(psi) = exp(-(|psi|/b_root)^c_root) as a spline,
// evaluated over the range [0, psi_max_root] where conductivity drops to 1%.
// This avoids repeated exp(pow(...)) calls inside E_from_Soil_to_Root_Collar.
void Leaf::setup_root_vulnerability(double resolution) {
  std::vector<double> x_psi_root, y_integral;
  build_cumulative_vulnerability_integral(root_b, root_c, resolution,
                                          x_psi_root, y_integral);

  // f_r conductivity knots on the same grid. f_r(0) = exp(-pow(0,root_c)) = 1.
  std::vector<double> y_f_r(x_psi_root.size());
  for (size_t i = 0; i < x_psi_root.size(); ++i) {
    y_f_r[i] = exp(-pow(x_psi_root[i] / root_b, root_c));
  }
  root_vuln_from_psi.init(x_psi_root, y_f_r);
  root_vuln_from_psi.set_extrapolate(true); // clamp to last value beyond range

  root_vuln_integral_from_psi.init(x_psi_root, y_integral);
  // linear extrapolation beyond range: slope ~= f_r at the tail (~1%), so the
  // integral keeps growing consistently with the clamped-conductivity tail.
  root_vuln_integral_from_psi.set_extrapolate(true);
}

// Same knot values setup_transpiration builds, at (b_at, c_at) rather than at
// the members, on the grid `x` the caller holds still (see the header).
void Leaf::set_transpiration_at(double b_at, double c_at,
                                const std::vector<double>& x) {
  std::vector<double> y(x.size(), 0.0);
  for (size_t i = 1; i < x.size(); ++i) {
    y[i] = (b_at / c_at) *
           boost::math::tgamma_lower(1.0 / c_at, pow(x[i] / b_at, c_at));
  }
  transpiration_from_psi.init(x, y);
  transpiration_from_psi.set_extrapolate(false);
  psi_from_transpiration.init(y, x);
  psi_from_transpiration.set_extrapolate(false);
  transpiration_cached_ = false;
}

// Same knot values setup_root_vulnerability builds, at (b_at, c_at) rather than
// at the members, on the grid `x` the caller holds still (see the header).
void Leaf::set_root_vulnerability_at(double b_at, double c_at,
                                     const std::vector<double>& x) {
  std::vector<double> y_integral(x.size(), 0.0), y_f_r(x.size());
  for (size_t i = 0; i < x.size(); ++i) {
    if (i > 0) {
      y_integral[i] =
          (b_at / c_at) *
          boost::math::tgamma_lower(1.0 / c_at, pow(x[i] / b_at, c_at));
    }
    y_f_r[i] = exp(-pow(x[i] / b_at, c_at));
  }
  root_vuln_from_psi.init(x, y_f_r);
  root_vuln_from_psi.set_extrapolate(true);
  root_vuln_integral_from_psi.init(x, y_integral);
  root_vuln_integral_from_psi.set_extrapolate(true);
}

// set spline for proportion of conductivity
void Leaf::setup_transpiration(double resolution) {
  std::vector<double> x_psi_, y_cumulative_transpiration_;
  build_cumulative_vulnerability_integral(b, c, resolution, x_psi_,
                                          y_cumulative_transpiration_);

  // setup interpolator
  transpiration_from_psi.init(x_psi_, y_cumulative_transpiration_);
  transpiration_from_psi.set_extrapolate(false);

  psi_from_transpiration.init(y_cumulative_transpiration_, x_psi_);
  psi_from_transpiration.set_extrapolate(false);
}

// replace f with some other function, returns E kg m^-2 s^-1

double Leaf::transpiration_full_integration(double psi_stem, double psi_upstream) {
  std::function<double(double)> f;
  f = [&](double psi) -> double { return proportion_of_conductivity(psi); };
  
  return leaf_specific_conductance_max_ * integrator.integrate(f, psi_upstream, psi_stem);
 }

//calculates supply-side transpiration from psi_stem and root_collar_psi_, returns kg h20 s^-1 m^-2 LA
// SIGN: psi_stem and psi_upstream are POSITIVE magnitudes here (passed straight
// to the spline). Contrast transpiration_to_psi_stem below. See the sign-
// conventions block above E_from_Soil_to_Root_Collar.
double Leaf::transpiration(double psi_stem, double psi_upstream) {

  // 1-entry memo: identical (psi_stem, psi_upstream) is requested several times
  // per profit evaluation; return the cached value (bit-identical) to skip the
  // redundant spline lookups. Cache invalidated in set_physiology.
  if (transpiration_cached_ &&
      psi_stem == transpiration_cache_psi_stem_ &&
      psi_upstream == transpiration_cache_psi_upstream_) {
    return transpiration_cache_value_;
  }

  // integration of proportion_of_conductivity over [root_collar_psi_, psi_stem]
  const double E = leaf_specific_conductance_max_ *
    (transpiration_from_psi.eval(psi_stem) - transpiration_from_psi.eval(psi_upstream));
  // return (transpiration_full_integration(psi_stem));

  transpiration_cache_psi_stem_ = psi_stem;
  transpiration_cache_psi_upstream_ = psi_upstream;
  transpiration_cache_value_ = E;
  transpiration_cached_ = true;
  return E;
}

// converts a known transpiration to its corresponding psi_stem, returns -MPa
// SIGN: unlike transpiration() above, psi_upstream here is a SIGNED (negative)
// potential, so it is flipped with a leading `-` before the spline lookup. The
// two functions are inverses called with opposite-sign psi_upstream.
double Leaf::transpiration_to_psi_stem(double transpiration_, double psi_upstream) {
  // integration of proportion_of_conductivity over [root_collar_psi_, psi_stem]


  double E_psi_stem = transpiration_/leaf_specific_conductance_max_ +  transpiration_from_psi.eval(-psi_upstream);


  return psi_from_transpiration.eval(E_psi_stem);
  }

// returns stomatal conductance to CO2, mol C m^-2 LA s^-1
double Leaf:: stom_cond_CO2(double psi_stem, double psi_upstream) {
  double transpiration_ = transpiration(psi_stem, psi_upstream);
  return atm_kpa_ * transpiration_ * kg_to_mol_h2o / atm_vpd_ / H2O_CO2_stom_diff_ratio;
}


// biochemical photosynthesis model equations
//ensure that units of PPFD_ actually correspond to something real.
// electron trnansport rate based on light availability and vcmax assuming co-limitation hypothesis
double Leaf::electron_transport() {



  double electron_transport_ = (a * PPFD_ + jmax_ - sqrt(pow(a * PPFD_ + jmax_, 2) - 
  4 * curv_fact_elec_trans * a * PPFD_ * jmax_)) / (2 * curv_fact_elec_trans); // check brackets are correct

  // double electron_transport_ = (4*a*PPFD_)/sqrt(pow(4*a*PPFD_/jmax_,2)+ 1);
    return electron_transport_;           
}

//calculate the rubisco-limited assimilation rate, returns umol m^-2 s^-1
double Leaf::assim_rubisco_limited(double ci_) {

  return (vcmax_ * (ci_ - gamma_ * umol_per_mol_to_Pa)) / (ci_ + km_);

}

//calculate the light-limited assimilation rate, returns umol m^-2 s^-1
double Leaf::assim_electron_limited(double ci_) {
  

  return electron_transport_ / 4 *
  ((ci_ - gamma_ * umol_per_mol_to_Pa) / (ci_ + 2 * gamma_ * umol_per_mol_to_Pa));
}

// returns co-limited assimilation umol m^-2 s^-1, NET of dark respiration
// (the trailing `- R_d_`), so gross assimilation is this value + R_d_.
double Leaf::assim_colimited(double ci_) {

  double assim_rubisco_limited_ = assim_rubisco_limited(ci_) ;
  double assim_electron_limited_ = assim_electron_limited(ci_);

  return (assim_rubisco_limited_ + assim_electron_limited_ - sqrt(pow(assim_rubisco_limited_ + assim_electron_limited_, 2) - 4 * curv_fact_colim * assim_rubisco_limited_ * assim_electron_limited_)) /
             (2 * curv_fact_colim)- R_d_;


}


// A - gc curves

// returns difference between co-limited assimilation and stom_cond_CO2, to be minimised (umol m^-2 s^-1)
double Leaf::assim_minus_stom_cond_CO2(double x, double psi_stem, double psi_upstream) {

  double assim_colimited_x_ = assim_colimited(x);

  double stom_cond_CO2_x_ = stom_cond_CO2(psi_stem, psi_upstream);
  return assim_colimited_x_ * umol_to_mol -
         (stom_cond_CO2_x_ * (ca_ - x) / (atm_kpa_ * kPa_to_Pa));
}

// converts psi stem to ci, used to find ci which makes A(ci) = gc(ca - ci)
double Leaf::psi_stem_to_ci(double psi_stem, double psi_upstream) {
  const double stom_cond_CO2_fixed = stom_cond_CO2(psi_stem, psi_upstream);

  // Propagate non-finite inputs as NA rather than entering the solver. A
  // non-finite psi_stem (e.g. NA from profit_psi_stem_TF(NA, .)) makes gc and
  // hence the whole target non-finite. The previous bisection returned NaN
  // silently in this case; the bracketing TOMS748 solver below instead throws
  // ("a and b do not bracket the root"), so guard explicitly to preserve the
  // NA-in -> NA-out contract (see test-leaf.r "Basic functions").
  if (!std::isfinite(stom_cond_CO2_fixed)) {
    return ci_ = NA_REAL;
  }

  auto target = [&](double x) mutable -> double {
    const double assim_colimited_x_ = assim_colimited(x);
    return assim_colimited_x_ * umol_to_mol -
      (stom_cond_CO2_fixed * (ca_ - x) / (atm_kpa_ * kPa_to_Pa));
  };

  // #486: this target (assim_colimited demand minus the linear gc supply) is
  // smooth and strictly monotone over (gamma*, ca] -- no singularity at the
  // bracket ends (Ar,Ae vanish linearly at gamma* so sqrt(disc) is linear, not
  // singular) and its only sharp feature is the colimitation elbow far below the
  // operating root. So it is a well-behaved case for a superlinear bracketing
  // solver: TOMS748 reaches the same root in ~9 evals vs bisection's ~29 at the
  // same 1e-7 tol. This is deliberately scoped to psi_stem_to_ci ONLY; the
  // hydraulic find_root_psi path keeps bisection (its target is not smooth -- see
  // the warning on util::uniroot_smooth).
  try {
    return ci_ = util::uniroot_smooth(target, gamma_ * umol_per_mol_to_Pa, ca_, 1e-7, ci_niter);
  } catch (const std::exception& e) {
    util::stop("psi_stem_to_ci failed: " + std::string(e.what()) +
               "; min=" + util::to_string(gamma_ * umol_per_mol_to_Pa) +
               "; max=" + util::to_string(ca_) +
               "; psi_stem=" + util::to_string(psi_stem) +
               "; psi_upstream=" + util::to_string(psi_upstream));
  }
}

// given psi_stem, find assimilation, transpiration and stomal conductance to c02
void Leaf::set_leaf_states_rates_from_psi_stem(double psi_stem, double psi_upstream) {

  if (psi_upstream >= psi_stem){
    ci_ = gamma_*umol_per_mol_to_Pa;
    transpiration_ = 0;
    stom_cond_CO2_ = 0;
    } else{
      if(assim_max_ < 0){
        ci_ = gamma_*umol_per_mol_to_Pa;
        transpiration_ = 0;
        stom_cond_CO2_ = 0;
        } else{
      ci_ = psi_stem_to_ci(psi_stem, psi_upstream);
      transpiration_ = transpiration(psi_stem, psi_upstream);
      stom_cond_CO2_ = atm_kpa_ * transpiration_ * kg_to_mol_h2o / atm_vpd_ / H2O_CO2_stom_diff_ratio;
      }
    }
  assim_colimited_ = assim_colimited(ci_);
}


// Hydraulic cost equations

// Sperry et al. 2017; Sabot et al. 2020 implementation

double Leaf::hydraulic_cost_Sperry(double psi_stem, double psi_upstream) {
  // Cost is definitionally zero when the potentials are equal. Returning it
  // explicitly avoids a tiny non-zero residual from FMA contraction of the
  // k_l_soil_ - k_l_stem_ subtraction (arch-dependent; see arm64 build, #468).
  if (psi_stem == psi_upstream) {
    hydraulic_cost_ = 0.0;
    return hydraulic_cost_;
  }
  double k_l_soil_ = leaf_specific_conductance_max_ * proportion_of_conductivity(psi_upstream);
  double k_l_stem_ = leaf_specific_conductance_max_ * proportion_of_conductivity(psi_stem);
  
  hydraulic_cost_ = k_l_soil_ - k_l_stem_;
  
  return hydraulic_cost_;
}

double Leaf::hydraulic_cost_TF(double psi_stem) {

  hydraulic_cost_ = g1_TF24 * pow((1 - proportion_of_conductivity(psi_stem)), beta2);

return hydraulic_cost_;
}

// Profit functions

double Leaf::profit_psi_stem_Sperry(double psi_stem, double psi_upstream) {

set_leaf_states_rates_from_psi_stem(psi_stem, psi_upstream);

  double benefit_ = assim_colimited_;
  double cost = hydraulic_cost_Sperry(psi_stem, psi_upstream);

  return benefit_ - lambda_ * cost;
}


double Leaf::profit_psi_stem_TF(double psi_stem, double psi_upstream) {
set_leaf_states_rates_from_psi_stem(psi_stem, psi_upstream);

double benefit_ = assim_colimited_;
  double cost = hydraulic_cost_TF(psi_stem);

  return benefit_ - cost;
}


//optimisation functions


// need docs on Golden Section Search.
void Leaf::optimise_psi_stem_Sperry() {

    if (!(psi_soil_.size() == 1)) {
    util::stop("psi soil must have only one value to use non-root-based profit optimisation methods");
  }

  opt_psi_stem_ = psi_soil_[0];


  if ((PPFD_ < 1.5e-8 )| (psi_soil_[0] > psi_crit)){
    profit_ = 0;
    transpiration_ = 0;
    stom_cond_CO2_ = 0;
    return;
  }

  // Maximise carbon profit over [psi_soil, psi_crit]. Brent's method (golden-
  // section + parabolic interpolation) converges super-linearly on this smooth
  // objective; we minimise -profit and recover the maximum from neg_profit_opt.
    double neg_profit_opt = 0.0;
    opt_psi_stem_ = util::brent_fmin(
        [&](double psi_stem) { return -profit_psi_stem_Sperry(psi_stem, psi_soil_[0]); },
        psi_soil_[0], psi_crit, GSS_tol_abs, &neg_profit_opt);
    profit_ = -neg_profit_opt;

  }
  

void Leaf::optimise_psi_stem_TF() {

  if (!(psi_soil_.size() == 1)) {
    util::stop("psi soil must have only one value to use non-root-based profit optimisation methods");
  }

  opt_psi_stem_ = psi_soil_[0];

  if (psi_soil_[0] > psi_crit){
    profit_ = profit_psi_stem_TF(psi_soil_[0], psi_soil_[0]);
    return;
  }

  // Maximise carbon profit over [psi_soil, psi_crit] via Brent's method
  // (minimise -profit), matching find_root_collar_psi's multi-layer solver.
    double neg_profit_opt = 0.0;
    opt_psi_stem_ = util::brent_fmin(
        [&](double psi_stem) { return -profit_psi_stem_TF(psi_stem, psi_soil_[0]); },
        psi_soil_[0], psi_crit, GSS_tol_abs, &neg_profit_opt);
    profit_ = -neg_profit_opt;

    return;
  }

// ===========================================================================
// MEDLYN STOMATAL-CONDUCTANCE MODEL (from develop #450)
// ---------------------------------------------------------------------------
// Standalone, R-callable coupling of the Medlyn (2011) optimal stomatal model
// to colimited photosynthesis. NOT invoked by the TF24 compute path, which
// optimises psi_stem directly via find_root_collar_psi; provided so the model
// remains available for R-level experimentation/comparison. beta_ is a soil-
// moisture stress factor in [0,1] from theta_/theta_w_/theta_fc_ (set in
// set_physiology from the default soil-moisture values).
// ===========================================================================
double Leaf::medlyn_model_gs(double assim_colimited_){

  double beta_ = (theta_ - theta_w_)/(theta_fc_ - theta_w_);

  if(atm_vpd == 0){
     medlyn_model_gs_ = g0;
  } else{
     medlyn_model_gs_ = g0 + 1.6*(1 + (g1*beta_)/sqrt(atm_vpd_))*(assim_colimited_/(ca_*(1/umol_per_mol_to_Pa)));
  }
  return medlyn_model_gs_;
}

// Supply==demand residual for the Medlyn solver, as a function of ci (x, Pa).
// It is the difference between the Medlyn optimal stomatal conductance and the
// diffusion-implied conductance, *multiplied through by (ca_ - x)* so it stays
// finite across the whole [gamma*, ca_] bracket -- the raw gs difference has a
// 1/(ca_-x) singularity at the upper end. The root (zero crossing) is identical
// to that of the raw difference for x < ca_:
//   gs_medlyn*(ca_-x) - gs_coupled*(ca_-x),  where
//   gs_coupled*(ca_-x) = assim * (atm_kpa_*kPa_to_Pa) * 1.6 / 1e6.
double Leaf::medlyn_stom_cond_minus_coupled_stom_cond(double x) {
  const double assim_colimited_x_ = assim_colimited(x);
  medlyn_model_gs_ = medlyn_model_gs(assim_colimited_x_);
  return medlyn_model_gs_ * (ca_ - x)
         - assim_colimited_x_ * (atm_kpa_ * kPa_to_Pa) * 1.6 / 1e6;
}

// Solve for the leaf-internal CO2 (ci) at which the Medlyn optimal stomatal
// conductance balances the diffusion-implied conductance, with Brent's method
// (util::uniroot). This replaces an earlier golden-section search on
// 1/|gs difference|, which could lock onto a spurious interior maximum
// (observed at high VPD).
//
// The residual is not monotone on [gamma*, ca_]: it has a hump and is negative
// at BOTH ends (near gamma* assimilation is below the compensation point;
// near ca_ the diffusion conductance diverges), so the interval can contain a
// spurious sub-compensation root as well as the meaningful Medlyn root. We
// therefore first locate the residual's maximum (Brent minimiser on -residual),
// then root-find on [argmax, ca_], which isolates the physically meaningful
// high-ci operating point in every case (including g0 == 0).
void Leaf::solve_medlyn_ci_numerical(){
  auto target = [&](double x) -> double {
    return medlyn_stom_cond_minus_coupled_stom_cond(x);
  };
  const double lo = gamma_ * umol_per_mol_to_Pa;
  const double hi = ca_;

  double neg_peak = 0.0;
  const double ci_peak =
      util::brent_fmin([&](double x) { return -target(x); }, lo, hi, ci_abs_tol,
                       &neg_peak);
  const double residual_peak = -neg_peak;

  if (residual_peak <= 0.0) {
    // Supply never reaches demand: no feasible Medlyn operating point. Report
    // the closest approach (the residual maximum) rather than failing.
    ci_ = ci_peak;
  } else {
    try {
      ci_ = util::uniroot(target, ci_peak, hi, ci_abs_tol, ci_niter);
    } catch (const std::exception& e) {
      util::stop("solve_medlyn_ci_numerical failed: " + std::string(e.what()) +
                 "; ci_peak=" + util::to_string(ci_peak) +
                 "; max=" + util::to_string(hi));
    }
  }
  assim_colimited_ = assim_colimited(ci_);
  stom_cond_CO2_ = medlyn_model_gs(assim_colimited_);
  return;
}

void Leaf::solve_medlyn_ci_analytical(){

  ci_ = ca_ * (g1/(g1 + sqrt(atm_vpd_)));
  assim_colimited_ = assim_colimited(ci_);
  stom_cond_CO2_ = medlyn_model_gs(assim_colimited_);
  return;
}

} // namespace plant
