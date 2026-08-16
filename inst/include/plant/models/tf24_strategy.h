// Built from  inst/include/plant/models/ff16_strategy.h on Mon Feb 12 09:52:27 2024 using the scaffolder, from the strategy:  FF16
// -*-c++-*-
#ifndef PLANT_PLANT_TF24_STRATEGY_H_
#define PLANT_PLANT_TF24_STRATEGY_H_

#include <plant/strategy.h>
#include <plant/models/tf24_environment.h>
#include <plant/qag.h>
#include <plant/leaf_model.h>
#include <plant/canopy_shape.h>
#include <odelia/ode_util.hpp>
#include <odelia/implicit_node.hpp>
#include <type_traits>

namespace plant {

// Biological (user-settable) parameters for the TF24 strategy. Held as a value
// member `pars` on TF24_Strategy and exposed to R as a nested RcppR6 list class
// (access as `s$pars$lma`). Only parameters that were previously exposed to R
// live here; derived quantities (eta_c, height_0, ...), the embedded Leaf
// model, solver tolerances and hard-coded hydraulic-root constants stay as
// plain members on the strategy.
template <typename S = double>
struct TF24_Pars {
  using value_type = S;

  // A default member initialiser has no block to declare the library pow in, so
  // the derived defaults below raise their base through here instead.
  static S power(const S& base, const S& exponent) {
    using std::pow;
    return pow(base, exponent);
  }

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
  // a_dG1 / a_dG2 now shape the *storage-dependent* growth mortality (#517):
  //   mortality_storage_dependent_dt(r) = a_dG1 * exp(-a_dG2 * r),  r = S/S_max.
  // Previously these fed the instantaneous productivity mortality
  // a_dG1*exp(-a_dG2*productivity_area), which overflowed to ~1e32 under deep
  // carbon deficit (the #550 blow-up). Buffering through the NSC pool bounds the
  // argument to r in [0,1], so this term is now bounded in [a_dG1*e^-a_dG2, a_dG1].
  S a_dG1  = 5.5;            // Max growth-related (low-reserve) mortality [/yr]
  S a_dG2  = 20.0;           // Sensitivity of mortality to relative reserves
  // * NSC storage pool (#517) -- buffers growth & mortality against short-term
  //   productivity swings. See Stefaniak et al. 2026 (plantNSC) for the design.
  //   Carbon first charges storage; growth/reproduction are mobilised from it,
  //   but gated so growth only proceeds when reserves are ample (a plant should
  //   not grow down its stores). Reserves-based mortality rises as they deplete.
  S a_st1  = 0.10;           // Storage capacity per unit sapwood mass [kg NSC / kg]
  S a_st2  = 0.1;            // Reserve fraction at which growth is half-on [0-1]
  S a_st3  = 0.8;            // Initial storage at birth [fraction of capacity]
  // * Light capture
  S k_I = 0.5;
  // * Leaf hydraulic / photosynthesis traits (default Eucalyptus saligna)
  S vcmax_25 = 96;
  S p_50 = 1.85;
  S K_s = 1;
  S c = log(log(1-0.5)/log(1-0.88))/(log(p_50) - log(5.16));
  S b = p_50 / power(-log(1 - 50.0 / 100.0), 1 / c);
  S psi_crit = b*power(log(1/0.05),1/c); // derived from b and c
  S beta1 = 20000;
  S beta2 = 1.5;
  S g1_TF24 = 7.5;
  S jmax_25 = vcmax_25*1.64;
  S a = 0.30; // effective quantum yield of electron transport
  S curv_fact_elec_trans = 0.7;
  S curv_fact_colim = 0.99;
  // Dark respiration at 25 C, the term net assimilation subtracts. The Leaf
  // constructor does not take it, so until this member existed the leaf ran at
  // its own default and no parameter could move it; the value here is that
  // default, which is what keeps the forward model unchanged.
  S R_d_25 = 1.44;
  S var_sapwood_volume_cost = 1;
  // nitrogen allocation traits (parameterised from Austraits 4.1.0)
  S nmass_l = 13e-3; // kg N kg^-1 mass
  S nmass_s = 1.98e-3; // kg N kg^-1 mass
  S nmass_b = 3.40e-3; // kg N kg^-1 mass
  S nmass_r = 3.35e-3; // kg N kg^-1 mass
  S dmass_dN = 0; // change in mass per change in kg kg^-1 N
  // shape exponent for the Q() root-fraction-with-depth profile
  S root_depth_shape_eta = 0.2;
  // * Root hydraulics
  // Root vulnerability curve, proportion of conductivity =
  // exp(-(psi/root_b)^root_c). Declared before root_psi_crit, which derives
  // from them. Previously fixed members of TF24_Strategy and unreachable from
  // R, which pinned root shutoff at ~5.87 MPa -- too conservative for taxa that
  // operate below that (e.g. Acacia aneura), and unavailable for calibration.
  S root_c = 2.680147;
  S root_b = 3.898245;
  // Potential at 5% remaining root conductivity [MPa]. Derived, exactly as
  // psi_crit is from b and c: if you set root_b or root_c directly, set this
  // too, or the vulnerability curve and the shutoff threshold disagree.
  S root_psi_crit = root_b*power(log(1/0.05),1/root_c);
  // Maximum rooting depth [m]. Rooting depth is min(height, rooting_depth_max),
  // so this also bounds the depth over which roots can draw water. Should not
  // exceed the soil column depth (TF24_Environment `depth`, default 1.5 m):
  // layers below the column do not exist, so deepening roots alone gains
  // nothing without deepening the soil as well.
  S rooting_depth_max = 1.5;
  // Germination
  S recruitment_decay = 0.0;
  // Penman-Monteith leaf energy balance (#523). use_energy_balance gates PM
  // (0 = off, today's Tleaf=Tair behaviour; != 0 = on); default off preserves
  // backward compatibility. d is the characteristic leaf dimension (m) for the
  // aerodynamic resistance ra = C_ra*sqrt(d/U0); inert while PM is off. Both stay
  // double: the gate is compared rather than differentiated, and d has no row in
  // the leaf's supplied Jacobian, so declaring it active would read an exact zero
  // that is indistinguishable from insensitivity.
  //
  // Both carry S so that field_ptrs() lists them and a rebind carries them: the
  // set field_ptrs() spans is the whole parameter set, not the differentiable
  // subset, so listing a parameter here does not make it a gradient target.
  // ad_parameter_names() decides that, and neither is in it -- d has no row in
  // the leaf's supplied Jacobian, and an unrouted target would read an exact
  // zero indistinguishable from insensitivity.
  S use_energy_balance = 0.0;
  S d = 0.05;

  // Every member above, in declaration order. Carries the whole parameter set
  // across a scalar change; ad_parameters() is the differentiable subset only.
  std::vector<S*> field_ptrs() {
    return {
      &lma, &rho, &hmat, &omega, &eta, &theta, &a_l1, &a_l2, &a_r1, &a_b1,
      &r_s, &r_b, &r_r, &r_l, &a_y, &a_bio, &k_l, &k_b, &k_s, &k_r,
      &a_p1, &a_p2, &a_f3, &a_f1, &a_f2, &S_D, &a_d0, &d_I, &a_dG1, &a_dG2,
      &a_st1, &a_st2, &a_st3, &k_I, &vcmax_25, &p_50, &K_s, &c, &b, &psi_crit,
      &beta1, &beta2, &g1_TF24, &jmax_25, &a, &curv_fact_elec_trans,
      &curv_fact_colim, &R_d_25,
      &var_sapwood_volume_cost, &nmass_l, &nmass_s, &nmass_b,
      &nmass_r, &dmass_dN, &root_depth_shape_eta, &root_c, &root_b,
      &root_psi_crit, &rooting_depth_max, &recruitment_decay,
      &use_energy_balance, &d
    };
  }
  static constexpr size_t field_count = 62;
};

// Every member of TF24_Pars is an S, so one added without extending
// field_ptrs() changes this size and is refused rather than dropped.
static_assert(sizeof(TF24_Pars<double>) ==
              TF24_Pars<double>::field_count * sizeof(double),
              "TF24_Pars has a member field_ptrs() does not list");

// Templated on the scalar S the state, the traits and everything derived from
// them carry; double is production. The embedded Leaf, the Control tolerances
// and the extrinsic drivers stay double.
template <typename S = double>
class TF24_Strategy: public Strategy<TF24_Environment<S>> {
public:
  using value_type = S;

  typedef std::shared_ptr<TF24_Strategy> ptr;
  TF24_Strategy();

  // Scientific version. Bump ONLY when equations or default parameters change
  // the simulation output for identical inputs. Do NOT bump for refactors,
  // performance, interface, or serialisation changes. Bumping invalidates
  // logpile's cache for this model (see plant::model_version() / model_id()).
  // Starts at 2: a published result exists using pre-versioning "v1" science.
  // v3 (#517/#554): NSC storage pool + reserve-gated growth and reserves-based
  // mortality change the simulation output for identical inputs; TF24f's
  // compound version auto-tracks this to 3.1.
  // v4: the hydraulic shut-down exits no longer leave the previous step's
  // transport state in place, so a shut-down plant stops drawing water from the
  // soil (it used to keep extracting its last wet-step uptake, because `Leaf`
  // is reused across steps and soil_consumption_ feeds the patch water
  // balance). Water-limited runs therefore change: on the scenario gateway,
  // offspring production moves by up to 5e-3 relative on 5 of 8 scenarios,
  // while every success/failure classification is unchanged. TF24f's compound
  // version auto-tracks this to 4.1.
  // v5: the leaf gas-exchange and hydraulics model is now the standalone
  // `phylloptim` package rather than a copy in this repo, and the swap carries four science
  // changes. Measured on the one-species SCM scenario of test-strategy-tf24.R
  // (max_patch_lifetime = 5), offspring production moves 81.9083 -> 83.9026,
  // i.e. **+2.4%**. Attributed by re-running both arms with the atm_kpa driver
  // forced to 101.3, which removes the pressure change and leaves the rest:
  //
  //     arm                  atm_kpa 100.5     atm_kpa 101.3
  //     this repo's leaf        81.9083           81.8201
  //     the leaf package        83.9026           81.8985
  //
  // So **the pressure fix is ~25x the rest of the swap put together** (+2.4%
  // against +0.10%), which was not the expectation going in. The leaf package's
  // ppm-to-Pa conversion is derived from atm_kpa (phylloptim #15 item 10c) instead
  // of hard-coded at 0.1013 = 1e-6 * 101300 Pa; TF24_Environment's atm_kpa driver
  // defaults to **100.5**, so Gamma*, Kc, Ko, Km and the ci root-find bounds all
  // move. Before the fix the conductance side of the model responded to atm_kpa
  // while the photosynthesis side silently assumed sea level -- visible in the
  // table as this repo's leaf moving only -0.11% across the same 0.8 kPa that
  // moves the package -2.4%.
  //
  // The remaining +0.10% is the two further stale-state exits (phylloptim #26,
  // ported from #585) plus the supply-path extraction (phylloptim #2). TF24f's
  // compound version auto-tracks this to 5.1.
  // v6: the leaf package moved to ONE representation for water potential --
  // positive magnitudes throughout (phylloptim #25). Two consequences, and the first
  // is why this is a version bump rather than a refactor:
  //   * the **`opt_root_psi` aux changes sign**. It is now the positive magnitude,
  //     which is what TF24f's `opt_root_psi_state` has always held. Before, the aux
  //     reported the signed potential while tf24f_strategy.cpp negated it back for
  //     the state -- an inconsistency in plant's own reported outputs, and the two
  //     compensating negations are deleted here. Any stored output or cached
  //     analysis reading that aux would silently change meaning, which is exactly
  //     what model_version() exists to catch.
  //   * outputs move slightly: one-species SCM offspring production 83.9026 ->
  //     83.8761, i.e. **-3.2e-4 relative**. This is NOT an equation change. The
  //     rewrite is exactly sign-symmetric in IEEE; what is not is boost's TOMS748,
  //     whose iterates depend on the bracket's orientation, and #25 reverses it.
  //     Measured there: 12 of 288 golden operating points differ by 1-3 ULP, the
  //     rest exactly. The SCM's adaptive stepper and node schedule amplify that to
  //     3e-4 -- within the ~GSS_tol_abs (1e-3) ceiling the leaf package documents,
  //     and small enough that every pinned test value and the exact stochastic
  //     TF24 counts (101/23) pass unchanged.
  // TF24f's compound version auto-tracks this to 6.1.
  // v7: the collar bracket is finally clamped to root_psi_crit, the potential at
  // which root conductivity is down to 5% (phylloptim #24, plant #584). The clamp was
  // written as a std::max against a *signed* root_psi_crit, so it could never bind
  // and the solver optimised over a collar the root system cannot supply. **The
  // window is 1.2 MPa wide at TF24's defaults** -- psi_crit = 7.085493 against
  // root_psi_crit = 5.870283 -- so this is a dry-corner correction, not a rounding
  // one. Two regimes inside it: the interval is tightened (the plant still
  // transpires, at a wetter collar), or root_psi_crit lands below the zero-uptake
  // collar and the plant shuts down because no operating point both moves water and
  // stays inside the root limit.
  //
  // No tested scenario moves: the standard SCM run stays at 83.8761 offspring, bit
  // for bit, because a mesic patch never drives the collar past 5.87 MPa. It is
  // still a version bump -- a user running a dry scenario gets different (correct)
  // numbers for identical inputs, and logpile's cache has to know.
  // TF24f's compound version auto-tracks this to 7.1.
  // v8: `TF24_Environment`'s `atm_kpa` driver default goes **100.5 -> 101.3**, and
  // this is the entry to read if you only read one. It largely CANCELS v5.
  //
  // v5 recorded +2.4% from deriving the leaf's ppm -> Pa conversion from `atm_kpa`
  // instead of hard-coding 0.1013. That constant *was* 101.3 kPa in disguise
  // (1e-6 * 101300 Pa), so the shift was not the fix doing damage -- it was this
  // driver disagreeing with the rest of the model. 100.5 arrived in `34d46ac2`
  // ("Simplify scm & environment interface", #446), an interface refactor that does
  // not mention atmospheric pressure, with no rationale recorded anywhere, while
  // every leaf-level test used 101.3. An artefact, not a site elevation.
  //
  // Pinning it to the value the model already assumed collapses the whole branch's
  // movement. Net effect of ALL of it (the swap, the #15 catch-up, the #26 ported
  // fixes, #25 and #24) against `develop`:
  //
  //     one-species SCM offspring   81.9083 -> 81.7426     -0.20%
  //     stochastic TF24 counts      103 / 28 -> 103 / 28   unchanged, exactly
  //
  // So **+2.43% became -0.20%**, and every pinned baseline reverts to develop's own
  // values: the two SCM offspring figures pass at their original 82.09077702 /
  // 67.54060383, and the seeded stochastic integers match bit for bit. That exact
  // match on discrete counts is a sharper statement than any tolerance-based check
  // that the swap preserves TF24's science.
  //
  // The fix itself is NOT undone -- an off-sea-level run still gets a self-consistent
  // Gamma*/Kc/Ko/Km and conductance side, which is the whole point of item 10c. Set
  // `atm_kpa` per site if you mean altitude; it just no longer defaults to an
  // altitude nobody chose.
  // v9: as FF16 v2 -- the light field and the census take their trapezium widths
  // from the coordinate the density is carried in. TF24 defaults to the
  // birth-date coordinate, so this moves its output; the water reduction already
  // integrated over birth dates and is unchanged.
  static constexpr int scientific_version = 9;

  S compute_average_light_environment(S z, S height,
                                      const TF24_Environment<S> &environment);

  // calculate the amount of water transpired relativised by leaf area index.

  S evapotranspiration_dt(S area_leaf_, int soil_layer);


  // Overrides ----------------------------------------------

  // update this when the length of state_names changes
  static size_t state_size () { return 6; }
  // update this when the length of aux_names changes
  size_t aux_size () { return aux_names().size(); }

  static std::vector<std::string> state_names() {
    return  std::vector<std::string>({
      "height",
      "mortality",
      "fecundity",
      "area_heartwood",
      "mass_heartwood",
      "storage"
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
      // Net CO2 assimilation at the optimal operating point, per unit leaf
      // area (umol CO2 m^-2 s^-1). Net, not gross: Leaf::assim_colimited()
      // subtracts dark respiration R_d_, so gross = assimilation + R_d_ with
      // R_d_ = 0.015 * vcmax_ at the acclimated vcmax_.
      "assimilation"
    });
    // add the associated computation to compute_rates and compute there
    if (this->collect_all_auxiliary) {
      ret.push_back("area_sapwood");
    }
    return ret;
  }

  // Addresses of the parameters a gradient can be taken with respect to, in the
  // order ad_parameter_names() gives. Both allocate, so take them once per
  // gradient evaluation and hold them for the run rather than per block; the
  // strategy is shared and the fields do not move. Index against .size().
  std::vector<S*> ad_parameters() {
    return {
      &pars.lma, &pars.rho, &pars.hmat, &pars.omega,
      &pars.theta, &pars.a_l1, &pars.a_l2, &pars.a_r1, &pars.a_b1,
      &pars.r_s, &pars.r_b, &pars.r_r, &pars.r_l, &pars.a_y, &pars.a_bio,
      &pars.k_l, &pars.k_b, &pars.k_s, &pars.k_r,
      &pars.a_f3, &pars.a_f1, &pars.a_f2,
      &pars.a_d0, &pars.d_I, &pars.a_dG1, &pars.a_dG2,
      &pars.a_st1, &pars.a_st2, &pars.a_st3,
      &pars.k_I,
      &pars.K_s, &pars.c, &pars.b, &pars.psi_crit,
      &pars.beta2, &pars.g1_TF24, &pars.a,
      &pars.curv_fact_elec_trans, &pars.curv_fact_colim,
      &pars.vcmax_25, &pars.jmax_25, &pars.R_d_25,
      &pars.root_c, &pars.root_b, &pars.root_psi_crit,
      &pars.rooting_depth_max, &pars.recruitment_decay
    };
  }

  // The TF24_Pars members ad_parameters() addresses, in the same order. Absent
  // from both: eta and root_depth_shape_eta, whose exponents reach a base of 0
  // in CanopyShape::Qp and the soil retention curves, where the recorded
  // derivative u^k * log(u) is a NaN; and
  // eleven that no equation on this path reads. Two shapes, both giving a gradient
  // row that is exactly zero for a reason no measurement reveals. a_p1 and a_p2
  // belong to the light-response curve the Farquhar leaf replaced; beta1, S_D,
  // var_sapwood_volume_cost, nmass_l, nmass_s, nmass_b, nmass_r and dmass_dN are
  // declared, carried and read by nothing here. p_50 is different and worse: the
  // only thing that reads it is TF24_Pars' own default initialisers for c and b, so
  // it is read once at construction and a value set afterwards reaches nothing.
  // A derivative with respect to it belongs to whatever computes c and b, which for
  // a run driven from traits is the R hyperparameter function.
  std::vector<std::string> ad_parameter_names() {
    return {
      "lma", "rho", "hmat", "omega",
      "theta", "a_l1", "a_l2", "a_r1", "a_b1",
      "r_s", "r_b", "r_r", "r_l", "a_y", "a_bio",
      "k_l", "k_b", "k_s", "k_r",
      "a_f3", "a_f1", "a_f2",
      "a_d0", "d_I", "a_dG1", "a_dG2",
      "a_st1", "a_st2", "a_st3",
      "k_I",
      "K_s", "c", "b", "psi_crit",
      "beta2", "g1_TF24", "a",
      "curv_fact_elec_trans", "curv_fact_colim",
      "vcmax_25", "jmax_25", "R_d_25",
      "root_c", "root_b", "root_psi_crit",
      "rooting_depth_max", "recruitment_decay"
    };
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

  void compute_rates(const TF24_Environment<S>& environment,
                Internals<S>& vars);
  
  void compute_roots(const TF24_Environment<S>& environment,
                Internals<S>& vars);

  void update_dependent_aux(const int index, Internals<S>& vars);

  // * Mass production
  // [eqn 12] Gross annual CO2 assimilation
  S assimilation(const TF24_Environment<S>& environment, S height,
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

  virtual S net_mass_production_dt(const TF24_Environment<S>& environment,
                                S height, S area_leaf_,
                                S height_inverse);

  // Resolve the leaf operating point on the already-set-up `leaf` (i.e. after
  // leaf.set_physiology(...)). Base TF24 optimises the root-collar psi via
  // golden-section search; the TF24f variant overrides this to make the optimum
  // chase a tracked ODE state (#525). Called per crown light point from
  // net_mass_production_dt, so it must be virtual to dispatch to the override
  // when net_mass_production_dt is reused unchanged by the subclass.
  virtual void solve_leaf();

  // Read how the leaf's outputs respond to what it was given, and record each
  // output that re-enters the active chain carrying that response.
  //
  // `drive` re-supplies the leaf's physiology at a radiation and a soil profile
  // handed to it, leaving the operating point unsolved. It is the caller's,
  // because the caller is what owns those inputs; the leaf never sees an active
  // value and this function never solves.
  template <typename Drive, typename Rebuild>
  void record_leaf_outputs(const S& radiation, const std::vector<S>& psi_soil,
                           const S& conductance_max,
                           Drive drive, Rebuild rebuild_roots);
  // Strategy-agnostic entry point used by Individual<TF24> (#266): reads the
  // height state and the cached aux slots itself, so the generic Individual
  // does not need to know TF24's state/aux layout.
  S net_mass_production_dt(const TF24_Environment<S>& environment,
                                const Internals<S>& vars) {
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


  S mortality_dt(S relative_reserves, S cumulative_mortality) const;
  S mortality_growth_independent_dt()const ;
  // Storage-dependent growth mortality (#517): rises smoothly as relative
  // reserves r = S/S_max deplete, bounded in [a_dG1*e^-a_dG2, a_dG1].
  S mortality_storage_dependent_dt(S relative_reserves) const;
  // NSC storage capacity S_max = a_st1 * mass_sapwood [kg NSC].
  S storage_capacity(S area_leaf, S height) const;
  // Seed the storage state for a newly germinated individual (#517).
  void set_initial_states(const TF24_Environment<S>& environment, Internals<S>& vars);
  // [eqn 20] Survival of seedlings during establishment, from the carbon a
  // seedling produces at birth size. This form works that carbon out.
  S establishment_probability(const TF24_Environment<S>& environment);
  // The same, for a newborn whose rates have just been computed. A newborn is
  // already at birth size, so compute_rates has left that carbon in aux and the
  // leaf need not be solved there twice.
  S establishment_probability(const TF24_Environment<S>& environment,
                              const Internals<S>& vars) {
    return establishment_probability(environment,
                                     vars.aux(aux_idx_net_mass_production_dt));
  }
  // The equation the two above share.
  S establishment_probability(const TF24_Environment<S>& environment,
                              S net_mass_production_dt_);

  // * Competitive environment
  // [eqn 11] total projected leaf area above height above height `z` for given plant
  S compute_competition(S z, S height) const;
  // Optimised overload called from Individual<TF24>::compute_competition with the
  // cached competition_effect (= area_leaf(height)) and height_inverse (= 1/height)
  // aux values, matching the shared individual.h interface (no recompute per call).
  S compute_competition(S z, S area_leaf_,
                        S height_inverse) const;
  // Strategy-agnostic entry point used by Individual<TF24> (#266): reads the
  // cached competition_effect and height_inverse aux slots itself.
  S compute_competition(S z, const Internals<S>& vars) const {
    return compute_competition(z, vars.aux(aux_idx_competition_effect),
                               vars.aux(aux_idx_height_inverse));
  }

  // The competition contribution and its vertical derivative from one pass, so
  // u^eta is evaluated once. The first entry is bit-for-bit the one
  // compute_competition() returns.
  std::pair<S, S> compute_competition_and_slope(S z, const Internals<S>& vars) const {
    const S area_leaf_ = vars.aux(aux_idx_competition_effect);
    const S height_inverse = vars.aux(aux_idx_height_inverse);
    const S scale = pars.k_I * area_leaf_;
    const std::pair<S, S> Qq =
      canopy_shape.Q_and_q(z * height_inverse, z, height_inverse);
    return {scale * Qq.first, -(scale * Qq.second)};
  }


  // The fraction of root mass below soil depth `z`, for a plant rooted to
  // `rooting_depth` with shape exponent `eta_x` (pars.root_depth_shape_eta). The
  // canopy's own cumulative form is CanopyShape::Q, at pars.eta.
  S Q(S z, S rooting_depth, S eta_x) const;

  // The inverse of dheight_darea_leaf, so the allometry has one source.
  S darea_leaf_dheight(S area_leaf) const {
    return 1.0 / dheight_darea_leaf(area_leaf);
  }

  // The aim is to find a plant height that gives the correct seed mass.
  double height_seed(void) const;

  // The seed's height and leaf area at the current scalar.
  //
  // Preparation solves the height in plain arithmetic and cannot run at an active
  // scalar, so on a differentiated path the height is declared by the residual that
  // defines it -- live mass equals seed mass -- and the leaf area is derived from
  // the height it returns. The two come from one call because the leaf area's own
  // partials in the allometric constants and its chain through the height are the
  // same channel: taking either against the other held fixed mixes them.
  struct SeedGeometry { S height; S area_leaf; };
  SeedGeometry seed_geometry() const {
    if constexpr (std::is_same_v<S, double>) {
      return {height_0, area_leaf_0};
    } else {
      const S h = odelia::implicit_value<S>(
        height_0,
        [&](S y) -> S { return mass_live_given_height(y) - pars.omega; });
      return {h, area_leaf(h)};
    }
  }

  // Set constants within TF24_Strategy
  void prepare_strategy();

  // The same strategy at scalar U.
  template <class U> using rebind = TF24_Strategy<U>;

  // This strategy, already prepared, copied onto scalar U. prepare_strategy()
  // is refused at an active scalar, so its results are carried, not rebuilt.
  template <class U>
  TF24_Strategy<U> rebind_from() const {
    TF24_Strategy<U> out;
    out.state_index = this->state_index;
    out.aux_index = this->aux_index;
    out.birth_rate_x = this->birth_rate_x;
    out.birth_rate_y = this->birth_rate_y;
    out.is_variable_birth_rate = this->is_variable_birth_rate;
    out.collect_all_auxiliary = this->collect_all_auxiliary;
    out.size_0 = this->size_0;
    out.control = this->control;
    out.name = this->name;
    out.extrinsic_drivers = this->extrinsic_drivers;

    TF24_Pars<S> src = pars;
    std::vector<S*> from = src.field_ptrs();
    std::vector<U*> to = out.pars.field_ptrs();
    for (size_t i = 0; i < from.size(); ++i) {
      *to[i] = U(*from[i]);
    }

    out.shading_model_ = shading_model_;
    out.eta_c = U(eta_c);
    out.canopy_shape.initialise(out.pars.eta, out.shading_model_);
    out.height_0 = height_0;
    out.area_leaf_0 = U(area_leaf_0);
    out.leaf = leaf;
    out.storage_gate_width = storage_gate_width;
    out.storage_prod_eps = storage_prod_eps;
    out.newton_tol_abs = newton_tol_abs;
    out.GSS_tol_abs = GSS_tol_abs;
    out.vulnerability_curve_ncontrol = vulnerability_curve_ncontrol;
    out.ci_abs_tol = ci_abs_tol;
    out.ci_niter = ci_niter;
    out.beta_R_H = beta_R_H;
    out.beta_R_V = beta_R_V;
    out.function_integrator = function_integrator;
    out.root_carbon_per_leaf_area_.assign(root_carbon_per_leaf_area_.size(),
                                          U(0.0));
    out.refresh_indices();
    return out;
  }

  // Birth height of a (germinated) seed. Strategy-agnostic accessor used by
  // the templated Individual; here height_0 is derived in prepare_strategy().
  double initial_height() const { return height_0; }

  // Crown shading model, resolved once from control.shading_model in
  // prepare_strategy(). TF24 supports deep-crown, mean-light (its default)
  // and crown-centre; PPA is not available for TF24.
  ShadingModel shading_model_ = ShadingModel::MeanLight;

  // Biological (user-settable) parameters; see TF24_Pars above.
  TF24_Pars<S> pars;

  // Derived / precomputed in prepare_strategy() (NOT user-set) -------------
  S eta_c     = NA_REAL; // crown shape factor, precomputed from pars.eta
  CanopyShape<S> canopy_shape;
  // Height and leaf area of a (germinated) seed
  double height_0  = NA_REAL;
  S area_leaf_0;

  // Embedded leaf hydraulic/photosynthesis sub-model, built in prepare_strategy()
  Leaf leaf;

  // Width of the smooth reserve gate G(r) on growth (#517); small relative to
  // [0,1] so the switch about the growth threshold a_st2 is fairly sharp but
  // differentiable. storage_prod_eps smooths the positive-part of net production
  // (replacing the old hard net>0 growth cutoff) for AD-readiness.
  double storage_gate_width = 0.1;
  double storage_prod_eps   = 1e-4;

  // Solver tolerances and other constants not currently exposed to R
  double newton_tol_abs = 0.001;
  double GSS_tol_abs = 1e-3;
  double vulnerability_curve_ncontrol = 100;
  double ci_abs_tol = 1e-6;
  double ci_niter = 1000;
  // The root-architecture model's two constants. They used to be handed to the
  // Leaf constructor; since phylloptim #33 the leaf takes the RESISTANCES and this
  // strategy owns the model that produces them, which is where they belong -- the
  // 1/3 : 2/3 root split and the dz^2 vertical scaling were never gas exchange.
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
  int aux_idx_assimilation = -1;
  int aux_idx_area_sapwood = -1;       // only present when collect_all_auxiliary
  int state_idx_area_heartwood = -1;
  int state_idx_mass_heartwood = -1;
  int state_idx_storage        = -1;

  // For integrating functions with using Gauss-Kronrod quadrature
  quadrature::QK function_integrator;

  // Reusable per-layer root-carbon buffer, refilled (not reallocated) each
  // net_mass_production_dt call to avoid a heap allocation per derivs eval.
  // Carries S: per-layer root carbon is one of the state directions the leaf's
  // supplied Jacobian has rows for, so it is a live gradient channel.
  std::vector<S> root_carbon_per_leaf_area_;

  // The same numbers with the derivative stripped, because the architecture
  // model and the leaf both take double. Held rather than made per call for the
  // reason the buffer above is.
  std::vector<double> root_carbon_value_;
  std::vector<double> psi_soil_value_;

  // And the resistances derived from it, held the same way and for the same
  // reason (phylloptim #33). Refilled through root_network_from_carbon's in-place
  // overload: building five fresh vectors per solve measured +0.074 us there,
  // about +2% of a whole solve. It must NOT be moved from -- phylloptim's
  // set_root_network takes it by const reference precisely so this buffer keeps
  // its capacity across calls. Stays double: the leaf is a black-box node.
  phylloptim::RootNetwork root_network_;

  // The leaf's two outputs on the active chain, carrying its supplied Jacobian.
  // Written by net_mass_production_dt before compute_rates reads either.
  S leaf_profit_;
  std::vector<S> leaf_soil_consumption_;
};

template <typename S>
typename TF24_Strategy<S>::ptr make_strategy_ptr(TF24_Strategy<S> s);

// --- Hard-coded root-distribution constants (review #9) ---------------------
// Named here for clarity; promotion to user-tunable traits (RcppR6) is a
// deliberate follow-up (see vignettes/models/code_review_leaf_tf24.qmd #9).
// rescales total fine-root mass into the per-layer carbon units expected by the
// root hydraulic network in Leaf::set_physiology.
inline const double root_mass_carbon_scale = 83.26 * 0.5;
// The rooting depth cap moved to TF24_Pars::rooting_depth_max so it can be set
// from R; it is no longer a file-static constant here.

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
template <typename S>
TF24_Strategy<S>::TF24_Strategy() {
  this->collect_all_auxiliary = false;
  // build the string state/aux name to index map
  refresh_indices();
  this->name = "TF24";
}

// not sure 'average' is the right term here..
template <typename S>
S TF24_Strategy<S>::compute_average_light_environment(
    S z, S height, const TF24_Environment<S> &environment) {
// NOTE: the light environment is clamped to a small positive floor (1e-4)
// rather than allowed to reach 0 (original rationale was never recorded;
// preserved as-is).

     using std::max;
     return max(environment.get_environment_at_height(z), S(0.0001)) *
       canopy_shape.q_from_height(z, height);
}

// assumes optimise_psi_stem_TF has been run for optimal psi_stem
template <typename S>
S TF24_Strategy<S>::evapotranspiration_dt(S area_leaf_, int soil_layer) {
  if constexpr (std::is_same_v<S, double>) {
    return leaf.soil_consumption_[soil_layer] * area_leaf_;
  } else {
    return leaf_soil_consumption_[soil_layer] * area_leaf_;
  }
}

// Two of the leaf's outputs re-enter the active chain, and they respond to the
// environment in two different ways.
//
// PROFIT is the objective at the operating point the leaf chose, so by the
// envelope theorem its response is the direct one, at a frozen operating point.
// The leaf supplies it analytically.
//
// PER-LAYER UPTAKE is set as a side effect AT that operating point. It consumes
// the choice rather than being it, so the choice's own movement is part of the
// answer:
//
//     dE_i/du = dE_i/du at a frozen collar  +  (dE_i/dcollar) * (dcollar/du)
//
// and dcollar/du comes from the condition that defines the collar rather than
// from the solve that found it: dprofit/dcollar is zero there, so
// dcollar/du = -(d2profit/dcollar du) / (d2profit/dcollar2).
//
// ⚠️ d2profit/dcollar du CANNOT be had by differencing the leaf's analytic
// dprofit/du in the collar, which is the obvious economy and is wrong by a
// factor of 1.40. That expression prices water at a lambda the collar solve
// equalises, so it is dprofit/du AT the operating point and not away from it.
// It has to come from the marginal profit itself, which is what the
// perturbations below read.
//
// Marginal profit reads the state only through total uptake and through uptake's
// own collar sensitivity, so a pair of scalars serves every direction. The second
// is closed form; the first is solved from ONE perturbation, in root carbon.
//
// The family that perturbation is taken from decides the answer, and a residual
// cannot say so. The soil potentials' pairs of sensitivities are collinear to
// about one part in 10^4, so a pair fitted from them fits every one of them and
// still carries the wrong split between the two scalars into every other
// direction. Root carbon is outside that span, which is why it anchors here.
template <typename S>
template <typename Drive, typename Rebuild>
void TF24_Strategy<S>::record_leaf_outputs(const S& radiation,
                                           const std::vector<S>& psi_soil,
                                           const S& conductance_max,
                                           Drive drive, Rebuild rebuild_roots) {
  using odelia::util::to_passive;
  namespace grad = phylloptim::gradient;

  grad::ProfitEnvDerivatives env;
  grad::profit_env_derivatives(leaf, env);
  if (!env.usable) {
    util::stop("TF24 gradient: the leaf cannot say how profit responds to its "
               "environment here -- " + env.message);
  }
  const size_t n_layer = env.dprofit_dpsi_soil.size();
  // The leaf answers only as deep as it is rooted, while psi_soil runs the whole
  // column. Deeper layers move no water and their response is a true zero.
  if (n_layer > psi_soil.size()) {
    util::stop("TF24 gradient: the leaf answered for " +
               util::to_string(static_cast<int>(n_layer)) + " layers out of " +
               util::to_string(static_cast<int>(psi_soil.size())));
  }

  const double collar = leaf.opt_root_psi_;
  const double radiation_value = to_passive(radiation);
  const std::vector<double>& psi_value = psi_soil_value_;

  // The marginal profit at the collar the solve left, and its slope there.
  const double h_collar = std::max(std::abs(collar), 1.0) * 1e-6;
  const double marginal = leaf.dprofit_droot_collar_psi(collar);
  const double marginal_hi = leaf.dprofit_droot_collar_psi(collar + h_collar);
  const double marginal_lo = leaf.dprofit_droot_collar_psi(collar - h_collar);
  const double curvature = (marginal_hi - marginal_lo) / (2.0 * h_collar);
  if (!util::is_finite(curvature) || curvature >= 0.0) {
    util::stop("TF24 gradient: the leaf's profit has no usable curvature at "
               "this operating point, so the collar's own response has nothing "
               "to stand on");
  }

  // How total uptake and its collar sensitivity respond to each layer, in closed
  // form. These are what the two scalars multiply.
  std::vector<double> dEup_dpsi, d2Eup_dcollar_dpsi, dE_dcollar;
  leaf.dE_from_soil_dpsi_soil(collar, psi_value, dEup_dpsi);
  leaf.roots_.d2uptake_dpsi_dpsi_soil(collar, psi_value, d2Eup_dcollar_dpsi);
  leaf.dE_from_soil_dpsi_collar_by_layer(collar, psi_value, dE_dcollar);
  for (size_t j = 0; j < n_layer; ++j) {
    if (!util::is_finite(dEup_dpsi[j]) || !util::is_finite(d2Eup_dcollar_dpsi[j]) ||
        !util::is_finite(dE_dcollar[j])) {
      util::stop("TF24 gradient: layer " + util::to_string(static_cast<int>(j)) +
                 " sits on a branch kink, so its water response does not exist");
    }
  }

  // One perturbed evaluation: re-supply the leaf, hold the collar where it was,
  // and read the marginal profit there. No solve, so the operating point does
  // not move and this is a partial derivative.
  const double kmax_value = to_passive(conductance_max);
  auto seat_at = [&](double rad, const std::vector<double>& psi,
                     double kmax) -> void {
    drive(rad, psi, kmax);
    leaf.evaluate_root_collar_psi(collar);
  };
  auto marginal_at = [&](double rad, const std::vector<double>& psi) -> double {
    seat_at(rad, psi, kmax_value);
    return leaf.dprofit_droot_collar_psi(collar);
  };

  // The two coefficients the state directions are read through. Both are
  // definite numbers -- the marginal profit's partials in total uptake and in
  // its collar slope -- so any two independent directions recover them in exact
  // arithmetic. Which directions are used decides only where the arithmetic's
  // error lands, and that is not a free choice: the soil potentials are
  // collinear with each other to about one part in 10^4, so a pair fitted from
  // them reproduces their own span and nothing outside it. Root carbon is
  // outside it, and its rows carried the whole of the compensation.
  //
  // So the second is read in closed form and the first is solved along a CARBON
  // direction, which is the family most exposed to the split. The potentials
  // then take whatever error is left, and they are the family that cannot see it.
  // A relative step of 1e-6 does not move total uptake above the solve's own
  // floor, so the largest step that is still local is the best-conditioned one.
  const double fit_step = 1e-3;
  // The leaf reports per-layer uptake in mol and every flux derivative in kg.
  const double to_mol_flux = 1.0 / phylloptim::kg_per_mol_h2o;

  const double b = leaf.dmarginal_profit_duptake_slope();
  if (!util::is_finite(b)) {
    util::stop("TF24 gradient: the stem's marginal profit per unit of collar "
               "conductance is not finite, so the water response has no second "
               "coefficient");
  }

  // The supply side answers every carbon direction analytically, so the only
  // thing this direction costs is the marginal profit either side of it.
  std::vector<std::vector<double>> dE_drc, dD_drc;
  leaf.roots_.duptake_droot_carbon(collar, psi_value, dE_drc, dD_drc);
  std::vector<double> dEup_drc(n_layer, 0.0), d2Eup_drc(n_layer, 0.0);
  for (size_t k = 0; k < n_layer; ++k) {
    if (!(root_carbon_value_[k] > 0.0)) { continue; }
    for (size_t i = 0; i < n_layer; ++i) {
      if (!util::is_finite(dE_drc[i][k]) || !util::is_finite(dD_drc[i][k])) {
        util::stop("TF24 gradient: layer " + util::to_string(static_cast<int>(i)) +
                   " sits on a branch kink, so the root carbon's route into its "
                   "flux does not exist");
      }
      dEup_drc[k] += dE_drc[i][k];
      d2Eup_drc[k] += dD_drc[i][k];
    }
  }
  // Whichever layer's carbon reaches total uptake hardest: the best-conditioned
  // direction of the family it is meant to serve.
  size_t anchor = n_layer;
  double anchor_reach = 0.0;
  for (size_t k = 0; k < n_layer; ++k) {
    if (root_carbon_value_[k] > 0.0 && std::abs(dEup_drc[k]) > anchor_reach) {
      anchor_reach = std::abs(dEup_drc[k]);
      anchor = k;
    }
  }
  if (anchor == n_layer) {
    util::stop("TF24 gradient: no layer's root carbon reaches total uptake, so "
               "the water response has no direction to be read along");
  }

  const double base_rc = root_carbon_value_[anchor];
  const double h_rc = base_rc * fit_step;
  double m_rc_up = 0.0, m_rc_dn = 0.0;
  for (int side = 0; side < 2; ++side) {
    root_carbon_value_[anchor] = base_rc + (side == 0 ? h_rc : -h_rc);
    rebuild_roots();
    seat_at(radiation_value, psi_value, kmax_value);
    (side == 0 ? m_rc_up : m_rc_dn) = leaf.dprofit_droot_collar_psi(collar);
  }
  root_carbon_value_[anchor] = base_rc;
  rebuild_roots();
  seat_at(radiation_value, psi_value, kmax_value);
  const double dR_drc_anchor = (m_rc_up - m_rc_dn) / (2.0 * h_rc);

  const double a = (dR_drc_anchor - b * d2Eup_drc[anchor]) / dEup_drc[anchor];
  if (!util::is_finite(a)) {
    util::stop("TF24 gradient: the first coefficient of the water response is "
               "not finite -- dR/drc=" + util::to_string(dR_drc_anchor) +
               " b=" + util::to_string(b) +
               " dEup/drc=" + util::to_string(dEup_drc[anchor]));
  }

  // The maximum leaf-specific conductance is not in the span the pair above
  // factors: it enters the marginal profit through the stem potential directly,
  // not through total uptake. So it is read along its own direction, at the same
  // step the pair was fitted at, and one pass gives all three of its partials.
  const double h_kmax = std::abs(kmax_value) * fit_step;
  std::vector<double> uptake_up(n_layer), uptake_dn(n_layer);
  double profit_up = 0.0, profit_dn = 0.0, m_kmax_up = 0.0, m_kmax_dn = 0.0;
  for (int side = 0; side < 2; ++side) {
    seat_at(radiation_value, psi_value,
            kmax_value + (side == 0 ? h_kmax : -h_kmax));
    (side == 0 ? profit_up : profit_dn) = leaf.profit_;
    (side == 0 ? m_kmax_up : m_kmax_dn) =
        leaf.dprofit_droot_collar_psi(collar);
    for (size_t i = 0; i < n_layer; ++i) {
      (side == 0 ? uptake_up : uptake_dn)[i] = leaf.soil_consumption_[i];
    }
  }
  seat_at(radiation_value, psi_value, kmax_value);
  const double dprofit_dkmax = (profit_up - profit_dn) / (2.0 * h_kmax);
  const double dR_dkmax = (m_kmax_up - m_kmax_dn) / (2.0 * h_kmax);
  if (!util::is_finite(dprofit_dkmax) || !util::is_finite(dR_dkmax)) {
    util::stop("TF24 gradient: the leaf's response to its maximum conductance "
               "is not finite, so the stem channel has nothing to stand on");
  }
  const double dcollar_dkmax = -dR_dkmax / curvature;

  // The root carbon each layer holds. Every quantity below is analytic: the
  // per-layer supply comes from the map the architecture model defines, which is
  // plain arithmetic with no interpolant in it, and the collar's own response
  // from the two coefficients above.
  std::vector<double> dprofit_drc(n_layer, 0.0), dcollar_drc(n_layer, 0.0);
  std::vector<std::vector<double>> dE_drc_frozen(
      n_layer, std::vector<double>(n_layer, 0.0));
  for (size_t k = 0; k < n_layer; ++k) {
    if (!(root_carbon_value_[k] > 0.0)) {
      continue;   // a layer with no roots moves no water, and that zero is the model
    }
    for (size_t i = 0; i < n_layer; ++i) {
      // soil_consumption_ is in mol and these are the kg quantities, so the
      // frozen row carries the same conversion the collar's does.
      dE_drc_frozen[i][k] = to_mol_flux * dE_drc[i][k];
    }
    const double dR_drc = a * dEup_drc[k] + b * d2Eup_drc[k];
    dcollar_drc[k] = -dR_drc / curvature;
    // At a frozen collar profit sees the carbon only through total uptake and
    // thence the stem potential, and that chain is the second coefficient.
    dprofit_drc[k] = b * dEup_drc[k];
  }

  // The leaf's own traits. It holds them, so the route to most of their rows is
  // to move one and re-solve: two evaluations each, which is what phylloptim's
  // own gradient module pays, and for its reason -- these have no closed form.
  // The order is set_traits(), and every trait goes through it -- including the
  // curve's position, which an identity could move without a rebuild. Measured
  // against a reference that rebuilds the strategy, the identity's row is 1.9e-04
  // off where the rebuild's is 2e-06: the grid is a function of the trait, so a
  // held grid differentiates a different model.
  const int n_leaf_trait = 14;
  double lt[n_leaf_trait] = {
      to_passive(pars.vcmax_25), to_passive(pars.c), to_passive(pars.b),
      to_passive(pars.psi_crit), to_passive(pars.root_c),
      to_passive(pars.root_b), to_passive(pars.root_psi_crit),
      to_passive(pars.beta2), to_passive(pars.jmax_25), to_passive(pars.a),
      to_passive(pars.curv_fact_elec_trans), to_passive(pars.curv_fact_colim),
      to_passive(pars.g1_TF24), to_passive(pars.R_d_25)};
  const S* lt_input[n_leaf_trait] = {
      &pars.vcmax_25, &pars.c, &pars.b, &pars.psi_crit, &pars.root_c,
      &pars.root_b, &pars.root_psi_crit, &pars.beta2, &pars.jmax_25, &pars.a,
      &pars.curv_fact_elec_trans, &pars.curv_fact_colim, &pars.g1_TF24,
      &pars.R_d_25};
  // Which are driven by moving the trait and re-solving. Ten are not, for two
  // different reasons.
  //
  // psi_crit and root_psi_crit set the dry bound of an interval the operating
  // point is strictly inside, so complementary slackness makes their rows zero
  // at an interior optimum -- which is the only kind of point this boundary
  // answers for. They are live at a pin, so a branch that begins answering
  // pinned points has to drive them again rather than read this list.
  //
  // Eight more the leaf answers for directly. beta2 and the cost scale reach
  // profit through the hydraulic cost and nothing else; the other six reach it
  // through assimilation and nothing else, three of them sharing the electron
  // transport as their only route. Every one of the eight has an exactly zero
  // frozen-collar uptake row for the same reason: at a fixed collar a
  // carbon-side trait moves no water.
  //
  // What is left driven is the four vulnerability-curve traits, which have no
  // analytic route at all -- their grid is a function of the trait and the
  // forward model rebuilds it.
  const bool lt_driven[n_leaf_trait] = {
      false, true,  true,  false, true,  true,  false,
      false, false, false, false, false, false, false};
  auto apply_leaf_traits = [&]() -> void {
    leaf.set_traits(lt[0], lt[1], lt[2], lt[3], lt[4], lt[5], lt[6], lt[7],
                    lt[8], lt[9], lt[10], lt[11], lt[12], lt[13]);
  };
  std::vector<double> dprofit_dlt(n_leaf_trait, 0.0),
      dcollar_dlt(n_leaf_trait, 0.0);
  std::vector<std::vector<double>> dE_dlt_frozen(
      n_layer, std::vector<double>(n_leaf_trait, 0.0));
  for (int k = 0; k < n_leaf_trait; ++k) {
    if (!lt_driven[k]) {
      continue;
    }
    const double base_t = lt[k];
    const double h_t = std::max(std::abs(base_t), 1.0) * fit_step;
    double p_up = 0.0, p_dn = 0.0, m_up = 0.0, m_dn = 0.0;
    std::vector<double> u_up(n_layer), u_dn(n_layer);
    for (int side = 0; side < 2; ++side) {
      lt[k] = base_t + (side == 0 ? h_t : -h_t);
      apply_leaf_traits();
      seat_at(radiation_value, psi_value, kmax_value);
      (side == 0 ? p_up : p_dn) = leaf.profit_;
      (side == 0 ? m_up : m_dn) = leaf.dprofit_droot_collar_psi(collar);
      for (size_t i = 0; i < n_layer; ++i) {
        (side == 0 ? u_up : u_dn)[i] = leaf.soil_consumption_[i];
      }
    }
    lt[k] = base_t;
    dprofit_dlt[k] = (p_up - p_dn) / (2.0 * h_t);
    const double dR_dlt = (m_up - m_dn) / (2.0 * h_t);
    if (!util::is_finite(dprofit_dlt[k]) || !util::is_finite(dR_dlt)) {
      util::stop("TF24 gradient: leaf trait " + util::to_string(k) +
                 " moves the leaf by an amount that is not finite");
    }
    dcollar_dlt[k] = -dR_dlt / curvature;
    for (size_t i = 0; i < n_layer; ++i) {
      dE_dlt_frozen[i][k] = (u_up[i] - u_dn[i]) / (2.0 * h_t);
    }
  }
  // Put the leaf back to the traits the value came from.
  apply_leaf_traits();
  seat_at(radiation_value, psi_value, kmax_value);

  // The eight the leaf answers for directly, read once at that operating point.
  {
    const int k_beta2 = 7, k_cost_scale = 12;
    const phylloptim::Leaf::CostTraitRows rows = leaf.cost_trait_rows();
    if (!util::is_finite(rows.dprofit_dbeta2) ||
        !util::is_finite(rows.dprofit_dcost_scale) ||
        !util::is_finite(rows.dmarginal_dbeta2) ||
        !util::is_finite(rows.dmarginal_dcost_scale)) {
      util::stop("TF24 gradient: the leaf's hydraulic cost does not respond "
                 "finitely to the two traits that set it");
    }
    dprofit_dlt[k_beta2] = rows.dprofit_dbeta2;
    dprofit_dlt[k_cost_scale] = rows.dprofit_dcost_scale;
    dcollar_dlt[k_beta2] = -rows.dmarginal_dbeta2 / curvature;
    dcollar_dlt[k_cost_scale] = -rows.dmarginal_dcost_scale / curvature;

    const int k_vcmax = 0, k_jmax = 8, k_a = 9, k_curv_elec = 10,
              k_curv_colim = 11, k_R_d = 13;
    const phylloptim::Leaf::PhotoTraitRows photo = leaf.photo_trait_rows();
    const double photo_rows[12] = {
        photo.dprofit_da, photo.dprofit_dcurv_elec, photo.dprofit_dcurv_colim,
        photo.dprofit_dvcmax_25, photo.dprofit_djmax_25, photo.dprofit_dR_d_25,
        photo.dmarginal_da, photo.dmarginal_dcurv_elec,
        photo.dmarginal_dcurv_colim, photo.dmarginal_dvcmax_25,
        photo.dmarginal_djmax_25, photo.dmarginal_dR_d_25};
    for (int i = 0; i < 12; ++i) {
      if (!util::is_finite(photo_rows[i])) {
        util::stop("TF24 gradient: the leaf's assimilation does not respond "
                   "finitely to the six traits that set it");
      }
    }
    dprofit_dlt[k_a] = photo.dprofit_da;
    dprofit_dlt[k_curv_elec] = photo.dprofit_dcurv_elec;
    dprofit_dlt[k_curv_colim] = photo.dprofit_dcurv_colim;
    dprofit_dlt[k_vcmax] = photo.dprofit_dvcmax_25;
    dprofit_dlt[k_jmax] = photo.dprofit_djmax_25;
    dprofit_dlt[k_R_d] = photo.dprofit_dR_d_25;
    dcollar_dlt[k_a] = -photo.dmarginal_da / curvature;
    dcollar_dlt[k_curv_elec] = -photo.dmarginal_dcurv_elec / curvature;
    dcollar_dlt[k_curv_colim] = -photo.dmarginal_dcurv_colim / curvature;
    dcollar_dlt[k_vcmax] = -photo.dmarginal_dvcmax_25 / curvature;
    dcollar_dlt[k_jmax] = -photo.dmarginal_djmax_25 / curvature;
    dcollar_dlt[k_R_d] = -photo.dmarginal_dR_d_25 / curvature;
  }

  // Profit, carrying the envelope response, and the two channels stay apart
  // because they are different objects: a soil term is a price times a supply
  // derivative, and radiation moves no water at all.
  std::vector<odelia::input_and_derivative<S>> against;
  against.reserve(2 + 2 * n_layer + n_leaf_trait);
  against.push_back({radiation, env.dprofit_dlight});
  against.push_back({conductance_max, dprofit_dkmax});
  for (size_t a = 0; a < n_layer; ++a) {
    against.push_back({root_carbon_per_leaf_area_[a], dprofit_drc[a]});
  }
  for (int k = 0; k < n_leaf_trait; ++k) {
    against.push_back({*lt_input[k], dprofit_dlt[k]});
  }
  for (size_t j = 0; j < n_layer; ++j) {
    against.push_back({psi_soil[j], env.dprofit_dpsi_soil[j]});
  }
  leaf_profit_ = odelia::record_with_derivatives<S>(leaf.profit_, against);

  // The collar's own response, then uptake. dE_i/d(radiation) at a frozen collar
  // is exactly zero -- light moves no water -- so the whole of the light channel
  // into the water arrives through the collar.
  const double dR_dlight =
      (marginal_at(radiation_value * (1.0 + 1e-6), psi_value) -
       marginal_at(radiation_value * (1.0 - 1e-6), psi_value)) /
      (2.0 * radiation_value * 1e-6);
  drive(radiation_value, psi_value, kmax_value);
  leaf.evaluate_root_collar_psi(collar);
  static_cast<void>(marginal);

  if (!util::is_finite(dR_dlight)) {
    util::stop("TF24 gradient: the marginal profit's response to radiation is "
               "not finite, so the collar's light response has nothing to "
               "stand on");
  }
  const double dcollar_dlight = -dR_dlight / curvature;
  std::vector<double> dcollar_dpsi(n_layer, 0.0);
  for (size_t j = 0; j < n_layer; ++j) {
    dcollar_dpsi[j] =
        -(a * dEup_dpsi[j] + b * d2Eup_dcollar_dpsi[j]) / curvature;
  }

  leaf_soil_consumption_.assign(psi_soil.size(), S(0.0));
  // The value below is soil_consumption_, in mol; every partial above it is
  // built from the leaf's kg-based flux accessors. The two are different units
  // by design, so the rows carry the same conversion the value already did.
  const double to_mol = to_mol_flux;
  for (size_t i = 0; i < n_layer; ++i) {
    against.clear();
    against.push_back({radiation, to_mol * dE_dcollar[i] * dcollar_dlight});
    // soil_consumption_ is already the mol quantity, so its own difference needs
    // no conversion; dE_dcollar is the kg one and does.
    const double dE_dkmax_frozen =
        (uptake_up[i] - uptake_dn[i]) / (2.0 * h_kmax);
    against.push_back({conductance_max,
                       dE_dkmax_frozen +
                         to_mol * dE_dcollar[i] * dcollar_dkmax});
    for (size_t a = 0; a < n_layer; ++a) {
      against.push_back({root_carbon_per_leaf_area_[a],
                         dE_drc_frozen[i][a] +
                           to_mol * dE_dcollar[i] * dcollar_drc[a]});
    }
    for (int k = 0; k < n_leaf_trait; ++k) {
      against.push_back({*lt_input[k],
                         dE_dlt_frozen[i][k] +
                           to_mol * dE_dcollar[i] * dcollar_dlt[k]});
    }
    for (size_t j = 0; j < n_layer; ++j) {
      const double frozen = (i == j) ? dEup_dpsi[i] : 0.0;
      against.push_back({psi_soil[j],
                         to_mol * (frozen + dE_dcollar[i] * dcollar_dpsi[j])});
    }
    // The leaf reports uptake in mol and E_up in kg; soil_consumption_ is the
    // mol one, which is what the patch water balance reads.
    leaf_soil_consumption_[i] =
        odelia::record_with_derivatives<S>(leaf.soil_consumption_[i], against);
  }
  // Layers below the deepest rooted one draw nothing and respond to nothing, and
  // that zero is the model rather than a missing answer.
  for (size_t i = n_layer; i < psi_soil.size(); ++i) {
    leaf_soil_consumption_[i] = S(0.0);
  }
}

template <typename S>
void TF24_Strategy<S>::refresh_indices () {
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
  aux_idx_assimilation          = this->aux_index.at("assimilation");
  // area_sapwood is only registered when collect_all_auxiliary is set.
  aux_idx_area_sapwood = this->aux_index.count("area_sapwood") ? this->aux_index.at("area_sapwood") : -1;
  state_idx_area_heartwood      = this->state_index.at("area_heartwood");
  state_idx_mass_heartwood      = this->state_index.at("mass_heartwood");
  state_idx_storage             = this->state_index.at("storage");
}

// [eqn 2] area_leaf (inverse of [eqn 3])
template <typename S>
S TF24_Strategy<S>::area_leaf(S height) const {
  return pow(height / pars.a_l1, 1.0 / pars.a_l2);
}

// [eqn 1] mass_leaf (inverse of [eqn 2])
template <typename S>
S TF24_Strategy<S>::mass_leaf(S area_leaf) const {
  return area_leaf * pars.lma;
}

// [eqn 4] area and mass of sapwood
template <typename S>
S TF24_Strategy<S>::area_sapwood(S area_leaf) const {
  return area_leaf * pars.theta;
}

template <typename S>
S TF24_Strategy<S>::mass_sapwood(S area_sapwood, S height) const {
  return area_sapwood * height * eta_c * pars.rho;
}

// [eqn 5] area and mass of bark
template <typename S>
S TF24_Strategy<S>::area_bark(S area_leaf) const {
  return pars.a_b1 * area_leaf * pars.theta;
}

template <typename S>
S TF24_Strategy<S>::mass_bark(S area_bark, S height) const {
  return area_bark * height * eta_c * pars.rho;
}

template <typename S>
S TF24_Strategy<S>::area_stem(S area_bark, S area_sapwood,
                            S area_heartwood) const {
  return area_bark + area_sapwood + area_heartwood;
}

template <typename S>
S TF24_Strategy<S>::diameter_stem(S area_stem) const {
  using std::sqrt;
  return sqrt(4 * area_stem / M_PI);
}

// [eqn 7] Mass of (fine) roots
template <typename S>
S TF24_Strategy<S>::mass_root(S area_leaf) const {
  return pars.a_r1 * area_leaf;
}

// [eqn 8] Total mass
template <typename S>
S TF24_Strategy<S>::mass_live(S mass_leaf, S mass_bark,
                           S mass_sapwood, S mass_root) const {
  return mass_leaf + mass_sapwood + mass_bark + mass_root;
}

template <typename S>
S TF24_Strategy<S>::mass_total(S mass_leaf, S mass_bark,
                            S mass_sapwood, S mass_heartwood,
                            S mass_root) const {
  return mass_leaf + mass_bark + mass_sapwood +  mass_heartwood + mass_root;
}

template <typename S>
S TF24_Strategy<S>::mass_above_ground(S mass_leaf, S mass_bark,
                            S mass_sapwood, S mass_heartwood) const {
  return mass_leaf + mass_bark + mass_sapwood + mass_heartwood;
}

// for updating auxiliary state
template <typename S>
void TF24_Strategy<S>::update_dependent_aux(const int index, Internals<S>& vars) {
  if (index == HEIGHT_INDEX) {
    S height = vars.state(HEIGHT_INDEX);
    vars.set_aux(aux_idx_competition_effect, area_leaf(height));
    vars.set_aux(aux_idx_height_inverse, 1.0 / height);
  }
}


// one-shot update of the scm variables
// i.e. setting rates of ode vars from the state and updating aux vars
template <typename S>
void TF24_Strategy<S>::compute_rates(const TF24_Environment<S>& environment,  Internals<S>& vars) {
  S height = vars.state(HEIGHT_INDEX);
  S area_leaf_ = vars.aux(aux_idx_competition_effect);

  const S net_mass_production_dt_ =
    net_mass_production_dt(environment, height, area_leaf_,
                           vars.aux(aux_idx_height_inverse));

  // store the aux sate
  vars.set_aux(aux_idx_net_mass_production_dt, net_mass_production_dt_);
  vars.set_aux(aux_idx_root_mass, mass_root(area_leaf_));
  vars.set_aux(aux_idx_opt_psi_stem, leaf.opt_psi_stem_);
  vars.set_aux(aux_idx_opt_root_psi, leaf.opt_root_psi_);
  vars.set_aux(aux_idx_transpiration, leaf.transpiration_);
  vars.set_aux(aux_idx_E_up, leaf.E_up_);
  vars.set_aux(aux_idx_profit, leaf.profit_);
  vars.set_aux(aux_idx_stom_cond_CO2, leaf.stom_cond_CO2_);
  vars.set_aux(aux_idx_assimilation, leaf.assim_colimited_);




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

  // --- NSC storage pool (#517) ---------------------------------------------
  // A storage pool buffers demography against short-term productivity swings:
  // growth is gated on having ample reserves, and mortality reads the buffered
  // relative reserves rather than instantaneous net production -- so a trough
  // draws reserves down and death is gradual, instead of the growth cutoff and
  // ~1e32 mortality spike that caused the #550 blow-up.
  const S storage     = std::max(vars.state(state_idx_storage), S(0.0));
  const S storage_max = storage_capacity(area_leaf_, height);
  using std::exp;
  using std::min;
  using std::sqrt;
  const S r           = storage_max > 0.0 ? min(storage / storage_max, S(1.0)) : S(0.0);
  // Reserve-gated growth (#517), following Daniel's intuition that a plant
  // should not grow unless it has ample carbon in storage. Growth and
  // reproduction proceed at the *production* rate, but scaled by a smooth gate
  // G(r) that is ~0 at low relative reserves and ~1 near capacity (logistic
  // centred on the growth threshold a_st2). Carbon not spent on growth refills
  // storage (dS/dt below), so: full reserves -> grow at production rate (healthy
  // dynamics preserved); low reserves -> redirect carbon to refill, growth
  // pauses; drought (net<0) -> reserves drain, growth halts and only resumes
  // once they refill (buffered growth AND mortality). Note growth is gated at
  // the production rate rather than metered out of the pool as a_st2*S: with
  // sapwood-scaled capacity, storage is tiny relative to seedling productivity,
  // so a discharge-rate limit would choke establishment.
  const S G = 1.0 / (1.0 + exp(-(r - pars.a_st2) / storage_gate_width));
  // Smooth positive part of net production (replaces the old hard net>0 cutoff).
  const S P = net_mass_production_dt_;
  const S Ppos =
    0.5 * (P + sqrt(P * P + storage_prod_eps * storage_prod_eps));
  const S growth_flux = Ppos * G;

  const S fraction_allocation_reproduction_ = fraction_allocation_reproduction(height);
  const S darea_leaf_dmass_live_ = darea_leaf_dmass_live(area_leaf_);
  const S fraction_allocation_growth_ = fraction_allocation_growth(height);
  const S area_leaf_dt = growth_flux * fraction_allocation_growth_ * darea_leaf_dmass_live_;

  vars.set_rate(HEIGHT_INDEX, dheight_darea_leaf(area_leaf_) * area_leaf_dt);
  vars.set_rate(FECUNDITY_INDEX,
    fecundity_dt(growth_flux, fraction_allocation_reproduction_));

  // Sapwood -> heartwood conversion is turnover-driven, so it proceeds
  // regardless of carbon status (previously gated behind net>0).
  vars.set_rate(state_idx_area_heartwood, area_heartwood_dt(area_leaf_));
  const S area_sapwood_ = area_sapwood(area_leaf_);
  const S mass_sapwood_ = mass_sapwood(area_sapwood_, height);
  vars.set_rate(state_idx_mass_heartwood, mass_heartwood_dt(mass_sapwood_));

  if (this->collect_all_auxiliary) {
    vars.set_aux(aux_idx_area_sapwood, area_sapwood_);
  }

  // Storage dynamics: dS/dt = net production - carbon spent on growth. When
  // production exceeds what the reserve gate lets through to growth, the surplus
  // charges storage; when it falls short (or net production is negative) storage
  // is drawn down. The net outflow is gated to vanish as S -> 0, flooring
  // storage at zero so relative reserves r stay in [0,1] and the storage-based
  // mortality stays bounded (the structural fix for #550). At the S~0 starvation
  // boundary the ungated part of the deficit is untracked (no worse than the
  // original model, which retained structure under net<0); by then the plant is
  // dying at the bounded maximum mortality anyway.
  const S net_flux = P - growth_flux;
  const S gate_ref = 1e-3 * storage_max;           // ~0.1% of capacity
  const S floor_gate =
    (storage + gate_ref) > 0.0 ? storage / (storage + gate_ref) : S(0.0);
  vars.set_rate(state_idx_storage,
                net_flux > 0.0 ? net_flux : floor_gate * net_flux);

  // [eqn 21] - Instantaneous mortality rate, now driven by relative reserves r.
  vars.set_rate(MORTALITY_INDEX,
      mortality_dt(r, vars.state(MORTALITY_INDEX)));

}

// [eqn 12] Gross annual CO2 assimilation (!!not in use for TF24 model!!)
template <typename S>
S TF24_Strategy<S>::assimilation(const TF24_Environment<S>& environment,
                                    S height,
                                    S area_leaf) {


  S A = 0.0;

  // Define an anonymous function to integrate
  // For given height in crown, take photosynthesis at depth multipled by 
  //   amount of leaf at that depth
  std::function<S(S)> f = [&](S z) -> S {
    return assimilation_leaf(environment.get_environment_at_height(z)) *
      canopy_shape.q_from_height(z, height);
  };

  // Integrate over crown depth using using Gauss-Kronrod quadrature.
  // The number of points used in the integration is determined by the control parameter
  // function_integration_rule. Rules defined in qk_rules.cpp
  A = function_integrator.integrate(f, S(0.0), height);

  return area_leaf * A;
}

// Photosynthetic rate per leaf area
// `x` is openness, ranging from 0 to 1.
template <typename S>
S TF24_Strategy<S>::assimilation_leaf(S x) const {
  return pars.a_p1 * x / (x + pars.a_p2);
}

// [eqn 13] Total maintenance respiration
// NOTE: In contrast with Falster ref model, we do not normalise by pars.a_y*pars.a_bio.
template <typename S>
S TF24_Strategy<S>::respiration(S mass_leaf, S mass_sapwood,
                             S mass_bark, S mass_root) const {
  return respiration_leaf(mass_leaf) +
         respiration_bark(mass_bark) +
         respiration_sapwood(mass_sapwood) +
         respiration_root(mass_root);
}

template <typename S>
S TF24_Strategy<S>::respiration_leaf(S mass) const {
  return pars.r_l * mass;
}

template <typename S>
S TF24_Strategy<S>::respiration_bark(S mass) const {
  return pars.r_b * mass;
}

template <typename S>
S TF24_Strategy<S>::respiration_sapwood(S mass) const {
  return pars.r_s * mass;
}

template <typename S>
S TF24_Strategy<S>::respiration_root(S mass) const {
  return pars.r_r * mass;
}

// [eqn 14] Total turnover
template <typename S>
S TF24_Strategy<S>::turnover(S mass_leaf, S mass_bark,
                          S mass_sapwood, S mass_root) const {
   return turnover_leaf(mass_leaf) +
          turnover_bark(mass_bark) +
          turnover_sapwood(mass_sapwood) +
          turnover_root(mass_root);
}

template <typename S>
S TF24_Strategy<S>::turnover_leaf(S mass) const {
  return pars.k_l * mass;
}

template <typename S>
S TF24_Strategy<S>::turnover_bark(S mass) const {
  return pars.k_b * mass;
}

template <typename S>
S TF24_Strategy<S>::turnover_sapwood(S mass) const {
  return pars.k_s * mass;
}

template <typename S>
S TF24_Strategy<S>::turnover_root(S mass) const {
  return pars.k_r * mass;
}

// [eqn 15] Net production
//
// NOTE: Translation of variable names from the Falster 2011.  Everything
// before the minus sign is SCM's N, our `net_mass_production_dt` is SCM's P.
template <typename S>
S TF24_Strategy<S>::net_mass_production_dt_A(S assimilation, S respiration,
                                S turnover) const {
  return pars.a_bio * pars.a_y * (assimilation - respiration) - turnover;
}

// One shot calculation of net_mass_production_dt
// Used by establishment_probability() and compute_rates().
template <typename S>
S TF24_Strategy<S>::net_mass_production_dt(const TF24_Environment<S>& environment,
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



  // The radiation that drives the leaf optimisation depends on the shading
  // model and is computed below (just before the optimisation), once the
  // depth-independent inputs are ready.

  // psi_soil (-MPa), computed once per soil state and cached in environment.
  const std::vector<S>& psi_soil = environment.get_soil_water_potential_state();
  
// find leaf specific max hydraulic conductance (kg m^-2 LA s^-1 MPa ^-1)
  // pars.K_s: max hydraulic conductivity (kg m^-2 s^-1 MPa^-1),
  // pars.theta: huber value
  // eta_c: accounts for average position of leaf mass
  // height: maximum plant height
  const S leaf_specific_conductance_max = pars.K_s * pars.theta / (height * eta_c);

  // Sapwood volume per leaf area (pars.theta * height * eta_c) used to be handed
  // to the leaf, which stored it and never read it. Recompute it here if a
  // caller ever needs it.

  // ----------------------------------------------------------------------
  // ROOT MASS DISTRIBUTION ACROSS SOIL LAYERS
  // ----------------------------------------------------------------------
  // Fine-root carbon is distributed over depth using the same cumulative shape
  // function Q() used for the leaf canopy, but parameterised over soil depth
  // instead of crown height. Q(z, rooting_depth, 0.2) gives the fraction of
  // roots *below* depth z, so layer a holds root_mass_scale * (Q(z_{a-1}) -
  // Q(z_a)). rooting_depth is capped at the soil column depth and the loop
  // breaks once Q reaches 0, so deeper layers stay at zero.
  //
  // The quantity is carbon PER UNIT LEAF AREA, and that is a hard requirement
  // rather than a convenience: the leaf is intensive, it takes resistances built
  // from this, and a network built from absolute carbon is five vectors of
  // positive numbers that no check on the far side can distinguish -- uptake
  // then comes back wrong by the leaf area with nothing raised. mass_root() is
  // strictly linear in area_leaf, so dividing it back out is exactly
  // root_mass_carbon_scale * pars.a_r1 and costs no arithmetic.
  //
  // Reuse the member buffer (assign refills + zeroes without reallocating when
  // the layer count is unchanged); zeroing matters because the loop below breaks
  // early below the rooting depth, leaving deep layers that must read as 0.
  // TODO (perf): the scale (83.26) is hard-coded and should become a trait.
  root_carbon_per_leaf_area_.assign(soil_number_of_depths_, 0.0);
  root_carbon_value_.assign(soil_number_of_depths_, 0.0);

  const S rooting_depth = std::min(height, pars.rooting_depth_max);
  const S root_mass_scale = root_mass_carbon_scale * pars.a_r1;

  {
    using odelia::util::to_passive;
    S prev_q = 1.0;
    for (int a = 0; a < soil_number_of_depths_; ++a) {
      if (prev_q == 0) {
        break;
      }
      const S q = Q(soil_depths_[a], rooting_depth, pars.root_depth_shape_eta);
      root_carbon_per_leaf_area_[a] = root_mass_scale * (prev_q - q);
      root_carbon_value_[a] = to_passive(root_carbon_per_leaf_area_[a]);
      prev_q = q;
    }
  }

  // Carbon -> resistance: the root architecture model runs here, on the plant
  // side. The leaf's supply solve reads r_R_H_min and r_R_V_sum and knows
  // nothing about root carbon, the horizontal/vertical split or the layer
  // thickness, so this strategy owns the model exactly as it already owns the
  // conductance-versus-height model above, and hands over the reduced quantity.
  // layer_thickness is the shared definition of dz -- do not open-code it, since
  // the vertical resistance scales with dz^2 and the two sides drifting apart
  // would be a silent squared factor neither side could detect.
  phylloptim::root_network_from_carbon(
      root_carbon_value_, phylloptim::layer_thickness(soil_depths_), beta_R_H,
      beta_R_V, root_network_);

  // Reuse geometry precomputed by environment; avoids rebuilding z midpoints each call.
  leaf.roots_.z_soil_mid_ = environment.get_soil_mid_depths();
  leaf.roots_.use_precomputed_z_soil_mid_ = true;

  // Per-timestep above-canopy wind, read only on the energy-balance path.
  leaf.wind_speed_ = environment.get_wind_speed();

  // Optimise the leaf at a given absorbed radiation: rebuilds physiology and
  // solves the root-collar water potential, leaving the leaf.* outputs
  // (profit_, transpiration_, soil_consumption_, opt_psi_stem_, ...) set. Only
  // the radiation argument varies between calls; every other input is
  // depth-independent and already computed above.
  //
  // Leaf carries double, so an active strategy hands it the values of its
  // inputs and reads its outputs back as constants: d(rates)/d(leaf inputs) is
  // exactly zero here, until the leaf's local Jacobian is injected across this
  // same boundary. Each leaf output enters the active chain at one place -- the
  // aux stores and leaf.profit_ in net_mass_production_dt, and
  // leaf.soil_consumption_ in evapotranspiration_dt -- so a partial attaches to
  // one expression per output.
  S radiation_used = 0.0;
  auto optimise_at = [&](const S& radiation) -> void {
    radiation_used = radiation;
    if constexpr (std::is_same_v<S, double>) {
      leaf.set_physiology(root_network_, radiation, psi_soil, soil_depths_,
                          leaf_specific_conductance_max,
                          environment.get_atm_vpd(), environment.get_ca(),
                          environment.get_leaf_temp(),
                          environment.get_atm_o2_kpa(),
                          environment.get_atm_kpa());
    } else {
      using odelia::util::to_passive;
      psi_soil_value_.resize(psi_soil.size());
      for (size_t a = 0; a < psi_soil.size(); ++a) {
        psi_soil_value_[a] = to_passive(psi_soil[a]);
      }
      leaf.set_physiology(root_network_, to_passive(radiation), psi_soil_value_,
                          soil_depths_,
                          to_passive(leaf_specific_conductance_max),
                          environment.get_atm_vpd(), environment.get_ca(),
                          environment.get_leaf_temp(),
                          environment.get_atm_o2_kpa(),
                          environment.get_atm_kpa());
    }
    solve_leaf();
  };

  // Convert canopy openness (0-1) into absorbed radiation: PPFD attenuated by
  // the self-shading coefficient pars.k_I. The light floor (1e-4) matches
  // compute_average_light_environment().
  const double PPFD = environment.get_PPFD();
  auto radiation_at = [&](S light) -> S {
    return pars.k_I * std::max(light, S(0.0001)) * PPFD;
  };

  // Aggregate the leaf submodel over the crown according to the shading model.
  // The expensive hydraulic optimisation is the unit of work here, so the model
  // choice is about how many times it runs and on what light:
  //  - crown-centre:  one optimisation at the crown-centre light.
  //  - mean-light:    one optimisation at the leaf-area-weighted mean light
  //                   (TF24's established default).
  //  - deep-crown:    one optimisation per crown-depth quadrature point, with
  //                   every leaf output integrated to a leaf-area-weighted mean.
  if (shading_model_ == ShadingModel::CrownCentre) {
    optimise_at(radiation_at(environment.get_environment_at_height(height * eta_c)));
  } else if (shading_model_ == ShadingModel::MeanLight) {
    // Leaf-area-weighted mean canopy openness = integral of (light * q) over the
    // crown (q integrates to one). radiation_at then applies pars.k_I * PPFD, exactly
    // reproducing TF24's established average_radiation.
    auto f = [&](S x) -> S {
      return compute_average_light_environment(x, height, environment);
    };
    optimise_at(radiation_at(function_integrator.integrate(f, S(0.0), height)));
  } else { // DeepCrown
    if constexpr (std::is_same_v<S, double>) {
      const std::vector<double> nodes =
        function_integrator.integrate_vector_x(0.0, height);
      const size_t nn = nodes.size();
      std::vector<double> profit_y(nn), trans_y(nn), eup_y(nn), psi_y(nn),
        root_psi_y(nn), gco2_y(nn), assim_y(nn);
      std::vector<std::vector<double>> soil_y(
        soil_number_of_depths_, std::vector<double>(nn));
      for (size_t i = 0; i < nn; ++i) {
        const S qi = canopy_shape.q_from_height(nodes[i], height);
        optimise_at(radiation_at(environment.get_environment_at_height(nodes[i])));
        profit_y[i]   = leaf.profit_ * qi;
        trans_y[i]    = leaf.transpiration_ * qi;
        eup_y[i]      = leaf.E_up_ * qi;
        psi_y[i]      = leaf.opt_psi_stem_ * qi;
        root_psi_y[i] = leaf.opt_root_psi_ * qi;
        gco2_y[i]     = leaf.stom_cond_CO2_ * qi;
        assim_y[i]    = leaf.assim_colimited_ * qi;
        for (int a = 0; a < soil_number_of_depths_; ++a) {
          soil_y[a][i] = leaf.soil_consumption_[a] * qi;
        }
      }
      // Integrate each leaf output to its leaf-area-weighted crown mean (q
      // integrates to one over the crown). soil_consumption_ feeds the patch
      // water balance, so it must be the depth-integrated total; the rest are
      // diagnostics reported through compute_rates.
      leaf.profit_          = function_integrator.integrate_vector(profit_y, 0.0, height);
      leaf.transpiration_   = function_integrator.integrate_vector(trans_y, 0.0, height);
      leaf.E_up_            = function_integrator.integrate_vector(eup_y, 0.0, height);
      leaf.opt_psi_stem_    = function_integrator.integrate_vector(psi_y, 0.0, height);
      leaf.opt_root_psi_    = function_integrator.integrate_vector(root_psi_y, 0.0, height);
      leaf.stom_cond_CO2_   = function_integrator.integrate_vector(gco2_y, 0.0, height);
      leaf.assim_colimited_ = function_integrator.integrate_vector(assim_y, 0.0, height);
      for (int a = 0; a < soil_number_of_depths_; ++a) {
        leaf.soil_consumption_[a] =
          function_integrator.integrate_vector(soil_y[a], 0.0, height);
      }
    } else {
      util::stop("shading_model '" + this->control.shading_model +
                 "' is not differentiable: its crown means pass through Leaf, "
                 "which carries double. Use crown-centre or mean-light "
                 "shading for a gradient.");
    }
  }


  //TODO: one point constant ratio and integral width for daylength
  // convert assimilation per leaf area per second (umol m^-2 s^-1) to canopy-level total yearly assimilation (mol yr^-1)
  // converts to canopy area, then years, then mols
  S profit_ = leaf.profit_;
  if constexpr (!std::is_same_v<S, double>) {
    // Re-supplying the leaf is the caller's job, because the caller is what owns
    // these inputs; every one of them crosses as a value.
    auto drive = [&](double rad, const std::vector<double>& psi,
                     double kmax) -> void {
      leaf.set_physiology(root_network_, rad, psi, soil_depths_, kmax,
                          environment.get_atm_vpd(), environment.get_ca(),
                          environment.get_leaf_temp(),
                          environment.get_atm_o2_kpa(),
                          environment.get_atm_kpa());
    };
    // The architecture model is this strategy's, so the leaf boundary is handed
    // a rebuild rather than the layer thickness and the two betas.
    auto rebuild_roots = [&]() -> void {
      phylloptim::root_network_from_carbon(
          root_carbon_value_, phylloptim::layer_thickness(soil_depths_),
          beta_R_H, beta_R_V, root_network_);
    };
    record_leaf_outputs(radiation_used, psi_soil,
                        leaf_specific_conductance_max, drive, rebuild_roots);
    profit_ = leaf_profit_;
  }
  const S assimilation_ = profit_ * area_leaf_* 60*60*12*365/1e6;
  // const double assimilation_ = assimilation(environment, height, area_leaf_);
  const S respiration_ =
    respiration(mass_leaf_, mass_sapwood_, mass_bark_, mass_root_);
  const S turnover_ =
    turnover(mass_leaf_, mass_bark_, mass_sapwood_, mass_root_);
  return net_mass_production_dt_A(assimilation_, respiration_, turnover_);
}

// Base TF24: optimise the root-collar water potential from scratch each call.
template <typename S>
void TF24_Strategy<S>::solve_leaf() {
  leaf.find_root_collar_psi();
}

// [eqn 16] Fraction of production allocated to reproduction
template <typename S>
S TF24_Strategy<S>::fraction_allocation_reproduction(S height) const {
  return pars.a_f1 / (1.0 + exp(pars.a_f2 * (1.0 - height / pars.hmat)));
}

// Fraction of production allocated to growth
template <typename S>
S TF24_Strategy<S>::fraction_allocation_growth(S height) const {
  return 1.0 - fraction_allocation_reproduction(height);
}

// [eqn 17] Rate of offspring production
template <typename S>
S TF24_Strategy<S>::fecundity_dt(S net_mass_production_dt,
                               S fraction_allocation_reproduction) const {
  return net_mass_production_dt * fraction_allocation_reproduction /
    (pars.omega + pars.a_f3);
}

template <typename S>
S TF24_Strategy<S>::darea_leaf_dmass_live(S area_leaf) const {
  return 1.0/(  dmass_leaf_darea_leaf(area_leaf)
              + dmass_sapwood_darea_leaf(area_leaf)
              + dmass_bark_darea_leaf(area_leaf)
              + dmass_root_darea_leaf(area_leaf));
}

template <typename S>
S TF24_Strategy<S>::dheight_darea_leaf(S area_leaf) const {
  return pars.a_l1 * pars.a_l2 * pow(area_leaf, pars.a_l2 - 1);
}

// Mass of leaf needed for new unit area leaf, d m_s / d a_l
template <typename S>
S TF24_Strategy<S>::dmass_leaf_darea_leaf(S /* area_leaf */) const {
  return pars.lma;
}

// Mass of stem needed for new unit area leaf, d m_s / d a_l
template <typename S>
S TF24_Strategy<S>::dmass_sapwood_darea_leaf(S area_leaf) const {
  return pars.rho * eta_c * pars.a_l1 * pars.theta * (pars.a_l2 + 1.0) * pow(area_leaf, pars.a_l2);
}

// Mass of bark needed for new unit area leaf, d m_b / d a_l
template <typename S>
S TF24_Strategy<S>::dmass_bark_darea_leaf(S area_leaf) const {
  return pars.a_b1 * dmass_sapwood_darea_leaf(area_leaf);
}

// Mass of root needed for new unit area leaf, d m_r / d a_l
template <typename S>
S TF24_Strategy<S>::dmass_root_darea_leaf(S /* area_leaf */) const {
  return pars.a_r1;
}

// Growth rate of basal diameter_stem per unit time
template <typename S>
S TF24_Strategy<S>::ddiameter_stem_darea_stem(S area_stem) const {
  return pow(M_PI * area_stem, -0.5);
}

// Growth rate of sapwood area at base per unit time
template <typename S>
S TF24_Strategy<S>::area_sapwood_dt(S area_leaf_dt) const {
  return area_leaf_dt * pars.theta;
}

// Note, unlike others, heartwood growth does not depend on leaf area growth, but
// rather existing sapwood
template <typename S>
S TF24_Strategy<S>::area_heartwood_dt(S area_leaf) const {
  return pars.k_s * area_sapwood(area_leaf);
}

// Growth rate of bark area at base per unit time
template <typename S>
S TF24_Strategy<S>::area_bark_dt(S area_leaf_dt) const {
  return pars.a_b1 * area_leaf_dt * pars.theta;
}

// Growth rate of stem basal area per unit time
template <typename S>
S TF24_Strategy<S>::area_stem_dt(S area_leaf,
                               S area_leaf_dt) const {
  return area_sapwood_dt(area_leaf_dt) +
    area_bark_dt(area_leaf_dt) +
    area_heartwood_dt(area_leaf);
}

// Growth rate of basal diameter_stem per unit time
template <typename S>
S TF24_Strategy<S>::diameter_stem_dt(S area_stem, S area_stem_dt) const {
  return ddiameter_stem_darea_stem(area_stem) * area_stem_dt;
}

// Growth rate of root mass per unit time
template <typename S>
S TF24_Strategy<S>::mass_root_dt(S area_leaf,
                               S area_leaf_dt) const {
  return area_leaf_dt * dmass_root_darea_leaf(area_leaf);
}

template <typename S>
S TF24_Strategy<S>::mass_live_dt(S fraction_allocation_reproduction,
                               S net_mass_production_dt) const {
  return (1 - fraction_allocation_reproduction) * net_mass_production_dt;
}

template <typename S>
S TF24_Strategy<S>::mass_total_dt(S fraction_allocation_reproduction,
                                     S net_mass_production_dt,
                                     S mass_heartwood_dt) const {
  return mass_live_dt(fraction_allocation_reproduction, net_mass_production_dt) +
    mass_heartwood_dt;
}

// TODO: Do we not track root mass change?
template <typename S>
S TF24_Strategy<S>::mass_above_ground_dt(S area_leaf,
                                       S fraction_allocation_reproduction,
                                       S net_mass_production_dt,
                                       S mass_heartwood_dt,
                                       S area_leaf_dt) const {
  const S mass_root_dt =
    area_leaf_dt * dmass_root_darea_leaf(area_leaf);
  return mass_total_dt(fraction_allocation_reproduction, net_mass_production_dt,
                        mass_heartwood_dt) - mass_root_dt;
}

template <typename S>
S TF24_Strategy<S>::mass_heartwood_dt(S mass_sapwood) const {
  return turnover_sapwood(mass_sapwood);
}


template <typename S>
S TF24_Strategy<S>::mass_live_given_height(S height) const {
  S area_leaf_ = area_leaf(height);
  return mass_leaf(area_leaf_) +
         mass_bark(area_bark(area_leaf_), height) +
         mass_sapwood(area_sapwood(area_leaf_), height) +
         mass_root(area_leaf_);
}

template <typename S>
S TF24_Strategy<S>::height_given_mass_leaf(S mass_leaf) const {
  return pars.a_l1 * pow(mass_leaf / pars.lma, pars.a_l2);
}

template <typename S>
S TF24_Strategy<S>::mortality_dt(S relative_reserves,
                              S cumulative_mortality) const {

  // Growth-dependent mortality is now driven by relative NSC reserves
  // r = S/S_max (in [0,1]) rather than instantaneous productivity, so the rate
  // is bounded (see mortality_storage_dependent_dt): death becomes gradual as
  // reserves deplete instead of spiking to ~1e32 under carbon deficit (#550).
  if (util::is_finite(cumulative_mortality)) {
    return
      mortality_growth_independent_dt() +
      mortality_storage_dependent_dt(relative_reserves);
 } else {
    // If mortality probability is 1 (latency = Inf) then the rate
    // calculations break.  Setting them to zero gives the correct
    // behaviour.
    return 0.0;
  }
}

template <typename S>
S TF24_Strategy<S>::mortality_growth_independent_dt() const {
  return pars.d_I;
}

// Storage-dependent growth mortality (#517), following Stefaniak et al. 2026
// (Eq 6). Bounded in [a_dG1*exp(-a_dG2), a_dG1] for r in [0,1]: full reserves
// (r=1) give near-zero excess mortality; empty reserves (r=0) give the finite
// maximum a_dG1. This boundedness is what removes the #550 ODE overflow.
template <typename S>
S TF24_Strategy<S>::mortality_storage_dependent_dt(S relative_reserves) const {
  return pars.a_dG1 * exp(-pars.a_dG2 * relative_reserves);
}

// NSC storage capacity: scales with sapwood mass (per Daniel, #517). mass_sapwood
// = area_sapwood(area_leaf) * height * eta_c * rho.
template <typename S>
S TF24_Strategy<S>::storage_capacity(S area_leaf_, S height) const {
  return pars.a_st1 * mass_sapwood(area_sapwood(area_leaf_), height);
}

// Seed the storage state for a newly germinated individual at a_st3 fraction of
// its capacity (Stefaniak et al. 2026, Eq 8), so seedlings are born with
// reserves rather than starting empty (which would kill them immediately).
template <typename S>
void TF24_Strategy<S>::set_initial_states(const TF24_Environment<S>& environment,
                                       Internals<S>& vars) {
  (void)environment;
  // The seed's height is written here, not inherited from the height an
  // Individual was constructed with. Construction runs in plain arithmetic and an
  // active strategy receives its results rather than deriving them, so a height
  // that arrives that way carries no trait derivative and every rate the newborn
  // goes on to have inherits the loss. Declaring it by its own residual at this
  // point puts it where the parameters are already differentiable inputs, and the
  // value written is the same one either way.
  const SeedGeometry seed = seed_geometry();
  vars.set_state(HEIGHT_INDEX, seed.height);
  // A height written into vars leaves the slots it determines holding what the
  // constructor's plain height derived, which is the same number carrying no
  // derivative. Re-deriving them costs one allometry and is the only thing that
  // puts the seed's row into the leaf area every rate at birth size is scaled by.
  update_dependent_aux(HEIGHT_INDEX, vars);
  vars.set_state(state_idx_storage,
                 pars.a_st3 * storage_capacity(seed.area_leaf, seed.height));
}

// [eqn 20] Survival of seedlings during establishment
template <typename S>
S TF24_Strategy<S>::establishment_probability(const TF24_Environment<S>& environment) {
  const SeedGeometry seed = seed_geometry();
  return establishment_probability(
    environment,
    net_mass_production_dt(environment, seed.height, seed.area_leaf,
                           1.0 / seed.height));
}

// Both forms above end here. The carbon is birth-size carbon either way, whatever
// height the caller's plant happens to be at.
template <typename S>
S TF24_Strategy<S>::establishment_probability(const TF24_Environment<S>& environment,
                                               S net_mass_production_dt_) {

  S decay_over_time = exp(-pars.recruitment_decay * environment.time);

  if (net_mass_production_dt_ > 0) {
    const S tmp = pars.a_d0 * seed_geometry().area_leaf / net_mass_production_dt_;
    return 1.0 / (tmp * tmp + 1.0) * decay_over_time;
  } else {
    return 0.0;
  }
}

template <typename S>
S TF24_Strategy<S>::compute_competition(S z, S height) const {
  return pars.k_I * area_leaf(height) * canopy_shape.Q_from_height(z, height);
}

// Ratio-first hot-path overload (see header): receives the cached
// competition_effect (= area_leaf(height)) and height_inverse (= 1/height), so the
// per-call area_leaf() evaluation and z/height division are hoisted out of the
// inner competition loop.
template <typename S>
S TF24_Strategy<S>::compute_competition(S z, S area_leaf_,
                                          S height_inverse) const {
  return pars.k_I * area_leaf_ * canopy_shape.Q(z * height_inverse);
}

// [eqn 10] Cumulative fraction of a quantity distributed over an extent with
//          shape exponent 'eta_x', above coordinate 'z' of a total 'height'.
//          Serves the root mass distribution over soil depth.
template <typename S>
S TF24_Strategy<S>::Q(S z, S height, S eta_x) const {
  if (z > height) {
    return S(0.0);
  }
  // u^eta_x. On double the plain pow; on an active scalar the recorded eta_x
  // derivative u^eta_x * log(u) is 0 * (-inf) -- a NaN -- at u = 0, where the
  // cumulative fraction is 1, so the guard supplies that value outright.
  S u_eta;
  if constexpr (std::is_same_v<S, double>) {
    u_eta = pow(z / height, eta_x);
  } else {
    const S u = z / height;
    u_eta = odelia::util::to_passive(u) <= 0.0 ? S(0.0) : pow(u, eta_x);
  }
  const S tmp = 1.0 - u_eta;
  return tmp * tmp;
}

// The aim is to find a plant height that gives the correct seed mass.
template <typename S>
double TF24_Strategy<S>::height_seed(void) const {

  // Note, these are not entirely correct bounds. Ideally we would use height
  // given *total* mass, not leaf mass, but that is difficult to calculate.
  // Using "height given leaf mass" will expand upper bound, but that's ok
  // most of time. Only issue is that could break with obscure parameter
  // values for LMA or height-leaf area scaling. Could instead use some
  // absolute maximum height for new seedling, e.g. 1m?
  const S
    h0 = height_given_mass_leaf(std::numeric_limits<double>::min()),
    h1 = height_given_mass_leaf(pars.omega);

  const double tol = this->control.offspring_production_tol;
  const size_t max_iterations = this->control.offspring_production_iterations;

  auto target = [&] (double x) mutable -> S {
    return mass_live_given_height(x) - pars.omega;
  };

  if constexpr (std::is_same_v<S, double>) {
    return util::uniroot(target, h0, h1, tol, max_iterations);
  } else {
    // Bisection is affine in its bracket and blind to the residual's values, so
    // recording the search would return d(h1)/d(trait) rather than the height at
    // which mass_live equals the seed mass. Declare it by the residual instead.
    static_assert(std::is_same_v<S, double>,
                  "height_seed() finds its root by iteration; an active scalar must "
                  "declare it through implicit_value on the residual.");
  }
}

template <typename S>
void TF24_Strategy<S>::prepare_strategy() {

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

  canopy_shape.initialise(pars.eta, shading_model_);

  eta_c = CanopyShape<S>::eta_c(pars.eta);
  // NOTE: Also pre-computing, though less trivial
  height_0 = height_seed();
  area_leaf_0 = area_leaf(height_0);

  if (this->is_variable_birth_rate) {
    this->extrinsic_drivers.set_variable("birth_rate", this->birth_rate_x, this->birth_rate_y);
  } else {
    this->extrinsic_drivers.set_constant("birth_rate", this->birth_rate_y[0]);
  }
  if constexpr (std::is_same_v<S, double>) {
    leaf = Leaf(pars.vcmax_25, pars.c, pars.b, pars.psi_crit,
                pars.root_c, pars.root_b, pars.root_psi_crit,
                pars.beta2, pars.jmax_25, pars.a,
                pars.curv_fact_elec_trans, pars.curv_fact_colim,
                this->control.GSS_tol_abs, this->control.vulnerability_curve_ncontrol,
                this->control.ci_abs_tol, this->control.ci_niter,
                pars.g1_TF24);
    // Not a constructor argument, and written before any physiology is set, so
    // the temperature block derives R_d_ from it on the first call.
    leaf.R_d_25 = pars.R_d_25;
  } else {
    static_assert(std::is_same_v<S, double>,
                  "Leaf carries double; an active strategy must supply the leaf's "
                  "local Jacobian across this boundary, not template Leaf.");
  }
}

template <typename S>
typename TF24_Strategy<S>::ptr make_strategy_ptr(TF24_Strategy<S> s) {
  s.prepare_strategy();
  return std::make_shared<TF24_Strategy<S> >(s);
}

}

#endif
