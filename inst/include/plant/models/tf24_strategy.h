// Built from  inst/include/plant/models/ff16_strategy.h on Mon Feb 12 09:52:27 2024 using the scaffolder, from the strategy:  FF16
// -*-c++-*-
#ifndef PLANT_PLANT_TF24_STRATEGY_H_
#define PLANT_PLANT_TF24_STRATEGY_H_

#include <plant/census.h>
#include <plant/strategy.h>
#include <plant/models/tf24_environment.h>
#include <plant/gradient_status.h>
#include <plant/qag.h>
#include <plant/leaf_model.h>
#include <plant/canopy_shape.h>
#include <odelia/ode_util.hpp>
#include <odelia/implicit_node.hpp>
#include <type_traits>

namespace plant {

// What a parameter is to the gradient. One table carries this beside the
// pointer and the name, so a parameter cannot be registered in one list and
// forgotten in another, and what an exact zero means is declared where the
// parameter is rather than recovered from its spelling.
enum class ad_role {
  differentiated,   // a column; an exact zero there is undeclared
  zero_slack,       // a column; an exact zero is complementary slackness
  zero_structural,  // a column; an exact zero is structural, on any trajectory
  refused,          // no column: the model cannot answer for it, and says so by name
  unread            // no column: reaches no equation on this path
};

constexpr bool ad_role_has_column(ad_role role) {
  return role == ad_role::differentiated || role == ad_role::zero_slack ||
         role == ad_role::zero_structural;
}

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
  // aerodynamic resistance ra = C_ra*sqrt(d/U0); inert while PM is off.
  S use_energy_balance = 0.0;
  S d = 0.05;

  struct ad_parameter {
    S* value;
    const char* name;
    ad_role role;
    // The leaf's own index for this parameter, or -1 where the leaf has no input
    // for it. One table, so a trait cannot be registered for the gradient here
    // and forgotten where the leaf's rows are asked for.
    int leaf_par;
  };

  // Every member above, in declaration order, with what each is to the
  // gradient. field_ptrs() and the ad_parameter_* projections all read this, so
  // a member reaches all of them or none.
#define PLANT_TF24_AD_PARAMETER(x, r) ad_parameter{&x, #x, ad_role::r, -1}
  // A parameter the leaf answers for. `p` is the leaf's spelling, which for four
  // of the fourteen is not this file's -- and writing it here is what makes a
  // rename on either side a compile error rather than a silently shifted column.
#define PLANT_TF24_LEAF_PARAMETER(x, r, p) \
  ad_parameter { &x, #x, ad_role::r, phylloptim::gradient::par_##p }
  std::vector<ad_parameter> ad_parameter_table() {
    return {
      PLANT_TF24_AD_PARAMETER(lma, differentiated),
      PLANT_TF24_AD_PARAMETER(rho, differentiated),
      PLANT_TF24_AD_PARAMETER(hmat, differentiated),
      PLANT_TF24_AD_PARAMETER(omega, differentiated),
      // The canopy profile's u^eta reaches u = 0, where the recorded derivative
      // u^eta*log(u) is 0*(-inf) and the guard returns a constant zero instead,
      // so a recorded row would be a silently wrong zero rather than a NaN.
      PLANT_TF24_AD_PARAMETER(eta, refused),
      PLANT_TF24_AD_PARAMETER(theta, differentiated),
      PLANT_TF24_AD_PARAMETER(a_l1, differentiated),
      PLANT_TF24_AD_PARAMETER(a_l2, differentiated),
      PLANT_TF24_AD_PARAMETER(a_r1, differentiated),
      PLANT_TF24_AD_PARAMETER(a_b1, differentiated),
      PLANT_TF24_AD_PARAMETER(r_s, differentiated),
      PLANT_TF24_AD_PARAMETER(r_b, differentiated),
      PLANT_TF24_AD_PARAMETER(r_r, differentiated),
      PLANT_TF24_AD_PARAMETER(r_l, differentiated),
      PLANT_TF24_AD_PARAMETER(a_y, differentiated),
      PLANT_TF24_AD_PARAMETER(a_bio, differentiated),
      PLANT_TF24_AD_PARAMETER(k_l, differentiated),
      PLANT_TF24_AD_PARAMETER(k_b, differentiated),
      PLANT_TF24_AD_PARAMETER(k_s, differentiated),
      PLANT_TF24_AD_PARAMETER(k_r, differentiated),
      // The light-response curve the Farquhar leaf replaced.
      PLANT_TF24_AD_PARAMETER(a_p1, unread),
      PLANT_TF24_AD_PARAMETER(a_p2, unread),
      // Occurs only in fecundity_dt's denominator, and no census metric reads
      // fecundity: zero on any trajectory rather than on this one.
      PLANT_TF24_AD_PARAMETER(a_f3, zero_structural),
      PLANT_TF24_AD_PARAMETER(a_f1, differentiated),
      PLANT_TF24_AD_PARAMETER(a_f2, differentiated),
      PLANT_TF24_AD_PARAMETER(S_D, unread),
      PLANT_TF24_AD_PARAMETER(a_d0, differentiated),
      PLANT_TF24_AD_PARAMETER(d_I, differentiated),
      PLANT_TF24_AD_PARAMETER(a_dG1, differentiated),
      PLANT_TF24_AD_PARAMETER(a_dG2, differentiated),
      PLANT_TF24_AD_PARAMETER(a_st1, differentiated),
      PLANT_TF24_AD_PARAMETER(a_st2, differentiated),
      PLANT_TF24_AD_PARAMETER(a_st3, differentiated),
      PLANT_TF24_AD_PARAMETER(k_I, differentiated),
      PLANT_TF24_LEAF_PARAMETER(vcmax_25, differentiated, vcmax_25),
      // Read only by c's and b's default initialisers, so a value set after
      // construction reaches nothing here; d/dp_50 belongs to whatever computes
      // c and b, which for a run driven from traits is the R hyperparameters.
      PLANT_TF24_AD_PARAMETER(p_50, unread),
      PLANT_TF24_AD_PARAMETER(K_s, differentiated),
      PLANT_TF24_LEAF_PARAMETER(c, differentiated, stem_c),
      PLANT_TF24_LEAF_PARAMETER(b, differentiated, stem_b),
      // Sets the dry bound of an interval the operating point is strictly
      // inside, so slackness makes its row zero at an interior optimum. Live at
      // a pin.
      PLANT_TF24_LEAF_PARAMETER(psi_crit, zero_slack, psi_crit),
      PLANT_TF24_AD_PARAMETER(beta1, unread),
      PLANT_TF24_LEAF_PARAMETER(beta2, differentiated, beta2),
      PLANT_TF24_LEAF_PARAMETER(g1_TF24, differentiated, cost_scale_TF24),
      PLANT_TF24_LEAF_PARAMETER(jmax_25, differentiated, jmax_25),
      PLANT_TF24_LEAF_PARAMETER(a, differentiated, a),
      PLANT_TF24_LEAF_PARAMETER(curv_fact_elec_trans, differentiated, curv_fact_elec_trans),
      PLANT_TF24_LEAF_PARAMETER(curv_fact_colim, differentiated, curv_fact_colim),
      PLANT_TF24_LEAF_PARAMETER(R_d_25, differentiated, R_d_25),
      // Declared and carried; no equation on this path reads them.
      PLANT_TF24_AD_PARAMETER(var_sapwood_volume_cost, unread),
      PLANT_TF24_AD_PARAMETER(nmass_l, unread),
      PLANT_TF24_AD_PARAMETER(nmass_s, unread),
      PLANT_TF24_AD_PARAMETER(nmass_b, unread),
      PLANT_TF24_AD_PARAMETER(nmass_r, unread),
      PLANT_TF24_AD_PARAMETER(dmass_dN, unread),
      // Q()'s u^eta_x reaches u = 0 with the same 0*(-inf) derivative as eta,
      // and the guard supplies the cumulative fraction there outright.
      PLANT_TF24_AD_PARAMETER(root_depth_shape_eta, refused),
      PLANT_TF24_LEAF_PARAMETER(root_c, differentiated, root_c),
      PLANT_TF24_LEAF_PARAMETER(root_b, differentiated, root_b),
      // The root curve's dry bound, slack for the same reason as psi_crit.
      PLANT_TF24_LEAF_PARAMETER(root_psi_crit, zero_slack, root_psi_crit),
      PLANT_TF24_AD_PARAMETER(rooting_depth_max, differentiated),
      PLANT_TF24_AD_PARAMETER(recruitment_decay, differentiated),
      // The gate is compared rather than differentiated.
      PLANT_TF24_AD_PARAMETER(use_energy_balance, unread),
      // No row in the leaf's supplied Jacobian.
      PLANT_TF24_AD_PARAMETER(d, unread)
    };
  }
#undef PLANT_TF24_AD_PARAMETER

  // Every entry of the table, which is the whole parameter set: a rebind carries
  // it across a scalar change. ad_parameters() is the subset with a column.
  std::vector<S*> field_ptrs() {
    const std::vector<ad_parameter> table = ad_parameter_table();
    std::vector<S*> ret;
    ret.reserve(table.size());
    for (const ad_parameter& p : table) {
      ret.push_back(p.value);
    }
    return ret;
  }
  static constexpr size_t field_count = 62;
};

// Every member of TF24_Pars is an S, so one added without extending
// ad_parameter_table() changes this size and is refused rather than dropped.
static_assert(sizeof(TF24_Pars<double>) ==
              TF24_Pars<double>::field_count * sizeof(double),
              "TF24_Pars has a member ad_parameter_table() does not list");

// Templated on the scalar S the state, the traits and everything derived from
// them carry; double is production. The embedded Leaf, the Control tolerances
// and the extrinsic drivers stay double.
// The operating points a run's leaf solves found, held against the rate evaluation
// that found them and in the order it found them -- which is the order the same
// evaluation finds them in again, because it runs at the same state.
//
// A rejected step's attempt writes the same slot as the retry that replaces it, so
// nothing here has to commit a step of its own.
class leaf_solved_points {
public:
  // The evaluation about to run. `keeping` is the run that fills the record; a pass
  // that is not it reads the slot instead, and finds nothing where the run addressed
  // no evaluation there.
  void begin_stage(odelia::ode::recorded_stage at, bool keeping) {
    solved = 0;
    slot = nullptr;
    filling = keeping;
    const size_t stage = static_cast<size_t>(at.stage);
    if (keeping) {
      if (kept.size() <= at.step) {
        kept.resize(at.step + 1);
      }
      if (kept[at.step].size() <= stage) {
        kept[at.step].resize(stage + 1);
      }
      slot = &kept[at.step][stage];
      slot->clear();
    } else if (at.step < kept.size() && stage < kept[at.step].size()) {
      slot = &kept[at.step][stage];
    }
  }

  // One evaluation is one address, so this ends where the rates are read. A solve
  // outside a step then keeps nothing and places nothing.
  void end_stage() {
    solved = 0;
    slot = nullptr;
  }

  // The next point in the open slot, or one holding nothing.
  Leaf::SolvedPoint next() {
    if (filling || slot == nullptr || solved >= slot->size()) {
      return Leaf::SolvedPoint();
    }
    ++placed;
    return (*slot)[solved++];
  }

  // How many points this record placed. Counted because a record that engages and
  // one that quietly does not are the same green suite: every number a placement
  // produces is the number a search produces, so nothing else can tell them apart.
  size_t placements() const { return placed; }

  void keep(const Leaf::SolvedPoint& point) {
    if (filling && slot != nullptr) {
      slot->push_back(point);
    }
  }

  void clear() {
    kept.clear();
    end_stage();
  }

private:
  std::vector<std::vector<std::vector<Leaf::SolvedPoint>>> kept;
  std::vector<Leaf::SolvedPoint>* slot = nullptr;
  bool filling = false;
  size_t solved = 0;
  size_t placed = 0;
};

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

  S compute_average_light_environment(const S& z, const S& height,
                                      const TF24_Environment<S> &environment);

  // calculate the amount of water transpired relativised by leaf area index.

  S evapotranspiration_dt(const S& area_leaf_, int soil_layer);


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

  // The metrics a census of this model sums, in the order it reports them. Beside
  // the state and aux names because all three say what this model calls its own
  // quantities -- and every slot is read by its cached index, which is how every
  // other reader of these slots reaches them.
  static std::vector<census_metric<TF24_Strategy<S>>> census_metrics() {
    using strategy = TF24_Strategy<S>;
    return {
      {"leaf_area",
       [](const strategy& p, const Internals<S>& vars) -> S {
         return p.area_leaf(vars.state(HEIGHT_INDEX));
       }},
      {"mass_above_ground",
       [](const strategy& p, const Internals<S>& vars) -> S {
         const S height = vars.state(HEIGHT_INDEX);
         const S area_leaf = p.area_leaf(height);
         return p.mass_above_ground(
             p.mass_leaf(area_leaf),
             p.mass_bark(p.area_bark(area_leaf), height),
             p.mass_sapwood(p.area_sapwood(area_leaf), height),
             vars.state(p.state_idx_mass_heartwood));
       }},
      {"area_stem",
       [](const strategy& p, const Internals<S>& vars) -> S {
         const S height = vars.state(HEIGHT_INDEX);
         const S area_leaf = p.area_leaf(height);
         return p.area_stem(p.area_bark(area_leaf), p.area_sapwood(area_leaf),
                            vars.state(p.state_idx_area_heartwood));
       }},
    };
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

  using ad_parameter = typename TF24_Pars<S>::ad_parameter;

  // Addresses of the parameters a gradient can be taken with respect to, in the
  // order ad_parameter_names() gives: the table entries whose role has a column.
  // Both allocate, so take them once per gradient evaluation and hold them for
  // the run rather than per block; the strategy is shared and the fields do not
  // move. Index against .size().
  std::vector<S*> ad_parameters() {
    const std::vector<ad_parameter> table = pars.ad_parameter_table();
    std::vector<S*> ret;
    ret.reserve(table.size());
    for (const ad_parameter& p : table) {
      if (ad_role_has_column(p.role)) {
        ret.push_back(p.value);
      }
    }
    return ret;
  }

  // The name each gradient column carries, one per ad_parameters() entry and in
  // that order. Which parameters are left out, and why, is on the table entries.
  std::vector<std::string> ad_parameter_names() {
    const std::vector<ad_parameter> table = pars.ad_parameter_table();
    std::vector<std::string> ret;
    ret.reserve(table.size());
    for (const ad_parameter& p : table) {
      if (ad_role_has_column(p.role)) {
        ret.push_back(p.name);
      }
    }
    return ret;
  }

  // What an exactly-zero entry in a column MEANS, one per ad_parameters() entry.
  // Consulted only where the number is exactly zero, which is what keeps it
  // honest as the census grows -- a parameter that becomes live stops being zero
  // and its role is never read.
  std::vector<gradient_status::Kind> ad_parameter_zero_classes() {
    const std::vector<ad_parameter> table = pars.ad_parameter_table();
    std::vector<gradient_status::Kind> ret;
    ret.reserve(table.size());
    for (const ad_parameter& p : table) {
      switch (p.role) {
      case ad_role::differentiated:
        ret.push_back(gradient_status::Kind::zero_undeclared);
        break;
      case ad_role::zero_slack:
        ret.push_back(gradient_status::Kind::zero_slack);
        break;
      case ad_role::zero_structural:
        ret.push_back(gradient_status::Kind::zero_structural);
        break;
      case ad_role::refused:
      case ad_role::unread:
        break;
      }
    }
    return ret;
  }

  // Translate generic methods to TF24 strategy leaf area methods

  S competition_effect(const S& height) const {
    return area_leaf(height);
  }

  void refresh_indices();


  // TF24 Methods  ----------------------------------------------

  // [eqn 2] area_leaf (inverse of [eqn 3])
  S area_leaf(const S& height) const;

  // [eqn 1] mass_leaf (inverse of [eqn 2])
  S mass_leaf(const S& area_leaf) const;

  // [eqn 4] area and mass of sapwood
  S area_sapwood(const S& area_leaf) const;
  S mass_sapwood(const S& area_sapwood, const S& height) const;

  // [eqn 5] area and mass of bark
  S area_bark(const S& area_leaf) const;
  S mass_bark (const S& area_bark, const S& height) const;

  S area_stem(const S& area_bark, const S& area_sapwood,
                            const S& area_heartwood) const;
  S diameter_stem(const S& area_stem) const;

  // [eqn 7] Mass of (fine) roots
  S mass_root(const S& area_leaf) const;

  // [eqn 8] Total Mass
  S mass_live(const S& mass_leaf, const S& mass_bark,
              const S& mass_sapwood, const S& mass_root) const;

  S mass_total(const S& mass_leaf, const S& mass_bark, const S& mass_sapwood,
               const S& mass_heartwood, const S& mass_root) const;

  // Above-ground mass = leaf + all stem components (bark + sapwood +
  // heartwood); excludes roots.
  S mass_above_ground(const S& mass_leaf, const S& mass_bark,
                      const S& mass_sapwood, const S& mass_heartwood) const;

  void compute_rates(const TF24_Environment<S>& environment,
                Internals<S>& vars);
  
  void compute_roots(const TF24_Environment<S>& environment,
                Internals<S>& vars);

  void update_dependent_aux(const int index, Internals<S>& vars);

  // * Mass production
  // [eqn 12] Gross annual CO2 assimilation
  S assimilation(const TF24_Environment<S>& environment, const S& height,
                 const S& area_leaf);
  // [Appendix S6] Per-leaf photosynthetic rate.
  S assimilation_leaf(const S& x) const;

  // [eqn 13] Total maintenance respiration
  S respiration(const S& mass_leaf, const S& mass_sapwood,
                const S& mass_bark, const S& mass_root) const;

  S respiration_leaf(const S& mass) const;
  S respiration_bark(const S& mass) const;
  S respiration_sapwood(const S& mass) const;
  S respiration_root(const S& mass) const;

  // [eqn 14] Total turnover
  S turnover(const S& mass_leaf, const S& mass_bark,
             const S& mass_sapwood, const S& mass_root) const;
  S turnover_leaf(const S& mass) const;
  S turnover_bark(const S& mass) const;
  S turnover_sapwood(const S& mass) const;
  S turnover_root(const S& mass) const;

  // [eqn 15] Net production
  S net_mass_production_dt_A(const S& assimilation, const S& respiration,
                             const S& turnover) const;

  virtual S net_mass_production_dt(const TF24_Environment<S>& environment,
                                const S& height, const S& area_leaf_,
                                const S& height_inverse);

  // Resolve the leaf operating point on the already-set-up `leaf` (i.e. after
  // leaf.set_physiology(...)). Base TF24 optimises the root-collar psi via
  // golden-section search; the TF24f variant overrides this to make the optimum
  // chase a tracked ODE state (#525). Called per crown light point from
  // net_mass_production_dt, so it must be virtual to dispatch to the override
  // when net_mass_production_dt is reused unchanged by the subclass.
  virtual void solve_leaf();

  // Read how the leaf's outputs respond to what it was given, and record each
  // output that re-enters the active chain carrying that response. The leaf
  // drives and solves itself from `leaf_drivers_`; nothing here re-supplies it.
  void record_leaf_outputs(const S& radiation, const std::vector<S>& psi_soil,
                           const S& conductance_max);
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
  virtual S fraction_allocation_reproduction(const S& height) const;
  S fraction_allocation_growth(const S& height) const;
  // [eqn 17] Rate of offspring production
  S fecundity_dt(const S& net_mass_production_dt,
                 const S& fraction_allocation_reproduction) const;

  // [eqn 18] Fraction of mass growth that is leaves
  S darea_leaf_dmass_live(const S& area_leaf) const;

  // change in height per change in leaf area
  S dheight_darea_leaf(const S& area_leaf) const;
  // Mass of leaf needed for new unit area leaf, d m_s / d a_l
  S dmass_leaf_darea_leaf(const S& area_leaf) const;
  // Mass of stem needed for new unit area leaf, d m_s / d a_l
  S dmass_sapwood_darea_leaf(const S& area_leaf) const;
  // Mass of bark needed for new unit area leaf, d m_b / d a_l
  S dmass_bark_darea_leaf(const S& area_leaf) const;
  // Mass of root needed for new unit area leaf, d m_r / d a_l
  S dmass_root_darea_leaf(const S& area_leaf) const;
  // Growth rate of basal diameter_stem per unit stem area
  S ddiameter_stem_darea_stem(const S& area_stem) const;
  // Growth rate of components per unit time:
  S area_leaf_dt(const S& area_leaf_dt) const;
  S area_sapwood_dt(const S& area_leaf_dt) const;
  S area_heartwood_dt(const S& area_leaf) const;
  S area_bark_dt(const S& area_leaf_dt) const;
  S area_stem_dt(const S& area_leaf, const S& area_leaf_dt) const;
  S diameter_stem_dt(const S& area_stem, const S& area_stem_dt) const;
  S mass_root_dt(const S& area_leaf,
                 const S& area_leaf_dt) const;
  S mass_live_dt(const S& fraction_allocation_reproduction,
                 const S& net_mass_production_dt) const;
  S mass_total_dt(const S& fraction_allocation_reproduction,
                  const S& net_mass_production_dt,
                  const S& mass_heartwood_dt) const;
  S mass_above_ground_dt(const S& area_leaf,
                         const S& fraction_allocation_reproduction,
                         const S& net_mass_production_dt,
                         const S& mass_heartwood_dt,
                         const S& area_leaf_dt) const;

  S mass_heartwood_dt(const S& mass_sapwood) const;

  S mass_live_given_height(const S& height) const;
  S height_given_mass_leaf(const S& mass_leaf_) const;


  S mortality_dt(const S& relative_reserves, const S& cumulative_mortality) const;
  S mortality_growth_independent_dt()const ;
  // Storage-dependent growth mortality (#517): rises smoothly as relative
  // reserves r = S/S_max deplete, bounded in [a_dG1*e^-a_dG2, a_dG1].
  S mortality_storage_dependent_dt(const S& relative_reserves) const;
  // NSC storage capacity S_max = a_st1 * mass_sapwood [kg NSC].
  S storage_capacity(const S& area_leaf, const S& height) const;
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
                              const S& net_mass_production_dt_);

  // * Competitive environment
  // [eqn 11] total projected leaf area above height above height `z` for given plant
  S compute_competition(const S& z, const S& height) const;
  // Optimised overload called from Individual<TF24>::compute_competition with the
  // cached competition_effect (= area_leaf(height)) and height_inverse (= 1/height)
  // aux values, matching the shared individual.h interface (no recompute per call).
  S compute_competition(const S& z, const S& area_leaf_,
                        const S& height_inverse) const;
  // Strategy-agnostic entry point used by Individual<TF24> (#266): reads the
  // cached competition_effect and height_inverse aux slots itself.
  S compute_competition(const S& z, const Internals<S>& vars) const {
    return compute_competition(z, vars.aux(aux_idx_competition_effect),
                               vars.aux(aux_idx_height_inverse));
  }

  // The competition contribution and its vertical derivative from one pass, so
  // u^eta is evaluated once. The first entry is bit-for-bit the one
  // compute_competition() returns, and both read the shading model's own profile:
  // the value used to read the smooth one directly, so under a flat-top profile
  // the two disagreed while a comment said they could not.
  std::pair<S, S> compute_competition_and_slope(const S& z, const Internals<S>& vars) const {
    const S& area_leaf_ = vars.aux(aux_idx_competition_effect);
    const S height_inverse = vars.aux(aux_idx_height_inverse);
    const S scale = pars.k_I * area_leaf_;
    const std::pair<S, S> Qq =
      canopy_shape.Q_and_q(z * height_inverse, z, height_inverse);
    return {scale * Qq.first, -(scale * Qq.second)};
  }


  // The fraction of root mass below soil depth `z`, for a plant rooted to
  // `rooting_depth` with shape exponent `eta_x` (pars.root_depth_shape_eta). The
  // canopy's own cumulative form is CanopyShape::Q, at pars.eta.
  S Q(const S& z, const S& rooting_depth, const S& eta_x) const;

  // The inverse of dheight_darea_leaf, so the allometry has one source.
  S darea_leaf_dheight(const S& area_leaf) const {
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
        [&](const S& y) -> S { return mass_live_given_height(y) - pars.omega; });
      return {h, area_leaf(h)};
    }
  }

  // Set constants within TF24_Strategy
  void prepare_strategy();

  // The same strategy at scalar U.

  // Another strategy's values, written into this one. prepare_strategy() is
  // refused at an active scalar, so what it produced is carried rather than
  // rebuilt: the leaf, and the canopy shape below.
  //
  // A rebind hands back a fresh strategy, so anything it leaves default has to be
  // set here rather than left alone -- an assignment writes into one that already
  // exists, and a member it does not write keeps a previous state's working.
  template <class S1>
  void assign_from(const TF24_Strategy<S1>& src) {
    // Qualified: these are the base's, and an unqualified name is not looked up
    // in a dependent base.
    this->birth_rate_x = src.birth_rate_x;
    this->birth_rate_y = src.birth_rate_y;
    this->is_variable_birth_rate = src.is_variable_birth_rate;
    this->collect_all_auxiliary = src.collect_all_auxiliary;
    this->size_0 = src.size_0;
    this->control = src.control;
    this->name = src.name;
    this->extrinsic_drivers = src.extrinsic_drivers;

    TF24_Pars<S1> from_pars = src.pars;
    std::vector<S1*> from = from_pars.field_ptrs();
    std::vector<S*> to = pars.field_ptrs();
    // field_count is a literal, so it is what a new member gets added without;
    // this is what says the two lists still describe the same struct.
    util::check_length(from.size(), to.size());
    util::check_length(to.size(), TF24_Pars<S>::field_count);
    for (size_t i = 0; i < from.size(); ++i) {
      *to[i] = S(odelia::util::to_passive(*from[i]));
    }

    shading_model_ = src.shading_model_;
    eta_c = S(odelia::util::to_passive(src.eta_c));
    canopy_shape.initialise(pars.eta, shading_model_);
    height_0 = src.height_0;
    area_leaf_0 = S(odelia::util::to_passive(src.area_leaf_0));
    leaf = src.leaf;
    storage_gate_width = src.storage_gate_width;
    storage_prod_eps = src.storage_prod_eps;
    beta_R_H = src.beta_R_H;
    beta_R_V = src.beta_R_V;
    function_integrator = src.function_integrator;
    // Shared, not copied: a clamp on this path fires inside a per-unit copy that
    // is then discarded, so the count has to land in storage the run still owns.
    // The missing-row flag is shared for the same reason.
    clamps.differentiated = src.clamps.differentiated;
    curvature_margin = src.curvature_margin;
    uptake_rows_unavailable = src.uptake_rows_unavailable;
    uptake_rows_reason = src.uptake_rows_reason;
    leaf_points = src.leaf_points;

    // Sized, not copied: these hold one right-hand side's working, and are
    // rewritten before they are read.
    root_carbon_per_leaf_area_.assign(src.root_carbon_per_leaf_area_.size(),
                                      S(0.0));

    // The index maps and the slot numbers are a function of the names, so they
    // are derived rather than carried.
    refresh_indices();
  }

  template <class U>
  TF24_Strategy<U> rebind_from() const {
    TF24_Strategy<U> out;
    out.assign_from(*this);
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

  // The clamp sites this strategy reaches; the list itself is shared with the
  // environment, which reaches the others.
  clamp_counter clamps;

  // Count a clamp against the path it fired on. Which path this is is a property
  // of the scalar, so it is decided at compile time.
  void note_clamp(int site) {
    if constexpr (std::is_same_v<S, double>) {
      ++clamps.forward[site];
    } else {
      ++(*clamps.differentiated)[site];
    }
  }

  // How close the interior derivation came to dividing by nothing. The guard
  // below refuses on a floor, and a floor nothing approaches and a floor nothing
  // reaches report the same green -- so the margin is carried out beside the
  // count rather than left to be assumed. Shared for the same reason the counts
  // are: the copy that measures it is discarded.
  std::shared_ptr<double> curvature_margin = std::make_shared<double>(-1.0);
  void note_curvature(double value) {
    const double m = std::abs(value);
    if (*curvature_margin < 0.0 || m < *curvature_margin) {
      *curvature_margin = m;
    }
  }

  // Below this the collar's own response is amplification rather than an answer.
  // It lives on Control because it moves which states answer, so two gradients
  // taken at different values are gradients of different functions and
  // stand_gradient() has to be able to say so.
  //
  // The units are the profit's, so a reparameterisation that rescales profit
  // needs it re-measured -- which is what curvature_margin is for.
  double curvature_floor() const { return this->control.gradient_curvature_floor; }

  // Set where a row the UPTAKE outputs need does not exist. The profit row is
  // built by the envelope theorem and reads no curvature, so it survives every
  // degeneracy the water rows do not -- and the values still do too, so the
  // recording carries them on and leaves the missing rows off the tape.
  //
  // Recorded rather than thrown because a throw from here is a refusal of the
  // whole recording, and this one is a refusal of two of its outputs. The sweep
  // reads it once it is over and turns it into the run's status.
  //
  // Shared for the reason the clamp counts are: the strategy that measures it is
  // a per-unit copy the sweep discards, so a flag on the copy is a flag nothing
  // can read. The first reason is kept and later ones dropped, and the sweep
  // clears it before it starts -- so this latches for one gradient call rather
  // than for one operating point.
  std::shared_ptr<bool> uptake_rows_unavailable = std::make_shared<bool>(false);
  std::shared_ptr<std::string> uptake_rows_reason =
      std::make_shared<std::string>();

  // How many leaf solves landed in each operating-point kind, indexed by the
  // enum. Diagnostic and deliberately outside rebind_from: a block copies the
  // strategy per unit and discards it, so carrying the tally across would count
  // the sweep's copies as well as the run.
  std::vector<size_t> operating_point_counts =
    std::vector<size_t>(Leaf::operating_point_kind_count, 0);

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

  // Everything the leaf is driven with, in one object, refilled rather than
  // rebuilt each call: the network is filled in place by the architecture model
  // and must NOT be moved from, because the leaf takes it by const reference
  // precisely so this buffer keeps its capacity across calls. All double -- the
  // leaf is a node this strategy supplies rows for, not a scalar it lifts.
  phylloptim::gradient::Drivers leaf_drivers_;

  // The input step is a relative 1e-3, which is what it takes to move total
  // uptake above the solve's own floor; the collar's stays small, because a step
  // crossing a narrow bracket is what detects a pinned point.
  static phylloptim::gradient::Settings leaf_row_settings() {
    phylloptim::gradient::Settings s;
    s.step = 1e-3;
    s.collar = 1e-6;
    return s;
  }

  // Set where a row the UPTAKE outputs need does not exist. The first reason is
  // kept: a later one is a consequence of the same degeneracy.
  void lose_uptake_rows(const std::string& why) {
    if (!*uptake_rows_unavailable) {
      *uptake_rows_unavailable = true;
      *uptake_rows_reason = why;
    }
  }

  // What the leaf clamped between two readings of its own tally. Everything
  // between re-supplies and re-solves it many times, so anything clamped there
  // is the gradient's and the forward run's share is the total less this.
  void note_leaf_clamps(const std::vector<std::size_t>& before) {
    const std::vector<std::size_t> after = leaf.clamp_counts();
    for (std::size_t s = 0; s < after.size(); ++s) {
      const std::size_t at = CLAMP_LEAF_FIRST + s;
      if (at < CLAMP_SITE_COUNT && after[s] > before[s]) {
        (*clamps.differentiated)[at] += after[s] - before[s];
      }
    }
  }

  // What this species' leaf solves found, kept against the rate evaluation that
  // found them so a pass re-running the model over the same states places the point
  // instead of searching for it again. The patch hands the address down.
  //
  // Shared, not copied, for the reason the clamp counts are: the pass that reads the
  // record runs on a rebound strategy.
  std::shared_ptr<leaf_solved_points> leaf_points =
      std::make_shared<leaf_solved_points>();
  void begin_stage(odelia::ode::recorded_stage at, bool keeping) {
    leaf_points->begin_stage(at, keeping);
  }
  void end_stage() { leaf_points->end_stage(); }
  size_t leaf_placements() const { return leaf_points->placements(); }

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
    const S& z, const S& height, const TF24_Environment<S> &environment) {
// NOTE: the light environment is clamped to a small positive floor (1e-4)
// rather than allowed to reach 0 (original rationale was never recorded;
// preserved as-is).
//
// Counted, because this is where the floor binds first: the mean this feeds is
// an integral of already-floored values against a shape that integrates to one,
// so a floored point here cannot lower the mean below the floor and the read
// downstream sees nothing.
     using std::max;
     const S light = environment.get_environment_at_height(z);
     if (light < S(0.0001)) {
       note_clamp(CLAMP_LIGHT_FLOOR_CROWN);
     }
     return max(light, S(0.0001)) * canopy_shape.q_from_height(z, height);
}

// assumes optimise_psi_stem_TF has been run for optimal psi_stem
template <typename S>
S TF24_Strategy<S>::evapotranspiration_dt(const S& area_leaf_, int soil_layer) {
  if constexpr (std::is_same_v<S, double>) {
    return leaf.soil_consumption_[soil_layer] * area_leaf_;
  } else {
    return leaf_soil_consumption_[soil_layer] * area_leaf_;
  }
}

// Two of the leaf's outputs re-enter the active chain, and they reach the
// operating point differently.
//
// PROFIT is the objective at the point the leaf chose, so by the envelope
// theorem its response at an interior optimum is the direct one and the point's
// own movement drops out. PER-LAYER UPTAKE is set as a side effect AT that
// point, so it consumes the choice rather than being it and the movement is part
// of its answer.
//
// The leaf returns both in PARTS -- the gradient of whatever condition defines
// the point, each output's response to the point, and the rest at a held one --
// so the point is recorded once here and every output hangs off it. The tape
// then carries one term per input plus one per output rather than their product,
// and the roles are already inside the numbers: nothing below asks which kind of
// point this is.
template <typename S>
void TF24_Strategy<S>::record_leaf_outputs(const S& radiation,
                                           const std::vector<S>& psi_soil,
                                           const S& conductance_max) {
  namespace grad = phylloptim::gradient;
  using odelia::util::to_passive;
  const int n_layer = static_cast<int>(psi_soil.size());

  // Every input the leaf answers for, paired with the active value its row is
  // recorded against. A layer holding no root carbon is left out: the
  // architecture model gave the network no slot for it, so the leaf has no row
  // there and asking would only carry a refusal further.
  //
  // theta[par_resistance] is neither written nor read -- it belongs to the
  // single-potential supply path and this leaf is on the multi-layer one.
  double theta[grad::n_pars] = {};
  const std::size_t room =
      static_cast<std::size_t>(grad::n_pars_total(n_layer));
  std::vector<int> input;
  std::vector<const S*> at;
  input.reserve(room);
  at.reserve(room);
  for (const auto& p : pars.ad_parameter_table()) {
    if (p.leaf_par < 0) {
      continue;
    }
    theta[p.leaf_par] = to_passive(*p.value);
    input.push_back(p.leaf_par);
    at.push_back(p.value);
  }
  theta[grad::par_kmax] = to_passive(conductance_max);
  input.push_back(grad::par_kmax);
  at.push_back(&conductance_max);
  input.push_back(grad::par_PPFD);
  at.push_back(&radiation);
  for (int j = 0; j < n_layer; ++j) {
    input.push_back(grad::par_psi_soil_first + j);
    at.push_back(&psi_soil[static_cast<std::size_t>(j)]);
  }
  for (int k = 0; k < n_layer; ++k) {
    if (root_carbon_value_[static_cast<std::size_t>(k)] > 0.0) {
      input.push_back(grad::par_root_carbon_first(n_layer) + k);
      at.push_back(&root_carbon_per_leaf_area_[static_cast<std::size_t>(k)]);
    }
  }
  const std::size_t n_input = input.size();

  // Which of these is the objective is the leaf's to say, not this caller's --
  // `role_of` reads it off the index. So this is a list of outputs and nothing
  // else, and there is no second list beside it to disagree with.
  std::vector<int> output;
  output.reserve(static_cast<std::size_t>(n_layer) + 1);
  output.push_back(grad::out_profit);
  for (int i = 0; i < n_layer; ++i) {
    output.push_back(grad::out_uptake_first + i);
  }

  // The leaf clamps in double on both paths, so which path a clamp fired on is a
  // question of when rather than of the scalar, and a delta answers it. Taken on
  // the throwing exit too, because a refusal is still a pass over the leaf.
  const std::vector<std::size_t> clamps_before = leaf.clamp_counts();
  const grad::RowRequest request{output.data(), output.size(), input.data(),
                                 n_input};
  grad::Rows rows;
  try {
    // `rows_at` is the read and takes the leaf and the request alone; this is the
    // entry that also differences whatever the read refuses, which is one branch --
    // a shade death whose wet bound has no derivative.
    rows = grad::rows_differenced(leaf, theta, leaf_drivers_, request,
                                  leaf_row_settings());
  } catch (const std::runtime_error& e) {
    note_leaf_clamps(clamps_before);
    throw gradient_refusal(std::string("TF24 gradient: ") + e.what());
  }
  note_leaf_clamps(clamps_before);

  // The interior derivation divides by the profit's curvature. At a bound the
  // slope is the bound's own and this floor is not about it.
  if (rows.kind == phylloptim::Leaf::OperatingPointKind::Interior) {
    note_curvature(rows.residual_slope);
    if (!(rows.residual_slope < 0.0) ||
        std::abs(rows.residual_slope) < curvature_floor()) {
      lose_uptake_rows(
          "TF24 gradient: the leaf's profit curvature at this operating point "
          "is " + util::to_string(rows.residual_slope) + ", against a floor of " +
          util::to_string(curvature_floor()) + ", so the collar's own response "
          "is amplification rather than an answer");
    }
  }

  // What a missing row costs depends on which output loses it: profit is what net
  // production reads and has no partial form, so losing its row ends the
  // gradient, where a water row can go missing and still leave the value the
  // patch balance needs. Which output is which is the leaf's to say rather than
  // this caller's -- the objective is the one that does not read the operating
  // point, so it is the one a lost point costs nothing.
  //
  // ⚠️ AND WHATEVER IS MISSING IS NAMED HERE. The row layer names an input whose
  // row it refused, but a number can also arrive non-finite from a channel that
  // reported no refusal, and the graft is what sees that. Reporting
  // `rows.message` alone lets it through as a refusal with an empty reason, which
  // is the one shape a refusal must never take: it tells a reader the gradient
  // stopped and nothing about where.
  const int n_layer_named = n_layer;
  auto missing = [&](std::size_t j, const odelia::graft_report& report,
                     bool carries_point) -> std::string {
    std::string out = std::string("TF24 gradient: `") +
                      grad::output_name(output[j], n_layer_named) + "`'s ";
    if (carries_point && report.at == 0) {
      out += "channel into the operating point";
    } else {
      const std::size_t i = report.at - (carries_point ? 1 : 0);
      out += "row in `" +
             grad::par_name(input[i], n_layer_named) + "`";
    }
    out += " cannot be recorded at a point the solve called " +
           std::string(Leaf::operating_point_kind_name(rows.kind)) + ": " +
           report.why;
    if (!rows.message.empty()) {
      out += "; " + rows.message;
    }
    return out;
  };

  leaf_soil_consumption_.assign(psi_soil.size(), S(0.0));
  std::vector<odelia::input_and_derivative<S>> against;
  against.reserve(n_input + 1);

  // Formed once however many outputs read it, which is the whole economy of an
  // implicit node: the theorem's quotient, and the reading where its slope cannot
  // be inverted, are taken once rather than per output. An output needs the point
  // only where it reads it, so the objective at an interior optimum is whole
  // without one -- which is the envelope theorem, and is why it survives the
  // degeneracy that costs the water rows theirs.
  S point;
  bool reads_point = false;
  for (std::size_t j = 0; j < output.size(); ++j) {
    reads_point = reads_point || rows.dy_dp[j] != 0.0;
  }
  bool point_whole = true;
  std::string point_why;
  if (reads_point) {
    for (std::size_t i = 0; i < n_input; ++i) {
      against.push_back({*at[i], rows.dresidual[i]});
    }
    const odelia::graft_report report =
        odelia::implicit_root<S>(rows.point, rows.residual_slope, against, point);
    point_whole = report.whole;
    if (!point_whole) {
      point_why = std::string("the operating point cannot be recorded at a ") +
                  Leaf::operating_point_kind_name(rows.kind) + " point: " +
                  report.why;
      if (!rows.message.empty()) {
        point_why += "; " + rows.message;
      }
    }
  }

  for (std::size_t j = 0; j < output.size(); ++j) {
    S& into = output[j] == grad::out_profit
                  ? leaf_profit_
                  : leaf_soil_consumption_[static_cast<std::size_t>(
                        output[j] - grad::out_uptake_first)];
    // The value is the LEAF's, not the row layer's: the read moves nothing, so the
    // solved leaf is still holding it and a second copy could only disagree.
    const double value =
        output[j] == grad::out_profit
            ? leaf.profit_
            : leaf.soil_consumption_[static_cast<std::size_t>(
                  output[j] - grad::out_uptake_first)];
    const bool objective = grad::role_of(output[j]) == grad::Role::Objective;
    // ⚠️ THE FLAG SUPPRESSES EVERY LATER PLANT'S WATER ROWS AND NOT ONLY THIS
    // ONE'S. A metric's gradient is a sum, so a term missing anywhere leaves the
    // whole water channel undefined and taping the rest would produce a gradient
    // that is partly an answer. The values go on flowing either way.
    //
    // Read BEFORE the graft rather than after it, because rows already on the
    // tape cannot be taken off it.
    if (*uptake_rows_unavailable && !objective) {
      into = S(value);
      continue;
    }
    // A lost point costs every output that reads it, and nothing is grafted for
    // those: a value carrying its held rows without the point's channel is the
    // channel gone missing with every number finite.
    if (rows.dy_dp[j] != 0.0 && !point_whole) {
      into = S(value);
      if (objective) {
        throw gradient_refusal("TF24 gradient: `" +
                               grad::output_name(output[j], n_layer_named) +
                               "` reads " + point_why);
      }
      lose_uptake_rows("TF24 gradient: `" +
                       grad::output_name(output[j], n_layer_named) +
                       "` reads " + point_why);
      continue;
    }
    against.clear();
    const bool carries_point = rows.dy_dp[j] != 0.0;
    if (carries_point) {
      against.push_back({point, rows.dy_dp[j]});
    }
    for (std::size_t i = 0; i < n_input; ++i) {
      against.push_back({*at[i], rows.held[j * n_input + i]});
    }
    // The VALUE survives a missing row and the patch balance still needs it, so
    // the graft leaves it carrying nothing rather than carrying part -- which is
    // what the flag below then says, against a zero row saying the opposite.
    const odelia::graft_report report =
        odelia::record_with_derivatives<S>(value, against, into);
    if (!report.whole) {
      if (objective) {
        throw gradient_refusal(missing(j, report, carries_point));
      }
      lose_uptake_rows(missing(j, report, carries_point));
    }
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
S TF24_Strategy<S>::area_leaf(const S& height) const {
  return pow(height / pars.a_l1, 1.0 / pars.a_l2);
}

// [eqn 1] mass_leaf (inverse of [eqn 2])
template <typename S>
S TF24_Strategy<S>::mass_leaf(const S& area_leaf) const {
  return area_leaf * pars.lma;
}

// [eqn 4] area and mass of sapwood
template <typename S>
S TF24_Strategy<S>::area_sapwood(const S& area_leaf) const {
  return area_leaf * pars.theta;
}

template <typename S>
S TF24_Strategy<S>::mass_sapwood(const S& area_sapwood, const S& height) const {
  return area_sapwood * height * eta_c * pars.rho;
}

// [eqn 5] area and mass of bark
template <typename S>
S TF24_Strategy<S>::area_bark(const S& area_leaf) const {
  return pars.a_b1 * area_leaf * pars.theta;
}

template <typename S>
S TF24_Strategy<S>::mass_bark(const S& area_bark, const S& height) const {
  return area_bark * height * eta_c * pars.rho;
}

template <typename S>
S TF24_Strategy<S>::area_stem(const S& area_bark, const S& area_sapwood,
                            const S& area_heartwood) const {
  return area_bark + area_sapwood + area_heartwood;
}

template <typename S>
S TF24_Strategy<S>::diameter_stem(const S& area_stem) const {
  using std::sqrt;
  return sqrt(4 * area_stem / M_PI);
}

// [eqn 7] Mass of (fine) roots
template <typename S>
S TF24_Strategy<S>::mass_root(const S& area_leaf) const {
  return pars.a_r1 * area_leaf;
}

// [eqn 8] Total mass
template <typename S>
S TF24_Strategy<S>::mass_live(const S& mass_leaf, const S& mass_bark,
                           const S& mass_sapwood, const S& mass_root) const {
  return mass_leaf + mass_sapwood + mass_bark + mass_root;
}

template <typename S>
S TF24_Strategy<S>::mass_total(const S& mass_leaf, const S& mass_bark,
                            const S& mass_sapwood, const S& mass_heartwood,
                            const S& mass_root) const {
  return mass_leaf + mass_bark + mass_sapwood +  mass_heartwood + mass_root;
}

template <typename S>
S TF24_Strategy<S>::mass_above_ground(const S& mass_leaf, const S& mass_bark,
                            const S& mass_sapwood, const S& mass_heartwood) const {
  return mass_leaf + mass_bark + mass_sapwood + mass_heartwood;
}

// for updating auxiliary state
template <typename S>
void TF24_Strategy<S>::update_dependent_aux(const int index, Internals<S>& vars) {
  if (index == HEIGHT_INDEX) {
    const S& height = vars.state(HEIGHT_INDEX);
    vars.set_aux(aux_idx_competition_effect, area_leaf(height));
    vars.set_aux(aux_idx_height_inverse, 1.0 / height);
  }
}


// one-shot update of the scm variables
// i.e. setting rates of ode vars from the state and updating aux vars
template <typename S>
void TF24_Strategy<S>::compute_rates(const TF24_Environment<S>& environment,  Internals<S>& vars) {
  const S& height = vars.state(HEIGHT_INDEX);
  const S& area_leaf_ = vars.aux(aux_idx_competition_effect);

  const S net_mass_production_dt_ =
    net_mass_production_dt(environment, height, area_leaf_,
                           vars.aux(aux_idx_height_inverse));

  // Read by establishment_probability for the boundary node, so this one is a
  // rate's input rather than a reading of it.
  vars.set_aux(aux_idx_net_mass_production_dt, net_mass_production_dt_);

  // The rest are readings: their only consumers are r_internals and ode_aux,
  // both R-facing. Writing them at double alone was tried and measured at no
  // difference, because a store into a slot the vector already holds pushes a
  // statement and registers nothing -- it is a newly constructed active scalar
  // that costs a slot, not a write to a live one.
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
  // Counted: below zero the pool stops being a function of its own state, so
  // every rate the gate scales loses its route to the reserve.
  if (vars.state(state_idx_storage) < S(0.0)) {
    note_clamp(CLAMP_STORAGE_FLOOR);
  }
  const S storage     = std::max(vars.state(state_idx_storage), S(0.0));
  const S storage_max = storage_capacity(area_leaf_, height);
  using std::exp;
  using std::min;
  using std::sqrt;
  // And at capacity the gate's argument stops reading the pool, which is the
  // same severance at the other end.
  if (storage_max > 0.0 && storage > storage_max) {
    note_clamp(CLAMP_RESERVE_CEILING);
  }
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
                                    const S& height,
                                    const S& area_leaf) {


  S A = 0.0;

  // Define an anonymous function to integrate
  // For given height in crown, take photosynthesis at depth multipled by 
  //   amount of leaf at that depth
  std::function<S(S)> f = [&](const S& z) -> S {
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
S TF24_Strategy<S>::assimilation_leaf(const S& x) const {
  return pars.a_p1 * x / (x + pars.a_p2);
}

// [eqn 13] Total maintenance respiration
// NOTE: In contrast with Falster ref model, we do not normalise by pars.a_y*pars.a_bio.
template <typename S>
S TF24_Strategy<S>::respiration(const S& mass_leaf, const S& mass_sapwood,
                             const S& mass_bark, const S& mass_root) const {
  return respiration_leaf(mass_leaf) +
         respiration_bark(mass_bark) +
         respiration_sapwood(mass_sapwood) +
         respiration_root(mass_root);
}

template <typename S>
S TF24_Strategy<S>::respiration_leaf(const S& mass) const {
  return pars.r_l * mass;
}

template <typename S>
S TF24_Strategy<S>::respiration_bark(const S& mass) const {
  return pars.r_b * mass;
}

template <typename S>
S TF24_Strategy<S>::respiration_sapwood(const S& mass) const {
  return pars.r_s * mass;
}

template <typename S>
S TF24_Strategy<S>::respiration_root(const S& mass) const {
  return pars.r_r * mass;
}

// [eqn 14] Total turnover
template <typename S>
S TF24_Strategy<S>::turnover(const S& mass_leaf, const S& mass_bark,
                          const S& mass_sapwood, const S& mass_root) const {
   return turnover_leaf(mass_leaf) +
          turnover_bark(mass_bark) +
          turnover_sapwood(mass_sapwood) +
          turnover_root(mass_root);
}

template <typename S>
S TF24_Strategy<S>::turnover_leaf(const S& mass) const {
  return pars.k_l * mass;
}

template <typename S>
S TF24_Strategy<S>::turnover_bark(const S& mass) const {
  return pars.k_b * mass;
}

template <typename S>
S TF24_Strategy<S>::turnover_sapwood(const S& mass) const {
  return pars.k_s * mass;
}

template <typename S>
S TF24_Strategy<S>::turnover_root(const S& mass) const {
  return pars.k_r * mass;
}

// [eqn 15] Net production
//
// NOTE: Translation of variable names from the Falster 2011.  Everything
// before the minus sign is SCM's N, our `net_mass_production_dt` is SCM's P.
template <typename S>
S TF24_Strategy<S>::net_mass_production_dt_A(const S& assimilation, const S& respiration,
                                const S& turnover) const {
  return pars.a_bio * pars.a_y * (assimilation - respiration) - turnover;
}

// One shot calculation of net_mass_production_dt
// Used by establishment_probability() and compute_rates().
template <typename S>
S TF24_Strategy<S>::net_mass_production_dt(const TF24_Environment<S>& environment,
                                const S& height, const S& area_leaf_,
                                const S& height_inverse) {
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

  // Counted, and this one severs a state row rather than smoothing a number:
  // above the cap the whole root profile stops depending on height, so every
  // layer's carbon -- and the uptake built on it -- reads a constant.
  if (height > pars.rooting_depth_max) {
    note_clamp(CLAMP_ROOTING_DEPTH);
  }
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
      beta_R_V, leaf_drivers_.root_network);

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
      leaf.set_physiology(leaf_drivers_.root_network, radiation, psi_soil,
                          soil_depths_, leaf_specific_conductance_max,
                          environment.get_atm_vpd(), environment.get_ca(),
                          environment.get_leaf_temp(),
                          environment.get_atm_o2_kpa(),
                          environment.get_atm_kpa());
    } else {
      using odelia::util::to_passive;
      // One object holds what the leaf is driven with, so this supply and the
      // perturbed ones its rows are read through cannot disagree about any of
      // it.
      leaf_drivers_.psi_soil.resize(psi_soil.size());
      for (size_t a = 0; a < psi_soil.size(); ++a) {
        leaf_drivers_.psi_soil[a] = to_passive(psi_soil[a]);
      }
      leaf_drivers_.PPFD = to_passive(radiation);
      leaf_drivers_.soil_depth = soil_depths_;
      leaf_drivers_.atm_vpd = environment.get_atm_vpd();
      leaf_drivers_.ca = environment.get_ca();
      leaf_drivers_.leaf_temp = environment.get_leaf_temp();
      leaf_drivers_.atm_o2_kpa = environment.get_atm_o2_kpa();
      leaf_drivers_.atm_kpa = environment.get_atm_kpa();
      leaf.set_physiology(leaf_drivers_.root_network, leaf_drivers_.PPFD,
                          leaf_drivers_.psi_soil, leaf_drivers_.soil_depth,
                          to_passive(leaf_specific_conductance_max),
                          leaf_drivers_.atm_vpd, leaf_drivers_.ca,
                          leaf_drivers_.leaf_temp, leaf_drivers_.atm_o2_kpa,
                          leaf_drivers_.atm_kpa);
    }
    solve_leaf();
  };

  // Convert canopy openness (0-1) into absorbed radiation: PPFD attenuated by
  // the self-shading coefficient pars.k_I. The light floor (1e-4) matches
  // compute_average_light_environment().
  const double PPFD = environment.get_PPFD();
  auto radiation_at = [&](const S& light) -> S {
    // Counted, not merely applied. Where this binds the cohort's radiation is a
    // constant with respect to every other cohort's height, so the row is
    // severed by the guard rather than by the model -- and the field is smooth
    // underneath. An uncounted severance is indistinguishable from a true zero.
    if (light < S(0.0001)) {
      // Below the floor the radiation this cohort receives stops depending on
      // any other cohort's height. The row is therefore exactly zero for the
      // model AS EVALUATED -- every light below the floor gives a bit-identical
      // census -- so the zero is the derivative of the function on the tape and
      // not an answer withheld. What makes it readable rather than silent is the
      // count, which is taken on this path because it is the only one a clamp
      // severs anything on.
      note_clamp(CLAMP_LIGHT_FLOOR);
    }
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
    auto f = [&](const S& x) -> S {
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
    record_leaf_outputs(radiation_used, psi_soil,
                        leaf_specific_conductance_max);
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

// Base TF24: place the operating point the run found at this rate evaluation, or
// optimise the root-collar water potential where it kept none -- a branch that exited
// on feasibility, or an evaluation the run addressed no record against.
template <typename S>
void TF24_Strategy<S>::solve_leaf() {
  if (!leaf.place_solved_point(leaf_points->next())) {
    leaf.find_root_collar_psi();
  }
  leaf_points->keep(leaf.solved_point());
  // The classification is decided by the branch taken and then overwritten by
  // the next plant, so without a tally the only route to its incidence is a
  // refusal message -- which reports the FIRST non-interior point and nothing
  // about how many followed it or what kinds they were.
  ++operating_point_counts[
      static_cast<size_t>(leaf.operating_point_kind())];
}

// [eqn 16] Fraction of production allocated to reproduction
template <typename S>
S TF24_Strategy<S>::fraction_allocation_reproduction(const S& height) const {
  return pars.a_f1 / (1.0 + exp(pars.a_f2 * (1.0 - height / pars.hmat)));
}

// Fraction of production allocated to growth
template <typename S>
S TF24_Strategy<S>::fraction_allocation_growth(const S& height) const {
  return 1.0 - fraction_allocation_reproduction(height);
}

// [eqn 17] Rate of offspring production
template <typename S>
S TF24_Strategy<S>::fecundity_dt(const S& net_mass_production_dt,
                               const S& fraction_allocation_reproduction) const {
  return net_mass_production_dt * fraction_allocation_reproduction /
    (pars.omega + pars.a_f3);
}

template <typename S>
S TF24_Strategy<S>::darea_leaf_dmass_live(const S& area_leaf) const {
  return 1.0/(  dmass_leaf_darea_leaf(area_leaf)
              + dmass_sapwood_darea_leaf(area_leaf)
              + dmass_bark_darea_leaf(area_leaf)
              + dmass_root_darea_leaf(area_leaf));
}

template <typename S>
S TF24_Strategy<S>::dheight_darea_leaf(const S& area_leaf) const {
  return pars.a_l1 * pars.a_l2 * pow(area_leaf, pars.a_l2 - 1);
}

// Mass of leaf needed for new unit area leaf, d m_s / d a_l
template <typename S>
S TF24_Strategy<S>::dmass_leaf_darea_leaf(const S& /* area_leaf */) const {
  return pars.lma;
}

// Mass of stem needed for new unit area leaf, d m_s / d a_l
template <typename S>
S TF24_Strategy<S>::dmass_sapwood_darea_leaf(const S& area_leaf) const {
  return pars.rho * eta_c * pars.a_l1 * pars.theta * (pars.a_l2 + 1.0) * pow(area_leaf, pars.a_l2);
}

// Mass of bark needed for new unit area leaf, d m_b / d a_l
template <typename S>
S TF24_Strategy<S>::dmass_bark_darea_leaf(const S& area_leaf) const {
  return pars.a_b1 * dmass_sapwood_darea_leaf(area_leaf);
}

// Mass of root needed for new unit area leaf, d m_r / d a_l
template <typename S>
S TF24_Strategy<S>::dmass_root_darea_leaf(const S& /* area_leaf */) const {
  return pars.a_r1;
}

// Growth rate of basal diameter_stem per unit time
template <typename S>
S TF24_Strategy<S>::ddiameter_stem_darea_stem(const S& area_stem) const {
  return pow(M_PI * area_stem, -0.5);
}

// Growth rate of sapwood area at base per unit time
template <typename S>
S TF24_Strategy<S>::area_sapwood_dt(const S& area_leaf_dt) const {
  return area_leaf_dt * pars.theta;
}

// Note, unlike others, heartwood growth does not depend on leaf area growth, but
// rather existing sapwood
template <typename S>
S TF24_Strategy<S>::area_heartwood_dt(const S& area_leaf) const {
  return pars.k_s * area_sapwood(area_leaf);
}

// Growth rate of bark area at base per unit time
template <typename S>
S TF24_Strategy<S>::area_bark_dt(const S& area_leaf_dt) const {
  return pars.a_b1 * area_leaf_dt * pars.theta;
}

// Growth rate of stem basal area per unit time
template <typename S>
S TF24_Strategy<S>::area_stem_dt(const S& area_leaf,
                               const S& area_leaf_dt) const {
  return area_sapwood_dt(area_leaf_dt) +
    area_bark_dt(area_leaf_dt) +
    area_heartwood_dt(area_leaf);
}

// Growth rate of basal diameter_stem per unit time
template <typename S>
S TF24_Strategy<S>::diameter_stem_dt(const S& area_stem, const S& area_stem_dt) const {
  return ddiameter_stem_darea_stem(area_stem) * area_stem_dt;
}

// Growth rate of root mass per unit time
template <typename S>
S TF24_Strategy<S>::mass_root_dt(const S& area_leaf,
                               const S& area_leaf_dt) const {
  return area_leaf_dt * dmass_root_darea_leaf(area_leaf);
}

template <typename S>
S TF24_Strategy<S>::mass_live_dt(const S& fraction_allocation_reproduction,
                               const S& net_mass_production_dt) const {
  return (1 - fraction_allocation_reproduction) * net_mass_production_dt;
}

template <typename S>
S TF24_Strategy<S>::mass_total_dt(const S& fraction_allocation_reproduction,
                                     const S& net_mass_production_dt,
                                     const S& mass_heartwood_dt) const {
  return mass_live_dt(fraction_allocation_reproduction, net_mass_production_dt) +
    mass_heartwood_dt;
}

// TODO: Do we not track root mass change?
template <typename S>
S TF24_Strategy<S>::mass_above_ground_dt(const S& area_leaf,
                                       const S& fraction_allocation_reproduction,
                                       const S& net_mass_production_dt,
                                       const S& mass_heartwood_dt,
                                       const S& area_leaf_dt) const {
  const S mass_root_dt =
    area_leaf_dt * dmass_root_darea_leaf(area_leaf);
  return mass_total_dt(fraction_allocation_reproduction, net_mass_production_dt,
                        mass_heartwood_dt) - mass_root_dt;
}

template <typename S>
S TF24_Strategy<S>::mass_heartwood_dt(const S& mass_sapwood) const {
  return turnover_sapwood(mass_sapwood);
}


template <typename S>
S TF24_Strategy<S>::mass_live_given_height(const S& height) const {
  S area_leaf_ = area_leaf(height);
  return mass_leaf(area_leaf_) +
         mass_bark(area_bark(area_leaf_), height) +
         mass_sapwood(area_sapwood(area_leaf_), height) +
         mass_root(area_leaf_);
}

template <typename S>
S TF24_Strategy<S>::height_given_mass_leaf(const S& mass_leaf) const {
  return pars.a_l1 * pow(mass_leaf / pars.lma, pars.a_l2);
}

template <typename S>
S TF24_Strategy<S>::mortality_dt(const S& relative_reserves,
                              const S& cumulative_mortality) const {

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
S TF24_Strategy<S>::mortality_storage_dependent_dt(const S& relative_reserves) const {
  return pars.a_dG1 * exp(-pars.a_dG2 * relative_reserves);
}

// NSC storage capacity: scales with sapwood mass (per Daniel, #517). mass_sapwood
// = area_sapwood(area_leaf) * height * eta_c * rho.
template <typename S>
S TF24_Strategy<S>::storage_capacity(const S& area_leaf_, const S& height) const {
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
                                               const S& net_mass_production_dt_) {

  S decay_over_time = exp(-pars.recruitment_decay * environment.time);

  if (net_mass_production_dt_ > 0) {
    const S tmp = pars.a_d0 * seed_geometry().area_leaf / net_mass_production_dt_;
    return 1.0 / (tmp * tmp + 1.0) * decay_over_time;
  } else {
    return 0.0;
  }
}

template <typename S>
S TF24_Strategy<S>::compute_competition(const S& z, const S& height) const {
  return pars.k_I * area_leaf(height) * canopy_shape.Q_from_height(z, height);
}

// Ratio-first hot-path overload (see header): receives the cached
// competition_effect (= area_leaf(height)) and height_inverse (= 1/height), so the
// per-call area_leaf() evaluation and z/height division are hoisted out of the
// inner competition loop.
template <typename S>
S TF24_Strategy<S>::compute_competition(const S& z, const S& area_leaf_,
                                          const S& height_inverse) const {
  return pars.k_I * area_leaf_ * canopy_shape.leaf_area_above(z * height_inverse);
}

// [eqn 10] Cumulative fraction of a quantity distributed over an extent with
//          shape exponent 'eta_x', above coordinate 'z' of a total 'height'.
//          Serves the root mass distribution over soil depth.
template <typename S>
S TF24_Strategy<S>::Q(const S& z, const S& height, const S& eta_x) const {
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
