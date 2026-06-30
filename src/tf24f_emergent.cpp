// Census-metric reconstruction for TF24f (#472 scope B, build-order step 1 -- the
// TF24f frozen-census R0 GATE). No AD: a pure double-precision replay that proves
// the collar-state trajectory + the census number-density (log_density) reproduce
// the SCM's stored stand before any reverse-mode tape is built (the R1 follow-up).
//
// Why TF24f, not TF24, gets the census gradient (the scope decision in
// notes/tf24-stand-gradient-scope.md): TF24's census number density needs the
// SECOND-order leaf-optimiser sensitivity d(growth-rate-gradient)/d(theta), which the
// linearised leaf-opt harvest of tf24_emergent.cpp does not differentiate faithfully.
// TF24f removes the per-step golden-section optimiser entirely: the optimal root-collar
// potential is a 6th ODE state (opt_root_psi_state) that relaxes toward its optimum by
// gradient ascent (dpsi/dt = k_acclim * dprofit_dpsi), and the leaf is evaluated AT the
// tracked collar by an analytic, IFT-able path (Leaf::dprofit_droot_collar_psi). So the
// trajectory is a clean 6-state ODE with differentiable rates and the census g' needs no
// curvature harvest -- carrying the collar as a replayed state propagates the height
// sensitivity directly (the R1 tape).
//
// This R0 gate establishes the prerequisite the scope note flags first (§7): that the
// collar-state replay is FAITHFUL. It re-evolves, per cohort, the 5 demographic states +
// the tracked collar + log_density over the SCM's frozen Cash-Karp schedule, reading the
// frozen resident environment and driving the REAL TF24f leaf evaluation at the tracked
// collar (crown-centre: one analytic leaf solve per stage, no optimiser). The census
// number density uses the SCM's own growth-rate-gradient scheme exactly: a backward
// finite difference of height_dt with step node_gradient_eps = 1e-6, holding the collar
// fixed and perturbing only height (Node::growth_rate_gradient with the default
// node_gradient_direction = -1, exact_ad off). It returns the reconstructed per-cohort
// heights / collar states / log-densities and the census metric values (LAI / biomass /
// size_moment); the R gate checks these against the SCM's stored stand and
// compute_competition(0).
#include <Rcpp.h>
#include <vector>
#include <string>
#include <algorithm>
#include <cmath>
#include <XAD/XAD.hpp>             // reverse-mode adjoint tape (the R1 census gradient)
#include <plant.h>                 // RcppR6 as<>/wrap for TF24_Environment
#include <plant/models/tf24f_strategy.h>
#include <plant/models/tf24_production_kernel.h>
#include <plant/models/ff16_production_kernel.h>   // ff16_cashkarp_replay, FF16State

namespace {

// Build a TF24f strategy from the SCM's parameter vector. Mirrors the TF24
// make_strategy in tf24_emergent.cpp but constructs the fast-acclimation variant and
// sets its acclimation knobs. crown-centre + GSS_tol_abs match the resident run.
plant::TF24f_Strategy make_tf24f(const Rcpp::NumericVector& pp, double k_acclim,
                                 bool use_ad_gradient,
                                 const std::string& shading = "crown-centre",
                                 double gss_tol = 1e-9) {
  plant::TF24f_Strategy s;
  s.control.shading_model = shading; s.control.GSS_tol_abs = gss_tol;
  s.k_acclim = k_acclim; s.use_ad_gradient = use_ad_gradient;
  auto& q = s.pars;
  q.lma=pp["lma"];q.rho=pp["rho"];q.hmat=pp["hmat"];q.omega=pp["omega"];q.eta=pp["eta"];
  q.theta=pp["theta"];q.a_l1=pp["a_l1"];q.a_l2=pp["a_l2"];q.a_r1=pp["a_r1"];q.a_b1=pp["a_b1"];
  q.r_s=pp["r_s"];q.r_b=pp["r_b"];q.r_r=pp["r_r"];q.r_l=pp["r_l"];q.a_y=pp["a_y"];q.a_bio=pp["a_bio"];
  q.k_l=pp["k_l"];q.k_b=pp["k_b"];q.k_s=pp["k_s"];q.k_r=pp["k_r"];
  q.a_f3=pp["a_f3"];q.a_f1=pp["a_f1"];q.a_f2=pp["a_f2"];q.d_I=pp["d_I"];
  q.a_dG1=pp["a_dG1"];q.a_dG2=pp["a_dG2"];q.k_I=pp["k_I"];q.a_d0=pp["a_d0"];
  q.recruitment_decay=pp["recruitment_decay"];
  q.vcmax_25=pp["vcmax_25"];q.K_s=pp["K_s"];q.b=pp["b"];q.c=pp["c"];q.beta2=pp["beta2"];
  q.jmax_25=pp["jmax_25"];q.a=pp["a"];q.curv_fact_elec_trans=pp["curv_fact_elec_trans"];
  q.curv_fact_colim=pp["curv_fact_colim"];
  s.prepare_strategy();
  return s;
}

// The 7-state cohort: TF24's 5 demographic states + the tracked root-collar potential
// (TF24f's 6th ODE state) + the census number density log_density.
struct St7 { plant::FF16State<double> demog; double collar; double log_density; };

// Evaluate the leaf at a GIVEN tracked collar potential, returning net production and
// (out) the profit gradient dprofit/dpsi that drives the collar's gradient-ascent rate.
// Drives the real TF24f leaf path: setting tracked_root_psi_ then calling the inherited
// net_mass_production_dt makes solve_leaf() evaluate at the tracked collar (clamped to
// the feasible interval) rather than re-optimise -- the analytic, IFT-able TF24f path.
// crown-centre runs exactly one leaf solve, so dprofit_dpsi_ is set cleanly.
double net_at_tracked(plant::TF24f_Strategy& s, const plant::TF24_Environment& e,
                      double h, double collar, double& dprofit_dpsi) {
  s.tracked_root_psi_ = collar;
  const double net = s.net_mass_production_dt(e, h, s.area_leaf(h), 1.0 / h);
  dprofit_dpsi = s.dprofit_dpsi_;
  return net;
}

} // namespace



// Compiled core of tf24f_census_recon(). Re-evolves every cohort's {5 demog, collar,
// log_density} over the harvested frozen schedule (per-RK-stage TF24_Environment in
// eh_list, step times sh, per-cohort birth steps), driving the real TF24f leaf at the
// tracked collar in double precision. Returns the reconstructed per-cohort final heights
// / collar states / log-densities and the census metric values (the height-trapezium of
// density*kI*area_leaf etc., matching Patch::compute_competition / FF16's census reduce).
// [[Rcpp::export]]
Rcpp::List tf24f_census_recon_impl(
    Rcpp::NumericVector pp, Rcpp::List eh_list, std::vector<double> sh,
    std::vector<int> birth, double birth_rate, double k_acclim,
    bool use_ad_gradient, std::vector<std::string> metrics, bool exact_ad_gprime) {
  plant::TF24f_Strategy s = make_tf24f(pp, k_acclim, use_ad_gradient);
  plant::TF24ProdPars<double> pd = s.prod_pars();
  const double h0v = s.initial_height(), a_d0 = s.pars.a_d0;
  const double kI = s.pars.k_I, recruitment_decay = s.pars.recruitment_decay;
  const double GEPS = 1e-6;                     // Control::node_gradient_eps (backward FD)

  // Validate metrics up front.
  for (auto& nm : metrics)
    if (nm != "LAI" && nm != "biomass" && nm != "size_moment")
      Rcpp::stop("unknown TF24f census metric: " + nm +
                 " (LAI / biomass / size_moment)");

  const std::size_t N = eh_list.size(), nC = birth.size();
  std::vector<std::vector<plant::TF24_Environment>> EH(N);
  for (std::size_t n = 0; n < N; ++n) {
    Rcpp::List st = eh_list[n];
    for (R_xlen_t k = 0; k < st.size(); ++k)
      EH[n].push_back(Rcpp::as<plant::TF24_Environment>(st[k]));
  }
  std::vector<double> step_h(N);
  for (std::size_t n = 0; n < N; ++n) step_h[n] = sh[n + 1] - sh[n];

  // Replay one cohort to its final state over the frozen schedule.
  auto replay_cohort = [&](std::size_t i) -> St7 {
    using std::exp; using std::log;
    const std::size_t b = (std::size_t)birth[i];
    const plant::TF24_Environment& eb = (b > 0) ? EH[b - 1][5] : EH[0][0];

    // --- Seed (Node::compute_initial_conditions) -----------------------------------
    // 1. collar_0 = the optimum at the birth env (set_initial_states runs the base
    //    golden-section optimiser ONCE via initializing_).
    s.initializing_ = true;
    s.net_mass_production_dt(eb, h0v, s.area_leaf(h0v), 1.0 / h0v);
    const double collar0 = -s.leaf.root_collar_psi_;
    s.initializing_ = false;
    // 2. Establishment + seedling growth rate, evaluated at the tracked collar_0 (the
    //    operating point the SCM uses: compute_rates set tracked_root_psi_ = collar_0
    //    before establishment_probability).
    const double decay = std::exp(-recruitment_decay * sh[b]);
    double dpsi0 = 0.0;
    const double area0 = plant::tf24_area_leaf<double>(pd.a_l1, pd.a_l2, h0v);
    const double net0 = net_at_tracked(s, eb, h0v, collar0, dpsi0);
    double pr_estab = 0.0;
    if (net0 > 0.0) { const double uu = a_d0 * area0 / net0; pr_estab = decay / (uu * uu + 1.0); }
    const double mort0 = (pr_estab > 0.0) ? -log(pr_estab)
                                          : std::numeric_limits<double>::infinity();
    const double g0 = plant::tf24_height_dt_from_net<double>(pd, h0v, area0, net0);
    const double logd0 = (g0 > 0.0 && pr_estab > 0.0)
                           ? log(birth_rate * pr_estab / g0)
                           : -std::numeric_limits<double>::infinity();

    // --- Cash-Karp replay over the frozen schedule ---------------------------------
    auto deriv = [&](const St7& y, std::size_t n, int stage) -> St7 {
      const plant::TF24_Environment& e =
        (stage == 0) ? ((n > 0) ? EH[n - 1][5] : EH[0][0]) : EH[n][stage - 1];
      const double h = y.demog.height, collar = y.collar;
      double dpsi = 0.0;
      const double net = net_at_tracked(s, e, h, collar, dpsi);
      const double al = plant::tf24_area_leaf<double>(pd.a_l1, pd.a_l2, h);
      plant::TF24Rates<double> r = plant::tf24_compute_rates_from_net<double>(pd, h, al, net, true);
      const double collar_dt = k_acclim * dpsi;
      // g' matches the SCM's Node::growth_rate_gradient, collar held fixed. Either the
      // exact-AD gradient (node_gradient_exact_ad = TRUE, the analytic kmax+light+E_up
      // channels) or the backward FD (dir = -1, step GEPS) -- whichever the resident
      // run used -- so the reconstructed density reproduces the SCM's stored density.
      double gprime;
      if (exact_ad_gprime) {
        s.tracked_root_psi_ = collar;
        gprime = s.growth_rate_gradient_height_ad(h, e);
        // restore the operating point at (h, collar) for the rate fill below
        net_at_tracked(s, e, h, collar, dpsi);
      } else {
        double dpsi_b = 0.0;
        const double net_b = net_at_tracked(s, e, h - GEPS, collar, dpsi_b);
        const double al_b = plant::tf24_area_leaf<double>(pd.a_l1, pd.a_l2, h - GEPS);
        const double g_back = plant::tf24_height_dt_from_net<double>(pd, h - GEPS, al_b, net_b);
        gprime = (r.height_dt - g_back) / GEPS;
      }
      const double log_density_dt = -gprime - r.mortality_dt;
      return St7{plant::FF16State<double>{r.height_dt, r.mortality_dt, r.fecundity_dt,
        r.area_heartwood_dt, r.mass_heartwood_dt}, collar_dt, log_density_dt};
    };
    auto axpy = [](const St7& a, double c, const St7& k) -> St7 {
      return St7{plant::FF16State<double>{
        a.demog.height + c * k.demog.height, a.demog.mortality + c * k.demog.mortality,
        a.demog.fecundity + c * k.demog.fecundity,
        a.demog.area_heartwood + c * k.demog.area_heartwood,
        a.demog.mass_heartwood + c * k.demog.mass_heartwood},
        a.collar + c * k.collar, a.log_density + c * k.log_density};
    };
    St7 y{plant::FF16State<double>{h0v, mort0, 0.0, 0.0, 0.0}, collar0, logd0};
    return plant::ff16_cashkarp_replay(y, step_h, b, deriv, axpy);
  };

  std::vector<St7> finals(nC);
  for (std::size_t i = 0; i < nC; ++i) finals[i] = replay_cohort(i);

  // The pending-seed (new_node) census tail term: a seed born into the FINAL-time env
  // (eh.back()[5]), giving density_new for the trapezium tail from h_last down to h0.
  double dens_new = 0.0;
  {
    const plant::TF24_Environment& ef = EH[N - 1][5];
    // collar seed at the final env
    s.initializing_ = true;
    s.net_mass_production_dt(ef, h0v, s.area_leaf(h0v), 1.0 / h0v);
    const double collar0 = -s.leaf.root_collar_psi_;
    s.initializing_ = false;
    const double decay = std::exp(-recruitment_decay * sh[N]);
    double dpsi0 = 0.0;
    const double area0 = plant::tf24_area_leaf<double>(pd.a_l1, pd.a_l2, h0v);
    const double net0 = net_at_tracked(s, ef, h0v, collar0, dpsi0);
    double pr_estab = 0.0;
    if (net0 > 0.0) { const double uu = a_d0 * area0 / net0; pr_estab = decay / (uu * uu + 1.0); }
    const double g0 = plant::tf24_height_dt_from_net<double>(pd, h0v, area0, net0);
    dens_new = (g0 > 0.0 && pr_estab > 0.0) ? birth_rate * pr_estab / g0 : 0.0;
  }

  // --- Census reduction (height-trapezium, descending height + pending-seed tail) ---
  // Mirrors ff16_emergent.cpp's census_reduce / Species::compute_competition exactly.
  std::vector<std::size_t> ord(nC);
  for (std::size_t i = 0; i < nC; ++i) ord[i] = i;
  std::sort(ord.begin(), ord.end(),
            [&](std::size_t a, std::size_t b){ return finals[a].demog.height > finals[b].demog.height; });

  auto census_reduce = [&](auto psi) -> double {
    std::vector<double> phi(nC);
    for (std::size_t i = 0; i < nC; ++i)
      phi[i] = psi(finals[i].demog.height, std::exp(finals[i].log_density),
                   finals[i].demog.mass_heartwood);
    double J = 0.0;
    for (std::size_t j = 0; j + 1 < nC; ++j) {
      const std::size_t a = ord[j], b = ord[j + 1];
      J += 0.5 * (finals[a].demog.height - finals[b].demog.height) * (phi[a] + phi[b]);
    }
    if (nC > 0) {
      const std::size_t last = ord[nC - 1];
      const double phi_new = psi(h0v, dens_new, 0.0);
      J += 0.5 * (finals[last].demog.height - h0v) * (phi[last] + phi_new);
    }
    return J;
  };

  Rcpp::NumericVector values(metrics.size());
  for (std::size_t m = 0; m < metrics.size(); ++m) {
    const std::string& nm = metrics[m];
    if (nm == "LAI") {
      values[m] = census_reduce([&](double h, double dens, double) -> double {
        return dens * kI * plant::tf24_area_leaf<double>(pd.a_l1, pd.a_l2, h); });
    } else if (nm == "biomass") {
      // live mass + heartwood (mirrors FF16 biomass); TF24 mass_live via the cascade.
      values[m] = census_reduce([&](double h, double dens, double mhw) -> double {
        const double al = plant::tf24_area_leaf<double>(pd.a_l1, pd.a_l2, h);
        const double mass_leaf    = al * pd.lma;
        const double area_sapwood = al * pd.theta;
        const double mass_sapwood = area_sapwood * h * pd.eta_c * pd.rho;
        const double area_bark    = pd.a_b1 * al * pd.theta;
        const double mass_bark    = area_bark * h * pd.eta_c * pd.rho;
        const double mass_root    = pd.a_r1 * al;
        return dens * (mass_leaf + mass_sapwood + mass_bark + mass_root + mhw); });
    } else { // size_moment
      values[m] = census_reduce([&](double h, double dens, double) -> double {
        return dens * h; });
    }
  }
  values.attr("names") = Rcpp::wrap(metrics);

  Rcpp::NumericVector heights(nC), collar(nC), log_density(nC);
  for (std::size_t i = 0; i < nC; ++i) {
    heights[i]     = finals[i].demog.height;
    collar[i]      = finals[i].collar;
    log_density[i] = finals[i].log_density;
  }
  return Rcpp::List::create(
    Rcpp::Named("values") = values,
    Rcpp::Named("heights") = heights,
    Rcpp::Named("collar") = collar,
    Rcpp::Named("log_density") = log_density,
    Rcpp::Named("density_new") = dens_new);
}


// ===========================================================================
// TF24f frozen-census reverse-mode AD trait gradient (#472 scope B, build-order
// step 2 -- the R1 tape, the "refine" step). One reverse sweep over the 7-state
// replay {5 demog, tracked collar, log_density} per metric, replacing the
// per-trait FD over the recon (tf24f_census_gradient_fd) it must reproduce.
//
// The structure extends the TF24 OFFSPRING tape (src/tf24_emergent.cpp): a
// double DISCOVERY pass runs the real TF24f leaf at the tracked collar along the
// frozen-schedule trajectory and HARVESTS the trait-independent operating point
// per RK stage; a second pass replays a leaf-opt-free, fully tapeable expression.
//
// The one ingredient beyond the offspring tape (notes/tf24f-census-tape-seed.md):
// TF24f's tracked collar is a strongly theta-dependent STATE that at k_acclim=1
// LAGS the optimum, so the envelope theorem does NOT zero its contribution. We
// therefore (a) carry the collar as a TAPED state whose gradient-ascent rate
// k_acclim*dprofit_dpsi is linearised with a CURVATURE harvest (the second
// derivatives d2profit/dpsi2, d2profit/dpsi dh, d2profit/dpsi dtheta_k), and
// (b) linearise profit itself in (h, collar, theta) -- the collar term
// dprofit_dpsi0*(collar-collar0) is what the offspring tape omits. The census
// number density's g' = d(height_dt)/d(height) reproduces the SCM's backward-FD
// scheme by harvesting a SECOND operating point at h-GEPS and differencing the
// two tapeable height_dt expressions (as ff16_emergent.cpp's census tape does).
//
// All harvested derivatives are taken by finite difference in the double pass
// (no leaf templating): profit/dprofit_dpsi are differenced w.r.t. collar, height
// and -- via independently-built perturbed strategies -- each requested trait.
// The trait set is a prototype subset (the seed: validate a few traits, not all
// 27); every requested trait must name a strategy `pars` entry (matching the FD
// gate, which perturbs the same parameter vector).
// ===========================================================================

namespace {

using ad   = xad::adj<double>;
using ad_t = ad::active_type;

double as_dbl(const ad_t& v)  { return xad::value(v); }

// Evaluate the real TF24f leaf at a GIVEN (height, tracked collar) in env e,
// returning the leaf profit and (out) the analytic dprofit/dpsi the acclimation
// rate uses. Drives the same path as net_at_tracked / the recon, but exposes
// leaf.profit_ (the quantity the tape linearises and feeds to the kernel).
void eval_leaf(plant::TF24f_Strategy& s, const plant::TF24_Environment& e,
               double h, double collar, double& profit, double& dpsi) {
  s.tracked_root_psi_ = collar;
  s.net_mass_production_dt(e, h, s.area_leaf(h), 1.0 / h);
  profit = s.leaf.profit_;
  dpsi   = s.dprofit_dpsi_;
}

// Trait-independent leaf-operating-point harvest at one (h, collar) point: the
// profit value, its first derivatives w.r.t. height (collar held) and collar, the
// per-trait profit sensitivity inj_k (the leaf channel theta -> profit, the TF24
// `inj` term done by FD over a perturbed strategy), and the curvature the collar
// gradient-ascent rate needs: d(dprofit_dpsi)/d{collar, height, theta_k}.
struct LH {
  double h, collar;                 // the harvested operating point
  double profit, dprofit_dh, dprofit_dpsi;
  double d2p_dpsi2, d2p_dpsidh;     // collar-rate curvature (collar / height)
  std::vector<double> inj;          // d(profit)/d(theta_k)            [per trait]
  std::vector<double> d2p_dpsidth;  // d(dprofit_dpsi)/d(theta_k)      [per trait]
};

// Harvest at (h, collar) in env e. s is the base strategy; sp[k]/sm[k] are the
// +/- trait-perturbed strategies (built once) and dk[k] the matching abs steps.
LH harvest_point(plant::TF24f_Strategy& s,
                 std::vector<plant::TF24f_Strategy>& sp,
                 std::vector<plant::TF24f_Strategy>& sm,
                 const std::vector<double>& dk,
                 const plant::TF24_Environment& e, double h, double collar) {
  const std::size_t T = sp.size();
  LH H; H.h = h; H.collar = collar; H.inj.resize(T); H.d2p_dpsidth.resize(T);
  double dpsi0;
  eval_leaf(s, e, h, collar, H.profit, dpsi0);
  H.dprofit_dpsi = dpsi0;

  // height channel (collar held): central FD of profit + dprofit_dpsi.
  const double dh = 1e-5 * h;
  double pph, dph, pmh, dmh;
  eval_leaf(s, e, h + dh, collar, pph, dph);
  eval_leaf(s, e, h - dh, collar, pmh, dmh);
  H.dprofit_dh = (pph - pmh) / (2.0 * dh);
  H.d2p_dpsidh = (dph - dmh) / (2.0 * dh);

  // collar channel (height held): d2profit/dpsi2 from FD of dprofit_dpsi.
  const double dc = 1e-5 * std::abs(collar) + 1e-8;
  double ppc, dpc, pmc, dmc;
  eval_leaf(s, e, h, collar + dc, ppc, dpc);
  eval_leaf(s, e, h, collar - dc, pmc, dmc);
  H.d2p_dpsi2 = (dpc - dmc) / (2.0 * dc);

  // trait channels: inj_k = d(profit)/d(theta_k); d2p_dpsidth_k = d(dpsi)/d(theta_k),
  // both at the SAME (h, collar) via the prebuilt perturbed strategies.
  for (std::size_t k = 0; k < T; ++k) {
    double pp_, dpp, pm_, dpm;
    eval_leaf(sp[k], e, h, collar, pp_, dpp);
    eval_leaf(sm[k], e, h, collar, pm_, dpm);
    H.inj[k]         = (pp_ - pm_) / (2.0 * dk[k]);
    H.d2p_dpsidth[k] = (dpp - dpm) / (2.0 * dk[k]);
  }
  return H;
}

// Per-stage harvest: the forward operating point (h, collar) and the backward
// point (h - GEPS, collar) the census g' (backward FD of height_dt) reads.
struct StageH { LH fwd, back; };

// Linearised leaf profit at (active height h, active collar) about a harvested
// point Hp, with the requested traits active: profit0 + dprofit_dh*(h - h0)
// + dprofit_dpsi*(collar - collar0) + sum_k inj_k*(theta_k - theta_k0). The
// collar term is the offspring tape's omission (collar is a taped state here).
// h is the active evaluation height; Hp.h the height the point was harvested at
// (for the backward point both already carry the -GEPS offset, so the difference
// h - Hp.h is the active displacement from the harvested trajectory either way).
ad_t profit_lin(const LH& Hp, ad_t h, ad_t collar,
                const std::vector<ad_t>& tr, const std::vector<double>& tr0) {
  ad_t p = ad_t(Hp.profit) + ad_t(Hp.dprofit_dh) * (h - ad_t(Hp.h))
         + ad_t(Hp.dprofit_dpsi) * (collar - ad_t(Hp.collar));
  for (std::size_t k = 0; k < tr.size(); ++k)
    if (Hp.inj[k] != 0.0) p += ad_t(Hp.inj[k]) * (tr[k] - ad_t(tr0[k]));
  return p;
}

// Linearised dprofit_dpsi (the gradient-ascent rate driver) about Hp:
// dpsi0 + d2p_dpsi2*(collar - collar0) + d2p_dpsidh*(h - h0)
//       + sum_k d2p_dpsidth_k*(theta_k - theta_k0).
ad_t dpsi_lin(const LH& Hp, ad_t h, ad_t collar,
              const std::vector<ad_t>& tr, const std::vector<double>& tr0) {
  ad_t d = ad_t(Hp.dprofit_dpsi) + ad_t(Hp.d2p_dpsi2) * (collar - ad_t(Hp.collar))
         + ad_t(Hp.d2p_dpsidh) * (h - ad_t(Hp.h));
  for (std::size_t k = 0; k < tr.size(); ++k)
    if (Hp.d2p_dpsidth[k] != 0.0) d += ad_t(Hp.d2p_dpsidth[k]) * (tr[k] - ad_t(tr0[k]));
  return d;
}

// Build TF24ProdPars<ad_t> from base doubles pd, overwriting any requested trait
// that names a kernel (cascade / area / allometry) field with its active scalar.
// Leaf-only traits (vcmax_25, K_s, b, c, ...) have no kernel field -- they enter
// solely through the harvested profit channel (inj), so they are skipped here.
void set_kernel_trait(plant::TF24ProdPars<ad_t>& p, const std::string& t, const ad_t& v) {
  if      (t=="lma")  p.lma=v;  else if (t=="rho")  p.rho=v;  else if (t=="theta") p.theta=v;
  else if (t=="a_b1") p.a_b1=v; else if (t=="a_r1") p.a_r1=v; else if (t=="r_l")   p.r_l=v;
  else if (t=="r_s")  p.r_s=v;  else if (t=="r_b")  p.r_b=v;  else if (t=="r_r")   p.r_r=v;
  else if (t=="k_l")  p.k_l=v;  else if (t=="k_b")  p.k_b=v;  else if (t=="k_s")   p.k_s=v;
  else if (t=="k_r")  p.k_r=v;  else if (t=="a_bio")p.a_bio=v;else if (t=="a_y")   p.a_y=v;
  else if (t=="a_l1") p.a_l1=v; else if (t=="a_l2") p.a_l2=v; else if (t=="a_f1")  p.a_f1=v;
  else if (t=="a_f2") p.a_f2=v; else if (t=="hmat") p.hmat=v; else if (t=="omega") p.omega=v;
  else if (t=="a_f3") p.a_f3=v; else if (t=="d_I")  p.d_I=v;  else if (t=="a_dG1") p.a_dG1=v;
  else if (t=="a_dG2")p.a_dG2=v;
  // (else: leaf-only trait -> profit channel only)
}

plant::TF24ProdPars<ad_t> pf_active(const plant::TF24ProdPars<double>& pd,
    const std::vector<std::string>& traits, const std::vector<ad_t>& tr) {
  plant::TF24ProdPars<ad_t> p;
  p.lma=pd.lma;p.rho=pd.rho;p.theta=pd.theta;p.a_b1=pd.a_b1;p.a_r1=pd.a_r1;p.eta_c=pd.eta_c;
  p.r_l=pd.r_l;p.r_s=pd.r_s;p.r_b=pd.r_b;p.r_r=pd.r_r;p.k_l=pd.k_l;p.k_b=pd.k_b;p.k_s=pd.k_s;
  p.k_r=pd.k_r;p.a_bio=pd.a_bio;p.a_y=pd.a_y;p.a_l1=pd.a_l1;p.a_l2=pd.a_l2;p.a_f1=pd.a_f1;
  p.a_f2=pd.a_f2;p.hmat=pd.hmat;p.omega=pd.omega;p.a_f3=pd.a_f3;p.d_I=pd.d_I;
  p.a_dG1=pd.a_dG1;p.a_dG2=pd.a_dG2;
  for (std::size_t k = 0; k < traits.size(); ++k) set_kernel_trait(p, traits[k], tr[k]);
  return p;
}

// The 7-state cohort on the tape: 5 demog + collar + log_density.
template <typename S> struct St7ad { plant::FF16State<S> demog; S collar; S log_density; };

} // namespace


// Compiled core of tf24f_census_gradient_ad(). Returns {jacobian = metrics x
// traits, values}. Harvests every cohort's per-stage leaf operating point (fwd +
// backward-for-g') in a double pass, then replays the 7-state system onto ONE
// tape per cohort and takes one reverse sweep per metric. The census metrics
// couple cohorts through the height-trapezium, so the reduction is taped over all
// cohorts' final states together (one tape spanning the stand, sweeping per
// metric). FROZEN resident light (the rare-mutant / invasion census gradient).
// [[Rcpp::export]]
Rcpp::List tf24f_census_gradient_ad_impl(
    Rcpp::NumericVector pp, Rcpp::List eh_list, std::vector<double> sh,
    std::vector<int> birth, double birth_rate, double k_acclim, bool use_ad_gradient,
    std::vector<std::string> traits, std::vector<std::string> metrics,
    double trait_rel_step) {
  using std::exp; using std::log;
  for (auto& nm : metrics)
    if (nm != "LAI" && nm != "biomass" && nm != "size_moment")
      Rcpp::stop("unknown TF24f census metric: " + nm);

  plant::TF24f_Strategy s = make_tf24f(pp, k_acclim, use_ad_gradient);
  plant::TF24ProdPars<double> pd = s.prod_pars();
  const double h0v = s.initial_height(), a_d0 = s.pars.a_d0;
  const double kI = s.pars.k_I, recruitment_decay = s.pars.recruitment_decay;
  const double GEPS = 1e-6;
  const std::size_t T = traits.size(), M = metrics.size();

  // Base trait values + per-trait FD steps; build the +/- perturbed strategies once
  // (reused across every stage harvest) and the IFT injection d(h0)/d(theta_k) and
  // d(collar0)/d(theta_k) (the latter taken at the seed operating point below).
  Rcpp::CharacterVector ppn = pp.names();
  std::vector<double> tr0(T), dk(T), dh0(T, 0.0);
  std::vector<plant::TF24f_Strategy> sp; sp.reserve(T);
  std::vector<plant::TF24f_Strategy> sm; sm.reserve(T);
  for (std::size_t k = 0; k < T; ++k) {
    bool found = false;
    for (R_xlen_t j = 0; j < ppn.size(); ++j)
      if (std::string(ppn[j]) == traits[k]) { found = true; break; }
    if (!found) Rcpp::stop("unknown TF24f trait (not in strategy pars): " + traits[k]);
    tr0[k] = pp[traits[k]];
    dk[k]  = trait_rel_step * std::max(std::abs(tr0[k]), 1e-8);
    Rcpp::NumericVector q1 = Rcpp::clone(pp); q1[traits[k]] = tr0[k] + dk[k];
    Rcpp::NumericVector q2 = Rcpp::clone(pp); q2[traits[k]] = tr0[k] - dk[k];
    sp.push_back(make_tf24f(q1, k_acclim, use_ad_gradient));
    sm.push_back(make_tf24f(q2, k_acclim, use_ad_gradient));
    dh0[k] = (sp[k].initial_height() - sm[k].initial_height()) / (2.0 * dk[k]);
  }

  const std::size_t N = eh_list.size(), nC = birth.size();
  std::vector<std::vector<plant::TF24_Environment>> EH(N);
  for (std::size_t n = 0; n < N; ++n) {
    Rcpp::List st = eh_list[n];
    for (R_xlen_t k = 0; k < st.size(); ++k)
      EH[n].push_back(Rcpp::as<plant::TF24_Environment>(st[k]));
  }
  std::vector<double> step_h(N);
  for (std::size_t n = 0; n < N; ++n) step_h[n] = sh[n + 1] - sh[n];

  auto env_at = [&](std::size_t n, int stage) -> const plant::TF24_Environment& {
    return (stage == 0) ? ((n > 0) ? EH[n - 1][5] : EH[0][0]) : EH[n][stage - 1];
  };

  // --- Per-cohort harvest of the double trajectory (mirrors tf24f_census_recon_impl's
  // replay, recording the per-stage operating point + the seed). ----------------
  struct CohortH {
    std::size_t b; double collar0; LH seed;            // birth seed (h0, collar0)
    std::vector<StageH> stages;                        // per RK stage
  };
  auto harvest_cohort = [&](std::size_t i) -> CohortH {
    CohortH C; const std::size_t b = (std::size_t)birth[i]; C.b = b;
    const plant::TF24_Environment& eb = (b > 0) ? EH[b - 1][5] : EH[0][0];
    // collar0 = optimum at birth (the base optimiser run once).
    s.initializing_ = true;
    s.net_mass_production_dt(eb, h0v, s.area_leaf(h0v), 1.0 / h0v);
    C.collar0 = -s.leaf.root_collar_psi_;
    s.initializing_ = false;
    C.seed = harvest_point(s, sp, sm, dk, eb, h0v, C.collar0);

    // Re-evolve the 7-state trajectory in double, harvesting per stage. We need only
    // the (h, collar) sequence the recon visits; reuse its deriv/axpy structure.
    auto deriv = [&](const St7& y, std::size_t n, int stage) -> St7 {
      const plant::TF24_Environment& e = env_at(n, stage);
      const double h = y.demog.height, collar = y.collar;
      C.stages.push_back(StageH{
        harvest_point(s, sp, sm, dk, e, h, collar),
        harvest_point(s, sp, sm, dk, e, h - GEPS, collar)});
      double dpsi = 0.0;
      const double net = net_at_tracked(s, e, h, collar, dpsi);
      const double al = plant::tf24_area_leaf<double>(pd.a_l1, pd.a_l2, h);
      plant::TF24Rates<double> r = plant::tf24_compute_rates_from_net<double>(pd, h, al, net, true);
      const double collar_dt = k_acclim * dpsi;
      double dpsi_b = 0.0;
      const double net_b = net_at_tracked(s, e, h - GEPS, collar, dpsi_b);
      const double al_b = plant::tf24_area_leaf<double>(pd.a_l1, pd.a_l2, h - GEPS);
      const double g_back = plant::tf24_height_dt_from_net<double>(pd, h - GEPS, al_b, net_b);
      const double gprime = (r.height_dt - g_back) / GEPS;
      const double log_density_dt = -gprime - r.mortality_dt;
      return St7{plant::FF16State<double>{r.height_dt, r.mortality_dt, r.fecundity_dt,
        r.area_heartwood_dt, r.mass_heartwood_dt}, collar_dt, log_density_dt};
    };
    auto axpy = [](const St7& a, double c, const St7& k) -> St7 {
      return St7{plant::FF16State<double>{
        a.demog.height+c*k.demog.height, a.demog.mortality+c*k.demog.mortality,
        a.demog.fecundity+c*k.demog.fecundity, a.demog.area_heartwood+c*k.demog.area_heartwood,
        a.demog.mass_heartwood+c*k.demog.mass_heartwood}, a.collar+c*k.collar,
        a.log_density+c*k.log_density};
    };
    // Seed the same way the recon does (so the harvested trajectory matches).
    const double decay = std::exp(-recruitment_decay * sh[b]);
    double dpsi0 = 0.0;
    const double area0 = plant::tf24_area_leaf<double>(pd.a_l1, pd.a_l2, h0v);
    const double net0 = net_at_tracked(s, eb, h0v, C.collar0, dpsi0);
    double pr_estab = 0.0;
    if (net0 > 0.0) { const double uu = a_d0 * area0 / net0; pr_estab = decay / (uu * uu + 1.0); }
    const double mort0 = (pr_estab > 0.0) ? -log(pr_estab)
                                          : std::numeric_limits<double>::infinity();
    const double g0 = plant::tf24_height_dt_from_net<double>(pd, h0v, area0, net0);
    const double logd0 = (g0 > 0.0 && pr_estab > 0.0) ? log(birth_rate * pr_estab / g0)
                                          : -std::numeric_limits<double>::infinity();
    St7 y{plant::FF16State<double>{h0v, mort0, 0.0, 0.0, 0.0}, C.collar0, logd0};
    plant::ff16_cashkarp_replay(y, step_h, b, deriv, axpy);
    return C;
  };

  std::vector<CohortH> CH(nC);
  for (std::size_t i = 0; i < nC; ++i) CH[i] = harvest_cohort(i);

  // Pending-seed (new_node) tail: seed harvested in the FINAL-time env.
  const plant::TF24_Environment& ef = EH[N - 1][5];
  double collar0_new;
  { s.initializing_ = true; s.net_mass_production_dt(ef, h0v, s.area_leaf(h0v), 1.0 / h0v);
    collar0_new = -s.leaf.root_collar_psi_; s.initializing_ = false; }
  LH seed_new = harvest_point(s, sp, sm, dk, ef, h0v, collar0_new);
  const double decay_new = std::exp(-recruitment_decay * sh[N]);

  // --- Tape: replay the 7-state system for cohort i with active traits, returning
  // the final (height, log_density, mass_heartwood) the census reduction needs. ----
  auto replay_ad = [&](std::size_t i, const std::vector<ad_t>& tr,
                       const plant::TF24ProdPars<ad_t>& pf,
                       ad_t& H_out, ad_t& logd_out, ad_t& mhw_out) {
    const CohortH& C = CH[i];
    const std::size_t b = C.b;
    // Active h0 (IFT) and collar0 (= optimum: d(collar0)/dtheta = -d2p_dpsidth/d2p_dpsi2).
    ad_t h0 = ad_t(h0v);
    for (std::size_t k = 0; k < T; ++k) if (dh0[k] != 0.0) h0 += ad_t(dh0[k]) * (tr[k] - ad_t(tr0[k]));
    ad_t collar0 = ad_t(C.collar0);
    const double inv_c2 = (C.seed.d2p_dpsi2 != 0.0) ? 1.0 / C.seed.d2p_dpsi2 : 0.0;
    for (std::size_t k = 0; k < T; ++k) {
      const double dcollar0 = -C.seed.d2p_dpsidth[k] * inv_c2;
      if (dcollar0 != 0.0) collar0 += ad_t(dcollar0) * (tr[k] - ad_t(tr0[k]));
    }
    // Establishment + seedling growth -> initial mortality, log_density.
    ad_t area0 = plant::tf24_area_leaf<ad_t>(pf.a_l1, pf.a_l2, h0);
    ad_t profit0 = profit_lin(C.seed, h0, collar0, tr, tr0);
    ad_t net0 = plant::tf24_net_mass_production<ad_t>(pf, h0, area0, profit0);
    ad_t uu = ad_t(a_d0) * area0 / net0;
    const double decay = std::exp(-recruitment_decay * sh[b]);
    ad_t pr_estab = ad_t(decay) / (uu * uu + ad_t(1.0));
    ad_t mort0 = -log(pr_estab);
    ad_t g0 = plant::tf24_height_dt_from_net<ad_t>(pf, h0, area0, net0);
    ad_t logd0 = log(ad_t(birth_rate) * pr_estab / g0);

    std::size_t idx = 0;
    auto deriv = [&](const St7ad<ad_t>& y, std::size_t n, int stage) -> St7ad<ad_t> {
      const StageH& S = C.stages[idx++];
      ad_t h = y.demog.height, collar = y.collar;
      ad_t profit = profit_lin(S.fwd, h, collar, tr, tr0);
      ad_t al = plant::tf24_area_leaf<ad_t>(pf.a_l1, pf.a_l2, h);
      ad_t net = plant::tf24_net_mass_production<ad_t>(pf, h, al, profit);
      plant::TF24Rates<ad_t> r = plant::tf24_compute_rates_from_net<ad_t>(pf, h, al, net, true);
      ad_t collar_dt = ad_t(k_acclim) * dpsi_lin(S.fwd, h, collar, tr, tr0);
      // g' = backward FD of height_dt (collar held), the SCM's scheme: the backward
      // profit is linearised about S.back (harvested at h-GEPS), height shift -GEPS.
      ad_t profit_b = profit_lin(S.back, h - ad_t(GEPS), collar, tr, tr0);
      ad_t al_b = plant::tf24_area_leaf<ad_t>(pf.a_l1, pf.a_l2, h - ad_t(GEPS));
      ad_t net_b = plant::tf24_net_mass_production<ad_t>(pf, h - ad_t(GEPS), al_b, profit_b);
      ad_t g_back = plant::tf24_height_dt_from_net<ad_t>(pf, h - ad_t(GEPS), al_b, net_b);
      ad_t gprime = (r.height_dt - g_back) / ad_t(GEPS);
      ad_t log_density_dt = -gprime - r.mortality_dt;
      return St7ad<ad_t>{plant::FF16State<ad_t>{r.height_dt, r.mortality_dt, r.fecundity_dt,
        r.area_heartwood_dt, r.mass_heartwood_dt}, collar_dt, log_density_dt};
    };
    auto axpy = [](const St7ad<ad_t>& a, double c, const St7ad<ad_t>& k) -> St7ad<ad_t> {
      return St7ad<ad_t>{plant::FF16State<ad_t>{
        a.demog.height+c*k.demog.height, a.demog.mortality+c*k.demog.mortality,
        a.demog.fecundity+c*k.demog.fecundity, a.demog.area_heartwood+c*k.demog.area_heartwood,
        a.demog.mass_heartwood+c*k.demog.mass_heartwood}, a.collar+c*k.collar,
        a.log_density+c*k.log_density};
    };
    St7ad<ad_t> y{plant::FF16State<ad_t>{h0, mort0, ad_t(0), ad_t(0), ad_t(0)}, collar0, logd0};
    St7ad<ad_t> f = plant::ff16_cashkarp_replay(y, step_h, b, deriv, axpy);
    H_out = f.demog.height; logd_out = f.log_density; mhw_out = f.demog.mass_heartwood;
  };

  // Active pending-seed density at the final env (the trapezium tail term).
  auto dens_new_ad = [&](const std::vector<ad_t>& tr,
                         const plant::TF24ProdPars<ad_t>& pf) -> ad_t {
    ad_t h0 = ad_t(h0v);
    for (std::size_t k = 0; k < T; ++k) if (dh0[k] != 0.0) h0 += ad_t(dh0[k]) * (tr[k] - ad_t(tr0[k]));
    ad_t collar0 = ad_t(collar0_new);
    const double inv_c2 = (seed_new.d2p_dpsi2 != 0.0) ? 1.0 / seed_new.d2p_dpsi2 : 0.0;
    for (std::size_t k = 0; k < T; ++k) {
      const double dc0 = -seed_new.d2p_dpsidth[k] * inv_c2;
      if (dc0 != 0.0) collar0 += ad_t(dc0) * (tr[k] - ad_t(tr0[k]));
    }
    ad_t area0 = plant::tf24_area_leaf<ad_t>(pf.a_l1, pf.a_l2, h0);
    ad_t profit0 = profit_lin(seed_new, h0, collar0, tr, tr0);
    ad_t net0 = plant::tf24_net_mass_production<ad_t>(pf, h0, area0, profit0);
    ad_t uu = ad_t(a_d0) * area0 / net0;
    ad_t pr_estab = ad_t(decay_new) / (uu * uu + ad_t(1.0));
    ad_t g0 = plant::tf24_height_dt_from_net<ad_t>(pf, h0, area0, net0);
    return ad_t(birth_rate) * pr_estab / g0;
  };

  // ONE tape over all cohorts (census couples them via the trapezium); one reverse
  // sweep per metric. Record the forward replay, then reduce each metric.
  ad::tape_type tape;
  std::vector<ad_t> tr(T);
  for (std::size_t k = 0; k < T; ++k) tr[k] = tr0[k];
  for (auto& x : tr) tape.registerInput(x);
  tape.newRecording();
  plant::TF24ProdPars<ad_t> pf = pf_active(pd, traits, tr);

  std::vector<ad_t> Hf(nC), Lf(nC), Mf(nC);
  for (std::size_t i = 0; i < nC; ++i) replay_ad(i, tr, pf, Hf[i], Lf[i], Mf[i]);
  ad_t h0a = ad_t(h0v);
  for (std::size_t k = 0; k < T; ++k) if (dh0[k] != 0.0) h0a += ad_t(dh0[k]) * (tr[k] - ad_t(tr0[k]));
  ad_t dens_new = dens_new_ad(tr, pf);

  // Descending-height order (frozen, on the double values) for the trapezium walk.
  std::vector<std::size_t> ord(nC);
  for (std::size_t i = 0; i < nC; ++i) ord[i] = i;
  std::sort(ord.begin(), ord.end(),
            [&](std::size_t a, std::size_t b){ return as_dbl(Hf[a]) > as_dbl(Hf[b]); });

  auto census_reduce = [&](auto psi) -> ad_t {
    std::vector<ad_t> phi(nC);
    for (std::size_t i = 0; i < nC; ++i) phi[i] = psi(Hf[i], exp(Lf[i]), Mf[i]);
    ad_t J = ad_t(0.0);
    for (std::size_t j = 0; j + 1 < nC; ++j) {
      const std::size_t a = ord[j], b = ord[j + 1];
      J += ad_t(0.5) * (Hf[a] - Hf[b]) * (phi[a] + phi[b]);
    }
    if (nC > 0) {
      const std::size_t last = ord[nC - 1];
      ad_t phi_new = psi(h0a, dens_new, ad_t(0.0));
      J += ad_t(0.5) * (Hf[last] - h0a) * (phi[last] + phi_new);
    }
    return J;
  };

  std::vector<ad_t> J(M);
  for (std::size_t m = 0; m < M; ++m) {
    const std::string& nm = metrics[m];
    if (nm == "LAI") {
      J[m] = census_reduce([&](ad_t h, ad_t dens, ad_t) -> ad_t {
        return dens * ad_t(kI) * plant::tf24_area_leaf<ad_t>(pf.a_l1, pf.a_l2, h); });
    } else if (nm == "biomass") {
      J[m] = census_reduce([&](ad_t h, ad_t dens, ad_t mhw) -> ad_t {
        ad_t al = plant::tf24_area_leaf<ad_t>(pf.a_l1, pf.a_l2, h);
        ad_t mass_leaf    = al * pf.lma;
        ad_t area_sapwood = al * pf.theta;
        ad_t mass_sapwood = area_sapwood * h * pf.eta_c * pf.rho;
        ad_t area_bark    = pf.a_b1 * al * pf.theta;
        ad_t mass_bark    = area_bark * h * pf.eta_c * pf.rho;
        ad_t mass_root    = pf.a_r1 * al;
        return dens * (mass_leaf + mass_sapwood + mass_bark + mass_root + mhw); });
    } else { // size_moment
      J[m] = census_reduce([&](ad_t h, ad_t dens, ad_t) -> ad_t { return dens * h; });
    }
  }

  // Reduce each metric: register, seed, sweep (clearing derivatives between metrics).
  Rcpp::NumericMatrix jac(M, T);
  Rcpp::NumericVector values(M);
  for (std::size_t m = 0; m < M; ++m) { values[m] = as_dbl(J[m]); tape.registerOutput(J[m]); }
  for (std::size_t m = 0; m < M; ++m) {
    tape.clearDerivatives();
    xad::derivative(J[m]) = 1.0;
    tape.computeAdjoints();
    for (std::size_t k = 0; k < T; ++k) jac(m, k) = xad::derivative(tr[k]);
  }
  jac.attr("dimnames") = Rcpp::List::create(Rcpp::wrap(metrics), Rcpp::wrap(traits));
  values.attr("names") = Rcpp::wrap(metrics);
  return Rcpp::List::create(Rcpp::Named("jacobian") = jac, Rcpp::Named("values") = values);
}


// ===========================================================================
// TF24f individual grow-to-size reverse-mode AD trait gradient (#472 scope B,
// build-order step 4 -- the "individuals" surface). A single TF24f plant grown
// in a FIXED environment to target size(s), differentiated w.r.t. traits: the
// AD refine of tf24f_grow_individual_to_size_gradient_fd. The TF24f analogue of
// ff16_grow_to_size_gradient_impl, but the trajectory carries the tracked collar
// (opt_root_psi_state) as a 6th state -- the SAME curvature-linearised collar
// machinery the census tape uses -- and there is NO canopy / density / g' (the
// env is given), so it is the lightest TF24f AD surface.
//
// No resident feedback, so the only machinery beyond the replay is the stopping-
// time implicit-function step. Pass 1 (R) runs the live grow (grow_individual_bracket)
// to harvest its adaptive Cash-Karp schedule + initial state. Pass 2 (here):
// discover t* (the size component crosses target) over the FROZEN schedule reading
// the FIXED env, harvest the per-RK-stage leaf operating point at the tracked
// collar along the schedule-to-t*, replay the 6-state system as a leaf-opt-free
// tapeable expression, and take one reverse sweep per state component. The collar
// starts at the individual's birth value (0, clamped into the feasible interval on
// the first eval) -- theta-independent, so unlike the census seed it needs no
// collar-IFT injection; only the seedling height h0 carries d(h0)/d(theta).
// The stopping time responds via the IFT on size(t*, theta) = target:
//   d(t*)/d(theta) = -(d size/d theta|t*) / size_dt(t*),
// and the TOTAL derivative of each component is the partial plus rate*d(t*)/d(theta).

namespace {

// 6-state grow cohort on the tape: 5 demog + tracked collar.
template <typename S> struct St6 { plant::FF16State<S> demog; S collar; };

template <typename S>
St6<S> st6_axpy(const St6<S>& a, double c, const St6<S>& k) {
  return St6<S>{plant::FF16State<S>{
    a.demog.height+c*k.demog.height, a.demog.mortality+c*k.demog.mortality,
    a.demog.fecundity+c*k.demog.fecundity, a.demog.area_heartwood+c*k.demog.area_heartwood,
    a.demog.mass_heartwood+c*k.demog.mass_heartwood}, a.collar+c*k.collar};
}
double st6_at(const St6<double>& s, int c) {
  switch (c) { case 0: return s.demog.height; case 1: return s.demog.mortality;
    case 2: return s.demog.fecundity; case 3: return s.demog.area_heartwood;
    case 4: return s.demog.mass_heartwood; default: return s.collar; }
}

} // namespace


// Compiled core of tf24f_grow_individual_to_size_gradient(). Returns, per target
// size: the stopping time t*, the ODE state at t*, d(t*)/d(theta) and the TOTAL
// d(state at t*)/d(theta) (a sizes x component x trait array). y0v is the live
// individual's initial ODE state {height, mortality, fecundity, area_heartwood,
// mass_heartwood, opt_root_psi_state}; sh the harvested adaptive step-time schedule
// (grow_individual_bracket$time); sidx the size component (0 = height). shading /
// gss_tol / k_acclim / use_ad_gradient mirror the live individual's strategy.
// [[Rcpp::export]]
Rcpp::List tf24f_grow_to_size_gradient_impl(
    Rcpp::NumericVector pp, plant::TF24_Environment env, Rcpp::NumericVector y0v,
    std::vector<double> sh, std::vector<double> targets, int sidx,
    std::vector<std::string> traits, double k_acclim, bool use_ad_gradient,
    std::string shading, double gss_tol, double trait_rel_step) {
  using std::exp; using std::log;
  if (sh.size() < 2) Rcpp::stop("schedule must have at least two step times");
  if (y0v.size() < 6) Rcpp::stop("y0 must have the 6 TF24f ODE states");

  plant::TF24f_Strategy s = make_tf24f(pp, k_acclim, use_ad_gradient, shading, gss_tol);
  plant::TF24ProdPars<double> pd = s.prod_pars();
  const std::size_t T = traits.size();

  // Base trait values + per-trait FD steps; build the +/- perturbed strategies once
  // (reused across targets/stages) and the seedling-size IFT d(h0)/d(theta).
  Rcpp::CharacterVector ppn = pp.names();
  std::vector<double> tr0(T), dk(T), dh0(T, 0.0);
  std::vector<plant::TF24f_Strategy> sp; sp.reserve(T);
  std::vector<plant::TF24f_Strategy> sm; sm.reserve(T);
  const double h0v = y0v[0];
  for (std::size_t k = 0; k < T; ++k) {
    bool found = false;
    for (R_xlen_t j = 0; j < ppn.size(); ++j)
      if (std::string(ppn[j]) == traits[k]) { found = true; break; }
    if (!found) Rcpp::stop("unknown TF24f trait (not in strategy pars): " + traits[k]);
    tr0[k] = pp[traits[k]];
    dk[k]  = trait_rel_step * std::max(std::abs(tr0[k]), 1e-8);
    Rcpp::NumericVector q1 = Rcpp::clone(pp); q1[traits[k]] = tr0[k] + dk[k];
    Rcpp::NumericVector q2 = Rcpp::clone(pp); q2[traits[k]] = tr0[k] - dk[k];
    sp.push_back(make_tf24f(q1, k_acclim, use_ad_gradient, shading, gss_tol));
    sm.push_back(make_tf24f(q2, k_acclim, use_ad_gradient, shading, gss_tol));
    dh0[k] = (sp[k].initial_height() - sm[k].initial_height()) / (2.0 * dk[k]);
  }

  std::vector<double> step_h(sh.size() - 1);
  for (std::size_t n = 0; n + 1 < sh.size(); ++n) step_h[n] = sh[n + 1] - sh[n];

  // Double 6-state derivative reading the fixed env (real leaf at tracked collar).
  auto deriv_d = [&](const St6<double>& y, std::size_t, int) -> St6<double> {
    double dpsi = 0.0;
    const double net = net_at_tracked(s, env, y.demog.height, y.collar, dpsi);
    const double al = plant::tf24_area_leaf<double>(pd.a_l1, pd.a_l2, y.demog.height);
    plant::TF24Rates<double> r = plant::tf24_compute_rates_from_net<double>(pd, y.demog.height, al, net, true);
    return St6<double>{plant::FF16State<double>{r.height_dt, r.mortality_dt, r.fecundity_dt,
      r.area_heartwood_dt, r.mass_heartwood_dt}, k_acclim * dpsi};
  };
  auto axpy_d = [](const St6<double>& a, double c, const St6<double>& k){ return st6_axpy<double>(a,c,k); };
  auto replay_one_d = [&](St6<double> y, const std::vector<double>& sched) -> St6<double> {
    return plant::ff16_cashkarp_replay(y, sched, (std::size_t)0, deriv_d, axpy_d);
  };

  St6<double> y0{plant::FF16State<double>{y0v[0], y0v[1], y0v[2], y0v[3], y0v[4]}, y0v[5]};

  const std::vector<std::string> comp =
    {"height","mortality","fecundity","area_heartwood","mass_heartwood","opt_root_psi_state"};
  const std::size_t nS = comp.size(), nG = targets.size();
  Rcpp::NumericVector tstar(nG);
  Rcpp::NumericMatrix state(nG, nS);
  Rcpp::NumericMatrix dtime(nG, T);
  Rcpp::NumericVector dstate(nG * nS * T);
  auto DS = [&](std::size_t g, std::size_t c, std::size_t k) -> double& {
    return dstate[g + nG * (c + nS * k)]; };

  for (std::size_t g = 0; g < nG; ++g) {
    // --- Discover t*: walk the schedule one RKCK step at a time until the size
    // component crosses target, then bisect the partial final step. -------------
    int nfull = -1; double dt_final = 0.0; St6<double> fin{};
    {
      St6<double> y = y0;
      for (std::size_t n = 0; n < step_h.size(); ++n) {
        St6<double> ynext = replay_one_d(y, {step_h[n]});
        if (st6_at(ynext, sidx) >= targets[g]) {
          nfull = (int)n;
          double lo = 0.0, hi = step_h[n];
          for (int it = 0; it < 80; ++it) {
            double mid = 0.5 * (lo + hi);
            double sm_ = st6_at(replay_one_d(y, {mid}), sidx);
            if (sm_ < targets[g]) lo = mid; else hi = mid;
          }
          dt_final = 0.5 * (lo + hi);
          fin = replay_one_d(y, {dt_final});
          break;
        }
        y = ynext;
      }
    }
    if (nfull < 0) Rcpp::stop("target size not reached within the schedule");
    tstar[g] = sh[nfull] + dt_final;
    for (std::size_t c = 0; c < nS; ++c) state(g, c) = st6_at(fin, (int)c);

    // Rates at t* (double) for the IFT correction.
    St6<double> rate = deriv_d(fin, 0, 0);
    double yd[6] = {rate.demog.height, rate.demog.mortality, rate.demog.fecundity,
                    rate.demog.area_heartwood, rate.demog.mass_heartwood, rate.collar};
    const double sdot = yd[sidx];

    // Frozen schedule to t*: nfull full steps + the single partial step dt_final.
    std::vector<double> sched(step_h.begin(), step_h.begin() + nfull);
    sched.push_back(dt_final);

    // --- Harvest the per-RK-stage leaf operating point along sched (double). ----
    std::vector<LH> stages;
    {
      auto deriv_h = [&](const St6<double>& y, std::size_t, int) -> St6<double> {
        stages.push_back(harvest_point(s, sp, sm, dk, env, y.demog.height, y.collar));
        return deriv_d(y, 0, 0);
      };
      St6<double> y = y0;
      plant::ff16_cashkarp_replay(y, sched, (std::size_t)0, deriv_h, axpy_d);
    }

    // --- Tape replay over sched: partial d(state)/d(theta) holding t* fixed. -----
    ad::tape_type tape;
    std::vector<ad_t> tr(T);
    for (std::size_t k = 0; k < T; ++k) tr[k] = tr0[k];
    for (auto& x : tr) tape.registerInput(x);
    tape.newRecording();
    plant::TF24ProdPars<ad_t> pf = pf_active(pd, traits, tr);
    ad_t h0 = ad_t(h0v);
    for (std::size_t k = 0; k < T; ++k) if (dh0[k] != 0.0) h0 += ad_t(dh0[k]) * (tr[k] - ad_t(tr0[k]));

    std::size_t idx = 0;
    auto deriv_ad = [&](const St6<ad_t>& y, std::size_t, int) -> St6<ad_t> {
      const LH& H = stages[idx++];
      ad_t h = y.demog.height, collar = y.collar;
      ad_t profit = profit_lin(H, h, collar, tr, tr0);
      ad_t al = plant::tf24_area_leaf<ad_t>(pf.a_l1, pf.a_l2, h);
      ad_t net = plant::tf24_net_mass_production<ad_t>(pf, h, al, profit);
      plant::TF24Rates<ad_t> r = plant::tf24_compute_rates_from_net<ad_t>(pf, h, al, net, true);
      ad_t collar_dt = ad_t(k_acclim) * dpsi_lin(H, h, collar, tr, tr0);
      return St6<ad_t>{plant::FF16State<ad_t>{r.height_dt, r.mortality_dt, r.fecundity_dt,
        r.area_heartwood_dt, r.mass_heartwood_dt}, collar_dt};
    };
    auto axpy_ad = [](const St6<ad_t>& a, double c, const St6<ad_t>& k){ return st6_axpy<ad_t>(a,c,k); };
    St6<ad_t> yad{plant::FF16State<ad_t>{h0, ad_t(y0v[1]), ad_t(y0v[2]), ad_t(y0v[3]),
      ad_t(y0v[4])}, ad_t(y0v[5])};
    St6<ad_t> out = plant::ff16_cashkarp_replay(yad, sched, (std::size_t)0, deriv_ad, axpy_ad);
    ad_t oc[6] = {out.demog.height, out.demog.mortality, out.demog.fecundity,
                  out.demog.area_heartwood, out.demog.mass_heartwood, out.collar};
    for (std::size_t c = 0; c < nS; ++c) tape.registerOutput(oc[c]);
    std::vector<std::vector<double>> P(nS, std::vector<double>(T, 0.0));
    for (std::size_t c = 0; c < nS; ++c) {
      tape.clearDerivatives();
      xad::derivative(oc[c]) = 1.0;
      tape.computeAdjoints();
      for (std::size_t k = 0; k < T; ++k) P[c][k] = xad::derivative(tr[k]);
    }
    // IFT: d(t*)/d(theta_k) = -(d size/d theta_k|t*) / size_dt(t*); total = P + ydot*dt*.
    for (std::size_t k = 0; k < T; ++k) {
      double dtk = (sdot != 0.0) ? -P[sidx][k] / sdot : 0.0;
      dtime(g, k) = dtk;
      for (std::size_t c = 0; c < nS; ++c) DS(g, c, k) = P[c][k] + yd[c] * dtk;
    }
  }

  state.attr("dimnames")  = Rcpp::List::create(R_NilValue, Rcpp::wrap(comp));
  dtime.attr("dimnames")  = Rcpp::List::create(R_NilValue, Rcpp::wrap(traits));
  dstate.attr("dim")      = Rcpp::IntegerVector::create((int)nG, (int)nS, (int)T);
  dstate.attr("dimnames") = Rcpp::List::create(R_NilValue, Rcpp::wrap(comp), Rcpp::wrap(traits));
  return Rcpp::List::create(Rcpp::Named("time") = tstar, Rcpp::Named("state") = state,
                            Rcpp::Named("d_time") = dtime, Rcpp::Named("d_state") = dstate);
}
