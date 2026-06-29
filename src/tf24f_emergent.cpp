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
#include <plant.h>                 // RcppR6 as<>/wrap for TF24_Environment
#include <plant/models/tf24f_strategy.h>
#include <plant/models/tf24_production_kernel.h>
#include <plant/models/ff16_production_kernel.h>   // ff16_cashkarp_replay, FF16State

namespace {

// Build a TF24f strategy from the SCM's parameter vector. Mirrors the TF24
// make_strategy in tf24_emergent.cpp but constructs the fast-acclimation variant and
// sets its acclimation knobs. crown-centre + GSS_tol_abs match the resident run.
plant::TF24f_Strategy make_tf24f(const Rcpp::NumericVector& pp, double k_acclim,
                                 bool use_ad_gradient) {
  plant::TF24f_Strategy s;
  s.control.shading_model = "crown-centre"; s.control.GSS_tol_abs = 1e-9;
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
