// Gates for Leaf::input_adjoints. Standalone: builds a Leaf, solves, and runs
// the three exact invariants plus the checks a finite difference can carry.
#include <plant/leaf_model.h>
#include <cstdio>
#include <cmath>
#include <vector>

using plant::Leaf;

static Leaf make_leaf() {
  double vcmax_25 = 100, jmax_25 = vcmax_25 * 1.67;
  double c = 2.04, b = 3.0, psi_crit = 5.0, beta2 = 1.0;
  double root_c = 2.65, root_b = 1.29;
  double root_psi_crit = root_b * std::pow(std::log(1.0 / 0.05), 1.0 / root_c);
  Leaf l(vcmax_25, c, b, psi_crit, root_c, root_b, root_psi_crit, beta2,
         jmax_25, 0.3, 0.7, 0.99, 1e-3, 100, 1e-6, 1000, 46.32995, 3.4e3, 9.4e4);
  l.initialize_integrator();
  l.setup_transpiration(100);
  l.setup_root_vulnerability(100);
  return l;
}

static void seat(Leaf& l, std::vector<double> psi_soil, double PPFD,
                 double area_leaf) {
  int n = psi_soil.size();
  std::vector<double> depth(n), mass(n);
  for (int i = 0; i < n; ++i) { depth[i] = (i + 1.0) / n; mass[i] = 1.0 / n; }
  double theta = 0.000157, h = 5.0;
  l.set_physiology(area_leaf, mass, 608, 0.0245, PPFD, psi_soil, depth,
                   1.0 * theta / h, 2.0, 40.0, theta * h, 25.0, 21.0, 101.3);
}

static double rel(double x, double y) {
  double s = std::max(std::abs(x), std::abs(y));
  return s == 0.0 ? 0.0 : std::abs(x - y) / s;
}

int main() {
  std::vector<std::vector<double>> states = {
      {0.10, 0.30, 0.60, 0.90, 1.20},
      {0.02, 0.05, 0.08, 0.12, 0.17},
      {0.50, 0.80, 1.10, 1.40, 1.70, 2.00, 2.30, 2.60},
      {1.60, 2.20, 2.80, 3.40, 4.00}};
  for (size_t st = 0; st < states.size(); ++st) {
    Leaf l = make_leaf();
    seat(l, states[st], 900.0, 0.05);
    l.find_root_collar_psi();
    const int n = l.max_soil_layer;
    const double p = -l.root_collar_psi_;
    printf("\n=== state %zu: n=%d p*=%.9g pinned=%d Pi_pp=%.6g\n",
           st, n, p, (int)l.collar_pinned_, l.dR_dcollar_at(p, 1e-6));
    printf("    |R| at p* = %.4e   profit_=%.9g\n",
           std::abs(l.dprofit_droot_collar_psi(p)), l.profit_);

    // dR_dcollar, and its own error under a halved step
    double d1 = l.dR_dcollar_at(p, 1e-6), d2 = l.dR_dcollar_at(p, 5e-7);
    printf("    Pi_pp = %.9g (h=1e-6) vs %.9g (h=5e-7), rel move %.3e\n",
           d1, d2, rel(d1, d2));

    // the leaf-side coefficients the waist needs
    l.refresh_soil_potentials();
    l.E_from_Soil_to_Root_Collar(l.root_collar_psi_, l.psi_soil_inverted_);
    const double E_up = l.E_up_;
    const double kappa = l.leaf_specific_conductance_max_;
    const double S_p = l.transpiration_from_psi.eval(p);
    const double E_ps = E_up / kappa + S_p;
    const double P_prime = l.psi_from_transpiration.deriv(E_ps);

    // CONTINUITY: E_up from the soil side against kappa*(S(psi_stem) - S(p)),
    // and the same identity differentiated (the interpolant-derivative check).
    const double psi_stem = l.opt_psi_stem_;
    const double E_stem = kappa * (l.transpiration_from_psi.eval(psi_stem) - S_p);
    const double S_prime_stem = l.transpiration_from_psi.deriv(psi_stem);
    printf("    CONTINUITY value  rel %.3e   (E_up=%.9g, stem side=%.9g)\n",
           rel(E_up, E_stem), E_up, E_stem);
    printf("    CONTINUITY deriv  S'(psi_stem)*P'(E_ps) - 1 = %.3e\n",
           S_prime_stem * P_prime - 1.0);

    Leaf::LayerFlux fx;
    l.layer_flux_partials(l.root_collar_psi_, l.psi_soil_inverted_, fx);

    // the two waist scalars
    const double ci = l.psi_stem_to_ci(psi_stem, p);
    // rebuild dprofit_dpsistem the way input_adjoints does, via R's own pieces
    const double dEup_dr = l.dE_from_soil_dpsi_collar(l.root_collar_psi_,
                                                      l.psi_soil_inverted_);
    double sum_dEdr = 0.0;
    for (int i = 0; i < n; ++i) sum_dEdr += fx.dE_dr[i];
    printf("    dE_up/dr: analytic sum %.9g vs dE_from_soil %.9g rel %.3e\n",
           sum_dEdr * plant::kg_per_mol_h2o, dEup_dr,
           rel(sum_dEdr * plant::kg_per_mol_h2o, dEup_dr));
    (void)ci;

    // b from its closed form, recovered by a residual pair in the slope
    // direction is not available, so it is taken from input_adjoints' own
    // expression and then tested inside the waist residual below.
    std::vector<double> lam_up(n, 0.0), adj;
    l.input_adjoints(1.0, lam_up, adj);
    std::vector<std::string> nm = l.inputs();
    printf("    inputs()=%zu adjoints=%zu\n", nm.size(), adj.size());

    // dR_dflux recovered from each layer in turn
    // recompute b exactly as input_adjoints does
    {
      const double gstar = l.gamma_ * plant::umol_per_mol_to_Pa;
      const double gc_const = l.atm_kpa_ * plant::kg_to_mol_h2o / l.atm_vpd_ /
                              plant::H2O_CO2_stom_diff_ratio;
      const double gc = gc_const * l.transpiration(psi_stem, p);
      const double inv_atm = 1.0 / (l.atm_kpa_ * plant::kPa_to_Pa);
      const double A_prime = (l.assim_colimited(ci + 1e-7) -
                              l.assim_colimited(ci - 1e-7)) / 2e-7;
      const double C_prime = (l.hydraulic_cost_TF(psi_stem + 1e-7) -
                              l.hydraulic_cost_TF(psi_stem - 1e-7)) / 2e-7;
      const double g_ci = A_prime * plant::umol_to_mol + gc * inv_atm;
      const double dgc_dps = gc_const * kappa * S_prime_stem;
      const double dci_dps = dgc_dps * (l.ca_ - ci) * inv_atm / g_ci;
      const double bcoef = -(A_prime * dci_dps - C_prime) * P_prime / kappa;
      (void)gstar;
      printf("    dR_dflux_slope (closed form) = %.9g\n", bcoef);
      double amin = 1e300, amax = -1e300;
      for (int j = 0; j < n; ++j) {
        double aj = l.dR_dflux_from_layer(j, 1e-6, bcoef);
        printf("      dR_dflux from layer %d = %.9g\n", j, aj);
        amin = std::min(amin, aj); amax = std::max(amax, aj);
      }
      printf("    dR_dflux spread over %d layers = %.3e relative\n", n,
             (amax - amin) / std::abs(0.5 * (amax + amin)));
      const double a_flux = l.dR_dflux_from_layer(0, 1e-6, bcoef);

      // WAIST: predicted dR/du against a residual pair, over 2n+1 directions
      double worst = 0.0;
      for (int j = 0; j < n; ++j) {
        const double dEup_du = -plant::kg_per_mol_h2o * fx.dE_dpsi[j];
        const double dslope_du = -plant::kg_per_mol_h2o * fx.d2E_dr_dpsi[j];
        const double pred = a_flux * dEup_du + bcoef * dslope_du;
        const double h = 1e-6, keep = l.psi_soil_[j];
        l.psi_soil_[j] = keep + h; l.transpiration_cached_ = false;
        double Rp = l.dprofit_droot_collar_psi(p);
        l.psi_soil_[j] = keep - h; l.transpiration_cached_ = false;
        double Rm = l.dprofit_droot_collar_psi(p);
        l.psi_soil_[j] = keep; l.transpiration_cached_ = false;
        l.refresh_soil_potentials();
        double meas = (Rp - Rm) / (2 * h);
        worst = std::max(worst, rel(pred, meas));
        if (j < 2) printf("      psi[%d]: pred %.6g meas %.6g rel %.2e\n", j,
                          pred, meas, rel(pred, meas));
      }
      // leaf area: both E_up and its collar slope carry a single 1/area factor
      {
        const double pred = a_flux * (-E_up / l.area_leaf_) +
                            bcoef * (-dEup_dr / l.area_leaf_);
        const double h = l.area_leaf_ * 1e-6, keep = l.area_leaf_;
        std::vector<double> mass(n, 1.0 / n), depth(n);
        for (int i = 0; i < n; ++i) depth[i] = (i + 1.0) / n;
        double theta = 0.000157, hh = 5.0;
        double Rpm[2];
        for (int k = 0; k < 2; ++k) {
          l.set_physiology(keep + (k ? -h : h), mass, 608, 0.0245, l.PPFD_,
                           l.psi_soil_, depth, 1.0 * theta / hh, 2.0, 40.0,
                           theta * hh, 25.0, 21.0, 101.3);
          Rpm[k] = l.dprofit_droot_collar_psi(p);
        }
        l.set_physiology(keep, mass, 608, 0.0245, l.PPFD_, l.psi_soil_, depth,
                         1.0 * theta / hh, 2.0, 40.0, theta * hh, 25.0, 21.0,
                         101.3);
        double meas = (Rpm[0] - Rpm[1]) / (2 * h);
        worst = std::max(worst, rel(pred, meas));
        printf("      area:   pred %.6g meas %.6g rel %.2e\n", pred, meas,
               rel(pred, meas));
      }
      // root masses
      {
        std::vector<double> depth(n);
        for (int i = 0; i < n; ++i) depth[i] = (i + 1.0) / n;
        double theta = 0.000157, hh = 5.0;
        for (int j = 0; j < n; ++j) {
          double dEup_dm = 0.0, dslope_dm = 0.0;
          const double m_j = 3.0 * l.c_r_V_[j];
          for (int i = j; i < n; ++i) {
            const double m_i = 3.0 * l.c_r_V_[i];
            const double r_R = fx.r_R[i];
            const double q = -l.r_R_V[j] / m_j +
                             (i == j ? -(r_R - l.r_R_V_sum[i]) / m_i : 0.0);
            const double d2 = (i == j) ? -fx.dr_R_dr[i] / m_i : 0.0;
            const double N = fx.dE_dr[i] * r_R * r_R;
            const double dN = (-1.0 / l.area_leaf_) * q - fx.num[i] * d2;
            dEup_dm += -fx.num[i] * q / (r_R * r_R);
            dslope_dm += (dN * r_R - 2.0 * N * q) / (r_R * r_R * r_R);
          }
          dEup_dm *= plant::kg_per_mol_h2o;
          dslope_dm *= plant::kg_per_mol_h2o;
          const double pred = a_flux * dEup_dm + bcoef * dslope_dm;
          const double hm = m_j * 1e-6;
          double Rpm[2];
          for (int k = 0; k < 2; ++k) {
            std::vector<double> mass(n, 1.0 / n);
            mass[j] += (k ? -hm : hm);
            l.set_physiology(l.area_leaf_, mass, 608, 0.0245, l.PPFD_,
                             l.psi_soil_, depth, 1.0 * theta / hh, 2.0, 40.0,
                             theta * hh, 25.0, 21.0, 101.3);
            Rpm[k] = l.dprofit_droot_collar_psi(p);
          }
          std::vector<double> mass(n, 1.0 / n);
          l.set_physiology(l.area_leaf_, mass, 608, 0.0245, l.PPFD_, l.psi_soil_,
                           depth, 1.0 * theta / hh, 2.0, 40.0, theta * hh, 25.0,
                           21.0, 101.3);
          double meas = (Rpm[0] - Rpm[1]) / (2 * hm);
          worst = std::max(worst, rel(pred, meas));
          if (j < 2) printf("      mass[%d]: pred %.6g meas %.6g rel %.2e\n", j,
                            pred, meas, rel(pred, meas));
        }
      }
      printf("    WAIST worst relative residual over %d directions = %.3e\n",
             2 * n + 1, worst);

      // The five inputs report 02 6.8 omits. root_b/root_c rebuild the root
      // vulnerability spline; beta_R_H/beta_R_V rescale the resistance network.
      // Each reaches R only through the soil->collar transport, so the waist
      // pair should predict its dR if it factors. Measured at the fixed
      // operating point (r, p), E_up and dE_up/dr by central difference.
      std::vector<double> depth2(n);
      for (int i = 0; i < n; ++i) depth2[i] = (i + 1.0) / n;
      std::vector<double> mass2(n, 1.0 / n);
      const double theta2 = 0.000157, hh2 = 5.0;
      auto factor_test = [&](const char* nm, double& member, double h,
                             bool rebuild_root) {
        const double keep = member;
        double E[2], slope[2], R[2];
        for (int k = 0; k < 2; ++k) {
          member = keep + (k ? -h : h);
          if (rebuild_root) l.setup_root_vulnerability(100);
          l.set_physiology(l.area_leaf_, mass2, 608, 0.0245, l.PPFD_, l.psi_soil_,
                           depth2, 1.0 * theta2 / hh2, 2.0, 40.0, theta2 * hh2,
                           25.0, 21.0, 101.3);
          l.transpiration_cached_ = false;
          R[k] = l.dprofit_droot_collar_psi(p);
          l.refresh_soil_potentials();
          l.E_from_Soil_to_Root_Collar(l.root_collar_psi_, l.psi_soil_inverted_);
          E[k] = l.E_up_;
          slope[k] = l.dE_from_soil_dpsi_collar(l.root_collar_psi_,
                                                l.psi_soil_inverted_);
        }
        member = keep;
        if (rebuild_root) l.setup_root_vulnerability(100);
        l.set_physiology(l.area_leaf_, mass2, 608, 0.0245, l.PPFD_, l.psi_soil_,
                         depth2, 1.0 * theta2 / hh2, 2.0, 40.0, theta2 * hh2,
                         25.0, 21.0, 101.3);
        const double dE = (E[0] - E[1]) / (2 * h);
        const double dslope = (slope[0] - slope[1]) / (2 * h);
        const double dR = (R[0] - R[1]) / (2 * h);
        const double pred = a_flux * dE + bcoef * dslope;
        printf("      %-10s dR pred %.6g meas %.6g  rel %.2e\n", nm, pred, dR,
               rel(pred, dR));
      };
      printf("    WAIST-EXT (the five omitted inputs):\n");
      factor_test("root_b", l.root_b, l.root_b * 1e-6, true);
      factor_test("root_c", l.root_c, l.root_c * 1e-6, true);
      factor_test("beta_R_H", l.beta_R_H, l.beta_R_H * 1e-6, false);
      factor_test("beta_R_V", l.beta_R_V, l.beta_R_V * 1e-6, false);
    }

    // TRANSLATION rows, against a difference along the translation
    {
      std::vector<double> dE_dd;
      double dpsistem_dd;
      l.find_root_collar_psi();
      l.translation_partials(dE_dd, dpsistem_dd);
      const double d = 1e-5;
      std::vector<double> keep = l.psi_soil_;
      std::vector<double> Ep(n), Em(n);
      double stem_p, stem_m;
      for (int k = 0; k < 2; ++k) {
        for (int i = 0; i < n; ++i) l.psi_soil_[i] = keep[i] + (k ? -d : d);
        l.refresh_soil_potentials();
        double rr = -(p + (k ? -d : d));
        l.E_from_Soil_to_Root_Collar(rr, l.psi_soil_inverted_);
        for (int i = 0; i < n; ++i) (k ? Em : Ep)[i] = l.soil_consumption_[i];
        double s = l.find_psi_stem_from_psi_root(rr, l.psi_soil_inverted_);
        (k ? stem_m : stem_p) = s;
      }
      l.psi_soil_ = keep; l.refresh_soil_potentials();
      double worst = 0.0;
      for (int i = 0; i < n; ++i) {
        double meas = (Ep[i] - Em[i]) / (2 * d);
        worst = std::max(worst, rel(dE_dd[i], meas));
        if (i < 3) printf("      dE[%d]/dd closed %.6g diff %.6g rel %.2e\n", i,
                          dE_dd[i], meas, rel(dE_dd[i], meas));
      }
      double meas_stem = (stem_p - stem_m) / (2 * d);
      printf("    TRANSLATION uptake worst rel %.3e; stem closed %.9g diff "
             "%.9g rel %.3e\n", worst, dpsistem_dd, meas_stem,
             rel(dpsistem_dd, meas_stem));
    }

    // STATIONARITY: dR/du + Pi_pp * dp*/du = 0, with dp*/du from a re-solve in
    // one layer (well resolved: a single layer is not the near-null direction)
    {
      l.find_root_collar_psi();
      const double bcoef_dummy = 0.0; (void)bcoef_dummy;
      for (int j = 0; j < n && j < 3; ++j) {
        const double h = 1e-5, keep = l.psi_soil_[j];
        double pstar[2];
        for (int k = 0; k < 2; ++k) {
          l.psi_soil_[j] = keep + (k ? -h : h);
          l.transpiration_cached_ = false;
          l.find_root_collar_psi();
          pstar[k] = -l.root_collar_psi_;
        }
        l.psi_soil_[j] = keep; l.transpiration_cached_ = false;
        l.find_root_collar_psi();
        const double dp_du = (pstar[0] - pstar[1]) / (2 * h);
        const double hR = 1e-6;
        l.psi_soil_[j] = keep + hR; l.transpiration_cached_ = false;
        double Rp = l.dprofit_droot_collar_psi(p);
        l.psi_soil_[j] = keep - hR; l.transpiration_cached_ = false;
        double Rm = l.dprofit_droot_collar_psi(p);
        l.psi_soil_[j] = keep; l.transpiration_cached_ = false;
        l.find_root_collar_psi();
        const double dR_du = (Rp - Rm) / (2 * hR);
        const double resid = dR_du + l.dR_dcollar_at(p, 1e-6) * dp_du;
        printf("    STATIONARITY psi[%d]: dR/du %.6g  Pi_pp*dp/du %.6g  "
               "residual %.3e (rel %.2e)\n", j, dR_du,
               l.dR_dcollar_at(p, 1e-6) * dp_du, resid,
               std::abs(resid) / std::abs(dR_du));
      }
    }

    // VACUITY: uptake's whole dependence on radiation is the argmax channel
    {
      l.find_root_collar_psi();
      std::vector<double> lam(n, 0.0); lam[0] = 1.0;
      std::vector<double> adj2;
      l.input_adjoints(0.0, lam, adj2);
      printf("    VACUITY d(uptake_0)/d(PPFD) = %.6g   d/d(psi_soil[0]) = %.6g\n",
             adj2[0], adj2[1]);
    }

    // A whole-solve difference, used only in the single-layer and radiation
    // directions where it resolves the response; not in the uniform direction.
    {
      std::vector<double> mass(n, 1.0 / n), depth(n);
      for (int i = 0; i < n; ++i) depth[i] = (i + 1.0) / n;
      const double theta = 0.000157, hh = 5.0;
      auto reseat = [&]() {
        l.transpiration_cached_ = false;
        l.set_physiology(0.05, mass, 608, 0.0245, l.PPFD_, l.psi_soil_, depth,
                         1.0 * theta / hh, 2.0, 40.0, theta * hh, 25.0, 21.0,
                         101.3);
        l.find_root_collar_psi();
      };
      reseat();
      std::vector<double> lamP(n, 0.0), adjP, lam0(n, 0.0), adj0;
      l.input_adjoints(1.0, lamP, adjP);
      lam0[0] = 1.0;
      l.input_adjoints(0.0, lam0, adj0);
      double* tgt[3] = {&l.PPFD_, &l.psi_soil_[0], &l.psi_soil_[1]};
      const char* nm2[3] = {"PPFD", "psi_soil[0]", "psi_soil[1]"};
      const double hs[3] = {l.PPFD_ * 1e-5, 1e-5, 1e-5};
      for (int d = 0; d < 3; ++d) {
        const double keep = *tgt[d];
        double prof[2], up0[2];
        for (int k = 0; k < 2; ++k) {
          *tgt[d] = keep + (k ? -hs[d] : hs[d]);
          reseat();
          prof[k] = l.profit_;
          up0[k] = l.soil_consumption_[0];
        }
        *tgt[d] = keep;
        reseat();
        const double dprof = (prof[0] - prof[1]) / (2 * hs[d]);
        const double dup = (up0[0] - up0[1]) / (2 * hs[d]);
        printf("    FULLSOLVE %-12s profit %.8g vs %.8g rel %.2e | uptake0 "
               "%.6g vs %.6g rel %.2e\n",
               nm2[d], adjP[d + (d ? 0 : 0)], dprof, rel(adjP[d], dprof),
               adj0[d], dup, rel(adj0[d], dup));
      }
    }

    // FINITE: no row of the supplied Jacobian may be NA, at either operating
    // point case. This is the gate the parameter and bound rows exist to turn
    // green, so it reads red on a tree where they are unfilled.
    {
      l.find_root_collar_psi();
      std::vector<double> lam(n, 1.0), adj3;
      l.input_adjoints(1.0, lam, adj3);
      std::vector<std::string> names = l.inputs();
      int bad = 0;
      for (size_t i = 0; i < adj3.size(); ++i) {
        if (!std::isfinite(adj3[i])) {
          if (bad < 6) printf("      non-finite: %s\n", names[i].c_str());
          ++bad;
        }
      }
      printf("    FINITE pinned=%d  non-finite rows %d of %zu\n",
             (int)l.collar_pinned_, bad, adj3.size());
    }

    // The interpolant rebuilders must put the splines back bit-for-bit, or
    // input_adjoints leaves the forward model on a different leaf than it found.
    {
      std::vector<double> xs, ys, xr, yr;
      l.build_cumulative_vulnerability_integral(l.b, l.c,
                                                l.vulnerability_curve_ncontrol,
                                                xs, ys);
      l.build_cumulative_vulnerability_integral(l.root_b, l.root_c,
                                                l.vulnerability_curve_ncontrol,
                                                xr, yr);
      const double probe = 0.5 * l.psi_crit;
      const double s0 = l.transpiration_from_psi.eval(probe);
      const double p0 = l.psi_from_transpiration.eval(s0);
      const double r0 = l.root_vuln_integral_from_psi.eval(probe);
      l.set_transpiration_at(l.b * 1.1, l.c, xs);
      l.set_root_vulnerability_at(l.root_b * 1.1, l.root_c, xr);
      l.set_transpiration_at(l.b, l.c, xs);
      l.set_root_vulnerability_at(l.root_b, l.root_c, xr);
      printf("    RESTORE bit-identical: S %d  P %d  G %d\n",
             (int)(l.transpiration_from_psi.eval(probe) == s0),
             (int)(l.psi_from_transpiration.eval(s0) == p0),
             (int)(l.root_vuln_integral_from_psi.eval(probe) == r0));
    }

    // The parameter rows against a central difference of the whole solve. Valid
    // for profit by the envelope theorem; for uptake it is an order-of-magnitude
    // check only, and for the four parameters that rebuild an interpolant it is
    // no reference at all: the builder's knot COUNT steps with b and root_b, so
    // the model's own output jumps there. Both the moving-grid difference and
    // the held-grid one are reported so the size of that is visible.
    {
      static double Leaf::*const slot[] = {
          &Leaf::vcmax_25, &Leaf::jmax_25, &Leaf::a, &Leaf::curv_fact_elec_trans,
          &Leaf::curv_fact_colim, &Leaf::b, &Leaf::c, &Leaf::psi_crit,
          &Leaf::beta2, &Leaf::g1_TF24, &Leaf::rho_, &Leaf::a_bio_,
          &Leaf::root_b, &Leaf::root_c, &Leaf::root_psi_crit};
      static const char* const pname[] = {
          "vcmax_25", "jmax_25", "a", "curv_fact_elec_trans", "curv_fact_colim",
          "b", "c", "psi_crit", "beta2", "g1_TF24", "rho", "a_bio",
          "root_b", "root_c", "root_psi_crit"};
      const int i_par0 = 3 + 2 * n;
      std::vector<double> mass(n, 1.0 / n), depth(n);
      for (int i = 0; i < n; ++i) depth[i] = (i + 1.0) / n;
      const double theta = 0.000157, hh = 5.0;
      std::vector<double> xs, ys, xr, yr;
      l.build_cumulative_vulnerability_integral(l.b, l.c,
                                                l.vulnerability_curve_ncontrol,
                                                xs, ys);
      l.build_cumulative_vulnerability_integral(l.root_b, l.root_c,
                                                l.vulnerability_curve_ncontrol,
                                                xr, yr);
      auto resolve = [&](int k, double v, bool hold_grid) {
        l.*slot[k] = v;
        if (k == 5 || k == 6) {
          if (hold_grid) l.set_transpiration_at(l.b, l.c, xs);
          else l.setup_transpiration(l.vulnerability_curve_ncontrol);
        }
        if (k == 12 || k == 13) {
          if (hold_grid) l.set_root_vulnerability_at(l.root_b, l.root_c, xr);
          else l.setup_root_vulnerability(l.vulnerability_curve_ncontrol);
        }
        l.transpiration_cached_ = false;
        // set_physiology keys the vcmax_/jmax_ block on (leaf_temp_,
        // atm_o2_kpa_) alone, so without this a perturbed vcmax_25 or jmax_25
        // never reaches the model and this reference reads exactly zero.
        l.photo_temp_cached_ = false;
        l.set_physiology(0.05, mass, 608, 0.0245, l.PPFD_, l.psi_soil_, depth,
                         1.0 * theta / hh, 2.0, 40.0, theta * hh, 25.0, 21.0,
                         101.3);
        l.find_root_collar_psi();
      };
      resolve(0, l.vcmax_25, false);
      std::vector<double> lamP(n, 0.0), adjP, lam0(n, 0.0), adj0;
      l.input_adjoints(1.0, lamP, adjP);
      lam0[0] = 1.0;
      l.input_adjoints(0.0, lam0, adj0);
      printf("    PARAMETER ROWS (adjoint | whole-solve difference)\n");
      for (int k = 0; k < 15; ++k) {
        const double keep = l.*slot[k];
        const double h = std::abs(keep) * 1e-6;
        if (!(h > 0.0)) { printf("      %-21s parameter is zero\n", pname[k]); continue; }
        double got_prof[2] = {0, 0}, got_up[2] = {0, 0};
        for (int grid = 0; grid < 2; ++grid) {
          double pr[2], up[2];
          for (int s = 0; s < 2; ++s) {
            resolve(k, keep + (s ? -h : h), grid == 1);
            pr[s] = l.profit_;
            up[s] = l.soil_consumption_[0];
          }
          got_prof[grid] = (pr[0] - pr[1]) / (2 * h);
          got_up[grid] = (up[0] - up[1]) / (2 * h);
          resolve(k, keep, grid == 1);
        }
        printf("      %-21s profit %13.6g | %13.6g moving %13.6g  rel %.2e\n",
               pname[k], adjP[i_par0 + k], got_prof[1], got_prof[0],
               rel(adjP[i_par0 + k], got_prof[1]));
        printf("      %-21s uptake %13.6g | %13.6g moving %13.6g  rel %.2e\n",
               "", adj0[i_par0 + k], got_up[1], got_up[0],
               rel(adj0[i_par0 + k], got_up[1]));
      }
      // vcmax_25 is the discriminator: uptake has no direct dependence on it, so
      // the whole row arrives through the operating point's movement and a
      // severed argmax channel returns exactly zero rather than a wrong number.
      printf("    DISCRIMINATOR d(uptake_0)/d(vcmax_25) = %.9g (nonzero=%d)\n",
             adj0[i_par0 + 0], (int)(adj0[i_par0 + 0] != 0.0));
    }

    // The bound rows, where the operating point is pinned. bound_a is where
    // uptake is zero; a tight bisection on E_up gives it to machine precision,
    // which prepare_collar_solve's own 1e-4 root-find does not, so differencing
    // that is a reference and differencing prepare_collar_solve is not.
    {
      l.find_root_collar_psi();
      if (l.collar_pinned_) {
        std::vector<double> lam(n, 0.0), adj4;
        lam[0] = 1.0;
        l.input_adjoints(0.0, lam, adj4);
        double ba, bb;
        l.prepare_collar_solve(ba, bb);
        l.find_root_collar_psi();
        const double p = -l.root_collar_psi_;
        printf("    BOUND p*=%.9g bound_a=%.9g bound_b=%.9g at_a=%d "
               "p*-bound=%.3e\n", p, ba, bb, (int)((p - ba) < (bb - p)),
               (p - ba) < (bb - p) ? p - ba : bb - p);
        auto zero_uptake_collar = [&]() -> double {
          double lo = 1e-6, hi = l.psi_crit;   // positive magnitudes
          for (int it = 0; it < 200; ++it) {
            const double m = 0.5 * (lo + hi);
            l.refresh_soil_potentials();
            l.E_from_Soil_to_Root_Collar(-m, l.psi_soil_inverted_);
            if (l.E_up_ > 0.0) hi = m; else lo = m;
          }
          return 0.5 * (lo + hi);
        };
        printf("      tight bound_a = %.12g vs prepare's %.12g\n",
               zero_uptake_collar(), ba);
        l.find_root_collar_psi();
        std::vector<double> dbound;
        l.bound_partials(dbound);
        l.find_root_collar_psi();
        for (int j = 0; j < n && j < 3; ++j) {
          const double keep = l.psi_soil_[j], h = 1e-6;
          double v[2];
          for (int s = 0; s < 2; ++s) {
            l.psi_soil_[j] = keep + (s ? -h : h);
            v[s] = zero_uptake_collar();
          }
          l.psi_soil_[j] = keep;
          l.refresh_soil_potentials();
          const double meas = (v[0] - v[1]) / (2 * h);
          printf("      d(bound_a)/d(psi_soil[%d]) rows %.9g tight difference "
                 "%.9g rel %.2e\n", j, dbound[1 + j], meas,
                 rel(dbound[1 + j], meas));
        }
        // Root mass and leaf area, against the same tight bisection.
        {
          std::vector<double> mass(n, 1.0 / n), depth(n);
          for (int i = 0; i < n; ++i) depth[i] = (i + 1.0) / n;
          const double theta = 0.000157, hh2 = 5.0;
          const double hm = (1.0 / n) * 1e-6;
          double v[2];
          for (int s = 0; s < 2; ++s) {
            std::vector<double> m2(n, 1.0 / n);
            m2[0] += (s ? -hm : hm);
            l.set_physiology(l.area_leaf_, m2, 608, 0.0245, l.PPFD_, l.psi_soil_,
                             depth, 1.0 * theta / hh2, 2.0, 40.0, theta * hh2,
                             25.0, 21.0, 101.3);
            v[s] = zero_uptake_collar();
          }
          l.set_physiology(l.area_leaf_, mass, 608, 0.0245, l.PPFD_, l.psi_soil_,
                           depth, 1.0 * theta / hh2, 2.0, 40.0, theta * hh2,
                           25.0, 21.0, 101.3);
          const double meas = (v[0] - v[1]) / (2 * hm);
          printf("      d(bound_a)/d(mass_root[0]) rows %.9g tight difference "
                 "%.9g rel %.2e\n", dbound[2 + n], meas,
                 rel(dbound[2 + n], meas));
          printf("      d(bound_a)/d(area_leaf)   rows %.9g (uptake is zero at "
                 "this endpoint, so leaf area cannot move it)\n", dbound[1 + n]);
        }
        // root_b and root_c reach the endpoint through the vulnerability curve.
        {
          const int idx[2] = {3 + 2 * n + 12, 3 + 2 * n + 13};
          double* mem[2] = {&l.root_b, &l.root_c};
          const char* nm3[2] = {"root_b", "root_c"};
          std::vector<double> xr, yr;
          l.build_cumulative_vulnerability_integral(
              l.root_b, l.root_c, l.vulnerability_curve_ncontrol, xr, yr);
          for (int t = 0; t < 2; ++t) {
            const double keep = *mem[t], hb = keep * 1e-6;
            double v[2][2];
            for (int grid = 0; grid < 2; ++grid) {
              for (int s = 0; s < 2; ++s) {
                *mem[t] = keep + (s ? -hb : hb);
                if (grid) l.set_root_vulnerability_at(l.root_b, l.root_c, xr);
                else l.setup_root_vulnerability(l.vulnerability_curve_ncontrol);
                v[grid][s] = zero_uptake_collar();
              }
              *mem[t] = keep;
              if (grid) l.set_root_vulnerability_at(l.root_b, l.root_c, xr);
              else l.setup_root_vulnerability(l.vulnerability_curve_ncontrol);
            }
            const double held = (v[1][0] - v[1][1]) / (2 * hb);
            printf("      d(bound_a)/d(%s) rows %.9g held-grid %.9g rel %.2e | "
                   "moving-grid %.9g\n", nm3[t], dbound[idx[t]], held,
                   rel(dbound[idx[t]], held), (v[0][0] - v[0][1]) / (2 * hb));
          }
        }
        l.find_root_collar_psi();
      } else {
        printf("    BOUND not pinned at this state\n");
      }
    }
  }
  return 0;
}
