# TF24 NET-PRODUCTION trait gradients for the HYDRAULIC leaf traits (#472 scope B,
# Phase F1-full). Extends scripts/ad_tf24_net_gradient.R (vcmax_25) to the traits
# that move the hydraulic COST and/or the transport (psi_stem), not just assim.
#
# This file covers:
#  - the COST-ONLY hydraulic traits g1_TF24 and beta2, which enter
#      C = g1_TF24 * (1 - exp(-(psi_stem/b)^c))^beta2
#    but NOT the transport (psi_stem) nor assimilation (ci). By the envelope
#    theorem (optimal collar frozen) their leaf gradient is just minus the
#    explicit cost derivative (Leaf::dprofit_dg1_TF24 / dprofit_dbeta2);
#  - the TRANSPORT trait K_s, which scales the supply-side conductance
#    k_max = K_s*theta/(h*eta_c) linearly; it moves psi_stem and ci (not the cost
#    explicitly), handled by the transport+IFT pattern (Leaf::dprofit_dkmax,
#    chained by k_max/K_s);
#  - the VULNERABILITY-SHAPE traits b and c (prop_cond = exp(-(psi/b)^c)), which
#    reshape the transpiration spline AND enter the cost explicitly. ci/benefit
#    are frozen (operating-point transpiration = the root-vulnerability uptake
#    E_up_), so dprofit/dt = -(C'(psi_stem)*dpsi_stem/dt + dC/dt|explicit), with
#    dpsi_stem/dt from the exact dS/dt of the cumulative curve (Leaf::dprofit_db /
#    dprofit_dc).
#
# Each hydraulic trait enters TF24 net production ONLY through the optimised leaf
# profit (the mass cascade, respiration, turnover and area_leaf are all
# hydraulic-trait-independent), so
#   d(net)/d(trait) = a_bio * a_y * area_leaf * conv * d(profit*)/d(trait).
# Validated end-to-end vs a finite difference of the live
# TF24_Strategy::net_mass_production_dt.
#
# Run from the package root after `R CMD INSTALL .`:
#   Rscript scripts/ad_tf24_hydraulic_gradient.R

suppressMessages({library(Rcpp); library(plant)})

plant_inc  <- system.file("include", package = "plant")
odelia_inc <- system.file("include", package = "odelia")
bh_inc     <- system.file("include", package = "BH")
plant_so   <- system.file("libs", "plant.so",  package = "plant")
odelia_so  <- system.file("libs", "odelia.so", package = "odelia")
if (!all(nzchar(c(plant_inc, odelia_inc, bh_inc))) ||
    !all(file.exists(c(plant_so, odelia_so))))
  stop("Need plant (installed from this branch), odelia and BH; run R CMD INSTALL .")
Sys.setenv(PKG_CPPFLAGS = paste(paste0("-I", shQuote(plant_inc)),
                                paste0("-I", shQuote(odelia_inc)),
                                paste0("-I", shQuote(bh_inc))))
Sys.setenv(PKG_LIBS = paste(shQuote(normalizePath(plant_so)),
                            shQuote(normalizePath(odelia_so))))

Rcpp::sourceCpp(code = '
#include <Rcpp.h>
#include <plant/models/tf24_strategy.h>
#include <plant/models/tf24_environment.h>
// [[Rcpp::plugins(cpp20)]]

// Live TF24 net production with a single hydraulic trait perturbed (rebuilds the
// strategy so the leaf is reconfigured exactly as prepare_strategy would).
static double tf24_net(const std::string& trait, double val, double light, double height) {
  plant::TF24_Strategy s; s.control.shading_model = "crown-centre";
  // Tighten the collar golden-section so the envelope theorem (collar frozen at
  // the optimum) is exact to high precision: the default GSS_tol_abs=1e-3 leaves
  // dprofit/dcollar ~ -4e-4, and for a transport trait (K_s) the collar moves
  // enough with the trait that this residual shows up as a ~1e-4 AD-vs-FD gap.
  s.control.GSS_tol_abs = 1e-9;
  // b/c reshape the transpiration spline; the AD uses the EXACT continuous dS/dt
  // (closed form for b, high-accuracy quadrature for c) while the FD rebuilds the
  // spline, so the AD-vs-FD gap is the spline interpolation error. At the default
  // ncontrol=100 that is ~5e-6; at 2000 it falls to ~1e-8 (confirming the AD is
  // the exact derivative). Use the dense spline for an unambiguous check.
  s.control.vulnerability_curve_ncontrol = 2000;
  if (trait == "g1_TF24")  s.g1_TF24  = val;
  else if (trait == "beta2") s.pars.beta2 = val;
  else if (trait == "K_s")   s.pars.K_s   = val;
  // b/c perturbed alone (psi_crit held at its stale default), matching the AD,
  // which differentiates the cost+transport at fixed psi_crit and c (resp. b).
  else if (trait == "b")     s.pars.b     = val;
  else if (trait == "c")     s.pars.c     = val;
  else Rcpp::stop("unknown trait");
  s.prepare_strategy();
  plant::TF24_Environment env; env.set_fixed_environment(light, 1e4);
  return s.net_mass_production_dt(env, height, s.area_leaf(height), 1.0 / height);
}

// [[Rcpp::export]]
Rcpp::NumericVector tf24_dnet_dhydraulic(std::string trait, double light, double height) {
  plant::TF24_Strategy s; s.control.shading_model = "crown-centre";
  s.control.GSS_tol_abs = 1e-9;
  // b/c reshape the transpiration spline; the AD uses the EXACT continuous dS/dt
  // (closed form for b, high-accuracy quadrature for c) while the FD rebuilds the
  // spline, so the AD-vs-FD gap is the spline interpolation error. At the default
  // ncontrol=100 that is ~5e-6; at 2000 it falls to ~1e-8 (confirming the AD is
  // the exact derivative). Use the dense spline for an unambiguous check.
  s.control.vulnerability_curve_ncontrol = 2000; s.prepare_strategy();
  plant::TF24_Environment env; env.set_fixed_environment(light, 1e4);
  const double al  = s.area_leaf(height);
  const double net = s.net_mass_production_dt(env, height, al, 1.0 / height);
  const double opt = -s.leaf.root_collar_psi_;                    // optimal collar
  const double dpdcollar = s.leaf.dprofit_droot_collar_psi(opt);  // ~0 => interior
  const double conv = 60.0 * 60.0 * 12.0 * 365.0 / 1e6;
  const double scale = s.pars.a_bio * s.pars.a_y * al * conv;

  double dprofit, v0;
  if (trait == "g1_TF24")    { dprofit = s.leaf.dprofit_dg1_TF24(opt); v0 = s.g1_TF24; }
  else if (trait == "beta2") { dprofit = s.leaf.dprofit_dbeta2(opt);   v0 = s.pars.beta2; }
  else if (trait == "K_s") {
    // k_max = K_s * theta / (h*eta_c): chain leaf dprofit/dkmax by k_max/K_s.
    const double kmax = s.leaf.leaf_specific_conductance_max_;
    dprofit = s.leaf.dprofit_dkmax(opt) * (kmax / s.pars.K_s);
    v0 = s.pars.K_s;
  }
  else if (trait == "b") { dprofit = s.leaf.dprofit_db(opt); v0 = s.pars.b; }
  else if (trait == "c") { dprofit = s.leaf.dprofit_dc(opt); v0 = s.pars.c; }
  else Rcpp::stop("unknown trait");

  const double ad = scale * dprofit;
  const double h  = 1e-6 * std::abs(v0);
  const double fd = (tf24_net(trait, v0 + h, light, height) -
                     tf24_net(trait, v0 - h, light, height)) / (2 * h);
  return Rcpp::NumericVector::create(Rcpp::_["net"] = net, Rcpp::_["profit"] = s.leaf.profit_,
    Rcpp::_["dprofit_dcollar"] = dpdcollar, Rcpp::_["AD"] = ad, Rcpp::_["FD"] = fd);
}')

ok <- TRUE
for (trait in c("g1_TF24", "beta2", "K_s", "b", "c")) {
  cat(sprintf("\nTF24 d(net_mass_production_dt)/d(%s): AD (envelope) vs strategy FD\n", trait))
  for (light in c(0.4, 0.7, 1.0)) {
    r <- tf24_dnet_dhydraulic(trait, light, 5.0)
    ae <- abs(r[["AD"]] - r[["FD"]])
    re <- ae / max(abs(r[["FD"]]), 1e-30)
    pass <- re < 1e-5 || ae < 1e-6   # small-value cases: judge on absolute error
    ok <- ok && pass && r[["profit"]] > 0 && abs(r[["dprofit_dcollar"]]) < 1e-2
    cat(sprintf("  light=%.1f  profit=%7.3f (interior: dp/dcollar=%.1e)  AD=%.7g FD=%.7g rel=%.1e %s\n",
                light, r[["profit"]], r[["dprofit_dcollar"]], r[["AD"]], r[["FD"]], re,
                if (pass) "OK" else "** MISMATCH **"))
  }
}
stopifnot(ok)
cat("\nTF24 hydraulic gradients (g1_TF24, beta2, K_s, b, c) validated vs net FD.\n")
