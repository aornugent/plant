# TF24 NET-PRODUCTION trait gradients for the PHOTOSYNTHESIS leaf traits (#472
# scope B, Phase F1-full). Companion to ad_tf24_net_gradient.R (vcmax_25) and
# ad_tf24_hydraulic_gradient.R (the hydraulic traits).
#
# jmax_25, a (quantum yield) and the two curvature factors (curv_fact_elec_trans,
# curv_fact_colim) are vcmax-like: they affect ONLY assimilation, not the
# transport (psi_stem) or the hydraulic cost. So by the envelope theorem (optimal
# collar frozen, dprofit/dcollar ~ 0) the leaf gradient follows the
# dprofit_dvcmax25 pattern:
#   dprofit/dt = A_t + A'(ci) * dci/dt,   dci/dt = -(A_t * umol_to_mol) / g_ci,
# where A_t = d(assim_colimited)/dt holding ci. jmax_25, a and curv_elec enter via
# the electron-transport rate et (A_t = A_et * det/dt; jmax_25 chains the linear
# jmax_/jmax_25); curv_colim enters the colimitation min directly
# (Leaf::dprofit_djmax25 / dprofit_da / dprofit_dcurv_elec / dprofit_dcurv_colim).
#
# Each enters TF24 net production ONLY through the optimised leaf profit, so
#   d(net)/d(trait) = a_bio * a_y * area_leaf * conv * d(profit*)/d(trait),
# validated end-to-end vs a finite difference of TF24_Strategy::net_mass_production_dt.
#
# Run from the package root after `R CMD INSTALL .`:
#   Rscript scripts/ad_tf24_photo_gradient.R

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

static double tf24_net(const std::string& trait, double val, double light, double height) {
  plant::TF24_Strategy s; s.control.shading_model = "crown-centre";
  s.control.GSS_tol_abs = 1e-9;
  if (trait == "jmax_25")        s.pars.jmax_25 = val;
  else if (trait == "a")         s.pars.a = val;
  else if (trait == "curv_elec") s.pars.curv_fact_elec_trans = val;
  else if (trait == "curv_colim") s.pars.curv_fact_colim = val;
  else Rcpp::stop("unknown trait");
  s.prepare_strategy();
  plant::TF24_Environment env; env.set_fixed_environment(light, 1e4);
  return s.net_mass_production_dt(env, height, s.area_leaf(height), 1.0 / height);
}

// [[Rcpp::export]]
Rcpp::NumericVector tf24_dnet_dphoto(std::string trait, double light, double height) {
  plant::TF24_Strategy s; s.control.shading_model = "crown-centre";
  s.control.GSS_tol_abs = 1e-9; s.prepare_strategy();
  plant::TF24_Environment env; env.set_fixed_environment(light, 1e4);
  const double al  = s.area_leaf(height);
  const double net = s.net_mass_production_dt(env, height, al, 1.0 / height);
  const double opt = -s.leaf.root_collar_psi_;
  const double dpdcollar = s.leaf.dprofit_droot_collar_psi(opt);
  const double conv = 60.0 * 60.0 * 12.0 * 365.0 / 1e6;
  const double scale = s.pars.a_bio * s.pars.a_y * al * conv;

  double dprofit, v0;
  if (trait == "jmax_25")    { dprofit = s.leaf.dprofit_djmax25(opt);    v0 = s.pars.jmax_25; }
  else if (trait == "a")     { dprofit = s.leaf.dprofit_da(opt);         v0 = s.pars.a; }
  else if (trait == "curv_elec")  { dprofit = s.leaf.dprofit_dcurv_elec(opt);  v0 = s.pars.curv_fact_elec_trans; }
  else if (trait == "curv_colim") { dprofit = s.leaf.dprofit_dcurv_colim(opt); v0 = s.pars.curv_fact_colim; }
  else Rcpp::stop("unknown trait");

  const double ad = scale * dprofit;
  const double h  = 1e-6 * std::abs(v0);
  const double fd = (tf24_net(trait, v0 + h, light, height) -
                     tf24_net(trait, v0 - h, light, height)) / (2 * h);
  return Rcpp::NumericVector::create(Rcpp::_["net"] = net, Rcpp::_["profit"] = s.leaf.profit_,
    Rcpp::_["dprofit_dcollar"] = dpdcollar, Rcpp::_["AD"] = ad, Rcpp::_["FD"] = fd);
}')

ok <- TRUE
for (trait in c("jmax_25", "a", "curv_elec", "curv_colim")) {
  cat(sprintf("\nTF24 d(net_mass_production_dt)/d(%s): AD (envelope) vs strategy FD\n", trait))
  for (light in c(0.4, 0.7, 1.0)) {
    r <- tf24_dnet_dphoto(trait, light, 5.0)
    ae <- abs(r[["AD"]] - r[["FD"]])
    re <- ae / max(abs(r[["FD"]]), 1e-30)
    pass <- re < 1e-5 || ae < 1e-7
    ok <- ok && pass && r[["profit"]] > 0 && abs(r[["dprofit_dcollar"]]) < 1e-2
    cat(sprintf("  light=%.1f  profit=%7.3f (interior: dp/dcollar=%.1e)  AD=%.7g FD=%.7g rel=%.1e %s\n",
                light, r[["profit"]], r[["dprofit_dcollar"]], r[["AD"]], r[["FD"]], re,
                if (pass) "OK" else "** MISMATCH **"))
  }
}
stopifnot(ok)
cat("\nTF24 photosynthesis gradients (jmax_25, a, curv_elec, curv_colim) validated vs net FD.\n")
