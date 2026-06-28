# First TF24 NET-PRODUCTION trait gradient (#472 scope B, Phase F1): exact
# d(net_mass_production_dt)/d(vcmax_25) for the TF24 strategy, via the leaf-level
# d(profit*)/d(vcmax_25) (envelope theorem + IFT, Leaf::dprofit_dvcmax25) carried
# through the net-production assembly.
#
# TF24 net production is
#   net = a_bio * a_y * (assimilation - respiration) - turnover,
#   assimilation = leaf.profit_ * area_leaf * (60*60*12*365/1e6),
# where profit_ is the OPTIMISED leaf profit (a max over collar potential nesting a
# psi_stem->ci root-find). vcmax_25 enters net ONLY through the leaf profit (the mass
# cascade and area_leaf are vcmax-independent), so
#   d(net)/d(vcmax_25) = a_bio * a_y * area_leaf * conv * d(profit*)/d(vcmax_25).
#
# Two things this confirms:
#  (a) the real TF24 strategy operates at an INTERIOR leaf optimum (positive profit,
#      dprofit/dcollar ~ 0) where the envelope theorem applies -- so freezing the
#      optimal collar and differentiating the partial is valid (it is NOT a stressed
#      boundary/shut-down optimum);
#  (b) the leaf-trait gradient plugs into the strategy net production exactly: AD
#      matches a finite difference of the live TF24_Strategy::net_mass_production_dt.
#
# This is the de-risked foundation for the full TF24 net-production kernel (all leaf
# + mass-cascade traits) and the TF24 emergent gradient via the two-pass replay.
#
# Run from the package root after `R CMD INSTALL .`:
#   Rscript scripts/ad_tf24_net_gradient.R

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

static double tf24_net(double vcmax_25, double light, double height) {
  plant::TF24_Strategy s; s.control.shading_model = "crown-centre";
  s.pars.vcmax_25 = vcmax_25; s.prepare_strategy();
  plant::TF24_Environment env; env.set_fixed_environment(light, 1e4);
  return s.net_mass_production_dt(env, height, s.area_leaf(height), 1.0 / height);
}

// [[Rcpp::export]]
Rcpp::NumericVector tf24_dnet_dvcmax(double light, double height) {
  plant::TF24_Strategy s; s.control.shading_model = "crown-centre"; s.prepare_strategy();
  plant::TF24_Environment env; env.set_fixed_environment(light, 1e4);
  const double al = s.area_leaf(height);
  const double net = s.net_mass_production_dt(env, height, al, 1.0 / height);
  const double opt = -s.leaf.root_collar_psi_;       // optimised collar potential
  const double dpdcollar = s.leaf.dprofit_droot_collar_psi(opt);  // ~0 => interior
  const double conv = 60.0 * 60.0 * 12.0 * 365.0 / 1e6;
  const double ad = s.pars.a_bio * s.pars.a_y * al * conv * s.leaf.dprofit_dvcmax25(opt);
  const double vc0 = s.pars.vcmax_25, h = 1e-5 * vc0;
  const double fd = (tf24_net(vc0 + h, light, height) - tf24_net(vc0 - h, light, height)) / (2 * h);
  return Rcpp::NumericVector::create(Rcpp::_["net"] = net, Rcpp::_["profit"] = s.leaf.profit_,
    Rcpp::_["dprofit_dcollar"] = dpdcollar, Rcpp::_["AD"] = ad, Rcpp::_["FD"] = fd);
}')

cat("TF24 d(net_mass_production_dt)/d(vcmax_25): AD (leaf envelope+IFT) vs strategy FD\n")
ok <- TRUE
for (light in c(0.4, 0.7, 1.0)) {
  r <- tf24_dnet_dvcmax(light, 5.0)
  re <- abs(r[["AD"]] - r[["FD"]]) / max(abs(r[["FD"]]), 1e-30)
  ok <- ok && re < 1e-5 && r[["profit"]] > 0 && abs(r[["dprofit_dcollar"]]) < 1e-2
  cat(sprintf("  light=%.1f  net=%8.5f profit=%7.3f (interior: dp/dcollar=%.1e)  AD=%.7g FD=%.7g rel=%.1e %s\n",
              light, r[["net"]], r[["profit"]], r[["dprofit_dcollar"]], r[["AD"]], r[["FD"]], re,
              if (re < 1e-5) "OK" else "** MISMATCH **"))
}
stopifnot(ok)
cat("\nTF24 net-production gradient validated (interior optimum; leaf gradient -> net).\n")
