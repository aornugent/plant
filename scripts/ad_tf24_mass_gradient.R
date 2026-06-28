# TF24 NET-PRODUCTION trait gradients for the MASS-CASCADE traits (#472 scope B,
# Phase F1-full). Companion to ad_tf24_net_gradient.R (vcmax_25),
# ad_tf24_hydraulic_gradient.R and ad_tf24_photo_gradient.R (the leaf traits).
#
# TF24 net production is
#   net = a_bio*a_y*(profit*area_leaf*conv - respiration) - turnover,
# and the respiration / turnover / mass-cascade algebra is IDENTICAL to FF16:
#   mass_leaf    = area_leaf * lma
#   area_sapwood = area_leaf * theta;  mass_sapwood = area_sapwood*height*eta_c*rho
#   area_bark    = a_b1*area_leaf*theta;  mass_bark = area_bark*height*eta_c*rho
#   mass_root    = a_r1 * area_leaf
#   respiration  = r_l*mass_leaf + r_b*mass_bark + r_s*mass_sapwood + r_r*mass_root
#   turnover     = k_l*mass_leaf + k_b*mass_bark + k_s*mass_sapwood + k_r*mass_root
# so the mass-cascade trait gradients come straight from a scalar-templated kernel
# (net_kernel<T> below), forward-AD per trait. Three pathways:
#
#  (1) PURE (no leaf coupling): lma, rho, a_b1, r_l, r_b, r_s, r_r, k_l, k_b, k_s,
#      k_r, a_bio, a_y -- profit and area_leaf are frozen; the trait moves only
#      respiration / turnover (exactly FF16). Exact to ~1e-10.
#  (2) area_leaf-active: a_l1, a_l2 set area_leaf = (height/a_l1)^(1/a_l2). The
#      profit PER LEAF AREA is area_leaf-independent (the root resistances scale
#      as 1/area_leaf, cancelling the 1/area_leaf in the soil->collar uptake), so
#      profit stays frozen and area_leaf is the only active input.
#  (3) leaf-coupled (inject the leaf-profit sensitivity, the Phase-D pattern):
#      - theta also sets k_max = K_s*theta/(h*eta_c): d profit/d theta =
#        dprofit_dkmax * (k_max/theta);
#      - a_r1 also scales every root hydraulic resistance by 1/a_r1, hence the
#        uptake E_up_ linearly: d profit/d a_r1 = dprofit_dEup * (E_up_/a_r1).
#
# Validated end-to-end vs a finite difference of TF24_Strategy::net_mass_production_dt.
# Run from the package root after `R CMD INSTALL .`:
#   Rscript scripts/ad_tf24_mass_gradient.R

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
#include <XAD/XAD.hpp>
// [[Rcpp::plugins(cpp20)]]
using AD = xad::fwd<double>::active_type;

// Scalar-templated TF24 net production. The mass cascade is FF16-identical; only
// assim = profit*area_leaf*conv is TF24-specific (profit = optimised leaf profit).
template <typename T>
T net_kernel(T profit, T area_leaf, double height, double eta_c,
             T lma, T rho, T theta, T a_b1, T a_r1,
             T r_l, T r_b, T r_s, T r_r, T k_l, T k_b, T k_s, T k_r,
             T a_bio, T a_y) {
  const double conv = 60.0*60.0*12.0*365.0/1e6;
  T mass_leaf    = area_leaf * lma;
  T area_sapwood = area_leaf * theta;
  T mass_sapwood = area_sapwood * height * eta_c * rho;
  T area_bark    = a_b1 * area_leaf * theta;
  T mass_bark    = area_bark * height * eta_c * rho;
  T mass_root    = a_r1 * area_leaf;
  T resp = r_l*mass_leaf + r_b*mass_bark + r_s*mass_sapwood + r_r*mass_root;
  T turn = k_l*mass_leaf + k_b*mass_bark + k_s*mass_sapwood + k_r*mass_root;
  return a_bio*a_y*(profit*area_leaf*conv - resp) - turn;
}

static double net_live(const std::string& t, double v, double light, double h) {
  plant::TF24_Strategy s; s.control.shading_model = "crown-centre";
  s.control.GSS_tol_abs = 1e-9;
  auto& p = s.pars;
  if (t=="lma")p.lma=v; else if(t=="rho")p.rho=v; else if(t=="theta")p.theta=v;
  else if(t=="a_b1")p.a_b1=v; else if(t=="a_r1")p.a_r1=v;
  else if(t=="a_l1")p.a_l1=v; else if(t=="a_l2")p.a_l2=v;
  else if(t=="r_l")p.r_l=v; else if(t=="r_b")p.r_b=v; else if(t=="r_s")p.r_s=v;
  else if(t=="r_r")p.r_r=v; else if(t=="k_l")p.k_l=v; else if(t=="k_b")p.k_b=v;
  else if(t=="k_s")p.k_s=v; else if(t=="k_r")p.k_r=v;
  else if(t=="a_bio")p.a_bio=v; else if(t=="a_y")p.a_y=v; else Rcpp::stop("?");
  s.prepare_strategy(); plant::TF24_Environment e; e.set_fixed_environment(light, 1e4);
  return s.net_mass_production_dt(e, h, s.area_leaf(h), 1.0/h);
}

// [[Rcpp::export]]
Rcpp::NumericVector tf24_dnet_dmass(std::string t, double light, double h) {
  plant::TF24_Strategy s; s.control.shading_model = "crown-centre";
  s.control.GSS_tol_abs = 1e-9; s.prepare_strategy();
  plant::TF24_Environment e; e.set_fixed_environment(light, 1e4);
  const double al = s.area_leaf(h);
  s.net_mass_production_dt(e, h, al, 1.0/h);
  const double prof = s.leaf.profit_;
  auto& p = s.pars;
  const double opt = -s.leaf.root_collar_psi_;

  auto val = [&](const std::string& n)->double{
    if(n=="lma")return p.lma; if(n=="rho")return p.rho; if(n=="theta")return p.theta;
    if(n=="a_b1")return p.a_b1; if(n=="a_r1")return p.a_r1; if(n=="a_l1")return p.a_l1;
    if(n=="a_l2")return p.a_l2; if(n=="r_l")return p.r_l; if(n=="r_b")return p.r_b;
    if(n=="r_s")return p.r_s; if(n=="r_r")return p.r_r; if(n=="k_l")return p.k_l;
    if(n=="k_b")return p.k_b; if(n=="k_s")return p.k_s; if(n=="k_r")return p.k_r;
    if(n=="a_bio")return p.a_bio; if(n=="a_y")return p.a_y; return 0; };

  // AD copies of every mass-cascade par; seed exactly one below.
  AD lma=p.lma,rho=p.rho,theta=p.theta,a_b1=p.a_b1,a_r1=p.a_r1,
     r_l=p.r_l,r_b=p.r_b,r_s=p.r_s,r_r=p.r_r,k_l=p.k_l,k_b=p.k_b,k_s=p.k_s,k_r=p.k_r,
     a_bio=p.a_bio,a_y=p.a_y;
  AD area_leaf = AD(al);   // frozen unless a_l1/a_l2
  AD profit    = AD(prof); // frozen unless theta/a_r1 (leaf-coupled)

  if (t=="a_l1" || t=="a_l2") {
    AD a_l1=p.a_l1, a_l2=p.a_l2;
    if (t=="a_l1") xad::derivative(a_l1)=1.0; else xad::derivative(a_l2)=1.0;
    area_leaf = pow(AD(h)/a_l1, 1.0/a_l2);          // profit per area is area_leaf-invariant
  } else if (t=="theta") {
    const double kmax = s.leaf.leaf_specific_conductance_max_;
    const double dprof = s.leaf.dprofit_dkmax(opt) * (kmax / p.theta);
    xad::derivative(theta)=1.0;
    profit = AD(prof) + AD(dprof)*(theta - AD(p.theta));   // inject leaf sensitivity
  } else if (t=="a_r1") {
    const double dprof = s.leaf.dprofit_dEup(opt) * (s.leaf.E_up_ / p.a_r1);
    xad::derivative(a_r1)=1.0;
    profit = AD(prof) + AD(dprof)*(a_r1 - AD(p.a_r1));      // inject leaf sensitivity
  } else {  // pure mass-cascade trait: seed it, profit + area_leaf frozen
    AD* tgt=nullptr;
    if(t=="lma")tgt=&lma; else if(t=="rho")tgt=&rho; else if(t=="a_b1")tgt=&a_b1;
    else if(t=="r_l")tgt=&r_l; else if(t=="r_b")tgt=&r_b; else if(t=="r_s")tgt=&r_s;
    else if(t=="r_r")tgt=&r_r; else if(t=="k_l")tgt=&k_l; else if(t=="k_b")tgt=&k_b;
    else if(t=="k_s")tgt=&k_s; else if(t=="k_r")tgt=&k_r; else if(t=="a_bio")tgt=&a_bio;
    else if(t=="a_y")tgt=&a_y; else Rcpp::stop("unknown trait");
    xad::derivative(*tgt)=1.0;
  }

  AD net = net_kernel<AD>(profit, area_leaf, h, s.eta_c, lma, rho, theta, a_b1,
    a_r1, r_l, r_b, r_s, r_r, k_l, k_b, k_s, k_r, a_bio, a_y);
  const double ad = xad::derivative(net);
  const double v0 = val(t), hh = 1e-6*std::abs(v0);
  const double fd = (net_live(t,v0+hh,light,h) - net_live(t,v0-hh,light,h))/(2*hh);
  return Rcpp::NumericVector::create(Rcpp::_["AD"]=ad, Rcpp::_["FD"]=fd);
}')

traits <- c("lma","rho","a_b1","r_l","r_b","r_s","r_r","k_l","k_b","k_s","k_r",
            "a_bio","a_y",          # pure (FF16-identical)
            "a_l1","a_l2",          # area_leaf-active
            "theta","a_r1")         # leaf-coupled (injection)
ok <- TRUE
cat("TF24 d(net_mass_production_dt)/d(mass-cascade trait): AD vs strategy FD (light=0.7, h=5)\n")
for (t in traits) {
  r <- tf24_dnet_dmass(t, 0.7, 5.0)
  ae <- abs(r[["AD"]] - r[["FD"]]); re <- ae / max(abs(r[["FD"]]), 1e-30)
  pass <- re < 1e-5 || ae < 1e-9
  ok <- ok && pass
  cat(sprintf("  %-6s AD=%- 14.7g FD=%- 14.7g rel=%.1e %s\n",
              t, r[["AD"]], r[["FD"]], re, if (pass) "OK" else "** MISMATCH **"))
}
stopifnot(ok)
cat("\nAll 17 TF24 mass-cascade trait gradients validated vs net FD.\n")
