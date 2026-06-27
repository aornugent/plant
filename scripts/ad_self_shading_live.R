# Self-shading gradient on the LIVE resident stand: the resident light RESPONDS to
# the trait (#472 scope B, Milestone C -- the active-knot / resident-reshaping path).
#
# All the earlier emergent-gradient scripts hold the resident light FROZEN (the
# mutant-through-frozen-canopy / invasion-fitness gradient: correct for a rare mutant
# that does not perturb the canopy). When the focal trait IS the resident's, an
# allometric trait (a_l1, a_l2) reshapes EVERY cohort's leaf area, hence the whole
# Beer's-law canopy and the light every plant reads. This script differentiates a
# focal output THROUGH that self-shaded light.
#
# The live SCM exposes, per node, the competition_effect ce_i and height h_i; the
# resident competition is competition(z) = trapezium_i( ce_i * Q(z/h_i) ), Q the
# Yokozawa leaf-area-above (1-u^eta)^2, and light(z) = exp(-competition(z)) (matching
# Patch::compute_competition + FF16_Environment Beer's law). Factor ce_i = C_i *
# area_leaf_i with C_i = ce_i / area_leaf_i frozen (density * survival weighting): then
#   competition(z; theta) = trapezium_i( C_i * area_leaf_i(theta) * Q(z/h_i) )
# is ACTIVE in the allometric trait through area_leaf_i, reconstructing the live light
# at the base trait and responding to perturbations -- the active-knot light.
#
# Validation: (a) the reconstructed light matches the live FF16_Environment (to the
# env spline's own interpolation tolerance); (b) d(focal net production)/d(a_l1) with
# the self-shaded light ACTIVE matches a two-pass FD over the same reconstruction; and
# the frozen-light value is reported to isolate the self-shading contribution.
#
# This is the static-census demonstration (final stand). The fully time-integrated
# version -- the active light at every replay step -- needs the per-step stand state
# (heights + ce per ODE step), a C++ harvest beyond environment_history; noted as the
# follow-up.
#
# Run from the package root after `R CMD INSTALL .`:
#   Rscript scripts/ad_self_shading_live.R

suppressMessages({library(Rcpp); library(plant)})

p <- scm_base_parameters("FF16")
p <- add_strategies(p, trait_matrix(0.0825, "lma"), hyperpar = FF16_hyperpar,
                    birth_rate = list(20))
p <- run_scm(p, Environment("FF16"), control(), refine_schedule = TRUE)$parameters
scm <- run_scm(p, Environment("FF16"), control(save_RK45_cache = TRUE),
               refine_schedule = FALSE)

sp <- scm$patch$species[[1]]
s  <- p$strategies[[1]]
h  <- sp$heights
ce <- sp$compute_competition_effect_by_nodes      # per-node competition_effect
a_l1 <- s$pars$a_l1; a_l2 <- s$pars$a_l2
# area_leaf is the allometry inverse, ff16_area_leaf = (height/a_l1)^(1/a_l2);
# C_i = ce_i / area_leaf_i is the frozen per-node weight (density * survival * k_I)
# so that C_i * ff16_area_leaf(a_l1,a_l2,h_i) reconstructs ce_i and responds to the trait.
C  <- ce / (h / a_l1)^(1 / a_l2)
# descending height order for the trapezium (matches Species::compute_competition)
o <- order(h, decreasing = TRUE)
h_desc <- h[o]; C_desc <- C[o]
pp <- unlist(s$pars)
cat(sprintf("Live stand: %d cohorts, heights %.2f..%.2f m\n", length(h), min(h), max(h)))

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
#include <vector>
#include <cmath>
#include <XAD/XAD.hpp>
#include <plant.h>
#include <plant/models/ff16_production_kernel.h>
using ad   = xad::adj<double>;
using ad_t = ad::active_type;
// [[Rcpp::plugins(cpp20)]]

static double as_double(double v)      { return v; }
static double as_double(const ad_t& v) { return xad::value(v); }

static plant::FF16_Strategy make_strategy(const Rcpp::NumericVector& pp) {
  plant::FF16_Strategy s; auto& q = s.pars;
  q.lma=pp["lma"];q.rho=pp["rho"];q.hmat=pp["hmat"];q.omega=pp["omega"];q.eta=pp["eta"];
  q.theta=pp["theta"];q.a_l1=pp["a_l1"];q.a_l2=pp["a_l2"];q.a_r1=pp["a_r1"];q.a_b1=pp["a_b1"];
  q.r_s=pp["r_s"];q.r_b=pp["r_b"];q.r_r=pp["r_r"];q.r_l=pp["r_l"];q.a_y=pp["a_y"];q.a_bio=pp["a_bio"];
  q.k_l=pp["k_l"];q.k_b=pp["k_b"];q.k_s=pp["k_s"];q.k_r=pp["k_r"];q.a_p1=pp["a_p1"];q.a_p2=pp["a_p2"];
  q.a_f3=pp["a_f3"];q.a_f1=pp["a_f1"];q.a_f2=pp["a_f2"];q.S_D=pp["S_D"];q.a_d0=pp["a_d0"];q.d_I=pp["d_I"];
  q.a_dG1=pp["a_dG1"];q.a_dG2=pp["a_dG2"];q.k_I=pp["k_I"];q.recruitment_decay=pp["recruitment_decay"];
  s.prepare_strategy(); return s;
}
template <typename S> static plant::FF16ProdPars<S> lift(const plant::FF16ProdPars<double>& d) {
  plant::FF16ProdPars<S> p;
  p.lma=d.lma;p.rho=d.rho;p.theta=d.theta;p.a_b1=d.a_b1;p.a_r1=d.a_r1;p.eta_c=d.eta_c;
  p.a_p1=d.a_p1;p.a_p2=d.a_p2;p.r_l=d.r_l;p.r_s=d.r_s;p.r_b=d.r_b;p.r_r=d.r_r;
  p.k_l=d.k_l;p.k_b=d.k_b;p.k_s=d.k_s;p.k_r=d.k_r;p.a_bio=d.a_bio;p.a_y=d.a_y;
  p.a_l1=d.a_l1;p.a_l2=d.a_l2;p.a_f1=d.a_f1;p.a_f2=d.a_f2;p.hmat=d.hmat;
  p.omega=d.omega;p.a_f3=d.a_f3;p.d_I=d.d_I;p.a_dG1=d.a_dG1;p.a_dG2=d.a_dG2; return p;
}

// Resident self-shaded light at z, ACTIVE in the allometric trait via area_leaf_i.
// competition(z) = (1/2) sum_adjacent (h_i - h_{i+1}) (g_i + g_{i+1}),
//   g_i = C_i * area_leaf(a_l1,a_l2,h_i) * Q(z/h_i),  Q = (1-u^eta)^2 (u=z/h_i<1).
// Heights, the trapezium spacing and C_i are frozen pass-1 doubles; light = exp(-comp).
template <typename S>
static S recon_light(double z, const plant::FF16ProdPars<S>& p, double eta,
                     const std::vector<double>& h, const std::vector<double>& C) {
  using std::pow; using std::exp;
  auto g = [&](std::size_t i) -> S {
    if (z >= h[i]) return S(0.0);
    double u = z / h[i]; double om = 1.0 - pow(u, eta);
    return S(C[i]) * plant::ff16_area_leaf(p.a_l1, p.a_l2, S(h[i])) * S(om * om);
  };
  S comp = S(0.0);
  S g_prev = g(0); double h_prev = h[0];
  for (std::size_t i = 1; i < h.size(); ++i) {
    S gi = g(i);
    comp = comp + S(h_prev - h[i]) * (g_prev + gi);
    h_prev = h[i]; g_prev = gi;
  }
  return exp(-S(0.5) * comp);
}

// Focal-cohort net production (crown-top) reading the self-shaded light at its crown.
template <typename S>
static S focal_net(const plant::FF16ProdPars<S>& p, double focal_h, double eta,
                   const std::vector<double>& hv, const std::vector<double>& C,
                   bool active_light) {
  S Ef = recon_light<S>(focal_h * as_double(p.eta_c), p, eta, hv, C);
  if (!active_light) Ef = S(as_double(Ef));   // freeze: strip the self-shading derivative
  S area_leaf = plant::ff16_area_leaf(p.a_l1, p.a_l2, S(focal_h));
  return plant::ff16_net_mass_production_crown_top(p, S(focal_h), area_leaf, Ef);
}

// [[Rcpp::export]]
Rcpp::List self_shading(Rcpp::NumericVector pp, std::vector<double> h,
                        std::vector<double> C, double focal_h) {
  auto s = make_strategy(pp); auto pd = s.prod_pars();
  const double eta = s.pars.eta;

  // (a) reconstructed light vs nothing here (checked in R); return at a few z.
  std::vector<double> zs = {1,3,5,8,12,15,17}, light;
  for (double z : zs) light.push_back(recon_light<double>(z, pd, eta, h, C));

  // (b) focal net + gradient w.r.t. a_l1, self-shading ACTIVE.
  const double Jd = focal_net<double>(pd, focal_h, eta, h, C, true);
  double dJ_active, dJ_frozen;
  { ad::tape_type t; ad_t a=pd.a_l1; t.registerInput(a); t.newRecording();
    auto pa=lift<ad_t>(pd); pa.a_l1=a;
    ad_t J=focal_net<ad_t>(pa, focal_h, eta, h, C, true);
    t.registerOutput(J); xad::derivative(J)=1.0; t.computeAdjoints(); dJ_active=xad::derivative(a); }
  { ad::tape_type t; ad_t a=pd.a_l1; t.registerInput(a); t.newRecording();
    auto pa=lift<ad_t>(pd); pa.a_l1=a;
    ad_t J=focal_net<ad_t>(pa, focal_h, eta, h, C, false);
    t.registerOutput(J); xad::derivative(J)=1.0; t.computeAdjoints(); dJ_frozen=xad::derivative(a); }
  std::vector<double> rel_h={1e-4,1e-5,1e-6,1e-7}, fd;
  for (double rh:rel_h){ double hh=rh*pd.a_l1; auto q1=pd; q1.a_l1+=hh; auto q2=pd; q2.a_l1-=hh;
    fd.push_back((focal_net<double>(q1,focal_h,eta,h,C,true)-focal_net<double>(q2,focal_h,eta,h,C,true))/(2*hh)); }
  return Rcpp::List::create(Rcpp::_["z"]=Rcpp::wrap(zs), Rcpp::_["light"]=Rcpp::wrap(light),
    Rcpp::_["J"]=Jd, Rcpp::_["dJ_active"]=dJ_active, Rcpp::_["dJ_frozen"]=dJ_frozen,
    Rcpp::_["fd"]=Rcpp::wrap(fd), Rcpp::_["rel_h"]=Rcpp::wrap(rel_h));
}')

focal_h <- 6.0
res <- self_shading(pp, h_desc, C_desc, focal_h)

## (a) reconstruction vs live env.
env <- scm$patch$environment
cat("\n(a) reconstructed self-shaded light vs live FF16_Environment:\n")
cat("    z    live_env      recon        abs.err\n")
maxe <- 0
for (k in seq_along(res$z)) {
  le <- env$get_environment_at_height(res$z[k]); maxe <- max(maxe, abs(le-res$light[k]))
  cat(sprintf("  %5.1f  %.8f  %.8f   %.1e\n", res$z[k], le, res$light[k], abs(le-res$light[k])))
}
cat(sprintf("    max abs err = %.1e (limited by the env light spline's interpolation)\n", maxe))

## (b) self-shading gradient.
cat(sprintf("\n(b) focal net production (h=%.1f m) = %.8g\n", focal_h, res$J))
cat("    d(focal net)/d(a_l1) with self-shaded light ACTIVE, vs two-pass FD:\n")
for (k in seq_along(res$rel_h))
  cat(sprintf("      h/a_l1 = %.0e   FD = %.9g   rel.err = %.2e\n",
              res$rel_h[k], res$fd[k], abs(res$fd[k]-res$dJ_active)/abs(res$dJ_active)))
best <- min(abs(res$fd - res$dJ_active) / abs(res$dJ_active))
cat(sprintf("\n      AD (light active) = %.9g   best FD = %.9g   min rel.err = %.2e %s\n",
            res$dJ_active, res$fd[which.min(abs(res$fd-res$dJ_active))], best,
            if (best < 1e-5) "OK" else "** MISMATCH **"))
cat(sprintf("      AD (light frozen) = %.9g   self-shading contribution = %.9g (%.1f%%)\n",
            res$dJ_frozen, res$dJ_active - res$dJ_frozen,
            100*(res$dJ_active-res$dJ_frozen)/res$dJ_active))
stopifnot(best < 1e-5)
cat("\nLive-stand self-shading (active-knot) gradient validated.\n")
