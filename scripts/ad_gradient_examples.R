# Reverse-mode AD trait gradients of FF16 outputs (#472 scope B / #537, Milestone C).
#
# A runnable demonstration of the scalar-templated FF16 kernels added for
# automatic differentiation: every kernel in inst/include/plant/models/
# ff16_production_kernel.h is templated on the scalar type S, so instantiating it
# with an XAD active type and running one reverse sweep yields exact derivatives
# of FF16 outputs w.r.t. traits. Each example below also computes a central finite
# difference of the SAME double computation and checks the two agree.
#
# A precedent first, then four FF16 examples increasing in scope:
#   0. Leaf gradient    -- the groundwork's first piece (TF24's leaf hydraulics):
#                          d(profit)/d(root-collar psi) by forward-mode AD + the
#                          implicit function theorem at the psi_stem->ci root-find.
#                          Pure R via the exposed Leaf class -- no compilation.
#   1. Rate fill        -- d(fecundity_dt)/d(a_p1) of the full demographic rate
#                          vector, and a bit-exact faithfulness check of the
#                          kernel against the live FF16_Strategy::compute_rates.
#   2. Emergent stand   -- d(stand LAI)/d{lma, a_p1} over a multi-cohort
#                          frozen-schedule replay (ff16_replay_cohort).
#   3. Self-shading     -- d(focal net production)/d(a_l1) THROUGH the resident
#                          light profile: the trait reshapes every cohort's leaf
#                          area -> competition -> Beer's law -> an active-value
#                          light spline (odelia basic_interpolator) -> focal light.
#   4. Whole gradient   -- d(net production)/d(ALL 19 production traits) from a
#                          SINGLE reverse sweep -- the headline reverse-mode
#                          advantage made concrete.
#
# Reverse-mode AD is the right tool here: many trait inputs -> one scalar output
# (a calibration objective) is differentiated in a single backward sweep,
# independent of the number of traits. Example 4 shows this directly -- 19
# derivatives from one backward pass, the same cost as one finite difference
# (which would instead need 19+ extra model evaluations).
#
# Requirements: the `plant` package must be INSTALLED from this branch (so its
# installed headers match its compiled .so -- a layout mismatch segfaults), plus
# `odelia` (the XAD tape) and `BH`. The AD path is C++-only; it links the live
# FF16 symbols from plant.so and the reverse-mode tape from odelia.so, exactly as
# the package's tape-linked tests do.
#
# Run from the package root, after `R CMD INSTALL .`:
#   Rscript scripts/ad_gradient_examples.R

suppressMessages({library(Rcpp); library(plant)})

plant_inc  <- system.file("include", package = "plant")
odelia_inc <- system.file("include", package = "odelia")
bh_inc     <- system.file("include", package = "BH")
plant_so   <- system.file("libs", "plant.so",  package = "plant")
odelia_so  <- system.file("libs", "odelia.so", package = "odelia")

ok <- nzchar(plant_inc) && nzchar(odelia_inc) && nzchar(bh_inc) &&
      file.exists(plant_so) && file.exists(odelia_so)
if (!ok) {
  stop("Need plant (installed from this branch), odelia, and BH available; ",
       "and plant.so + odelia.so on disk. Install with `R CMD INSTALL .` first.")
}

Sys.setenv(PKG_CPPFLAGS = paste(paste0("-I", shQuote(plant_inc)),
                                paste0("-I", shQuote(odelia_inc)),
                                paste0("-I", shQuote(bh_inc))))
Sys.setenv(PKG_LIBS = paste(shQuote(normalizePath(plant_so)),
                            shQuote(normalizePath(odelia_so))))

Rcpp::sourceCpp(code = '
#include <Rcpp.h>
#include <optional>
#include <vector>
#include <cmath>
#include <XAD/XAD.hpp>
#include <odelia/interpolator.hpp>
#include <plant/models/ff16_strategy.h>
#include <plant/models/ff16_environment.h>
#include <plant/individual.h>
using ad   = xad::adj<double>;
using ad_t = ad::active_type;
namespace oi = odelia::interpolator;
// [[Rcpp::plugins(cpp20)]]

// Lift a prepared FF16ProdPars<double> to <ad_t> (caller then registers the one
// trait of interest as a tape input and overwrites it).
static plant::FF16ProdPars<ad_t> lift(const plant::FF16ProdPars<double>& d) {
  plant::FF16ProdPars<ad_t> p;
  p.lma=d.lma;p.rho=d.rho;p.theta=d.theta;p.a_b1=d.a_b1;p.a_r1=d.a_r1;p.eta_c=d.eta_c;
  p.a_p1=d.a_p1;p.a_p2=d.a_p2;p.r_l=d.r_l;p.r_s=d.r_s;p.r_b=d.r_b;p.r_r=d.r_r;
  p.k_l=d.k_l;p.k_b=d.k_b;p.k_s=d.k_s;p.k_r=d.k_r;p.a_bio=d.a_bio;p.a_y=d.a_y;
  p.a_l1=d.a_l1;p.a_l2=d.a_l2;p.a_f1=d.a_f1;p.a_f2=d.a_f2;p.hmat=d.hmat;
  p.omega=d.omega;p.a_f3=d.a_f3;p.d_I=d.d_I;p.a_dG1=d.a_dG1;p.a_dG2=d.a_dG2;
  return p;
}
static double as_double(double v)        { return v; }
static double as_double(const ad_t& v)   { return xad::value(v); }

// ---- Example 1: demographic rate fill -------------------------------------
// Faithfulness (live crown-top compute_rates vs the kernel) + d(fecundity_dt)/d(a_p1).
// [[Rcpp::export]]
Rcpp::List ex1_rates(double height, double light_E) {
  plant::FF16_Strategy s; s.control.shading_model = "crown-centre"; s.prepare_strategy();
  auto sp = plant::make_strategy_ptr(s);
  plant::Individual<plant::FF16_Strategy, plant::FF16_Environment> ind(sp);
  ind.set_state("height", height);
  plant::FF16_Environment env; env.set_fixed_environment(light_E, 1e4);
  ind.compute_rates(env);

  auto pd = s.prod_pars();
  auto r  = plant::ff16_compute_rates_crown_top<double>(pd, height, light_E, true);
  Rcpp::NumericVector live = Rcpp::NumericVector::create(
    ind.rate("height"), ind.rate("fecundity"), ind.rate("area_heartwood"),
    ind.rate("mass_heartwood"), ind.rate("mortality"));
  Rcpp::NumericVector kern = Rcpp::NumericVector::create(
    r.height_dt, r.fecundity_dt, r.area_heartwood_dt, r.mass_heartwood_dt, r.mortality_dt);

  double dfec;
  { ad::tape_type tape; ad_t a_p1 = pd.a_p1; tape.registerInput(a_p1); tape.newRecording();
    auto p = lift(pd); p.a_p1 = a_p1;
    ad_t f = plant::ff16_compute_rates_crown_top<ad_t>(p, ad_t(height), ad_t(light_E), true).fecundity_dt;
    tape.registerOutput(f); xad::derivative(f) = 1.0; tape.computeAdjoints();
    dfec = xad::derivative(a_p1); }
  auto kf = [&](double v){ auto q = pd; q.a_p1 = v;
    return plant::ff16_compute_rates_crown_top<double>(q, height, light_E, true).fecundity_dt; };
  double h = 1e-4 * pd.a_p1, dfec_fd = (kf(pd.a_p1+h) - kf(pd.a_p1-h)) / (2*h);
  return Rcpp::List::create(Rcpp::_["live"]=live, Rcpp::_["kernel"]=kern,
                            Rcpp::_["d_ap1_ad"]=dfec, Rcpp::_["d_ap1_fd"]=dfec_fd);
}

// ---- Example 2: emergent multi-cohort stand LAI ---------------------------
template <typename S>
S stand_LAI(const plant::FF16ProdPars<S>& p, S h0, double dt,
            const std::vector<double>& light, const std::vector<int>& intro,
            const std::vector<double>& w) {
  S LAI = S(0.0);
  for (size_t i = 0; i < intro.size(); ++i) {
    plant::FF16State<S> y{h0, S(0), S(0), S(0), S(0)};
    y = plant::ff16_replay_cohort<S>(p, y, dt, light, (size_t)intro[i], true);
    LAI += w[i] * plant::ff16_area_leaf(p.a_l1, p.a_l2, y.height);
  }
  return LAI;
}
// [[Rcpp::export]]
Rcpp::NumericVector ex2_stand(double h0, double dt) {
  plant::FF16_Strategy s; s.control.shading_model = "crown-centre"; s.prepare_strategy();
  auto pd = s.prod_pars();
  const int M = 120; std::vector<double> light(M);
  for (int t = 0; t < M; ++t) light[t] = 0.97 - 0.5 * ((double)t / (M-1));
  std::vector<int> intro = {0,15,30,45,60,75};
  std::vector<double> w = {1.0,0.85,0.7,0.55,0.4,0.25};
  double Ld = stand_LAI<double>(pd, h0, dt, light, intro, w);
  double dlma, dap1;
  { ad::tape_type tp; ad_t lma=pd.lma; tp.registerInput(lma); tp.newRecording();
    auto p=lift(pd); p.lma=lma; ad_t L=stand_LAI<ad_t>(p,ad_t(h0),dt,light,intro,w);
    tp.registerOutput(L); xad::derivative(L)=1.0; tp.computeAdjoints(); dlma=xad::derivative(lma); }
  { ad::tape_type tp; ad_t ap1=pd.a_p1; tp.registerInput(ap1); tp.newRecording();
    auto p=lift(pd); p.a_p1=ap1; ad_t L=stand_LAI<ad_t>(p,ad_t(h0),dt,light,intro,w);
    tp.registerOutput(L); xad::derivative(L)=1.0; tp.computeAdjoints(); dap1=xad::derivative(ap1); }
  auto kl=[&](double v){auto q=pd;q.lma=v;return stand_LAI<double>(q,h0,dt,light,intro,w);};
  auto ka=[&](double v){auto q=pd;q.a_p1=v;return stand_LAI<double>(q,h0,dt,light,intro,w);};
  double hl=1e-5*pd.lma, ha=1e-5*pd.a_p1;
  return Rcpp::NumericVector::create(Ld, dlma,(kl(pd.lma+hl)-kl(pd.lma-hl))/(2*hl),
                                          dap1,(ka(pd.a_p1+ha)-ka(pd.a_p1-ha))/(2*ha));
}

// ---- Example 3: full self-shading gradient --------------------------------
template <typename S>
S focal_net(const plant::FF16ProdPars<S>& p, double k_I, double eta,
            const std::vector<double>& h, const std::vector<double>& dens,
            const std::vector<double>& zk, int focal) {
  std::vector<S> Ek(zk.size());
  for (size_t j = 0; j < zk.size(); ++j)
    Ek[j] = plant::ff16_resident_light_at<S>(zk[j], p.a_l1, p.a_l2, k_I, eta, h, dens);
  oi::basic_interpolator<S> light; light.init(zk, Ek);     // frozen x, active y
  double zf = h[focal] * as_double(p.eta_c);               // focal crown height
  S Ef  = light.eval(zf);
  S aLf = plant::ff16_area_leaf(p.a_l1, p.a_l2, S(h[focal]));
  return plant::ff16_net_mass_production_crown_top(p, S(h[focal]), aLf, Ef);
}
// [[Rcpp::export]]
Rcpp::NumericVector ex3_selfshade() {
  plant::FF16_Strategy s; s.control.shading_model = "crown-centre"; s.prepare_strategy();
  auto pd = s.prod_pars(); double k_I = s.pars.k_I, eta = s.pars.eta;
  std::vector<double> h = {12,9,6,4,2.5}, dens = {0.2,0.4,0.7,1.0,1.5};
  std::vector<double> zk; for (int j = 0; j <= 40; ++j) zk.push_back(12.0*j/40.0);
  int focal = 2;
  double Jd = focal_net<double>(pd, k_I, eta, h, dens, zk, focal);
  double dJ;
  { ad::tape_type tp; ad_t a_l1=pd.a_l1; tp.registerInput(a_l1); tp.newRecording();
    auto p=lift(pd); p.a_l1=a_l1; ad_t J=focal_net<ad_t>(p,k_I,eta,h,dens,zk,focal);
    tp.registerOutput(J); xad::derivative(J)=1.0; tp.computeAdjoints(); dJ=xad::derivative(a_l1); }
  auto kf=[&](double v){auto q=pd;q.a_l1=v;return focal_net<double>(q,k_I,eta,h,dens,zk,focal);};
  double hh=1e-6*pd.a_l1;
  return Rcpp::NumericVector::create(Jd, dJ, (kf(pd.a_l1+hh)-kf(pd.a_l1-hh))/(2*hh));
}

// ---- Example 4: the whole trait-gradient vector in ONE reverse sweep -------
// FF16 net production (crown-top) as the scalar objective; differentiate it
// w.r.t. all 19 production-relevant traits at once. a_l1/a_l2 also flow through
// area_leaf, so the gradient covers allometry as well as physiology.
template <typename S>
S ff16_netprod(const plant::FF16ProdPars<S>& p, double height, double light_E) {
  S al = plant::ff16_area_leaf(p.a_l1, p.a_l2, S(height));
  return plant::ff16_net_mass_production_crown_top(p, S(height), al, S(light_E));
}
// [[Rcpp::export]]
Rcpp::List ex4_all_traits(double height, double light_E) {
  plant::FF16_Strategy s; s.control.shading_model = "crown-centre"; s.prepare_strategy();
  auto pd = s.prod_pars();

  // ONE tape, one forward eval, one backward sweep -> every trait derivative.
  ad::tape_type tape;
  auto p = lift(pd);
  std::vector<ad_t*> in = {&p.lma,&p.rho,&p.theta,&p.a_b1,&p.a_r1,&p.a_p1,&p.a_p2,
    &p.r_l,&p.r_s,&p.r_b,&p.r_r,&p.k_l,&p.k_b,&p.k_s,&p.k_r,&p.a_bio,&p.a_y,&p.a_l1,&p.a_l2};
  for (auto* x : in) tape.registerInput(*x);
  tape.newRecording();
  ad_t J = ff16_netprod<ad_t>(p, height, light_E);
  tape.registerOutput(J); xad::derivative(J) = 1.0; tape.computeAdjoints();
  Rcpp::NumericVector grad(in.size());
  for (size_t i = 0; i < in.size(); ++i) grad[i] = xad::derivative(*in[i]);

  // Per-trait central FD for comparison (one extra pair of evals per trait --
  // the cost reverse mode avoids).
  auto dd = pd;
  std::vector<double*> dp = {&dd.lma,&dd.rho,&dd.theta,&dd.a_b1,&dd.a_r1,&dd.a_p1,&dd.a_p2,
    &dd.r_l,&dd.r_s,&dd.r_b,&dd.r_r,&dd.k_l,&dd.k_b,&dd.k_s,&dd.k_r,&dd.a_bio,&dd.a_y,&dd.a_l1,&dd.a_l2};
  Rcpp::NumericVector fd(dp.size());
  for (size_t i = 0; i < dp.size(); ++i) {
    double b = *dp[i], hh = 1e-6 * std::max(1.0, std::abs(b));
    *dp[i] = b + hh; double jp = ff16_netprod<double>(dd, height, light_E);
    *dp[i] = b - hh; double jm = ff16_netprod<double>(dd, height, light_E);
    *dp[i] = b; fd[i] = (jp - jm) / (2 * hh);
  }
  return Rcpp::List::create(Rcpp::_["J"]=xad::value(J), Rcpp::_["grad"]=grad, Rcpp::_["fd"]=fd);
}')

rel <- function(a, b) abs(a - b) / pmax(abs(b), 1e-30)
chk <- function(name, ad, fd, tol = 1e-6) {
  r <- rel(ad, fd)
  cat(sprintf("  %-28s AD = % .8g   FD = % .8g   rel.err = %.1e %s\n",
              name, ad, fd, r, if (r < tol) "OK" else "** MISMATCH **"))
  stopifnot(r < tol)
}

cat("== Example 0 (precedent): leaf-level gradient -- forward-mode AD + IFT ==\n")
# TF24's leaf hydraulics, the first exact AD gradient in plant (#531/#539). All
# methods used here are exposed on the Leaf R class, so this needs no compilation.
local({
  root_c <- 2.65; root_b <- 1.29; theta <- 0.000157; h <- 5
  l <- Leaf(vcmax_25 = 100, jmax_25 = 100 * 167, c = 2.04, b = 3, psi_crit = 5,
            root_c = root_c, root_b = root_b,
            root_psi_crit = root_b * (log(1 / 0.05))^(1 / root_c), beta2 = 1,
            hk_s = 75, a = 0.3, curv_fact_elec_trans = 0.7, curv_fact_colim = 0.99,
            GSS_tol_abs = 1e-8, vulnerability_curve_ncontrol = 100,
            ci_abs_tol = 1e-6, ci_niter = 1000, g1_TF24 = 46.32995,
            beta_R_H = 3.4e3, beta_R_V = 9.4e4)
  l$set_physiology(area_leaf = 0.05, mass_root_prop = 1, rho = 608, a_bio = 0.0245,
                   PPFD = 900, psi_soil = 2, soil_depth = 1,
                   leaf_specific_conductance_max = theta / h, atm_vpd = 2, ca = 40,
                   sapwood_volume_per_leaf_area = theta * h, leaf_temp = 25,
                   atm_o2_kpa = 21, atm_kpa = 101.3)
  l$find_root_collar_psi()
  opt <- -l$root_collar_psi_   # operating root-collar potential (positive magnitude)
  # profit(psi) via evaluate_root_collar_psi; the exact gradient via AD + IFT.
  fd <- function(psi, e = 1e-5)
    (l$evaluate_root_collar_psi(psi + e) - l$evaluate_root_collar_psi(psi - e)) / (2 * e)
  for (psi in c(opt + 0.1, opt + 0.2)) {        # strictly-interior, feasible points
    l$evaluate_root_collar_psi(psi)
    if (abs(-l$root_collar_psi_ - psi) > 1e-8) next   # skip if clamped to the boundary
    chk(sprintf("d(profit)/d(psi) @ %.2f", psi),
        l$dprofit_droot_collar_psi(psi), fd(psi), tol = 1e-4)
  }
})

cat("\n== Example 1: FF16 demographic rate fill (height=3.7 m, light=0.92) ==\n")
e1 <- ex1_rates(3.7, 0.92)
rn <- c("height_dt","fecundity_dt","area_heartwood_dt","mass_heartwood_dt","mortality_dt")
cat("  faithfulness of ff16_compute_rates_crown_top vs live FF16_Strategy::compute_rates:\n")
for (i in seq_along(rn))
  cat(sprintf("    %-20s live = % .10g   kernel = % .10g   rel = %.1e\n",
              rn[i], e1$live[i], e1$kernel[i], rel(e1$live[i], e1$kernel[i])))
stopifnot(max(rel(e1$live, e1$kernel)) < 1e-12)
chk("d(fecundity_dt)/d(a_p1)", e1$d_ap1_ad, e1$d_ap1_fd)

cat("\n== Example 2: emergent multi-cohort stand LAI (6-cohort frozen-schedule replay) ==\n")
e2 <- ex2_stand(0.4, 0.05)
cat(sprintf("  stand LAI = %.8g\n", e2[1]))
chk("d(LAI)/d(lma)",  e2[2], e2[3])
chk("d(LAI)/d(a_p1)", e2[4], e2[5])

cat("\n== Example 3: full self-shading gradient (resident light responds to trait) ==\n")
e3 <- ex3_selfshade()
cat(sprintf("  focal net production = %.8g\n", e3[1]))
chk("d(focal net)/d(a_l1)", e3[2], e3[3])

cat("\n== Example 4: whole trait-gradient vector in ONE reverse sweep ==\n")
e4 <- ex4_all_traits(3.7, 0.92)
traits <- c("lma","rho","theta","a_b1","a_r1","a_p1","a_p2","r_l","r_s","r_b",
            "r_r","k_l","k_b","k_s","k_r","a_bio","a_y","a_l1","a_l2")
cat(sprintf("  net production = %.8g  (%d trait derivatives from one backward sweep)\n",
            e4$J, length(e4$grad)))
for (i in seq_along(traits)) chk(sprintf("d(net)/d(%s)", traits[i]), e4$grad[i], e4$fd[i])

cat("\nAll reverse-mode AD gradients match finite differences. ",
    "Kernels: inst/include/plant/models/ff16_production_kernel.h\n", sep = "")


