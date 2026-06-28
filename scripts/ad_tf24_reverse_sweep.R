# TF24 NET-PRODUCTION whole-gradient in ONE reverse sweep (#472 scope B, Phase
# F1-full) -- the headline reverse-mode advantage, made concrete for TF24.
#
# The four per-class TF24 scripts (ad_tf24_net_gradient.R [vcmax], ad_tf24_photo_
# gradient.R, ad_tf24_hydraulic_gradient.R, ad_tf24_mass_gradient.R) each validate
# ONE trait's d(net_mass_production_dt)/d(trait) by FORWARD-mode AD vs a live FD.
# This script assembles all 27 into a SINGLE reverse-mode pass:
#
#   net = a_bio*a_y*(profit*area_leaf*conv - respiration) - turnover
#
# carried through the committed scalar-templated kernel (tf24_net_mass_production,
# tf24_production_kernel.h). The leaf optimisation profit*(theta) is NOT taped
# (it nests a root-find/optimiser); instead its trait sensitivities -- the
# validated Leaf::dprofit_d* numbers -- are INJECTED first-order into an active
# `profit` (the #539 IFT / FF16 height_0 injection pattern):
#
#   profit_ad = profit_v + sum_k dprofit_dk * (trait_k - trait_k_v)
#
# over the 12 profit-coupled traits (10 leaf + theta via k_max + a_r1 via E_up_),
# while a_l1/a_l2 drive area_leaf actively and the 13 pure mass-cascade traits flow
# through the cascade. ONE backward pass then yields the FULL 27-vector
# d(net)/d(theta_k) -- at the SAME cost as one trait (the reverse tape size is
# input-count-independent), whereas the per-trait forward scripts need 27 sweeps.
#
# Validation (the contract): the reverse 27-vector is checked against (a) the
# per-trait FORWARD-mode value built from the IDENTICAL injection (reverse ==
# forward, the tape-machinery check, ~1e-10) and (b) a live two-sided FD of
# TF24_Strategy::net_mass_production_dt (the ground truth, ~1e-5..1e-8; b/c use a
# dense vulnerability spline so the FD resolves the exact continuous dS/dt).
#
# Run from the package root after `R CMD INSTALL .`:
#   Rscript scripts/ad_tf24_reverse_sweep.R

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
#include <vector>
#include <string>
#include <cmath>
#include <XAD/XAD.hpp>
#include <plant/models/tf24_strategy.h>
#include <plant/models/tf24_environment.h>
#include <plant/models/tf24_production_kernel.h>
// [[Rcpp::plugins(cpp20)]]
using radj = xad::adj<double>;   using ad_t  = radj::active_type;   // reverse
using rfwd = xad::fwd<double>;   using fad_t = rfwd::active_type;   // forward

// The 27 net-production traits, in a fixed order. 10 leaf (profit-coupled),
// 13 pure mass-cascade, 2 area_leaf-active, 2 leaf-coupled cascade.
static const std::vector<std::string> TRAITS = {
  "vcmax_25","g1_TF24","beta2","K_s","b","c","jmax_25","a","curv_elec","curv_colim",
  "lma","rho","a_b1","r_l","r_b","r_s","r_r","k_l","k_b","k_s","k_r","a_bio","a_y",
  "a_l1","a_l2","theta","a_r1"};

// Configure a crown-centre TF24 strategy with a tight collar optimum (envelope
// theorem exact) and a dense vulnerability spline (so the b/c FD resolves the
// exact dS/dt the AD computes). One named trait optionally overridden.
static plant::TF24_Strategy make_strategy(const std::string& over = "", double v = 0) {
  plant::TF24_Strategy s;
  s.control.shading_model = "crown-centre";
  s.control.GSS_tol_abs = 1e-9;
  s.control.vulnerability_curve_ncontrol = 2000;
  if (over == "g1_TF24")          s.g1_TF24 = v;
  else if (over == "curv_elec")   s.pars.curv_fact_elec_trans = v;
  else if (over == "curv_colim")  s.pars.curv_fact_colim = v;
  else if (over == "vcmax_25")    s.pars.vcmax_25 = v;
  else if (over == "beta2")       s.pars.beta2 = v;
  else if (over == "K_s")         s.pars.K_s = v;
  else if (over == "b")           s.pars.b = v;
  else if (over == "c")           s.pars.c = v;
  else if (over == "jmax_25")     s.pars.jmax_25 = v;
  else if (over == "a")           s.pars.a = v;
  else if (over == "lma")         s.pars.lma = v;
  else if (over == "rho")         s.pars.rho = v;
  else if (over == "a_b1")        s.pars.a_b1 = v;
  else if (over == "r_l")         s.pars.r_l = v;
  else if (over == "r_b")         s.pars.r_b = v;
  else if (over == "r_s")         s.pars.r_s = v;
  else if (over == "r_r")         s.pars.r_r = v;
  else if (over == "k_l")         s.pars.k_l = v;
  else if (over == "k_b")         s.pars.k_b = v;
  else if (over == "k_s")         s.pars.k_s = v;
  else if (over == "k_r")         s.pars.k_r = v;
  else if (over == "a_bio")       s.pars.a_bio = v;
  else if (over == "a_y")         s.pars.a_y = v;
  else if (over == "a_l1")        s.pars.a_l1 = v;
  else if (over == "a_l2")        s.pars.a_l2 = v;
  else if (over == "theta")       s.pars.theta = v;
  else if (over == "a_r1")        s.pars.a_r1 = v;
  else if (!over.empty())         Rcpp::stop("unknown trait " + over);
  s.prepare_strategy();
  return s;
}

// The live value of a trait on a configured strategy.
static double trait_value(const plant::TF24_Strategy& s, const std::string& t) {
  if (t == "g1_TF24")     return s.g1_TF24;
  if (t == "curv_elec")   return s.pars.curv_fact_elec_trans;
  if (t == "curv_colim")  return s.pars.curv_fact_colim;
  const auto& p = s.pars;
  if (t=="vcmax_25")return p.vcmax_25; if(t=="beta2")return p.beta2; if(t=="K_s")return p.K_s;
  if (t=="b")return p.b; if(t=="c")return p.c; if(t=="jmax_25")return p.jmax_25; if(t=="a")return p.a;
  if (t=="lma")return p.lma; if(t=="rho")return p.rho; if(t=="a_b1")return p.a_b1;
  if (t=="r_l")return p.r_l; if(t=="r_b")return p.r_b; if(t=="r_s")return p.r_s; if(t=="r_r")return p.r_r;
  if (t=="k_l")return p.k_l; if(t=="k_b")return p.k_b; if(t=="k_s")return p.k_s; if(t=="k_r")return p.k_r;
  if (t=="a_bio")return p.a_bio; if(t=="a_y")return p.a_y; if(t=="a_l1")return p.a_l1;
  if (t=="a_l2")return p.a_l2; if(t=="theta")return p.theta; if(t=="a_r1")return p.a_r1;
  Rcpp::stop("unknown trait " + t);
}

// Live net production at the operating point (height, light), trait overridden.
static double net_live(const std::string& t, double v, double light, double h) {
  plant::TF24_Strategy s = make_strategy(t, v);
  plant::TF24_Environment e; e.set_fixed_environment(light, 1e4);
  return s.net_mass_production_dt(e, h, s.area_leaf(h), 1.0 / h);
}

// Build a TF24ProdPars<S> from the live (double) prod_pars, with the per-trait S
// variables substituted for the cascade fields. The 13 pure-cascade, the 2 area
// and the 2 leaf-coupled cascade traits live in prod_pars; the 10 leaf traits and
// the demographic params do not (net does not depend on the latter).
template <typename S>
static plant::TF24ProdPars<S> make_prodpars(const plant::TF24ProdPars<double>& d,
                                            const std::vector<S>& tr) {
  // index helper into TRAITS
  auto IX = [](const std::string& n){ for (std::size_t i=0;i<TRAITS.size();++i) if(TRAITS[i]==n) return i; return (std::size_t)-1; };
  plant::TF24ProdPars<S> p;
  p.lma=tr[IX("lma")]; p.rho=tr[IX("rho")]; p.theta=tr[IX("theta")]; p.a_b1=tr[IX("a_b1")];
  p.a_r1=tr[IX("a_r1")]; p.eta_c=S(d.eta_c);
  p.r_l=tr[IX("r_l")]; p.r_s=tr[IX("r_s")]; p.r_b=tr[IX("r_b")]; p.r_r=tr[IX("r_r")];
  p.k_l=tr[IX("k_l")]; p.k_b=tr[IX("k_b")]; p.k_s=tr[IX("k_s")]; p.k_r=tr[IX("k_r")];
  p.a_bio=tr[IX("a_bio")]; p.a_y=tr[IX("a_y")];
  p.a_l1=tr[IX("a_l1")]; p.a_l2=tr[IX("a_l2")];
  // demographic params irrelevant to net; set to the (frozen) double values.
  p.a_f1=S(d.a_f1); p.a_f2=S(d.a_f2); p.hmat=S(d.hmat);
  p.omega=S(d.omega); p.a_f3=S(d.a_f3);
  p.d_I=S(d.d_I); p.a_dG1=S(d.a_dG1); p.a_dG2=S(d.a_dG2);
  return p;
}

// Net production as a function of the 27 active trait scalars, given the frozen
// leaf optimisation (profit value + injected dprofit/dtrait sensitivities). This
// is the single expression both the reverse and the forward sweeps differentiate.
template <typename S>
static S net_expr(const std::vector<S>& tr, double h,
                  const plant::TF24ProdPars<double>& pd,
                  double profit_v, const std::vector<double>& dprofit) {
  auto IX = [](const std::string& n){ for (std::size_t i=0;i<TRAITS.size();++i) if(TRAITS[i]==n) return i; return (std::size_t)-1; };
  plant::TF24ProdPars<S> p = make_prodpars<S>(pd, tr);
  // area_leaf active in a_l1, a_l2.
  S area_leaf = plant::tf24_area_leaf<S>(tr[IX("a_l1")], tr[IX("a_l2")], S(h));
  // profit: value + first-order injection of the leaf sensitivities.
  S profit = S(profit_v);
  for (std::size_t i = 0; i < TRAITS.size(); ++i)
    if (dprofit[i] != 0.0) profit += S(dprofit[i]) * (tr[i] - S(xad::value(tr[i])));
  return plant::tf24_net_mass_production<S>(p, S(h), area_leaf, profit);
}

// [[Rcpp::export]]
Rcpp::List tf24_net_reverse_sweep(double light, double h) {
  plant::TF24_Strategy s = make_strategy();
  plant::TF24_Environment env; env.set_fixed_environment(light, 1e4);
  const double al  = s.area_leaf(h);
  const double net = s.net_mass_production_dt(env, h, al, 1.0 / h);
  const double opt = -s.leaf.root_collar_psi_;
  const double conv = plant::tf24_assimilation_conv;
  const double scale = s.pars.a_bio * s.pars.a_y * al * conv;   // d(net)/d(profit)
  const double profit_v = s.leaf.profit_;
  const plant::TF24ProdPars<double> pd = s.prod_pars();
  const std::size_t N = TRAITS.size();
  auto IX = [](const std::string& n){ for (std::size_t i=0;i<TRAITS.size();++i) if(TRAITS[i]==n) return i; return (std::size_t)-1; };

  // Leaf-profit sensitivity d(profit)/d(trait) for the 12 profit-coupled traits
  // (the validated Leaf::dprofit_d* numbers); 0 for the rest. The reverse/forward
  // sweeps multiply these by scale = d(net)/d(profit) through the kernel.
  std::vector<double> dprofit(N, 0.0);
  dprofit[IX("vcmax_25")]   = s.leaf.dprofit_dvcmax25(opt);
  dprofit[IX("g1_TF24")]    = s.leaf.dprofit_dg1_TF24(opt);
  dprofit[IX("beta2")]      = s.leaf.dprofit_dbeta2(opt);
  dprofit[IX("b")]          = s.leaf.dprofit_db(opt);
  dprofit[IX("c")]          = s.leaf.dprofit_dc(opt);
  dprofit[IX("jmax_25")]    = s.leaf.dprofit_djmax25(opt);
  dprofit[IX("a")]          = s.leaf.dprofit_da(opt);
  dprofit[IX("curv_elec")]  = s.leaf.dprofit_dcurv_elec(opt);
  dprofit[IX("curv_colim")] = s.leaf.dprofit_dcurv_colim(opt);
  const double kmax = s.leaf.leaf_specific_conductance_max_;
  dprofit[IX("K_s")]   = s.leaf.dprofit_dkmax(opt) * (kmax / s.pars.K_s);
  dprofit[IX("theta")] = s.leaf.dprofit_dkmax(opt) * (kmax / s.pars.theta);
  dprofit[IX("a_r1")]  = s.leaf.dprofit_dEup(opt)  * (s.leaf.E_up_ / s.pars.a_r1);

  std::vector<double> v0(N);
  for (std::size_t i = 0; i < N; ++i) v0[i] = trait_value(s, TRAITS[i]);

  // ---- ONE reverse sweep -> the whole 27-vector. -------------------------
  std::vector<ad_t> tr(N);
  for (std::size_t i = 0; i < N; ++i) tr[i] = v0[i];
  radj::tape_type tape;
  for (auto& x : tr) tape.registerInput(x);
  tape.newRecording();
  ad_t net_ad = net_expr<ad_t>(tr, h, pd, profit_v, dprofit);
  tape.registerOutput(net_ad);
  xad::derivative(net_ad) = 1.0;
  tape.computeAdjoints();
  std::vector<double> rev(N);
  for (std::size_t i = 0; i < N; ++i) rev[i] = xad::derivative(tr[i]);

  // ---- Per-trait forward sweep (same injected expression) ----------------
  std::vector<double> fwd(N);
  for (std::size_t j = 0; j < N; ++j) {
    std::vector<fad_t> trf(N);
    for (std::size_t i = 0; i < N; ++i) trf[i] = v0[i];
    xad::derivative(trf[j]) = 1.0;
    fad_t nf = net_expr<fad_t>(trf, h, pd, profit_v, dprofit);
    fwd[j] = xad::derivative(nf);
  }

  // ---- Live two-sided FD of net_mass_production_dt (ground truth) ---------
  // The FD step that resolves each trait differs by orders of magnitude: the 12
  // profit-coupled traits go through the leaf hydraulic optimisation, whose
  // root-find has an absolute noise floor (~1e-9 in net), so a too-small step has
  // noise swamp the signal; yet a trait like jmax_25 (nonlinear through electron
  // transport) is truncation-limited and wants a SMALL step; the 15 pure cascade
  // traits give a smooth closed-form net (profit frozen). No single step works
  // for all -- exactly the asymmetry reverse-mode AD sidesteps (the AD needs no
  // step tuning: reverse == forward to machine eps for ALL 27). So the FD here
  // uses a ROBUST PLATEAU picker (non-circular w.r.t. the AD): evaluate the
  // central difference over a ladder of relative steps and report the value on
  // the most self-consistent (smallest adjacent-difference) rung.
  const std::vector<double> rel_steps = {1e-4, 3e-5, 1e-5, 3e-6, 1e-6, 3e-7};
  std::vector<double> fd(N);
  for (std::size_t i = 0; i < N; ++i) {
    const double b0 = v0[i], scl = (std::abs(b0) > 0 ? std::abs(b0) : 1.0);
    std::vector<double> cand(rel_steps.size());
    for (std::size_t k = 0; k < rel_steps.size(); ++k) {
      const double step = rel_steps[k] * scl;
      cand[k] = (net_live(TRAITS[i], b0 + step, light, h) -
                 net_live(TRAITS[i], b0 - step, light, h)) / (2 * step);
    }
    std::size_t best = 0; double best_gap = std::abs(cand[1] - cand[0]);
    for (std::size_t k = 1; k + 1 < cand.size(); ++k) {
      const double gap = std::abs(cand[k + 1] - cand[k]);
      if (gap < best_gap) { best_gap = gap; best = k; }
    }
    fd[i] = cand[best];   // the plateau value
  }

  return Rcpp::List::create(
    Rcpp::_["trait"] = TRAITS, Rcpp::_["net"] = net,
    Rcpp::_["reverse"] = rev, Rcpp::_["forward"] = fwd, Rcpp::_["fd"] = fd,
    Rcpp::_["dprofit_dcollar"] = s.leaf.dprofit_droot_collar_psi(opt),
    Rcpp::_["profit"] = profit_v);
}')

cat("TF24 d(net_mass_production_dt)/d(trait): ONE reverse sweep vs forward + live FD\n")
cat("(crown-centre, light=0.7, h=12; interior optimum)\n\n")
r <- tf24_net_reverse_sweep(0.7, 12.0)
stopifnot(r$profit > 0, abs(r$dprofit_dcollar) < 1e-2)
cat(sprintf("net = %.6g  (profit=%.4g, interior dp/dcollar=%.1e)\n\n",
            r$net, r$profit, r$dprofit_dcollar))

tr <- r$trait
df <- data.frame(trait = tr, reverse = r$reverse, forward = r$forward, fd = r$fd)
# reverse vs forward: identical expression differentiated two ways -> ~machine eps
df$rel_rf <- abs(df$reverse - df$forward) / pmax(abs(df$forward), 1e-30)
# reverse vs live FD: the ground-truth contract.
df$rel_fd <- abs(df$reverse - df$fd) / pmax(abs(df$fd), 1e-30)
df$abs_fd <- abs(df$reverse - df$fd)

# reverse vs forward is the PRIMARY contract (identical injected expression, two
# AD modes -> machine eps). The live FD is the ground-truth sanity check; 26/27
# match to ~1e-7, and the lone abs-floor case is c, the vulnerability-SHAPE trait:
# the AD computes the EXACT continuous d(transpiration integral)/dc while the FD
# rebuilds the discretised vulnerability spline, so the residual is the spline
# interpolation error (~4e-6, shrinks with vulnerability_curve_ncontrol -- the
# convergence that proves the AD is the exact derivative; see the dedicated
# ad_tf24_hydraulic_gradient.R). Hence an abs-OR-rel criterion.
ok_rf <- all(df$rel_rf < 1e-9)
ok_fd <- all(df$rel_fd < 1e-5 | df$abs_fd < 1e-5)
cat(sprintf("%-11s %15s %15s %15s %9s %9s\n",
            "trait","reverse","forward","live FD","rel(r-f)","rel(r-FD)"))
for (i in seq_along(tr)) {
  pass <- (df$rel_fd[i] < 1e-5 || df$abs_fd[i] < 1e-5) && df$rel_rf[i] < 1e-9
  cat(sprintf("%-11s %15.7g %15.7g %15.7g %9.1e %9.1e %s\n",
              df$trait[i], df$reverse[i], df$forward[i], df$fd[i],
              df$rel_rf[i], df$rel_fd[i], if (pass) "OK" else "** MISMATCH **"))
}
cat(sprintf("\nreverse == forward (all 27): %s   reverse == live FD (all 27): %s\n",
            if (ok_rf) "YES" else "NO", if (ok_fd) "YES" else "NO"))
stopifnot(ok_rf, ok_fd)
cat("\nONE reverse sweep reproduced all 27 per-trait forward gradients AND the live FD.\n")
cat("Reverse-tape cost is input-count-independent: 27 derivatives for the price of 1.\n")
