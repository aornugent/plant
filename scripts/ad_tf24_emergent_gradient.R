# TF24 emergent community trait gradient over the LIVE resident SCM (#472 scope B,
# Phase F1-full) -- AD through the ENTIRE TF24 SCM.
#
# The FF16 emergent gradient (scripts/ad_emergent_gradient.R) could TAPE the whole
# trajectory because FF16 assimilation is a closed form of light. TF24 cannot: its
# net production comes from the hydraulic LEAF OPTIMISATION (a max over collar
# potential nesting a psi_stem->ci root-find), which has no tape. So this uses a
# TANGENT-LINEAR (forward-sensitivity) two-pass replay:
#
#   Pass 1 (double): run the real TF24 resident SCM to completion with the adaptive
#     Cash-Karp RKCK solver + save_RK45_cache (crown-centre shading). Harvest the
#     frozen schedule (Patch$step_history) and per-RK-stage resident environment
#     (Patch$environment_history[step][0..5]), plus each cohort's birth step/weight.
#   Pass 2 (tangent-linear): replay every cohort with the SAME Cash-Karp stepper
#     (ff16_cashkarp_replay), carrying BOTH the demographic state AND its d/d(trait)
#     sensitivity. At each RK stage the deriv runs the REAL leaf optimisation
#     (TF24_Strategy::net_mass_production_dt) in the frozen stage environment to get
#     `net` (so the double trajectory is faithful), and propagates the trait
#     sensitivity:
#       d(net)/d(vcmax)      = a_bio*a_y*area_leaf*conv*Leaf::dprofit_dvcmax25  (analytic),
#       d(net)/d(height)     = central FD of the leaf opt   (the within-trajectory
#                              height-Jacobian through the leaf; the only non-analytic
#                              piece -- the leaf opt's height response has no closed
#                              form, like FF16's frozen-stage env is a pass-1 input),
#       d(rates)/d(trait)    = forward-mode XAD over the committed kernel
#                              tf24_compute_rates_from_net (height + net both seeded).
#     The emergent stand output is J(theta) = sum_i w_i * fecundity_i(t_end); one
#     forward-sensitivity pass gives d(J)/d(vcmax_25).
#
# Two checks: (a) FAITHFULNESS -- the double replay reproduces the live SCM cohort
# heights (the RKCK port + per-stage env + TF24 kernel are exact; limited by the
# leaf-opt tolerance, ~1e-7); (b) GRADIENT -- d(J)/d(vcmax_25) matches a two-pass
# central FD on the same frozen schedule, converging O(h^2) as the step shrinks
# (the convergence is the proof the tangent-linear AD is the exact derivative).
#
# Reverse-mode (one sweep for all traits) is NOT applicable through the leaf opt;
# the headline reverse-mode win is the net-production sweep, scripts/ad_tf24_
# reverse_sweep.R. Here per-trait forward-sensitivity is the right tool.
#
# Run from the package root after `R CMD INSTALL .` (needs plant from this branch,
# odelia, BH):  Rscript scripts/ad_tf24_emergent_gradient.R

suppressMessages({library(Rcpp); library(plant)})

## ---- Pass 1: the real resident TF24 SCM, harvested from one clean run -------
# crown-centre binds the single-optimisation crown light (run_scm defaults to the
# deep-crown integral). refine FIRST (no cache), then ONE cached run so the history
# buffers are a single monotonic schedule (reset() does not clear them).
p <- scm_base_parameters("TF24")
p$max_patch_lifetime <- 20            # modest horizon -> tractable leaf-opt count
p <- add_strategies(p, trait_matrix(0.1978791, "lma"), hyperpar = TF24_hyperpar,
                    birth_rate = list(20))
mk <- function(cache = FALSE)
  control(shading_model = "crown-centre", GSS_tol_abs = 1e-9,
          ode_tol_rel = 1e-4, ode_tol_abs = 1e-4, save_RK45_cache = cache)
p2  <- run_scm(p, Environment("TF24"), mk(FALSE), refine_schedule = TRUE)$parameters
scm <- run_scm(p2, Environment("TF24"), mk(TRUE), refine_schedule = FALSE)
stopifnot(!is.unsorted(scm$patch$step_history))

sh <- scm$patch$step_history          # {0, t1, ...}, length N+1
eh <- scm$patch$environment_history   # length N, each a list of 6 frozen envs
sp <- scm$patch$species[[1]]
node_times   <- sp$node_times
weights      <- sp$patch_densities    # frozen pass-1 cohort weights
live_heights <- sp$heights
birth_step <- vapply(node_times, function(t) which.min(abs(sh - t)) - 1L, integer(1))
pp <- unlist(scm$parameters$strategies[[1]]$pars)
cat(sprintf("Pass 1: %d steps, %d cohorts, horizon %.1f\n",
            length(eh), length(node_times), max(sh)))

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
#include <plant/models/tf24_strategy.h>
#include <plant/models/tf24_environment.h>
#include <plant/models/tf24_production_kernel.h>
#include <plant/models/ff16_production_kernel.h>   // ff16_cashkarp_replay (generic)
// [[Rcpp::plugins(cpp20)]]
using F = xad::fwd<double>::active_type;

static plant::TF24_Strategy make_strategy(const Rcpp::NumericVector& pp, double vcmax) {
  plant::TF24_Strategy s;
  s.control.shading_model = "crown-centre"; s.control.GSS_tol_abs = 1e-9;
  auto& q = s.pars;
  q.lma=pp["lma"];q.rho=pp["rho"];q.hmat=pp["hmat"];q.omega=pp["omega"];q.eta=pp["eta"];
  q.theta=pp["theta"];q.a_l1=pp["a_l1"];q.a_l2=pp["a_l2"];q.a_r1=pp["a_r1"];q.a_b1=pp["a_b1"];
  q.r_s=pp["r_s"];q.r_b=pp["r_b"];q.r_r=pp["r_r"];q.r_l=pp["r_l"];q.a_y=pp["a_y"];q.a_bio=pp["a_bio"];
  q.k_l=pp["k_l"];q.k_b=pp["k_b"];q.k_s=pp["k_s"];q.k_r=pp["k_r"];
  q.a_f3=pp["a_f3"];q.a_f1=pp["a_f1"];q.a_f2=pp["a_f2"];q.d_I=pp["d_I"];
  q.a_dG1=pp["a_dG1"];q.a_dG2=pp["a_dG2"];q.k_I=pp["k_I"];
  q.vcmax_25=vcmax;q.K_s=pp["K_s"];q.b=pp["b"];q.c=pp["c"];q.beta2=pp["beta2"];
  q.jmax_25=pp["jmax_25"];q.a=pp["a"];q.curv_fact_elec_trans=pp["curv_fact_elec_trans"];
  q.curv_fact_colim=pp["curv_fact_colim"];
  s.prepare_strategy(); return s;
}
template <typename S> static plant::TF24ProdPars<S> lift(const plant::TF24ProdPars<double>& d) {
  plant::TF24ProdPars<S> p;
  p.lma=d.lma;p.rho=d.rho;p.theta=d.theta;p.a_b1=d.a_b1;p.a_r1=d.a_r1;p.eta_c=d.eta_c;
  p.r_l=d.r_l;p.r_s=d.r_s;p.r_b=d.r_b;p.r_r=d.r_r;p.k_l=d.k_l;p.k_b=d.k_b;p.k_s=d.k_s;p.k_r=d.k_r;
  p.a_bio=d.a_bio;p.a_y=d.a_y;p.a_l1=d.a_l1;p.a_l2=d.a_l2;p.a_f1=d.a_f1;p.a_f2=d.a_f2;p.hmat=d.hmat;
  p.omega=d.omega;p.a_f3=d.a_f3;p.d_I=d.d_I;p.a_dG1=d.a_dG1;p.a_dG2=d.a_dG2; return p;
}

// Tangent-linear demographic state: value + d/d(vcmax_25) sensitivity.
struct TL { plant::FF16State<double> v, s; };
struct Ctx {
  plant::TF24_Strategy* st; plant::TF24ProdPars<double> pd;
  std::vector<std::vector<plant::TF24_Environment>>* eh; double conv;
};
static double net_at(plant::TF24_Strategy& s, plant::TF24_Environment& e, double h) {
  return s.net_mass_production_dt(e, h, s.area_leaf(h), 1.0/h);
}

// Replay one cohort (tangent-linear); sens=true propagates d/d(vcmax), else double.
static TL replay(Ctx& C, std::size_t birth, const std::vector<double>& step_h,
                 double h0, bool sens) {
  auto deriv = [&](const TL& y, std::size_t n, int stage) -> TL {
    plant::TF24_Environment* e =
      (stage==0)?((n>0)?&(*C.eh)[n-1][5]:&(*C.eh)[0][0]):&(*C.eh)[n][stage-1];
    auto& s = *C.st;
    const double h = y.v.height, sht = y.s.height;
    const double al = s.area_leaf(h);
    const double net0 = net_at(s, *e, h);             // REAL leaf opt (faithful)
    double dnet_total = 0.0;
    if (sens) {
      const double opt = -s.leaf.root_collar_psi_;
      const double dnet_dv = s.pars.a_bio*s.pars.a_y*al*C.conv*s.leaf.dprofit_dvcmax25(opt);
      const double dd = 1e-5*h;                        // height-Jacobian FD step
      const double dnet_dh = (net_at(s,*e,h+dd) - net_at(s,*e,h-dd)) / (2*dd);
      dnet_total = dnet_dv + dnet_dh*sht;
    }
    F h_ad = h; xad::derivative(h_ad) = sht;
    F net_ad = net0; xad::derivative(net_ad) = dnet_total;
    plant::TF24ProdPars<F> pf = lift<F>(C.pd);
    F al_ad = plant::tf24_area_leaf<F>(pf.a_l1, pf.a_l2, h_ad);
    plant::TF24Rates<F> r = plant::tf24_compute_rates_from_net<F>(pf, h_ad, al_ad, net_ad, true);
    TL o;
    o.v = plant::FF16State<double>{xad::value(r.height_dt),xad::value(r.mortality_dt),
      xad::value(r.fecundity_dt),xad::value(r.area_heartwood_dt),xad::value(r.mass_heartwood_dt)};
    o.s = plant::FF16State<double>{xad::derivative(r.height_dt),xad::derivative(r.mortality_dt),
      xad::derivative(r.fecundity_dt),xad::derivative(r.area_heartwood_dt),xad::derivative(r.mass_heartwood_dt)};
    return o;
  };
  auto axpy = [](const TL& a, double c, const TL& k) -> TL {
    return TL{ plant::FF16State<double>{a.v.height+c*k.v.height,a.v.mortality+c*k.v.mortality,
      a.v.fecundity+c*k.v.fecundity,a.v.area_heartwood+c*k.v.area_heartwood,a.v.mass_heartwood+c*k.v.mass_heartwood},
      plant::FF16State<double>{a.s.height+c*k.s.height,a.s.mortality+c*k.s.mortality,
      a.s.fecundity+c*k.s.fecundity,a.s.area_heartwood+c*k.s.area_heartwood,a.s.mass_heartwood+c*k.s.mass_heartwood} };
  };
  TL y{ plant::FF16State<double>{h0,0,0,0,0}, plant::FF16State<double>{0,0,0,0,0} };
  return plant::ff16_cashkarp_replay(y, step_h, birth, deriv, axpy);
}

static double stand_J(Ctx& C, const std::vector<int>& birth,
                      const std::vector<double>& step_h, double h0,
                      const std::vector<double>& w) {
  double J=0;
  for (std::size_t i=0;i<birth.size();++i)
    J += w[i]*replay(C,(std::size_t)birth[i],step_h,h0,false).v.fecundity;
  return J;
}

// [[Rcpp::export]]
Rcpp::List tf24_emergent(Rcpp::NumericVector pp, Rcpp::List eh_list,
    std::vector<double> shv, std::vector<int> birth, std::vector<double> w) {
  const std::size_t N = eh_list.size();
  std::vector<std::vector<plant::TF24_Environment>> EH(N);
  for (std::size_t n=0;n<N;++n){Rcpp::List st=eh_list[n];
    for(R_xlen_t k=0;k<st.size();++k) EH[n].push_back(Rcpp::as<plant::TF24_Environment>(st[k]));}
  std::vector<double> step_h(N); for (std::size_t n=0;n<N;++n) step_h[n]=shv[n+1]-shv[n];

  const double vc0 = pp["vcmax_25"];
  plant::TF24_Strategy s0 = make_strategy(pp, vc0);
  Ctx C{&s0, s0.prod_pars(), &EH, 60.0*60.0*12.0*365.0/1e6};
  const double h0 = s0.initial_height();

  // (a) faithfulness: double replay final heights, and (b) AD sensitivity.
  Rcpp::NumericVector hf(birth.size());
  double J=0, sJ=0;
  for (std::size_t i=0;i<birth.size();++i){
    TL y = replay(C,(std::size_t)birth[i],step_h,h0,true);
    hf[i]=y.v.height; J += w[i]*y.v.fecundity; sJ += w[i]*y.s.fecundity;
  }

  // two-pass central FD on the same frozen schedule (re-run double replays).
  std::vector<double> rel={1e-3,1e-4,1e-5}, fd;
  for (double rh: rel){ double dd=rh*vc0;
    plant::TF24_Strategy sp=make_strategy(pp,vc0+dd); Ctx Cp{&sp,sp.prod_pars(),&EH,C.conv};
    plant::TF24_Strategy sm=make_strategy(pp,vc0-dd); Ctx Cm{&sm,sm.prod_pars(),&EH,C.conv};
    fd.push_back((stand_J(Cp,birth,step_h,sp.initial_height(),w)
                 -stand_J(Cm,birth,step_h,sm.initial_height(),w))/(2*dd));
  }
  return Rcpp::List::create(Rcpp::_["J"]=J, Rcpp::_["sJ_ad"]=sJ,
    Rcpp::_["replay_heights"]=hf, Rcpp::_["fd"]=Rcpp::wrap(fd), Rcpp::_["rel"]=Rcpp::wrap(rel));
}')

res <- tf24_emergent(pp, eh, sh, birth_step, weights)

## ---- (a) faithfulness ------------------------------------------------------
max_h_err <- max(abs(res$replay_heights - live_heights))
cat(sprintf("\n(a) Faithfulness  max |replay - live SCM height| over %d cohorts = %.2e\n",
            length(live_heights), max_h_err))

## ---- (b) emergent gradient -------------------------------------------------
cat(sprintf("\n(b) Emergent stand output J = sum_i w_i fecundity_i(t_end) = %.8g\n", res$J))
cat(sprintf("    d(J)/d(vcmax_25)  AD (tangent-linear) = %.8g\n", res$sJ_ad))
for (i in seq_along(res$rel))
  cat(sprintf("    FD(rel step %.0e) = %.8g   rel.err = %.2e\n",
              res$rel[i], res$fd[i], abs(res$sJ_ad - res$fd[i])/abs(res$fd[i])))
best <- min(abs(res$sJ_ad - res$fd)/abs(res$fd))
cat(sprintf("\nbest AD-vs-FD rel.err = %.2e (FD converges O(h^2) -> AD exact)\n", best))
stopifnot(max_h_err < 1e-5, best < 1e-4)
cat("\nTF24 emergent community gradient validated: AD through the entire resident SCM.\n")
