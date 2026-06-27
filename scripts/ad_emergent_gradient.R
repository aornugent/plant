# Live-SCM two-pass emergent trait gradient (#472 scope B, Milestone C / #537).
#
# The headline scope-B result: a reverse-mode trait gradient of an EMERGENT,
# community-level FF16 output, computed over the schedule and resident light of a
# REAL Solver-for-Characteristics-Method run -- not a hand-built stand-in.
#
#   Pass 1 (double): run the real FF16 resident SCM to completion with the
#     adaptive Cash-Karp RKCK solver and save_RK45_cache. Harvest the frozen
#     schedule (Patch$step_history -> the actual adaptive step sizes) and the
#     per-RK-stage resident light (Patch$environment_history[step][0..5]), plus
#     each cohort's birth step and weight. (step_history/environment_history are
#     exposed on Patch for exactly this; see inst/RcppR6_classes.yml.)
#   Pass 2 (AD): replay every cohort with ff16_replay_cohort_rkck -- the SAME
#     Cash-Karp stepper the SCM used (NOT forward Euler), reading the FROZEN
#     per-stage resident light actively at the cohort's crown height (value +
#     slope, so the within-cohort self-shading feedback flows). Form an emergent
#     stand output J(theta) = sum_i w_i * fecundity_i(t_end) and take ONE reverse
#     sweep for d(J)/d(trait).
#
# This is the faithful counterpart of the mutant-fitness replay (run_mutant ->
# advance_fixed + cached environment), lifted to an XAD active scalar. Two checks:
#   (a) faithfulness  -- the double replay reproduces the live SCM cohort heights
#                        to machine precision (the RKCK port + per-stage env wiring
#                        are exact);
#   (b) gradient      -- d(J)/d(a_p1) by AD matches a two-pass central finite
#                        difference on the same frozen schedule (h -> 0 limit).
#
# Requirements: `plant` INSTALLED from this branch (installed headers must match
# its compiled .so), plus `odelia` (XAD tape) and `BH`. Run from the package root
# after `R CMD INSTALL .`:
#   Rscript scripts/ad_emergent_gradient.R

suppressMessages({library(Rcpp); library(plant)})

## ---- Pass 1: the real resident SCM, harvested from a single clean run -----
# crown-centre shading binds FF16_Strategy::assimilation_crown_top, matching the
# crown-top kernel exactly (run_scm otherwise defaults to the deep-crown integral).
# refine_schedule() runs the SCM repeatedly and reset() does NOT clear the history
# buffers, so refine FIRST (no cache) then take ONE clean cached run.
p <- scm_base_parameters("FF16")
p <- add_strategies(p, trait_matrix(0.0825, "lma"), hyperpar = FF16_hyperpar,
                    birth_rate = list(20))
ctrl_refine <- control(); ctrl_refine$shading_model <- "crown-centre"
p <- run_scm(p, Environment("FF16"), ctrl_refine, refine_schedule = TRUE)$parameters
ctrl <- control(save_RK45_cache = TRUE); ctrl$shading_model <- "crown-centre"
scm  <- run_scm(p, Environment("FF16"), ctrl, refine_schedule = FALSE)
stopifnot(!is.unsorted(scm$patch$step_history))   # single clean run => monotonic

sh <- scm$patch$step_history          # {0, t1, ...}, length N+1
eh <- scm$patch$environment_history   # length N, each a list of 6 frozen envs
sp <- scm$patch$species[[1]]
node_times   <- sp$node_times
live_heights <- sp$heights
weights      <- sp$patch_densities    # frozen pass-1 cohort weights
pp           <- unlist(scm$parameters$strategies[[1]]$pars)
birth_step   <- vapply(node_times, function(t) which.min(abs(sh - t)) - 1L, integer(1))
cat(sprintf("Pass 1: %d ODE steps, %d cohorts, t_end = %.2f\n",
            length(eh), length(node_times), max(sh)))

## ---- Pass 2: Cash-Karp RKCK replay carried in an XAD active type ----------
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
  q.lma=pp["lma"]; q.rho=pp["rho"]; q.hmat=pp["hmat"]; q.omega=pp["omega"];
  q.eta=pp["eta"]; q.theta=pp["theta"]; q.a_l1=pp["a_l1"]; q.a_l2=pp["a_l2"];
  q.a_r1=pp["a_r1"]; q.a_b1=pp["a_b1"]; q.r_s=pp["r_s"]; q.r_b=pp["r_b"];
  q.r_r=pp["r_r"]; q.r_l=pp["r_l"]; q.a_y=pp["a_y"]; q.a_bio=pp["a_bio"];
  q.k_l=pp["k_l"]; q.k_b=pp["k_b"]; q.k_s=pp["k_s"]; q.k_r=pp["k_r"];
  q.a_p1=pp["a_p1"]; q.a_p2=pp["a_p2"]; q.a_f3=pp["a_f3"]; q.a_f1=pp["a_f1"];
  q.a_f2=pp["a_f2"]; q.S_D=pp["S_D"]; q.a_d0=pp["a_d0"]; q.d_I=pp["d_I"];
  q.a_dG1=pp["a_dG1"]; q.a_dG2=pp["a_dG2"]; q.k_I=pp["k_I"];
  q.recruitment_decay=pp["recruitment_decay"];
  s.prepare_strategy();
  return s;
}
template <typename S>
static plant::FF16ProdPars<S> lift(const plant::FF16ProdPars<double>& d) {
  plant::FF16ProdPars<S> p;
  p.lma=d.lma;p.rho=d.rho;p.theta=d.theta;p.a_b1=d.a_b1;p.a_r1=d.a_r1;p.eta_c=d.eta_c;
  p.a_p1=d.a_p1;p.a_p2=d.a_p2;p.r_l=d.r_l;p.r_s=d.r_s;p.r_b=d.r_b;p.r_r=d.r_r;
  p.k_l=d.k_l;p.k_b=d.k_b;p.k_s=d.k_s;p.k_r=d.k_r;p.a_bio=d.a_bio;p.a_y=d.a_y;
  p.a_l1=d.a_l1;p.a_l2=d.a_l2;p.a_f1=d.a_f1;p.a_f2=d.a_f2;p.hmat=d.hmat;
  p.omega=d.omega;p.a_f3=d.a_f3;p.d_I=d.d_I;p.a_dG1=d.a_dG1;p.a_dG2=d.a_dG2;
  return p;
}

// Materialised frozen resident environment trajectory (built once).
struct Harvest {
  std::vector<std::vector<plant::FF16_Environment>> eh;  // [step][0..5]
  std::vector<double> step_h;                            // adaptive step sizes
  double eta_c, h0;
};
static Harvest build(const plant::FF16_Strategy& s, Rcpp::List eh_list,
                     const std::vector<double>& sh) {
  Harvest H; H.eta_c = s.prod_pars().eta_c; H.h0 = s.initial_height();
  const std::size_t N = eh_list.size();
  H.eh.resize(N);
  for (std::size_t n = 0; n < N; ++n) {
    Rcpp::List st = eh_list[n];
    for (R_xlen_t k = 0; k < st.size(); ++k)
      H.eh[n].push_back(Rcpp::as<plant::FF16_Environment>(st[k]));
  }
  H.step_h.resize(N);
  for (std::size_t n = 0; n < N; ++n) H.step_h[n] = sh[n + 1] - sh[n];
  return H;
}

// stage 0 -> step-start env (prev step final stage; birth env for n==0);
// stage 1..5 -> environment_history[n][0..4] (the k2..k6 derivs envs). The crown
// reads light at height*eta_c; for an active (ad) height we seed value + slope
// from the FROZEN spline so d(light)/d(height) flows (frozen-knot self-shading).
template <typename S>
static S crown_light(const Harvest& H, std::size_t n, int stage, S height) {
  const plant::FF16_Environment* e =
    (stage == 0) ? ((n > 0) ? &H.eh[n - 1][5] : &H.eh[0][0]) : &H.eh[n][stage - 1];
  const double hd = as_double(height), z = hd * H.eta_c;
  const double Lv = e->get_environment_at_height(z);
  const double Ld = e->get_environment_deriv_at_height(z) * H.eta_c;
  return S(Lv) + S(Ld) * (height - S(hd));   // value + slope (slope*0 for double)
}

// Emergent stand output J(theta) = sum_i w_i * fecundity_i(t_end).
template <typename S>
static S stand_J(const plant::FF16ProdPars<S>& pd, const Harvest& H,
                 const std::vector<int>& birth, const std::vector<double>& w) {
  auto cl = [&](std::size_t n, int stage, S h){ return crown_light<S>(H, n, stage, h); };
  S J = S(0.0);
  for (std::size_t i = 0; i < birth.size(); ++i) {
    plant::FF16State<S> y{S(H.h0), S(0), S(0), S(0), S(0)};
    y = plant::ff16_replay_cohort_rkck<S>(pd, y, H.step_h, (std::size_t)birth[i],
                                          cl, true);
    J += S(w[i]) * y.fecundity;
  }
  return J;
}

// [[Rcpp::export]]
Rcpp::List emergent(Rcpp::NumericVector pp, Rcpp::List eh_list,
                    std::vector<double> sh, std::vector<int> birth,
                    std::vector<double> w) {
  auto s  = make_strategy(pp);
  auto pd = s.prod_pars();
  Harvest H = build(s, eh_list, sh);

  // (a) faithfulness: double replay final heights.
  Rcpp::NumericVector hf(birth.size());
  for (std::size_t i = 0; i < birth.size(); ++i) {
    plant::FF16State<double> y{H.h0,0,0,0,0};
    auto cl = [&](std::size_t n, int st, double h){ return crown_light<double>(H,n,st,h); };
    y = plant::ff16_replay_cohort_rkck<double>(pd, y, H.step_h, (std::size_t)birth[i], cl, true);
    hf[i] = y.height;
  }

  const double Jd = stand_J<double>(pd, H, birth, w);

  // (b) reverse-mode AD: d(J)/d(a_p1) in one sweep.
  double dJ_ad;
  {
    ad::tape_type tape;
    ad_t a_p1 = pd.a_p1; tape.registerInput(a_p1); tape.newRecording();
    auto pa = lift<ad_t>(pd); pa.a_p1 = a_p1;
    ad_t J = stand_J<ad_t>(pa, H, birth, w);
    tape.registerOutput(J); xad::derivative(J) = 1.0; tape.computeAdjoints();
    dJ_ad = xad::derivative(a_p1);
  }

  // Two-pass central FD on the same frozen schedule, swept over step sizes.
  auto Jof = [&](double v){ auto q = pd; q.a_p1 = v; return stand_J<double>(q, H, birth, w); };
  std::vector<double> rel_h = {1e-3,1e-4,1e-5,1e-6,1e-7}, fd;
  for (double rh : rel_h) { double h = rh*pd.a_p1; fd.push_back((Jof(pd.a_p1+h)-Jof(pd.a_p1-h))/(2*h)); }

  return Rcpp::List::create(Rcpp::_["J"]=Jd, Rcpp::_["dJ_ad"]=dJ_ad,
                            Rcpp::_["replay_heights"]=hf,
                            Rcpp::_["fd"]=Rcpp::wrap(fd),
                            Rcpp::_["rel_h"]=Rcpp::wrap(rel_h));
}')

res <- emergent(pp, eh, sh, birth_step, weights)

## ---- (a) faithfulness: replay vs live SCM heights ------------------------
max_h_err <- max(abs(res$replay_heights - live_heights))
cat(sprintf("\n(a) Faithfulness  max |replay - live SCM height| over %d cohorts = %.2e\n",
            length(live_heights), max_h_err))
stopifnot(max_h_err < 1e-8)

## ---- (b) emergent gradient: AD vs two-pass FD (h -> 0 limit) --------------
cat(sprintf("\n(b) Emergent stand fecundity  J = sum_i w_i * fecundity_i(t_end) = %.8g\n", res$J))
cat("    two-pass central FD vs AD (AD is the h->0 limit; O(h^2) convergence):\n")
for (k in seq_along(res$rel_h))
  cat(sprintf("      h/a_p1 = %.0e   FD = %.9g   rel.err = %.2e\n",
              res$rel_h[k], res$fd[k], abs(res$fd[k]-res$dJ_ad)/abs(res$dJ_ad)))
best <- min(abs(res$fd - res$dJ_ad) / abs(res$dJ_ad))
cat(sprintf("\n    d(J)/d(a_p1):  AD = %.9g   best FD = %.9g   min rel.err = %.2e %s\n",
            res$dJ_ad, res$fd[which.min(abs(res$fd-res$dJ_ad))], best,
            if (best < 1e-5) "OK" else "** MISMATCH **"))
stopifnot(best < 1e-5)

cat("\nLive-SCM two-pass emergent trait gradient validated",
    "(faithful RKCK replay + reverse-mode dJ/dtrait).\n")
