# Deep-crown gradient of the REAL SCM offspring_production (#472 scope B, Milestone C).
#
# The capstone: differentiate the SCM's actual emergent fitness output,
# offspring_production, for the DEFAULT FF16 shading model (deep-crown crown integral).
# Combines the two earlier pieces:
#   - the survival-weighted-offspring 6th state + node-spacing trapezium of
#     scripts/ad_offspring_gradient.R (crown-top), and
#   - the moving-node Gauss-Kronrod crown integral of scripts/ad_deep_crown_gradient.R.
#
# Pass 2 replays each cohort as a 6-state FF16LifeState (5 FF16 states + survival-
# weighted offspring) through the shared ff16_cashkarp_replay stepper, with a deriv
# that (i) forms `net` via the deep-crown GK integral over the FROZEN per-RK-stage
# resident light, (ii) fills the demographic rates via the shared
# ff16_compute_rates_from_net, and (iii) accumulates
#   d(offspring)/dt = fecundity_dt * exp(-mortality) * pr_patch_survival(t)/ppsab.
# offspring_production is the node-spacing trapezium of offspring * patch_density *
# S_D * birth_rate -- a frozen linear post-weighting, so one reverse sweep gives
# d(offspring_production)/d(trait).
#
# Validation: (a) the double reconstruction matches the SCM scalar; (b) the AD
# gradient matches a two-pass central FD (establishment frozen in both).
#
# Run from the package root after `R CMD INSTALL .`:
#   Rscript scripts/ad_deep_crown_offspring_gradient.R

suppressMessages({library(Rcpp); library(plant)})

## Pass 1: real resident SCM, DEFAULT deep-crown shading, single clean cached run.
p <- scm_base_parameters("FF16")
p <- add_strategies(p, trait_matrix(0.0825, "lma"), hyperpar = FF16_hyperpar,
                    birth_rate = list(20))
p <- run_scm(p, Environment("FF16"), control(), refine_schedule = TRUE)$parameters
scm <- run_scm(p, Environment("FF16"), control(save_RK45_cache = TRUE),
               refine_schedule = FALSE)
stopifnot(!is.unsorted(scm$patch$step_history))

sh <- scm$patch$step_history
eh <- scm$patch$environment_history
sp <- scm$patch$species[[1]]
node_times <- sp$node_times
pdens <- sp$patch_densities
ppsab <- sp$pr_patch_survival_at_birth
S_D   <- p$strategies[[1]]$pars$S_D
br    <- 20
pp    <- unlist(scm$parameters$strategies[[1]]$pars)
birth_step <- vapply(node_times, function(t) which.min(abs(sh - t)) - 1L, integer(1))
N <- length(eh)

tcoef <- numeric(length(node_times)); x <- node_times; nn <- length(x)
tcoef[1] <- 0.5 * (x[2] - x[1]); tcoef[nn] <- 0.5 * (x[nn] - x[nn - 1])
if (nn > 2) tcoef[2:(nn - 1)] <- 0.5 * (x[3:nn] - x[1:(nn - 2)])
tw <- tcoef * pdens * S_D * br

ah <- c(0.0, 0.2, 0.3, 0.6, 1.0, 0.875)
hN <- diff(sh)
ppsurv <- matrix(0.0, N, 6)
for (k in seq_len(N)) for (s in 1:6) ppsurv[k, s] <- scm$patch$pr_survival(sh[k] + ah[s] * hN[k])

cat(sprintf("Pass 1 (deep-crown): %d steps, %d cohorts, SCM offspring_production = %.8g\n",
            N, length(node_times), scm$offspring_production))

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

struct Frozen {
  std::vector<std::vector<plant::FF16_Environment>> eh;
  std::vector<double> step_h; double eta, h0;
  std::vector<int> birth; std::vector<double> mort0, ppsab, tw;
  Rcpp::NumericMatrix ppsurv;
  const plant::quadrature::QK* integ;
};

// J(theta) = sum_i tw_i * offspring_weighted_i, deep-crown 6-state replay.
template <typename S>
static S stand_offspring_deep(const plant::FF16ProdPars<S>& pd, const Frozen& F) {
  S J = S(0.0);
  for (std::size_t i = 0; i < F.birth.size(); ++i) {
    const double ppsab = F.ppsab[i];
    auto deriv = [&](const plant::FF16LifeState<S>& s, std::size_t n, int stage)
        -> plant::FF16LifeState<S> {
      const plant::FF16_Environment* e =
        (stage==0)?((n>0)?&F.eh[n-1][5]:&F.eh[0][0]):&F.eh[n][stage-1];
      const double canopy_top = e->max_environment_height();
      const S height = s.demog.height;
      auto integrand = [&](S z) -> S {
        double zv = as_double(z);
        double lv = e->get_environment_at_height(zv, canopy_top);
        double ld = e->get_environment_deriv_at_height(zv);
        S light = S(lv) + S(std::isfinite(ld)?ld:0.0) * (z - S(zv));
        return plant::ff16_assimilation_leaf<S>(pd.a_p1, pd.a_p2, light) *
               plant::ff16_canopy_q<S>(F.eta, z / height, z);
      };
      S area_leaf = plant::ff16_area_leaf(pd.a_l1, pd.a_l2, height);
      S assim = area_leaf * F.integ->integrate_ad<S>(integrand, S(0.0), height);
      S net = plant::ff16_net_from_components(pd, height, area_leaf, assim);
      plant::FF16Rates<S> r = plant::ff16_compute_rates_from_net(pd, height, area_leaf, net, true);
      using std::exp;
      S off_dt = r.fecundity_dt * exp(-s.demog.mortality) * S(F.ppsurv(n, stage) / ppsab);
      return plant::FF16LifeState<S>{plant::FF16State<S>{r.height_dt, r.mortality_dt,
        r.fecundity_dt, r.area_heartwood_dt, r.mass_heartwood_dt}, off_dt};
    };
    auto axpy = [](const plant::FF16LifeState<S>& a, double c, const plant::FF16LifeState<S>& k)
        -> plant::FF16LifeState<S> {
      return plant::FF16LifeState<S>{plant::FF16State<S>{
        a.demog.height+c*k.demog.height, a.demog.mortality+c*k.demog.mortality,
        a.demog.fecundity+c*k.demog.fecundity, a.demog.area_heartwood+c*k.demog.area_heartwood,
        a.demog.mass_heartwood+c*k.demog.mass_heartwood}, a.offspring+c*k.offspring};
    };
    plant::FF16LifeState<S> y{plant::FF16State<S>{S(F.h0), S(F.mort0[i]), S(0), S(0), S(0)}, S(0)};
    y = plant::ff16_cashkarp_replay(y, F.step_h, (std::size_t)F.birth[i], deriv, axpy);
    J += S(F.tw[i]) * y.offspring;
  }
  return J;
}

// [[Rcpp::export]]
Rcpp::List deep_offspring(Rcpp::NumericVector pp, Rcpp::List eh_list,
    std::vector<double> sh, std::vector<int> birth, Rcpp::NumericMatrix ppsurv,
    std::vector<double> ppsab, std::vector<double> tw) {
  auto s = make_strategy(pp); auto pd = s.prod_pars();
  Frozen F; F.eta=s.pars.eta; F.h0=s.initial_height(); F.birth=birth; F.ppsab=ppsab;
  F.tw=tw; F.ppsurv=ppsurv; F.integ=&s.function_integrator;
  const std::size_t N = eh_list.size(); F.eh.resize(N); F.step_h.resize(N);
  for (std::size_t n=0;n<N;++n){Rcpp::List st=eh_list[n]; for(R_xlen_t k=0;k<st.size();++k) F.eh[n].push_back(Rcpp::as<plant::FF16_Environment>(st[k]));}
  for (std::size_t n=0;n<N;++n) F.step_h[n]=sh[n+1]-sh[n];

  auto sp = plant::make_strategy_ptr(s);
  F.mort0.resize(birth.size());
  for (std::size_t i=0;i<birth.size();++i){
    const std::size_t b=(std::size_t)birth[i];
    const plant::FF16_Environment& eb=(b>0)?F.eh[b-1][5]:F.eh[0][0];
    plant::Individual<plant::FF16_Strategy,plant::FF16_Environment> ind(sp);
    ind.set_state("height", F.h0);
    F.mort0[i] = -std::log(ind.establishment_probability(eb));
  }

  const double Jd = stand_offspring_deep<double>(pd, F);
  double dJ_ad;
  { ad::tape_type tape; ad_t a_p1=pd.a_p1; tape.registerInput(a_p1); tape.newRecording();
    auto pa=lift<ad_t>(pd); pa.a_p1=a_p1;
    ad_t J=stand_offspring_deep<ad_t>(pa, F);
    tape.registerOutput(J); xad::derivative(J)=1.0; tape.computeAdjoints(); dJ_ad=xad::derivative(a_p1); }
  std::vector<double> rel_h={1e-4,1e-5,1e-6,1e-7}, fd;
  for (double rh:rel_h){ double h=rh*pd.a_p1; auto q1=pd; q1.a_p1+=h; auto q2=pd; q2.a_p1-=h;
    fd.push_back((stand_offspring_deep<double>(q1,F)-stand_offspring_deep<double>(q2,F))/(2*h)); }
  return Rcpp::List::create(Rcpp::_["J"]=Jd, Rcpp::_["dJ_ad"]=dJ_ad,
                            Rcpp::_["fd"]=Rcpp::wrap(fd), Rcpp::_["rel_h"]=Rcpp::wrap(rel_h));
}')

res <- deep_offspring(pp, eh, sh, birth_step, ppsurv, ppsab, tw)
re_J <- abs(res$J - scm$offspring_production) / scm$offspring_production
cat(sprintf("\n(a) Reconstructed deep-crown offspring_production = %.8g  (SCM = %.8g, rel.err = %.2e)\n",
            res$J, scm$offspring_production, re_J))
cat("\n(b) d(offspring_production)/d(a_p1):  AD vs two-pass FD (AD = h->0 limit):\n")
for (k in seq_along(res$rel_h))
  cat(sprintf("      h/a_p1 = %.0e   FD = %.9g   rel.err = %.2e\n",
              res$rel_h[k], res$fd[k], abs(res$fd[k]-res$dJ_ad)/abs(res$dJ_ad)))
best <- min(abs(res$fd - res$dJ_ad) / abs(res$dJ_ad))
cat(sprintf("\n      AD = %.9g   best FD = %.9g   min rel.err = %.2e  %s\n",
            res$dJ_ad, res$fd[which.min(abs(res$fd-res$dJ_ad))], best,
            if (best < 1e-5) "OK" else "** MISMATCH **"))
stopifnot(re_J < 1e-4, best < 1e-5)
cat("\nDeep-crown gradient of the real SCM offspring_production validated.\n")
