# Time-integrated self-shading gradient (#472 scope B, Milestone C): a focal cohort
# replayed over its WHOLE lifetime through a resident light that RESPONDS to the trait
# at EVERY step -- the time-integrated active-knot path, extending the single-census
# scripts/ad_self_shading_live.R.
#
# Enabled by the new per-ODE-step stand harvest: Patch$stand_height_history /
# Patch$stand_competition_history record, alongside environment_history during a
# save_RK45_cache run, the species-0 node heights and per-node competition effects
# (ce_i = node.compute_competition(0) = k_I*area_leaf, and Q(0)=1, so the resident
# competition(z) = trapezium_i( ce_i * Q(z/h_i) )). Factoring ce_i = C_i * area_leaf_i
# with C_i = ce_i / area_leaf_i frozen (density * survival weighting) makes the per-step
# resident light differentiable in an allometric trait.
#
# Pass 2 feeds that reconstruction to the committed ff16_replay_cohort_rkck via a
# crown_light callable: at each step the focal crown reads the resident light
# reconstructed from the step-start stand, ACTIVE in a_l1 (every cohort's area_leaf)
# AND in the focal's own crown height (the within-cohort feedback). One reverse sweep
# gives d(focal lifetime fecundity)/d(a_l1) including the self-shading response.
#
# Validation: (a) the per-step reconstruction matches the live FF16_Environment across
# the run; (b) AD vs a two-pass FD over the same reconstruction. The constant-light
# value is reported to isolate the self-shading contribution.
#
# Run from the package root after `R CMD INSTALL .`:
#   Rscript scripts/ad_self_shading_timeint.R
suppressMessages({library(Rcpp); library(plant)})

p <- scm_base_parameters("FF16")
p <- add_strategies(p, trait_matrix(0.0825,"lma"), hyperpar=FF16_hyperpar, birth_rate=list(20))
p <- run_scm(p, Environment("FF16"), control(), refine_schedule=TRUE)$parameters
scm <- run_scm(p, Environment("FF16"), control(save_RK45_cache=TRUE), refine_schedule=FALSE)
stopifnot(!is.unsorted(scm$patch$step_history))

sh <- scm$patch$step_history
eh <- scm$patch$environment_history
shist <- scm$patch$stand_height_history
chist <- scm$patch$stand_competition_history
N <- length(eh)
stopifnot(length(shist) == N, length(chist) == N)
sp <- scm$patch$species[[1]]
pp <- unlist(scm$parameters$strategies[[1]]$pars)
a_l1 <- pp[["a_l1"]]; a_l2 <- pp[["a_l2"]]; eta <- pp[["eta"]]
birth_step <- vapply(sp$node_times, function(t) which.min(abs(sh - t)) - 1L, integer(1))
cat(sprintf("Pass 1: %d steps; stand sizes %d..%d cohorts\n",
            N, length(shist[[1]]), length(shist[[N]])))

## ---- (a) validate per-step reconstruction vs the live (step-end) env ----------
area_leaf <- function(h) (h / a_l1)^(1 / a_l2)
recon_light_R <- function(hv, cv, z) {           # cv = ce_i; C_i*area_leaf == ce_i
  if (length(hv) < 2) return(1.0)
  o <- order(hv, decreasing=TRUE); hh <- hv[o]; cc <- cv[o]
  Q <- ifelse(z/hh < 1, (1 - (z/hh)^eta)^2, 0)
  f <- cc * Q
  comp <- sum(diff(-hh) * (head(f,-1) + tail(f,-1))) / 2
  exp(-comp)
}
chk_steps <- unique(round(seq(2, N, length.out = 8)))
errs <- c()
for (n in chk_steps) {
  envn <- eh[[n]][[6]]                            # step-end env (== stand at step n)
  z <- 5
  errs <- c(errs, abs(recon_light_R(shist[[n]], chist[[n]], z) - envn$get_environment_at_height(z)))
}
cat(sprintf("(a) per-step recon vs live env @z=5: max abs err over %d steps = %.2e\n",
            length(chk_steps), max(errs)))

## ---- (b) time-integrated focal self-shading gradient (C++/AD) -----------------
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
using ad=xad::adj<double>; using ad_t=ad::active_type;
// [[Rcpp::plugins(cpp20)]]
static double as_double(double v){return v;} static double as_double(const ad_t&v){return xad::value(v);}
static plant::FF16_Strategy make_strategy(const Rcpp::NumericVector& pp){
  plant::FF16_Strategy s; auto& q=s.pars;
  q.lma=pp["lma"];q.rho=pp["rho"];q.hmat=pp["hmat"];q.omega=pp["omega"];q.eta=pp["eta"];
  q.theta=pp["theta"];q.a_l1=pp["a_l1"];q.a_l2=pp["a_l2"];q.a_r1=pp["a_r1"];q.a_b1=pp["a_b1"];
  q.r_s=pp["r_s"];q.r_b=pp["r_b"];q.r_r=pp["r_r"];q.r_l=pp["r_l"];q.a_y=pp["a_y"];q.a_bio=pp["a_bio"];
  q.k_l=pp["k_l"];q.k_b=pp["k_b"];q.k_s=pp["k_s"];q.k_r=pp["k_r"];q.a_p1=pp["a_p1"];q.a_p2=pp["a_p2"];
  q.a_f3=pp["a_f3"];q.a_f1=pp["a_f1"];q.a_f2=pp["a_f2"];q.S_D=pp["S_D"];q.a_d0=pp["a_d0"];q.d_I=pp["d_I"];
  q.a_dG1=pp["a_dG1"];q.a_dG2=pp["a_dG2"];q.k_I=pp["k_I"];q.recruitment_decay=pp["recruitment_decay"];
  s.prepare_strategy(); return s;
}
template <typename S> static plant::FF16ProdPars<S> lift(const plant::FF16ProdPars<double>& d){
  plant::FF16ProdPars<S> p;
  p.lma=d.lma;p.rho=d.rho;p.theta=d.theta;p.a_b1=d.a_b1;p.a_r1=d.a_r1;p.eta_c=d.eta_c;
  p.a_p1=d.a_p1;p.a_p2=d.a_p2;p.r_l=d.r_l;p.r_s=d.r_s;p.r_b=d.r_b;p.r_r=d.r_r;
  p.k_l=d.k_l;p.k_b=d.k_b;p.k_s=d.k_s;p.k_r=d.k_r;p.a_bio=d.a_bio;p.a_y=d.a_y;
  p.a_l1=d.a_l1;p.a_l2=d.a_l2;p.a_f1=d.a_f1;p.a_f2=d.a_f2;p.hmat=d.hmat;
  p.omega=d.omega;p.a_f3=d.a_f3;p.d_I=d.d_I;p.a_dG1=d.a_dG1;p.a_dG2=d.a_dG2; return p;
}

// per-step resident stand: heights (desc) + frozen weights C_i = ce_i/area_leaf_i(base).
struct Stand { std::vector<std::vector<double>> h, C; double eta; };

// Reconstruct resident light at z from the stand at step n. ACTIVE in a_l1 via each
// resident area_leaf (the self-shading reshaping) AND in z (the focal crown height,
// so d(light)/d(focal height) -- the within-cohort feedback -- also flows). Resident
// heights hv[i] and the weights Cv[i] stay frozen doubles.
template <typename S>
static S recon_light(const plant::FF16ProdPars<S>& p, const Stand& st, std::size_t n, S z) {
  using std::pow; using std::exp;
  const auto& hv = st.h[n]; const auto& Cv = st.C[n];
  if (hv.size() < 2) return S(1.0);
  auto g = [&](std::size_t i) -> S {
    if (as_double(z) >= hv[i]) return S(0.0);
    S u = z / S(hv[i]); S om = S(1.0) - pow(u, st.eta);
    return S(Cv[i]) * plant::ff16_area_leaf(p.a_l1, p.a_l2, S(hv[i])) * (om * om);
  };
  S comp = S(0.0); S gp = g(0); double hp = hv[0];
  for (std::size_t i=1;i<hv.size();++i){ S gi=g(i); comp = comp + S(hp-hv[i])*(gp+gi); hp=hv[i]; gp=gi; }
  return exp(-S(0.5)*comp);
}

// [[Rcpp::export]]
Rcpp::List timeint(Rcpp::NumericVector pp, std::vector<double> sh,
    Rcpp::List shist, Rcpp::List chist, int focal_birth) {
  auto s = make_strategy(pp); auto pd = s.prod_pars();
  const double eta_c = pd.eta_c, h0 = s.initial_height(), eta = s.pars.eta;
  const std::size_t N = sh.size()-1;
  Stand st; st.eta=eta; st.h.resize(N); st.C.resize(N);
  for (std::size_t n=0;n<N;++n){
    std::vector<double> hv = Rcpp::as<std::vector<double>>(shist[n]);
    std::vector<double> cv = Rcpp::as<std::vector<double>>(chist[n]);
    st.h[n]=hv; st.C[n].resize(hv.size());
    for (std::size_t i=0;i<hv.size();++i){
      double al = std::pow(hv[i]/pd.a_l1, 1.0/pd.a_l2);   // base area_leaf
      st.C[n][i] = (al>0)? cv[i]/al : 0.0;                // frozen weight
    }
  }
  std::vector<double> step_h(N); for(std::size_t n=0;n<N;++n) step_h[n]=sh[n+1]-sh[n];

  // crown_light: resident light reconstructed from the step-START stand (held across
  // the RK stages of a step), read at the focal crown height*eta_c, ACTIVE in a_l1.
  auto replay = [&](const auto& p, bool active){
    using S = std::decay_t<decltype(p.a_p1)>;
    auto cl = [&](std::size_t n, int /*stage*/, S height) -> S {
      std::size_t sn = (n>0)? n-1 : 0;            // step-start stand
      S z = height * S(eta_c);
      S L = recon_light<S>(p, st, sn, z);         // active in a_l1 AND in focal height z
      if (!active) L = S(as_double(L));
      return L;
    };
    plant::FF16State<S> y{S(h0),S(0),S(0),S(0),S(0)};
    return plant::ff16_replay_cohort_rkck<S>(p, y, step_h, (std::size_t)focal_birth, cl, true);
  };

  double Jd = replay(pd, true).fecundity;
  double dJ_active, dJ_frozen;
  { ad::tape_type t; ad_t a=pd.a_l1; t.registerInput(a); t.newRecording();
    auto pa=lift<ad_t>(pd); pa.a_l1=a; ad_t J=replay(pa, true).fecundity;
    t.registerOutput(J); xad::derivative(J)=1.0; t.computeAdjoints(); dJ_active=xad::derivative(a); }
  { ad::tape_type t; ad_t a=pd.a_l1; t.registerInput(a); t.newRecording();
    auto pa=lift<ad_t>(pd); pa.a_l1=a; ad_t J=replay(pa, false).fecundity;
    t.registerOutput(J); xad::derivative(J)=1.0; t.computeAdjoints(); dJ_frozen=xad::derivative(a); }
  std::vector<double> rel_h={1e-4,1e-5,1e-6,1e-7}, fd;
  for(double rh:rel_h){ double hh=rh*pd.a_l1; auto q1=pd;q1.a_l1+=hh; auto q2=pd;q2.a_l1-=hh;
    fd.push_back((replay(q1,true).fecundity - replay(q2,true).fecundity)/(2*hh)); }
  return Rcpp::List::create(Rcpp::_["J"]=Jd, Rcpp::_["dJ_active"]=dJ_active,
    Rcpp::_["dJ_frozen"]=dJ_frozen, Rcpp::_["fd"]=Rcpp::wrap(fd), Rcpp::_["rel_h"]=Rcpp::wrap(rel_h));
}')

res <- timeint(pp, sh, shist, chist, birth_step[1])
cat(sprintf("\n(b) focal (born step %d) lifetime fecundity = %.8g\n", birth_step[1], res$J))
for (k in seq_along(res$rel_h))
  cat(sprintf("    h/a_l1=%.0e  FD=%.9g  rel.err=%.2e\n", res$rel_h[k], res$fd[k],
              abs(res$fd[k]-res$dJ_active)/abs(res$dJ_active)))
best <- min(abs(res$fd-res$dJ_active)/max(abs(res$dJ_active),1e-30))
cat(sprintf("\n    d(fecundity)/d(a_l1) active=%.9g  frozen=%.9g  self-shading=%.4g (%.1f%%)  best rel.err=%.2e %s\n",
    res$dJ_active, res$dJ_frozen, res$dJ_active-res$dJ_frozen,
    100*(res$dJ_active-res$dJ_frozen)/res$dJ_active, best, if (best<1e-5) "OK" else "**MISMATCH**"))
stopifnot(max(errs) < 5e-3, best < 1e-5)
cat("\nTime-integrated self-shading gradient validated.\n")
