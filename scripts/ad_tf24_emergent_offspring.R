# TF24 EMERGENT gradient of the REAL SCM offspring_production (#472 scope B,
# Phase F1-full) -- the canonical emergent scalar, vs the stand-fecundity proxy of
# scripts/ad_tf24_emergent_gradient.R. AD through the entire TF24 SCM.
#
# offspring_production = trapezium over node times of
#   offspring_weighted_i * patch_density_i * S_D * birth_rate,
# where offspring_weighted_i is a 6th survival-weighted ODE state (mirrors
# Node::compute_rates):
#   d(offspring)/dt = fecundity_dt * exp(-mortality) * pr_patch_survival(t)/ppsab_i,
# the node's mortality initialised to -log(establishment_probability(birth env)).
#
# Same two-pass tangent-linear machinery as ad_tf24_emergent_gradient.R: pass 1 runs
# the live crown-centre TF24 resident SCM (RKCK + save_RK45_cache) and harvests the
# frozen schedule, per-RK-stage env, cohort birth steps + survival weights; pass 2
# replays each cohort with the shared Cash-Karp stepper over a 6-state tangent-linear
# value+sensitivity vector. The deriv runs the REAL leaf opt for the faithful net,
# injects d(net)/d(vcmax_25) = a_bio*a_y*area_leaf*conv*dprofit_dvcmax25 and the
# FD height-Jacobian, and carries the offspring accumulator's sensitivity
#   d(off_dt)/dvcmax = sw * exp(-mortality) * (s_fecundity_dt - fecundity_dt*s_mortality).
# ESTABLISHMENT is FROZEN here (mortality_0 from the resident vcmax, in AD AND FD --
# a clean separable partial; differentiating it is the #539/C-26 follow-up).
#
# Checks: (a) reconstruction -- the double replay's trapezium offspring_production
# matches the live SCM; (b) gradient -- d(offspring_production)/d(vcmax_25) AD vs a
# two-pass central FD on the same frozen schedule (O(h^2) -> AD exact).
#
# Run from the package root after `R CMD INSTALL .`:
#   Rscript scripts/ad_tf24_emergent_offspring.R

suppressMessages({library(Rcpp); library(plant)})

p <- scm_base_parameters("TF24")
p$max_patch_lifetime <- 20
p <- add_strategies(p, trait_matrix(0.1978791, "lma"), hyperpar = TF24_hyperpar,
                    birth_rate = list(20))
mk <- function(cache = FALSE)
  control(shading_model = "crown-centre", GSS_tol_abs = 1e-9,
          ode_tol_rel = 1e-4, ode_tol_abs = 1e-4, save_RK45_cache = cache)
p2  <- run_scm(p, Environment("TF24"), mk(FALSE), refine_schedule = TRUE)$parameters
scm <- run_scm(p2, Environment("TF24"), mk(TRUE), refine_schedule = FALSE)
stopifnot(!is.unsorted(scm$patch$step_history))

sh <- scm$patch$step_history
eh <- scm$patch$environment_history
sp <- scm$patch$species[[1]]
node_times <- sp$node_times
pdens      <- sp$patch_densities
ppsab      <- sp$pr_patch_survival_at_birth
S_D <- scm$parameters$strategies[[1]]$pars[["S_D"]]; br <- 20
pp  <- unlist(scm$parameters$strategies[[1]]$pars)
N   <- length(eh)
birth_step <- vapply(node_times, function(t) which.min(abs(sh - t)) - 1L, integer(1))

# Trapezoid coefficients over node times -> emergent post-weighting (frozen).
x <- node_times; nn <- length(x); tcoef <- numeric(nn)
tcoef[1] <- 0.5*(x[2]-x[1]); tcoef[nn] <- 0.5*(x[nn]-x[nn-1])
if (nn > 2) tcoef[2:(nn-1)] <- 0.5*(x[3:nn] - x[1:(nn-2)])
tw <- tcoef * pdens * S_D * br

# Per-RK-stage frozen patch survival (Cash-Karp node fractions).
ah <- c(0.0, 0.2, 0.3, 0.6, 1.0, 0.875); hN <- diff(sh)
ppsurv <- matrix(0.0, N, 6)
for (k in seq_len(N)) for (s in 1:6) ppsurv[k, s] <- scm$patch$pr_survival(sh[k] + ah[s]*hN[k])
cat(sprintf("Pass 1: %d steps, %d cohorts, SCM offspring_production = %.8g\n",
            N, nn, scm$offspring_production))

plant_inc<-system.file("include",package="plant");odelia_inc<-system.file("include",package="odelia");bh_inc<-system.file("include",package="BH")
plant_so<-system.file("libs","plant.so",package="plant");odelia_so<-system.file("libs","odelia.so",package="odelia")
if (!all(nzchar(c(plant_inc,odelia_inc,bh_inc))) || !all(file.exists(c(plant_so,odelia_so))))
  stop("Need plant (installed from this branch), odelia and BH; run R CMD INSTALL .")
Sys.setenv(PKG_CPPFLAGS=paste(paste0("-I",shQuote(plant_inc)),paste0("-I",shQuote(odelia_inc)),paste0("-I",shQuote(bh_inc))))
Sys.setenv(PKG_LIBS=paste(shQuote(normalizePath(plant_so)),shQuote(normalizePath(odelia_so))))

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
  s.control.shading_model="crown-centre"; s.control.GSS_tol_abs=1e-9;
  auto& q=s.pars;
  q.lma=pp["lma"];q.rho=pp["rho"];q.hmat=pp["hmat"];q.omega=pp["omega"];q.eta=pp["eta"];
  q.theta=pp["theta"];q.a_l1=pp["a_l1"];q.a_l2=pp["a_l2"];q.a_r1=pp["a_r1"];q.a_b1=pp["a_b1"];
  q.r_s=pp["r_s"];q.r_b=pp["r_b"];q.r_r=pp["r_r"];q.r_l=pp["r_l"];q.a_y=pp["a_y"];q.a_bio=pp["a_bio"];
  q.k_l=pp["k_l"];q.k_b=pp["k_b"];q.k_s=pp["k_s"];q.k_r=pp["k_r"];
  q.a_f3=pp["a_f3"];q.a_f1=pp["a_f1"];q.a_f2=pp["a_f2"];q.d_I=pp["d_I"];
  q.a_dG1=pp["a_dG1"];q.a_dG2=pp["a_dG2"];q.k_I=pp["k_I"];q.a_d0=pp["a_d0"];
  q.recruitment_decay=pp["recruitment_decay"];
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

// 6-state tangent-linear life state: value + d/d(vcmax) for {h,m,f,ah,mh,off}.
struct L6 { double h,m,f,ah,mh,off; };
struct TL { L6 v, s; };
struct Ctx {
  plant::TF24_Strategy* st; plant::TF24ProdPars<double> pd;
  std::vector<std::vector<plant::TF24_Environment>>* eh; double conv;
  const Rcpp::NumericMatrix* ppsurv; double ppsab;
};
static double net_at(plant::TF24_Strategy& s, plant::TF24_Environment& e, double h) {
  return s.net_mass_production_dt(e, h, s.area_leaf(h), 1.0/h);
}

static TL replay(Ctx& C, std::size_t birth, const std::vector<double>& step_h,
                 double h0, double mort0, bool sens) {
  auto deriv = [&](const TL& y, std::size_t n, int stage) -> TL {
    plant::TF24_Environment* e =
      (stage==0)?((n>0)?&(*C.eh)[n-1][5]:&(*C.eh)[0][0]):&(*C.eh)[n][stage-1];
    auto& s = *C.st;
    const double h = y.v.h, sht = y.s.h;
    const double al = s.area_leaf(h);
    const double net0 = net_at(s, *e, h);
    double dnet_total = 0.0;
    if (sens) {
      const double opt = -s.leaf.root_collar_psi_;
      const double dnet_dv = s.pars.a_bio*s.pars.a_y*al*C.conv*s.leaf.dprofit_dvcmax25(opt);
      const double dd = 1e-5*h;
      const double dnet_dh = (net_at(s,*e,h+dd) - net_at(s,*e,h-dd)) / (2*dd);
      dnet_total = dnet_dv + dnet_dh*sht;
    }
    F h_ad=h; xad::derivative(h_ad)=sht;
    F net_ad=net0; xad::derivative(net_ad)=dnet_total;
    plant::TF24ProdPars<F> pf = lift<F>(C.pd);
    F al_ad = plant::tf24_area_leaf<F>(pf.a_l1, pf.a_l2, h_ad);
    plant::TF24Rates<F> r = plant::tf24_compute_rates_from_net<F>(pf, h_ad, al_ad, net_ad, true);
    // 6th state: survival-weighted offspring.
    const double sw = (*C.ppsurv)(n, stage) / C.ppsab;
    using std::exp;
    const double fv = xad::value(r.fecundity_dt), fs = xad::derivative(r.fecundity_dt);
    const double em = exp(-y.v.m);
    const double off_v = fv * em * sw;
    const double off_s = sw * em * (fs - fv * y.s.m);   // d/dvcmax (mort state sens)
    TL o;
    o.v = L6{xad::value(r.height_dt),xad::value(r.mortality_dt),fv,
             xad::value(r.area_heartwood_dt),xad::value(r.mass_heartwood_dt),off_v};
    o.s = L6{xad::derivative(r.height_dt),xad::derivative(r.mortality_dt),fs,
             xad::derivative(r.area_heartwood_dt),xad::derivative(r.mass_heartwood_dt),off_s};
    return o;
  };
  auto axpy = [](const TL& a, double c, const TL& k) -> TL {
    return TL{ L6{a.v.h+c*k.v.h,a.v.m+c*k.v.m,a.v.f+c*k.v.f,a.v.ah+c*k.v.ah,a.v.mh+c*k.v.mh,a.v.off+c*k.v.off},
               L6{a.s.h+c*k.s.h,a.s.m+c*k.s.m,a.s.f+c*k.s.f,a.s.ah+c*k.s.ah,a.s.mh+c*k.s.mh,a.s.off+c*k.s.off} };
  };
  TL y{ L6{h0,mort0,0,0,0,0}, L6{0,0,0,0,0,0} };
  return plant::ff16_cashkarp_replay(y, step_h, birth, deriv, axpy);
}

// mortality at birth = -log(establishment_probability) in the frozen birth env.
static double mort0_for(plant::TF24_Strategy& s, plant::TF24_Environment& eb) {
  return -std::log(s.establishment_probability(eb));
}

static double offspring_production(Ctx& C, const std::vector<int>& birth,
    const std::vector<double>& step_h, double h0, const std::vector<double>& tw,
    std::vector<std::vector<plant::TF24_Environment>>& EH) {
  double J=0;
  for (std::size_t i=0;i<birth.size();++i){
    const std::size_t b=(std::size_t)birth[i];
    plant::TF24_Environment& eb = (b>0)?EH[b-1][5]:EH[0][0];
    double m0 = mort0_for(*C.st, eb);
    J += tw[i]*replay(C,b,step_h,h0,m0,false).v.off;
  }
  return J;
}

// [[Rcpp::export]]
Rcpp::List tf24_emergent_offspring(Rcpp::NumericVector pp, Rcpp::List eh_list,
    std::vector<double> shv, std::vector<int> birth, std::vector<double> tw,
    Rcpp::NumericMatrix ppsurv, std::vector<double> ppsab) {
  const std::size_t N = eh_list.size();
  std::vector<std::vector<plant::TF24_Environment>> EH(N);
  for (std::size_t n=0;n<N;++n){Rcpp::List st=eh_list[n];
    for(R_xlen_t k=0;k<st.size();++k) EH[n].push_back(Rcpp::as<plant::TF24_Environment>(st[k]));}
  std::vector<double> step_h(N); for (std::size_t n=0;n<N;++n) step_h[n]=shv[n+1]-shv[n];

  const double vc0 = pp["vcmax_25"];
  plant::TF24_Strategy s0 = make_strategy(pp, vc0);
  const double h0 = s0.initial_height();

  // AD: reconstruct offspring_production + d/d(vcmax) in one tangent-linear pass.
  // Establishment frozen: mort0_i computed at the RESIDENT vcmax, reused everywhere.
  std::vector<double> mort0(birth.size());
  { Ctx C0{&s0, s0.prod_pars(), &EH, 60.0*60.0*12.0*365.0/1e6, &ppsurv, 1.0};
    for (std::size_t i=0;i<birth.size();++i){
      const std::size_t b=(std::size_t)birth[i];
      plant::TF24_Environment& eb=(b>0)?EH[b-1][5]:EH[0][0];
      mort0[i]=mort0_for(*C0.st, eb);
    }
  }
  double J=0, sJ=0;
  for (std::size_t i=0;i<birth.size();++i){
    Ctx C{&s0, s0.prod_pars(), &EH, 60.0*60.0*12.0*365.0/1e6, &ppsurv, ppsab[i]};
    TL y = replay(C,(std::size_t)birth[i],step_h,h0,mort0[i],true);
    J += tw[i]*y.v.off; sJ += tw[i]*y.s.off;
  }

  // two-pass FD (perturb vcmax; establishment mort0 stays FROZEN at resident value).
  std::vector<double> rel={1e-3,1e-4,1e-5}, fd;
  for (double rh: rel){ double dd=rh*vc0;
    double Jp=0, Jm=0;
    plant::TF24_Strategy sp=make_strategy(pp,vc0+dd);
    plant::TF24_Strategy sm=make_strategy(pp,vc0-dd);
    for (std::size_t i=0;i<birth.size();++i){
      Ctx Cp{&sp,sp.prod_pars(),&EH,60.0*60.0*12.0*365.0/1e6,&ppsurv,ppsab[i]};
      Ctx Cm{&sm,sm.prod_pars(),&EH,60.0*60.0*12.0*365.0/1e6,&ppsurv,ppsab[i]};
      Jp += tw[i]*replay(Cp,(std::size_t)birth[i],step_h,sp.initial_height(),mort0[i],false).v.off;
      Jm += tw[i]*replay(Cm,(std::size_t)birth[i],step_h,sm.initial_height(),mort0[i],false).v.off;
    }
    fd.push_back((Jp-Jm)/(2*dd));
  }
  return Rcpp::List::create(Rcpp::_["offspring"]=J, Rcpp::_["grad_ad"]=sJ,
    Rcpp::_["fd"]=Rcpp::wrap(fd), Rcpp::_["rel"]=Rcpp::wrap(rel));
}')

res <- tf24_emergent_offspring(pp, eh, sh, birth_step, tw, ppsurv, ppsab)

## ---- (a) reconstruction ----------------------------------------------------
rel_recon <- abs(res$offspring - scm$offspring_production) /
             abs(scm$offspring_production)
cat(sprintf("\n(a) Reconstruction  offspring_production replay = %.8g  SCM = %.8g  rel = %.2e\n",
            res$offspring, scm$offspring_production, rel_recon))

## ---- (b) emergent gradient -------------------------------------------------
cat(sprintf("\n(b) d(offspring_production)/d(vcmax_25)  AD = %.8g\n", res$grad_ad))
for (i in seq_along(res$rel))
  cat(sprintf("    FD(rel step %.0e) = %.8g   rel.err = %.2e\n",
              res$rel[i], res$fd[i], abs(res$grad_ad - res$fd[i])/abs(res$fd[i])))
best <- min(abs(res$grad_ad - res$fd)/abs(res$fd))
cat(sprintf("\nbest AD-vs-FD rel.err = %.2e (FD converges O(h^2) -> AD exact)\n", best))
stopifnot(rel_recon < 5e-3, best < 1e-4)
cat("\nTF24 offspring_production gradient validated: AD through the entire resident SCM.\n")
