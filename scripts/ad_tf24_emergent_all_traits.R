# TF24 emergent community gradient over the LIVE resident SCM for ALL 27 traits
# (#472 scope B, Phase F1-full) -- the full d(J)/d(theta_k) vector through the SCM.
#
# Builds on ad_tf24_emergent_gradient.R (which did vcmax_25 + lma). The expensive
# per-RK-stage leaf optimisation harvest -- the optimised profit, the 10 leaf
# sensitivities Leaf::dprofit_d*, k_max, E_up_, and the height-Jacobian d(profit)/
# d(height) -- is TRAIT-INDEPENDENT. So each cohort is integrated ONCE in double
# (running the real leaf opt, RECORDING that harvest per stage), then all 27 trait
# sensitivities are propagated through the SAME recorded harvest by cheap
# forward-mode XAD over the committed kernel (tf24_net_mass_production ->
# tf24_compute_rates_from_net) -- NO further leaf opts. Cost ~ one tangent-linear
# pass for the whole 27-vector.
#
# Each trait's injection (the only per-trait difference):
#   - 10 leaf traits (vcmax_25,g1_TF24,beta2,K_s,b,c,jmax_25,a,curv_elec,curv_colim):
#     d(profit)/d(trait) from the recorded dprofit_d* (K_s via dprofit_dkmax*kmax/K_s);
#   - 13 pure cascade (lma,rho,a_b1,r_l,r_b,r_s,r_r,k_l,k_b,k_s,k_r,a_bio,a_y) and the
#     2 area traits (a_l1,a_l2): seeded in TF24ProdPars, kernel differentiates the
#     cascade/area analytically; lma,rho,a_b1,a_l1,a_l2 also shift the seedling height_0
#     (d(h0)/d(trait) by IFT = FD of initial_height);
#   - 2 leaf-coupled cascade (theta via k_max, a_r1 via E_up_): BOTH a cascade seed AND
#     a profit injection.
#   The within-trajectory height feedback rides the active height + recorded
#   d(profit)/d(height) for every trait.
#
# Validated: faithfulness (double replay heights == live SCM) + d(J)/d(trait) vs a
# two-pass central FD for representatives across the three classes (a leaf, a pure
# cascade, an area and a leaf-coupled trait). Full per-trait FD would be 27 stand
# re-runs; the representatives plus the machinery shared with the single-trait script
# (validated there) cover all classes.
#
# Run from the package root after `R CMD INSTALL .`:
#   Rscript scripts/ad_tf24_emergent_all_traits.R

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
node_times <- sp$node_times; weights <- sp$patch_densities
live_heights <- sp$heights
birth_step <- vapply(node_times, function(t) which.min(abs(sh - t)) - 1L, integer(1))
pp <- unlist(scm$parameters$strategies[[1]]$pars)
cat(sprintf("Pass 1: %d steps, %d cohorts, horizon %.1f\n",
            length(eh), length(node_times), max(sh)))

plant_inc<-system.file("include",package="plant");odelia_inc<-system.file("include",package="odelia");bh_inc<-system.file("include",package="BH")
plant_so<-system.file("libs","plant.so",package="plant");odelia_so<-system.file("libs","odelia.so",package="odelia")
if (!all(nzchar(c(plant_inc,odelia_inc,bh_inc))) || !all(file.exists(c(plant_so,odelia_so))))
  stop("Need plant (installed from this branch), odelia and BH; run R CMD INSTALL .")
Sys.setenv(PKG_CPPFLAGS=paste(paste0("-I",shQuote(plant_inc)),paste0("-I",shQuote(odelia_inc)),paste0("-I",shQuote(bh_inc))))
Sys.setenv(PKG_LIBS=paste(shQuote(normalizePath(plant_so)),shQuote(normalizePath(odelia_so))))

Rcpp::sourceCpp(code = '
#include <Rcpp.h>
#include <vector>
#include <string>
#include <cmath>
#include <XAD/XAD.hpp>
#include <plant.h>
#include <plant/models/tf24_strategy.h>
#include <plant/models/tf24_environment.h>
#include <plant/models/tf24_production_kernel.h>
#include <plant/models/ff16_production_kernel.h>
// [[Rcpp::plugins(cpp20)]]
using F = xad::fwd<double>::active_type;

static const std::vector<std::string> TRAITS = {
  "vcmax_25","g1_TF24","beta2","K_s","b","c","jmax_25","a","curv_elec","curv_colim",
  "lma","rho","a_b1","r_l","r_b","r_s","r_r","k_l","k_b","k_s","k_r","a_bio","a_y",
  "a_l1","a_l2","theta","a_r1"};

static plant::TF24_Strategy make_strategy(const Rcpp::NumericVector& pp,
                                          const std::string& over="", double v=0) {
  plant::TF24_Strategy s;
  s.control.shading_model="crown-centre"; s.control.GSS_tol_abs=1e-9;
  auto& q=s.pars;
  q.lma=pp["lma"];q.rho=pp["rho"];q.hmat=pp["hmat"];q.omega=pp["omega"];q.eta=pp["eta"];
  q.theta=pp["theta"];q.a_l1=pp["a_l1"];q.a_l2=pp["a_l2"];q.a_r1=pp["a_r1"];q.a_b1=pp["a_b1"];
  q.r_s=pp["r_s"];q.r_b=pp["r_b"];q.r_r=pp["r_r"];q.r_l=pp["r_l"];q.a_y=pp["a_y"];q.a_bio=pp["a_bio"];
  q.k_l=pp["k_l"];q.k_b=pp["k_b"];q.k_s=pp["k_s"];q.k_r=pp["k_r"];
  q.a_f3=pp["a_f3"];q.a_f1=pp["a_f1"];q.a_f2=pp["a_f2"];q.d_I=pp["d_I"];
  q.a_dG1=pp["a_dG1"];q.a_dG2=pp["a_dG2"];q.k_I=pp["k_I"];
  q.vcmax_25=pp["vcmax_25"];q.K_s=pp["K_s"];q.b=pp["b"];q.c=pp["c"];q.beta2=pp["beta2"];
  q.jmax_25=pp["jmax_25"];q.a=pp["a"];q.curv_fact_elec_trans=pp["curv_fact_elec_trans"];
  q.curv_fact_colim=pp["curv_fact_colim"];
  if (over=="g1_TF24") s.g1_TF24=v;
  else if (over=="curv_elec") q.curv_fact_elec_trans=v;
  else if (over=="curv_colim") q.curv_fact_colim=v;
  else if (over=="vcmax_25") q.vcmax_25=v; else if (over=="beta2") q.beta2=v;
  else if (over=="K_s") q.K_s=v; else if (over=="b") q.b=v; else if (over=="c") q.c=v;
  else if (over=="jmax_25") q.jmax_25=v; else if (over=="a") q.a=v;
  else if (over=="lma") q.lma=v; else if (over=="rho") q.rho=v; else if (over=="a_b1") q.a_b1=v;
  else if (over=="r_l") q.r_l=v; else if (over=="r_b") q.r_b=v; else if (over=="r_s") q.r_s=v;
  else if (over=="r_r") q.r_r=v; else if (over=="k_l") q.k_l=v; else if (over=="k_b") q.k_b=v;
  else if (over=="k_s") q.k_s=v; else if (over=="k_r") q.k_r=v; else if (over=="a_bio") q.a_bio=v;
  else if (over=="a_y") q.a_y=v; else if (over=="a_l1") q.a_l1=v; else if (over=="a_l2") q.a_l2=v;
  else if (over=="theta") q.theta=v; else if (over=="a_r1") q.a_r1=v;
  else if (!over.empty()) Rcpp::stop("unknown trait "+over);
  s.prepare_strategy(); return s;
}
static double trait_value(const plant::TF24_Strategy& s, const std::string& t) {
  if (t=="g1_TF24") return s.g1_TF24;
  if (t=="curv_elec") return s.pars.curv_fact_elec_trans;
  if (t=="curv_colim") return s.pars.curv_fact_colim;
  const auto& p=s.pars;
  if(t=="vcmax_25")return p.vcmax_25;if(t=="beta2")return p.beta2;if(t=="K_s")return p.K_s;
  if(t=="b")return p.b;if(t=="c")return p.c;if(t=="jmax_25")return p.jmax_25;if(t=="a")return p.a;
  if(t=="lma")return p.lma;if(t=="rho")return p.rho;if(t=="a_b1")return p.a_b1;if(t=="r_l")return p.r_l;
  if(t=="r_b")return p.r_b;if(t=="r_s")return p.r_s;if(t=="r_r")return p.r_r;if(t=="k_l")return p.k_l;
  if(t=="k_b")return p.k_b;if(t=="k_s")return p.k_s;if(t=="k_r")return p.k_r;if(t=="a_bio")return p.a_bio;
  if(t=="a_y")return p.a_y;if(t=="a_l1")return p.a_l1;if(t=="a_l2")return p.a_l2;if(t=="theta")return p.theta;
  if(t=="a_r1")return p.a_r1; Rcpp::stop("?"); return 0;
}
template <typename S> static plant::TF24ProdPars<S> lift(const plant::TF24ProdPars<double>& d) {
  plant::TF24ProdPars<S> p;
  p.lma=d.lma;p.rho=d.rho;p.theta=d.theta;p.a_b1=d.a_b1;p.a_r1=d.a_r1;p.eta_c=d.eta_c;
  p.r_l=d.r_l;p.r_s=d.r_s;p.r_b=d.r_b;p.r_r=d.r_r;p.k_l=d.k_l;p.k_b=d.k_b;p.k_s=d.k_s;p.k_r=d.k_r;
  p.a_bio=d.a_bio;p.a_y=d.a_y;p.a_l1=d.a_l1;p.a_l2=d.a_l2;p.a_f1=d.a_f1;p.a_f2=d.a_f2;p.hmat=d.hmat;
  p.omega=d.omega;p.a_f3=d.a_f3;p.d_I=d.d_I;p.a_dG1=d.a_dG1;p.a_dG2=d.a_dG2; return p;
}
// Seed the (cascade/area/leaf-coupled) ProdPars field for a trait; returns false
// for a pure leaf trait (no ProdPars field).
static bool seed_field(plant::TF24ProdPars<F>& p, const std::string& t) {
  F* f=nullptr;
  if(t=="lma")f=&p.lma; else if(t=="rho")f=&p.rho; else if(t=="a_b1")f=&p.a_b1;
  else if(t=="r_l")f=&p.r_l; else if(t=="r_b")f=&p.r_b; else if(t=="r_s")f=&p.r_s; else if(t=="r_r")f=&p.r_r;
  else if(t=="k_l")f=&p.k_l; else if(t=="k_b")f=&p.k_b; else if(t=="k_s")f=&p.k_s; else if(t=="k_r")f=&p.k_r;
  else if(t=="a_bio")f=&p.a_bio; else if(t=="a_y")f=&p.a_y; else if(t=="a_l1")f=&p.a_l1;
  else if(t=="a_l2")f=&p.a_l2; else if(t=="theta")f=&p.theta; else if(t=="a_r1")f=&p.a_r1;
  if(!f) return false; xad::derivative(*f)=1.0; return true;
}

// Per-RK-stage leaf-opt harvest (trait-independent).
struct H {
  double profit, dprofit_dh, kmax, Eup;
  double dvcmax,dg1,dbeta2,dkmax,db,dc,djmax,da,dcelec,dccolim,dEup;
};
struct TL { plant::FF16State<double> v, s; };
static double profit_at(plant::TF24_Strategy& s, plant::TF24_Environment& e, double h) {
  s.net_mass_production_dt(e, h, s.area_leaf(h), 1.0/h); return s.leaf.profit_;
}
// d(profit)/d(trait) direct term from the recorded harvest + the strategy pars.
static double direct(const std::string& t, const H& h, const plant::TF24_Pars& p) {
  if(t=="vcmax_25")return h.dvcmax; if(t=="g1_TF24")return h.dg1; if(t=="beta2")return h.dbeta2;
  if(t=="b")return h.db; if(t=="c")return h.dc; if(t=="jmax_25")return h.djmax; if(t=="a")return h.da;
  if(t=="curv_elec")return h.dcelec; if(t=="curv_colim")return h.dccolim;
  if(t=="K_s")return h.dkmax*(h.kmax/p.K_s);
  if(t=="theta")return h.dkmax*(h.kmax/p.theta);
  if(t=="a_r1")return h.dEup*(h.Eup/p.a_r1);
  return 0.0;
}

// Record pass: integrate one cohort in double, recording the harvest per stage.
static plant::FF16State<double> record_cohort(plant::TF24_Strategy& s,
    const plant::TF24ProdPars<double>& pd,
    std::vector<std::vector<plant::TF24_Environment>>& EH,
    std::size_t birth, const std::vector<double>& step_h, double h0,
    std::vector<H>& rec) {
  rec.clear();
  auto deriv = [&](const plant::FF16State<double>& y, std::size_t n, int stage)
      -> plant::FF16State<double> {
    plant::TF24_Environment* e =
      (stage==0)?((n>0)?&EH[n-1][5]:&EH[0][0]):&EH[n][stage-1];
    const double h = y.height;
    const double profit_v = profit_at(s, *e, h);
    const double opt = -s.leaf.root_collar_psi_;
    H hh;
    hh.profit = profit_v; hh.kmax = s.leaf.leaf_specific_conductance_max_; hh.Eup = s.leaf.E_up_;
    hh.dvcmax=s.leaf.dprofit_dvcmax25(opt); hh.dg1=s.leaf.dprofit_dg1_TF24(opt);
    hh.dbeta2=s.leaf.dprofit_dbeta2(opt); hh.dkmax=s.leaf.dprofit_dkmax(opt);
    hh.db=s.leaf.dprofit_db(opt); hh.dc=s.leaf.dprofit_dc(opt);
    hh.djmax=s.leaf.dprofit_djmax25(opt); hh.da=s.leaf.dprofit_da(opt);
    hh.dcelec=s.leaf.dprofit_dcurv_elec(opt); hh.dccolim=s.leaf.dprofit_dcurv_colim(opt);
    hh.dEup=s.leaf.dprofit_dEup(opt);
    const double dd=1e-5*h;
    hh.dprofit_dh = (profit_at(s,*e,h+dd)-profit_at(s,*e,h-dd))/(2*dd);
    rec.push_back(hh);
    // value rates via the kernel from the harvested profit (== live net).
    plant::TF24ProdPars<F> pf=lift<F>(pd); F h_ad=h; F prof=profit_v;
    F al=plant::tf24_area_leaf<F>(pf.a_l1,pf.a_l2,h_ad);
    F net=plant::tf24_net_mass_production<F>(pf,h_ad,al,prof);
    plant::TF24Rates<F> r=plant::tf24_compute_rates_from_net<F>(pf,h_ad,al,net,true);
    return plant::FF16State<double>{xad::value(r.height_dt),xad::value(r.mortality_dt),
      xad::value(r.fecundity_dt),xad::value(r.area_heartwood_dt),xad::value(r.mass_heartwood_dt)};
  };
  auto axpy=[](const plant::FF16State<double>&a,double c,const plant::FF16State<double>&k){
    return plant::FF16State<double>{a.height+c*k.height,a.mortality+c*k.mortality,
      a.fecundity+c*k.fecundity,a.area_heartwood+c*k.area_heartwood,a.mass_heartwood+c*k.mass_heartwood};};
  plant::FF16State<double> y{h0,0,0,0,0};
  return plant::ff16_cashkarp_replay(y, step_h, birth, deriv, axpy);
}

// Sensitivity pass for one trait: reads the recorded harvest (NO leaf opts).
static double sens_cohort(const plant::TF24ProdPars<double>& pd, const plant::TF24_Pars& pars,
    const std::string& trait, std::size_t birth, const std::vector<double>& step_h,
    double h0, double dh0, const std::vector<H>& rec) {
  std::size_t idx=0;
  auto deriv = [&](const TL& y, std::size_t, int) -> TL {
    const H& hh = rec[idx++];
    const double h=y.v.height, sh=y.s.height;
    const double dprofit_d = direct(trait,hh,pars) + hh.dprofit_dh*sh;
    plant::TF24ProdPars<F> pf=lift<F>(pd); seed_field(pf,trait);
    F h_ad=h; xad::derivative(h_ad)=sh;
    F prof=hh.profit; xad::derivative(prof)=dprofit_d;
    F al=plant::tf24_area_leaf<F>(pf.a_l1,pf.a_l2,h_ad);
    F net=plant::tf24_net_mass_production<F>(pf,h_ad,al,prof);
    plant::TF24Rates<F> r=plant::tf24_compute_rates_from_net<F>(pf,h_ad,al,net,true);
    TL o;
    o.v=plant::FF16State<double>{xad::value(r.height_dt),xad::value(r.mortality_dt),xad::value(r.fecundity_dt),xad::value(r.area_heartwood_dt),xad::value(r.mass_heartwood_dt)};
    o.s=plant::FF16State<double>{xad::derivative(r.height_dt),xad::derivative(r.mortality_dt),xad::derivative(r.fecundity_dt),xad::derivative(r.area_heartwood_dt),xad::derivative(r.mass_heartwood_dt)};
    return o;
  };
  auto axpy=[](const TL&a,double c,const TL&k)->TL{return TL{
    plant::FF16State<double>{a.v.height+c*k.v.height,a.v.mortality+c*k.v.mortality,a.v.fecundity+c*k.v.fecundity,a.v.area_heartwood+c*k.v.area_heartwood,a.v.mass_heartwood+c*k.v.mass_heartwood},
    plant::FF16State<double>{a.s.height+c*k.s.height,a.s.mortality+c*k.s.mortality,a.s.fecundity+c*k.s.fecundity,a.s.area_heartwood+c*k.s.area_heartwood,a.s.mass_heartwood+c*k.s.mass_heartwood}};};
  TL y{ plant::FF16State<double>{h0,0,0,0,0}, plant::FF16State<double>{dh0,0,0,0,0} };
  return plant::ff16_cashkarp_replay(y, step_h, birth, deriv, axpy).s.fecundity;
}

// double stand J (for FD of a single trait).
static double stand_J(plant::TF24_Strategy& s, const plant::TF24ProdPars<double>& pd,
    std::vector<std::vector<plant::TF24_Environment>>& EH, const std::vector<int>& birth,
    const std::vector<double>& step_h, double h0, const std::vector<double>& w) {
  std::vector<H> rec; double J=0;
  for (std::size_t i=0;i<birth.size();++i)
    J += w[i]*record_cohort(s,pd,EH,(std::size_t)birth[i],step_h,h0,rec).fecundity;
  return J;
}

static std::vector<std::vector<plant::TF24_Environment>> build_EH(Rcpp::List eh_list) {
  const std::size_t N=eh_list.size();
  std::vector<std::vector<plant::TF24_Environment>> EH(N);
  for (std::size_t n=0;n<N;++n){Rcpp::List st=eh_list[n];
    for(R_xlen_t k=0;k<st.size();++k) EH[n].push_back(Rcpp::as<plant::TF24_Environment>(st[k]));}
  return EH;
}

// The CHEAP part: the full 27-trait AD gradient in one shared-harvest pass (the
// leaf opts run once per cohort, all 27 sensitivities propagated through them).
// [[Rcpp::export]]
Rcpp::List tf24_emergent_ad(Rcpp::NumericVector pp, Rcpp::List eh_list,
    std::vector<double> shv, std::vector<int> birth, std::vector<double> w) {
  auto EH = build_EH(eh_list);
  const std::size_t N=eh_list.size();
  std::vector<double> step_h(N); for(std::size_t n=0;n<N;++n) step_h[n]=shv[n+1]-shv[n];
  plant::TF24_Strategy s0 = make_strategy(pp);
  plant::TF24ProdPars<double> pd = s0.prod_pars();
  const double h0 = s0.initial_height();
  const std::size_t T=TRAITS.size();
  std::vector<double> dh0(T,0.0);
  for (std::size_t k=0;k<T;++k){ double v0=trait_value(s0,TRAITS[k]); double dd=1e-4*std::abs(v0);
    dh0[k]=(make_strategy(pp,TRAITS[k],v0+dd).initial_height()
           -make_strategy(pp,TRAITS[k],v0-dd).initial_height())/(2*dd); }
  std::vector<double> sJ(T,0.0); double J=0;
  Rcpp::NumericVector hf(birth.size());
  std::vector<H> rec;
  for (std::size_t i=0;i<birth.size();++i){
    plant::FF16State<double> yf =
      record_cohort(s0,pd,EH,(std::size_t)birth[i],step_h,h0,rec);
    hf[i]=yf.height; J += w[i]*yf.fecundity;
    for (std::size_t k=0;k<T;++k)
      sJ[k] += w[i]*sens_cohort(pd,s0.pars,TRAITS[k],(std::size_t)birth[i],step_h,h0,dh0[k],rec);
  }
  return Rcpp::List::create(Rcpp::_["trait"]=TRAITS, Rcpp::_["J"]=J, Rcpp::_["grad"]=Rcpp::wrap(sJ),
    Rcpp::_["dh0"]=Rcpp::wrap(dh0), Rcpp::_["replay_heights"]=hf);
}

// The EXPENSIVE part for ONE trait: a two-pass live FD (re-runs the stand replay at
// the perturbed trait, with leaf opts). Independent per trait -> the R driver runs
// these across cores. Most traits resolve at a 1e-4 step; the photosynthesis
// COLIMITATION/electron-transport curvatures (curv_*) sit near theta~0.99 where the
// second derivative is large and amplified along the trajectory, so their emergent
// FD truncation needs a finer 1e-5 step (their AD is independently exact -- the
// net-level FD converges to ~1e-8). This per-trait step is exactly the asymmetry
// reverse-mode AD sidesteps.
// [[Rcpp::export]]
double tf24_emergent_fd_one(std::string trait, Rcpp::NumericVector pp, Rcpp::List eh_list,
    std::vector<double> shv, std::vector<int> birth, std::vector<double> w) {
  auto EH = build_EH(eh_list);
  const std::size_t N=eh_list.size();
  std::vector<double> step_h(N); for(std::size_t n=0;n<N;++n) step_h[n]=shv[n+1]-shv[n];
  plant::TF24_Strategy s0 = make_strategy(pp);
  const double rel = (trait=="curv_colim"||trait=="curv_elec") ? 1e-5 : 1e-4;
  double v0=trait_value(s0,trait), dd=rel*std::abs(v0);
  plant::TF24_Strategy sp=make_strategy(pp,trait,v0+dd);
  plant::TF24_Strategy sm=make_strategy(pp,trait,v0-dd);
  double Jp=stand_J(sp,sp.prod_pars(),EH,birth,step_h,sp.initial_height(),w);
  double Jm=stand_J(sm,sm.prod_pars(),EH,birth,step_h,sm.initial_height(),w);
  return (Jp-Jm)/(2*dd);
}')

# The AD 27-vector is a single cheap shared-harvest pass. The two-pass FD ground
# truth is one independent stand replay per trait, so it is run ACROSS CORES with
# parallel::mclapply (the sourceCpp functions are compiled once in the parent; forks
# inherit the loaded library). By default FD-validate one trait per class (fast);
# `Rscript ... full` FD-validates ALL 27 traits.
suppressMessages(library(parallel))
all27 <- c("vcmax_25","g1_TF24","beta2","K_s","b","c","jmax_25","a","curv_elec",
           "curv_colim","lma","rho","a_b1","r_l","r_b","r_s","r_r","k_l","k_b","k_s",
           "k_r","a_bio","a_y","a_l1","a_l2","theta","a_r1")
full_fd <- "full" %in% commandArgs(trailingOnly = TRUE)
fd_reps <- if (full_fd) all27 else c("vcmax_25", "lma", "a_l1", "theta", "a_y")

res <- tf24_emergent_ad(pp, eh, sh, birth_step, weights)
max_h_err <- max(abs(res$replay_heights - live_heights))
cat(sprintf("\nFaithfulness  max |replay - live SCM height| = %.2e\n", max_h_err))
cat(sprintf("J = sum_i w_i fecundity_i(t_end) = %.8g\n\n", res$J))

g <- setNames(res$grad, res$trait)
cat("Full 27-trait emergent gradient  d(J)/d(trait)  (ONE shared-harvest pass):\n")
for (t in res$trait) cat(sprintf("  %-11s % .7g\n", t, g[[t]]))

ncores <- max(1L, min(length(fd_reps), detectCores() - 2L))
cat(sprintf("\nFD validation (%s) across %d cores:\n",
            if (full_fd) "ALL 27 traits" else "representatives, one per class", ncores))
fd_vals <- unlist(mclapply(fd_reps, function(t)
  tf24_emergent_fd_one(t, pp, eh, sh, birth_step, weights), mc.cores = ncores))
fd <- setNames(fd_vals, fd_reps); ok <- TRUE
for (t in fd_reps) {
  re <- abs(g[[t]] - fd[[t]]) / max(abs(fd[[t]]), 1e-30)
  # cascade traits that also shift height_0 have a noisier two-pass FD ground truth
  # (stiff leaf-opt + seedling root-find); a leaf trait is clean. AD is exact either way.
  pass <- re < 1e-3; ok <- ok && pass
  cat(sprintf("  %-11s AD=% .7g  FD=% .7g  rel=%.2e %s\n",
              t, g[[t]], fd[[t]], re, if (pass) "OK" else "** CHECK **"))
}
stopifnot(max_h_err < 1e-5, ok)
cat(sprintf("\nFull 27-trait emergent gradient through the live TF24 SCM in one shared pass;\n%s validated vs two-pass FD.\n",
            if (full_fd) "all 27 traits" else "representatives of every class"))
