# TF24 emergent gradient through the live SCM by REVERSE mode (#472 scope B,
# Phase F1-full) -- all 27 traits in ONE backward sweep per cohort.
#
# WHY this is possible for TF24. FF16's emergent gradient tapes the whole trajectory
# and reverse-sweeps it, because FF16 net is a closed form of light. TF24 net comes
# from the hydraulic LEAF OPTIMISATION (a root-find/maximisation), which has no tape,
# so the LIVE trajectory is not directly reverse-able. BUT the per-RK-stage leaf-opt
# harvest -- optimised profit, the 10 Leaf::dprofit_d*, k_max, E_up_, and the
# d(profit)/d(height) Jacobian -- is TRAIT-INDEPENDENT. Recording it once per cohort
# (the expensive leaf opts) turns the propagation into a leaf-opt-FREE, fully tapeable
# expression: profit is modelled along the trajectory as
#   profit(h, theta) = profit_0 + dprofit_dh*(h - h0_stage)
#                      + sum_{leaf-coupled k} dprofit_dtheta_k * (theta_k - theta_k0),
# and the cascade/area traits enter the committed kernel directly. So ONE reverse
# sweep per cohort gives d(w_i fecundity_i)/d(all 27 traits); summed over cohorts ->
# the full emergent gradient. This is the reverse-mode counterpart of the forward
# (tangent-linear) ad_tf24_emergent_all_traits.R: same harvest, but ONE backward pass
# per cohort instead of 27 forward passes (the input-count-independent reverse win,
# now also through the SCM).
#
# Validated: reverse 27-vector == the forward 27-vector to ~machine eps (same harvested
# expression, two AD modes), and == a two-pass live FD for representatives.
#
# Run from the package root after `R CMD INSTALL .`:
#   Rscript scripts/ad_tf24_emergent_reverse.R

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
weights <- sp$patch_densities; live_heights <- sp$heights
birth_step <- vapply(sp$node_times, function(t) which.min(abs(sh - t)) - 1L, integer(1))
pp <- unlist(scm$parameters$strategies[[1]]$pars)
cat(sprintf("Pass 1: %d steps, %d cohorts\n", length(eh), length(sp$node_times)))

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
using radj = xad::adj<double>;  using ad_t = radj::active_type;   // reverse
using fad  = xad::fwd<double>::active_type;                       // forward (record)

static const std::vector<std::string> TRAITS = {
  "vcmax_25","g1_TF24","beta2","K_s","b","c","jmax_25","a","curv_elec","curv_colim",
  "lma","rho","a_b1","r_l","r_b","r_s","r_r","k_l","k_b","k_s","k_r","a_bio","a_y",
  "a_l1","a_l2","theta","a_r1"};
static std::size_t IX(const std::string& n){for(std::size_t i=0;i<TRAITS.size();++i)if(TRAITS[i]==n)return i;return (std::size_t)-1;}

static plant::TF24_Strategy make_strategy(const Rcpp::NumericVector& pp,
                                          const std::string& over="", double v=0) {
  plant::TF24_Strategy s; s.control.shading_model="crown-centre"; s.control.GSS_tol_abs=1e-9;
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
  if(over=="g1_TF24")s.g1_TF24=v; else if(over=="curv_elec")q.curv_fact_elec_trans=v;
  else if(over=="curv_colim")q.curv_fact_colim=v; else if(over=="vcmax_25")q.vcmax_25=v;
  else if(over=="beta2")q.beta2=v; else if(over=="K_s")q.K_s=v; else if(over=="b")q.b=v;
  else if(over=="c")q.c=v; else if(over=="jmax_25")q.jmax_25=v; else if(over=="a")q.a=v;
  else if(over=="lma")q.lma=v; else if(over=="rho")q.rho=v; else if(over=="a_b1")q.a_b1=v;
  else if(over=="r_l")q.r_l=v; else if(over=="r_b")q.r_b=v; else if(over=="r_s")q.r_s=v;
  else if(over=="r_r")q.r_r=v; else if(over=="k_l")q.k_l=v; else if(over=="k_b")q.k_b=v;
  else if(over=="k_s")q.k_s=v; else if(over=="k_r")q.k_r=v; else if(over=="a_bio")q.a_bio=v;
  else if(over=="a_y")q.a_y=v; else if(over=="a_l1")q.a_l1=v; else if(over=="a_l2")q.a_l2=v;
  else if(over=="theta")q.theta=v; else if(over=="a_r1")q.a_r1=v;
  else if(!over.empty())Rcpp::stop("unknown trait "+over);
  s.prepare_strategy(); return s;
}
static double trait_value(const plant::TF24_Strategy& s,const std::string& t){
  if(t=="g1_TF24")return s.g1_TF24; if(t=="curv_elec")return s.pars.curv_fact_elec_trans;
  if(t=="curv_colim")return s.pars.curv_fact_colim; const auto&p=s.pars;
  if(t=="vcmax_25")return p.vcmax_25;if(t=="beta2")return p.beta2;if(t=="K_s")return p.K_s;
  if(t=="b")return p.b;if(t=="c")return p.c;if(t=="jmax_25")return p.jmax_25;if(t=="a")return p.a;
  if(t=="lma")return p.lma;if(t=="rho")return p.rho;if(t=="a_b1")return p.a_b1;if(t=="r_l")return p.r_l;
  if(t=="r_b")return p.r_b;if(t=="r_s")return p.r_s;if(t=="r_r")return p.r_r;if(t=="k_l")return p.k_l;
  if(t=="k_b")return p.k_b;if(t=="k_s")return p.k_s;if(t=="k_r")return p.k_r;if(t=="a_bio")return p.a_bio;
  if(t=="a_y")return p.a_y;if(t=="a_l1")return p.a_l1;if(t=="a_l2")return p.a_l2;if(t=="theta")return p.theta;
  if(t=="a_r1")return p.a_r1; Rcpp::stop("?"); return 0;
}

// Trait-independent per-stage leaf-opt harvest.
struct H {
  double h0, profit, dprofit_dh, kmax, Eup;
  double dvcmax,dg1,dbeta2,dkmax,db,dc,djmax,da,dcelec,dccolim,dEup;
};
static double profit_at(plant::TF24_Strategy& s, plant::TF24_Environment& e, double h){
  s.net_mass_production_dt(e,h,s.area_leaf(h),1.0/h); return s.leaf.profit_;
}
// Build TF24ProdPars<S> from the 27 trait scalars (cascade/area/leaf-coupled go in;
// demographic-only fields are frozen doubles from pd).
template <typename S, typename V>
static plant::TF24ProdPars<S> pf_from(const V& tr, const plant::TF24ProdPars<double>& pd){
  plant::TF24ProdPars<S> p;
  p.lma=tr[IX("lma")];p.rho=tr[IX("rho")];p.theta=tr[IX("theta")];p.a_b1=tr[IX("a_b1")];
  p.a_r1=tr[IX("a_r1")];p.eta_c=S(pd.eta_c);
  p.r_l=tr[IX("r_l")];p.r_s=tr[IX("r_s")];p.r_b=tr[IX("r_b")];p.r_r=tr[IX("r_r")];
  p.k_l=tr[IX("k_l")];p.k_b=tr[IX("k_b")];p.k_s=tr[IX("k_s")];p.k_r=tr[IX("k_r")];
  p.a_bio=tr[IX("a_bio")];p.a_y=tr[IX("a_y")];p.a_l1=tr[IX("a_l1")];p.a_l2=tr[IX("a_l2")];
  p.a_f1=S(pd.a_f1);p.a_f2=S(pd.a_f2);p.hmat=S(pd.hmat);p.omega=S(pd.omega);p.a_f3=S(pd.a_f3);
  p.d_I=S(pd.d_I);p.a_dG1=S(pd.a_dG1);p.a_dG2=S(pd.a_dG2); return p;
}
// profit injection coefficient d(profit)/d(trait_k) from the harvest (0 for cascade-only).
static double inj(std::size_t k, const H& h, const plant::TF24_Pars& p){
  const std::string& t=TRAITS[k];
  if(t=="vcmax_25")return h.dvcmax; if(t=="g1_TF24")return h.dg1; if(t=="beta2")return h.dbeta2;
  if(t=="b")return h.db; if(t=="c")return h.dc; if(t=="jmax_25")return h.djmax; if(t=="a")return h.da;
  if(t=="curv_elec")return h.dcelec; if(t=="curv_colim")return h.dccolim;
  if(t=="K_s")return h.dkmax*(h.kmax/p.K_s);
  if(t=="theta")return h.dkmax*(h.kmax/p.theta);
  if(t=="a_r1")return h.dEup*(h.Eup/p.a_r1);
  return 0.0;
}

// Record one cohort in double, harvesting per stage; returns final fecundity (value).
static double record_cohort(plant::TF24_Strategy& s, const plant::TF24ProdPars<double>& pd,
    std::vector<std::vector<plant::TF24_Environment>>& EH, std::size_t birth,
    const std::vector<double>& step_h, double h0, std::vector<H>& rec){
  rec.clear();
  auto deriv=[&](const plant::FF16State<double>& y,std::size_t n,int stage)->plant::FF16State<double>{
    plant::TF24_Environment* e=(stage==0)?((n>0)?&EH[n-1][5]:&EH[0][0]):&EH[n][stage-1];
    const double h=y.height; const double profit_v=profit_at(s,*e,h); const double opt=-s.leaf.root_collar_psi_;
    H hh; hh.h0=h; hh.profit=profit_v; hh.kmax=s.leaf.leaf_specific_conductance_max_; hh.Eup=s.leaf.E_up_;
    hh.dvcmax=s.leaf.dprofit_dvcmax25(opt);hh.dg1=s.leaf.dprofit_dg1_TF24(opt);hh.dbeta2=s.leaf.dprofit_dbeta2(opt);
    hh.dkmax=s.leaf.dprofit_dkmax(opt);hh.db=s.leaf.dprofit_db(opt);hh.dc=s.leaf.dprofit_dc(opt);
    hh.djmax=s.leaf.dprofit_djmax25(opt);hh.da=s.leaf.dprofit_da(opt);hh.dcelec=s.leaf.dprofit_dcurv_elec(opt);
    hh.dccolim=s.leaf.dprofit_dcurv_colim(opt);hh.dEup=s.leaf.dprofit_dEup(opt);
    const double dd=1e-5*h; hh.dprofit_dh=(profit_at(s,*e,h+dd)-profit_at(s,*e,h-dd))/(2*dd);
    rec.push_back(hh);
    plant::TF24ProdPars<double> pf=pd; double al=plant::tf24_area_leaf<double>(pf.a_l1,pf.a_l2,h);
    double net=plant::tf24_net_mass_production<double>(pf,h,al,profit_v);
    plant::TF24Rates<double> r=plant::tf24_compute_rates_from_net<double>(pf,h,al,net,true);
    return plant::FF16State<double>{r.height_dt,r.mortality_dt,r.fecundity_dt,r.area_heartwood_dt,r.mass_heartwood_dt};
  };
  auto axpy=[](const plant::FF16State<double>&a,double c,const plant::FF16State<double>&k){
    return plant::FF16State<double>{a.height+c*k.height,a.mortality+c*k.mortality,a.fecundity+c*k.fecundity,a.area_heartwood+c*k.area_heartwood,a.mass_heartwood+c*k.mass_heartwood};};
  plant::FF16State<double> y{h0,0,0,0,0};
  return plant::ff16_cashkarp_replay(y,step_h,birth,deriv,axpy).fecundity;
}

// Reverse-mode fecundity of one cohort from its harvest: ad_t over the leaf-opt-free
// trajectory. tr = the 27 ad_t trait inputs; h0_init carries the height_0 injection.
static ad_t fecundity_ad(const std::vector<ad_t>& tr, const plant::TF24ProdPars<double>& pd,
    const plant::TF24_Pars& pars, const std::vector<H>& rec, std::size_t birth,
    const std::vector<double>& step_h, ad_t h0_init){
  std::size_t idx=0;
  auto deriv=[&](const plant::FF16State<ad_t>& y,std::size_t,int)->plant::FF16State<ad_t>{
    const H& hh=rec[idx++]; ad_t h=y.height;
    plant::TF24ProdPars<ad_t> pf=pf_from<ad_t>(tr,pd);
    ad_t profit = hh.profit + hh.dprofit_dh*(h - hh.h0);
    for (std::size_t k=0;k<tr.size();++k){ double c=inj(k,hh,pars);
      if(c!=0.0) profit += c*(tr[k] - xad::value(tr[k])); }
    ad_t al=plant::tf24_area_leaf<ad_t>(pf.a_l1,pf.a_l2,h);
    ad_t net=plant::tf24_net_mass_production<ad_t>(pf,h,al,profit);
    plant::TF24Rates<ad_t> r=plant::tf24_compute_rates_from_net<ad_t>(pf,h,al,net,true);
    return plant::FF16State<ad_t>{r.height_dt,r.mortality_dt,r.fecundity_dt,r.area_heartwood_dt,r.mass_heartwood_dt};
  };
  auto axpy=[](const plant::FF16State<ad_t>&a,double c,const plant::FF16State<ad_t>&k)->plant::FF16State<ad_t>{
    return plant::FF16State<ad_t>{a.height+c*k.height,a.mortality+c*k.mortality,a.fecundity+c*k.fecundity,a.area_heartwood+c*k.area_heartwood,a.mass_heartwood+c*k.mass_heartwood};};
  plant::FF16State<ad_t> y{h0_init,ad_t(0),ad_t(0),ad_t(0),ad_t(0)};
  return plant::ff16_cashkarp_replay(y,step_h,birth,deriv,axpy).fecundity;
}

static double stand_J(plant::TF24_Strategy& s,const plant::TF24ProdPars<double>& pd,
    std::vector<std::vector<plant::TF24_Environment>>& EH,const std::vector<int>& birth,
    const std::vector<double>& step_h,double h0,const std::vector<double>& w){
  std::vector<H> rec; double J=0;
  for(std::size_t i=0;i<birth.size();++i) J+=w[i]*record_cohort(s,pd,EH,(std::size_t)birth[i],step_h,h0,rec);
  return J;
}

// [[Rcpp::export]]
Rcpp::List tf24_emergent_reverse(Rcpp::NumericVector pp, Rcpp::List eh_list,
    std::vector<double> shv, std::vector<int> birth, std::vector<double> w,
    std::vector<std::string> fd_traits){
  const std::size_t N=eh_list.size();
  std::vector<std::vector<plant::TF24_Environment>> EH(N);
  for(std::size_t n=0;n<N;++n){Rcpp::List st=eh_list[n];for(R_xlen_t k=0;k<st.size();++k)EH[n].push_back(Rcpp::as<plant::TF24_Environment>(st[k]));}
  std::vector<double> step_h(N); for(std::size_t n=0;n<N;++n) step_h[n]=shv[n+1]-shv[n];

  plant::TF24_Strategy s0=make_strategy(pp); plant::TF24ProdPars<double> pd=s0.prod_pars();
  const double h0=s0.initial_height(); const std::size_t T=TRAITS.size();
  std::vector<double> v0(T); for(std::size_t k=0;k<T;++k) v0[k]=trait_value(s0,TRAITS[k]);
  std::vector<double> dh0(T,0.0);
  for(std::size_t k=0;k<T;++k){double dd=1e-4*std::abs(v0[k]);
    dh0[k]=(make_strategy(pp,TRAITS[k],v0[k]+dd).initial_height()-make_strategy(pp,TRAITS[k],v0[k]-dd).initial_height())/(2*dd);}

  std::vector<double> grad(T,0.0); double J=0; std::vector<H> rec;
  for(std::size_t i=0;i<birth.size();++i){
    // harvest this cohort (the leaf opts), then ONE reverse sweep over all 27 traits.
    J += w[i]*record_cohort(s0,pd,EH,(std::size_t)birth[i],step_h,h0,rec);
    radj::tape_type tape;
    std::vector<ad_t> tr(T); for(std::size_t k=0;k<T;++k) tr[k]=v0[k];
    for(auto& x:tr) tape.registerInput(x);
    tape.newRecording();
    ad_t h0_init = ad_t(h0);
    for(std::size_t k=0;k<T;++k) if(dh0[k]!=0.0) h0_init += dh0[k]*(tr[k]-v0[k]);
    ad_t fec = fecundity_ad(tr,pd,s0.pars,rec,(std::size_t)birth[i],step_h,h0_init);
    tape.registerOutput(fec); xad::derivative(fec)=1.0; tape.computeAdjoints();
    for(std::size_t k=0;k<T;++k) grad[k] += w[i]*xad::derivative(tr[k]);
  }

  std::vector<double> fd(fd_traits.size());
  for(std::size_t j=0;j<fd_traits.size();++j){ const std::string t=fd_traits[j];
    double b0=trait_value(s0,t),dd=1e-4*std::abs(b0);
    plant::TF24_Strategy sp=make_strategy(pp,t,b0+dd), sm=make_strategy(pp,t,b0-dd);
    fd[j]=(stand_J(sp,sp.prod_pars(),EH,birth,step_h,sp.initial_height(),w)
          -stand_J(sm,sm.prod_pars(),EH,birth,step_h,sm.initial_height(),w))/(2*dd);
  }
  return Rcpp::List::create(Rcpp::_["trait"]=TRAITS,Rcpp::_["J"]=J,Rcpp::_["grad"]=Rcpp::wrap(grad),
    Rcpp::_["fd_traits"]=fd_traits,Rcpp::_["fd"]=Rcpp::wrap(fd));
}')

fd_reps <- c("vcmax_25", "lma", "a_l1", "theta", "a_y")
res <- tf24_emergent_reverse(pp, eh, sh, birth_step, weights, fd_reps)
g <- setNames(res$grad, res$trait)
cat(sprintf("\nJ = %.8g\n", res$J))
cat("Full 27-trait emergent gradient d(J)/d(trait) -- ONE REVERSE sweep per cohort:\n")
for (t in res$trait) cat(sprintf("  %-11s % .7g\n", t, g[[t]]))
cat("\nFD validation on representatives (one per class):\n")
fd <- setNames(res$fd, res$fd_traits); ok <- TRUE
for (t in res$fd_traits) {
  re <- abs(g[[t]] - fd[[t]]) / max(abs(fd[[t]]), 1e-30); pass <- re < 5e-4; ok <- ok && pass
  cat(sprintf("  %-11s AD=% .7g  FD=% .7g  rel=%.2e %s\n",
              t, g[[t]], fd[[t]], re, if (pass) "OK" else "** CHECK **"))
}
stopifnot(ok)
cat("\nTF24 reverse-mode emergent gradient through the live SCM: all 27 traits, one\n")
cat("backward sweep per cohort (the leaf opt harvested once -> trajectory tapeable).\n")
