# Deep-crown live-SCM two-pass emergent trait gradient (#472 scope B, Milestone C).
#
# The DEFAULT FF16 shading model: the resident SCM integrates photosynthesis over
# crown depth with adaptive Gauss-Kronrod (FF16_Strategy::assimilation_deep_crown),
# not the single crown-top light read of scripts/ad_emergent_gradient.R.
#
#   Pass 1 (double): run the real resident SCM (default deep-crown shading) with
#     save_RK45_cache; harvest the frozen schedule (step_history) + per-RK-stage
#     resident light (environment_history).
#   Pass 2 (AD): replay each cohort with the MOVING-NODE Gauss-Kronrod crown integral
#     -- QK::integrate_ad over the active bounds [0, height], reading the FROZEN
#     per-stage resident light at each (moving) node with value + slope so d(light)/dz
#     flows -- exactly the deep-crown path of
#     FF16_Strategy::growth_rate_gradient_height_ad (#537 A1), carried through the
#     whole trajectory. The net feeds the SHARED rate-fill tail
#     (ff16_compute_rates_from_net), identical to crown-top downstream.
#
# Validation: (a) the double replay reproduces the live SCM cohort heights to machine
# precision; (b) d(J)/d(a_p1) for J = sum_i w_i * fecundity_i(t_end) matches a
# two-pass central finite difference.
#
# Run from the package root after `R CMD INSTALL .`:
#   Rscript scripts/ad_deep_crown_gradient.R
suppressMessages({library(Rcpp); library(plant)})

## Pass 1: real resident SCM with the DEFAULT deep-crown shading + cache.
p <- scm_base_parameters("FF16")
p <- add_strategies(p, trait_matrix(0.0825,"lma"), hyperpar=FF16_hyperpar, birth_rate=list(20))
p <- run_scm(p, Environment("FF16"), control(), refine_schedule=TRUE)$parameters
scm <- run_scm(p, Environment("FF16"), control(save_RK45_cache=TRUE), refine_schedule=FALSE)
stopifnot(!is.unsorted(scm$patch$step_history))

sh <- scm$patch$step_history
eh <- scm$patch$environment_history
sp <- scm$patch$species[[1]]
node_times <- sp$node_times; live_heights <- sp$heights; pdens <- sp$patch_densities
pp <- unlist(scm$parameters$strategies[[1]]$pars)
birth_step <- vapply(node_times, function(t) which.min(abs(sh - t)) - 1L, integer(1))
# trapezoid weights for an emergent J = sum_i tw_i * fecundity_i
tc <- numeric(length(node_times)); x <- node_times; nn <- length(x)
tc[1] <- 0.5*(x[2]-x[1]); tc[nn] <- 0.5*(x[nn]-x[nn-1])
if (nn>2) tc[2:(nn-1)] <- 0.5*(x[3:nn]-x[1:(nn-2)])
tw <- tc * pdens * 0.25 * 20            # * S_D * birth_rate (frozen)
cat(sprintf("Pass 1 (deep-crown): %d steps, %d cohorts, t_end=%.1f\n",
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

struct Frozen {
  std::vector<std::vector<plant::FF16_Environment>> eh;
  std::vector<double> step_h; double eta, h0;
  const plant::quadrature::QK* integ;
};

// Deep-crown derivative: net = area_leaf * GK_integral over [0,h] of
// assimilation_leaf(light(z)) * q(z/h, z), light read ACTIVELY from the frozen
// per-stage resident env (value + slope). Reuses ff16_compute_rates_from_net.
template <typename S>
static plant::FF16State<S> deep_crown_deriv(const plant::FF16ProdPars<S>& p,
    const Frozen& F, const plant::FF16State<S>& st, std::size_t n, int stage,
    bool mortality_finite) {
  const plant::FF16_Environment* e =
    (stage==0)?((n>0)?&F.eh[n-1][5]:&F.eh[0][0]):&F.eh[n][stage-1];
  const double canopy_top = e->max_environment_height();
  const S height = st.height;
  auto integrand = [&](S z) -> S {
    double zv = as_double(z);
    double lv = e->get_environment_at_height(zv, canopy_top);
    double ld = e->get_environment_deriv_at_height(zv);
    S light = S(lv) + S(std::isfinite(ld)?ld:0.0) * (z - S(zv));
    S u = z / height;
    return plant::ff16_assimilation_leaf<S>(p.a_p1, p.a_p2, light) *
           plant::ff16_canopy_q<S>(F.eta, u, z);
  };
  S area_leaf = plant::ff16_area_leaf(p.a_l1, p.a_l2, height);
  S assim = area_leaf * F.integ->integrate_ad<S>(integrand, S(0.0), height);
  S net = plant::ff16_net_from_components(p, height, area_leaf, assim);
  plant::FF16Rates<S> r = plant::ff16_compute_rates_from_net(p, height, area_leaf, net, mortality_finite);
  return plant::FF16State<S>{r.height_dt, r.mortality_dt, r.fecundity_dt,
                             r.area_heartwood_dt, r.mass_heartwood_dt};
}

template <typename S>
static plant::FF16State<S> replay_deep(const plant::FF16ProdPars<S>& p, const Frozen& F,
                                       std::size_t step0) {
  auto deriv = [&](const plant::FF16State<S>& st, std::size_t n, int stage){
    return deep_crown_deriv<S>(p, F, st, n, stage, true); };
  auto axpy = [](const plant::FF16State<S>& a, double c, const plant::FF16State<S>& k){
    return plant::FF16State<S>{a.height+c*k.height,a.mortality+c*k.mortality,a.fecundity+c*k.fecundity,
      a.area_heartwood+c*k.area_heartwood,a.mass_heartwood+c*k.mass_heartwood}; };
  plant::FF16State<S> y{S(F.h0), S(0),S(0),S(0),S(0)};
  return plant::ff16_cashkarp_replay(y, F.step_h, step0, deriv, axpy);
}

// [[Rcpp::export]]
Rcpp::List deep_crown(Rcpp::NumericVector pp, Rcpp::List eh_list, std::vector<double> sh,
                      std::vector<int> birth, std::vector<double> tw) {
  auto s = make_strategy(pp); auto pd = s.prod_pars();
  Frozen F; F.eta = s.pars.eta; F.h0 = s.initial_height(); F.integ = &s.function_integrator;
  const std::size_t N = eh_list.size(); F.eh.resize(N); F.step_h.resize(N);
  for(std::size_t n=0;n<N;++n){Rcpp::List st=eh_list[n]; for(R_xlen_t k=0;k<st.size();++k) F.eh[n].push_back(Rcpp::as<plant::FF16_Environment>(st[k]));}
  for(std::size_t n=0;n<N;++n) F.step_h[n]=sh[n+1]-sh[n];

  // (a) faithfulness: double replay final heights.
  Rcpp::NumericVector hf(birth.size());
  for(std::size_t i=0;i<birth.size();++i) hf[i]=replay_deep<double>(pd, F, (std::size_t)birth[i]).height;

  // emergent J = sum_i tw_i * fecundity_i, double + AD + FD.
  auto standJ=[&](const plant::FF16ProdPars<double>& q)->double{ double J=0; for(std::size_t i=0;i<birth.size();++i) J+=tw[i]*replay_deep<double>(q,F,(std::size_t)birth[i]).fecundity; return J; };
  double Jd = standJ(pd);
  double dJ_ad;
  { ad::tape_type tape; ad_t a_p1=pd.a_p1; tape.registerInput(a_p1); tape.newRecording();
    auto pa=lift<ad_t>(pd); pa.a_p1=a_p1;
    ad_t J=ad_t(0.0);
    for(std::size_t i=0;i<birth.size();++i) J += tw[i]*replay_deep<ad_t>(pa,F,(std::size_t)birth[i]).fecundity;
    tape.registerOutput(J); xad::derivative(J)=1.0; tape.computeAdjoints(); dJ_ad=xad::derivative(a_p1); }
  std::vector<double> rel_h={1e-4,1e-5,1e-6,1e-7}, fd;
  for(double rh:rel_h){ double h=rh*pd.a_p1; auto q1=pd; q1.a_p1+=h; auto q2=pd; q2.a_p1-=h; fd.push_back((standJ(q1)-standJ(q2))/(2*h)); }
  return Rcpp::List::create(Rcpp::_["replay_heights"]=hf, Rcpp::_["J"]=Jd,
                            Rcpp::_["dJ_ad"]=dJ_ad, Rcpp::_["fd"]=Rcpp::wrap(fd), Rcpp::_["rel_h"]=Rcpp::wrap(rel_h));
}')

res <- deep_crown(pp, eh, sh, birth_step, tw)
max_h_err <- max(abs(res$replay_heights - live_heights))
cat(sprintf("\n(a) deep-crown faithfulness: max |replay - live SCM height| over %d cohorts = %.2e\n",
            length(live_heights), max_h_err))
cat(sprintf("\n(b) emergent J (sum w_i*fecundity_i) = %.8g\n", res$J))
for (k in seq_along(res$rel_h))
  cat(sprintf("    h/a_p1=%.0e  FD=%.9g  rel.err=%.2e\n", res$rel_h[k], res$fd[k], abs(res$fd[k]-res$dJ_ad)/abs(res$dJ_ad)))
best <- min(abs(res$fd-res$dJ_ad)/abs(res$dJ_ad))
cat(sprintf("\n    d(J)/d(a_p1): AD=%.9g  best FD=%.9g  min rel.err=%.2e  %s\n",
            res$dJ_ad, res$fd[which.min(abs(res$fd-res$dJ_ad))], best, if (best<1e-5) "OK" else "** MISMATCH **"))
stopifnot(max_h_err < 1e-7, best < 1e-5)
cat("\nDeep-crown two-pass emergent gradient validated.\n")
