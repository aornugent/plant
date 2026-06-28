// Reverse-mode AD emergent-gradient routines for FF16, compiled into plant.so
// (#472 scope B, Phase C). Unlike growth_rate_gradient_height_ad (forward mode,
// header-only), these use the XAD adjoint TAPE (xad::adj). plant.so has no tape of
// its own: the tape symbols are odelia's single compiled copy (src/Tape.cpp),
// resolved at load against odelia's globally-loaded DLL -- the same mechanism the
// odelia ODE Solver already relies on (see src/Makevars; odelia is imported first
// via importFrom(odelia, odelia_load_dll), so its DLL loads before plant's).
//
// The headline routine differentiates the SCM's emergent offspring_production
// w.r.t. a set of FF16 traits in ONE reverse sweep, over the frozen resident
// schedule + per-RK-stage resident light harvested by a save_RK45_cache run
// (deep-crown / default shading). Establishment is frozen (a separable partial).
// It takes the harvested data as plain arrays, so it carries no templating into the
// SCM class; the R-facing offspring_production_gradient() gathers these from a run
// SCM and calls it.
#include <Rcpp.h>
#include <vector>
#include <string>
#include <algorithm>
#include <cmath>
#include <XAD/XAD.hpp>
#include <plant.h>                 // RcppR6 as<>/wrap for FF16_Environment etc.
#include <plant/models/ff16_production_kernel.h>

using ad   = xad::adj<double>;
using ad_t = ad::active_type;

namespace {

double as_double(double v)      { return v; }
double as_double(const ad_t& v) { return xad::value(v); }

plant::FF16_Strategy make_strategy(const Rcpp::NumericVector& pp) {
  plant::FF16_Strategy s; auto& q = s.pars;
  q.lma=pp["lma"];q.rho=pp["rho"];q.hmat=pp["hmat"];q.omega=pp["omega"];q.eta=pp["eta"];
  q.theta=pp["theta"];q.a_l1=pp["a_l1"];q.a_l2=pp["a_l2"];q.a_r1=pp["a_r1"];q.a_b1=pp["a_b1"];
  q.r_s=pp["r_s"];q.r_b=pp["r_b"];q.r_r=pp["r_r"];q.r_l=pp["r_l"];q.a_y=pp["a_y"];q.a_bio=pp["a_bio"];
  q.k_l=pp["k_l"];q.k_b=pp["k_b"];q.k_s=pp["k_s"];q.k_r=pp["k_r"];q.a_p1=pp["a_p1"];q.a_p2=pp["a_p2"];
  q.a_f3=pp["a_f3"];q.a_f1=pp["a_f1"];q.a_f2=pp["a_f2"];q.S_D=pp["S_D"];q.a_d0=pp["a_d0"];q.d_I=pp["d_I"];
  q.a_dG1=pp["a_dG1"];q.a_dG2=pp["a_dG2"];q.k_I=pp["k_I"];q.recruitment_decay=pp["recruitment_decay"];
  s.prepare_strategy(); return s;
}
template <typename S> plant::FF16ProdPars<S> lift(const plant::FF16ProdPars<double>& d) {
  plant::FF16ProdPars<S> p;
  p.lma=d.lma;p.rho=d.rho;p.theta=d.theta;p.a_b1=d.a_b1;p.a_r1=d.a_r1;p.eta_c=d.eta_c;
  p.a_p1=d.a_p1;p.a_p2=d.a_p2;p.r_l=d.r_l;p.r_s=d.r_s;p.r_b=d.r_b;p.r_r=d.r_r;
  p.k_l=d.k_l;p.k_b=d.k_b;p.k_s=d.k_s;p.k_r=d.k_r;p.a_bio=d.a_bio;p.a_y=d.a_y;
  p.a_l1=d.a_l1;p.a_l2=d.a_l2;p.a_f1=d.a_f1;p.a_f2=d.a_f2;p.hmat=d.hmat;
  p.omega=d.omega;p.a_f3=d.a_f3;p.d_I=d.d_I;p.a_dG1=d.a_dG1;p.a_dG2=d.a_dG2; return p;
}
template <typename S> std::vector<S*> field_ptrs(plant::FF16ProdPars<S>& p) {
  return {&p.lma,&p.rho,&p.theta,&p.a_b1,&p.a_r1,&p.eta_c,&p.a_p1,&p.a_p2,
          &p.r_l,&p.r_s,&p.r_b,&p.r_r,&p.k_l,&p.k_b,&p.k_s,&p.k_r,&p.a_bio,&p.a_y,
          &p.a_l1,&p.a_l2,&p.a_f1,&p.a_f2,&p.hmat,&p.omega,&p.a_f3,&p.d_I,&p.a_dG1,&p.a_dG2};
}
std::vector<std::string> field_names() {
  return {"lma","rho","theta","a_b1","a_r1","eta_c","a_p1","a_p2","r_l","r_s","r_b","r_r",
          "k_l","k_b","k_s","k_r","a_bio","a_y","a_l1","a_l2","a_f1","a_f2","hmat","omega",
          "a_f3","d_I","a_dG1","a_dG2"};
}

struct Frozen {
  std::vector<std::vector<plant::FF16_Environment>> eh;  // [step][0..5]
  std::vector<double> step_h, ppsab, tw, decay;          // decay = exp(-recr_decay*t_birth)
  Rcpp::NumericMatrix ppsurv;                            // [step][0..5] stage survival
  std::vector<int> birth;
  double eta, h0, a_d0;
  const plant::quadrature::QK* integ;
};

// Deep-crown net at `height` reading the frozen env `e` (moving-node GK integral).
template <typename S>
S deep_net(const plant::FF16ProdPars<S>& pd, const plant::quadrature::QK* integ,
           double eta, const plant::FF16_Environment* e, S height) {
  const double canopy_top = e->max_environment_height();
  auto integrand = [&](S z) -> S {
    double zv = as_double(z);
    double lv = e->get_environment_at_height(zv, canopy_top);
    double ld = e->get_environment_deriv_at_height(zv);
    S light = S(lv) + S(std::isfinite(ld)?ld:0.0) * (z - S(zv));
    return plant::ff16_assimilation_leaf<S>(pd.a_p1, pd.a_p2, light) *
           plant::ff16_canopy_q<S>(eta, z / height, z);
  };
  S area_leaf = plant::ff16_area_leaf(pd.a_l1, pd.a_l2, height);
  S assim = area_leaf * integ->integrate_ad<S>(integrand, S(0.0), height);
  return plant::ff16_net_from_components(pd, height, area_leaf, assim);
}

// Emergent offspring_production = sum_i tw_i * offspring_weighted_i (deep-crown
// 6-state replay through ff16_cashkarp_replay); establishment frozen via mort0.
template <typename S>
S stand_offspring(const plant::FF16ProdPars<S>& pd, const Frozen& F, S h0) {
  using std::exp; using std::log;
  S J = S(0.0);
  for (std::size_t i = 0; i < F.birth.size(); ++i) {
    const std::size_t b = (std::size_t)F.birth[i];
    const double ppsab = F.ppsab[i];
    // Establishment (recruitment filter), DIFFERENTIATED: mortality_0 =
    // -log(pr_estab), pr_estab from the seedling net production (deep-crown) in the
    // frozen birth env -> the trait flows through net0 and area_leaf_0. h0 (seedling
    // height) carries its own d/d(trait) via the IFT injection (see the caller).
    const plant::FF16_Environment* eb = (b > 0) ? &F.eh[b - 1][5] : &F.eh[0][0];
    S area_leaf_0 = plant::ff16_area_leaf(pd.a_l1, pd.a_l2, h0);
    S net0 = deep_net<S>(pd, F.integ, F.eta, eb, h0);
    S pr_estab = plant::ff16_establishment_probability<S>(area_leaf_0, net0, F.a_d0, F.decay[i]);
    S mort0 = -log(pr_estab);
    auto deriv = [&](const plant::FF16LifeState<S>& s, std::size_t n, int stage)
        -> plant::FF16LifeState<S> {
      const plant::FF16_Environment* e =
        (stage==0)?((n>0)?&F.eh[n-1][5]:&F.eh[0][0]):&F.eh[n][stage-1];
      S net = deep_net<S>(pd, F.integ, F.eta, e, s.demog.height);
      S area_leaf = plant::ff16_area_leaf(pd.a_l1, pd.a_l2, s.demog.height);
      plant::FF16Rates<S> r = plant::ff16_compute_rates_from_net(pd, s.demog.height, area_leaf, net, true);
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
    plant::FF16LifeState<S> y{plant::FF16State<S>{h0, mort0, S(0), S(0), S(0)}, S(0)};
    y = plant::ff16_cashkarp_replay(y, F.step_h, b, deriv, axpy);
    J += S(F.tw[i]) * y.offspring;
  }
  return J;
}

} // namespace

// Reverse-mode probe (CI smoke test of the tape-at-load): d(fecundity_dt)/d(a_p1)
// of a crown-top plant of the given height/light, via one backward sweep.
// [[Rcpp::export]]
double ff16_reverse_tape_probe(double height, double light_E) {
  plant::FF16_Strategy s; s.control.shading_model = "crown-centre"; s.prepare_strategy();
  auto pd = s.prod_pars();
  double d_ap1;
  { ad::tape_type tape; ad_t a_p1 = pd.a_p1; tape.registerInput(a_p1); tape.newRecording();
    auto p = lift<ad_t>(pd); p.a_p1 = a_p1;
    ad_t f = plant::ff16_compute_rates_crown_top<ad_t>(p, ad_t(height), ad_t(light_E), true).fecundity_dt;
    tape.registerOutput(f); xad::derivative(f) = 1.0; tape.computeAdjoints();
    d_ap1 = xad::derivative(a_p1); }
  return d_ap1;
}

// Compiled core of offspring_production_gradient(). Takes the harvested resident
// schedule (env per RK stage, step sizes), the per-cohort birth steps / weights /
// survival, and the trait names to differentiate. Returns d(offspring_production)/
// d(trait), establishment frozen. eh_list is steps x 6 FF16_Environment objects.
// [[Rcpp::export]]
Rcpp::NumericVector ff16_offspring_production_gradient_impl(
    Rcpp::NumericVector pp, Rcpp::List eh_list, std::vector<double> sh,
    std::vector<int> birth, Rcpp::NumericMatrix ppsurv, std::vector<double> ppsab,
    std::vector<double> tw, std::vector<std::string> traits) {
  auto s = make_strategy(pp);
  auto pd = s.prod_pars();
  Frozen F; F.eta = s.pars.eta; F.h0 = s.initial_height(); F.birth = birth;
  F.ppsab = ppsab; F.tw = tw; F.ppsurv = ppsurv; F.integ = &s.function_integrator;
  F.a_d0 = s.pars.a_d0;
  const std::size_t N = eh_list.size(); F.eh.resize(N); F.step_h.resize(N);
  for (std::size_t n=0;n<N;++n){Rcpp::List st=eh_list[n]; for(R_xlen_t k=0;k<st.size();++k) F.eh[n].push_back(Rcpp::as<plant::FF16_Environment>(st[k]));}
  for (std::size_t n=0;n<N;++n) F.step_h[n]=sh[n+1]-sh[n];

  // Per-cohort establishment decay factor exp(-recruitment_decay * birth_time);
  // birth_time = step_history at the birth step. recruitment_decay / a_d0 are not
  // among the differentiated traits, so they fold to doubles; the establishment
  // gradient flows through the (active) seedling net production inside stand_offspring.
  F.decay.resize(birth.size());
  for (std::size_t i=0;i<birth.size();++i)
    F.decay[i] = std::exp(-s.pars.recruitment_decay * sh[(std::size_t)birth[i]]);

  // Map requested trait names to FF16ProdPars field indices.
  auto names = field_names();
  std::vector<std::size_t> idx;
  for (auto& t : traits) {
    auto it = std::find(names.begin(), names.end(), t);
    if (it == names.end()) Rcpp::stop("unknown FF16 trait: " + t);
    idx.push_back(std::distance(names.begin(), it));
  }

  // d(height_0)/d(trait) by the implicit function theorem at the height_seed root
  // (mass_live_given_height(h0) == omega). A separate reverse sweep of mass_live at
  // h0 gives d(mass_live)/d(theta_k) and d(mass_live)/d(height); IFT then gives
  // d(h0)/d(theta_k) (the seedling-size response, the #539 pattern). Scoped in its
  // own block so its tape is torn down before the main recording.
  const double h0v = s.initial_height();
  std::vector<double> dh0(idx.size(), 0.0);
  {
    ad::tape_type tape0;
    auto pm = lift<ad_t>(pd);
    auto fm = field_ptrs<ad_t>(pm);
    ad_t hin = h0v;
    for (auto i : idx) tape0.registerInput(*fm[i]);
    tape0.registerInput(hin);
    tape0.newRecording();
    ad_t m = plant::ff16_mass_live_given_height<ad_t>(pm, hin);
    tape0.registerOutput(m); xad::derivative(m) = 1.0; tape0.computeAdjoints();
    const double dm_dh = xad::derivative(hin);
    auto names_o = field_names();
    for (std::size_t k = 0; k < idx.size(); ++k) {
      double dm_dtheta = xad::derivative(*fm[idx[k]]);
      double dF_dtheta = dm_dtheta - (names_o[idx[k]] == "omega" ? 1.0 : 0.0);
      dh0[k] = (dm_dh != 0.0) ? -dF_dtheta / dm_dh : 0.0;
    }
  }

  // ONE reverse sweep over the requested traits.
  ad::tape_type tape;
  auto pa = lift<ad_t>(pd);
  auto fp = field_ptrs<ad_t>(pa);
  for (auto i : idx) tape.registerInput(*fp[i]);
  tape.newRecording();
  // h0 active: value h0v + the IFT first-order injection so the tape carries
  // d(h0)/d(theta_k) for each registered trait (zero for traits not in mass_live).
  ad_t h0 = h0v;
  for (std::size_t k = 0; k < idx.size(); ++k)
    h0 = h0 + ad_t(dh0[k]) * (*fp[idx[k]] - ad_t(xad::value(*fp[idx[k]])));
  ad_t J = stand_offspring<ad_t>(pa, F, h0);
  tape.registerOutput(J); xad::derivative(J) = 1.0; tape.computeAdjoints();

  Rcpp::NumericVector grad(idx.size());
  for (std::size_t k=0;k<idx.size();++k) grad[k] = xad::derivative(*fp[idx[k]]);
  grad.attr("names") = Rcpp::wrap(traits);
  grad.attr("offspring_production") = as_double(J);
  return grad;
}
