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
#include <functional>
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

// Replay ONE cohort (index i) over the frozen schedule, returning its FINAL
// FF16LifeState<S> (the 5 demographic states + survival-weighted offspring
// accumulator at the end of the run, i.e. at the final census). Establishment is
// differentiated via the seedling net production -> mortality_0 = -log(pr_estab);
// h0 carries its own d/d(trait) via the IFT injection done by the caller. This is
// the single per-cohort replay that the generic reduction engine sums over: every
// emergent metric is a weighted reduction over these final cohort states, so the
// replay is recorded ONCE and reused for every metric (one tape, one reverse sweep
// per metric).
template <typename S>
plant::FF16LifeState<S> replay_cohort_final(const plant::FF16ProdPars<S>& pd,
                                            const Frozen& F, std::size_t i, S h0) {
  using std::exp; using std::log;
  const std::size_t b = (std::size_t)F.birth[i];
  const double ppsab = F.ppsab[i];
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
  return plant::ff16_cashkarp_replay(y, F.step_h, b, deriv, axpy);
}

// Emergent offspring_production = sum_i tw_i * offspring_weighted_i (deep-crown
// 6-state replay); establishment frozen via mort0. Kept as the dedicated routine
// behind offspring_production_gradient(); the generic engine below reproduces it
// as one symmetric (w, f) registry entry.
template <typename S>
S stand_offspring(const plant::FF16ProdPars<S>& pd, const Frozen& F, S h0) {
  S J = S(0.0);
  for (std::size_t i = 0; i < F.birth.size(); ++i)
    J += S(F.tw[i]) * replay_cohort_final<S>(pd, F, i, h0).offspring;
  return J;
}

// ---------------------------------------------------------------------------
// Generic (w, f) reduction engine (#472 scope B, build-order step 1). A stand
// metric is a weighted reduction over the replayed cohort final states,
//   metric = sum_i w_i * f(state_i),
// where w_i is a frozen per-cohort weight (from pass 1) and f maps a cohort's
// final FF16LifeState to a scalar contribution. offspring_production, LAI, biomass
// and size-distribution moments are all symmetric instantiations -- NONE privileged.
// The engine records ONE forward replay of every cohort onto a single tape, then
// takes one cheap reverse (adjoint) sweep PER metric, so M metrics cost
// replay + M*sweep, not M replays. This is the calibration-facing core: a
// metrics x traits Jacobian + the metric values, out of one resident baseline.
// ---------------------------------------------------------------------------

// A registered metric: its frozen per-cohort weight w_i and the contribution
// functor f(state) -> S. Census-time per-cohort number density (for LAI/biomass/
// size moments) is supplied to f as `density`, evolved as a replayed state.
template <typename S>
struct Metric {
  std::string name;
  std::vector<double> w;                                 // frozen per-cohort weight
  std::function<S(const plant::FF16LifeState<S>&)> f;    // contribution from state
};

// Assemble the Frozen harvest (resident schedule + per-RK-stage env + per-cohort
// survival/weights/decay) from the R-side arrays. Shared by every entry point.
Frozen build_frozen(const plant::FF16_Strategy& s, Rcpp::List eh_list,
                    const std::vector<double>& sh, std::vector<int> birth,
                    Rcpp::NumericMatrix ppsurv, std::vector<double> ppsab,
                    std::vector<double> tw) {
  Frozen F; F.eta = s.pars.eta; F.h0 = s.initial_height(); F.birth = birth;
  F.ppsab = ppsab; F.tw = tw; F.ppsurv = ppsurv; F.integ = &s.function_integrator;
  F.a_d0 = s.pars.a_d0;
  const std::size_t N = eh_list.size(); F.eh.resize(N); F.step_h.resize(N);
  for (std::size_t n=0;n<N;++n){Rcpp::List st=eh_list[n]; for(R_xlen_t k=0;k<st.size();++k) F.eh[n].push_back(Rcpp::as<plant::FF16_Environment>(st[k]));}
  for (std::size_t n=0;n<N;++n) F.step_h[n]=sh[n+1]-sh[n];
  // Per-cohort establishment decay exp(-recruitment_decay * birth_time).
  F.decay.resize(birth.size());
  for (std::size_t i=0;i<birth.size();++i)
    F.decay[i] = std::exp(-s.pars.recruitment_decay * sh[(std::size_t)birth[i]]);
  return F;
}

// Map requested trait names to FF16ProdPars field indices (errors on unknown).
std::vector<std::size_t> resolve_traits(const std::vector<std::string>& traits) {
  auto names = field_names();
  std::vector<std::size_t> idx;
  for (auto& t : traits) {
    auto it = std::find(names.begin(), names.end(), t);
    if (it == names.end()) Rcpp::stop("unknown FF16 trait: " + t);
    idx.push_back(std::distance(names.begin(), it));
  }
  return idx;
}

// d(height_0)/d(trait_k) by the implicit function theorem at the height_seed root
// (mass_live_given_height(h0) == omega), in a SCOPED reverse sweep of mass_live at
// h0 (the #539 seedling-size pattern). Returns the per-requested-trait dh0.
std::vector<double> compute_dh0(const plant::FF16ProdPars<double>& pd, double h0v,
                                const std::vector<std::size_t>& idx) {
  std::vector<double> dh0(idx.size(), 0.0);
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
  return dh0;
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
  Frozen F = build_frozen(s, eh_list, sh, birth, ppsurv, ppsab, tw);
  std::vector<std::size_t> idx = resolve_traits(traits);
  const double h0v = s.initial_height();
  std::vector<double> dh0 = compute_dh0(pd, h0v, idx);

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

// Compiled core of the generic stand-gradient engine (#472 scope B, build-order
// step 1). Records ONE forward replay of every cohort's final state onto a single
// adjoint tape, then for EACH requested metric m = sum_i w_i * f_m(state_i) takes
// one reverse sweep, giving the metrics x traits Jacobian (+ the metric values) out
// of one resident baseline. No metric is privileged: offspring_production is just
// one registered (w, f) entry. Census metrics (LAI/biomass/size moments, which need
// the per-cohort census number density) are added in a later increment; this first
// increment carries offspring_production to prove the engine reproduces the
// dedicated routine as one symmetric entry.
// [[Rcpp::export]]
Rcpp::List ff16_stand_gradient_impl(
    Rcpp::NumericVector pp, Rcpp::List eh_list, std::vector<double> sh,
    std::vector<int> birth, Rcpp::NumericMatrix ppsurv, std::vector<double> ppsab,
    std::vector<double> tw, std::vector<std::string> traits,
    std::vector<std::string> metrics) {
  auto s  = make_strategy(pp);
  auto pd = s.prod_pars();
  Frozen F = build_frozen(s, eh_list, sh, birth, ppsurv, ppsab, tw);
  std::vector<std::size_t> idx = resolve_traits(traits);
  const double h0v = s.initial_height();
  std::vector<double> dh0 = compute_dh0(pd, h0v, idx);
  const std::size_t M = metrics.size(), nC = F.birth.size(), nT = idx.size();

  // Build the requested metrics as (w, f) reductions. offspring_production:
  // w_i = tw_i (the node-spacing trapezoid * patch_density * S_D * birth_rate from
  // pass 1), f = the survival-weighted offspring accumulator at the final census.
  std::vector<Metric<ad_t>> mets; mets.reserve(M);
  for (auto& nm : metrics) {
    if (nm == "offspring_production") {
      mets.push_back({nm, F.tw,
        [](const plant::FF16LifeState<ad_t>& st) -> ad_t { return st.offspring; }});
    } else {
      Rcpp::stop("unknown stand metric: " + nm);
    }
  }

  // ONE forward recording: replay every cohort, then form each metric scalar.
  ad::tape_type tape;
  auto pa = lift<ad_t>(pd);
  auto fp = field_ptrs<ad_t>(pa);
  for (auto i : idx) tape.registerInput(*fp[i]);
  tape.newRecording();
  // h0 active via the IFT first-order injection (see compute_dh0).
  ad_t h0 = h0v;
  for (std::size_t k = 0; k < nT; ++k)
    h0 = h0 + ad_t(dh0[k]) * (*fp[idx[k]] - ad_t(xad::value(*fp[idx[k]])));

  std::vector<plant::FF16LifeState<ad_t>> finals; finals.reserve(nC);
  for (std::size_t i = 0; i < nC; ++i)
    finals.push_back(replay_cohort_final<ad_t>(pa, F, i, h0));

  std::vector<ad_t> J(M);
  for (std::size_t m = 0; m < M; ++m) {
    ad_t acc = ad_t(0.0);
    for (std::size_t i = 0; i < nC; ++i) acc = acc + ad_t(mets[m].w[i]) * mets[m].f(finals[i]);
    J[m] = acc;
    tape.registerOutput(J[m]);
  }

  // One cheap reverse sweep PER metric over the single recording (clear adjoints
  // between, the XAD multi-output Jacobian pattern).
  Rcpp::NumericMatrix jac(M, nT);
  Rcpp::NumericVector values(M);
  for (std::size_t m = 0; m < M; ++m) {
    tape.clearDerivatives();
    xad::derivative(J[m]) = 1.0;
    tape.computeAdjoints();
    for (std::size_t k = 0; k < nT; ++k) jac(m, k) = xad::derivative(*fp[idx[k]]);
    values[m] = as_double(J[m]);
  }

  jac.attr("dimnames") = Rcpp::List::create(Rcpp::wrap(metrics), Rcpp::wrap(traits));
  values.attr("names") = Rcpp::wrap(metrics);
  return Rcpp::List::create(Rcpp::Named("jacobian") = jac,
                            Rcpp::Named("values")   = values);
}

// Escape hatch (#472 scope B, build-order step 1): the per-cohort state x trait
// Jacobian. For ANY emergent metric that is NOT a simple reduction -- quantiles,
// ratios, bespoke statistics a downstream package invents -- expose d(state_i,c)/
// d(theta_k) for every cohort i, every final-state component c, and let the metric
// gradient compose downstream by chain rule (the same boundary as "likelihoods
// live downstream"). Each cohort's final state is INDEPENDENT (it depends only on
// its own replay), so we tape ONE cohort at a time and take one reverse sweep per
// state component -- nC small sweeps, not one giant tape. Returns the cohort final
// states + the [cohort, component, trait] Jacobian.
// [[Rcpp::export]]
Rcpp::List ff16_state_jacobian_impl(
    Rcpp::NumericVector pp, Rcpp::List eh_list, std::vector<double> sh,
    std::vector<int> birth, Rcpp::NumericMatrix ppsurv, std::vector<double> ppsab,
    std::vector<double> tw, std::vector<std::string> traits) {
  auto s  = make_strategy(pp);
  auto pd = s.prod_pars();
  Frozen F = build_frozen(s, eh_list, sh, birth, ppsurv, ppsab, tw);
  std::vector<std::size_t> idx = resolve_traits(traits);
  const double h0v = s.initial_height();
  std::vector<double> dh0 = compute_dh0(pd, h0v, idx);
  const std::size_t nC = F.birth.size(), nT = idx.size();

  const std::vector<std::string> comp =
    {"height","mortality","fecundity","area_heartwood","mass_heartwood","offspring"};
  const std::size_t nS = comp.size();

  Rcpp::NumericMatrix states(nC, nS);
  // 3D Jacobian, R column-major dim [nC, nS, nT].
  Rcpp::NumericVector jac(nC * nS * nT);
  auto JAC = [&](std::size_t i, std::size_t c, std::size_t k) -> double& {
    return jac[i + nC * (c + nS * k)];
  };

  for (std::size_t i = 0; i < nC; ++i) {
    ad::tape_type tape;
    auto pa = lift<ad_t>(pd);
    auto fp = field_ptrs<ad_t>(pa);
    for (auto j : idx) tape.registerInput(*fp[j]);
    tape.newRecording();
    ad_t h0 = h0v;
    for (std::size_t k = 0; k < nT; ++k)
      h0 = h0 + ad_t(dh0[k]) * (*fp[idx[k]] - ad_t(xad::value(*fp[idx[k]])));

    plant::FF16LifeState<ad_t> y = replay_cohort_final<ad_t>(pa, F, i, h0);
    ad_t out[6] = {y.demog.height, y.demog.mortality, y.demog.fecundity,
                   y.demog.area_heartwood, y.demog.mass_heartwood, y.offspring};
    for (std::size_t c = 0; c < nS; ++c) {
      states(i, c) = as_double(out[c]);
      tape.registerOutput(out[c]);
    }
    for (std::size_t c = 0; c < nS; ++c) {
      tape.clearDerivatives();
      xad::derivative(out[c]) = 1.0;
      tape.computeAdjoints();
      for (std::size_t k = 0; k < nT; ++k) JAC(i, c, k) = xad::derivative(*fp[idx[k]]);
    }
  }

  states.attr("dimnames") = Rcpp::List::create(R_NilValue, Rcpp::wrap(comp));
  jac.attr("dim") = Rcpp::IntegerVector::create((int)nC, (int)nS, (int)nT);
  jac.attr("dimnames") = Rcpp::List::create(R_NilValue, Rcpp::wrap(comp),
                                            Rcpp::wrap(traits));
  return Rcpp::List::create(Rcpp::Named("states") = states,
                            Rcpp::Named("jacobian") = jac);
}
