// -*-c++-*-
#ifndef PLANT_PLANT_CENSUS_H_
#define PLANT_PLANT_CENSUS_H_

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <string>
#include <utility>
#include <vector>

#include <plant/canopy_shape.h>
#include <plant/util.h>

// A stand census and the FIRST OF ITS TWO SENSITIVITY TERMS.
//
// A census summarises a stand by integrating a per-individual quantity over the
// size distribution:
//
//   C = sum_species sum_k w_k n_k m(h_k; phi)
//
// with w_k the trapezium weight on the height grid, n_k the stem number density
// of node k, and m the per-individual quantity the metric sums (leaf area, stem
// area, above-ground mass).
//
// Differentiating a census with respect to a strategy parameter phi gives two
// terms:
//
//   dC/dphi = sum_k w_k n_k dm/dphi        (1) DIRECT, at fixed state
//           + sum_k (d(w_k n_k)/dphi) m_k  (2) the response of the size
//                                              distribution itself
//
// Only (1) is implemented here. It is pure algebra on a frozen state: freeze the
// cohorts -- same heights, same densities, same weights -- and the census still
// moves with phi, because what each individual *is* at a given height is a
// function of the parameter. Term (2) requires the stand's development and is
// not computed anywhere in this package; a census gradient is not complete
// without it.
//
// The two terms are separable, and (1) is the one a trajectory/adjoint sweep of
// the size distribution can never produce, so it is worth having on its own and
// worth naming so nobody mistakes it for the whole derivative.

namespace plant {
namespace census {

// The census metrics, in output order. Each is a quadrature of a different
// per-individual quantity over the same size distribution.
inline std::vector<std::string> metric_names() {
  return {"area_leaf", "area_stem", "mass_above_ground"};
}

// The registered parameters: the union over metrics of the strategy parameters
// the census algebra actually reads. Derived from the allometry (see
// per_individual() below), not from the metrics' names:
//
//   area_leaf          A                       -> a_l1, a_l2
//   area_stem          theta (1+a_b1) A + Ah   -> a_l1, a_l2, theta, a_b1
//   mass_above_ground  lma A
//                      + eta_c rho theta (1+a_b1) h A
//                      + Mh                    -> all seven
//
// where A = area_leaf(h) = (h / a_l1)^(1 / a_l2) and eta_c = eta_c(eta).
//
// These are STRATEGY parameters, not traits. A trait sensitivity additionally
// needs the hyperpar Jacobian d(pars)/d(trait); that chain rule is not applied
// here.
inline std::vector<std::string> parameter_names() {
  return {"a_l1", "a_l2", "lma", "rho", "theta", "a_b1", "eta"};
}

inline std::size_t n_metrics() { return 3; }
inline std::size_t n_parameters() { return 7; }

// Indices into parameter_names().
enum ParIndex { PAR_A_L1 = 0, PAR_A_L2, PAR_LMA, PAR_RHO, PAR_THETA, PAR_A_B1,
                PAR_ETA };
enum MetricIndex { MET_AREA_LEAF = 0, MET_AREA_STEM, MET_MASS_ABOVE_GROUND };

// Which (metric, parameter) pairs the algebra above reaches. A pair that is
// false is a STRUCTURAL zero: no accumulator runs for it, and the reported
// derivative is exactly 0. A pair that is true had its accumulator run, so a
// zero there is a real zero (or a bug). Distinguishing the two is the point:
// without the flag an absent term and a vanishing term look identical, and a
// missing accumulator reads as a correct answer.
inline bool reaches(std::size_t metric, std::size_t par) {
  switch (metric) {
  case MET_AREA_LEAF:
    return par == PAR_A_L1 || par == PAR_A_L2;
  case MET_AREA_STEM:
    return par == PAR_A_L1 || par == PAR_A_L2 || par == PAR_THETA ||
           par == PAR_A_B1;
  case MET_MASS_ABOVE_GROUND:
    return true;
  default:
    return false;
  }
}

// d(eta_c)/d(eta) for eta_c = 1 - 2/(1+eta) + 1/(1+2 eta) (CanopyShape::eta_c).
inline double deta_c_deta(double eta) {
  const double a = 1.0 + eta, b = 1.0 + 2.0 * eta;
  return 2.0 / (a * a) - 2.0 / (b * b);
}

// The per-individual quantities and their parameter derivatives, at one height
// and heartwood state, for one strategy. `Pars` is any FF16-family pars struct
// (FF16, TF24, TF24f share this allometry exactly) and `eta_c` is the strategy's
// prepared crown constant.
//
// m[i]    = metric i's per-individual quantity
// dm[i*P+j] = d m[i] / d parameter j, zero (and never written) where
//             reaches(i, j) is false.
template <typename Pars>
void per_individual(const Pars& pars, double eta_c, double height,
                    double area_heartwood, double mass_heartwood,
                    double* m, double* dm) {
  const std::size_t P = n_parameters();

  // Leaf area [eqn 2] and its allometric derivatives.
  const double A = std::pow(height / pars.a_l1, 1.0 / pars.a_l2);
  const double dA_da_l1 = -A / (pars.a_l1 * pars.a_l2);
  // log(A) = log(height/a_l1)/a_l2, so d A/d a_l2 = -A log(A)/a_l2.
  const double dA_da_l2 = (A > 0.0) ? -A * std::log(A) / pars.a_l2 : 0.0;

  const double bark_plus_sap = 1.0 + pars.a_b1; // sapwood + bark, per theta*A
  // Mass of sapwood + bark per unit leaf area [eqn 4, 5].
  const double woody = eta_c * pars.rho * pars.theta * bark_plus_sap * height;

  // -- leaf area -----------------------------------------------------------
  m[MET_AREA_LEAF] = A;
  dm[MET_AREA_LEAF * P + PAR_A_L1] = dA_da_l1;
  dm[MET_AREA_LEAF * P + PAR_A_L2] = dA_da_l2;

  // -- stem (basal) area ---------------------------------------------------
  // Bark and sapwood are each an area per leaf area; heartwood area is state
  // and carries no parameter.
  m[MET_AREA_STEM] = pars.theta * bark_plus_sap * A + area_heartwood;
  dm[MET_AREA_STEM * P + PAR_A_L1]  = pars.theta * bark_plus_sap * dA_da_l1;
  dm[MET_AREA_STEM * P + PAR_A_L2]  = pars.theta * bark_plus_sap * dA_da_l2;
  dm[MET_AREA_STEM * P + PAR_THETA] = bark_plus_sap * A;
  dm[MET_AREA_STEM * P + PAR_A_B1]  = pars.theta * A;

  // -- above-ground mass ---------------------------------------------------
  // leaf + bark + sapwood + heartwood; heartwood mass is state.
  m[MET_MASS_ABOVE_GROUND] = (pars.lma + woody) * A + mass_heartwood;
  const std::size_t o = MET_MASS_ABOVE_GROUND * P;
  dm[o + PAR_A_L1]  = (pars.lma + woody) * dA_da_l1;
  dm[o + PAR_A_L2]  = (pars.lma + woody) * dA_da_l2;
  dm[o + PAR_LMA]   = A;
  dm[o + PAR_RHO]   = eta_c * pars.theta * bark_plus_sap * height * A;
  dm[o + PAR_THETA] = eta_c * pars.rho * bark_plus_sap * height * A;
  dm[o + PAR_A_B1]  = eta_c * pars.rho * pars.theta * height * A;
  dm[o + PAR_ETA]   = deta_c_deta(pars.eta) * pars.rho * pars.theta *
                      bark_plus_sap * height * A;
}

// Trapezium weights on the height grid, with the monotone guard.
//
// A census is a quadrature and it must be taken on a monotone grid. The node
// list is the grid, and reserve-gated growth (#517) lets a younger node overtake
// an older one, so the grid can invert (#571). On an inverted grid neighbouring
// trapezia have opposite-signed widths and CANCEL instead of accumulating -- the
// census is then wrong before any derivative of it exists. Sort by height, as
// Species::compute_competition_unordered and Species::consumption_rate already
// do for the other two field reductions; an already-monotone grid is untouched
// and its weights are bit-identical to the unsorted ones.
//
// Returns weights aligned with the input order, so the caller can multiply them
// straight onto per-node quantities. Every weight is non-negative.
inline std::vector<double> quadrature_weights(const std::vector<double>& height) {
  const std::size_t n = height.size();
  std::vector<double> w(n, 0.0);
  if (n < 2) {
    return w;
  }
  std::vector<std::size_t> ord(n);
  for (std::size_t i = 0; i < n; ++i) {
    ord[i] = i;
  }
  std::sort(ord.begin(), ord.end(),
            [&height](std::size_t a, std::size_t b) {
              return height[a] < height[b];
            });
  for (std::size_t j = 0; j + 1 < n; ++j) {
    const double dh = height[ord[j + 1]] - height[ord[j]];
    w[ord[j]]     += 0.5 * dh;
    w[ord[j + 1]] += 0.5 * dh;
  }
  return w;
}

// Is the grid already monotone in height (either direction)? Reported alongside
// the census so a caller can see when the guard above did something.
inline bool grid_is_monotone(const std::vector<double>& height) {
  return std::is_sorted(height.begin(), height.end()) ||
         std::is_sorted(height.rbegin(), height.rend());
}

// One species' census and its direct sensitivity term, on a frozen state.
//
// `height`, `density`, `area_heartwood`, `mass_heartwood` are the node grid,
// in any order; they are plain values, so nothing the strategy does can move
// them. `census` is filled with n_metrics() totals, `gradient` with
// n_metrics() x n_parameters() row-major entries. Both are accumulated into, so
// a caller can sum species by passing the same buffers.
template <typename Strategy>
void accumulate(const Strategy& s, const std::vector<double>& height,
                const std::vector<double>& density,
                const std::vector<double>& area_heartwood,
                const std::vector<double>& mass_heartwood,
                double* census, double* gradient) {
  const std::size_t n = height.size(), M = n_metrics(), P = n_parameters();
  util::check_length(density.size(), n);
  util::check_length(area_heartwood.size(), n);
  util::check_length(mass_heartwood.size(), n);

  const std::vector<double> w = quadrature_weights(height);
  std::vector<double> m(M), dm(M * P, 0.0);

  for (std::size_t k = 0; k < n; ++k) {
    per_individual(s.pars, s.eta_c, height[k], area_heartwood[k],
                   mass_heartwood[k], m.data(), dm.data());
    const double wn = w[k] * density[k];
    for (std::size_t i = 0; i < M; ++i) {
      census[i] += wn * m[i];
      for (std::size_t j = 0; j < P; ++j) {
        // Structural zeros are never accumulated, so they stay exactly 0 and
        // stay distinguishable from a term that ran and came out small.
        if (reaches(i, j)) {
          gradient[i * P + j] += wn * dm[i * P + j];
        }
      }
    }
  }
}

}
}

#endif
