// -*-c++-*-
#ifndef PLANT_PLANT_INTERNALS_MINIMAL_H_
#define PLANT_PLANT_INTERNALS_MINIMAL_H_

#define HEIGHT_INDEX 0
#define MORTALITY_INDEX 1
#define FECUNDITY_INDEX 2

#include <memory> // std::shared_ptr
#include <odelia/ode_interface.hpp>
#include <vector>
// #include <plant/plant_internals.h>

// TODO(#483): extra_state bounds, upper and lower limits
namespace plant {

class Internals {
public:
  Internals(size_t s_size=0, size_t a_size=0, size_t r_size=0)
      :
      state_size(s_size),
      aux_size(a_size),
      resource_size(r_size),
      states(s_size, 0.0),
      rates(s_size, NA_REAL) ,
      auxs(a_size, 0.0),
      consumption_rates(r_size, NA_REAL)
    {}
  size_t state_size;
  size_t aux_size;
  size_t resource_size;


  // Perhaps make these private so the () overloads below have some use
  std::vector<double> states;
  std::vector<double> rates;
  std::vector<double> auxs;
  std::vector<double> consumption_rates;  // not quite as pithy

  // T6 Slice 3b: per-individual uptake Jacobian d(consumption_rate[i])/d(theta_k),
  // row-major (i*n_soil + k), already in soil-moisture (ODE-state) space -- the
  // retention chain is folded in by the strategy fill, so all aggregation above
  // (Individual/Node/Species/Patch) stays env-agnostic. Empty unless the gated
  // fill runs (control.compute_uptake_jacobian); zero overhead and bit-identical
  // when off. Lazily sized by the strategy fill, not resize_consumption_rates.
  std::vector<double> duptake_jacobian;

  double state(int i) const { return states[i]; }
  double rate(int i) const { return rates[i]; }
  double aux(int i) const { return auxs[i]; }
  double consumption_rate(int i) const { return consumption_rates[i]; }
  // (i,k) into the row-major n_soil x n_soil Jacobian; ns = soil layer count.
  double duptake_jacobian_entry(int i, int k, int ns) const {
    return duptake_jacobian.empty() ? 0.0 : duptake_jacobian[i * ns + k];
  }

  void set_state(int i, double v) { states[i] = v; }
  void set_rate(int i, double v) { rates[i] = v; }
  void set_aux(int i, double v) { auxs[i] = v; }
  void set_consumption_rate(int i, double v) { consumption_rates[i] = v; }

  void resize(size_t new_size, size_t new_aux_size) {
    state_size = new_size;
    aux_size = new_aux_size;
    states.resize(new_size, 0.0);
    rates.resize(new_size, NA_REAL);
    auxs.resize(new_aux_size, 0.0);
  }

  void resize_consumption_rates(size_t new_resource_size) {
    resource_size = new_resource_size;
    consumption_rates.resize(new_resource_size, NA_REAL);
  }
};

} // namespace plant

#endif
