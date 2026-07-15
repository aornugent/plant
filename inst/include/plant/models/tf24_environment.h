// Built from  inst/include/plant/models/ff16_environment.h on Mon Feb 12 09:52:27 2024 using the scaffolder, from the strategy:  FF16
// -*-c++-*-
#ifndef PLANT_PLANT_TF24_ENVIRONMENT_H_
#define PLANT_PLANT_TF24_ENVIRONMENT_H_

#include <plant/environment.h>
#include <plant/resource_spline.h>
#include <odelia/interpolator.hpp>
#include <odelia/ode_util.hpp> // odelia::util::to_passive (full AD strip, nesting-safe)
#include <limits>
#include <cmath>
#include <type_traits> // std::is_same_v

using namespace Rcpp;

namespace plant {

// TF24's environment: a light field plus a multi-layer soil-water state.
// Templated on the scalar S carried by the light knot values and the soil ODE
// state; S = double is the production path (the `TF24_Environment` alias below).
// The soil retention/conductivity curves and the water-balance rate arithmetic
// carry S so a trait derivative flows through the resident soil coupling; the
// fixed grid geometry (z / z_mid / dz) and the R-facing narrowings stay double.
// Base members are reached through this-> because Environment_<S> is a dependent
// base.
template <class S = double>
class TF24_Environment_ : public Environment_<S> {
public:
  // constructor for R interface - default settings can be modified
  // except for soil_number_of_depths
  // which are only updated on construction

  TF24_Environment_(bool light_availability_spline_rescale_usually = true,
                   int soil_number_of_depths = 5,
                   double delta_z = 9999, // not using this
                   double soil_moist_sat = 0.428, // saturated soil moisture content (m3 water m^-3 soil)
                   double K_sat = 163.0411, //saturated hydraulic conductivity of soil
                   double a_psi = 1.78e3, // not currently being used
                   double n_psi = 6.57, // not currently being used
                   double a_infil = 1, // infiltration switch (0-1), 0 no runoff, 1 runoff
                   double b_infil = 8, // unitless, determines infiltration rate
                   double depth = 1.5)  // total depth of soil (m)
      : delta_z(delta_z),
      soil_moist_sat(soil_moist_sat),
      a_psi(a_psi),
      n_psi(n_psi),
      K_sat(K_sat),
      a_infil(a_infil),
      b_infil(b_infil),
      depth(depth)
  {
    this->time = 0.0;



    // Shading defaults have lower tolerance which are overwritten for speed
    light_availability = ResourceSpline_<S>(
                   1e-4,  // light_availability_spline_tol,
                   17,    // light_availability_spline_nbase,
                   16,    // light_availability_spline_max_depth,
                   true //light_availability_spline_rescale_usually)
                  );

    this->extrinsic_drivers_set_constant("PPFD",1800);
    this->extrinsic_drivers_set_constant("rainfall",1);
    this->extrinsic_drivers_set_constant("atm_vpd",1);
    this->extrinsic_drivers_set_constant("ca",40);
    this->extrinsic_drivers_set_constant("leaf_temp",25);
    this->extrinsic_drivers_set_constant("atm_o2_kpa",21);
    this->extrinsic_drivers_set_constant("atm_kpa",100.5);

    set_soil_number_of_depths(soil_number_of_depths);
    set_soil_water_state(std::vector<double>(soil_number_of_depths, soil_moist_sat*0.5));
  };

  // Number of cumulative auxilliary variables to track in soil moisture model
  static constexpr size_t aux_num = 4;

  // Setup soil water distribtuion
  void set_soil_number_of_depths(int n) {
    soil_number_of_depths = n;

    this->vars = Internals_<S>(soil_number_of_depths + aux_num);

    z.resize(soil_number_of_depths);
    z_mid.resize(soil_number_of_depths);
    dz.resize(soil_number_of_depths);
    // positive downwards
    water_flux.resize(soil_number_of_depths);

    delta_z = depth / soil_number_of_depths;

    for (int i = 0; i < soil_number_of_depths; i++)
    {
      z[i] = (i + 1) * delta_z;
      if (i == 0) {
        z_mid[i] = z[i] / 2.0;
      } else {
        z_mid[i] = (z[i - 1] + z[i]) / 2.0;
      }
    }

    for (int i = 0; i < soil_number_of_depths; i++)
    {
      dz[i] = delta_z;
    }

    psi_soil_cache_.resize(soil_number_of_depths);
    psi_soil_cache_state_.resize(soil_number_of_depths);
    psi_soil_cache_valid_ = false;
  }
  int get_soil_number_of_depths() const {return soil_number_of_depths;}
  std::vector<double> get_soil_mid_depths() const { return z_mid; }

  // TODO: should we use auxilliary in internals
  // Per-layer gravitational drainage; carries S because it reads the active soil
  // state through soil_K_from_soil_theta.
  std::vector<S> water_flux;
  std::vector<double> z;
  std::vector<double> z_mid;
  std::vector<double> dz;
  // Cached soil water potential per layer (carries S). On the active path the
  // exact-compare memo is skipped (a value-only compare is blind to a changed
  // derivative), so it is recomputed each read; see get_soil_water_potential_state.
  mutable std::vector<S> psi_soil_cache_;
  mutable std::vector<S> psi_soil_cache_state_;
  mutable bool psi_soil_cache_valid_ = false;

  // Per-driver memo for the time-varying extrinsic drivers (see get_* above).
  // cache_time NaN-initialised so the first call always misses (NaN != time).
  double cached_driver_(const std::string &name, double &cache_val,
                        double &cache_time) const {
    if (cache_time != this->time) {
      cache_val = this->extrinsic_drivers.evaluate(name, this->time);
      cache_time = this->time;
    }
    return cache_val;
  }
  static constexpr double NAN_TIME_ = std::numeric_limits<double>::quiet_NaN();
  mutable double ppfd_cache_ = 0, atm_vpd_cache_ = 0, ca_cache_ = 0,
                 leaf_temp_cache_ = 0, atm_o2_kpa_cache_ = 0, atm_kpa_cache_ = 0;
  mutable double ppfd_cache_time_ = NAN_TIME_, atm_vpd_cache_time_ = NAN_TIME_,
                 ca_cache_time_ = NAN_TIME_, leaf_temp_cache_time_ = NAN_TIME_,
                 atm_o2_kpa_cache_time_ = NAN_TIME_, atm_kpa_cache_time_ = NAN_TIME_;

  // A ResourceSpline used for storing light availbility (0-1)
  ResourceSpline_<S> light_availability;

  // Light interface
  bool canopy_rescale_usually;
  //distance between layers
  int soil_number_of_depths;
  double delta_z;

  double depth;
  //saturated soil moisture
  double soil_moist_sat;
  //Saturated soil hydraulic conductivity
  double K_sat;
  double a_psi;
  double n_psi;
  double a_infil;
  double b_infil;

  // Residual soil moisture floor (m3 m^-3). The water balance cannot dry a
  // layer below this (see the positivity guard in compute_rates), and
  // psi_from_soil_moist is evaluated at this floor so the retention curve
  // (which diverges as theta->0) stays finite. Well below any realistic
  // operating moisture, so it does not perturb non-drought runs. See issue
  // #485.
  double soil_moist_residual = 1e-2;

  // Ability to prescribe a fixed value
  // TODO: add setting to set other variables like water
  void set_fixed_environment(double value, double height_max) {
    light_availability.set_fixed_value(value, height_max);
  }

  void set_fixed_environment(double value) {
    double height_max = 150.0;
    set_fixed_environment(value, height_max);
  }

  // Fix the light field to a scalar that may carry AD derivatives (unlike
  // set_fixed_environment(double), which strips them). Used to seed the resident
  // light as an active input when exercising the leaf seam's light channel.
  void set_fixed_environment_scalar(S value, double height_max) {
    light_availability.set_fixed_value_scalar(value, height_max);
  }

  S get_environment_at_height(S height) const {
    return light_availability.get_value_at_height(height);
  }

  // Highest height covered by the light spline; hoist out of hot per-point
  // loops and feed back into the capped get_environment_at_height() overload.
  double max_environment_height() const {
    return light_availability.max_height();
  }

  S get_environment_at_height(S height, S cap) const {
    return light_availability.get_value_at_height(height, cap);
  }

  virtual void r_init_interpolators(const std::vector<double> &state)
  {
    light_availability.r_init_interpolators(state);
  }

  // ------------------------------------------------------------------
  // SOIL WATER BALANCE (multi-layer bucket model)
  // ------------------------------------------------------------------
  // Each soil layer i is a bucket holding volumetric moisture vars.state(i)
  // (m3 water m^-3 soil). The rate of change is a simple mass balance:
  //
  //   d(theta_i)/dt = (water_in_i - water_out_i - root_uptake_i) / dz[i]
  //
  // where:
  //   * water_in_0   = infiltration (rainfall reduced by a saturation-excess
  //                    runoff term controlled by a_infil/b_infil);
  //   * water_in_i>0 = drainage out of the layer above (water_flux[i-1]);
  //   * water_out_i  = gravitational drainage = soil_K_from_soil_theta(theta_i),
  //                    a Clapp & Hornberger (1978) / Zeng & Decker (2009)
  //                    unsaturated hydraulic conductivity;
  //   * root_uptake_i= resource_depletion[i], supplied by the plants via the
  //                    strategy's evapotranspiration_dt (m yr^-1).
  //
  // The final `aux_num` state slots accumulate diagnostic cumulative fluxes
  // (rainfall, infiltration, deep drainage, total root uptake).
  //
  // This is an explicit, first-order representation; drainage is instantaneous
  // single-direction (no upward capillary flux between layers - that is handled
  // hydraulically inside the plant via E_from_Soil_to_Root_Collar).
  virtual void compute_rates(std::vector<S> const &resource_depletion)
  {

    S water_input;
    double rainfall = this->extrinsic_drivers.evaluate("rainfall", this->time);
    // saturation-excess runoff: 1 - a_infil*(theta_0/theta_sat)^b_infil, floored
    // at 0. The floor is a max against a scalar, written as a select so it also
    // holds at an active S.
    const S runoff_factor =
        1 - a_infil * pow(this->vars.state(0) / soil_moist_sat, b_infil);
    S infiltration = rainfall * ((runoff_factor > 0.0) ? runoff_factor : S(0.0));
    // Carries S: resource_depletion (plant water uptake) is active on the
    // resident soil coupling, so the cumulative-uptake diagnostic state is too.
    S total_resource_depletion = 0;


    // treat each soil layer as a separate resource pool
    for (size_t i = 0; i < soil_number_of_depths; i++)
    {

      // initial representation of drainage; to be improved
      if (i == 0)
      {
        water_input = infiltration;
      }
      else
      {
        // m3 m^-2
        water_input = water_flux[i-1];
      }
        // TODO: m3 m^-2
      water_flux[i] = soil_K_from_soil_theta(this->vars.state(i));
      // this function does runoff

      // Positivity guard (issue #485): a layer at or below the residual
      // moisture theta_r is not dried further (only rewetting is allowed). This
      // keeps the explicit fixed-step solver from driving a drought-stressed
      // layer to theta <= 0, where the retention curve psi_from_soil_moist and
      // the conductivity curve soil_K_from_soil_theta go non-finite. Wetter
      // layers are unaffected, so non-drought runs are unchanged.
      const S theta = this->vars.state(i);
      S rate = (water_input - water_flux[i] - resource_depletion[i]) / dz[i];
      if (theta <= soil_moist_residual && rate < 0.0) {
        rate = 0.0;
      }
      this->vars.set_rate(i, rate);
      total_resource_depletion += resource_depletion[i];
    }
      this->vars.set_rate(soil_number_of_depths, rainfall);
      this->vars.set_rate(soil_number_of_depths + 1, infiltration);
      this->vars.set_rate(soil_number_of_depths + 2, water_flux[soil_number_of_depths - 1]);
      this->vars.set_rate(soil_number_of_depths + 3, total_resource_depletion);

  }

  // calculate K from K_sat based on theta
  S soil_K_from_soil_theta(S theta) {
    //Eq. 5 Zeng and Decker (2009), ref Clapp and Hornberger (1978)
    // Floor at 0: an intermediate explicit-RK stage can probe theta < 0, and
    // std::pow(negative, non-integer) is NaN. A non-positive layer simply
    // drains nothing (K = 0). See issue #485.
    const S t = (theta > 0.0) ? theta : S(0.0);
    return K_sat * pow(t/soil_moist_sat, 2*n_psi + 3);
  }


  // convert soil moisture to soil water potential
  S psi_from_soil_moist(S soil_moist_) const {
    // Floor at the residual moisture: the retention curve (negative exponent)
    // diverges to +inf as theta->0, so an empty layer would otherwise yield a
    // non-finite potential. At/below theta_r the potential is large but finite
    // and the plant's root vulnerability curve has already shut uptake to ~0.
    const S t = (soil_moist_ > soil_moist_residual) ? soil_moist_
                                                    : S(soil_moist_residual);
    return a_psi * pow(t/soil_moist_sat, -n_psi)/1e6; // convert from Pa to MPa
  }

  // convert soil water potential to soil moisture
  S soil_moist_from_psi(S psi_soil_) const {
    return pow((psi_soil_/a_psi), (-1/n_psi))*soil_moist_sat;
  }

  // Easy wrappers. Cn also use `extrinsic_drivers_evaluate("PPFD", time)
  //
  // Each is read once per individual per ODE derivs evaluation (in the leaf
  // physiology setup), always at the current `time`. Memoise per driver keyed
  // on `time` so the repeated unordered_map<string,...> lookups across
  // individuals at the same time collapse to one lookup per distinct time.
  // The memo is per-driver (not eager-refresh-all) so an unset driver is only
  // looked up if it is actually requested — preserving the existing throw-if-
  // absent behaviour. The returned value is bit-identical to a direct evaluate.
  double get_PPFD()       const { return cached_driver_("PPFD", ppfd_cache_, ppfd_cache_time_); }
  double get_atm_vpd()    const { return cached_driver_("atm_vpd", atm_vpd_cache_, atm_vpd_cache_time_); }
  double get_ca()         const { return cached_driver_("ca", ca_cache_, ca_cache_time_); }
  double get_leaf_temp()  const { return cached_driver_("leaf_temp", leaf_temp_cache_, leaf_temp_cache_time_); }
  double get_atm_o2_kpa() const { return cached_driver_("atm_o2_kpa", atm_o2_kpa_cache_, atm_o2_kpa_cache_time_); }
  double get_atm_kpa()    const { return cached_driver_("atm_kpa", atm_kpa_cache_, atm_kpa_cache_time_); }


  std::vector<S> get_soil_water_state() const { return {this->vars.states.begin(), this->vars.states.end() - aux_num}; }
  const std::vector<S>& get_soil_water_potential_state() const {
    bool cache_stale = !psi_soil_cache_valid_ ||
      psi_soil_cache_state_.size() != static_cast<size_t>(soil_number_of_depths);

    // On the active path a value-only exact compare cannot see a changed
    // derivative, so never serve a cached psi: recompute at S each read (the
    // curve is closed-form and cheap). The memo is a double-run-only speedup.
    if constexpr (!std::is_same_v<S, double>) {
      cache_stale = true;
    }

    if (!cache_stale) {
      for (int i = 0; i < soil_number_of_depths; ++i) {
        if (psi_soil_cache_state_[i] != this->vars.state(i)) {
          cache_stale = true;
          break;
        }
      }
    }

    if (cache_stale) {
      psi_soil_cache_.resize(soil_number_of_depths);
      psi_soil_cache_state_.resize(soil_number_of_depths);
      for (int i = 0; i < soil_number_of_depths; ++i) {
        const S soil_moist = this->vars.state(i);
        psi_soil_cache_state_[i] = soil_moist;
        psi_soil_cache_[i] = psi_from_soil_moist(soil_moist);
      }
      psi_soil_cache_valid_ = true;
    }

    return psi_soil_cache_;
  }
  std::vector<S> get_soil_water_state_cumulative_flux() const { return {this->vars.states.end()-aux_num, this->vars.states.end()}; }
  std::vector<double> get_soil_depths() const { return z; }
  // double get_soil_depth(int layer) const { return z[layer]; }


  // TODO: I wonder if this needs a better name? See also environment.h
  Internals_<S> r_internals() const { return this->vars; }

  // R interface
  void set_soil_water_state(std::vector<double> state) {
    if(state.size() != (this->vars.state_size- aux_num)) {
      throw std::invalid_argument("Input vector size does not match soil state size.");
    }
    for (size_t i = 0; i < (this->vars.state_size); i++) {
      if(i < soil_number_of_depths){
        this->vars.set_state(i, state[i]);
      } else {
        this->vars.set_state(i, 0);
      }
  }
    psi_soil_cache_valid_ = false;
}

  // Pre-compute resources available in the environment, as a function of height
  template <typename Function>
  void compute_environment(Function f_compute_competition, double height_max, bool rescale) {

    // Define an anonymous function to use in creation of light_availability spline
    // Note: extinction coefficient was already applied in strategy, so
    // f_compute_competition gives sum of projected leaf area (k L) across species. Just need to apply Beer's law, E = exp(- (k L))
    auto f_light_availability = [&](double height) -> double
    { return exp(-f_compute_competition(height)); };

    // Calculates the light_availability spline, by fitting to the function
    // `f_compute_competition` as a function of height
    light_availability.compute_environment(f_light_availability, height_max, rescale);
  }

  virtual void clear_environment() {
    light_availability.clear();
  }

  virtual Rcpp::List r_get_state() const
  {

    // Surely an easier way?
    auto const &soil_depth_list = get_soil_depths();
    auto rcpp_soil_depth_vec = Rcpp::NumericVector(soil_depth_list.begin(), soil_depth_list.end());

    // Only doubles cross to R; strip every AD layer off the (possibly active)
    // soil state via to_passive (a no-op on the double path).
    auto const soil_moist_list = get_soil_water_state();
    Rcpp::NumericVector rcpp_soil_moist_vec(soil_moist_list.size());
    for (size_t i = 0; i < soil_moist_list.size(); ++i)
      rcpp_soil_moist_vec[i] = odelia::util::to_passive(soil_moist_list[i]);

    auto const soil_moist_cumulative_flux_list = get_soil_water_state_cumulative_flux();
    Rcpp::NumericVector rcpp_soil_moist_vec_cumulative_flux(soil_moist_cumulative_flux_list.size());
    for (size_t i = 0; i < soil_moist_cumulative_flux_list.size(); ++i)
      rcpp_soil_moist_vec_cumulative_flux[i] = odelia::util::to_passive(soil_moist_cumulative_flux_list[i]);

    return Rcpp::List::create(
        // auto ret = get_state(environment.extrinsic_drivers, time);

        _["light_availability"] = light_availability.r_get_state(),
        _["soil_moist"] = rcpp_soil_moist_vec,
        _["soil_depth"] = rcpp_soil_depth_vec,
        _["soil_moist_cumulative_flux"] = rcpp_soil_moist_vec_cumulative_flux
    );
  }
  };

// The double instantiation is the production path bound by RcppR6 as
// `plant::TF24_Environment`.
using TF24_Environment = TF24_Environment_<double>;
}

#endif
