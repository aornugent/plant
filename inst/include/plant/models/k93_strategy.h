// Built from  inst/include/plant/models/ff16_strategy.h on Fri Jul 24 10:23:19 2020 using the scaffolder, from the strategy:  FF16
// -*-c++-*-
#ifndef PLANT_PLANT_K93_STRATEGY_H_
#define PLANT_PLANT_K93_STRATEGY_H_

#include <plant/strategy.h>
#include <plant/models/k93_environment.h>
#include <plant/canopy_shape.h>

namespace plant {

// Biological (user-settable) parameters for the K93 strategy. Held as a value
// member `pars` on K93_Strategy and exposed to R as a nested RcppR6 list class
// (access as `s$pars$b_0`). Templated on the scalar S so a trait derivative
// flows when a field is seeded; S = double is the production path (the `K93_Pars`
// alias below). Derived quantities (canopy_shape) stay as plain members on the
// strategy.
template <class S = double>
struct K93_Pars_ {
  // Initial seedling size (dbh cm)
  S height_0 = 2.0;
  // * Growth
  S b_0 = 0.059;    // Growth intercept (yr-1)
  S b_1 = 0.012;    // Growth asymptote (yr-1.(ln cm)-1)
  S b_2 = 0.00041;  // Growth suppression rate (m2.cm-2.yr-1)
  // * Mortality
  S c_0 = 0.008;    // Intercept (yr-1)
  S c_1 = 0.00044;  // Suppression rate (m2.cm-2.yr-1)
  // * Reproduction
  S d_0 = 0.00073;  // Recruitment rate (cm2.yr-1)
  S d_1 = 0.044;    // Reduction from suppression (m2.cm-2.yr-1)
  // Probability of survival during dispersal (required by scm.h)
  S S_D = 1.0;
  // Smoothing / canopy shape parameter
  S eta = 12;
  // Light capture parameter
  S k_I = 0.01;
};

using K93_Pars = K93_Pars_<double>;

// The K93 (Kohyama 1993) strategy. Templated on the scalar S carried by its
// physiology; S = double is the production path (the `K93_Strategy` alias
// below). Closed-form rates -- no iteration, no quadrature -- so this is the
// smallest strategy to carry an active scalar. Base members of Strategy<E> are
// reached through this-> because the base is dependent on S.
template <class S = double>
class K93_Strategy_ : public Strategy<K93_Environment_<S>> {
public:
  using environment_type = K93_Environment_<S>;
  using value_type = S;
  typedef std::shared_ptr<K93_Strategy_<S>> ptr;

  K93_Strategy_() {
    // Empirical parameter defaults (Table 1) live in the K93_Pars struct
    // member initialisers.

    // build the string state/aux name to index map
    refresh_indices();
    this->name = "K93";
  }

  // Direct aux indices for the hot path, avoiding aux_index.at("...") string-map
  // lookups (these showed up in profiling; see #466). MUST stay in sync with the
  // order of aux_names() below. refresh_indices() still fills the named maps used
  // by the R-facing paths.
  static constexpr int COMPETITION_EFFECT_AUX_INDEX = 0;
  static constexpr int HEIGHT_INVERSE_AUX_INDEX = 1;

  // update this when the length of state_names changes
  static size_t state_size () { return 3; }
  // update this when the length of aux_names changes
  size_t aux_size () { return aux_names().size(); }

  static std::vector<std::string> state_names() {
    return  std::vector<std::string>({
      "height",
      "mortality",
      "fecundity"
      });
  }

  std::vector<std::string> aux_names() {
    return std::vector<std::string>({"competition_effect", "height_inverse"});
  }

  // Signatures fixed in plant.h
  void update_dependent_aux(const int index, Internals_<S>& vars) {
    if (index == HEIGHT_INDEX) {
      S height = vars.state(HEIGHT_INDEX);
      vars.set_aux(COMPETITION_EFFECT_AUX_INDEX,
                   compute_competition_by_ratio(
                     0.0, pars.k_I * size_to_basal_area(height)));
      vars.set_aux(HEIGHT_INVERSE_AUX_INDEX, 1.0 / height);
    }
  }

  void refresh_indices() {
    // Create and fill the name to state index maps
    this->state_index = std::map<std::string,int>();
    this->aux_index   = std::map<std::string,int>();
    std::vector<std::string> aux_names_vec = aux_names();
    std::vector<std::string> state_names_vec = state_names();
    for (size_t i = 0; i < state_names_vec.size(); i++) {
      this->state_index[state_names_vec[i]] = i;
    }
    for (size_t i = 0; i < aux_names_vec.size(); i++) {
      this->aux_index[aux_names_vec[i]] = i;
    }
  }

  double establishment_probability(const environment_type& environment) {
    (void) environment;
    // K93 (Kohyama 1993) has no carbon-budget establishment filter:
    // establishment is deterministic, so every dispersed seed establishes.
    // Demographic filtering instead happens via growth suppression (size_dt is
    // clamped to >= 0 under competition) and the g > 0 density guard in
    // Node::compute_initial_conditions. This value is consumed by the SCM
    // (initial mortality / density) and the stochastic germination test.
    // TODO(#480): may want to make this dependent on achieving positive growth rate.
    return 1.0;
  }

  double net_mass_production_dt(const environment_type& environment,
                                S height, S area_leaf_) {
    (void) environment;
    (void) height;
    (void) area_leaf_;
    // K93 models growth directly (size_dt); it has no carbon mass-production
    // budget, so net mass production is undefined for this strategy. The value
    // is never consumed by compute_rates -- it is reachable only via the R
    // diagnostic accessor and resource_compensation_point(). Return NA to make
    // "not applicable" explicit rather than returning a misleading number.
    return NA_REAL;
  }

  double net_mass_production_dt(const environment_type& environment,
                                S height, S area_leaf_,
                                S height_inverse) {
    (void) height_inverse;
    return net_mass_production_dt(environment, height, area_leaf_);
  }
  // Strategy-agnostic entry point used by Individual<K93> (#266). K93 has no
  // carbon budget, so the worker ignores these arguments and returns NA; the
  // wrapper exists to keep Individual's interface uniform across strategies.
  double net_mass_production_dt(const environment_type& environment,
                                const Internals_<S>& vars) {
    return net_mass_production_dt(environment, vars.state(HEIGHT_INDEX),
                                  vars.aux(COMPETITION_EFFECT_AUX_INDEX),
                                  vars.aux(HEIGHT_INVERSE_AUX_INDEX));
  }

  // i.e. setting rates of ode vars from the state and updating aux vars
  void compute_rates(const environment_type& environment, Internals_<S>& vars) {

    S height = vars.state(HEIGHT_INDEX);

    // suppression integral mapped [0, 1] using adaptive spline
    // back transform to basal area and add suppression from self
    S competition = environment.get_environment_at_height(height);

    S cumulative_basal_area = -log(competition) / pars.k_I;

    if (!util::is_finite(cumulative_basal_area))
    {
      util::stop("Environmental interpolator has gone out of bounds, try lowering the extinction coefficient k_I");
    }

    vars.set_rate(HEIGHT_INDEX,
      size_dt(height, cumulative_basal_area));

    vars.set_rate(FECUNDITY_INDEX,
      fecundity_dt(height, cumulative_basal_area));

    vars.set_rate(MORTALITY_INDEX,
      mortality_dt(cumulative_basal_area, vars.state(MORTALITY_INDEX)));
  }

  S compute_competition(double z, S size) const {
    return compute_competition(z, pars.k_I * size_to_basal_area(size), 1.0 / size);
  }
  // Hot path (called per node from Species::compute_competition): defined inline
  // here so the no-LTO build folds it into the loop instead of paying a cross-TU
  // call (mirrors FF16_Strategy; see the profile-plant skill).
  S compute_competition(double z, S whole_plant_competition,
                        S height_inverse) const {
    return compute_competition_by_ratio(z * height_inverse, whole_plant_competition);
  }
  // Templated on the ratio scalar R: in update_dependent_aux the ratio is a
  // double literal (0.0), so canopy_shape.Q stays a double read there; only the
  // resident field-build hot path (deferred) drives an active ratio, which will
  // need a templated CanopyShape.
  template <typename R>
  S compute_competition_by_ratio(R z_over_size,
                                 S whole_plant_competition) const {
    // Competition only felt if plant bigger than target size z.
    return whole_plant_competition * canopy_shape.Q(z_over_size);
  }
  // Strategy-agnostic entry point used by Individual<K93> (#266): reads the
  // cached competition_effect and height_inverse aux slots itself.
  S compute_competition(double z, const Internals_<S>& vars) const {
    return compute_competition(z, vars.aux(COMPETITION_EFFECT_AUX_INDEX),
                               vars.aux(HEIGHT_INVERSE_AUX_INDEX));
  }


  // K93 Methods  ----------------------------------------------
  S size_to_basal_area(S size) const {
    return M_PI / 4 * pow(size, 2);
  }

  // [eqn 10] Growth
  S size_dt(S size, S cumulative_basal_area) const {

    S growth = size * (pars.b_0 - pars.b_1 * log(size) - pars.b_2 * cumulative_basal_area);

    if(growth < 0.0) {
      growth = 0.0;
    }

    return growth;
  }

  // [eqn 12] Reproduction
  S fecundity_dt(S size, S cumulative_basal_area) const {
    S basal_area = size_to_basal_area(size);
    return pars.d_0 * basal_area * exp(-pars.d_1 * cumulative_basal_area);
  }

  // [eqn 11] Mortality
  S mortality_dt(S cumulative_basal_area,
                 S cumulative_mortality) const {
    // If mortality probability is 1 (latency = Inf) then the rate
    // calculations break.  Setting them to zero gives the correct
    // behaviour.
    if (util::is_finite(cumulative_mortality)) {
      S mu = -pars.c_0 + pars.c_1 * cumulative_basal_area;
      return (mu > 0)? mu : static_cast<S>(0.0);
   } else {
      return 0.0;
    }
  }

  // useful for pre-computing expensive objects
  void prepare_strategy() {
    canopy_shape.initialise(xad::value(pars.eta));

    if (this->is_variable_birth_rate) {
      this->extrinsic_drivers.set_variable("birth_rate", this->birth_rate_x, this->birth_rate_y);
    } else {
      this->extrinsic_drivers.set_constant("birth_rate", this->birth_rate_y[0]);
    }
  }

  // Birth height (dbh, cm) of a seedling. Strategy-agnostic accessor used by
  // the templated Individual; for K93 this is the user-set pars.height_0.
  S initial_height() const { return pars.height_0; }

  // Biological (user-settable) parameters; see K93_Pars above.
  K93_Pars_<S> pars;

  // Derived / precomputed in prepare_strategy() (NOT user-set)
  CanopyShape canopy_shape;
};

using K93_Strategy = K93_Strategy_<double>;

}

#endif
