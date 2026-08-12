// -*-c++-*-
#ifndef SPECIES
#define SPECIES

#include <vector>
#include <algorithm>
#include <limits>
#include <tuple>
#include <utility>
#include <plant/util.h>
#include <plant/environment.h>
#include <odelia/ode_interface.hpp>
#include <plant/node.h>
#include <plant/species_base.h>
#include <plant/transport_census.h>
#include <odelia/drivers.hpp>

namespace plant {

// One census metric maps one cohort to the quantity summed over the size
// distribution, at the scalar the strategy carries. A metric reads the strategy
// and one cohort and nothing else, so a metric with a height cut of its own
// carries that cut. Adding a metric is a struct of this shape plus its name in
// the tuple below.
namespace census_metric {

struct leaf_area {
  template <class Strategy, class Individual>
  typename Individual::value_type
  operator()(const Strategy& strategy, const Individual& individual) const {
    return strategy.area_leaf(individual.state(HEIGHT_INDEX));
  }
  static const char* name() { return "leaf_area"; }
};

struct mass_above_ground {
  template <class Strategy, class Individual>
  typename Individual::value_type
  operator()(const Strategy& strategy, const Individual& individual) const {
    using S = typename Individual::value_type;
    const S height = individual.state(HEIGHT_INDEX);
    const S area_leaf = strategy.area_leaf(height);
    return strategy.mass_above_ground(
        strategy.mass_leaf(area_leaf),
        strategy.mass_bark(strategy.area_bark(area_leaf), height),
        strategy.mass_sapwood(strategy.area_sapwood(area_leaf), height),
        individual.state("mass_heartwood"));
  }
  static const char* name() { return "mass_above_ground"; }
};

struct area_stem {
  template <class Strategy, class Individual>
  typename Individual::value_type
  operator()(const Strategy& strategy, const Individual& individual) const {
    using S = typename Individual::value_type;
    const S height = individual.state(HEIGHT_INDEX);
    const S area_leaf = strategy.area_leaf(height);
    return strategy.area_stem(strategy.area_bark(area_leaf),
                              strategy.area_sapwood(area_leaf),
                              individual.state("area_heartwood"));
  }
  static const char* name() { return "area_stem"; }
};

}

// The metrics TF24 is censused on. A census's codomain is the tuple's size, so a
// fourth metric is one name added here.
using tf24_census = std::tuple<census_metric::leaf_area,
                               census_metric::mass_above_ground,
                               census_metric::area_stem>;

// This is purely for running the deterministic model. It shares its storage and
// ODE plumbing with the stochastic species through SpeciesBase (species_base.h);
// the size-density-specific machinery below (density-weighted competition,
// survival-weighted rates, lifetime fitness, schedule-refinement error) stays
// here.

template <typename T, typename E>
class Species : public SpeciesBase<Species<T, E>, T, E, Node<T, E>> {
  typedef SpeciesBase<Species<T, E>, T, E, Node<T, E>> base_type;
public:
  using value_type = typename T::value_type;

  typedef T         strategy_type;
  typedef E         environment_type;
  typedef Individual<T,E>  individual_type;
  typedef Node<T,E> node_type;
  typedef typename strategy_type::ptr strategy_type_ptr;
  Species(strategy_type s);
  // Build on a strategy that is already prepared; see SpeciesBase.
  explicit Species(strategy_type_ptr s);

  // ODE plumbing and the per-element serialisers are inherited from SpeciesBase
  // and iterate all nodes (the deterministic model has no notion of "dead").
  using base_type::ode_size;
  using base_type::ode_state;
  using base_type::ode_rates;
  using base_type::get_node_state;
  using base_type::get_node_aux;
  typename std::vector<node_type>::iterator node_begin() { return nodes.begin(); }
  typename std::vector<node_type>::iterator node_end() { return nodes.end(); }
  typename std::vector<node_type>::const_iterator node_begin() const { return nodes.begin(); }
  typename std::vector<node_type>::const_iterator node_end() const { return nodes.end(); }

  size_t size() const;
  void clear();
  void introduce_new_node();
  // Introduce a node, stamping it with the introduction time and patch-age
  // density at birth (called by Patch, which knows the time and disturbance).
  void introduce_new_node(double time, double patch_density);
  // Drop the node introduce_new_node pushed last, which is the newest: the width
  // a reverse sweep needs before an introduction.
  void remove_newest_node();

  value_type height_max() const;
  // The query height is a knot position on the interpolant's own grid, which is
  // double (see ResourceSpline::rebuild_spline); what the cohorts put into the
  // field arrives through the node contributions summed over below.
  value_type compute_competition(double height) const;

  // The reduction and its vertical derivative from one traversal, so each node's
  // u^eta is evaluated once and the two sums add their terms in the same order.
  // The first entry equals compute_competition(height) bit for bit.
  std::pair<value_type, value_type>
  compute_competition_and_slope(double height) const;

  // The same reduction with the inflow boundary interval left off, so it is a
  // function of the ODE state alone: it never reads new_node. The boundary
  // condition n_b = birth_rate * pr_estab / g needs a field to be evaluated in,
  // and evaluating it in this one is what breaks the cycle. See
  // Patch::compute_environment.
  value_type compute_competition_excl_boundary(double height) const;

  // compute_competition_and_slope() with that same interval left off.
  std::pair<value_type, value_type>
  compute_competition_and_slope_excl_boundary(double height) const;

  // Transpose of compute_competition_and_slope at `height`, closing boundary
  // trapezium and all. `out` points at this species' first node.
  void compute_competition_and_slope_adjoint(double height,
                                             double lambda_value,
                                             double lambda_slope,
                                             node_size_adjoints* out) const;

  // Transpose of consumption_rate for one resource. The grid is the quadrature
  // abscissa, so the weights reach the heights only where that abscissa is one.
  void consumption_rate_adjoint(int resource, double lambda_uptake,
                                node_uptake_adjoints* out) const;

  // Evaluate the inflow boundary condition in the environment passed. Split out
  // of compute_rates() so the field build owns it and the field stops reading a
  // density carried from the previous evaluation.
  void compute_boundary_node(const environment_type& environment,
                             double pr_patch_survival, double birth_rate) {
    new_node.compute_initial_conditions(environment, pr_patch_survival, birth_rate);
  }

  // The trapezium integral of n_k psi(state_k) over the size distribution, with
  // n_k = exp(l_k). A census is a quadrature of a density, so the grid is the
  // coordinate that density is carried in; taking gaps in any other variable
  // integrates one density against another's spacing. The inflow boundary node
  // closes the grid, being the lower limit of the distribution.
  template <class Psi> value_type census(Psi psi) const;

  // Whether the decreasing-height node ordering still holds (see height_max()).
  bool heights_are_decreasing() const;

  // The tallest node height and whether the heights are ordered, from a single
  // pass. compute_competition() needs both on every call, and walking the heights
  // twice was measurably slower on FF16 (~5% on the SCM benchmark) than the
  // O(1) nodes.front() it replaced.
  // h_max is the canopy top, a position; decreasing is a comparison outcome and
  // structural, so it stays bool whatever the heights are made of.
  struct HeightScan { value_type h_max; bool decreasing; };
  // Cached: heights change only when the ODE state is set or a node is
  // introduced/cleared, whereas compute_competition() is called once per spline
  // knot, so this is hundreds of calls per change. Every mutator invalidates.
  HeightScan scan_heights() const;

  // Setting the ODE state rewrites every node's height, so the cached scan goes
  // with it. Shadows (rather than uses) the SpeciesBase version for that reason.
  template <typename It> It set_ode_state(It it) {
    invalidate_height_scan();
    return base_type::set_ode_state(it);
  }
  void compute_rates(const environment_type& environment, double pr_patch_survival, double birth_rate);

  // -dg/dh on the cohort grid: node i spans the interval down to its lower
  // neighbour, the lowest down to new_node. The spacing between two
  // characteristics has an exact rate, d(dh)/dt = g_i - g_below, so this is
  // exactly d(log dh)/dt rather than an estimate of dg/dh. Where the interval has
  // no width the node takes the compression of the one above.
  value_type growth_rate_gradient(std::size_t i) const;

  std::vector<double> net_reproduction_ratio_by_node() const;
  // Per-node lifetime offspring, weighted by patch-age density and S_D.
  std::vector<double> net_reproduction_ratio_by_node_weighted() const;
  // Introduction times of each node (the integration x-axis for fitness).
  std::vector<double> node_times() const;

  // The boundary node's birth date is *now*, but compute_initial_conditions()
  // (which stamps it) runs inside compute_rates(), which the ODE stepper calls
  // after set_ode_state() has already rebuilt the environment. Reading the stamp
  // during that rebuild would therefore pick up the previous derivs call's time
  // and shorten the boundary trapezium by one Runge-Kutta stage. Patch refreshes
  // it before building the profile. Harmless on the height path, where the
  // boundary abscissa is the constant initial height.
  void set_new_node_birth_date(double time) {
    new_node.set_introduction_time(time);
  }

  // Two nodes sharing a birth date give a zero-width trapezium interval, so the
  // birth-date quadrature silently loses them. Cannot happen for a scheduled
  // run (introduction times are distinct by construction) but can for a patch
  // seeded or imported without per-node times.
  bool birth_dates_are_distinct() const;

  // Which coordinate this species' size distribution is carried in. Exposed so
  // Patch can check every species agrees before summing their contributions.
  bool density_in_birth_date() const {
    return control().node_density_in_birth_date;
  }

  // * ODE interface
  // NOTE: We are a time-independent model here so no need to pass
  // time in as an argument.  All the bits involving time are taken
  // care of by Environment for us.
  // (ode_size/set_ode_state/ode_state/ode_rates come from SpeciesBase.)
  size_t aux_size() const;

  void resize_consumption_rates(int i);
  value_type consumption_rate(int i) const;
  std::vector<value_type> consumption_rate_by_node_rev(int i) const;
  std::vector<value_type> consumption_rate_by_node(int i) const;

  template <typename It> It ode_aux(It it) const;

  Rcpp::NumericMatrix r_get_state() const;

  // * R interface
  std::vector<double> r_heights() const;
  std::vector<double> r_heights_rev() const;
  void r_set_heights(std::vector<double> heights);
  const node_type& node_at(size_t i) const {return nodes[i];}
  const node_type& r_new_node() const {return new_node;}
  std::vector<node_type> r_nodes() const {return nodes;}
  const node_type& r_node_at(util::index idx) const {
    return nodes[idx.check_bounds(size())];
  }

  // Do this with set_ode_state, using an iterator?
  /* double state(int i) const { return vars.state(i); } */

  /* double rate(int i) const { return vars.rate(i); } */

  /* void set_state(int i, double v) { */
  /*   vars.set_state(i, v); */
  /* } */


  // These are used to determine the degree of node refinement.
  std::vector<double> r_compute_competition_effect_by_nodes() const;
  std::vector<double> r_compute_competition_effect_by_nodes_error(double scal) const;

  // Per-node size density, **always** as a density in height whichever
  // coordinate the solver carried it in, so downstream code (tidy_outputs.R's
  // `density`, interpolate_to_heights(), the plots) keeps its meaning. In
  // birth-date coordinates that means dividing the carried quantity by the
  // Jacobian; see height_jacobian(). NA for a node where the Jacobian vanishes,
  // which is where the height density genuinely does not exist.
  std::vector<double> r_log_densities() const;
  // The quantity actually integrated: the density in height on the height path,
  // and the density in birth date (nu) on the birth-date one. Reported
  // alongside r_log_densities() rather than instead of it, because it is the
  // thing whose ODE the solver solves and the thing that cannot go non-finite.
  std::vector<double> r_log_densities_state() const;
  // Per-node rate of change of log density; used to guard against initial
  // conditions whose densities would explode to non-finite values. This is the
  // rate of the *carried* quantity, so on the birth-date path it is -mortality.
  std::vector<double> r_log_density_rates() const;

  // |dh/dtau| per node, with the boundary node appended last: the Jacobian of
  // the change of variables between the two coordinates, N = nu / |dh/dtau|.
  std::vector<double> height_jacobian() const;

  // Per-node birth bookkeeping, exposed so an exported patch state can be
  // re-imported faithfully (see node.h::set_birth_state). node_times() above
  // already returns the per-node introduction times.
  std::vector<double> r_patch_densities() const;
  std::vector<double> r_pr_patch_survival_at_birth() const;
  // Restore birth bookkeeping for imported nodes (resume); the argument lengths
  // must each match the current node count.
  void set_birth_state(const std::vector<double>& times,
                       const std::vector<double>& patch_density,
                       const std::vector<double>& pr_patch_survival);

  ExtrinsicDrivers extrinsic_drivers() const {return strategy->extrinsic_drivers;}

  // The prepared strategy this species and its nodes share.
  strategy_type_ptr strategy_ptr() const {return this->strategy;}

private:
  // compute_competition() for the case where the node heights are no longer
  // ordered, so the node list cannot be used directly as the quadrature grid.
  // Height coordinate only -- it integrates in height, and the birth-date
  // abscissa cannot invert (see compute_competition).
  value_type compute_competition_unordered(double height,
                                           bool include_boundary) const;

  // compute_competition_and_slope() over a height-sorted view, for the same
  // broken-ordering case compute_competition_unordered() handles.
  std::pair<value_type, value_type>
  compute_competition_and_slope_unordered(double height,
                                          bool include_boundary) const;

  // The fused reduction, with the closing boundary trapezium included or not.
  std::pair<value_type, value_type>
  compute_competition_and_slope_impl(double height, bool include_boundary) const;

  // The reduction, with the closing boundary trapezium included or not. The
  // included case is the arithmetic compute_competition() has always done, in one
  // accumulator, so that path keeps its rounding exactly.
  value_type compute_competition_impl(double height,
                                      bool include_boundary) const;

  // Cache for scan_heights(). Every path that can change a node height must call
  // invalidate_height_scan(); a stale cache here would silently reintroduce the
  // wrong competition profile of #571, so the coverage of these calls was checked
  // by asserting cache == freshly-computed on every call across the whole suite
  // and the scenario gateway.
  HeightScan compute_height_scan() const;
  void invalidate_height_scan() { height_scan_valid = false; }
  mutable HeightScan height_scan_cache{value_type(0.0), true};
  mutable bool height_scan_valid = false;

  // Storage (strategy, nodes) and control() live in SpeciesBase; the
  // using-declarations let the unqualified references below resolve through the
  // dependent base.
  using base_type::nodes;
  using base_type::strategy;
  using base_type::control;
  node_type new_node;

  // The abscissa every reduction over the size distribution is taken over,
  // increasing as the node list is walked from the tallest down. Heights are
  // negated so that both coordinates increase in the same direction; negation is
  // exact, so the height branch's trapezium widths are bit-identical to
  // differencing the heights themselves. Callers in hot loops read the
  // coordinate once and pass it in.
  //
  // It is a position, so it is read at its value even where the height it comes
  // from is active: a quadrature grid is structure. On the birth-date coordinate
  // that is exact, the date being fixed at birth. On the height coordinate it
  // drops the weights' own channel, which the transposes below supply by hand.
  static double abscissa_of(const node_type& n, bool birth_date) {
    using odelia::util::to_passive;
    return birth_date ? to_passive(n.introduction_time())
                      : -to_passive(n.height());
  }
  double quadrature_abscissa(const node_type& n) const {
    return abscissa_of(n, control().node_density_in_birth_date);
  }
  std::vector<double> quadrature_abscissae() const {
    std::vector<double> ret;
    ret.reserve(size());
    const bool birth_date = control().node_density_in_birth_date;
    for (auto& c : nodes) {
      ret.push_back(abscissa_of(c, birth_date));
    }
    return ret;
  }

  typedef typename std::vector<node_type>::iterator nodes_iterator;
  typedef typename std::vector<node_type>::const_iterator nodes_const_iterator;
};

template <typename T, typename E>
Species<T,E>::Species(strategy_type s)
  : base_type(s),
    new_node(this->strategy) {
}

template <typename T, typename E>
Species<T,E>::Species(strategy_type_ptr s)
  : base_type(s),
    new_node(this->strategy) {
}

template <typename T, typename E>
size_t Species<T,E>::size() const {
  return nodes.size();
}

template <typename T, typename E>
void Species<T,E>::clear() {
  invalidate_height_scan();
  nodes.clear();
  // Reset the new_node to a blank new_node, too.
  new_node = node_type(strategy);
}

template <typename T, typename E>
void Species<T,E>::introduce_new_node() {
  invalidate_height_scan();
  // new_node already holds the initial conditions computed against the current
  // environment by the most recent compute_rates() call (see compute_rates ->
  // new_node.compute_initial_conditions above), and the member is refreshed
  // again on the next compute_rates() ready for the following introduction.
  // Recomputing it here would be redundant, and would (wrongly) re-seed against
  // the post-introduction environment rather than the environment at the
  // node's introduction time (resolves the recompute question in #478).
  nodes.push_back(new_node);
}

// If a species contains no individuals, we return the height of a
// seed of the species.  Otherwise we return the height of the largest
// individual, which will be at least as tall as a seed.
//
// This used to return nodes.front(), relying on the decreasing-height ordering
// asserted below. That ordering is guaranteed only while height growth is a
// function of height and the shared environment, which TF24 broke: its
// reserve-gated growth (#517) makes dh/dt depend on a cohort's own storage, so
// two cohorts born moments apart into a rapidly changing environment can cross
// in height. When they had, this returned a height 0.1 m *below* the tallest and
// only living cohort, truncating the light spline's domain (#571). Scanning is
// O(n) in heights only -- negligible against the crown integrals in
// compute_competition -- and returns exactly nodes.front() whenever the ordering
// does hold, so results are unchanged in that case.
template <typename T, typename E>
typename Species<T,E>::value_type Species<T,E>::height_max() const {
  if (nodes.empty()) {
    return new_node.height();
  }
  value_type ret = -std::numeric_limits<double>::infinity();
  for (nodes_const_iterator it = nodes.begin(); it != nodes.end(); ++it) {
    const value_type h = it->height();
    if (h > ret) {
      ret = h;
    }
  }
  return ret;
}

// Are the node heights still ordered largest to smallest? See height_max() above
// for why this can no longer be assumed. Heights only, so this is cheap relative
// to the per-node crown integrals it guards.
template <typename T, typename E>
bool Species<T,E>::heights_are_decreasing() const {
  return scan_heights().decreasing;
}

template <typename T, typename E>
typename Species<T,E>::HeightScan Species<T,E>::scan_heights() const {
  if (!height_scan_valid) {
    height_scan_cache = compute_height_scan();
    height_scan_valid = true;
  }
  return height_scan_cache;
}

// Tallest height and orderedness in one pass over the heights.
template <typename T, typename E>
typename Species<T,E>::HeightScan Species<T,E>::compute_height_scan() const {
  HeightScan ret{value_type(-std::numeric_limits<double>::infinity()), true};
  value_type h_prev = std::numeric_limits<double>::infinity();
  for (nodes_const_iterator it = nodes.begin(); it != nodes.end(); ++it) {
    const value_type h = it->height();
    if (h > h_prev) {
      ret.decreasing = false;
    }
    if (h > ret.h_max) {
      ret.h_max = h;
    }
    h_prev = h;
  }
  return ret;
}

// Because of nodes are always ordered from largest to smallest, we
// need not continue down the list once the leaf area above a certain
// height is zero, because it will be zero for all nodes further down
// the list.
//
// NOTE: This is simply performing numerical integration,  via the
// trapezium rule, of the compute_competition with respect to plant
// height.  You'd think that this would be nicer to do in terms of a
// call to an external trapezium integration function, but building
// and discarding the intermediate storage ends up being a nontrivial
// cost.  A more general iterator version might be possible, but with
// the fiddliness around the boundary conditions that won't likely be
// useful.
//
// NOTE: In the cases where there is no individuals, we return 0 for
// all heights.  The integral is not defined, but an empty light
// environment seems appropriate.
//
// NOTE: A similar early-exit condition to the Plant version is used;
// once the lower bound of the trazpeium is zero, we stop including
// individuals.  Working with the boundary node is tricky here,
// because we might need to include that, too: always in the case of a
// single node (needed to be the second half of the trapezium) and
// also needed if the last looked at plant was still contributing to
// the integral).
template <typename T, typename E>
typename Species<T,E>::value_type
Species<T,E>::compute_competition(double height) const {
  return compute_competition_impl(height, true);
}

// The interior sum alone: never touches new_node, so it is a function of the ODE
// state and the strategy only.
template <typename T, typename E>
typename Species<T,E>::value_type
Species<T,E>::compute_competition_excl_boundary(double height) const {
  return compute_competition_impl(height, false);
}

template <typename T, typename E>
typename Species<T,E>::value_type
Species<T,E>::compute_competition_impl(double height,
                                       bool include_boundary) const {
  if (size() == 0) {
    return value_type(0.0);
  }
  const HeightScan scan = scan_heights();
  if (scan.h_max < height) {
    return value_type(0.0);
  }
  // Read the coordinate once: this is the hottest loop in the solver (one pass
  // per spline knot per Runge-Kutta stage), so the control lookup does not
  // belong inside it.
  const bool birth_date = control().node_density_in_birth_date;
  // The loop below uses the node list itself as the quadrature grid, and the
  // early exit is valid only if that grid is monotone. When it is not, the exit
  // fires at the first node below `height` and silently drops every node beyond
  // it -- including, in #571, the only cohort with non-zero density, which put a
  // fictitious step in the competition profile. Take the ordered path instead.
  // Heights only, so the usual (ordered) case keeps this loop and its results
  // exactly.
  //
  // Only the height abscissa can invert. Introduction times are fixed at birth
  // and nodes are appended in that order, so the birth-date grid is monotone
  // whatever the heights do -- and compute_competition_unordered integrates in
  // height, so sending the birth-date coordinate down it would silently swap
  // coordinates mid-run.
  if (!birth_date && !scan.decreasing) {
    return compute_competition_unordered(height, include_boundary);
  }
  value_type tot = 0.0;
  nodes_const_iterator it = nodes.begin();
  // x1/f1 are the taller end of the interval and its contribution, x0/f0 the
  // lower end and its; each pair comes from one node, and the width multiplies
  // the sum of the two contributions.
  value_type x1 = abscissa_of(*it, birth_date),
             f1 = it->compute_competition(height);

  // Loop over nodes
  for (++it; it != nodes.end(); ++it) {
    const value_type x0 = abscissa_of(*it, birth_date), h0 = it->height(),
                     f0 = it->compute_competition(height);
    if (!util::is_finite(f0)) {
      util::stop("Detected non-finite contribution");
    }
    // Integration
    tot += (x0 - x1) * (f1 + f0);
    // Upper point moves for next time:
    x1 = x0;
    f1 = f0;
    // It is the decreasing height ordering, not the abscissa, that licenses
    // stopping here: every later node is then shorter than `height` and
    // contributes nothing. On the birth-date axis that ordering can break while
    // the abscissa stays monotone, and a node below `height` may be followed by
    // a taller one, so walk the whole list instead. Always true on the height
    // path, which returned above otherwise.
    if (scan.decreasing && h0 < height) {
      break;
    }
  }

  // On the birth-date axis this segment is zero-width at the moment of
  // introduction and contributes nothing once the boundary node is below
  // `height`, so it is always safe to include; f1 can legitimately be zero here
  // when the walk ran to the end.
  if (include_boundary && (size() == 1 || birth_date || f1 > 0)) {
    const value_type x0 = abscissa_of(new_node, birth_date),
                     f0 = new_node.compute_competition(height);
    tot += (x0 - x1) * (f1 + f0);
  }

  return tot / 2;
}

// The same trapezium integral as compute_competition(), and alongside it the
// integral of the vertical derivative, from one traversal of the nodes. Both
// sums visit the same nodes in the same order and associate identically, so the
// first entry is compute_competition(height) bit for bit; a check that it is
// lives in test-canopy-methods.R. The early exit and the closing boundary
// trapezium are driven by the value, as they are there.
template <typename T, typename E>
std::pair<typename Species<T,E>::value_type,
          typename Species<T,E>::value_type>
Species<T,E>::compute_competition_and_slope(double height) const {
  return compute_competition_and_slope_impl(height, true);
}

// The fused pair for the interior sum alone, so the field the boundary condition
// is evaluated in is a function of the ODE state only.
template <typename T, typename E>
std::pair<typename Species<T,E>::value_type,
          typename Species<T,E>::value_type>
Species<T,E>::compute_competition_and_slope_excl_boundary(double height) const {
  return compute_competition_and_slope_impl(height, false);
}

template <typename T, typename E>
std::pair<typename Species<T,E>::value_type,
          typename Species<T,E>::value_type>
Species<T,E>::compute_competition_and_slope_impl(double height,
                                                 bool include_boundary) const {
  if (size() == 0) {
    return {value_type(0.0), value_type(0.0)};
  }
  const HeightScan scan = scan_heights();
  if (scan.h_max < height) {
    return {value_type(0.0), value_type(0.0)};
  }
  const bool birth_date = control().node_density_in_birth_date;
  if (!birth_date && !scan.decreasing) {
    return compute_competition_and_slope_unordered(height, include_boundary);
  }
  value_type tot = 0.0, tot_slope = 0.0;
  nodes_const_iterator it = nodes.begin();
  std::pair<value_type, value_type> fs1 =
    it->compute_competition_and_slope(height);
  value_type x1 = abscissa_of(*it, birth_date), f_h1 = fs1.first,
             s_h1 = fs1.second;

  for (++it; it != nodes.end(); ++it) {
    const std::pair<value_type, value_type> fs0 =
      it->compute_competition_and_slope(height);
    const value_type x0 = abscissa_of(*it, birth_date), h0 = it->height(),
                     f_h0 = fs0.first, s_h0 = fs0.second;
    if (!util::is_finite(f_h0) || !util::is_finite(s_h0)) {
      util::stop("Detected non-finite contribution");
    }
    tot       += (x0 - x1) * (f_h1 + f_h0);
    tot_slope += (x0 - x1) * (s_h1 + s_h0);
    x1   = x0;
    f_h1 = f_h0;
    s_h1 = s_h0;
    if (scan.decreasing && h0 < height) {
      break;
    }
  }

  if (include_boundary && (size() == 1 || birth_date || f_h1 > 0)) {
    const std::pair<value_type, value_type> fs0 =
      new_node.compute_competition_and_slope(height);
    const value_type x0 = abscissa_of(new_node, birth_date);
    tot       += (x0 - x1) * (f_h1 + fs0.first);
    tot_slope += (x0 - x1) * (s_h1 + fs0.second);
  }

  return {tot / 2, tot_slope / 2};
}

// The same trapezium integral as compute_competition(), but over a height-sorted
// view of the nodes rather than the node list in place. Used only when the
// ordering has broken (#571): it agrees with the in-place version whenever the
// ordering holds, so this is a fallback rather than a change of method.
//
// Dropping the zero-density nodes instead would be wrong. A node whose density
// has collapsed to exactly zero contributes f = 0, and that zero is meaningful --
// it is the reconstruction saying density vanishes at that size. Removing those
// grid points would interpolate live density straight across the band and
// overestimate it, so they stay in and the grid gets sorted.
//
// No early exit here: it would need the same monotonicity that is missing. Nodes
// below `height` contribute f = 0 at both ends, so including them costs time but
// changes nothing. The scratch buffer is thread_local and reused, so the repeated
// calls that build one spline do not each allocate.
template <typename T, typename E>
typename Species<T,E>::value_type
Species<T,E>::compute_competition_unordered(double height,
                                            bool include_boundary) const {
  // Cleared on entry, so no element outlives the call that made it.
  thread_local std::vector<std::pair<value_type, value_type>> hf;
  hf.clear();
  hf.reserve(size());

  for (nodes_const_iterator it = nodes.begin(); it != nodes.end(); ++it) {
    const value_type f = it->compute_competition(height);
    if (!util::is_finite(f)) {
      util::stop("Detected non-finite contribution");
    }
    hf.push_back({it->height(), f});
  }
  // Ordering the quadrature grid is structural, so the key is compared and
  // nothing here is differentiated.
  std::sort(hf.begin(), hf.end(),
            [](std::pair<value_type, value_type> const& a,
               std::pair<value_type, value_type> const& b) {
              return a.first > b.first;
            });

  value_type tot = 0.0;
  value_type h1 = hf.front().first, f_h1 = hf.front().second;
  for (size_t j = 1; j < hf.size(); ++j) {
    const value_type h0 = hf[j].first, f_h0 = hf[j].second;
    tot += (h1 - h0) * (f_h1 + f_h0);
    h1   = h0;
    f_h1 = f_h0;
  }

  if (include_boundary && (size() == 1 || f_h1 > 0)) {
    const value_type h0 = new_node.height(),
                     f_h0 = new_node.compute_competition(height);
    tot += (h1 - h0) * (f_h1 + f_h0);
  }

  return tot / 2;
}

// compute_competition_unordered() with the vertical derivative alongside it. The
// sort key and its tie-break are the same, so the two sums here also visit the
// same nodes in the same order.
template <typename T, typename E>
std::pair<typename Species<T,E>::value_type,
          typename Species<T,E>::value_type>
Species<T,E>::compute_competition_and_slope_unordered(double height,
                                                      bool include_boundary) const {
  // Cleared on entry, so no element outlives the call that made it.
  thread_local std::vector<
    std::pair<value_type, std::pair<value_type, value_type>>> hfs;
  hfs.clear();
  hfs.reserve(size());

  for (nodes_const_iterator it = nodes.begin(); it != nodes.end(); ++it) {
    const std::pair<value_type, value_type> fs =
      it->compute_competition_and_slope(height);
    if (!util::is_finite(fs.first) || !util::is_finite(fs.second)) {
      util::stop("Detected non-finite contribution");
    }
    hfs.push_back({it->height(), fs});
  }
  std::sort(hfs.begin(), hfs.end(),
            [](std::pair<value_type, std::pair<value_type, value_type>> const& a,
               std::pair<value_type, std::pair<value_type, value_type>> const& b) {
              return a.first > b.first;
            });

  value_type tot = 0.0, tot_slope = 0.0;
  value_type h1 = hfs.front().first;
  value_type f_h1 = hfs.front().second.first, s_h1 = hfs.front().second.second;
  for (size_t j = 1; j < hfs.size(); ++j) {
    const value_type h0 = hfs[j].first;
    const value_type f_h0 = hfs[j].second.first, s_h0 = hfs[j].second.second;
    tot       += (h1 - h0) * (f_h1 + f_h0);
    tot_slope += (h1 - h0) * (s_h1 + s_h0);
    h1   = h0;
    f_h1 = f_h0;
    s_h1 = s_h0;
  }

  if (include_boundary && (size() == 1 || f_h1 > 0)) {
    const std::pair<value_type, value_type> fs0 =
      new_node.compute_competition_and_slope(height);
    const value_type h0 = new_node.height();
    tot       += (h1 - h0) * (f_h1 + fs0.first);
    tot_slope += (h1 - h0) * (s_h1 + fs0.second);
  }

  return {tot / 2, tot_slope / 2};
}

template <typename T, typename E>
void Species<T,E>::compute_rates(const E& environment, double pr_patch_survival, double birth_rate) {
  for (auto& c : nodes) {
    c.compute_rates(environment, pr_patch_survival);
  }
  // The boundary condition, evaluated in the field the nodes above were just
  // rated in. This is not the evaluation the field itself reads -- that one is in
  // a field excluding the boundary interval, and Patch::compute_environment owns
  // it -- so the two are the same function at different arguments rather than one
  // computed twice. This value is the one an introduced node inherits.
  new_node.compute_initial_conditions(environment, pr_patch_survival, birth_rate);
  if (internals::transport_census_active()) {
    // The sub-grid value is recovered from the rate each node has just written,
    // log_density_dt = -growth_rate_gradient - mortality, so the census adds no
    // rate evaluation of its own.
    internals::transport_census& census = internals::the_transport_census();
    // A census is read rather than differentiated, so each column is recorded at
    // its value.
    using odelia::util::to_passive;
    for (std::size_t i = 0; i < nodes.size(); ++i) {
      const node_type& below = i + 1 < nodes.size() ? nodes[i + 1] : new_node;
      census.add(environment.time, to_passive(nodes[i].height()),
                 to_passive(nodes[i].height() - below.height()),
                 to_passive(nodes[i].growth_rate()),
                 to_passive(below.growth_rate()),
                 to_passive(-(nodes[i].get_log_density_rate() +
                              nodes[i].mortality_rate())),
                 to_passive(growth_rate_gradient(i)));
    }
  }
}

template <typename T, typename E>
typename Species<T,E>::value_type
Species<T,E>::growth_rate_gradient(std::size_t i) const {
  const node_type& below = i + 1 < size() ? nodes[i + 1] : new_node;
  const value_type dh = nodes[i].height() - below.height();
  if (dh == 0.0) {
    return i > 0 ? growth_rate_gradient(i - 1) : value_type(0.0);
  }
  return (nodes[i].growth_rate() - below.growth_rate()) / dh;
}

template <typename T, typename E>
void Species<T,E>::introduce_new_node(double time, double patch_density) {
  invalidate_height_scan();
  // Stamp the pushed copy (not new_node) so the member stays pristine for
  // the no-arg introduction paths.
  nodes.push_back(new_node);
  nodes.back().set_introduction(time, patch_density);
}

template <typename T, typename E>
void Species<T,E>::remove_newest_node() {
  if (nodes.empty()) {
    util::stop("no node to remove from this species");
  }
  invalidate_height_scan();
  nodes.pop_back();
}

template <typename T, typename E>
std::vector<double> Species<T,E>::net_reproduction_ratio_by_node() const {
  std::vector<double> ret;
  ret.reserve(size());
  for (auto& c : nodes) {
    ret.push_back(odelia::util::to_passive(c.fecundity()));
  }
  return ret;
}

template <typename T, typename E>
std::vector<double> Species<T,E>::net_reproduction_ratio_by_node_weighted() const {
  std::vector<double> ret;
  ret.reserve(size());
  for (auto& c : nodes) {
    ret.push_back(
      odelia::util::to_passive(c.weighted_fecundity(strategy->pars.S_D)));
  }
  return ret;
}

template <typename T, typename E>
std::vector<double> Species<T,E>::node_times() const {
  std::vector<double> ret;
  ret.reserve(size());
  for (auto& c : nodes) {
    ret.push_back(c.introduction_time());
  }
  return ret;
}

template <typename T, typename E>
bool Species<T,E>::birth_dates_are_distinct() const {
  for (size_t i = 1; i < nodes.size(); ++i) {
    if (nodes[i].introduction_time() == nodes[i - 1].introduction_time()) {
      return false;
    }
  }
  return true;
}

template <typename T, typename E>
void Species<T,E>::resize_consumption_rates(int r) {
  new_node.resize_consumption_rates(r);
}

template <typename T, typename E>
typename Species<T,E>::value_type
Species<T,E>::consumption_rate(int i) const {
  if (size() == 0) {
    return value_type(0.0);
  }
  if (control().node_density_in_birth_date) {
    // Introduction times are fixed at birth and nodes are appended in that
    // order, so this grid is ascending however the heights behave -- there is no
    // inverted case to sort. new_node's birth date is the current time, which is
    // the newest, so it goes on the end rather than the front.
    // The abscissa is double and the integrand is not: birth dates are fixed at
    // birth and carry no derivative, so the quadrature weights are constants,
    // while the rates are what the uptake's sensitivity runs through.
    std::vector<double> times = node_times();
    times.push_back(new_node.introduction_time());
    std::vector<value_type> rates = consumption_rate_by_node(i);
    rates.push_back(new_node.consumption_rate(i));
    return util::trapezium(times, rates);
  }
  // node heights are in descending order - we need ascending for integration,
  // starting at new_node, which is where the size distribution starts. The
  // heights are the integration grid; the rates integrated over it are what the
  // uptake depends on.
  std::vector<value_type> heights;
  heights.reserve(size() + 1);
  heights.push_back(new_node.height());
  for (auto it = nodes.rbegin(); it != nodes.rend(); ++it) {
    heights.push_back(it->height());
  }
  std::vector<value_type> rates = consumption_rate_by_node_rev(i);

  // The node list is the quadrature grid here as it is in compute_competition,
  // so an inverted grid (#571) makes neighbouring trapezia cancel rather than
  // accumulate. Sort the pairs when the ordering has broken, as
  // compute_competition_unordered does; an already-ascending grid is untouched.
  // The order is decided on the passive value: an ordering chosen on an active
  // key would make the recorded computation depend on the state.
  auto below = [](const value_type& a, const value_type& b) -> bool {
    return odelia::util::to_passive(a) < odelia::util::to_passive(b);
  };
  if (!std::is_sorted(heights.begin(), heights.end(), below)) {
    std::vector<size_t> order(heights.size());
    for (size_t j = 0; j < order.size(); ++j) {
      order[j] = j;
    }
    std::sort(order.begin(), order.end(), [&](size_t a, size_t b) -> bool {
      return below(heights[a], heights[b]);
    });
    std::vector<value_type> h_sorted, r_sorted;
    h_sorted.reserve(heights.size());
    r_sorted.reserve(rates.size());
    for (size_t j : order) {
      h_sorted.push_back(heights[j]);
      r_sorted.push_back(rates[j]);
    }
    heights.swap(h_sorted);
    rates.swap(r_sorted);
  }
  return util::trapezium(heights, rates);
}

template <typename T, typename E>
std::vector<typename Species<T,E>::value_type>
Species<T,E>::consumption_rate_by_node_rev(int i) const {
  std::vector<value_type> ret;
  ret.reserve(size() + 1);
  ret.push_back(new_node.consumption_rate(i));
  for(auto it = nodes.rbegin(); it != nodes.rend(); ++it) {
    ret.push_back(it->consumption_rate(i));
  }
  return ret;
}



template <typename T, typename E>
std::vector<typename Species<T,E>::value_type>
Species<T,E>::consumption_rate_by_node(int i) const {
  std::vector<value_type> ret;
  ret.reserve(size());
  for(auto& c : nodes) {
    ret.push_back(c.consumption_rate(i));
  }
  return ret;
}


template <typename T, typename E>
template <class Psi>
typename Species<T,E>::value_type
Species<T,E>::census(Psi psi) const {
  if (size() == 0) {
    return value_type(0.0);
  }
  // The abscissa ascends as the nodes are walked in storage order and the
  // boundary node closes the grid: it is the youngest, and the shortest.
  const bool birth_date = control().node_density_in_birth_date;
  if (!birth_date && !heights_are_decreasing()) {
    util::stop("The census needs the node heights in decreasing order; on a"
               " crossed grid neighbouring trapezia cancel instead of"
               " accumulating");
  }
  std::vector<double> x;
  std::vector<value_type> weighted;
  x.reserve(size() + 1);
  weighted.reserve(size() + 1);
  for (auto& c : nodes) {
    x.push_back(abscissa_of(c, birth_date));
    weighted.push_back(c.get_density() * psi(*strategy, c.individual));
  }
  x.push_back(abscissa_of(new_node, birth_date));
  weighted.push_back(new_node.get_density() *
                     psi(*strategy, new_node.individual));
  return util::trapezium(x, weighted);
}

// Mirrors compute_competition_and_slope_impl(height, true) term by term: the
// same early exit, the same closing trapezium, the same node order.
template <typename T, typename E>
void Species<T,E>::compute_competition_and_slope_adjoint(
    double height, double lambda_value, double lambda_slope,
    node_size_adjoints* out) const {
  using odelia::util::to_passive;
  if (size() == 0) {
    return;
  }
  const HeightScan scan = scan_heights();
  if (to_passive(scan.h_max) < height) {
    return;
  }
  const bool birth_date = control().node_density_in_birth_date;
  if (!birth_date && !scan.decreasing) {
    util::stop("The competition adjoint needs the node heights in decreasing"
               " order; the sorted-view reduction has no transpose here");
  }
  // The forward halves both sums once at the end. One slot per grid point, so
  // the boundary node -- the distribution's lower point -- carries its own row
  // rather than only lending its abscissa to its neighbour's weight.
  const double lv = lambda_value * 0.5, ls = lambda_slope * 0.5;
  const size_t n_slot = size() + 1;
  std::vector<double> lambda_f(n_slot, 0.0), lambda_s(n_slot, 0.0);

  std::pair<value_type, value_type> fs1 =
    nodes[0].compute_competition_and_slope(height);
  size_t upper = 0, last = 0;
  double x1 = abscissa_of(nodes[0], birth_date);
  double f_h1 = to_passive(fs1.first), s_h1 = to_passive(fs1.second);

  for (size_t k = 1; k < size(); ++k) {
    const std::pair<value_type, value_type> fs0 =
      nodes[k].compute_competition_and_slope(height);
    const double x0 = abscissa_of(nodes[k], birth_date);
    const double h0 = to_passive(nodes[k].height());
    const double f_h0 = to_passive(fs0.first), s_h0 = to_passive(fs0.second);
    const double width = x0 - x1;
    if (!birth_date) {
      // On the height abscissa the interval's width is built from two node
      // heights, so the quadrature moves with them and not only the integrand
      // does. A birth date is fixed at birth, so there the width is a constant
      // and this term does not exist.
      const double edge = lv * (f_h1 + f_h0) + ls * (s_h1 + s_h0);
      out[upper].height += edge;
      out[k].height     -= edge;
    }
    lambda_f[upper] += lv * width;
    lambda_s[upper] += ls * width;
    lambda_f[k]     += lv * width;
    lambda_s[k]     += ls * width;
    upper = k; last = k;
    x1 = x0; f_h1 = f_h0; s_h1 = s_h0;
    if (scan.decreasing && h0 < height) {
      break;
    }
  }

  if (size() == 1 || birth_date || f_h1 > 0) {
    // The boundary node's own height and density are not ODE state.
    const std::pair<value_type, value_type> fs0 =
      new_node.compute_competition_and_slope(height);
    const double x0 = abscissa_of(new_node, birth_date);
    if (!birth_date) {
      out[upper].height += lv * (f_h1 + to_passive(fs0.first)) +
                           ls * (s_h1 + to_passive(fs0.second));
    }
    lambda_f[upper] += lv * (x0 - x1);
    lambda_s[upper] += ls * (x0 - x1);
    lambda_f[size()] += lv * (x0 - x1);
    lambda_s[size()] += ls * (x0 - x1);
  }

  for (size_t k = 0; k < n_slot; ++k) {
    if (lambda_f[k] == 0.0 && lambda_s[k] == 0.0) {
      continue;
    }
    const node_type& node = k + 1 < n_slot ? nodes[k] : new_node;
    const typename node_type::competition_partials p =
      node.compute_competition_and_slope_partials(height);
    out[k].area_leaf += lambda_f[k] * to_passive(p.value_darea_leaf) +
                        lambda_s[k] * to_passive(p.slope_darea_leaf);
    out[k].height    += lambda_f[k] * to_passive(p.value_dheight) +
                        lambda_s[k] * to_passive(p.slope_dheight);
    out[k].log_density += lambda_f[k] * to_passive(p.value_dlog_density) +
                          lambda_s[k] * to_passive(p.slope_dlog_density);
    out[k].extinction += lambda_f[k] * to_passive(p.value_dk_I) +
                         lambda_s[k] * to_passive(p.slope_dk_I);
  }
}

// Mirrors consumption_rate: the grid runs from new_node upwards, so array slot
// j holds node size() - j and slot 0 the boundary node, which is not state.
template <typename T, typename E>
void Species<T,E>::consumption_rate_adjoint(int resource, double lambda_uptake,
                                            node_uptake_adjoints* out) const {
  using odelia::util::to_passive;
  if (size() == 0) {
    return;
  }
  // The grid is the quadrature abscissa, which ascends as the nodes are walked
  // in storage order and puts the boundary node last on either coordinate: it is
  // the youngest, and the shortest. So slot j holds node j, slot n-1 holds the
  // boundary node, and that last one is not ODE state.
  const bool birth_date = control().node_density_in_birth_date;
  if (!birth_date && !heights_are_decreasing()) {
    util::stop("The uptake adjoint needs the node heights in decreasing order;"
               " the sorted-view reduction has no transpose here");
  }
  const size_t n = size() + 1;
  std::vector<double> x(n), y(n);
  for (size_t j = 0; j + 1 < n; ++j) {
    x[j] = abscissa_of(nodes[j], birth_date);
    y[j] = to_passive(nodes[j].consumption_rate(resource));
  }
  x[n - 1] = abscissa_of(new_node, birth_date);
  y[n - 1] = to_passive(new_node.consumption_rate(resource));

  std::vector<double> lambda_x(n, 0.0), lambda_y(n, 0.0);
  for (size_t j = 0; j + 1 < n; ++j) {
    const double half = 0.5 * lambda_uptake;
    if (!birth_date) {
      // The height abscissa is state, so the width carries a derivative. A
      // birth date is fixed at birth and the width is a constant.
      lambda_x[j]     -= half * (y[j + 1] + y[j]);
      lambda_x[j + 1] += half * (y[j + 1] + y[j]);
    }
    lambda_y[j]     += half * (x[j + 1] - x[j]);
    lambda_y[j + 1] += half * (x[j + 1] - x[j]);
  }
  // Every slot, boundary node included. It is the distribution's lower grid
  // point, so a transpose stopping at the introduced nodes is the transpose of a
  // reduction the forward model is not computing -- and a quadrature node is not
  // a term that can be dropped, because its weight is shared with its neighbour.
  for (size_t j = 0; j < n; ++j) {
    const node_type& node = j + 1 < n ? nodes[j] : new_node;
    const double density = to_passive(node.get_density());
    out[j].uptake      += lambda_y[j] * density;
    out[j].log_density += lambda_y[j] * y[j];
    // lambda_x is taken in the abscissa, and on the height coordinate the
    // abscissa is minus the height.
    out[j].height      -= lambda_x[j];
  }
}

// bit clunky...
template <typename T, typename E>
size_t Species<T,E>::aux_size() const {
  return size() * strategy->aux_size();
}

template <typename T, typename E>
template <typename It>
It Species<T,E>::ode_aux(It it) const {
  return odelia::ode::ode_aux(nodes.begin(), nodes.end(), it);
}

template <typename T, typename E>
Rcpp::NumericMatrix Species<T, E>::r_get_state() const {

  size_t ode_size = node_type::ode_size(), n_nodes = size();
  size_t aux_size = strategy->aux_size();

  // On the birth-date path `log_density` is converted to the density in height
  // before it leaves C++, so every downstream consumer of this matrix keeps its
  // meaning, and the quantity actually integrated is reported alongside it as
  // `log_density_state`. Note export_patch_state() resumes from patch$ode_state,
  // not from here, so the raw state is what a resume reloads.
  const size_t extra = control().node_density_in_birth_date ? 1 : 0;

  // Set output size. // +1 is seed
  Rcpp::NumericMatrix ret(static_cast<int>(ode_size + aux_size + extra), n_nodes + 1);
  Rcpp::NumericMatrix::iterator it = ret.begin();

  for (size_t i = 0; i < n_nodes; ++i)
  {
    it = get_node_state(nodes[i], it);
    it = get_node_aux(nodes[i], it);
    it += extra;
  }

  it = get_node_state(new_node, it);
  it = get_node_aux(new_node, it);
  it += extra;

  // Combine ode_names and aux_names into a single vector for dimnames
  std::vector<std::string> names = node_type::ode_names();
  std::vector<std::string> aux = strategy->aux_names();
  names.insert(names.end(), aux.begin(), aux.end());

  if (extra > 0) {
    const int ld = static_cast<int>(
      std::find(names.begin(), names.end(), "log_density") - names.begin());
    const int st = static_cast<int>(ode_size + aux_size);
    names.push_back("log_density_state");
    // This matrix includes the boundary node as its last column, so the
    // Jacobian's trailing entry is used here (unlike r_log_densities()).
    const std::vector<double> jac = height_jacobian();
    for (int col = 0; col <= static_cast<int>(n_nodes); ++col) {
      const double nu = ret(ld, col);
      ret(st, col) = nu;
      ret(ld, col) = util::is_finite(jac[col]) ? nu - std::log(jac[col]) : NA_REAL;
    }
  }

  ret.attr("dimnames") = Rcpp::List::create(names, R_NilValue);

  return ret;
}

template <typename T, typename E>
std::vector<double> Species<T,E>::r_heights() const {
  std::vector<double> ret;
  ret.reserve(size());
  for (nodes_const_iterator it = nodes.begin();
       it != nodes.end(); ++it) {
    ret.push_back(odelia::util::to_passive(it->height()));
  }
  return ret;
}

template <typename T, typename E>
std::vector<double> Species<T,E>::r_heights_rev() const {
  std::vector<double> ret;
  ret.reserve(size());
  for (nodes_const_iterator it = nodes.begin();
       it != nodes.end(); ++it) {
    ret.push_back(odelia::util::to_passive(it->height()));
  }
  std::reverse(ret.begin(), ret.end());
  return ret;
}

template <typename T, typename E>
void Species<T,E>::r_set_heights(std::vector<double> heights) {
  invalidate_height_scan();
  util::check_length(heights.size(), size());
  if (!util::is_decreasing(heights.begin(), heights.end())) {
    util::stop("height must be decreasing (ties allowed)");
  }
  size_t i = 0;
  for (nodes_iterator it = nodes.begin(); it != nodes.end(); ++it, ++i) {
    it->individual.set_state("height", heights[i]);
  }
}

template <typename T, typename E>
std::vector<double> Species<T,E>::r_compute_competition_effect_by_nodes() const {
  std::vector<double> ret;
  ret.reserve(size());
  for (auto& c : nodes) {
    ret.push_back(odelia::util::to_passive(c.compute_competition(0.0)));
  }
  return ret;
}

template <typename T, typename E>
std::vector<double> Species<T,E>::r_compute_competition_effect_by_nodes_error(double scal) const {
  // Over the same abscissa the competition integral uses, so schedule
  // refinement measures the error of the quadrature actually being taken.
  // local_error_integration takes absolute differences, so the height branch's
  // sign flip leaves it unchanged.
  return util::local_error_integration(quadrature_abscissae(),
                                       r_compute_competition_effect_by_nodes(), scal);
}

// Central differences in the interior, one-sided at the two ends. The boundary
// node supplies the extra point at the young end, so the youngest real node --
// the one whose Jacobian a forward difference gets worst -- gets a genuine
// central difference. NA where the Jacobian vanishes: two cohorts at the same
// height is exactly where the height density is undefined (it is the multivalued
// case), and reporting +Inf there would be worse than reporting nothing.
template <typename T, typename E>
std::vector<double> Species<T,E>::height_jacobian() const {
  const size_t n = size();
  std::vector<double> h(n + 1), t(n + 1);
  for (size_t i = 0; i < n; ++i) {
    h[i] = nodes[i].height();
    t[i] = nodes[i].introduction_time();
  }
  h[n] = new_node.height();
  t[n] = new_node.introduction_time();

  std::vector<double> ret(n + 1, NA_REAL);

  // The boundary node needs no difference at all: dh/dtau = -g(H_0) at birth,
  // and compute_initial_conditions() re-evaluates that every step, so the value
  // it recorded is both exact and current. Only for *this* node -- an introduced
  // one has aged since, and its recorded rate is frozen at its own birth.
  const double g0 = new_node.growth_rate_at_birth();
  if (util::is_finite(g0) && g0 > 0.0) {
    ret[n] = g0;
  }

  if (n < 1) {
    return ret; // no interior to difference
  }
  for (size_t j = 0; j < n; ++j) {
    const size_t lo = (j == 0) ? 0 : j - 1;
    const size_t hi = j + 1;
    const double jac = std::abs((h[hi] - h[lo]) / (t[hi] - t[lo]));
    if (util::is_finite(jac) && jac > 0.0) {
      ret[j] = jac;
    }
  }
  if (!util::is_finite(ret[n])) {
    // Fall back to a one-sided difference if the birth rate is unavailable
    // (a node loaded from an exported state records zero).
    const double jac = std::abs((h[n] - h[n - 1]) / (t[n] - t[n - 1]));
    if (util::is_finite(jac) && jac > 0.0) {
      ret[n] = jac;
    }
  }
  return ret;
}

template <typename T, typename E>
std::vector<double> Species<T,E>::r_log_densities_state() const {
  std::vector<double> ret;
  ret.reserve(size());
  for (nodes_const_iterator it = nodes.begin();
       it != nodes.end(); ++it) {
    ret.push_back(odelia::util::to_passive(it->get_log_density()));
  }
  return ret;
}

template <typename T, typename E>
std::vector<double> Species<T,E>::r_log_densities() const {
  std::vector<double> ret = r_log_densities_state();
  if (!control().node_density_in_birth_date) {
    return ret;
  }
  // N = nu / |dh/dtau|. jac carries one extra trailing entry for the boundary
  // node, which this accessor does not report.
  const std::vector<double> jac = height_jacobian();
  for (size_t i = 0; i < ret.size(); ++i) {
    ret[i] = util::is_finite(jac[i]) ? ret[i] - std::log(jac[i]) : NA_REAL;
  }
  return ret;
}

template <typename T, typename E>
std::vector<double> Species<T,E>::r_log_density_rates() const {
  std::vector<double> ret;
  ret.reserve(size());
  for (nodes_const_iterator it = nodes.begin(); it != nodes.end(); ++it) {
    ret.push_back(odelia::util::to_passive(it->get_log_density_rate()));
  }
  return ret;
}

template <typename T, typename E>
std::vector<double> Species<T,E>::r_patch_densities() const {
  std::vector<double> ret;
  ret.reserve(size());
  for (nodes_const_iterator it = nodes.begin(); it != nodes.end(); ++it) {
    ret.push_back(it->patch_density());
  }
  return ret;
}

template <typename T, typename E>
std::vector<double> Species<T,E>::r_pr_patch_survival_at_birth() const {
  std::vector<double> ret;
  ret.reserve(size());
  for (nodes_const_iterator it = nodes.begin(); it != nodes.end(); ++it) {
    ret.push_back(it->get_pr_patch_survival_at_birth());
  }
  return ret;
}

template <typename T, typename E>
void Species<T,E>::set_birth_state(const std::vector<double>& times,
                                   const std::vector<double>& patch_density,
                                   const std::vector<double>& pr_patch_survival) {
  util::check_length(times.size(), size());
  util::check_length(patch_density.size(), size());
  util::check_length(pr_patch_survival.size(), size());
  for (size_t i = 0; i < size(); ++i) {
    nodes[i].set_birth_state(times[i], patch_density[i], pr_patch_survival[i]);
  }
}

}

#endif /* SPECIES */
