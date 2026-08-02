// -*-c++-*-
#ifndef SPECIES
#define SPECIES

#include <vector>
#include <algorithm>
#include <limits>
#include <utility>
#include <plant/util.h>
#include <plant/environment.h>
#include <odelia/ode_interface.hpp>
#include <plant/node.h>
#include <plant/species_base.h>
#include <plant/transport_census.h>
#include <odelia/drivers.hpp>

namespace plant {

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

  // Transpose of compute_competition_and_slope at `height`, with the closing
  // boundary trapezium the field build includes. `out` points at this species'
  // first node. The boundary node's own height and density are not ODE state.
  void compute_competition_and_slope_adjoint(double height,
                                             double lambda_value,
                                             double lambda_slope,
                                             node_size_adjoints* out) const;

  // Transpose of consumption_rate for one resource. The trapezium's grid is the
  // node heights, so the weights reach the heights as well as the rates.
  void consumption_rate_adjoint(int resource, double lambda_uptake,
                                node_uptake_adjoints* out) const;

  // Evaluate the inflow boundary condition in the environment passed. Split out
  // of compute_rates() so the field build owns it and the field stops reading a
  // density carried from the previous evaluation.
  void compute_boundary_node(const environment_type& environment,
                             double pr_patch_survival, double birth_rate) {
    new_node.compute_initial_conditions(environment, pr_patch_survival, birth_rate);
  }

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

  // * ODE interface
  // NOTE: We are a time-independent model here so no need to pass
  // time in as an argument.  All the bits involving time are taken
  // care of by Environment for us.
  // (ode_size/set_ode_state/ode_state/ode_rates come from SpeciesBase.)
  size_t aux_size() const;

  void resize_consumption_rates(int i);
  value_type consumption_rate(int i) const;
  std::vector<value_type> consumption_rate_by_node_rev(int i) const;

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

  // This is just kind of useful
  std::vector<double> r_log_densities() const;
  // Per-node rate of change of log density; used to guard against initial
  // conditions whose densities would explode to non-finite values.
  std::vector<double> r_log_density_rates() const;

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

private:
  // compute_competition() for the case where the node heights are no longer
  // ordered, so the node list cannot be used directly as the quadrature grid.
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

  typedef typename std::vector<node_type>::iterator nodes_iterator;
  typedef typename std::vector<node_type>::const_iterator nodes_const_iterator;
};

template <typename T, typename E>
Species<T,E>::Species(strategy_type s)
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
  // The loop below uses the node list itself as the quadrature grid, and the
  // early exit is valid only if that grid is monotone. When it is not, the exit
  // fires at the first node below `height` and silently drops every node beyond
  // it -- including, in #571, the only cohort with non-zero density, which put a
  // fictitious step in the competition profile. Take the ordered path instead.
  // Heights only, so the usual (ordered) case keeps this loop and its results
  // exactly.
  if (!scan.decreasing) {
    return compute_competition_unordered(height, include_boundary);
  }
  value_type tot = 0.0;
  nodes_const_iterator it = nodes.begin();
  // h1/f_h1 are the taller end of the interval and its contribution, h0/f_h0 the
  // lower end and its; each pair comes from one node, and the width multiplies
  // the sum of the two contributions.
  value_type h1 = it->height(), f_h1 = it->compute_competition(height);

  // Loop over nodes
  for (++it; it != nodes.end(); ++it) {
    const value_type h0 = it->height(), f_h0 = it->compute_competition(height);
    if (!util::is_finite(f_h0)) {
      util::stop("Detected non-finite contribution");
    }
    // Integration
    tot += (h1 - h0) * (f_h1 + f_h0);
    // Upper point moves for next time:
    h1   = h0;
    f_h1 = f_h0;
    if (h0 < height) {
      break;
    }
  }

  if (include_boundary && (size() == 1 || f_h1 > 0)) {
    const value_type h0 = new_node.height(),
                     f_h0 = new_node.compute_competition(height);
    tot += (h1 - h0) * (f_h1 + f_h0);
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
  if (!scan.decreasing) {
    return compute_competition_and_slope_unordered(height, include_boundary);
  }
  value_type tot = 0.0, tot_slope = 0.0;
  nodes_const_iterator it = nodes.begin();
  std::pair<value_type, value_type> fs1 =
    it->compute_competition_and_slope(height);
  value_type h1 = it->height(), f_h1 = fs1.first, s_h1 = fs1.second;

  for (++it; it != nodes.end(); ++it) {
    const std::pair<value_type, value_type> fs0 =
      it->compute_competition_and_slope(height);
    const value_type h0 = it->height(), f_h0 = fs0.first, s_h0 = fs0.second;
    if (!util::is_finite(f_h0) || !util::is_finite(s_h0)) {
      util::stop("Detected non-finite contribution");
    }
    tot       += (h1 - h0) * (f_h1 + f_h0);
    tot_slope += (h1 - h0) * (s_h1 + s_h0);
    h1   = h0;
    f_h1 = f_h0;
    s_h1 = s_h0;
    if (h0 < height) {
      break;
    }
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
void Species<T,E>::resize_consumption_rates(int r) {
  new_node.resize_consumption_rates(r);
}

template <typename T, typename E>
typename Species<T,E>::value_type
Species<T,E>::consumption_rate(int i) const {
  if (size() == 0) {
    return value_type(0.0);
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
  return util::trapezium(heights, consumption_rate_by_node_rev(i));
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
  if (!scan.decreasing) {
    util::stop("The competition adjoint needs the node heights in decreasing"
               " order; the sorted-view reduction has no transpose here");
  }
  // The forward halves both sums once at the end.
  const double lv = lambda_value * 0.5, ls = lambda_slope * 0.5;
  std::vector<double> lambda_f(size(), 0.0), lambda_s(size(), 0.0);

  std::pair<value_type, value_type> fs1 =
    nodes[0].compute_competition_and_slope(height);
  size_t upper = 0, last = 0;
  double h1 = to_passive(nodes[0].height());
  double f_h1 = to_passive(fs1.first), s_h1 = to_passive(fs1.second);

  for (size_t k = 1; k < size(); ++k) {
    const std::pair<value_type, value_type> fs0 =
      nodes[k].compute_competition_and_slope(height);
    const double h0 = to_passive(nodes[k].height());
    const double f_h0 = to_passive(fs0.first), s_h0 = to_passive(fs0.second);
    const double width = h1 - h0;
    // The interval's width is built from two node heights, so the quadrature
    // moves with them and not only the integrand does.
    const double edge = lv * (f_h1 + f_h0) + ls * (s_h1 + s_h0);
    out[upper].height += edge;
    out[k].height     -= edge;
    lambda_f[upper] += lv * width;
    lambda_s[upper] += ls * width;
    lambda_f[k]     += lv * width;
    lambda_s[k]     += ls * width;
    upper = k; last = k;
    h1 = h0; f_h1 = f_h0; s_h1 = s_h0;
    if (h0 < height) {
      break;
    }
  }

  if (size() == 1 || f_h1 > 0) {
    const std::pair<value_type, value_type> fs0 =
      new_node.compute_competition_and_slope(height);
    const double h0 = to_passive(new_node.height());
    out[upper].height += lv * (f_h1 + to_passive(fs0.first)) +
                         ls * (s_h1 + to_passive(fs0.second));
    lambda_f[upper] += lv * (h1 - h0);
    lambda_s[upper] += ls * (h1 - h0);
  }

  for (size_t k = 0; k <= last; ++k) {
    if (lambda_f[k] == 0.0 && lambda_s[k] == 0.0) {
      continue;
    }
    const typename node_type::competition_partials p =
      nodes[k].compute_competition_and_slope_partials(height);
    out[k].area_leaf += lambda_f[k] * to_passive(p.value_darea_leaf) +
                        lambda_s[k] * to_passive(p.slope_darea_leaf);
    out[k].height    += lambda_f[k] * to_passive(p.value_dheight) +
                        lambda_s[k] * to_passive(p.slope_dheight);
    out[k].log_density += lambda_f[k] * to_passive(p.value_dlog_density) +
                          lambda_s[k] * to_passive(p.slope_dlog_density);
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
  const size_t n = size() + 1;
  std::vector<double> x(n), y(n);
  x[0] = to_passive(new_node.height());
  y[0] = to_passive(new_node.consumption_rate(resource));
  for (size_t j = 1; j < n; ++j) {
    x[j] = to_passive(nodes[size() - j].height());
    y[j] = to_passive(nodes[size() - j].consumption_rate(resource));
  }
  std::vector<double> lambda_x(n, 0.0), lambda_y(n, 0.0);
  for (size_t j = 0; j + 1 < n; ++j) {
    const double half = 0.5 * lambda_uptake;
    lambda_x[j]     -= half * (y[j + 1] + y[j]);
    lambda_x[j + 1] += half * (y[j + 1] + y[j]);
    lambda_y[j]     += half * (x[j + 1] - x[j]);
    lambda_y[j + 1] += half * (x[j + 1] - x[j]);
  }
  for (size_t j = 1; j < n; ++j) {
    const size_t k = size() - j;
    const double density = to_passive(nodes[k].get_density());
    out[k].uptake      += lambda_y[j] * density;
    out[k].log_density += lambda_y[j] * y[j];
    out[k].height      += lambda_x[j];
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

  // Set output size. // +1 is seed
  Rcpp::NumericMatrix ret(static_cast<int>(ode_size + aux_size), n_nodes + 1); 
  Rcpp::NumericMatrix::iterator it = ret.begin();
  
  for (size_t i = 0; i < n_nodes; ++i)
  {
    it = get_node_state(nodes[i], it);
    it = get_node_aux(nodes[i], it);
  }

  it = get_node_state(new_node, it);
  it = get_node_aux(new_node, it);

  // Combine ode_names and aux_names into a single vector for dimnames
  std::vector<std::string> names = node_type::ode_names();
  std::vector<std::string> aux = strategy->aux_names();
  names.insert(names.end(), aux.begin(), aux.end());

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
  return util::local_error_integration(r_heights(), r_compute_competition_effect_by_nodes(), scal);
}

template <typename T, typename E>
std::vector<double> Species<T,E>::r_log_densities() const {
  std::vector<double> ret;
  ret.reserve(size());
  for (nodes_const_iterator it = nodes.begin();
       it != nodes.end(); ++it) {
    ret.push_back(odelia::util::to_passive(it->get_log_density()));
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
