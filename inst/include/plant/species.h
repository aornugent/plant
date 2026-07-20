// -*-c++-*-
#ifndef SPECIES
#define SPECIES

#include <vector>
#include <plant/util.h>
#include <plant/environment.h>
#include <odelia/ode_interface.hpp>
#include <odelia/ode_util.hpp> // odelia::util::diagnostic (intent-named AD strip)
#include <odelia/mass_transport.hpp>
#include <plant/node.h>
#include <plant/species_base.h>
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
  typedef T         strategy_type;
  typedef E         environment_type;
  typedef Individual<T,E>  individual_type;
  typedef Node<T,E> node_type;
  typedef typename strategy_type::ptr strategy_type_ptr;
  // The scalar this species' demographic state carries (double / active).
  using value_type = typename base_type::value_type;
  Species(strategy_type s);

  // ODE plumbing and the per-element serialisers are inherited from SpeciesBase
  // and iterate all nodes (the deterministic model has no notion of "dead").
  using base_type::ode_size;
  using base_type::set_ode_state;
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

  double height_max() const;
  value_type compute_competition(double height) const;
  // Population census: the mass-weighted reduction Sum_i n_i * psi(individual_i)
  // over the stand, with the per-cohort abundance n_i consumed as transported mass
  // exp(lambda_i) on the chart (never the reconstructed density, which overflows at
  // a tiny-but-nonzero spacing, odelia#46). `query` skips cohorts below that height
  // (a plant shorter than it shades nothing there); 0 censuses the whole stand.
  // compute_competition is the self-shading member (psi = per-individual shade).
  template <class Psi> value_type census(Psi psi, double query = 0.0) const;
  void compute_rates(const environment_type& environment, double pr_patch_survival, double birth_rate);
  // Transport-log-mass chart: reconstruct the log_density/density view from the
  // transported lambda (before any density is read), and seed a newborn's lambda
  // from its birth log_density. No-ops off the chart.
  void reconstruct_densities();
  void seed_newborn_log_mass();
  std::vector<value_type> net_reproduction_ratio_by_node() const;
  // Per-node lifetime offspring, weighted by patch-age density and S_D.
  std::vector<value_type> net_reproduction_ratio_by_node_weighted() const;
  // Introduction times of each node (the integration x-axis for fitness).
  std::vector<double> node_times() const;

  // * ODE interface
  // NOTE: We are a time-independent model here so no need to pass
  // time in as an argument.  All the bits involving time are taken
  // care of by Environment for us.
  // (ode_size/set_ode_state/ode_state/ode_rates come from SpeciesBase.)
  size_t aux_size() const;

  void resize_consumption_rates(int i);
  // Stand-total resource consumption for ODE channel i (e.g. TF24 soil-layer
  // water uptake). Carries value_type so the resident soil coupling
  // differentiates; the *_by_node_rev sibling stays double for the R interface.
  value_type consumption_rate(int i) const;
  std::vector<double> consumption_rate_by_node_rev(int i) const;

  template <typename It> It ode_aux(It it) const;

  Rcpp::NumericMatrix r_get_state() const;

  // * R interface
  std::vector<double> r_heights() const;
  std::vector<double> r_heights_rev() const;
  void r_set_heights(std::vector<double> heights);
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

  // The low-level strategy parameters the gradient is taken with respect to
  // (§8.1). All cohorts of this species share one strategy, so a species
  // contributes a single parameter block regardless of node count.
  std::vector<typename base_type::value_type*> ad_parameters() {
    return strategy->field_ptrs();
  }

  // Re-derive the strategy's precomputed quantities (canopy shape, birth size,
  // eta_c) from its current parameters. The gradient driver seeds ad_parameters()
  // and then reset()s, so this must run after the seed or a parameter whose effect
  // is mediated by a precomputed quantity (initial_height_, area_leaf_0, eta_c)
  // loses that part of its derivative -- the same reason IndividualRunner::reset()
  // re-prepares. All cohorts share one strategy, so one call covers the species.
  void prepare_strategy() { strategy->prepare_strategy(); }

private:
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
  nodes.clear();
  // Reset the new_node to a blank new_node, too.
  new_node = node_type(strategy);
}

template <typename T, typename E>
void Species<T,E>::introduce_new_node() {
  // new_node already holds the initial conditions computed against the current
  // environment by the most recent compute_rates() call (see compute_rates ->
  // new_node.compute_initial_conditions above), and the member is refreshed
  // again on the next compute_rates() ready for the following introduction.
  // Recomputing it here would be redundant, and would (wrongly) re-seed against
  // the post-introduction environment rather than the environment at the
  // node's introduction time (resolves the recompute question in #478).
  nodes.push_back(new_node);
  seed_newborn_log_mass();
}

// If a species contains no individuals, we return the height of a
// seed of the species.  Otherwise we return the height of the largest
// individual (always the first in the list) which will be at least
// tall as a seed.
template <typename T, typename E>
double Species<T,E>::height_max() const {
  // The stand's field-domain top: a frozen double L2 position (§5.5), never a
  // differentiated quantity, so narrow the tallest node's height.
  return xad::value(nodes.empty() ? new_node.height() : nodes.front().height());
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
// The population census Sum_i n_i * psi(individual_i), reduced as the same
// trapezium-in-mass every field reduction uses. On the chart the per-cohort
// abundance n_i is the bounded transported mass exp(lambda_i), consumed directly:
// a trapezium edge (h1-h0)*(d_i*psi_i + d_j*psi_j) is split into node-attributed
// halves ((h1-h0)/dx_i)*exp(lambda_i)*psi_i -- a moderate spacing ratio times a
// bounded mass -- so the reconstructed density exp(lambda)/dx (which overflows at
// a tiny-but-nonzero spacing near a growth stall, odelia#46) is never formed. dx_i
// is the SAME odelia::cohort_spacing that defines lambda, so /dx cancels the
// chart's *dx exactly on the interior. The new_node boundary carries no lambda and
// keeps its finite birth density; a coincident cohort (dx == 0) contributes 0
// (odelia's -inf convention). Off the chart log_density is the state, density is
// finite, and the plain density*psi trapezium is used. `query` is the self-shading
// early-exit height (a plant shorter than it shades nothing above it); 0 for a
// whole-stand census. Compare the design's Species::census<Psi>.
template <typename T, typename E>
template <class Psi>
typename Species<T,E>::value_type
Species<T,E>::census(Psi psi, double query) const {
  if (size() == 0) return value_type(0.0);

  const bool on_chart = nodes.front().on_mass_chart();
  std::vector<value_type> dx;
  if (on_chart) {
    std::vector<value_type> h(nodes.size());
    for (std::size_t i = 0; i < nodes.size(); ++i) h[i] = nodes[i].height();
    dx = odelia::cohort_spacing(h);
  }
  // One endpoint's contribution to an edge of width `gap`, before the final /2.
  // On the chart (and not the lambda-less new_node, idx == npos): the overflow-free
  // ((gap/dx_i)*exp(lambda_i))*psi_i. Otherwise the plain density*psi the edge used.
  const std::size_t npos = static_cast<std::size_t>(-1);
  auto edge_end = [&](std::size_t i, const value_type& gap,
                      const node_type& node) -> value_type {
    if (on_chart && i != npos) {
      if (!(dx[i] > 0.0)) return value_type(0.0);  // coincident cohort: zero mass
      return (gap / dx[i]) * exp(node.get_log_mass()) * psi(node.individual);
    }
    return gap * node.get_density() * psi(node.individual);
  };

  value_type tot = 0.0;
  std::size_t i_prev = 0;
  value_type h1 = nodes[0].height();
  bool broke = false;
  for (std::size_t i = 1; i < nodes.size(); ++i) {
    const value_type h0 = nodes[i].height();
    const value_type gap = h1 - h0;
    const value_type e =
        edge_end(i_prev, gap, nodes[i_prev]) + edge_end(i, gap, nodes[i]);
    if (!util::is_finite(e)) {
      util::stop("Detected non-finite contribution");
    }
    tot += e;
    h1 = h0;
    i_prev = i;
    if (h0 < query) { broke = true; break; }
  }

  // Include the new_node bottom boundary exactly when the old edge loop did: a
  // single node (it is the trapezium's second point), or the loop reached the
  // stand's base without early-exit and the last node still contributes. f_last
  // mirrors the old `f_h1 > 0` test on the bounded mass, so it agrees on sign.
  value_type mass_last;
  if (on_chart) mass_last = exp(nodes[i_prev].get_log_mass());
  else          mass_last = nodes[i_prev].get_density();
  const value_type f_last = mass_last * psi(nodes[i_prev].individual);
  if (size() == 1 || (!broke && f_last > 0)) {
    const value_type gap = h1 - new_node.height();
    tot += edge_end(i_prev, gap, nodes[i_prev]) + edge_end(npos, gap, new_node);
  }

  return tot / 2;
}

// The resident self-shading integral: a census whose per-individual weight is the
// shade that individual casts at `height`, with the same query-height early-exit.
template <typename T, typename E>
typename Species<T,E>::value_type
Species<T,E>::compute_competition(double height) const {
  // The query height and the domain guard stay double (height_max is a frozen L2
  // domain position, §5.5); the node heights and per-node shade carry value_type,
  // so a trait re-shades the stand on a gradient pass.
  if (size() == 0 || height_max() < height) {
    return value_type(0.0);
  }
  return census(
      [height](const individual_type& ind) { return ind.compute_competition(height); },
      height);
}

// NOTE: We should probably prefer to rescale when this is called
// through the ode stepper.
template <typename T, typename E>
void Species<T,E>::compute_rates(const E& environment, double pr_patch_survival, double birth_rate) {
  for (auto& c : nodes) {
    c.compute_rates(environment, pr_patch_survival);
  }
  new_node.compute_initial_conditions(environment, pr_patch_survival, birth_rate);

  // Transport-log-mass chart (geometric strategies with the flag on): the
  // transported quantity is lambda = log_density + log(cohort_spacing) and its
  // rate is simply -mortality -- the compression term -d(g)/d(size) cancels
  // against the log-spacing's own evolution and is never formed (no numerical
  // d(g)/d(size) on the tape, stable through a growth stall). Node::compute_rates
  // has already set log_density_dt = -mortality for the geometric branch, so
  // there is nothing to add here; log_density/density are reconstructed as a
  // read-side view (reconstruct_densities) before the field is assembled.
}

// Reconstruct the log_density/density view from the transported log mass and the
// canonical cohort spacing (odelia). Called once per environment build, before
// any density is read. A lone cohort has no spacing, so lambda == log_density.
template <typename T, typename E>
void Species<T,E>::reconstruct_densities() {
  const std::size_t n = nodes.size();
  if (n == 0) return;
  if (!nodes[0].on_mass_chart()) return;  // off the chart: log_density is the state
  if (n == 1) { nodes[0].set_log_density(nodes[0].get_log_mass()); return; }
  std::vector<value_type> heights(n), log_mass(n);
  for (std::size_t i = 0; i < n; ++i) {
    heights[i]  = nodes[i].height();
    log_mass[i] = nodes[i].get_log_mass();
  }
  const std::vector<value_type> ld =
      odelia::log_density_from_log_mass(heights, log_mass);
  for (std::size_t i = 0; i < n; ++i) nodes[i].set_log_density(ld[i]);
}

// Rebuild the transported log mass from the current log_density view whenever a
// cohort is introduced. Introducing a node changes its neighbours' spacing, so
// every cell's mass (lambda = log_density + log(cohort_spacing)) is re-derived
// from the density it holds now with the new spacing -- a remesh that preserves
// the physical density across the discretisation change. For cells whose spacing
// is unchanged this is the exact round-trip identity; for the newborn (and a
// lone cohort promoted to a pair) it sets the correct initial lambda. Must run at
// introduction, before the stepper reads the initial condition off ode_state.
template <typename T, typename E>
void Species<T,E>::seed_newborn_log_mass() {
  const std::size_t n = nodes.size();
  if (n == 0 || !nodes[0].on_mass_chart()) return;
  if (n == 1) { nodes[0].set_log_mass(nodes[0].get_log_density()); return; }
  std::vector<value_type> heights(n), log_density(n);
  for (std::size_t i = 0; i < n; ++i) {
    heights[i]     = nodes[i].height();
    log_density[i] = nodes[i].get_log_density();
  }
  const std::vector<value_type> lm =
      odelia::log_mass_from_log_density(heights, log_density);
  for (std::size_t i = 0; i < n; ++i) nodes[i].set_log_mass(lm[i]);
}

template <typename T, typename E>
void Species<T,E>::introduce_new_node(double time, double patch_density) {
  // Stamp the pushed copy (not new_node) so the member stays pristine for
  // the no-arg introduction paths.
  nodes.push_back(new_node);
  nodes.back().set_introduction(time, patch_density);
  // Seed the newborn's transported lambda from its birth log_density now that it
  // has a neighbour (needed before the stepper reads its initial condition).
  seed_newborn_log_mass();
}

template <typename T, typename E>
std::vector<typename Species<T,E>::value_type>
Species<T,E>::net_reproduction_ratio_by_node() const {
  std::vector<value_type> ret;
  ret.reserve(size());
  for (auto& c : nodes) {
    ret.push_back(c.fecundity());
  }
  return ret;
}

template <typename T, typename E>
std::vector<typename Species<T,E>::value_type>
Species<T,E>::net_reproduction_ratio_by_node_weighted() const {
  std::vector<value_type> ret;
  ret.reserve(size());
  for (auto& c : nodes) {
    ret.push_back(c.weighted_fecundity(strategy->pars.S_D));
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
typename Species<T,E>::value_type Species<T,E>::consumption_rate(int i) const {
  // can't determine density for one node
  if(size() < 2) {
    return 0.0;
  } else {
    // node heights are in descending order - we need ascending for integration.
    // Trapezium rule accumulated in value_type (util::trapezium narrows to
    // double); the derivative w.r.t. traits flows through the node heights and
    // per-node consumption rates. Reproduces the double path exactly.
    std::vector<value_type> h, c;
    h.reserve(size());
    c.reserve(size());
    for (auto it = nodes.rbegin(); it != nodes.rend(); ++it) {
      h.push_back(it->height());
      c.push_back(it->consumption_rate(i));
    }
    value_type tot = 0.0;
    for (std::size_t k = 1; k < h.size(); ++k) {
      tot += (h[k] - h[k - 1]) * (c[k] + c[k - 1]);
    }
    return tot * 0.5;
  }
}

template <typename T, typename E>
std::vector<double> Species<T,E>::consumption_rate_by_node_rev(int i) const {
  std::vector<double> ret;
  ret.reserve(size());
  for(auto it = nodes.rbegin(); it != nodes.rend(); ++it) {
    ret.push_back(odelia::util::diagnostic(it->consumption_rate(i)));  // R-facing: double only
  }
  return ret;
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
    ret.push_back(odelia::util::diagnostic(it->height()));  // R-facing: double only
  }
  return ret;
}

template <typename T, typename E>
std::vector<double> Species<T,E>::r_heights_rev() const {
  std::vector<double> ret;
  ret.reserve(size());
  for (nodes_const_iterator it = nodes.begin();
       it != nodes.end(); ++it) {
    ret.push_back(odelia::util::diagnostic(it->height()));  // R-facing: double only
  }
  std::reverse(ret.begin(), ret.end());
  return ret;
}

template <typename T, typename E>
void Species<T,E>::r_set_heights(std::vector<double> heights) {
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
    ret.push_back(odelia::util::diagnostic(c.compute_competition(0.0)));  // R-facing: double
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
    ret.push_back(odelia::util::diagnostic(it->get_log_density()));  // R-facing: double only
  }
  return ret;
}

template <typename T, typename E>
std::vector<double> Species<T,E>::r_log_density_rates() const {
  std::vector<double> ret;
  ret.reserve(size());
  for (nodes_const_iterator it = nodes.begin(); it != nodes.end(); ++it) {
    // R-facing diagnostic (only double crosses to R): narrow the density rate,
    // which is value_type on a gradient pass. check_initial_density_rates reads
    // it as a threshold guard -- off the differentiated value.
    ret.push_back(odelia::util::diagnostic(it->get_log_density_rate()));
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
