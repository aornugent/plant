// -*-c++-*-
#ifndef PLANT_PLANT_COUPLING_H_
#define PLANT_PLANT_COUPLING_H_

#include <concepts>
#include <vector>
#include <cstddef>

namespace plant {

// The contract Patch<T,E> requires of its environment E. Patch static_asserts
// Coupling<environment_type>, so E is checked *at the scalar Patch is
// instantiated with* -- on the production path that is E::value_type = double;
// on a gradient pass it is the active scalar, checked against the active
// value_type iterator. A conforming E therefore cannot then fail deep in the
// templated <It> pass (plan §6.1 [S4]).
//
// This is the pre-gradient form, keyed to the environment interface as it stands
// (the physiology point-read, the ODE-state size, the templated
// ODE-serialisation seam on the value_type iterator, and the environment's own
// state-rate hook). It deliberately does NOT yet require the recording-era
// surface (at(Query) / set_field_active / light_knots() / compute_state_rates /
// rebind_from<S>, plan §6.1): EnvironmentRecording and the active replay are
// deferred, and the raw record/replay members still live on Patch. Tighten this
// toward the §6.1 shape when that machinery lands.
template <class E>
concept Coupling =
    requires(E e, const E ce,
             typename E::value_type q,
             const std::vector<double>& depletion,
             typename std::vector<typename E::value_type>::iterator it,
             typename std::vector<typename E::value_type>::const_iterator cit) {
  typename E::value_type;

  // The physiology point-read: light at an (active) height, returning value_type.
  { ce.get_environment_at_height(q) } -> std::same_as<typename E::value_type>;

  // ODE-state size, and the ODE-serialisation seam templated on the value_type
  // iterator (NOT odelia's double-only free helpers) -- plan §3/§6.1. FF16/K93
  // carry no environment ODE state (ode_size() == 0), so these are
  // pass-throughs; TF24's soil state rides them on the double alias.
  { ce.ode_size() } -> std::convertible_to<std::size_t>;
  { e.set_ode_state(cit) };
  { ce.ode_state(it) };
  { ce.ode_rates(it) };

  // The environment's own state rates (soil water for TF24; a no-op for the
  // light-only FF16/K93 environments).
  e.compute_rates(depletion);
};

} // namespace plant

#endif
