// -*-c++-*-
#ifndef PLANT_PLANT_WITH_SLOPE_H_
#define PLANT_PLANT_WITH_SLOPE_H_

namespace plant {

// A quantity and its vertical slope, which is what every frame of the competition
// path carries: an individual's contribution, a species' sum over its nodes, a
// patch's sum over its species, and a light-environment knot.
//
// `slope` is the SIGNED d(value)/d(height), so it is negative wherever the profile
// falls with height. The canopy shape forms the pair the other way round -- q is
// exactly -dQ/dz -- and each strategy's compute_competition_and_slope negates it on
// the way out, which is the one place the convention turns over. Above there it is
// signed, which is what build_extinction_field's chain rule assumes.
//
// The shape's own pair stays a pair of Q and q: those are the model's names for two
// quantities it takes from one u^eta, not a value and a slope.
//
// No includes of its own, so the headers on this path can take it without taking
// R's headers with it.
template <typename T>
struct with_slope {
  T value;
  T slope;

  // ⚠️ WITHOUT THIS THE TAPE AUDIT SKIPS THE PAIR IN SILENCE. odelia's
  // visit_active passes over any shape it does not open, without refusing it, and
  // it does not open an aggregate of two scalars -- so this type has to say what
  // it holds or its two members leave the walk. Both are scalars, so the visitor
  // takes each directly. active_system::release is what reports a miss.
  template <class F>
  void for_each_active(F&& f) {
    f(value);
    f(slope);
  }
};

}

#endif
