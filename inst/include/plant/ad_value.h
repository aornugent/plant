// -*-c++-*-
#ifndef PLANT_PLANT_AD_VALUE_H_
#define PLANT_PLANT_AD_VALUE_H_

#include <type_traits>

namespace plant {

// Passive value of a model scalar at an active->double boundary. For the resident
// numerics (S = double) it is the identity, so the double path is bit-unchanged.
// For an active scalar it returns the underlying double via xad::value, resolved
// by ADL in the translation unit that instantiates the active model (XAD is
// included there, not in the double build). Use only where a genuine double is
// required -- competition heights feeding the double light spline, cohort
// densities and survival, the finite-difference growth gradient -- never to store
// a physiological value that must keep its derivative.
template <class S>
double ad_value(const S& x) {
  if constexpr (std::is_same_v<S, double>) {
    return x;
  } else {
    return value(x);
  }
}

}

#endif
