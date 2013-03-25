// -*-c++-*-
#ifndef TREE_UTIL_H_
#define TREE_UTIL_H_

#include <Rcpp.h>
#include <gsl/gsl_nan.h>

namespace util {

template <typename T>
bool is_finite(T x) {
  // TODO: Get the finite check in here!
  // Rf_warning("Requesting finite check, but not yet implemented");
  return true;
}

void check_bounds(int idx, int max);

// Based on C++11's is_sorted
template <class ForwardIterator> 
bool is_decreasing(ForwardIterator first, ForwardIterator last) {
  if ( first == last ) 
    return true;

  ForwardIterator next = first;
  while ( ++next != last ) {
    if ( *next > *first )
      return false;
    ++first;
  }
  return true;
}


}

#endif
