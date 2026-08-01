// Compile-only reproducer for TF24 at the adjoint active scalar: compilation is
// the whole check, and what it reports is the list of sites the active build still
// has to reach. Run it by hand; nothing runs it for you, so the count below is a
// record of one run rather than a gate.
//
//   x86_64-linux-gnu-g++ -std=gnu++20 -fsyntax-only -fmax-errors=200 \
//     -I$(R.home include) -I$(Rcpp include) -I$(BH include) \
//     -I$(odelia include) -isystem inst/include -DNDEBUG \
//     scripts/tf24-active-probe.cpp
//
// At plant p1/phase-1 against odelia p1/odelia-integration: 41 errors at 33 sites,
// in the six groups implementation-notes.md records.

#include <plant/models/tf24_strategy.h>
#include <plant/individual_runner.h>
#include <odelia/ode_solver.hpp>

// The active scalar comes from odelia's own alias, on a Solver plant already
// instantiates, so plant keeps naming the AD library to odelia.
using tf24_runner =
    plant::tools::IndividualRunner<plant::TF24_Strategy<double>,
                                   plant::TF24_Environment<double> >;
using active_scalar = odelia::ode::Solver<tf24_runner>::active_scalar;

// The class bodies instantiate; the member bodies below are what still fails.
static_assert(sizeof(plant::TF24_Environment<active_scalar>) > 0);
static_assert(sizeof(plant::TF24_Strategy<active_scalar>) > 0);

template class plant::TF24_Strategy<active_scalar>;
