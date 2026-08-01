// Compile-only reproducer for TF24 at the adjoint active scalar. Compilation is
// the whole check, and it currently reports 41 errors; testthat never runs this.
//
//   x86_64-linux-gnu-g++ -std=gnu++20 -fsyntax-only -fmax-errors=200 \
//     -I$(R.home include) -I$(Rcpp include) -I$(BH include) \
//     -I$(odelia include) -isystem inst/include -DNDEBUG \
//     tests/testthat/tf24-active-probe.cpp

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
