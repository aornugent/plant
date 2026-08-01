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
// 19 errors, all of them the sites the design refuses: 17 in the DeepCrown
// branch, which TF24 does not default to, and 2 named static assertions inside
// prepare_strategy(), which never runs inside a recorded block.

#include <plant/models/tf24_strategy.h>
#include <plant/individual_runner.h>
#include <odelia/ode_solver.hpp>

// The active scalar comes from odelia's own alias, on a Solver plant already
// instantiates, so plant keeps naming the AD library to odelia.
using tf24_runner =
    plant::tools::IndividualRunner<plant::TF24_Strategy<double>,
                                   plant::TF24_Environment<double> >;
using active_scalar = odelia::ode::Solver<tf24_runner>::active_scalar;

// A sizeof assertion passes on a class body alone, so what is checked here is
// the strategy's member bodies.
static_assert(sizeof(plant::TF24_Environment<active_scalar>) > 0);
static_assert(sizeof(plant::TF24_Strategy<active_scalar>) > 0);

template class plant::TF24_Strategy<active_scalar>;

// Instantiating plant::Individual on top of this needs its state store to carry
// the scalar as well: it holds Internals<double>, so every place it hands vars to
// the strategy is a passive seam. So the strategy is what this probe claims.
