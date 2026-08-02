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
// 2 errors, both of them named static assertions inside prepare_strategy(),
// which never runs inside a recorded block; the DeepCrown branch is now
// refused at the active scalar with if constexpr rather than instantiated.
//
// Patch and StochasticPatch are instantiated below, and with them Species,
// Node, StochasticSpecies, StochasticNode and Individual, so those 2 are also
// the whole count: a new line in any other file is a passive seam. Instantiate
// the outermost consumer -- a census taken through Individual alone reports
// neither the containers that hold it nor the numerics they reach.
//
// What this leaves out: TF24 only, so FF16 and K93 are unmeasured; SCM and
// StochasticPatchRunner, which no line below names; and every member template
// a call from these instantiations does not reach.

#include <plant/models/tf24_strategy.h>
#include <plant/individual_runner.h>
#include <plant/node.h>
#include <plant/patch.h>
#include <plant/stochastic_patch.h>
#include <odelia/ode_solver.hpp>
#include <type_traits>
#include <vector>

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

using active_environment = plant::TF24_Environment<active_scalar>;
using active_node = plant::Node<plant::TF24_Strategy<active_scalar>,
                                active_environment>;

template class plant::Individual<plant::TF24_Strategy<active_scalar>,
                                 active_environment>;

// Node::growth_rate_gradient is private; r_growth_rate_gradient forwards to it
// and returns what it returns. A double on either declaration still compiles,
// because the value is taken through a return type written as double, so the
// transport term would read a derivative of zero with nothing raised. The
// assertion is the only check for that.
static_assert(std::is_same_v<
    decltype(std::declval<active_node&>().r_growth_rate_gradient(
               std::declval<const active_environment&>())),
    active_scalar>,
  "the transport probe returns the strategy's scalar");

template class plant::Patch<plant::TF24_Strategy<active_scalar>,
                            plant::TF24_Environment<active_scalar> >;

template class plant::StochasticPatch<plant::TF24_Strategy<active_scalar>,
                                      plant::TF24_Environment<active_scalar> >;

// An explicit class instantiation leaves member templates uninstantiated, so
// the serialisers need naming one by one at the double iterator R hands them.
using active_individual = plant::Individual<plant::TF24_Strategy<active_scalar>,
                                            active_environment>;
using state_iterator = std::vector<double>::iterator;
template state_iterator active_individual::ode_state(state_iterator) const;
template state_iterator active_individual::ode_rates(state_iterator) const;
template state_iterator active_individual::ode_aux(state_iterator) const;

// Instantiating the class does not instantiate its member templates, so the
// five the solver drives are called here. A census reaches only what it calls.
using active_patch = plant::Patch<plant::TF24_Strategy<active_scalar>,
                                  plant::TF24_Environment<active_scalar> >;

void solver_driven_members(active_patch& patch,
                           std::vector<active_scalar>& y) {
  patch.set_ode_state(y.begin(), 0.0);
  patch.ode_state(y.begin());
  patch.ode_rates(y.begin());
  patch.ode_aux(y.begin());
  patch.set_ode_aux(y.begin());
}

// rebind_from is a member template, so the class instantiation above does not
// reach it; name the double->active crossing the adjoint stepper takes.
template plant::Patch<plant::TF24_Strategy<active_scalar>,
                      plant::TF24_Environment<active_scalar> >
plant::Patch<plant::TF24_Strategy<double>,
             plant::TF24_Environment<double> >::rebind_from<active_scalar>() const;
