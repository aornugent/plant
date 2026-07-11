#include <plant.h>

// The FF16 Patch must satisfy odelia's Replayable concept: it is what wires the
// stepper's record/replay hooks and the frozen-environment dispatch. A rename
// that broke a hook would silently drop mutant runs back to recompute; assert
// it at compile time on the concrete FF16 Patch instead.
static_assert(
    odelia::ode::Replayable<plant::Patch<plant::FF16_Strategy,
                                         plant::FF16_Environment>>,
    "Patch must satisfy odelia::ode::Replayable (record_stage / record_ode_step "
    "/ replay_step / has_recorded_field)");

// Helpers for FF16 model

// [[Rcpp::export]]
plant::NodeSchedule node_schedule_default__Parameters___FF16__FF16_Env(const plant::Parameters<plant::FF16_Strategy,plant::FF16_Environment>& p) {
   return plant::node_schedule_default<plant::Parameters<plant::FF16_Strategy,plant::FF16_Environment> >(p);
}

// [[Rcpp::export]]
plant::NodeSchedule make_node_schedule__Parameters___FF16__FF16_Env(const plant::Parameters<plant::FF16_Strategy,plant::FF16_Environment>& p) {
   return plant::make_node_schedule<plant::Parameters<plant::FF16_Strategy,plant::FF16_Environment> >(p);
}
