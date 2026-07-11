#include <plant.h>

// Patch satisfies odelia's Replayable concept, so derivs compiles the frozen-field
// dispatch the mutant replay needs (record/replay hooks + has_recorded_field).
static_assert(odelia::ode::Replayable<
              plant::Patch<plant::FF16_Strategy, plant::FF16_Environment>>);

// Helpers for FF16 model

// [[Rcpp::export]]
plant::NodeSchedule node_schedule_default__Parameters___FF16__FF16_Env(const plant::Parameters<plant::FF16_Strategy,plant::FF16_Environment>& p) {
   return plant::node_schedule_default<plant::Parameters<plant::FF16_Strategy,plant::FF16_Environment> >(p);
}

// [[Rcpp::export]]
plant::NodeSchedule make_node_schedule__Parameters___FF16__FF16_Env(const plant::Parameters<plant::FF16_Strategy,plant::FF16_Environment>& p) {
   return plant::make_node_schedule<plant::Parameters<plant::FF16_Strategy,plant::FF16_Environment> >(p);
}
