// Generated by RcppR6 (0.2.4): do not edit by hand
#ifndef _PLANT_RCPPR6_PRE_HPP_
#define _PLANT_RCPPR6_PRE_HPP_

#include <RcppCommon.h>


namespace plant {
namespace RcppR6 {
template <typename T> class RcppR6;
}
}

namespace plant { namespace ode { namespace test { class OdeR; } } }

namespace Rcpp {
template <typename T> SEXP wrap(const plant::RcppR6::RcppR6<T>&);
namespace traits {
template <typename T> class Exporter<plant::RcppR6::RcppR6<T> >;
}

template <> SEXP wrap(const plant::ode::test::Lorenz&);
template <> plant::ode::test::Lorenz as(SEXP);
template <> SEXP wrap(const plant::ode::test::OdeR&);
template <> plant::ode::test::OdeR as(SEXP);
template <> SEXP wrap(const plant::ode::Runner<plant::ode::test::Lorenz>&);
template <> plant::ode::Runner<plant::ode::test::Lorenz> as(SEXP);

template <> SEXP wrap(const plant::ode::Runner<plant::ode::test::OdeR>&);
template <> plant::ode::Runner<plant::ode::test::OdeR> as(SEXP);

template <> SEXP wrap(const plant::ode::Runner<plant::tools::PlantRunner<plant::FF16_Strategy, plant::FF16_Environment> >&);
template <> plant::ode::Runner<plant::tools::PlantRunner<plant::FF16_Strategy, plant::FF16_Environment> > as(SEXP);

template <> SEXP wrap(const plant::ode::Runner<plant::tools::PlantRunner<plant::FF16r_Strategy, plant::FF16r_Environment> >&);
template <> plant::ode::Runner<plant::tools::PlantRunner<plant::FF16r_Strategy, plant::FF16r_Environment> > as(SEXP);
template <> SEXP wrap(const plant::CohortScheduleEvent&);
template <> plant::CohortScheduleEvent as(SEXP);
template <> SEXP wrap(const plant::CohortSchedule&);
template <> plant::CohortSchedule as(SEXP);
template <> SEXP wrap(const plant::Disturbance&);
template <> plant::Disturbance as(SEXP);
template <> SEXP wrap(const plant::Control&);
template <> plant::Control as(SEXP);
template <> SEXP wrap(const plant::ode::OdeControl&);
template <> plant::ode::OdeControl as(SEXP);
template <> SEXP wrap(const plant::quadrature::QK&);
template <> plant::quadrature::QK as(SEXP);
template <> SEXP wrap(const plant::quadrature::QAG&);
template <> plant::quadrature::QAG as(SEXP);
template <> SEXP wrap(const plant::interpolator::Interpolator&);
template <> plant::interpolator::Interpolator as(SEXP);
template <> SEXP wrap(const plant::Plant<plant::FF16_Strategy,plant::FF16_Environment>&);
template <> plant::Plant<plant::FF16_Strategy,plant::FF16_Environment> as(SEXP);

template <> SEXP wrap(const plant::Plant<plant::FF16r_Strategy,plant::FF16r_Environment>&);
template <> plant::Plant<plant::FF16r_Strategy,plant::FF16r_Environment> as(SEXP);
template <> SEXP wrap(const plant::tools::PlantRunner<plant::FF16_Strategy,plant::FF16_Environment>&);
template <> plant::tools::PlantRunner<plant::FF16_Strategy,plant::FF16_Environment> as(SEXP);

template <> SEXP wrap(const plant::tools::PlantRunner<plant::FF16r_Strategy,plant::FF16r_Environment>&);
template <> plant::tools::PlantRunner<plant::FF16r_Strategy,plant::FF16r_Environment> as(SEXP);
template <> SEXP wrap(const plant::Internals&);
template <> plant::Internals as(SEXP);
template <> SEXP wrap(const plant::Parameters<plant::FF16_Strategy,plant::FF16_Environment>&);
template <> plant::Parameters<plant::FF16_Strategy,plant::FF16_Environment> as(SEXP);

template <> SEXP wrap(const plant::Parameters<plant::FF16r_Strategy,plant::FF16r_Environment>&);
template <> plant::Parameters<plant::FF16r_Strategy,plant::FF16r_Environment> as(SEXP);
template <> SEXP wrap(const plant::Cohort<plant::FF16_Strategy,plant::FF16_Environment>&);
template <> plant::Cohort<plant::FF16_Strategy,plant::FF16_Environment> as(SEXP);

template <> SEXP wrap(const plant::Cohort<plant::FF16r_Strategy,plant::FF16r_Environment>&);
template <> plant::Cohort<plant::FF16r_Strategy,plant::FF16r_Environment> as(SEXP);
template <> SEXP wrap(const plant::Species<plant::FF16_Strategy,plant::FF16_Environment>&);
template <> plant::Species<plant::FF16_Strategy,plant::FF16_Environment> as(SEXP);

template <> SEXP wrap(const plant::Species<plant::FF16r_Strategy,plant::FF16r_Environment>&);
template <> plant::Species<plant::FF16r_Strategy,plant::FF16r_Environment> as(SEXP);
template <> SEXP wrap(const plant::Patch<plant::FF16_Strategy,plant::FF16_Environment>&);
template <> plant::Patch<plant::FF16_Strategy,plant::FF16_Environment> as(SEXP);

template <> SEXP wrap(const plant::Patch<plant::FF16r_Strategy,plant::FF16r_Environment>&);
template <> plant::Patch<plant::FF16r_Strategy,plant::FF16r_Environment> as(SEXP);
template <> SEXP wrap(const plant::SCM<plant::FF16_Strategy,plant::FF16_Environment>&);
template <> plant::SCM<plant::FF16_Strategy,plant::FF16_Environment> as(SEXP);

template <> SEXP wrap(const plant::SCM<plant::FF16r_Strategy,plant::FF16r_Environment>&);
template <> plant::SCM<plant::FF16r_Strategy,plant::FF16r_Environment> as(SEXP);
template <> SEXP wrap(const plant::StochasticSpecies<plant::FF16_Strategy,plant::FF16_Environment>&);
template <> plant::StochasticSpecies<plant::FF16_Strategy,plant::FF16_Environment> as(SEXP);

template <> SEXP wrap(const plant::StochasticSpecies<plant::FF16r_Strategy,plant::FF16r_Environment>&);
template <> plant::StochasticSpecies<plant::FF16r_Strategy,plant::FF16r_Environment> as(SEXP);
template <> SEXP wrap(const plant::StochasticPatch<plant::FF16_Strategy,plant::FF16_Environment>&);
template <> plant::StochasticPatch<plant::FF16_Strategy,plant::FF16_Environment> as(SEXP);

template <> SEXP wrap(const plant::StochasticPatch<plant::FF16r_Strategy,plant::FF16r_Environment>&);
template <> plant::StochasticPatch<plant::FF16r_Strategy,plant::FF16r_Environment> as(SEXP);
template <> SEXP wrap(const plant::StochasticPatchRunner<plant::FF16_Strategy,plant::FF16_Environment>&);
template <> plant::StochasticPatchRunner<plant::FF16_Strategy,plant::FF16_Environment> as(SEXP);

template <> SEXP wrap(const plant::StochasticPatchRunner<plant::FF16r_Strategy,plant::FF16r_Environment>&);
template <> plant::StochasticPatchRunner<plant::FF16r_Strategy,plant::FF16r_Environment> as(SEXP);
template <> SEXP wrap(const plant::FF16_Strategy&);
template <> plant::FF16_Strategy as(SEXP);
template <> SEXP wrap(const plant::FF16_Environment&);
template <> plant::FF16_Environment as(SEXP);
template <> SEXP wrap(const plant::FF16r_Strategy&);
template <> plant::FF16r_Strategy as(SEXP);
template <> SEXP wrap(const plant::FF16r_Environment&);
template <> plant::FF16r_Environment as(SEXP);
}

#endif
