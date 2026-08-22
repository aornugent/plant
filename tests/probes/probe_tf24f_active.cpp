// Report 07 §12's last falsifier, fired at compile time.
//
// The tracked-collar variant needs no implicit solve, no stationarity condition
// and no curvature: its collar is an ODE state whose derivative arrives from the
// adjoint like any other. So it is the case that says whether this interface was
// carved around the argmax, and the symptom of a bad carve is that the EASY case
// cannot be written in it.
//
// Never run. Instantiating it is the whole test, and the errors are the
// specification.
#include <plant/models/tf24f_strategy.h>
#include <plant/models/tf24_environment.h>

using tangent = plant::tangent;

void probe_tf24f_at_an_active_scalar() {
  plant::TF24f_Strategy<tangent> s;
  plant::TF24_Environment<tangent> env;
  plant::Internals<tangent> vars(plant::TF24f_Strategy<tangent>::state_size());
  s.refresh_indices();
  s.compute_rates(env, vars);
}

// Compile it, do not link it:
//
//   g++ -std=c++20 -fsyntax-only -DXAD_NO_THREADLOCAL -DXAD_USE_STRONG_INLINE \
//       -I../../inst/include $(Rscript -e 'cat(paste0("-I", c(
//         system.file("include", package = "Rcpp"),
//         system.file("include", package = "BH"),
//         system.file("include", package = "odelia"),
//         system.file("include", package = "phylloptim"),
//         system.file("include", package = "RcppR6")), collapse = " "))') \
//       $(R CMD config --cppflags) probe_tf24f_active.cpp
//
// What it reports today, and the third one is the one to care about:
//
//   161: cannot convert xad::FReal<double,1> to double
//        -- the tracked collar is handed to evaluate_root_collar_psi, which takes
//        a double because the leaf is a double-only model.
//   194: no matching call to max(FReal&, double&)
//        -- the same, through the finite-difference branch's clamp.
//   AND NOTHING AT LINE 162, which is the finding. `dprofit_dpsi_` is S and
//        `dprofit_droot_collar_psi` returns double, so the assignment is a legal
//        conversion that drops the derivative. The tracked state's rate is
//        k_acclim * dprofit_dpsi_, so every trait row reaching the collar is a
//        structural zero -- finite, plausible, and silent.
