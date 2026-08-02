# Building and running the gate harnesses

`leaf_jac_gate.cpp` is a standalone program against `plant/leaf_model.h`. It needs
`Leaf`'s own translation unit and the quadrature and `util::stop` ones it calls, and
`util::stop` reaches R, so the link takes `-lR`:

```sh
cd <plant worktree>
RCPP=$(Rscript -e 'cat(system.file("include", package = "Rcpp"))')
BH=$(Rscript -e 'cat(system.file("include", package = "BH"))')
RINC=$(Rscript -e 'cat(R.home("include"))')
RLIB=$(Rscript -e 'cat(R.home("lib"))')
g++ -std=gnu++20 -O2 -DNDEBUG -I inst/include \
  -isystem <odelia library>/odelia/include -isystem "$RCPP" -isystem "$BH" -isystem "$RINC" \
  scratch/leaf_jac_gate.cpp src/leaf_model.cpp src/qk.cpp src/qk_rules.cpp \
  src/qag.cpp src/qag_internals.cpp src/util.cpp -L"$RLIB" -lR -o leaf_jac_gate
./leaf_jac_gate
```

It prints one block per hand-built soil state: the collar residual and curvature, the
continuity, waist, translation and stationarity residuals, the `inputs()`/`input_adjoints`
size identity, and the `FULLSOLVE` comparisons. The third state is collar-pinned, so its
stationarity rows read order 1 and its `FULLSOLVE` rows read `nan`, which is the unbuilt
`d(bound)/du`. The link is about 17 s and the run is under 0.01 s.

`wire_gates.cpp` is `Rcpp::sourceCpp`'d and its exports take a `Patch` R6 object through
`SEXP`, so it needs plant loaded and an R driver that hands it a patch. The exports resolve
plant's own compiled symbols only if the build links the package's `.so`, which
`sourceCpp` will not do on its own:

```r
library(odelia); pkgload::load_all("<plant worktree>")
Sys.setenv(PKG_CPPFLAGS = "-I<plant worktree>/inst/include -DNDEBUG",
           PKG_LIBS = "<plant worktree>/src/plant.so")
Rcpp::sourceCpp("<plant worktree>/scratch/wire_gates.cpp")
```

The driver then builds TF24 parameters, runs an `SCM` so the patch is at a state the model
visits, takes `scm$patch`, and calls the exports with an output-adjoint seed:
`patch_adjoint_partial(patch, lam, upto)` against `patch_recording_vjp(patch, lam)` for the
decomposition, `block_vjp` against a central difference of `block_value` for one cohort's
block, and `colSums(knot_contributions(patch, lam))` against the accumulated
`knot_value` for the knots. `block_vjp`'s output-adjoint argument must be at least as long
as the block's output count (11 for TF24); a shorter one overruns and corrupts the heap.
A driver at `max_patch_lifetime = 10`, compile included, is about 45 s, dominated by the
per-output block sweeps.
