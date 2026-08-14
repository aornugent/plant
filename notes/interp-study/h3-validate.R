## Validate the harness before measuring anything with it:
##  1. the R Hermite reproduces the shipped interpolant at arbitrary heights;
##  2. exact_field reproduces the shipped build's knot values and slopes.
## Resolve siblings relative to this file, so the study runs from the repo.
STUDY <- tryCatch(dirname(normalizePath(sys.frame(1)$ofile)), error = function(e) ".")
if (!file.exists(file.path(STUDY, "lib-field.R"))) STUDY <- "notes/interp-study"
source(file.path(STUDY, "lib-field.R"))

p <- ladder_patch_two_by_two()
p$compute_environment()

k <- ladder_field_knots_tf24(p)
f <- herm(k$height, k$value, k$slope)

hmax <- max(k$height)
eta <- ladder_strategy_parameter(p, 1, "eta")
cat("eta =", eta, "\n")

## 1. the R Hermite against the shipped one, away from the knots
z <- seq(1e-6, hmax * 0.999, length.out = 997)
cpp <- vapply(z, function(zz) p$environment$get_environment_at_height(zz), numeric(1))
mine <- herm_eval(f, z)$value
cat(sprintf("R Hermite vs shipped interpolant : max abs diff %.3e\n",
            max(abs(mine - cpp))))

## 2. exact_field against the shipped build's knots
ex <- exact_field(p, k$height, "L")
cat(sprintf("exact_field vs shipped knot value: max abs diff %.3e\n",
            max(abs(ex$value - k$value))))
cat(sprintf("exact_field vs shipped knot slope: max abs diff %.3e\n",
            max(abs(ex$slope - k$slope))))

## 3. the interpolant is exact at its own knots
at <- herm_eval(f, k$height)
cat(sprintf("interpolant at its knots         : value %.3e  slope %.3e\n",
            max(abs(at$value - k$value)), max(abs(at$slope - k$slope))))

## 4. the ground slope really is zero, so the flat first span is a fact
cat(sprintf("A'(0) = %.3e   L'(0) = %.3e\n",
            p$compute_competition_and_slope(0)[[2]], ex$slope[[1]]))
