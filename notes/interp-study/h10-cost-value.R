## Cost at a canopy tall enough to stress the fixed grid's knot count, and what
## switching placement costs in forward values -- which is the re-blessing.
## Resolve siblings relative to this file, so the study runs from the repo.
STUDY <- tryCatch(dirname(normalizePath(sys.frame(1)$ofile)), error = function(e) ".")
if (!file.exists(file.path(STUDY, "lib-field.R"))) STUDY <- "notes/interp-study"
source(file.path(STUDY, "lib-field.R"))

bench <- function(f, reps) { f(); system.time(for (i in seq_len(reps)) f())[["elapsed"]] / reps }

run_long <- function(lifetime = 20) {
  p <- ladder_parameters(c("fast", "slow"), lifetime = lifetime)
  p$node_schedule_times <- list(seq(0, lifetime * 0.9, length.out = 20),
                                seq(0, lifetime * 0.9, length.out = 20))
  ladder_run(p)
}

settings <- list(
  list("canopy", 0.10,  65, 1),
  list("canopy", 0.10, 129, 1),
  list("fixed",  0.50,  65, 1),
  list("fixed",  0.25,  65, 1),
  list("fixed",  0.10,  65, 1))

cat("=== a 20-year run: canopy height, knots, cost, and the census ===\n")
cat(sprintf("%-20s %7s %6s %11s %8s %14s %14s\n",
            "policy", "h_max", "knots", "sec/run", "rel", "leaf area", "h_max chk"))
base <- NA; ref <- NULL
for (s in settings) {
  do.call(interp_policy_set, s)
  st <- try(run_long(), silent = TRUE)
  if (inherits(st, "try-error")) {
    cat(sprintf("%-20s   FAILED: %s\n", s[[1]], conditionMessage(attr(st, "condition"))))
    next
  }
  pa <- ladder_as_patch(st)
  hs <- unlist(lapply(pa$species, function(sp)
    vapply(sp$nodes, function(n) n$height, numeric(1))))
  nk <- length(ladder_field_knots_tf24(pa)$height)
  la <- ladder_census_leaf_area(st)
  off <- max(hs)
  tm <- bench(function() run_long(), 3)
  if (is.na(base)) { base <- tm; ref <- c(la, off) }
  cat(sprintf("%-20s %7.2f %6d %11.3e %7.2fx %14.7e %14.7e\n",
              sprintf("%s d=%.2f n=%d", s[[1]], s[[2]], s[[3]]),
              max(hs), nk, tm, tm / base, la, off))
}
cat(sprintf("\nreference (first row): leaf area %.7e  offspring %.7e\n", ref[1], ref[2]))

## ---- what the forward answer moves by --------------------------------------
cat("\n=== forward value shift, relative to the current placement ===\n")
interp_policy_set("canopy", 0.10, 65, 1)
top_of <- function(st) {
  pa <- ladder_as_patch(st)
  max(unlist(lapply(pa$species, function(sp)
    vapply(sp$nodes, function(n) n$height, numeric(1)))))
}
st0 <- run_long()
la0 <- ladder_census_leaf_area(st0); h0 <- top_of(st0)
for (s in list(list("canopy", 0.10, 257, 1), list("fixed", 0.50, 65, 1),
               list("fixed", 0.25, 65, 1), list("fixed", 0.10, 65, 1))) {
  do.call(interp_policy_set, s)
  st <- run_long()
  cat(sprintf("%-20s leaf area %+.3e   canopy top %+.3e\n",
              sprintf("%s d=%.2f n=%d", s[[1]], s[[2]], s[[3]]),
              ladder_census_leaf_area(st) / la0 - 1, top_of(st) / h0 - 1))
}
