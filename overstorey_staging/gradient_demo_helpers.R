## Shared driving and interpretation code for the reverse-mode gradient demos.
##
## Kept out of the .qmd so tests/testthat/test-gradient-demo.R exercises the same
## code the documents do. No plotting here.
##
## Assumes `plant` is loaded (library(plant) or pkgload::load_all()).
##
## Two of these exist because a raw gradient column is not the number an
## ecologist means, for two separate reasons, and every demo hits both:
##
##   gd_elasticity()        a partial derivative carries the units of its
##                          parameter, so d(mass)/d(omega) and d(mass)/d(lma) are
##                          not comparable and the larger one is usually just the
##                          one with the smaller units.
##   gd_hyperpar_jacobian() TF24_hyperpar derives k_l, r_l, nmass_l from lma;
##                          d_I, g1_TF24, k_s, r_s, r_b from rho; a_f3 from
##                          omega; and c, p_50, b, psi_crit from K_s. The
##                          gradient's `lma` column holds those fixed, so it is
##                          the partial and not the trait derivative.

## The trait point the demos vary around. Five of the parameters TF24 carries,
## chosen because trait_matrix + TF24_hyperpar drive them and a trait database
## reports them.
gd_traits <- function(...) {
  utils::modifyList(
    list(lma = 0.0825, hmat = 5.13, k_I = 0.5, a_l1 = 5.44, a_l2 = 0.306),
    list(...))
}

## Seven south-eastern Australian woody species, spanning the wet sclerophyll to
## cool-temperate-rainforest gradient and the dry sclerophyll end of it. They
## co-occur across the Victorian Central Highlands, the Otways and Tasmanian wet
## forests, so the set is an assemblage rather than a list of extremes.
##
## ⚠️ THESE ARE REPRESENTATIVE VALUES ASSEMBLED FOR A DEMONSTRATION, not measured
## data. A study takes them from AusTraits, which is the curated Australian trait
## database and reports exactly these four for these species. They are here so
## the demo has recognisable rows; DO NOT cite them.
##
## Four traits, because four is what a trait database reports and three of them
## drive further parameters through TF24_hyperpar:
##   lma    leaf mass per area   [kg/m2]  -> k_l, r_l, nmass_l
##   rho    wood density         [kg/m3]  -> d_I, g1_TF24, k_s, r_s, r_b
##   hmat   height at MATURATION [m]      (not maximum height: a mountain ash
##                                        first reproduces near 16 m, not 90)
##   omega  seed mass            [kg]     -> a_f3
gd_species <- function() {
  data.frame(
    species = c("Pomaderris aspera", "Acacia dealbata",
                "Atherosperma moschatum", "Nothofagus cunninghamii",
                "Eucalyptus regnans", "Banksia serrata",
                "Allocasuarina littoralis"),
    common = c("hazel pomaderris", "silver wattle", "southern sassafras",
               "myrtle beech", "mountain ash", "old man banksia",
               "black she-oak"),
    niche = c("wet-forest understorey", "nitrogen-fixing pioneer",
              "shade-tolerant subcanopy", "rainforest canopy",
              "fire-killed obligate seeder", "serotinous sclerophyll",
              "dry sclerophyll, dense wood"),
    lma   = c(0.070, 0.090, 0.110, 0.145, 0.130, 0.230, 0.190),
    rho   = c(450,   550,   500,   620,   480,   650,   800),
    hmat  = c(4.0,   6.0,   8.0,   12.0,  16.0,  5.0,   5.0),
    omega = c(3e-7,  1.1e-5, 5e-6, 2e-6,  1.5e-6, 6e-5, 4e-6),
    stringsAsFactors = FALSE)
}

## One species' trait vector, ready for gd_point(). The traits the demo does not
## vary keep their TF24 defaults.
gd_species_traits <- function(row) {
  c(lma = row$lma, rho = row$rho, hmat = row$hmat, omega = row$omega)
}

## Parameters at a trait point. `traits` is a named numeric vector, or a matrix
## with one row per species.
##
## `schedule` transplants a node-introduction schedule from another point instead
## of refining here. DO NOT transplant `ode_times` with it: pinning the ODE grid
## chosen at another trait point makes the leaf fail to place an operating point
## and the gradient refuses every metric, on about half of a +-30% scatter.
gd_parameters <- function(traits, lifetime = 40, birth_rate = 1.10,
                          schedule = NULL) {
  m <- if (is.matrix(traits)) traits else
    matrix(traits, nrow = 1, dimnames = list(NULL, names(traits)))
  p <- scm_base_parameters("TF24")
  p$max_patch_lifetime <- lifetime
  p <- add_strategies(p, m, hyperpar = TF24_hyperpar,
                      birth_rate = as.list(rep(birth_rate, nrow(m))))
  if (!is.null(schedule)) {
    p$node_schedule_times <- schedule
  }
  p
}

## One point: parameters in, a gradient and everything needed to read it out.
##
## `refine` costs about twice what the gradient does and is what makes the
## schedule this point's own. Passing `schedule` and refine = FALSE reuses
## another point's, which is cheaper and is a different function -- see the
## note in the demos on which the workflow wants.
gd_point <- function(traits, lifetime = 40, schedule = NULL, refine = TRUE,
                     ctrl = Control(node_density_in_birth_date = TRUE)) {
  p <- gd_parameters(traits, lifetime = lifetime, schedule = schedule)
  scm <- run_scm(p, Environment("TF24"), ctrl, refine_schedule = refine,
                 collect = FALSE, record_trajectory = TRUE)
  g <- stand_gradient(scm)
  list(traits = traits,
       value = g$value,
       gradient = g$gradient,
       refused = stand_gradient_refused(g),
       x = gd_param_values(scm, colnames(g$gradient)),
       schedule = scm$parameters$node_schedule_times,
       steps = length(scm$ode_times) - 1L)
}

## The parameter value behind every gradient column, in the same order. Columns
## are named "<species>.<parameter>".
gd_param_values <- function(scm, cols) {
  strategies <- scm$parameters$strategies
  vapply(cols, function(cn) {
    sp <- as.integer(sub("\\..*$", "", cn))
    nm <- sub("^[0-9]+\\.", "", cn)
    as.numeric(strategies[[sp]]$pars[[nm]])
  }, numeric(1))
}

## Elasticity: the percentage change in a metric per percentage change in a
## parameter, (x/y) * dy/dx. Dimensionless, so columns are comparable to each
## other and rows to each other, which raw partials are not.
##
## A metric at zero has no elasticity and comes back NaN rather than Inf.
gd_elasticity <- function(gradient, value, x) {
  stopifnot(ncol(gradient) == length(x), nrow(gradient) == length(value))
  e <- sweep(gradient, 2L, x, `*`)
  e <- sweep(e, 1L, ifelse(value == 0, NA_real_, value), `/`)
  dimnames(e) <- dimnames(gradient)
  e
}

## d(parameter)/d(trait), by differencing TF24_hyperpar itself rather than
## restating its algebra here -- the derived forms are its own and a second copy
## would drift from it silently.
##
## filter = FALSE deliberately: with the filter on, a derived column equal to the
## default is dropped, so the set of columns changes between the two evaluations
## of a difference and the Jacobian comes out ragged.
gd_hyperpar_jacobian <- function(traits, strategy = NULL, rel_h = 1e-6) {
  if (is.null(strategy)) strategy <- scm_base_parameters("TF24")$strategy_default
  m0 <- matrix(traits, nrow = 1, dimnames = list(NULL, names(traits)))
  at <- function(m) {
    out <- TF24_hyperpar(m, strategy, filter = FALSE)
    stats::setNames(as.numeric(out[1, ]), colnames(out))
  }
  base <- at(m0)
  J <- matrix(0.0, nrow = length(base), ncol = length(traits),
              dimnames = list(names(base), names(traits)))
  for (k in seq_along(traits)) {
    h <- rel_h * max(abs(traits[[k]]), 1e-8)
    mp <- m0; mp[1, k] <- mp[1, k] + h
    mm <- m0; mm[1, k] <- mm[1, k] - h
    J[, k] <- (at(mp) - at(mm)) / (2 * h)
  }
  J
}

## The trait derivative: the gradient chained through the hyperpar map, so a
## column is d(metric)/d(trait) with every parameter the trait drives moving with
## it. Columns the gradient does not carry contribute nothing, which is right for
## the ones the model declares undifferentiable (nmass_l among them).
##
## Single species: the gradient's columns are matched on the parameter name with
## the species prefix stripped.
gd_trait_gradient <- function(gradient, J, species = 1L) {
  cols <- colnames(gradient)
  keep <- as.integer(sub("\\..*$", "", cols)) == species
  nm <- sub("^[0-9]+\\.", "", cols[keep])
  shared <- intersect(nm, rownames(J))
  g <- gradient[, keep, drop = FALSE][, match(shared, nm), drop = FALSE]
  out <- g %*% J[shared, , drop = FALSE]
  dimnames(out) <- list(rownames(gradient), colnames(J))
  out
}

## The trait-space elasticity, which is what compares against a trait-database
## correlation: percentage change in a metric per percentage change in a trait,
## with the hyperpar-derived parameters carried along.
gd_trait_elasticity <- function(trait_gradient, value, traits) {
  gd_elasticity(trait_gradient, value,
                as.numeric(traits[colnames(trait_gradient)]))
}

## The iso-surface trade-off: holding a metric fixed, the percentage change in
## `b` that offsets a one-percent rise in `a`. Read off one gradient, because one
## gradient is the whole tangent plane.
##
## In elasticities the ratio is unit-free, which is the form a trait-database
## correlation is stated in. A near-zero denominator means `b` does not move this
## metric and the offset is unbounded; those come back NA rather than a large
## number that looks like an answer.
gd_tradeoff <- function(elasticity, metric, tol = 1e-8) {
  e <- elasticity[metric, ]
  n <- length(e)
  out <- matrix(NA_real_, n, n, dimnames = list(names(e), names(e)))
  for (i in seq_len(n)) {
    for (j in seq_len(n)) {
      if (i != j && is.finite(e[[j]]) && abs(e[[j]]) > tol) {
        out[i, j] <- -e[[i]] / e[[j]]
      }
    }
  }
  out
}

## Which parameter is the largest lever on a metric here, and by how much over
## the runner-up. The margin is what says whether a boundary in a phase diagram
## is a real change of regime or two levers swapping inside the noise.
gd_largest_lever <- function(elasticity, metric) {
  e <- abs(elasticity[metric, ])
  e <- e[is.finite(e)]
  if (!length(e)) return(list(name = NA_character_, margin = NA_real_))
  o <- order(e, decreasing = TRUE)
  list(name = names(e)[o[1]],
       margin = if (length(e) > 1) e[o[1]] / max(e[o[2]], 1e-30) else Inf)
}

## The angle between two gradients, in degrees, over shared columns. Taken on
## elasticities: an angle between raw partials is an angle in a space whose axes
## have different units, which is not a property of the model.
gd_angle <- function(a, b) {
  ok <- is.finite(a) & is.finite(b)
  if (!any(ok)) return(NA_real_)
  ca <- sum(a[ok] * b[ok]) / sqrt(sum(a[ok]^2) * sum(b[ok]^2))
  as.numeric(acos(max(-1, min(1, ca))) * 180 / pi)
}
