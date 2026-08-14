## Shared harness for the develop-vs-branch comparison.
##
## Every quantity is measured WITHIN one build, against that build's own exact
## reduction, so nothing depends on transferring state between two plant
## versions with different interfaces.

BUILD <- Sys.getenv("PLANT_BUILD")
stopifnot(nzchar(BUILD))
suppressMessages({
  library(odelia)
  pkgload::load_all(BUILD, quiet = TRUE)
})

## Gauss-Legendre on [0,1], for the crown integral.
gl_nodes <- function(n) {
  i <- 1:(n - 1); b <- i / sqrt(4 * i^2 - 1)
  J <- matrix(0, n, n); J[cbind(i, i + 1)] <- b; J[cbind(i + 1, i)] <- b
  e <- eigen(J, symmetric = TRUE)
  x <- rev(e$values); w <- 2 * rev(e$vectors[1, ])^2
  list(x = (x + 1) / 2, w = w / 2)
}
GL <- gl_nodes(240)

## Crown-mean light for a cohort of height h, given a function returning L(z).
##   M(h) = int_0^h L(z) q(z,h) dz, and with t = u^eta the weight is 2(1-t),
## so the sampling sits where the crown's leaf area is.
crown_mean <- function(Lfun, h, eta) {
  sum(GL$w * 2 * (1 - GL$x) * Lfun(h * GL$x^(1 / eta)))
}

heights_of <- function(patch)
  unlist(lapply(patch$species, function(s)
    vapply(s$nodes, function(n) n$height, numeric(1))))

densities_of <- function(patch)
  unlist(lapply(patch$species, function(s)
    vapply(s$nodes, function(n) exp(n$log_density), numeric(1))))

## The field this build actually uses, and the field it should be.
field_error <- function(patch, eta) {
  hs <- heights_of(patch)
  hs <- hs[is.finite(hs) & hs > 0]
  if (!length(hs)) return(NULL)
  fit   <- function(z) vapply(z, function(q)
    patch$environment$get_environment_at_height(q), numeric(1))
  exact <- function(z) exp(-vapply(z, function(q)
    patch$compute_competition(q), numeric(1)))
  mf <- vapply(hs, function(h) crown_mean(fit, h, eta), numeric(1))
  me <- vapply(hs, function(h) crown_mean(exact, h, eta), numeric(1))
  list(heights = hs, rel = abs(mf - me) / me)
}

## An FF16 stand run to `lifetime`, collected so the whole trajectory is
## available. Defaults everywhere else, so this is the shipped configuration.
ff16_scm <- function(lifetime, collect = TRUE) {
  p0 <- scm_base_parameters("FF16")
  p0$max_patch_lifetime <- lifetime
  p1 <- add_strategies(p0, trait_matrix(0.0825, "lma"),
                       hyperpar = FF16_hyperpar, birth_rate = list(20))
  scm <- SCM("FF16", "FF16_Env")(p1, Environment("FF16"), Control())
  scm$collect <- collect
  scm$run()
  scm
}
