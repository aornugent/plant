# AD-6 (#11): the FF16 invasion (offspring) gradient end-to-end, against a
# re-optimising central finite difference. offspring_production_gradient
# differentiates the existing run_mutant at S = active with the canopy read frozen
# (off-tape) and self-competition suppressed; the reverse-mode result must match a
# central FD of the same frozen-canopy run -- perturb the raw strategy parameter the
# AD seeds, re-run run_mutant, read offspring. This is the AD-vs-FD oracle that gates
# the port. It exercises the compiled stand_gradient path in plant.so.

test_that("FF16 invasion offspring gradient matches re-optimising central FD", {
  p0 <- scm_base_parameters("FF16")
  p0$max_patch_lifetime <- 50
  e <- Environment("FF16")
  ctrl <- Control()
  ctrl$save_RK45_cache <- TRUE

  pr1 <- add_strategies(p0, trait_matrix(0.1, "lma"), birth_rate = 1)
  scm <- run_scm(pr1, e, ctrl)

  traits <- c("lma", "hmat")
  ad <- offspring_production_gradient(scm, traits)
  expect_true(all(is.finite(ad)))

  # Re-optimising central FD: perturb the same raw strategy parameter the AD seeds
  # and re-run the frozen-canopy mutant. Richardson-extrapolate the central
  # difference to remove the O(h^2) truncation bias so the comparison isolates the
  # AD's exactness rather than the differencing error.
  fd_offspring <- function(param, delta) {
    p <- pr1
    p$strategies[[1]]$pars[[param]] <- p$strategies[[1]]$pars[[param]] + delta
    scm$run_mutant(p)
    scm$offspring_production[1]
  }
  central <- function(param, x0, h) {
    (fd_offspring(param, h * x0) - fd_offspring(param, -h * x0)) / (2 * h * x0)
  }
  richardson <- function(param, x0, h) {
    (4 * central(param, x0, h / 2) - central(param, x0, h)) / 3
  }

  for (param in traits) {
    x0 <- pr1$strategies[[1]]$pars[[param]]
    fd <- richardson(param, x0, 1e-4)
    expect_equal(unname(ad[param]), fd, tolerance = 1e-4)
  }
})

test_that("FF16 invasion gradient stays finite for a final-step boundary cohort", {
  # A cohort introduced on the final step (birth == N) is established at the seed
  # height and must never be stepped past the last recorded canopy; the active crown
  # derivative would otherwise read one slice beyond environment_history. Force such a
  # boundary cohort (an introduction exactly at max_patch_lifetime) and assert every
  # trait column of the gradient is finite (the zero-height fix).
  life <- 40
  p0 <- scm_base_parameters("FF16")
  p0$max_patch_lifetime <- life
  e <- Environment("FF16")
  ctrl <- Control()
  ctrl$save_RK45_cache <- TRUE

  pr1 <- add_strategies(p0, trait_matrix(0.1, "lma"), birth_rate = 1)
  scm0 <- run_scm(pr1, e, ctrl)
  times0 <- scm0$parameters$node_schedule_times[[1]]
  times_boundary <- sort(unique(c(times0[times0 < life], life)))

  types <- extract_RcppR6_template_types(pr1, "Parameters")
  scm <- do.call("SCM", types)(pr1, e, ctrl)
  scm$set_node_schedule_times(list(times_boundary))
  scm$run()

  # A genuine boundary cohort exists at the final step.
  expect_equal(max(scm$patch$species[[1]]$node_times), life)

  g <- offspring_production_gradient(scm, c("lma", "hmat", "a_l1", "a_l2"))
  expect_length(g, 4)
  expect_true(all(is.finite(g)))
})
