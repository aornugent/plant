# Generic calibration-facing stand-gradient engine (#472 scope B, build-order
# step 1). CI-runnable in plain R: the engine (ff16_stand_gradient_impl,
# ff16_state_jacobian_impl) is compiled into plant.so with the XAD adjoint tape
# resolved at load against odelia, so no on-the-fly compilation is needed.

test_that("stand_gradient reproduces offspring_production_gradient as one entry", {
  p <- scm_base_parameters("FF16")
  p <- add_strategies(p, trait_matrix(0.0825, "lma"), hyperpar = FF16_hyperpar,
                      birth_rate = list(20))
  p <- run_scm(p, Environment("FF16"), control(), refine_schedule = TRUE)$parameters
  scm <- run_scm(p, Environment("FF16"), control(save_RK45_cache = TRUE),
                 refine_schedule = FALSE)

  tr <- c("a_p1", "lma", "theta", "omega")
  g_ded <- as.numeric(offspring_production_gradient(scm, traits = tr))
  g_eng <- stand_gradient(scm, metrics = "offspring_production", traits = tr)
  # The generic engine is bit-for-bit the dedicated routine for this metric.
  expect_equal(as.numeric(g_eng$jacobian["offspring_production", ]), g_ded,
               tolerance = 1e-10)
  expect_equal(g_eng$values[["offspring_production"]], scm$offspring_production[[1]],
               tolerance = 1e-4)
})

test_that("stand_gradient census metrics reconstruct + match a frozen-resident FD", {
  p <- scm_base_parameters("FF16")
  p <- add_strategies(p, trait_matrix(0.0825, "lma"), hyperpar = FF16_hyperpar,
                      birth_rate = list(20))
  p <- run_scm(p, Environment("FF16"), control(), refine_schedule = TRUE)$parameters
  scm <- run_scm(p, Environment("FF16"), control(save_RK45_cache = TRUE),
                 refine_schedule = FALSE)

  mets <- c("offspring_production", "LAI", "biomass", "size_moment")
  tr   <- c("a_p1", "lma", "a_l1")
  g <- stand_gradient(scm, metrics = mets, traits = tr)

  # LAI reconstructs the SCM's own compute_competition(0) (the size-distribution
  # leaf-area integral, including the pending-seed tail term).
  expect_equal(g$values[["LAI"]], scm$patch$compute_competition(0), tolerance = 1e-6)
  expect_equal(g$values[["offspring_production"]], scm$offspring_production[[1]],
               tolerance = 1e-4)
  expect_equal(dim(g$jacobian), c(length(mets), length(tr)))

  # Frozen-resident two-pass FD of each metric value (perturb a trait in the
  # parameter vector, re-reduce over the SAME harvested schedule + resident light).
  h <- plant:::ff16_harvest(scm, 1L, NULL)
  val_at <- function(q) plant:::ff16_stand_gradient_impl(
    q, h$eh, h$sh, h$birth_step, h$ppsurv, h$ppsab, h$tw, tr, mets, h$birth_rate)$values
  fd_best <- function(trait, mm, rel_h) {
    fds <- vapply(rel_h, function(rh) {
      hh <- rh * abs(h$pp[[trait]]); q1 <- q2 <- h$pp
      q1[[trait]] <- q1[[trait]] + hh; q2[[trait]] <- q2[[trait]] - hh
      (val_at(q1)[[mm]] - val_at(q2)[[mm]]) / (2 * hh)
    }, numeric(1))
    fds[which.min(abs(fds - g$jacobian[mm, trait]))]
  }
  for (mm in mets) {
    expect_equal(g$jacobian[mm, "a_p1"], fd_best("a_p1", mm, c(1e-5, 1e-6)),
                 tolerance = 1e-3)
    # lma exercises the IFT seedling-size (height_0) path through every metric.
    expect_equal(g$jacobian[mm, "lma"], fd_best("lma", mm, c(1e-4, 1e-5)),
                 tolerance = 1e-3)
  }
})

test_that("stand_state_jacobian matches a frozen-resident FD per cohort", {
  p <- scm_base_parameters("FF16")
  p <- add_strategies(p, trait_matrix(0.0825, "lma"), hyperpar = FF16_hyperpar,
                      birth_rate = list(20))
  p$node_schedule_times <- list(seq(0, 60, length.out = 9))
  p$max_patch_lifetime <- 60
  scm <- run_scm(p, Environment("FF16"), control(save_RK45_cache = TRUE),
                 refine_schedule = FALSE)

  tr <- c("a_p1", "lma")
  J <- stand_state_jacobian(scm, traits = tr)
  nC <- nrow(J$states)
  expect_equal(dim(J$jacobian), c(nC, 6L, length(tr)))

  # Frozen-resident FD of the cohort final states, perturbing a_p1.
  h <- plant:::ff16_harvest(scm, 1L, NULL)
  states_at <- function(q) plant:::ff16_state_jacobian_impl(
    q, h$eh, h$sh, h$birth_step, h$ppsurv, h$ppsab, h$tw, tr)$states
  hh <- 1e-6 * abs(h$pp[["a_p1"]]); q1 <- q2 <- h$pp
  q1[["a_p1"]] <- q1[["a_p1"]] + hh; q2[["a_p1"]] <- q2[["a_p1"]] - hh
  fd <- (states_at(q1) - states_at(q2)) / (2 * hh)
  ti <- match("a_p1", tr)
  # Check a mid and a late cohort's height + offspring sensitivity.
  for (i in c(2L, nC - 1L)) for (cc in c("height", "offspring")) {
    ci <- match(cc, dimnames(J$states)[[2]])
    expect_equal(unname(J$jacobian[i, ci, ti]), unname(fd[i, ci]),
                 tolerance = 1e-3)
  }
})
