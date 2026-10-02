## What an invasion run has to reproduce. Every number here was measured on
## develop, against a recorder reached through three hooks the ODE solver called
## into the patch. odelia's rewrite deleted those hooks; what replaces them is
## odelia's own store/load channel -- the run keeps the field in the same
## per-(step, stage) row it already keeps what a rate evaluation solved for, and
## the invasion pass reads it there. The numbers did not move, which is the point of
## keeping them: they were the specification the replacement was written to, not
## a re-pin taken from it.

test_that("mutant method works", {
  # basic setup
  p0 <- scm_base_parameters("FF16")
  p0$max_patch_lifetime <- 50

  e <- Environment("FF16")
  ctrl <- Control()

  tol <- 1e-4
  
  # We'll run tests with 1 and 3 residents, each with different numbers of mutants
  
  lma <- c(0.05, 0.1, 0.2)
  birth_rate <- 1

  # 1 resident strategies
  pr1 <- add_strategies(p0, trait_matrix(lma[2], "lma"), birth_rate = rep(birth_rate, 1))

  pr1m1 <- add_strategies(pr1, trait_matrix(lma[3], "lma"), birth_rate = rep(birth_rate, 1))

  pr1m3 <- add_strategies(pr1, trait_matrix(lma, "lma"), birth_rate = rep(birth_rate, 3))

  pr1m10 <- add_strategies(pr1, trait_matrix(seq(lma[1], lma[3], length.out=10), "lma"), birth_rate = rep(birth_rate, 10))

  # test error handling
  # scm object but not yet run
  types <- extract_RcppR6_template_types(pr1, "Parameters")
  scm <- do.call("SCM", types)(pr1, e, empty_events(), ctrl)

  expect_error(scm$run_mutant(p0), "Run a resident first")

  # check mutant fitness against resindet and expected values
  scm <- run_scm(pr1, e, ctrl)
  pr1_rr <- scm$net_reproduction_ratios
  expected <- 2.77322
  expect_equal(pr1_rr, expected, tolerance = tol)

  scm$run_mutant(pr1m1)
  pr1m1_rr <- scm$net_reproduction_ratios
  expected <- c(2.77322, 3.707605)
  expect_equal(pr1m1_rr, expected, tolerance = tol)
  expect_equal(pr1m1_rr[1], pr1_rr, tolerance = tol)

  scm$run_mutant(pr1m3)
  pr1m3_rr <- scm$net_reproduction_ratios
  expected <- c(2.77322, 3.7429e-10, 2.77322, 3.70753)
  expect_equal(pr1m3_rr, expected, tolerance = tol)
  expect_equal(pr1m3_rr[1], pr1_rr, tolerance = tol)

  scm$run_mutant(pr1m10)
  pr1m10_rr <- scm$net_reproduction_ratios
  expected <- c(2.773222, 3.742935e-10, 9.308944e-07, 0.1363641, 2.773222, 3.890554, 1.524582, 1.160212, 1.871261, 2.765328, 3.707372)
  expect_equal(pr1m10_rr, expected, tolerance = tol)
  expect_equal(pr1m10_rr[1], pr1_rr, tolerance = tol)

  # 3 resident strategies
  pr3 <- add_strategies(p0, trait_matrix(lma, "lma"), birth_rate = rep(birth_rate, 3))
  
  pr3m1 <- add_strategies(pr3, trait_matrix(lma[3], "lma"), birth_rate = rep(birth_rate, 1))

  pr3m3 <- add_strategies(pr3, trait_matrix(lma, "lma"), birth_rate = rep(birth_rate, 3))

  pr3m10 <- add_strategies(pr3, trait_matrix(seq(lma[1], lma[3], length.out = 10), "lma"), birth_rate = rep(birth_rate, 10))

  scm <- run_scm(pr3, e, ctrl)
  pr3_rr <- scm$net_reproduction_ratios
  expected <- c(4.265e-10, 2.831741, 0.09125339)
  expect_equal(pr3_rr, expected, tolerance = tol)


  scm$run_mutant(pr3m1)
  pr3m1_rr <- scm$net_reproduction_ratios
  expected <- c(4.265e-10, 2.831741, 0.09125339, 0.09125339)
  expect_equal(pr3m1_rr, expected, tolerance = tol)
  expect_equal(pr3m1_rr[1:3], pr3_rr, tolerance = tol)

  scm$run_mutant(pr3m3)
  pr3m3_rr <- scm$net_reproduction_ratios
  expected <- c(4.265e-10, 2.831741, 0.09125339, 4.265e-10, 2.831741, 0.09125339)
  expect_equal(pr3m3_rr, expected, tolerance = tol)
  expect_equal(pr3m3_rr[1:3], pr3_rr, tolerance = tol)

  scm$run_mutant(pr3m10)
  pr3m10_rr <- scm$net_reproduction_ratios
  expected <- c(4.265011e-10, 2.831741, 0.09125377, 4.265011e-10, 5.587752e-06, 0.266188, 2.831741, 2.690585, 0.3796333, 0.07098642, 0.07226859, 0.08342181, 0.09125377)
  expect_equal(pr3m10_rr, expected, tolerance = tol)
  expect_equal(pr3m3_rr[1:3], pr3_rr, tolerance = tol)
})

test_that("mutant method densities", {
  # For a mutant strategy identical to the resident, the mutant method must
  # reproduce exactly the fitness that strategy attains when run as a resident.
  # This is an identity of the machinery rather than a near-equilibrium
  # approximation -- it holds at any birth rate and any patch lifetime (the two
  # sides agree to ~1e-13 below, far inside the 1e-3 tolerance). We therefore
  # check the invariant across a spread of birth rates (including the degenerate
  # zero-birth case) at two patch lifetimes.
  #
  # Short patch lifetimes are used on purpose: because the agreement is
  # lifetime-independent, a shorter patch retains the full strength of the check
  # while running several times faster than the model's default lifetime. (The
  # earlier versions sampled birth rates around a hard-coded equilibrium at the
  # default lifetime, but the test never asserted anything *about* that
  # equilibrium -- only the resident-vs-mutant identity -- so the long, costly
  # patch bought no extra coverage.)
  ctrl <- Control()

  traits <- trait_matrix(0.0825, c("lma"))
  tol <- 1e-3

  # fitness at birth rate x computed two ways: as a resident, and as a mutant
  # of the resident -- which must agree. The identity holds against whatever
  # competitive landscape the resident run produces, so we do NOT refine the
  # cohort schedule first: refinement is irrelevant to the invariant but was the
  # dominant cost (it roughly tripled this block's run time).
  f_test <- function(p, x) {
    p1 <- p
    p1$strategies[[1]]$birth_rate_y <- x

    scm <- run_scm(p1, ctrl = ctrl)
    r_rr <- scm$net_reproduction_ratios

    scm$run_mutant(p1)
    m_rr <- scm$net_reproduction_ratios

    dplyr::tibble(birth_rate = x, resident_f = log(r_rr), mutant_f = log(m_rr))
  }

  run_case <- function(life, birth_rates) {
    p0 <- scm_base_parameters("FF16")
    p0$max_patch_lifetime <- life
    pr1 <- add_strategies(p0, traits, birth_rate = 1)

    outputs <- purrr::map_df(birth_rates, ~ f_test(pr1, .x))

    expect_equal(birth_rates, outputs$birth_rate, tolerance = tol)
    expect_equal(outputs$resident_f, outputs$mutant_f, tolerance = tol)
  }

  run_case(30, c(0, 5, 10, 20))
  run_case(20, c(0, 5, 10, 20))
})

# The TF24 stand the invasion tests below replay: lifetime 6, twenty of the default
# schedule's introductions for every species, and constant rain.
tf24_invasion_fixture <- function(traits = trait_matrix(0, "TF24_floor_lambda_o")) {
  p0 <- scm_base_parameters("TF24")
  p0$max_patch_lifetime <- 6
  p1 <- add_strategies(p0, traits, hyperpar = TF24_hyperpar,
                       birth_rate = rep(1, nrow(traits)))
  full <- p1$node_schedule_times[[1]]
  p1$node_schedule_times <-
    rep(list(full[round(seq(1, length(full), length.out = 20))]), nrow(traits))

  env <- Environment("TF24")
  env$set_soil_water_state(rep(0.428 * 0.5, env$get_soil_number_of_depths()))
  env$extrinsic_drivers_set_constant("rainfall", 1)
  list(p = p1, env = env)
}

test_that("an identical invader repeats the run's fitness exactly, and a recorded run is not repeated, TF24", {
  # An invasion walks the run's accepted steps in the field each evaluation was
  # taken in, applying the same entries: the same events in the same order, then
  # the same introductions. A copy of the run's own strategy therefore makes every
  # evaluation the run made, and its fitness is the run's to every digit -- with
  # no events, and under a pulse, a harvest where it meets the schedule's 19th
  # introduction, a climate extreme and a second harvest.
  #
  # TF24's storage pool refuses a state below empty by throwing, which the run's
  # stepper answers by shrinking and retrying. The walk takes only accepted
  # steps and would fail rather than shrink, so passing without an error is the
  # statement that it meets no refusal.
  #
  # A run that kept its states kept those fields beside them, so the invasion
  # walks its recording and repeats nothing; a run that kept none is repeated
  # first. The identity holds either way.
  #
  # Lifetime 6 with twenty of the default schedule's introductions is 194 steps
  # a run; more introductions buy only more of the same regime.
  introduction <- tf24_invasion_fixture()$p$node_schedule_times[[1]][19]
  schedules <- list(
    none = NULL,
    events = list(rainfall_pulse(time = 1, depth = 0.05),
                  harvest(time = introduction, fraction = 0.5),
                  climate_extreme(time = 3.5, intensity = 5, threshold = 1,
                                  sensitivity = 20),
                  harvest(time = 4, fraction = 0.5)))
  for (name in names(schedules)) {
    for (record in c(FALSE, TRUE)) {
      fixture <- tf24_invasion_fixture()
      p1 <- fixture$p
      ev <- do.call(events, c(list(events_default(p1)), schedules[[name]]))
      scm <- run_scm(p1, env = fixture$env, ctrl = Control(), events = ev,
                     record_trajectory = record)
      run_rr <- scm$net_reproduction_ratios
      label <- paste("with", name, if (record) "after a recorded run")

      # Guards that the walk has something to repeat: a stand that died out
      # would make the identity trivial, and a short recording would make it
      # cheap in the wrong way.
      expect_true(all(is.finite(run_rr)) && all(run_rr > 0))
      expect_gt(length(scm$ode_times), 150)
      expect_equal(scm$patch$species[[1]]$size, 20L)

      run_log <- scm$event_log
      expect_equal(scm$runs, 1)
      expect_no_error(scm$run_mutant(p1))
      expect_equal(scm$runs, if (record) 2 else 3,
                   label = paste("the runs", label))
      expect_identical(scm$net_reproduction_ratios, run_rr,
                       label = paste("the invader's fitness", label))
      # Its log is the run's: the same events, each taking out what it took out.
      expect_identical(scm$event_log$time, run_log$time)
      expect_identical(scm$event_log$applied, run_log$applied)
    }
  }
})

test_that("an invader a little costlier in leaf runs on the run's steps, TF24", {
  # The storage relaxation offset slows every pool, so on the run's steps these
  # invaders' pools stay non-negative. At a zero offset the farther invader's
  # stages go so far below empty that its density overflows on this height stand.
  lma <- scm_base_parameters("TF24")$strategy_default$pars[["lma"]]
  with_offset <- function(p, offset) {
    for (i in seq_along(p$strategies)) {
      s <- p$strategies[[i]]
      s$pars$storage_relaxation_offset <- offset
      p$strategies[[i]] <- s
    }
    p
  }
  invader <- function(m, offset) {
    traits <- trait_matrix(c(0, m * lma), c("TF24_floor_lambda_o", "lma"))
    with_offset(tf24_invasion_fixture(traits)$p, offset)
  }
  stand <- tf24_invasion_fixture()

  scm <- run_scm(stand$p, env = stand$env, ctrl = Control())
  run_rr <- scm$net_reproduction_ratios
  for (m in c(1.001, 1.05)) {
    scm$run_mutant(invader(m, 7 / 365))
    expect_lt(scm$net_reproduction_ratios, run_rr)
  }

  scm <- run_scm(with_offset(stand$p, 0), env = stand$env, ctrl = Control())
  expect_no_error(scm$run_mutant(invader(1.001, 0)))
  expect_error(scm$run_mutant(invader(1.05, 0)), "density")
})

test_that("an invader introduced where the run introduced nothing is refused", {
  # A walk applies an entry wherever the recording inserts, so the invaders'
  # introductions have to be the run's; one elsewhere would be skipped with every
  # number finite.
  p <- add_strategies(scm_base_parameters("FF16"), trait_matrix(0.0825, "lma"))
  p$max_patch_lifetime <- 10
  times <- p$node_schedule_times[[1]]
  p$node_schedule_times <- list(times[times <= 10])
  scm <- run_scm(p)
  run_times <- p$node_schedule_times[[1]]
  invader <- p
  invader$node_schedule_times <- list(sort(c(run_times, 0.25)))
  expect_error(scm$run_mutant(invader), "where the schedule's next is at t=0.25")
  # Past the run's last introduction, where no insertion comes to refuse it.
  invader$node_schedule_times <- list(c(run_times, (max(run_times) + 10) / 2))
  expect_error(scm$run_mutant(invader), "past the rows' last insertion")
})

test_that("two invaders together each have the fitness they have alone, TF24", {
  # Each invader is evaluated in the recorded field and solves for its own leaf
  # operating points, so invading beside another invader changes neither of them.
  lma <- scm_base_parameters("TF24")$strategy_default$pars[["lma"]]
  traits <- trait_matrix(c(0, 0, lma, 0.95 * lma),
                         c("TF24_floor_lambda_o", "lma"))
  stand <- tf24_invasion_fixture(traits[1, , drop = FALSE])
  scm <- run_scm(stand$p, env = stand$env, ctrl = Control())
  run_rr <- scm$net_reproduction_ratios

  scm$run_mutant(tf24_invasion_fixture(traits)$p)
  together <- scm$net_reproduction_ratios
  scm$run_mutant(tf24_invasion_fixture(traits[2, , drop = FALSE])$p)
  alone <- scm$net_reproduction_ratios

  expect_identical(together[1], run_rr)
  expect_identical(together[2], alone)
  # The second invader's fitness is not the run's.
  expect_gt(abs(log(alone) - log(run_rr)), 1)
})

test_that("an invader's sweep across a harvest agrees with a difference and a tangent", {
  # On the birth-date coordinate, which the sweep runs on, under a harvest where it
  # meets the stand's second introduction. The gradient holds the recorded field
  # fixed, and each side of the difference invades against the same recording, so
  # both take its steps and apply its harvest.
  p <- ladder_parameters("fast")
  p$node_schedule_times <- list(c(0, 0.63))
  ev <- events(events_default(p), harvest(time = 0.63, fraction = 0.5))
  run <- function(q) run_scm(q, Environment("TF24"), ladder_control(), events = ev)
  with_parameter <- function(q, name, value) {
    strategies <- q$strategies
    pars <- strategies[[1]]$pars
    pars[[name]] <- value
    strategies[[1]]$pars <- pars
    q$strategies <- strategies
    q
  }
  invader <- with_parameter(p, "lma", 0.95 * p$strategies[[1]]$pars[["lma"]])

  scm <- run(p)
  scm$run_mutant(invader)
  columns <- c("1.hmat", "1.k_I")
  # The invasion kept no states, so the sweep repeats it to keep them.
  swept <- stand_gradient(scm, traits = columns)$gradient
  # The tangent walks the same recording forward, through the same entries.
  every <- census_trait_names_tf24(scm)
  for (column in columns) {
    tangent <- ladder_trajectory_tangent(scm, as.numeric(every == column))$tangent
    expect_equal(tangent, unname(swept[, column]), tolerance = 1e-12)
  }
  scm$record_trajectory <- FALSE

  # lma is left out: its difference has a floor near 1e-5 without an invasion too.
  for (name in c("hmat", "k_I")) {
    value <- invader$strategies[[1]]$pars[[name]]
    h <- 1e-5 * value
    census_at <- function(x) {
      scm$run_mutant(with_parameter(invader, name, x))
      stand_census(scm)
    }
    differenced <- (census_at(value + h) - census_at(value - h)) / (2 * h)
    expect_equal(swept[, paste0("1.", name)], differenced, tolerance = 1e-7)
  }

  # Run on its own, the same strategy moves the field it is evaluated in, and its
  # gradient differs by up to a fifth.
  moving <- stand_gradient(run(invader), traits = columns)$gradient
  expect_gt(max(abs(moving / swept - 1)), 0.1)
})
