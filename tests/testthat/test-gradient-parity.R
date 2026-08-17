# The parity gate: for every state the forward model returns a number for, the
# reverse returns a row or names a violated constraint.
#
# ⚠️ READ AS WRITTEN THIS IS SATISFIED BY A SWEEP THAT REFUSES EVERYTHING, which
# is why it is two checks and not one.
#
# The GATE is the first: nothing escapes unnamed. A raw error, a status outside
# the declared set, a finite-looking row under a refusal, an undefined metric
# that reads as a zero one -- each is a failure however few drivers reach it.
#
# The COVERAGE is the second, and it is a measurement rather than a threshold:
# which drivers answer, and for those that do not, the branch named. A regime
# that stops answering fails here; a regime that has never answered is listed
# below by name, so the difference between the two is a line of code rather than
# a reading of the output.

parity_stand <- function(rain, lifetime, k_I = 0.5, amplitude = 0) {
  p <- scm_base_parameters("TF24")
  p$max_patch_lifetime <- lifetime
  tr <- c(lma = 0.0825, hmat = 5.13, k_I = k_I, a_l1 = 5.44, a_l2 = 0.306)
  p <- add_strategies(p, trait_matrix(unname(tr), names(tr)),
                      hyperpar = TF24_hyperpar, birth_rate = list(1.10))
  env <- Environment("TF24")
  env$set_soil_water_state(rep(0.428 * 0.5, 5))
  if (amplitude == 0) {
    env$extrinsic_drivers_set_constant("rainfall", rain)
  } else {
    # Knots per year rather than per run, so the seasonality a driver realises
    # does not depend on the lifetime it is run for.
    t <- seq(0, lifetime,
             length.out = max(200L, as.integer(ceiling(lifetime * 48))))
    env$extrinsic_drivers_set_variable("rainfall", t,
                                       pmax(rain * (1 + amplitude * sin(2 * pi * t)), 0))
  }
  ctrl <- Control()
  ctrl$node_density_in_birth_date <- TRUE
  scm <- SCM("TF24", "TF24_Env")(p, env, ctrl)
  census_clear_operating_point_counts_tf24(scm)
  scm$run()
  scm
}

# The branches that have never returned a row, by the name the refusal carries.
# A refusal naming anything else is a regime that USED to answer, which is what
# this file exists to catch.
#
# ⚠️ THIS LIST IS EMPTY, so any refusal at all fails the coverage check below.
# `shade-death` left it when the shut branches answered; the LIGHT FLOOR left it
# when the row it severs was found to be exactly zero for the model as evaluated
# -- every light below the floor gives a bit-identical census, so the honest row
# is the zero rather than a refusal, and what makes the zero readable is the
# clamp counter rather than a name in this list.
parity_known_gaps <- character(0)

parity_kinds <- c("answered", "zero-slack", "zero-structural", "zero-undeclared",
                  "refused")

# One driver, reduced to what the two directions each said.
parity_of <- function(scm) {
  g <- stand_gradient(scm)
  st <- as.vector(g$status)
  counts <- census_operating_point_counts_tf24(scm)[[1]]
  names(counts) <- census_operating_point_names_tf24()
  list(status = st, gradient = g, refused = any(st == "refused"),
       reason = if (is.null(g$refusal[[1]])) NA_character_ else g$refusal[[1]]$reason,
       refusal = g$refusal[[1]],
       # Read after the sweep, so it is the sweep's own severances and not the
       # forward run's -- the two differ, and the forward one cannot say whether
       # a gradient carries a declared zero.
       clamp = census_clamp_counts_differentiated_tf24(scm)[[1]],
       kinds = counts[counts > 0])
}

parity_drivers <- list(
  list(name = "wet",      rain = 2.00, lifetime = 5),
  list(name = "drought",  rain = 0.10, lifetime = 5),
  list(name = "seasonal", rain = 1.00, lifetime = 5, amplitude = 1.0),
  list(name = "shaded",   rain = 2.00, lifetime = 5, k_I = 20),
  list(name = "clamped",  rain = 2.00, lifetime = 5, k_I = 40)
)

# Each driver is run and swept ONCE. The sweep is the whole cost here -- the run
# is free -- and three checks reading one sweep is the difference between a file
# that costs a minute and one that costs four.
parity_cache <- new.env(parent = emptyenv())
parity_shared <- function() {
  if (is.null(parity_cache$all)) {
    parity_cache$all <- lapply(parity_drivers, function(d) {
      scm <- parity_stand(d$rain, d$lifetime,
                          k_I = if (is.null(d$k_I)) 0.5 else d$k_I,
                          amplitude = if (is.null(d$amplitude)) 0 else d$amplitude)
      c(list(name = d$name, census = stand_census(scm)), parity_of(scm))
    })
  }
  parity_cache$all
}

test_that("no driver the forward model answers leaves the reverse unnamed", {
  # The gate. Every clause fails on one driver, so none of them is a statement
  # about how many drivers there are.
  for (r in parity_shared()) {
    # The forward model has to have answered, or this driver says nothing about
    # the reverse and is not evidence either way.
    expect_true(all(is.finite(r$census)))
    expect_true(all(r$status %in% parity_kinds))

    entries <- unlist(r$gradient$gradient)
    if (r$refused) {
      # Named, always. A refusal a reader cannot name is the raw error the
      # status channel replaced.
      expect_true(is.character(r$reason) && nchar(r$reason) > 0)

      # Located, or honestly unlocated -- never half. A refusal raised in the
      # SWEEP carries species, node and a step range; one raised while forming
      # the census SEEDS has no node loop to be caught in and carries none of
      # them. Both are readable; a partly-filled location is the shape that
      # reads as an answer, so that is what fails here.
      where <- c(r$refusal$species, r$refusal$node,
                 r$refusal$step_first, r$refusal$step_last)
      expect_true(all(where >= 0) || all(where == -1))
      if (all(where >= 0)) {
        expect_gte(r$refusal$step_last, r$refusal$step_first)
      } else {
        message(sprintf("      unlocated: raised forming the census seeds, not in the sweep"))
      }
      # And an undefined metric must not read as a zero one -- the distinction
      # the status channel exists for, asserted at the boundary rather than
      # inside it.
      expect_true(all(is.na(entries[r$status == "refused"])))
    } else {
      expect_true(all(is.finite(entries[r$status == "answered"])))
      expect_gt(sum(r$status == "answered"), 0)
    }
  }
})

test_that("every refusal names a branch that has never answered", {
  # The coverage half, and the reason it is not a count. A regime that stops
  # answering and one that never answered both come back "refused"; only the name
  # tells them apart, so the name is what is asserted.
  answered <- character(0)
  refused <- character(0)
  for (r in parity_shared()) {
    message(sprintf("  %-9s %-9s %s", r$name,
                    if (r$refused) "refused" else "answered",
                    paste(sprintf("%s=%.0f", names(r$kinds), r$kinds), collapse = " ")))
    if (r$refused) {
      refused <- c(refused, r$name)
      message(sprintf("      %s", substr(r$reason, 1, 140)))
      expect_true(any(vapply(parity_known_gaps,
                             function(g) grepl(g, r$reason, fixed = TRUE),
                             logical(1))))
    } else {
      answered <- c(answered, r$name)
    }
  }
  message(sprintf("  parity: %d of %d drivers answer; %s refuse",
                  length(answered), length(parity_shared()),
                  if (length(refused)) paste(refused, collapse = ", ") else "none"))

  # Non-vacuity. Every driver answering is the end point this file was built to
  # reach, so "something still refuses" can no longer be the guard against a
  # vacuous pass -- and dropping it without a replacement would leave a sweep
  # that answers by doing nothing indistinguishable from one that answers by
  # covering the regimes.
  expect_equal(length(answered), length(parity_shared()))
  expect_length(refused, 0)
})

test_that("the driver that answers through a clamp says so", {
  # What replaces the refusal as this file's non-vacuity guard. The `clamped`
  # driver answers, and the reason that is not vacuous is that it reaches the
  # light floor and reports the severance rather than being quietly unaffected by
  # it: an answered gradient carrying a declared zero and one carrying no clamp
  # at all are the same numbers, and only the count separates them.
  by_name <- stats::setNames(parity_shared(), vapply(parity_shared(),
                                                     function(r) r$name, ""))
  nm <- census_clamp_names_tf24()
  at <- function(driver, site) by_name[[driver]]$clamp[[match(site, nm)]]
  for (d in names(by_name)) {
    c_d <- by_name[[d]]$clamp
    message(sprintf("  %-9s %s", d,
                    paste(sprintf("%s=%.0f", nm[c_d > 0], c_d[c_d > 0]),
                          collapse = " ")))
  }

  # Every driver, including the control, severs its light-independent root rows:
  # above the rooting-depth cap the root profile stops reading height, and that
  # is the model rather than a guard. It is the non-vacuity guard because it
  # holds everywhere -- a build reporting zero here has lost the counter, not
  # found a cleaner stand.
  for (d in names(by_name)) {
    expect_gt(at(d, "rooting_depth"), 0)
  }

  # The light floor is the driver-specific one: the shaded driver reaches it and
  # the control does not, so these counts measure the driver rather than the
  # machinery. Both sites, because the crown one is where it binds first.
  expect_gt(at("clamped", "light_floor"), 0)
  expect_gt(at("clamped", "light_floor_crown"), 0)
  expect_gt(at("clamped", "light_floor_crown"), at("clamped", "light_floor"))
  expect_equal(at("wet", "light_floor"), 0)
  expect_equal(at("wet", "light_floor_crown"), 0)

  # And three sites fire on no driver at all, which is the guard census's own
  # entry and a different statement from "it held".
  for (s in c("soil_positivity", "rainfall", "infiltration")) {
    for (d in names(by_name)) {
      expect_equal(at(d, s), 0)
    }
  }
})

test_that("each driver reaches the branch it is here for", {
  # A fixture must be shown to reach what it tests. The drought and seasonal
  # drivers exist to put PINNED operating points on the gradient's path, and a
  # run that answered without ever reaching one would pass the gate while
  # testing nothing; the shaded driver exists to reach the one branch that has
  # no rows.
  by_name <- stats::setNames(parity_shared(), vapply(parity_shared(),
                                                     function(r) r$name, ""))
  reach <- function(nm, kind) {
    k <- by_name[[nm]]$kinds
    if (is.na(k[kind])) 0 else k[[kind]]
  }
  expect_gt(reach("drought", "pinned-dry-root-crit"), 0)
  expect_gt(reach("seasonal", "pinned-dry-root-crit"), 0)
  expect_gt(reach("shaded", "shade-death"), 0)
  expect_gt(reach("clamped", "shade-death"), 0)
  # And the control has to be a control: the wet driver never leaves the branch
  # the gradient was first built for, which is what makes it the one fixture a
  # regression shows up against cleanly.
  expect_equal(unname(reach("wet", "interior")), unname(sum(by_name[["wet"]]$kinds)))
})
