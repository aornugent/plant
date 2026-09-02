# Smoke test for the staging reverse-mode gradient demos. Exercises the exact
# helpers the two vignettes use, so they cannot silently rot.
# overstorey_staging/ is .Rbuildignore'd (not installed), so this is dev-only:
# skip when the helper file is absent (installed/CRAN checks).
#
# A short stand on purpose: this is a smoke test of the helpers, not of the
# model, and the demos run at 30 years.

helpers_or_skip <- function() {
  skip_on_cran()
  helpers <- test_path("..", "..", "overstorey_staging",
                       "gradient_demo_helpers.R")
  skip_if_not(file.exists(helpers),
              "overstorey_staging/ not present (built package)")
  helpers
}

test_that("one point gives a gradient with a value behind every column", {
  source(helpers_or_skip(), local = TRUE)
  tr <- gd_traits()
  pt <- gd_point(unlist(tr), lifetime = 15)

  expect_equal(nrow(pt$gradient), 3L)
  expect_equal(ncol(pt$gradient), length(pt$x))
  expect_identical(colnames(pt$gradient), names(pt$x))
  expect_true(all(is.finite(pt$x)))
  expect_false(any(pt$refused))
  expect_true(all(is.finite(pt$gradient)))
  # Every column carries its species index, which is what keeps two species'
  # like-named parameters apart.
  expect_true(all(grepl("^[0-9]+\\.", colnames(pt$gradient))))
})

test_that("elasticity is dimensionless and rescaling a parameter cannot change it", {
  source(helpers_or_skip(), local = TRUE)
  g <- matrix(c(2, 4, 6, 8), nrow = 2,
              dimnames = list(c("m1", "m2"), c("1.a", "1.b")))
  v <- c(m1 = 10, m2 = 20)
  x <- c(`1.a` = 5, `1.b` = 0.5)
  e <- gd_elasticity(g, v, x)
  expect_equal(unname(e["m1", "1.a"]), 2 * 5 / 10)
  expect_equal(unname(e["m2", "1.b"]), 8 * 0.5 / 20)

  # The property that makes it the right summary: state a parameter in different
  # units -- x scaled by k, so dy/dx scaled by 1/k -- and the elasticity is the
  # same number. A raw partial is not.
  k <- 1000
  e2 <- gd_elasticity(g * c(1 / k, 1), v, x * c(k, 1))
  expect_equal(unname(e2["m1", "1.a"]), unname(e["m1", "1.a"]))

  # A metric at zero has no elasticity, and must not come back as Inf.
  expect_true(all(is.na(gd_elasticity(g, c(m1 = 0, m2 = 20), x)["m1", ])))
})

test_that("the hyperpar Jacobian finds the parameters a trait actually drives", {
  source(helpers_or_skip(), local = TRUE)
  tr <- unlist(gd_traits())
  J <- gd_hyperpar_jacobian(tr)

  expect_equal(colnames(J), names(tr))
  # A trait is its own parameter, so the diagonal is one.
  for (k in names(tr)) expect_equal(unname(J[k, k]), 1, tolerance = 1e-6)

  # TF24_hyperpar derives k_l, r_l and nmass_l from lma and nothing from hmat.
  # Getting this wrong is the whole reason the chain exists, so it is asserted
  # rather than trusted.
  expect_true(all(c("k_l", "r_l", "nmass_l") %in% rownames(J)))
  expect_true(all(J[c("k_l", "r_l", "nmass_l"), "lma"] != 0))
  expect_true(all(J[c("k_l", "r_l", "nmass_l"), "hmat"] == 0))
})

test_that("chaining changes the trait lma carries and leaves the others alone", {
  source(helpers_or_skip(), local = TRUE)
  tr <- unlist(gd_traits())
  pt <- gd_point(tr, lifetime = 15)
  J <- gd_hyperpar_jacobian(tr)
  gt <- gd_trait_gradient(pt$gradient, J)

  expect_equal(colnames(gt), names(tr))
  expect_equal(rownames(gt), rownames(pt$gradient))

  # hmat drives no derived parameter, so its chained column must be its own
  # partial. lma drives three, so its column must move. This pair is the control:
  # if both changed, the chain is touching columns it should not.
  raw <- pt$gradient[, paste0("1.", names(tr)), drop = FALSE]
  expect_equal(unname(gt[, "hmat"]), unname(raw[, "1.hmat"]), tolerance = 1e-8)
  expect_false(isTRUE(all.equal(unname(gt[, "lma"]), unname(raw[, "1.lma"]))))
})

test_that("a trade-off is the elasticity ratio, and an unmoved parameter gives NA", {
  source(helpers_or_skip(), local = TRUE)
  e <- matrix(c(2, 1, 0), nrow = 1,
              dimnames = list("m", c("1.a", "1.b", "1.c")))
  tr <- gd_tradeoff(e, "m")
  expect_equal(unname(tr["1.a", "1.b"]), -2)
  expect_equal(unname(tr["1.b", "1.a"]), -0.5)
  # c moves the metric not at all, so nothing offsets a change in it.
  expect_true(all(is.na(tr[, "1.c"])))
  expect_true(all(is.na(diag(tr))))
})

test_that("largest lever reports its margin, and the angle is in degrees", {
  source(helpers_or_skip(), local = TRUE)
  e <- matrix(c(4, 2, 1), nrow = 1,
              dimnames = list("m", c("1.a", "1.b", "1.c")))
  lev <- gd_largest_lever(e, "m")
  expect_equal(lev$name, "1.a")
  expect_equal(unname(lev$margin), 2)

  expect_equal(gd_angle(c(1, 0), c(1, 0)), 0)
  expect_equal(gd_angle(c(1, 0), c(0, 1)), 90)
  expect_equal(gd_angle(c(1, 0), c(-1, 0)), 180)
})

test_that("a shared node schedule holds the node set fixed across trait points", {
  source(helpers_or_skip(), local = TRUE)
  tr <- unlist(gd_traits())
  a <- gd_point(tr, lifetime = 15)

  moved <- tr; moved[["lma"]] <- moved[["lma"]] * 1.05
  b <- gd_point(moved, lifetime = 15, schedule = a$schedule, refine = FALSE)

  # The point of sharing: the gradient is O(1)-sensitive to which node set it
  # got, so two points being comparable requires the same one.
  expect_equal(lengths(b$schedule), lengths(a$schedule))
  expect_equal(b$schedule[[1]], a$schedule[[1]])
  expect_false(any(b$refused))
  expect_true(all(is.finite(b$gradient)))
})
