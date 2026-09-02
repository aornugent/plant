## Generate the per-species gradients the fingerprints demo reads.
##
## Run deliberately; about 25 minutes. Writes
## overstorey_staging/gradient_fingerprints.rds.
##
##   Rscript overstorey_staging/gradient_fingerprints.R
##
## SEED MASS IS EXCLUDED, and that is a finding rather than a simplification.
## fecundity_dt divides by (omega + a_f3) with a_f3 = 3*omega, so the denominator
## is about 4*omega and offspring production goes as 1/omega. TF24 is calibrated
## at omega = 3.8e-5 kg; real Australian seeds run one to two orders below that.
## Given the seven species' own seed masses, six of seven produced gradients from
## 1e+05 to 1e+128 or would not run at all, ordered exactly by how far below the
## default they sat -- and the only well-behaved one was the only species with a
## seed HEAVIER than the default. So seed mass stays at its default here and the
## demo says so.
##
## THE NUDGE LADDER. The maturation switch is a near-step in height, and a cohort
## crosses it inside one or two ODE steps; which step catches the crossing sets
## the gradient and flips under changes of 1e-4 in a trait. So a single gradient
## per species is not trustworthy: each species is evaluated at a ladder of tiny
## lma offsets, and an answer is kept only where THREE consecutive rungs agree.
## Two is not enough -- one species had a pair agreeing within 2.5x sitting
## between rungs of 826137 and 5144.

library(odelia)
library(plant)

pkg_root <- if (file.exists("DESCRIPTION")) "." else ".."
source(file.path(pkg_root, "overstorey_staging", "gradient_demo_helpers.R"))

LIFETIME <- as.numeric(Sys.getenv("FP_LIFETIME", "40"))
OUT <- Sys.getenv("FP_OUT", file.path(pkg_root, "overstorey_staging",
                                      "gradient_fingerprints.rds"))
RUNGS <- c(0, 1e-3, 2e-3, 3e-3, 4e-3, 5e-3, 6e-3)
AGREE <- 1.5   # three consecutive rungs must sit inside this factor
METRIC <- "mass_above_ground"

ctrl <- Control(node_density_in_birth_date = TRUE)
spp <- gd_species()

# One schedule for every species, refined at the assemblage median, so the
# fingerprints are derivatives of the same discrete model.
mid <- c(lma = stats::median(spp$lma), rho = stats::median(spp$rho),
         hmat = stats::median(spp$hmat))
message("refining the shared schedule once")
S <- run_scm(gd_parameters(mid, LIFETIME), Environment("TF24"), ctrl,
             refine_schedule = TRUE, collect = FALSE)$parameters$node_schedule_times
message("  ", length(S[[1]]), " nodes")

probe <- function(tr) {
  tryCatch({
    p <- gd_parameters(tr, LIFETIME, schedule = S)
    scm <- run_scm(p, Environment("TF24"), ctrl, refine_schedule = FALSE,
                   collect = FALSE, record_trajectory = TRUE)
    g <- stand_gradient(scm)
    if (any(stand_gradient_refused(g))) return(NULL)
    x <- gd_param_values(scm, colnames(g$gradient))
    list(gradient = g$gradient, value = g$value,
         elasticity = gd_elasticity(g$gradient, g$value, x))
  }, error = function(e) NULL)
}

rows <- vector("list", nrow(spp))
for (i in seq_len(nrow(spp))) {
  tr0 <- c(lma = spp$lma[i], rho = spp$rho[i], hmat = spp$hmat[i])
  rungs <- lapply(RUNGS, function(d) {
    tr <- tr0; tr[["lma"]] <- tr[["lma"]] * (1 + d); probe(tr)
  })
  norms <- vapply(rungs, function(r) if (is.null(r)) NA_real_ else {
    e <- r$elasticity[METRIC, ]; sqrt(sum(e[is.finite(e)]^2)) }, numeric(1))

  # The middle of the first run of three consecutive rungs inside AGREE.
  pick <- NA_integer_
  for (j in seq_len(length(RUNGS) - 2)) {
    w <- norms[j + 0:2]
    if (!anyNA(w) && max(w) / min(w) < AGREE) { pick <- j + 1L; break }
  }
  message(sprintf("  %-26s %s -> %s", spp$species[i],
                  paste(sprintf("%.2f", norms), collapse = " "),
                  if (is.na(pick)) "no stable run of three" else
                    sprintf("rung %d", pick)))
  rows[[i]] <- list(species = spp$species[i], common = spp$common[i],
                    niche = spp$niche[i], traits = tr0, rungs = RUNGS,
                    norms = norms, pick = pick,
                    chosen = if (is.na(pick)) NULL else rungs[[pick]])
}

saveRDS(list(rows = rows, species = spp, lifetime = LIFETIME,
             nodes = length(S[[1]]), rungs = RUNGS, agree = AGREE,
             metric = METRIC, generated = Sys.time()), OUT)
message("wrote ", OUT)
