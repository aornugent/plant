## Generate the trait-space grid the phase-diagram demo reads.
##
## Run this deliberately; it is about 45 minutes and writes
## overstorey_staging/gradient_phase_grid.rds. The demo does not regenerate it,
## because a document that recomputes an hour of gradients is a document nobody
## renders.
##
##   Rscript overstorey_staging/gradient_phase_grid.R
##
## ONE node schedule, refined once at the grid centre and shared by every cell.
## That is not only the cheaper half -- it is what makes the cells comparable.
## Refinement chooses the node set from the forward leaf-area and seed-rain
## errors, and the gradient is O(1)-sensitive to which node set it got: over a
## schedule_eps sweep the census moved by 4e-4 while the gradient moved by more
## than its own norm. Sharing the schedule holds the node set constant across the
## grid, so a boundary in the map is the model changing its mind rather than the
## mesh changing under it.
##
## DO NOT also transplant ode_times. Pinning the ODE grid chosen at the centre
## makes the leaf fail to place an operating point out at the corners, and the
## gradient then refuses every metric on about a third of the cells.

library(odelia)
library(plant)

# Run from the package root or from overstorey_staging/; both resolve.
pkg_root <- if (file.exists("DESCRIPTION")) "." else ".."
source(file.path(pkg_root, "overstorey_staging", "gradient_demo_helpers.R"))

LIFETIME <- as.numeric(Sys.getenv("GRID_LIFETIME", "30"))
N        <- as.integer(Sys.getenv("GRID_N", "8"))
SPAN     <- as.numeric(Sys.getenv("GRID_SPAN", "0.25"))
OUT      <- Sys.getenv("GRID_OUT",
                       file.path(pkg_root, "overstorey_staging",
                                 "gradient_phase_grid.rds"))

ctrl <- Control(node_density_in_birth_date = TRUE)
tr0 <- unlist(gd_traits())

# The two axes: the leaf economic trait and the size at maturity. Both are set
# through trait_matrix, both are in a trait database, and lma drives three
# further parameters through TF24_hyperpar while hmat drives none -- so the pair
# also shows the hyperpar chain doing work on one axis and nothing on the other.
AX <- c("lma", "hmat")
grid_at <- function(fx, fy) { tr <- tr0; tr[[AX[1]]] <- tr0[[AX[1]]] * fx
                              tr[[AX[2]]] <- tr0[[AX[2]]] * fy; tr }
fx <- exp(seq(log(1 - SPAN), log(1 + SPAN), length.out = N))
fy <- exp(seq(log(1 - SPAN), log(1 + SPAN), length.out = N))

message("refining the centre schedule once")
centre <- run_scm(gd_parameters(tr0, LIFETIME), Environment("TF24"), ctrl,
                  refine_schedule = TRUE, collect = FALSE)
S <- centre$parameters$node_schedule_times
message("  ", length(S[[1]]), " nodes")

cells <- vector("list", N * N)
i <- 0L
t0 <- Sys.time()
for (a in seq_len(N)) for (b in seq_len(N)) {
  i <- i + 1L
  tr <- grid_at(fx[a], fy[b])
  cell <- tryCatch({
    pt <- gd_point(tr, LIFETIME, schedule = S, refine = FALSE, ctrl = ctrl)
    e <- gd_elasticity(pt$gradient, pt$value, pt$x)
    J <- gd_hyperpar_jacobian(tr)
    gt <- gd_trait_gradient(pt$gradient, J)
    list(fx = fx[a], fy = fy[b], traits = tr, value = pt$value,
         gradient = pt$gradient, elasticity = e,
         trait_gradient = gt,
         trait_elasticity = gd_trait_elasticity(gt, pt$value, tr),
         refused = pt$refused, steps = pt$steps, error = NULL)
  }, error = function(e) list(fx = fx[a], fy = fy[b], traits = tr,
                              error = conditionMessage(e)))
  cells[[i]] <- cell
  if (i %% 8L == 0L)
    message(sprintf("  %d/%d  %.1f min elapsed", i, N * N,
                    as.numeric(difftime(Sys.time(), t0, units = "mins"))))
}

saveRDS(list(axes = AX, fx = fx, fy = fy, centre = tr0, lifetime = LIFETIME,
             nodes = length(S[[1]]), span = SPAN, cells = cells,
             generated = Sys.time()), OUT)
message("wrote ", OUT, " in ",
        sprintf("%.1f", as.numeric(difftime(Sys.time(), t0, units = "mins"))),
        " min")
