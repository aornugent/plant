# Collar-census probe. One production lifetime per arm; the census header is
# gated on PLANT_COLLAR_CENSUS, so an unset variable reproduces the tree exactly.
#
# Rscript scripts/collar_census.R <max_patch_lifetime> <rain> [census-output-path]
library(odelia)
pkgload::load_all("/home/user/wt-p3-pinned", quiet = TRUE)

args <- commandArgs(trailingOnly = TRUE)
lifetime <- as.numeric(args[[1]])
rain <- as.numeric(args[[2]])

p0 <- scm_base_parameters("TF24", "TF24_Env")
p0$max_patch_lifetime <- lifetime
p <- add_strategies(p0, trait_matrix(0.1978791, "lma"))

env <- Environment("TF24")
if (is.finite(rain)) env$extrinsic_drivers_set_constant("rainfall", rain)

scm <- run_scm(p, env, Control(), collect = FALSE, refine_schedule = FALSE)
cat(sprintf("%.17g", scm$offspring_production), length(scm$ode_times), "\n")
