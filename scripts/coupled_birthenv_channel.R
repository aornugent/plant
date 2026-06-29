# Step 2 measurement (#472 scope B): the birth-env channel = AD(active birth-env) -
# AD(frozen birth-env), and whether it closes part of the full-SCM-FD gap.
suppressMessages({library(devtools); load_all(".", compile = FALSE, quiet = TRUE)})
traits <- c("lma", "a_p1", "a_l1")

p <- scm_base_parameters("FF16")
p <- add_strategies(p, trait_matrix(0.0825, "lma"), hyperpar = FF16_hyperpar,
                    birth_rate = list(20))
pr <- run_scm(p, Environment("FF16"), control(), refine_schedule = TRUE)$parameters
scm <- run_scm(pr, Environment("FF16"), control(save_RK45_cache = TRUE),
               refine_schedule = FALSE)
h <- plant:::ff16_harvest(scm, 1L, NULL); patch <- scm$patch
nn_h <- patch$stand_newnode_height_stage_history
nn_c <- patch$stand_newnode_competition_stage_history

ad_at <- function(active_be) plant:::ff16_coupled_gradient_impl(
  h$pp, h$eh, h$sh, h$birth_step, h$ppsurv, h$ppsab, h$tw, traits, "LAI",
  h$birth_rate, nn_h, nn_c, h$patch_area, active_be)$jacobian["LAI", ]
ad_active <- ad_at(TRUE)
ad_frozen <- ad_at(FALSE)

lai_full <- function(trait, delta) {
  q <- pr; q$strategies[[1]]$pars[[trait]] <- q$strategies[[1]]$pars[[trait]] + delta
  s2 <- run_scm(q, Environment("FF16"), control(save_RK45_cache = TRUE),
                refine_schedule = FALSE)
  e <- s2$patch$environment_history[[length(s2$patch$environment_history)]][[6]]
  -log(e$get_environment_at_height(0))
}

cat(sprintf("%-5s  AD(frozenBE)  AD(activeBE)  birthEnv-chan  full-SCM-FD  ",
            "trait"))
cat("gap(frz)  gap(act)\n")
for (tr in traits) {
  d <- 1e-3 * abs(pr$strategies[[1]]$pars[[tr]])
  fd <- (lai_full(tr, +d) - lai_full(tr, -d)) / (2 * d)
  gfrz <- abs(ad_frozen[tr] - fd) / max(abs(fd), 1e-30)
  gact <- abs(ad_active[tr] - fd) / max(abs(fd), 1e-30)
  cat(sprintf("%-5s  % .4e  % .4e  % .4e  % .4e  %.1f%%   %.1f%%\n",
              tr, ad_frozen[tr], ad_active[tr], ad_active[tr] - ad_frozen[tr],
              fd, 100 * gfrz, 100 * gact))
}
