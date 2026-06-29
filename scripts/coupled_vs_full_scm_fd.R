suppressMessages({library(devtools); load_all(".", compile=FALSE, quiet=TRUE)})
p <- scm_base_parameters("FF16"); p <- add_strategies(p, trait_matrix(0.0825,"lma"), hyperpar=FF16_hyperpar, birth_rate=list(20))
pr <- run_scm(p, Environment("FF16"), control(), refine_schedule=TRUE)$parameters
scm <- run_scm(pr, Environment("FF16"), control(save_RK45_cache=TRUE), refine_schedule=FALSE)
h <- plant:::ff16_harvest(scm,1L,NULL); patch<-scm$patch

# coupled AD dLAI/dtheta (LAI == competition(0) at final census)
g <- plant:::ff16_coupled_gradient_impl(h$pp,h$eh,h$sh,h$birth_step,h$ppsurv,h$ppsab,h$tw,
       c("lma","a_p1","a_l1"),"LAI",h$birth_rate,
       patch$stand_newnode_height_stage_history,patch$stand_newnode_competition_stage_history,h$patch_area)

# LAI from a FULL re-run on the SAME (frozen) node schedule -- run_scm lets the birth
# canopy respond, so FD here = AD + the dropped birth-env/boundary channels.
lai_full <- function(trait, delta) {
  pp <- pr$strategies[[1]]$pars; pp[[trait]] <- pp[[trait]] + delta
  q <- pr; q$strategies[[1]]$pars[[trait]] <- pp[[trait]]
  s2 <- run_scm(q, Environment("FF16"), control(save_RK45_cache=TRUE), refine_schedule=FALSE)
  e <- s2$patch$environment_history[[length(s2$patch$environment_history)]][[6]]
  list(lai = -log(e$get_environment_at_height(0)),
       used = s2$parameters$strategies[[1]]$pars[[trait]])  # confirm perturbation kept
}
cat(sprintf("%-6s  coupled-AD      full-SCM FD     rel.diff\n","trait"))
for (tr in c("lma","a_p1","a_l1")) {
  d <- 1e-3*abs(pr$strategies[[1]]$pars[[tr]])
  hp <- lai_full(tr,+d); hm <- lai_full(tr,-d)
  # sanity: perturbation survived prepare_strategy (raw input pars only)
  fd <- (hp$lai - hm$lai)/(2*d)
  ad <- g$jacobian["LAI",tr]
  cat(sprintf("%-6s  % .6e  % .6e  %.2e\n", tr, ad, fd, abs(ad-fd)/pmax(abs(fd),1e-30)))
}
