## ---- echo = FALSE-------------------------------------------------------
library(razimuth)

## ------------------------------------------------------------------------
set.seed(3141)
grouse <- sim_atd(n_loc = 100, sq_km = 4, n_azimuth = 3, dist_vec = c(200,300), 
                  kappa_sim = 50, prior_r = 1000, plot = FALSE)$sim_df
head(x = grouse, n = 3)

## ------------------------------------------------------------------------
grouse_atm <- convert_atm(df = grouse)

## ---- out.width = '80%', fig.align = 'center', fig.height = 7, fig.width = 7----
visualize_atm(atm_df = grouse_atm, obs_id = 93, add_prior = T)

## ---- hide = T-----------------------------------------------------------
atm_fit <- atm_mcmc(atm_df = grouse_atm, n_mcmc = 11000, n_burn = 1000)

## ---- out.width = '80%', fig.align = 'center', fig.height = 7, fig.width = 7----
id_tmp <- which_pid(atm_df = grouse_atm, obs_id = 93)
visualize_atm(atm_df = grouse_atm, obs_id = id_tmp, add_prior = TRUE)
p_isopleth(df = atm_fit$mu_ls[[id_tmp]]$pdraws, prob_lvls = c(0.5,0.9), range_extend = 0,
           kde_n = 50, col_vec = c(4,4))
points(matrix(atm_fit$pmode[id_tmp, 2:3], ncol = 2), pch = 21, bg = 4)
legend("topleft", c("Posterior Mode"), pch = 21, pt.bg = 4, bty = "n")

## ---- out.width = '80%', fig.align = 'center', fig.height = 7, fig.width = 7----
visualize_atm(atm_df = grouse_atm, obs_id = id_tmp, add_prior = T, add_legend = F)
points(atm_fit$mu_ls[[id_tmp]]$pdraws, pch = 20,
       col = viridis::viridis(n = nrow(atm_fit$mu_ls[[id_tmp]]$pdraws)))

## ---- out.width = '80%', fig.align = 'center', fig.height = 7, fig.width = 10----
plot_kappa(atm_obj = atm_fit$kappa_ls, item = "traceplot")
plot_kappa(atm_obj = atm_fit$kappa_ls, item = "run_mean")

