#' @title ATM master function
#' @description Master function implementing an MCMC algorithm for fitting the ATM model with single kappa parameter.
#' @param atm_df object returned by convert_atm()
#' @param n_mcmc number of mcmc iterations
#' @param n_burn number of burn-in iterations
#' @export
#' @useDynLib razimuth
#' @importFrom Rcpp evalCpp
atm_mcmc <- function(atm_df, n_mcmc, n_burn){

  cat("[[--START--]]\n")
  ## pointers ----
  I_n <- length(atm_df$pid)

  ## generate grid for direct approximation ----
  exp_grid_vals <- c(rev(-1*(cumsum(round(exp(1:150 * 0.07)) + 4))), 0, cumsum(round(exp(1:150 * 0.07)) + 4))

  cat("\n[[--generating grids--]]\n")
  grid.list <- list()
  for(i in 1:I_n){
    grid.pts <- as.matrix(sapply(expand.grid(exp_grid_vals + atm_df$loc_pts[i,1],
                                             exp_grid_vals + atm_df$loc_pts[i,2]), as.numeric))
    grid.list[[i]] <- grid.pts[which(c(priorUpts(gridpts = grid.pts,
                                                 obspt = t(atm_df$obsvr_pts[[i]]),
                                                 radius = atm_df$radius_ls[[i]])) == 1),]
    cat(paste0("[",i,"]"))
  }

  ## compute theta tilde for each grid point ----
  # cat("\n")
  theta.tmp <- list()
  for(i in 1:I_n){
    theta.tmp[[i]] <- matrix(data = NA, nrow = nrow(grid.list[[i]]), ncol = ncol(atm_df$obsvr_pts[[i]]))
    for(j in 1:ncol(atm_df$obsvr_pts[[i]])){
      for(m in 1:nrow(grid.list[[i]])){
        theta.tmp[[i]][m,j] <- atan2(y = grid.list[[i]][m,2] - atm_df$obsvr_pts[[i]][2,j],
                                     x = grid.list[[i]][m,1] - atm_df$obsvr_pts[[i]][1,j])
      }
    }
    # cat(paste0("[",i,"]"))
  }

  ## mcmc ----
  cat("\n\n[[--mcmc--]]\n")
  atmfit_out <- atm_vM_1kappa(z.list = atm_df$obsvr_pts, theta.list = atm_df$theta_vals,
                              grid.list = grid.list, theta.tmp = theta.tmp, I_n = I_n,
                              n.mcmc = n_mcmc, n.burn = n_burn, cat_n = 100)

  ## post-processing mcmc output ----
  cat("\n\n[[--post-processing mcmc output--]]\n")
  ## compute location posterior mode ----
  pmode_mat <- matrix(data = NA, nrow = I_n, ncol = 2)
  for(i in 1:I_n){
    mu_tmp <- t(atmfit_out$mu.save[i,,])
    ux <- unique(mu_tmp)
    pmode_mat[i,] <- ux[which.max(colSums(pmode_cpp(x = mu_tmp, ux = ux)$modeI)),]
  }
  pmode_mat <- cbind(atm_df$pid, pmode_mat)
  colnames(pmode_mat) <- c("pid","utm_x","utm_y")

  mu_ls <- list()
  for(i in 1:I_n){
    mu_ls[[i]] <- list()
    mu_ls[[i]][[1]] <- atm_df$pid[i]
    mu_ls[[i]][[2]] <- unique(atm_df$dt_ls[[i]])
    mu_ls[[i]][[3]] <- pmode_mat[i,c("utm_x","utm_y")]
    mu_ls[[i]][[4]] <- t(atmfit_out$mu.save[i,,])
    names(mu_ls[[i]]) <- c("pid","date","pmode","pdraws")
  }

  kappa_ls <- list()
  kappa_ls[[1]] <- atmfit_out$kappa.accept
  kappa_ls[[2]] <- atmfit_out$kappa.save
  kappa_ls[[3]] <- atmfit_out$ktune.save
  names(kappa_ls) <- c("acceptance","pdraws","tuning")
  class(kappa_ls) <- "atm_kappa"

  ## write output ----
  cat("\n[[--END--]]\n")
  list(mu_ls = mu_ls,
       kappa_ls = kappa_ls,
       pmode_mat = pmode_mat)
}

#' @title ATM model MCMC sampler
#' @description MCMC sampler for ATM model with single kappa parameter.
#' @param z.list list of observer coordinates
#' @param theta.list list of observed bearings
#' @param grid.list list of direct approximation grid points
#' @param theta.tmp list of theta values corresponding to grid points
#' @param I_n number of locations
#' @param n.mcmc number of total mcmc iterations
#' @param n.burn number of burn-in iterations
#' @param cat_n mcmc iteration printing interval
#' @useDynLib razimuth
atm_vM_1kappa <- function(z.list, theta.list, grid.list, theta.tmp, I_n, n.mcmc, n.burn, cat_n){

  ## starting values ----
  kappa <- stats::runif(n = 1, min = 2.355644, max = 100)
  mu.mat <- matrix(data = NA, nrow = I_n, ncol = 2)

  ## saved parameter containers ----
  mu.save <- array(data = NA, dim = c(I_n, 2, (n.mcmc - n.burn)))
  kappa.save <- rep(x = NA, times = (n.mcmc - n.burn))
  ktune.save <- matrix(data = NA, ncol = 1, nrow = n.burn/50)

  ## tuning parameter containers ----
  kappa_tune <- 1
  kappa.accept <- rep(x = 0, times = n.mcmc)
  kappa_accept_batch <- 0

  ## begin mcmc loop ----
  for(k in 1:n.mcmc){

    ## sample mu ----
    for(i in 1:I_n){
      mu.mat[i,] <- grid.list[[i]][sample(x = 1:nrow(grid.list[[i]]),
                                          size = 1,
                                          prob = vhf_cc_mu(cp = kappa,
                                                           thetalist = theta.list[[i]],
                                                           thetatmp = theta.tmp[[i]],
                                                           distid = 1)),]
    }

    theta.old <- theta.list
    for(i in 1:I_n){
      for(j in 1:ncol(z.list[[i]])){
        theta.old[[i]][j] <- atan2(y = mu.mat[i,2] - z.list[[i]][2,j],
                                   x = mu.mat[i,1] - z.list[[i]][1,j])
      }
    }

    ## sample kappa ----
    kappa.star <- truncdist::rtrunc(n = 1, spec = "norm", a = 0, b = 1e5, mean = kappa, sd = kappa_tune)

    tmp.sum1 <- rep(x = 0, times = I_n)
    tmp.sum2 <- rep(x = 0, times = I_n)
    for(i in 1:I_n){
      tmp.sum1[i] <- sum(kappa.star * cos(theta.list[[i]] - theta.old[[i]]) - log(2*pi*besselI(x = kappa.star, nu = 0, expon.scaled = T)) - kappa.star)
      tmp.sum2[i] <- sum(kappa * cos(theta.list[[i]] - theta.old[[i]]) - log(2*pi*besselI(x = kappa, nu = 0, expon.scaled = T)) - kappa)
    }
    mh.1 <- sum(tmp.sum1) + truncdist::dtrunc(kappa, spec = "norm", a = 0, b = 1e5, mean = kappa.star, sd = kappa_tune, log = T)
    mh.2 <- sum(tmp.sum2) + truncdist::dtrunc(kappa.star, spec = "norm", a = 0, b = 1e5, mean = kappa, sd = kappa_tune, log = T)
    mh.k <- exp(mh.1 - mh.2)

    if(mh.k > stats::runif(n = 1, min = 0, max = 1)){
      kappa <- kappa.star
      kappa.accept[k] <- 1
      kappa_accept_batch <- kappa_accept_batch + (1/50)
    }

    ## update tuning ----
    if (k %% 50 == 0 & k <= n.burn) {
      update_k <- update_tuning(k = k, accept_tmp = kappa_accept_batch, tune = kappa_tune)
      kappa_accept_batch <- 0
      kappa_tune <- update_k$tune_update
      ktune.save[(k/50),] <- update_k$tune_update
    }

    ## save samples ----
    if(k > n.burn){
      mu.save[,,(k-n.burn)] <- mu.mat
      kappa.save[(k-n.burn)] <- kappa
    }

    ## print mcmc progress ----
    if(k %% cat_n == 0){
      cat(paste0("[",k,"]"))
    }
  }

  ## compute acceptance rates ----
  kappa.accept <-  mean(kappa.accept[-c(1:n.burn)])

  ## write output to list ----
  list(mu.save = mu.save,
       kappa.save = kappa.save,
       kappa.accept = kappa.accept,
       ktune.save = ktune.save)
}


#' @title Update tuning parameters
#' @description Tuning parameter update during burn-in.
#' @param k iteration
#' @param accept_tmp acceptance rate during past 50 iterations
#' @param tune previous batch tuning parameter values
update_tuning <- function(k, accept_tmp, tune){
  delta <- 1 / sqrt(k)
  n <- length(tune)
  for(i in 1:n){
    if(accept_tmp[i] > 0.44){
      tune[i] <- exp(log(tune[i]) + delta)
    } else {
      tune[i] <- exp(log(tune[i]) - delta)
    }
  }
  list(tune_update = tune)
}
