#' @title Visualize azimuthal telemetry data
#' @description Plotting function to visualize user supplied azimuthal telemetry data.
#' @param atm_df object returned by convert_atm()
#' @param obs_id observation id
#' @param add_prior logical specifying whether to plot circular prior
#' @param add_legend logical specifying whether to add legend to plot
#' @importFrom stats dist
#' @importFrom grDevices rgb
#' @importFrom graphics lines points legend abline
#' @export
visualize_atm <- function(atm_df, obs_id, add_prior = FALSE, add_legend = TRUE){

  ## pointer ----
  uid <- which(unique(atm_df$pid) == obs_id)

  ## plotting ----
  if(add_prior){
    circle_tmp_pts <- cbind(atm_df$radius_ls[[uid]][1] * cos(seq(from = 0, to = 2*pi, length.out = 5)) + atm_df$obsvr_pts[[uid]][1,1],
                            atm_df$radius_ls[[uid]][1] * cos(seq(from = 0, to = 2*pi, length.out = 5)) + atm_df$obsvr_pts[[uid]][2,1])
    for(j in 2:length(atm_df$radius_ls[[uid]])){
      circle_tmp_pts <- rbind(circle_tmp_pts,
                              cbind(atm_df$radius_ls[[uid]][j] * cos(seq(from = 0, to = 2*pi, length.out = 5)) + atm_df$obsvr_pts[[uid]][1,j],
                                    atm_df$radius_ls[[uid]][j] * cos(seq(from = 0, to = 2*pi, length.out = 5)) + atm_df$obsvr_pts[[uid]][2,j]))
    }
    plot(x = circle_tmp_pts[,1], y = circle_tmp_pts[,2], type = "n",
         xlab = "Easting (m)", ylab = "Northing (m)", main = paste0("\n Observation ID = ", obs_id), asp = T)
    for(j in 1:length(atm_df$radius_ls[[uid]])){
      lines(x = atm_df$radius_ls[[uid]][j] * cos(seq(from = 0, to = 2*pi, length.out = 360)) + atm_df$obsvr_pts[[uid]][1,j],
            y = atm_df$radius_ls[[uid]][j] * sin(seq(from = 0, to = 2*pi, length.out = 360)) + atm_df$obsvr_pts[[uid]][2,j],
            lwd = 2, col = rgb(0,0,0, alpha = 0.2))
    }
    for(i in 1:length(atm_df$theta_vals[[uid]])){
      segments(x0 = atm_df$obsvr_pts[[uid]][1,i],
               y0 = atm_df$obsvr_pts[[uid]][2,i],
               x1 = atm_df$obsvr_pts[[uid]][1,i] + cos(atm_df$theta_vals[[uid]][i]) * 1e6,
               y1 = atm_df$obsvr_pts[[uid]][2,i] + sin(atm_df$theta_vals[[uid]][i]) * 1e6)
    }
    points(x = atm_df$obsvr_pts[[uid]][1,], y = atm_df$obsvr_pts[[uid]][2,], pch = 21, bg = "gray35")
    points(x = atm_df$loc_pts[uid,1], y = atm_df$loc_pts[uid,2], pch = 21, bg = 2)
    if(add_legend){
      legend("topright", c("observer location", "transmitter estimate"), pch = c(21,21),
             pt.bg = c("gray35",2), bty = "n")
    }
  } else{
    center_pt <- bearing_int_avg(obs_pts = atm_df$obsvr_pts[[uid]], theta = atm_df$theta_vals[[uid]])$int_avg
    tmp_r <- max(as.matrix(dist(rbind(center_pt, t(atm_df$obsvr_pts[[uid]]))))[,1])
    plot(x = tmp_r * cos(seq(from = 0, to = 2*pi, length.out = 5)) + center_pt[1],
         y = tmp_r * sin(seq(from = 0, to = 2*pi, length.out = 5)) + center_pt[2], type = "n",
         xlab = "Easting (m)", ylab = "Northing (m)", main = paste0("\n Observation ID = ", obs_id), asp = T)
    for(i in 1:length(atm_df$theta_vals[[uid]])){
      segments(x0 = atm_df$obsvr_pts[[uid]][1,i],
               y0 = atm_df$obsvr_pts[[uid]][2,i],
               x1 = atm_df$obsvr_pts[[uid]][1,i] + cos(atm_df$theta_vals[[uid]][i]) * 1e6,
               y1 = atm_df$obsvr_pts[[uid]][2,i] + sin(atm_df$theta_vals[[uid]][i]) * 1e6)
    }
    points(x = atm_df$obsvr_pts[[uid]][1,], y = atm_df$obsvr_pts[[uid]][2,], pch = 21, bg = 2)
    points(x = atm_df$loc_pts[uid,1], y = atm_df$loc_pts[uid,2], pch = 21, bg = 4)
    if(add_legend){
      legend("topright", c("Observer Location", "sigloc Estimate"), pch = c(21,21),
             pt.bg = c("gray35",2), bty = "n")
    }
  }
}

#' @title Diagnostic plots for kappa parameter
#' @description Plotting function to produce diagnostic plots of the kappa parameter MCMC output.
#' @param atm_obj object returned by atm_mcmc()
#' @param item specify which diagnostic to plot
#' @importFrom graphics plot segments
#' @importFrom stats density
#' @export
plot_kappa <- function(atm_obj, item = c("traceplot","density","tuning","run_mean")){
  if(item == "traceplot"){
    plot(x = 1:length(atm_obj$pdraws), y = atm_obj$pdraws, type = "l",
         xlab = "iteration (post burn-in)", ylab = "parameter value", main = "\n kappa diagnostics: traceplot")
  }
  if(item == "density"){
    plot(density(atm_obj$pdraws), xlab = "parameter value", ylab = "density",
         main = "\n kappa diagnostics: posterior density")
    segments(x0 = mean(atm_obj$pdraws), y0 = 0, x1 = mean(atm_obj$pdraws), y1 = 1e6, lty = 2)
  }
  if(item == "tuning"){
    plot(1:length(atm_obj$tuning), atm_obj$tuning, type = "b", xlab = "batch",
         ylab = "tuning parameter value", main = "\n kappa diagnostics: tuning parameter updates",
         pch = 20)
  }
  if(item == "run_mean"){
    plot(x = seq(along = atm_obj$pdraws), y = cumsum(atm_obj$pdraws)/seq(along = atm_obj$pdraws),
         type = "l", xlab = "iteration (post burn-in)", ylab = "running mean",
         main = "\n kappa diagnostic: running mean plot", lwd = 2)
    abline(h = mean(atm_obj$pdraws), col = rgb(0,0,0, alpha = 0.5))
  }
}

#' @title Plot transmitter location posterior credible isopleth
#' @description Plotting function to plot posterior credible isopleth(s) of transmitter location.
#' @param df matrix containing posterior draws
#' @param prob_lvls credible isopleth probability value(s)
#' @param range_extend buffer distance to add to range of df to prevent isopleth discontinuity
#' @param kde_n number of grid points in each direction - scalar n argument from kde2d()
#' @param col_vec vector of length prob_lvls
#' @importFrom grDevices contourLines
#' @importFrom stats approx runif var
#' @export
p_isopleth <- function(df, prob_lvls, range_extend, kde_n, col_vec){

  if(0 == min(apply(df,2,var))){
    dist_tmp <- dist(df)
    dist_tmp <- ifelse(dist_tmp == 0, NA, dist_tmp)
    min_d <- min(dist_tmp, na.rm = TRUE)
    mu_tmp <- matrix(NA, nrow = nrow(df), ncol = 2)
    for(i in 1:nrow(mu_tmp)){
      mu_tmp[i,] <- c(runif(1, min = df[i,1] - min_d/2, max = df[i,1] + min_d/2),
                      runif(1, min = df[i,2] - min_d/2, max = df[i,2] + min_d/2))
    }
    tmp_kernel <- MASS::kde2d(x = mu_tmp[,1], y = mu_tmp[,2], n = kde_n,
                              lims = c(range(mu_tmp[,1]) + c(-1,1)*range_extend,
                                       range(mu_tmp[,2]) + c(-1,1)*range_extend))
  }else{
    tmp_kernel <- MASS::kde2d(x = df[,1], y = df[,2], n = kde_n,
                              lims = c(range(df[,1]) + c(-1,1)*range_extend,
                                       range(df[,2]) + c(-1,1)*range_extend))
  }

  dx <- diff(tmp_kernel$x[1:2])
  dy <- diff(tmp_kernel$y[1:2])
  sz <- sort(tmp_kernel$z)
  c1 <- cumsum(sz) * dx * dy
  lvls <- sapply(prob_lvls, function(x) {
    approx(c1, sz, xout = 1 - x)$y
  })
  out <- contourLines(x = tmp_kernel$x, y = tmp_kernel$y, z = tmp_kernel$z, levels = lvls)

  for(jj in 1:length(out)){
    lines(x = out[[jj]]$x, y = out[[jj]]$y, lwd = 2, col = col_vec[jj])
  }
}
