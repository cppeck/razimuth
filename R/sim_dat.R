#' @title Simulate azimuthal telemetry data
#' @description Simulate azimuthal telemetry data under random design for demonstration purposes
#' @param n_loc number of transmitter locations to simulate uniformly in square of area sq_km
#' @param sq_km area of square spatial domain in km
#' @param n_azimuth number of azimuths taken to triangulate transmitter location
#' @param dist_vec min and max distance between observer and transmitter location (drawn uniformly)
#' @param kappa_sim value for kappa governing uncertainty in which we observe azimuths
#' @param prior_r prior radius
#' @param plot logical on whether or not to plot each simulated transmitter location with observed azimuths
#' @importFrom graphics par
#' @export
sim_atd <- function(n_loc = 100, sq_km = 4, n_azimuth = 3, dist_vec = c(200,500), kappa_sim = 50, prior_r = 1000, plot = FALSE){

  ## simulate true locations ----
  tmp_r <- (sq_km/4)*1000
  x_vals <- runif(n = n_loc, min = 493392 - tmp_r, max = 493392 + tmp_r)
  y_vals <- runif(n = n_loc, min = 4489825 - tmp_r, max = 4489825 + tmp_r)
  spt_mat <- cbind(x_vals, y_vals)

  ## number of bearings per location ----
  jList <- rep(x = n_azimuth, times = n_loc)

  ## azimuths ----
  theta_ls <- vector("list",n_loc)
  for(i in 1:n_loc){
    theta_ls[[i]] <- runif(n = n_azimuth, min = 0, max = 2*pi)
  }

  ## distance between transmitter and observer ----
  d_ls <- vector("list",n_loc)
  for(i in 1:n_loc){
    d_ls[[i]] <- runif(n = n_azimuth, min = dist_vec[1], max = dist_vec[2])
  }

  ## observer locations ----
  tpt_ls <- vector("list",n_loc)
  for(i in 1:n_loc){
    tpt_ls[[i]] <- matrix(NA, nrow = 2, ncol = jList[i])
    for(j in 1:jList[i]){
      tpt_ls[[i]][,j] <- c(spt_mat[i,1] + cos(theta_ls[[i]][j]) * d_ls[[i]][j],
                           spt_mat[i,2] + sin(theta_ls[[i]][j]) * d_ls[[i]][j])
    }
  }

  ## azimuth uncertainty ----
  kappa.val <- kappa_sim
  theta_adj <- theta_ls
  for(i in 1:n_loc){
    for(j in 1:length(theta_ls[[i]])){
      theta_adj[[i]][j] <- CircStats::rvm(n = 1, mean = (theta_ls[[i]][j] - pi), k = kappa.val)
    }
  }

  ## construct df ----
  z.tmp <- matrix(unlist(tpt_ls), ncol = 2, byrow = T)
  df <- cbind(z.tmp, unlist(theta_adj))
  df[,3] <- abs((df[,3] * (180/pi)) - 360) + 90
  df[,3] <- ifelse(df[,3] > 360, df[,3] - 360, df[,3])
  tmp <- rep(1, jList[1])
  for(i in 2:n_loc){
    tmp <- c(tmp, rep(i,jList[i]))
  }
  df <- cbind("sim_hen", tmp, df)
  df <- data.frame(df)
  df$date <- as.POSIXct(Sys.Date(), format = "%Y-%m-%d", tz = "UTC")
  df$prior_r <- prior_r
  colnames(df) <- c("indiv","obs_id","utm_x","utm_y","azimuth","date","prior_r")

  df$utm_x <- as.numeric(as.character(df$utm_x))
  df$utm_y <- as.numeric(as.character(df$utm_y))
  df$azimuth <- round(as.numeric(as.character(df$azimuth)))
  df$obs_id <- as.integer(as.character(df$obs_id))
  df$prior_r <- as.numeric(as.character(df$prior_r))
  df$date <- as.POSIXct(df$date)

  if(plot){
    ## visual check ----
    par(mfrow = c(2,2))
    for(i in 1:n_loc){
      tmp <- t(tpt_ls[[i]])
      ext_val <- dist_vec[2]*2
      plot(spt_mat[i,1], spt_mat[i,2], xlim = c(spt_mat[i,1] - ext_val, spt_mat[i,1] + ext_val),
           ylim = c(spt_mat[i,2] - ext_val, spt_mat[i,2] + ext_val), pch = 20, col = 2, main = paste0("i = ",i),
           xlab = "Easting (m)", ylab = "Northing (m)", asp = TRUE)
      points(tmp, pch = 20)
      for(j in 1:nrow(tmp)){
        segments(x0 = tmp[j,1], y0 = tmp[j,2], x1 = tmp[j,1] + cos(theta_adj[[i]][j])*1e10,
                 y1 = tmp[j,2] + sin(theta_adj[[i]][j])*1e10)
      }
      legend("topright", c("observer location(s)","true transmitter location"), pch = c(20,20), col = c(1,2),
             bty = "n")
    }
    par(mfrow = c(1,1))
  }

  ## return df ----
  list(sim_df = df)
}
