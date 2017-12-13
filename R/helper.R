#' @title Generate ATM dataframe
#' @description Generate objects in format usable by razimuth package functions.
#' @param df user supplied dataframe
#' @export
convert_atm <- function(df){

  ## pointers ----
  pid <- unique(df$obs_id)

  ## containers ----
  obsvr_pts <- list()
  theta_vals <- list()
  radius_ls <- list()
  dt_ls <- list()

  ## deconstruct df into containers for use in atm_mcmc() ----
  for(i in 1:length(pid)){
    tmp_df <- subset(x = df,
                     subset = df$obs_id == pid[i])
    obsvr_pts[[i]] <- matrix(data = t(tmp_df[,c("utm_x","utm_y")]),
                             nrow = 2,
                             byrow = F)
    tmp_vec <- circular::circular(x = tmp_df[,"azimuth"],
                                  units = "degrees",
                                  zero = pi/2,
                                  rotation = "clock")
    theta_vals[[i]] <- as.numeric(circular::conversion.circular(x = tmp_vec,
                                                                units = "radians",
                                                                zero = 0,
                                                                rotation = "counter"))
    radius_ls[[i]] <- tmp_df[,"prior_r"]
    dt_ls[[i]] <- tmp_df[,c("date")]
  }

  ## grid anchor points (based on pairwise average of intersections) ----
  loc_pts <- matrix(NA, nrow = length(pid), ncol = 2)
  for(i in 1:length(pid)){
    loc_pts[i,] <- bearing_int_avg(obs_pts = obsvr_pts[[i]], theta = theta_vals[[i]])$int_avg
  }
  # loc_pts <- loc_signal(df = df)$loc_out[,c("utm_x","utm_y")]

  ## write output ----
  list(obsvr_pts = obsvr_pts,
       theta_vals = theta_vals,
       radius_ls = radius_ls,
       dt_ls = dt_ls,
       loc_pts = loc_pts,
       pid = pid)
}

#' @title Modified locate()
#' @description Modified version of locate() from the R package sigloc.
#' @param x user supplied dataframe
#' @param int_avg average intersection
#' @param int_mat matrix containing all unique intersection points
locate_est <- function(x, int_avg, int_mat){

  ## Define function to solve system of MLE equations
  solver <- function(par){

    ## Isolate portion of data corresponding to the unique grouping
    subdata <- data.frame(X = x$utm_x[x$obs_id == group],
                          Y = x$utm_y[x$obs_id == group],
                          SIN = x$sin[x$obs_id == group],
                          COS = x$cos[x$obs_id == group])

    ## Create variables to hold transmitter location estimates
    loc_est <- numeric(length(par))
    transmitter_x <- par[1]
    transmitter_y <- par[2]

    ## Define system of MLE equations (Lenth 1981)
    loc_est[1] = -sum((transmitter_y-subdata$Y)*(subdata$SIN*(transmitter_x-subdata$X)-subdata$COS*(transmitter_y-subdata$Y))/(((transmitter_x-subdata$X)^2+(transmitter_y-subdata$Y)^2)^0.5)^3)
    loc_est[2] = sum((transmitter_x-subdata$X)*(subdata$SIN*(transmitter_x-subdata$X)-subdata$COS*(transmitter_y-subdata$Y))/(((transmitter_x-subdata$X)^2+(transmitter_y-subdata$Y)^2)^0.5)^3)

    ## Return transmitter location estimate
    loc_est
  }

  ## Convert compass bearing to standard radian measure (Lenth 1981)
  x$theta <- (pi/180*(90-x$azimuth))

  ## Calculate necessary variables
  x$sin = sin(x$theta)
  x$cos = cos(x$theta)
  x$tan = tan(x$theta)

  ## Create data frame to store calculated transmitter locations
  transmitter <- data.frame(X = NA, Y = NA, BadPoint = NA)

  ## Declare buffer (in meters) for ensuring location validity
  buffer <- 5

  ## Calculate transmitter location for each grouping
  for(group in 1:length(unique(x$obs_id))){

    ## Calculate starting point for solver function based on the mean of bearing intersections
    start_point_x <- int_avg[group, 2]
    start_point_y <- int_avg[group, 3]

    ## Temporarily store results of solver function
    triangulation_results <- nleqslv::nleqslv(c(start_point_x, start_point_y), solver)

    ## Withdraws transmitter location from list of solver function results
    location<-triangulation_results$x

    ## Ensure location validity prior to recording
    valid <- TRUE
    if(location[1] < min(int_mat[[group]][,1] - buffer)) valid <- FALSE
    if(location[1] > max(int_mat[[group]][,1] + buffer)) valid <- FALSE
    if(location[2] < min(int_mat[[group]][,2] - buffer)) valid <- FALSE
    if(location[2] > max(int_mat[[group]][,2] + buffer)) valid <- FALSE

    ## Set default color scheme for transmitter ploting
    transmitter[group,3] <- 0

    ## Transfer calculated locations into results variable
    if(valid) {
      transmitter[group,1] <- location[1]
      transmitter[group,2] <- location[2]
    } else {
      # warning("Bad point detected in Grouping ", group)
      transmitter[group,1] <- start_point_x
      transmitter[group,2] <- start_point_y
      transmitter[group,3] <- 1
    }
  }

  ## Return transmitter locations
  return(transmitter)
}

#' @title pid identifier helper
#' @description Function to identify pid for a given obs_id
#' @param atm_df object returned by convert_atm()
#' @param obs_id observation id
#' @export
which_pid <- function(atm_df, obs_id){
  return(which(unique(atm_df$pid) == obs_id))
}

#' @title Compute average intersection
#' @description Function that computes average intersection of bearings
#' @param obs_pts observer locations
#' @param theta bearing taken from each observer location
bearing_int_avg <- function(obs_pts, theta){
  ## pointers ----
  npt <- length(theta)

  ## initialize containers ----
  int_mat <- matrix(NA, ncol = 2, nrow = factorial(npt))
  keep <- as.numeric()

  for(i in 1:npt){
    theta[i] <- ifelse(theta[i] < -pi, theta[i] + 2*pi, ifelse(theta[i] > pi, theta[i] - 2*pi, theta[i]))
  }

  ## loop
  ct <- 1
  for(i in 1:npt){
    for(j in 1:npt){
      if(i != j){
        ## compute intersections
        int_mat[ct,2] <- ((obs_pts[1,j] - obs_pts[1,i])*tan(theta[i])*tan(theta[j]) - obs_pts[2,j]*tan(theta[i]) + obs_pts[2,i]*tan(theta[j]))/(tan(theta[j]) - tan(theta[i]))
        int_mat[ct,1] <- (int_mat[ct,2] - obs_pts[2,j])/tan(theta[j]) + obs_pts[1,j]

        ## direction of intersection indicator
        theta_i <- atan2(int_mat[ct,2] - obs_pts[2,i], int_mat[ct,1] - obs_pts[1,i])
        theta_j <- atan2(int_mat[ct,2] - obs_pts[2,j], int_mat[ct,1] - obs_pts[1,j])
        keep[ct] <- ifelse(sum(round(c(theta_i - theta[i], theta_j - theta[j]), digits = 4)) == 0, 1, 0)

        ct <- ct + 1
      }
    }
  }

  ui <- which(keep == 1)
  if(sum(ui) != 0){
    int_mat <- int_mat[ui,]
    int_avg <- apply(int_mat, 2, mean)
  } else{
    int_avg <- apply(t(obs_pts), 2, mean)
  }

  list(int_avg = int_avg,
       int_mat = int_mat)
}
