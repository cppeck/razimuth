#' @title Locate transmitter signal - Lenth's method
#' @description Estimation of transmitter location using azimuthal telemetry data and the algorithms outlined in <http://www.jstor.org/stable/1268030>
#' @param df user supplied dataframe
#' @param max_iter maximum number of iterations
#' @param param_tol tolerance after which algorithm terminates iterations
#' @param robust logical for M-estimation
#' @param psi_fun psi function
#' @param psi_c psi function tuning constant
#' @param buffer buffer distance (m) extending bbox secondary location filter
#' @importFrom stats complete.cases
#' @importFrom utils tail
#' @export
loc_signal <- function(df, max_iter = 100, param_tol = 1e-3, robust = FALSE, psi_fun = c("Huber","Andrews"), psi_c = 1.5, buffer = 100){

  ## pointers ----
  pid <- unique(x = df$obs_id)

  ## containers ----
  obsvr_pts <- theta_vals <- vector(mode = "list", length = length(pid))

  ## deconstruct df into containers ----
  for(i in 1:length(pid)){
    tmp_df <- subset(x = df, subset = df$obs_id == pid[i])
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
  }

  ## construct loc_out ----
  loc_out <- matrix(NA, ncol = 9, nrow = length(pid))
  colnames(loc_out) <- c("obs_id","utm_x","utm_y","var_x","var_y","cov_xy","kappa","iter","error")

  ## implement algorithms specified in Lenth (1981) - <http://www.jstor.org/stable/1268030>
  if(!robust){
    for(i in 1:length(pid)){
      error_code <- 0
      ## initialize loop container ----
      loc_save <- matrix(NA, ncol = 2, nrow = max_iter)

      ## compute location starting value ----
      theta_i <- theta_vals[[i]]
      x_i <- obsvr_pts[[i]][1,]
      y_i <- obsvr_pts[[i]][2,]
      s_i <- sin(theta_i)
      c_i <- cos(theta_i)
      z_i <- s_i*x_i - c_i*y_i
      loc_0 <- solve(a = matrix(c(sum(s_i*s_i), -sum(s_i*c_i),
                                  -sum(c_i*s_i), sum(c_i*c_i)), ncol = 2),
                     b = c(sum(s_i*z_i),
                           -sum(c_i*z_i)))
      loc_save[1,] <- loc_est <- loc_0

      ## iteratively solve linear system given in eq.(2.6) ----
      iter_uid <- 1
      repeat{
        d_i <- sqrt((loc_est[1] - x_i)^2 + (loc_est[2] - y_i)^2)
        s_star <- (loc_est[2] - y_i)/(d_i^3)
        c_star <- (loc_est[1] - x_i)/(d_i^3)

        loc_est <- solve(a = matrix(c(sum(s_i*s_star), -sum(s_i*c_star),
                                      -sum(c_i*s_star), sum(c_i*c_star)), ncol = 2),
                         b = c(sum(s_star*z_i),
                               -sum(c_star*z_i)))

        loc_save[iter_uid + 1,] <- loc_est

        ## break repeat if tolerance reached ----
        if(max(abs(loc_save[iter_uid + 1,] - loc_save[iter_uid,])) < param_tol) {break}

        ## break repeat if max_iter reached ----
        if((iter_uid + 1) == max_iter) {break}

        iter_uid <- iter_uid + 1
      }

      ## compute quantities of interest ----
      loc_tmp <- tail(loc_save[complete.cases(loc_save),], n = 1)
      mu_hat <- atan2(y = (loc_tmp[2] - y_i), x = (loc_tmp[1] - x_i))
      C_bar <- sum(cos(theta_i - mu_hat))/length(theta_i)
      kappa_hat_inv <- 2*(1 - C_bar) + (1 - C_bar)^2 * (0.48794 - 0.82905*C_bar - 1.3915*C_bar^2)/C_bar  ## approx. given in eq.(2.10)
      d_hat_i <- sqrt((loc_tmp[1] - x_i)^2 + (loc_tmp[2] - y_i)^2)
      s_star_hat <- (loc_tmp[2] - y_i)/(d_hat_i^3)
      c_star_hat <- (loc_tmp[1] - x_i)/(d_hat_i^3)
      Q_hat <- kappa_hat_inv * solve(matrix(c(sum(s_star_hat*s_i), (-0.5)*sum(s_star_hat*c_i + c_star_hat*s_i),
                                              (-0.5)*sum(s_star_hat*c_i + c_star_hat*s_i), sum(c_star_hat*c_i)), ncol = 2))  ## approx. given in eq.(2.9)

      ## save loop results ----
      if(kappa_hat_inv <= 0){
        error_code <- 1
        loc_out[i,] <- c(pid[i],rep(NA,6),(iter_uid+1),error_code)
        next()
      }else{
        if((Q_hat[1,1] < 0) + (Q_hat[2,2] < 0) > 0){
          error_code <- 2
          loc_out[i,] <- c(pid[i],rep(NA,6),(iter_uid+1),error_code)
          next()
        }else{
          if((iter_uid + 1) == max_iter){
            error_code <- 3
            loc_out[i,] <- c(pid[i],rep(NA,6),(iter_uid+1),error_code)
            next()
          }else{
            ## bounding box check ----
            box_tmp <- rbind(bearing_int_avg(obs_pts = obsvr_pts[[i]], theta = theta_vals[[i]])$int_mat, t(obsvr_pts[[i]]))
            box_x <- c(min(box_tmp[,1]),max(box_tmp[,1]),max(box_tmp[,1]),min(box_tmp[,1])) + c(-1,1,1,-1)*buffer
            box_y <- c(min(box_tmp[,2]),min(box_tmp[,2]),max(box_tmp[,2]),max(box_tmp[,2])) + c(-1,-1,1,1)*buffer

            if(mgcv::in.out(bnd = cbind(box_x,box_y), x = loc_tmp)){
              loc_out[i,] <- c(pid[i],loc_tmp,Q_hat[1,1],Q_hat[2,2],Q_hat[1,2],(1/kappa_hat_inv),(iter_uid+1),error_code)
            }else{
              error_code <- 4
              loc_out[i,] <- c(pid[i],rep(NA,6),(iter_uid+1),error_code)
            }
          }
        }
      }
    }
  }

  if(robust){
    if(!(psi_fun %in% c("Huber","Andrews"))) stop("'psi_fun' option not available; available options: 'Huber' or 'Andrews'")
    weights_ls <- vector("list", length(pid))
    for(i in 1:length(pid)){
      error_code <- 0
      ## initialize loop containers ----
      loc_save <- matrix(NA, ncol = 2, nrow = max_iter)
      w_save <- matrix(NA, ncol = length(theta_vals[[i]]), nrow = max_iter)

      ## compute location starting value and initial kappa ----
      theta_i <- theta_vals[[i]]
      x_i <- obsvr_pts[[i]][1,]
      y_i <- obsvr_pts[[i]][2,]
      s_i <- sin(theta_i)
      c_i <- cos(theta_i)
      z_i <- s_i*x_i - c_i*y_i
      loc_0 <- solve(a = matrix(c(sum(s_i*s_i), -sum(s_i*c_i),
                                  -sum(c_i*s_i), sum(c_i*c_i)), ncol = 2),
                     b = c(sum(s_i*z_i),
                           -sum(c_i*z_i)))

      loc_save[1,] <- loc_est <- loc_0
      w_save[1,] <- rep(1, length(theta_i))

      mu_hat <- atan2(y = (loc_0[2] - y_i), x = (loc_0[1] - x_i))
      Cw_bar <- sum(cos(theta_i - mu_hat))/length(theta_i)  ## w_i = 1
      wkappa_hat_inv <- 2*(1 - Cw_bar) + (1 - Cw_bar)^2 * (0.48794 - 0.82905*Cw_bar - 1.3915*Cw_bar^2)/Cw_bar
      if(wkappa_hat_inv <= 0){
        error_code <- 1
        loc_out[i,] <- c(pid[i],rep(NA,6),1,error_code)
        next()
      }

      ## iteratively solve linear system given in Equation (4.7) ----
      iter_uid <- 1
      repeat{
        t_fun <- sqrt(2*(1/wkappa_hat_inv)*(1 - cos(theta_i - mu_hat)))  ## standardization function - eq.(4.2)
        # t_fun <- t_fun * ifelse(t_fun %% 2*pi >= 0 & t_fun %% 2*pi <= pi, 1, -1)  ## eq.(3.5) from Lenth (1981) <http://www.jstor.org/stable/1267979>
        if(psi_fun == "Huber"){
          psi_huber <- as.numeric(apply(matrix(t_fun, ncol = 1), 1, function(x) sign(x) * min(abs(x), psi_c)))  ## psi function - eq.(4.5)
          w_i <- ifelse(t_fun == 0, 1, psi_huber/t_fun)
        } else{
          if(psi_fun == "Andrews"){
            psi_andrews <- psi_c * sin(t_fun/psi_c) * ifelse(abs(t_fun) < psi_c * pi, 1, 0)  ## psi function - eq.(4.6)
            w_i <- ifelse(t_fun == 0, 1, psi_andrews/t_fun)
          }
        }

        d_i <- sqrt((loc_est[1] - x_i)^2 + (loc_est[2] - y_i)^2)
        s_star <- (loc_est[2] - y_i)/(d_i^3)
        c_star <- (loc_est[1] - x_i)/(d_i^3)

        loc_est <- solve(a = matrix(c(sum(w_i*(s_i*s_star)), -sum(w_i*(s_i*c_star)),
                                      -sum(w_i*(c_i*s_star)), sum(w_i*(c_i*c_star))), ncol = 2),
                         b = c(sum(w_i*(s_star*z_i)),
                               -sum(w_i*(c_star*z_i))))

        loc_save[iter_uid + 1,] <- loc_est
        w_save[iter_uid + 1,] <- w_i

        ## break repeat if tolerance reached for both w_i and location estimates ----
        if(max(abs(w_save[iter_uid + 1,] - w_save[iter_uid,])) < param_tol &
           max(abs(loc_save[iter_uid + 1,] - loc_save[iter_uid,])) < param_tol) {break}

        ## break repeat if max_iter reached ----
        if((iter_uid + 1) == max_iter) {break}

        iter_uid <- iter_uid + 1
        mu_hat <- atan2(y = (loc_est[2] - y_i), x = (loc_est[1] - x_i))
        Cw_bar <- sum(w_i*(cos(theta_i - mu_hat)))/sum(w_i)
        wkappa_hat_inv <- 2*(1 - Cw_bar) + (1 - Cw_bar)^2 * (0.48794 - 0.82905*Cw_bar - 1.3915*Cw_bar^2)/Cw_bar
        if(wkappa_hat_inv <= 0){break}
      }

      ## compute quantities of interest ----
      loc_tmp <- tail(loc_save[complete.cases(loc_save),], n = 1)
      mu_hat <- atan2(y = (loc_tmp[2] - y_i), x = (loc_tmp[1] - x_i))
      Cw_bar <- sum(w_i*(cos(theta_i - mu_hat)))/sum(w_i)
      wkappa_hat_inv <- 2*(1 - Cw_bar) + (1 - Cw_bar)^2 * (0.48794 - 0.82905*Cw_bar - 1.3915*Cw_bar^2)/Cw_bar

      d_hat_i <- sqrt((loc_tmp[1] - x_i)^2 + (loc_tmp[2] - y_i)^2)
      s_star_hat <- (loc_tmp[2] - y_i)/(d_hat_i^3)
      c_star_hat <- (loc_tmp[1] - x_i)/(d_hat_i^3)
      Q_hat <- wkappa_hat_inv * solve(matrix(c(sum(w_i*(s_star_hat*s_i)), (-0.5)*sum(w_i*(s_star_hat*c_i + c_star_hat*s_i)),
                                              (-0.5)*sum(w_i*(s_star_hat*c_i + c_star_hat*s_i)), sum(w_i*(c_star_hat*c_i))), ncol = 2))  ## approx. given in eq.(2.9)

      ## save loop results ----
      if(wkappa_hat_inv <= 0){
        error_code <- 1
        loc_out[i,] <- c(pid[i],rep(NA,6),(iter_uid+1),error_code)
        next()
      }else{
        if((Q_hat[1,1] < 0) + (Q_hat[2,2] < 0) > 0){
          error_code <- 2
          loc_out[i,] <- c(pid[i],rep(NA,6),(iter_uid+1),error_code)
          next()
        }else{
          if((iter_uid + 1) == max_iter){
            error_code <- 3
            loc_out[i,] <- c(pid[i],rep(NA,6),(iter_uid+1),error_code)
            next()
          }else{
            ## bounding box check ----
            box_tmp <- rbind(bearing_int_avg(obs_pts = obsvr_pts[[i]], theta = theta_vals[[i]])$int_mat, t(obsvr_pts[[i]]))
            box_x <- c(min(box_tmp[,1]),max(box_tmp[,1]),max(box_tmp[,1]),min(box_tmp[,1])) + c(-1,1,1,-1)*buffer
            box_y <- c(min(box_tmp[,2]),min(box_tmp[,2]),max(box_tmp[,2]),max(box_tmp[,2])) + c(-1,-1,1,1)*buffer

            if(mgcv::in.out(bnd = cbind(box_x,box_y), x = loc_tmp)){
              loc_out[i,] <- c(pid[i],loc_tmp,Q_hat[1,1],Q_hat[2,2],Q_hat[1,2],(1/wkappa_hat_inv),(iter_uid+1),error_code)
            }else{
              error_code <- 4
              loc_out[i,] <- c(pid[i],rep(NA,6),(iter_uid+1),error_code)
            }
          }
        }
      }
      weights_ls[[i]] <- tail(w_save[complete.cases(w_save),], n = 1)
    }
  }

  ## write output ----
  list(loc_out = loc_out,
       max_iter = max_iter,
       param_tol = param_tol,
       buffer = buffer,
       psi_fun = if(robust){psi_fun} else{NULL},
       psi_c = if(robust){psi_c} else{NULL},
       weights_ls = if(robust){weights_ls} else{NULL})
}
