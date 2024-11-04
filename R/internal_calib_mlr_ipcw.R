### Internal functions to allow estimation of calibration using the MLR-IPCW method.

#' Estimate data for calibration plots using MLR-IPCW.
#' @description
#' Called in calib_msm::calib_msm to apply the MLR-IPCW method.
#'
#' @details
#' Observed event probabilities are already returned only for the landmarked cohort of
#' individuals also uncesored at time ``t`. There is no option to plot over a different
#' set of predicted transition probabilities (unlike in `calib_blr` and `calib_pv`) as
#' that option is provided purely to allow bootstrapping, which makes no sense with the
#' calibration scatter plots produced by MLR-IPCW.
#'
#' @returns A list of datasets for each calibration plot.
#'
#' @noRd
calib_mlr_ipcw <- function(data_raw,
                           data_ms,
                           j,
                           s,
                           t,
                           weights_provided,
                           w_covs,
                           w_function,
                           w_landmark_type,
                           w_max,
                           w_stabilised,
                           w_max_follow,
                           mlr_smoother_type,
                           mlr_ps_int,
                           mlr_degree,
                           mlr_s_df,
                           mlr_niknots,
                           CI,
                           CI_R_boot,
                           CI_seed,
                           valid_transitions,
                           assess_moderate,
                           assess_mean, ...){

  ## Create landmarked dataset
  data_raw_lmk_js_uncens <-  apply_landmark(data_raw = data_raw, data_ms = data_ms, j = j, s = s, t = t, exclude_cens_t = TRUE)

  ## Calculate weights
  ## Note this is done in the entire dataset data_boot, which has its own functionality (w_landmark_type) to landmark on j and s, or just s, before
  ## calculating the weights
  if (weights_provided == FALSE){

    ## If custom function for estimating weights has been inputted ("w_function"), replace "calc_weights" with this function
    if (!is.null(w_function)){
      calc_weights <- w_function
    }

    ## Calculate the weights
    weights <- calc_weights(data_ms = data_ms,
                            data_raw = data_raw,
                            covs = w_covs,
                            t = t,
                            s = s,
                            landmark_type = w_landmark_type,
                            j = j,
                            max_weight = w_max,
                            stabilised = w_stabilised,
                            max_follow = w_max_follow,
                            ...)

    ## Add weights to data_boot
    data_raw_lmk_js_uncens <- dplyr::left_join(data_raw_lmk_js_uncens, dplyr::distinct(weights), by = dplyr::join_by(id))

  }

  ### Assign reference category
  ref_cat <- paste(valid_transitions[1])

  ###
  ### Calibration plots/Moderate Calibration
  ###
  if (assess_moderate == TRUE){

    ### Define equation
    eq_LHS <- paste("state_poly_fac ~ ")
    if (mlr_smoother_type == "s"){
      eq_RHS <- paste("s(mlr_lp", 1:(length(valid_transitions) - 1), ", df = mlr_s_df)", sep = "", collapse = "+")
    } else if (mlr_smoother_type == "sm.ps"){
      eq_RHS <- paste("sm.ps(mlr_lp", 1:(length(valid_transitions) - 1), ", ps.int = mlr_ps_int, degree = mlr_degree)", sep = "", collapse = "+")
    } else if (mlr_smoother_type == "sm.os"){
      eq_RHS <- paste("sm.os(mlr_lp", 1:(length(valid_transitions) - 1), ", niknots = mlr_niknots)", sep = "", collapse = "+")
    }
    eq_mlr <- stats::as.formula(paste(eq_LHS, eq_RHS, sep =""))

    ### Apply nominal recalibration framework of van Hoorde et al., (2014)
    calib_model <- VGAM::vgam(eq_mlr, weights = data_raw_lmk_js_uncens[, "ipcw"],
                              data = data_raw_lmk_js_uncens, family = VGAM::multinomial(refLevel = ref_cat))

    ###
    ### Generate predicted-observed risks and add to data_raw_lmk_js_uncens

    ## Create dataframe to store
    dat_mlr_pred_obs <- data.frame(matrix(NA, ncol = length(valid_transitions), nrow = nrow(data_raw_lmk_js_uncens)))
    ## Assign colnames
    colnames(dat_mlr_pred_obs) <- paste("mlr_pred_obs", valid_transitions, sep = "")
    ## Calc pred_obs for those who are uncensored
    mlr_pred_obs <- VGAM::predictvglm(calib_model, newdata = data_raw_lmk_js_uncens, type = "response")
    ## Assign to appropriate individuals
    dat_mlr_pred_obs[, ] <- mlr_pred_obs

    ### Then add it to data_raw_lmk_js_uncens
    data_raw_lmk_js_uncens <- cbind(data_raw_lmk_js_uncens, dat_mlr_pred_obs)

    ### Assign output
    output_object <- dplyr::select(data_raw_lmk_js_uncens, id, paste("tp_pred", valid_transitions, sep = ""), paste("mlr_pred_obs", valid_transitions, sep = ""))

    ### Get plotdata in same format as calib_blr
    ## Start by creating new output object
    output_object_plots <- vector("list", length(valid_transitions))
    names(output_object_plots) <- paste("state", valid_transitions, sep = "")

    ## Loop through and create output object for each valid transition (same format as for BLR-IPCW and PV)
    for (k in 1:length(valid_transitions)){

      ## Assign state of interest
      state_k <- valid_transitions[k]

      ## Create output object
      output_object_plots[[k]] <- data.frame("id" = output_object[, "id"],
                                             "pred" = output_object[, paste("tp_pred", valid_transitions[k], sep = "")],
                                             "obs" = output_object[, paste("mlr_pred_obs", valid_transitions[k], sep = "")])
    }

  }

  ###
  ### Calibration intercept
  ###
  if (assess_mean == TRUE){

    if (CI != FALSE){

      ### Create object to store plot data
      output_object_mean <- vector("list", length(valid_transitions))
      names(output_object_mean) <- paste("state", valid_transitions, sep = "")

      ### Define alpha for CI's
      alpha <- (1-CI/100)/2

      ### Set seed for bootstrapping
      if (!is.null(CI_seed)){
        set.seed(CI_seed)
      }

      ### Apply bootstrapping
      boot_mean <- boot::boot(data_raw, calc_mean_mlr_boot, R = CI_R_boot,
                              data_ms = data_ms,
                              state_k = state_k,
                              j = j,
                              s2 = s,
                              t = t,
                              weights_provided = weights_provided,
                              w_function = w_function,
                              w_covs = w_covs,
                              w_landmark_type = w_landmark_type,
                              w_max = w_max,
                              w_stabilised = w_stabilised,
                              w_max_follow = w_max_follow,
                              valid_transitions = valid_transitions, ...)

      ### Extract confidence bands
      lower <- apply(boot_mean$t, 2, stats::quantile, probs = alpha, na.rm = TRUE)
      upper <- apply(boot_mean$t, 2, stats::quantile, probs = 1-alpha, na.rm = TRUE)

      ### Cycle through states and assign to correct part of output object
      for (state in 1:length(valid_transitions)){

        ### Put into output object
        output_object_mean[[state]] <- c("mean" = as.numeric(boot_mean$t0[state]),
                                         "mean_lower" = as.numeric(lower[state]),
                                         "mean_upper" = as.numeric(upper[state]))

      }

    } else if (CI == FALSE){

      ### Assess mean calibration
      output_object_mean <- calc_mean_mlr_boot(data_raw = data_raw,
                                               indices = 1:nrow(data_raw),
                                               data_ms = data_ms,
                                               state_k = state_k,
                                               j = j,
                                               s2 = s,
                                               t = t,
                                               weights_provided = weights_provided,
                                               w_function = w_function,
                                               w_covs = w_covs,
                                               w_landmark_type = w_landmark_type,
                                               w_max = w_max,
                                               w_stabilised = w_stabilised,
                                               w_max_follow = w_max_follow,
                                               valid_transitions = valid_transitions, ...)

      ### Put into output object
      names(output_object_mean) <- paste("state", valid_transitions, sep = "")

    }

  }

  # ###
  # ### PLACEHOLDER FOR Calibration slope
  # ###
  #
  # ## Add constraints
  # i <- diag(4)
  # i1 <- rbind(1, 0, 0, 0)
  # i2 <- rbind(0, 1, 0, 0)
  # i3 <- rbind(0, 0, 1, 0)
  # i4 <- rbind(0, 0, 0, 1)
  # clist <- list("(Intercept)" = i, "mlr_lp1" = i1, "mlr_lp2" = i2, "mlr_lp3" = i3, "mlr_lp4" = i4)
  # clist
  #
  # ### Apply nominal recalibration framework
  # ###
  #
  # ### Define equation
  # eq_LHS <- paste("state_poly_fac ~ ")
  # eq_RHS <- paste("mlr_lp", 1:(length(valid_transitions) - 1), sep = "", collapse = "+")
  # eq_mlr <- stats::as.formula(paste(eq_LHS, eq_RHS, sep =""))
  #
  # ### Fit model to estimate slopes
  # mlr_slope_model <- VGAM::vgam(state_poly_fac ~ mlr_lp1 + mlr_lp2 + mlr_lp3 + mlr_lp4, constraints = clist, weights = ipcw,
  #                           data = data_raw_lmk_js_uncens, family = VGAM::multinomial(refLevel = ref_cat))
  #
  # ### Extract slopes and standard errors
  # mlr_slopes <- mlr_slope_model@coefficients[paste("mlr_lp", 1:(length(valid_transitions) - 1), sep = "")]
  # mlr_slopes_se <- sqrt(diag(stats::vcov(mlr_slope_model))[paste("mlr_lp", 1:(length(valid_transitions) - 1), sep = "")])
  #
  # ### Define output object
  # output_object_weak <- list("slopes" = mlr_slopes, "slopes_se" = mlr_slopes_se)

  ### Define combined output object
  output_object_comb <- vector("list")

  if (assess_moderate == TRUE){
    output_object_comb[["plotdata"]] <- output_object_plots
  }
  if (assess_mean == TRUE){
    output_object_comb[["mean"]] <- output_object_mean
  }
  # if (assess_weak == TRUE){
  #   output_object_comb[["weak"]] <- output_object_weak
  # }

  return(output_object_comb)

}



#' Estimate mean calibration using MLR-IPCW.
#' @description
#' Estimate mean calibration using MLR-IPCW.
#'
#' @details
#' Function written in a format so that it can be used in combination with \code{\link[boot]{boot}}
#' for bootstrapping. Specifying `indices = 1:nrow(data_raw)` will return the calibration
#' of interest.
#'
#' @returns Mean calibration and calibration slope for a given state.
#'
#' @noRd
calc_mean_mlr_boot <- function(data_raw,
                               indices,
                               data_ms,
                               state_k,
                               j,
                               s2, # can't use 's' because it matches an argument for the boot function
                               t,
                               weights_provided,
                               w_function,
                               w_covs,
                               w_landmark_type,
                               w_max,
                               w_stabilised,
                               w_max_follow,
                               valid_transitions, ...){

  ## Create bootstrapped dataset
  data_boot <- data_raw[indices, ]

  ## Create landmarked dataset
  data_boot_lmk_js_uncens <-  apply_landmark(data_raw = data_boot, data_ms = data_ms, j = j, s = s2, t = t, exclude_cens_t = TRUE)

  ## Calculate weights
  ## Note this is done in the entire dataset data_boot, which has its own functionality (w_landmark_type) to landmark on j and s, or just s, before
  ## calculating the weights
  if (weights_provided == FALSE){

    ## If custom function for estimating weights has been inputted ("w_function"), replace "calc_weights" with this function
    if (!is.null(w_function)){
      calc_weights <- w_function
    }

    ## Calculate the weights
    weights <- calc_weights(data_ms = data_ms,
                            data_raw = data_boot,
                            covs = w_covs,
                            t = t,
                            s = s2,
                            landmark_type = w_landmark_type,
                            j = j,
                            max_weight = w_max,
                            stabilised = w_stabilised,
                            max_follow = w_max_follow,
                            ...)

    ## Add weights to data_boot
    data_boot_lmk_js_uncens <- dplyr::left_join(data_boot_lmk_js_uncens, dplyr::distinct(weights), by = dplyr::join_by(id))

  }

  ### Assign reference category
  ref_cat <- paste(valid_transitions[1])

  ### Define equation
  eq_LHS <- paste("state_poly_fac ~ ")
  eq_RHS <- paste("mlr_lp", 1:(length(valid_transitions) - 1), sep = "", collapse = "+")
  eq_mlr <- stats::as.formula(paste(eq_LHS, eq_RHS, sep =""))

  ### Fit model to estimate intercepts
  mlr_int_model <- VGAM::vgam(data_boot_lmk_js_uncens$state_poly_fac ~ 1,
                              offset = as.matrix(data_boot_lmk_js_uncens[, paste("mlr_lp", 1:(length(valid_transitions) - 1), sep = "")]),
                              weights = data_boot_lmk_js_uncens$ipcw,
                              family = VGAM::multinomial(refLevel = ref_cat))

  ###  Now need to calculate predicted-observed values for each state using this model, and calculate mean difference from predicted risks

  ### Create a new dataset with linear predictors from the model (intercept plus offset)
  data_pred_obs_lp <- matrix( nrow = nrow(data_boot_lmk_js_uncens), ncol = length(valid_transitions) - 1)
  for (state in 1:(length(valid_transitions) - 1)){
    data_pred_obs_lp[,state] <-
      exp(stats::coefficients(mlr_int_model)[state] + data_boot_lmk_js_uncens[,paste("mlr_lp", state, sep = "")])
  }

  ### Now calculate the predicted probabilities
  p1 <- 1/(1 + rowSums(data_pred_obs_lp))
  p_rest <- data_pred_obs_lp/(1 + rowSums(data_pred_obs_lp))

  ### Put together into dataset of predicted-observed values from the intercept model
  int_pred_obs <- data.frame(cbind(p1, p_rest))
  colnames(int_pred_obs) <- paste("pred_obs_", valid_transitions, sep = "")

  ### Check rows sum to 1 as expected
  # rowSums(int_pred_obs)

  ### Get difference between predicted-observed and predicted risks
  int_diff <- colMeans(int_pred_obs - data_boot_lmk_js_uncens[,paste("tp_pred", valid_transitions, sep = "")])
  names(int_diff) <- paste("state", valid_transitions, sep = "")

  return(int_diff)

}
