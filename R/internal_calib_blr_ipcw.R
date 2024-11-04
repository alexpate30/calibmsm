### Internal functions to allow estimation of calibration using the BLR-IPCW method.

#' Estimate data for calibration plots using BLR-IPCW.
#' @description
#' Called in calibmsm::calib_msm to apply the BLR-IPCW method.
#'
#' @details
#' Calls heavily on calc_obs_blr_loess_boot or calc_obs_blr_rcs_boot to estimate the observed
#' event probabilities for a given set of predicted transition probabilities. Bootstrapping
#' may be applied depending on user input.
#'
#' @returns A list of datasets for each calibration plot.
#'
#' @noRd
calib_blr_ipcw <- function(data_raw,
                           data_ms,
                           tp_pred_plot,
                           j,
                           s,
                           t,
                           curve_type,
                           rcs_nk,
                           loess_span,
                           loess_degree,
                           loess_surface,
                           loess_trace_hat,
                           loess_cell,
                           loess_iterations,
                           loess_iterTrace,
                           weights_provided,
                           w_function,
                           w_covs,
                           w_landmark_type,
                           w_max,
                           w_stabilised,
                           w_max_follow,
                           CI,
                           CI_type,
                           CI_R_boot,
                           CI_seed,
                           transitions_out,
                           assess_moderate,
                           assess_mean, ...){

  ###
  ### Create the object of predicted risks over which the calibration plots will be plotted

  ### For calib_blr_ipcw, this is the landmarked cohort of individuals who are also uncensored at time t,
  ### or tp_pred_plot if specified
  if (is.null(tp_pred_plot)){
    ## Note that tp_pred has been added to data_raw and the predicted transition probabilities (and relevant transformations)
    ## are contained within this dataset.
    data_to_plot <- apply_landmark(data_raw = data_raw, data_ms = data_ms, j = j, s = s, t = t, exclude_cens_t = TRUE)
  } else if (!is.null(tp_pred_plot)){
    data_to_plot <- tp_pred_plot
  }

  ### Create object to store output
  output_object_plots <- vector("list", length(transitions_out))
  names(output_object_plots) <- paste("state", transitions_out, sep = "")

  output_object_mean <- vector("list", length(transitions_out))
  names(output_object_mean) <- paste("state", transitions_out, sep = "")

  # output_object_weak <- vector("list", length(transitions_out))
  # names(output_object_weak) <- paste("state", transitions_out, sep = "")

  ### Loop through and fit models if BLR selected
  for (k in 1:length(transitions_out)){

    ### Assign state of interest
    state_k <- transitions_out[k]

    ### Calculate confidence intervals
    if (CI != FALSE & CI_type == "bootstrap"){

      ### Define alpha for CI's
      alpha <- (1-CI/100)/2

      ### Set seed for bootstrapping
      if (!is.null(CI_seed)){
        set.seed(CI_seed)
      }

      ###
      ### Calibration plots
      ###
      if (assess_moderate == TRUE){
        if (curve_type == "loess"){
          boot_obs <- boot::boot(data_raw, calc_obs_blr_loess_boot, R = CI_R_boot,
                                 data_ms = data_ms,
                                 data_to_plot = data_to_plot,
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
                                 loess_span = loess_span,
                                 loess_degree = loess_degree,
                                 loess_surface = loess_surface,
                                 loess_trace_hat = loess_trace_hat,
                                 loess_cell = loess_cell,
                                 loess_iterations = loess_iterations,
                                 loess_iterTrace = loess_iterTrace, ...)
        } else if (curve_type == "rcs"){
          boot_obs <- boot::boot(data_raw, calc_obs_blr_rcs_boot, R = CI_R_boot,
                                 data_ms = data_ms,
                                 data_to_plot = data_to_plot,
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
                                 rcs_nk = rcs_nk, ...)
        }

        ### Extract confidence bands
        lower <- apply(boot_obs$t, 2, stats::quantile, probs = alpha, na.rm = TRUE)
        upper <- apply(boot_obs$t, 2, stats::quantile, probs = 1-alpha, na.rm = TRUE)

        ### Produce a warning if any NA values
        if(sum(is.na(boot_obs$t)) > 0){
          warning(paste("WARNING, SOME BOOTSTRAPPED OBSERVED EVENT PROBABILITIES WERE NA FOR STATE", state_k, "\n",
                        "THERE ARE ", sum(apply(boot_obs$t, 1, function(x) {sum(is.na(x)) > 0})), " ITERATIONS WITH NA's \n",
                        "THE MEAN NUMBER OF NA's IN EACH ITERATION IS", mean(apply(boot_obs$t, 1, function(x) {sum(is.na(x))}))
          ))
        }

        ### Assign output
        if ("id" %in% colnames(data_to_plot)){
          output_object_plots[[k]] <- data.frame(
            "id" = data_to_plot[, "id"],
            "pred" = data_to_plot[, paste("tp_pred", state_k, sep = "")],
            "obs" = boot_obs$t0,
            "obs_lower" = lower,
            "obs_upper" = upper)
        } else {
          output_object_plots[[k]] <- data.frame(
            "pred" = data_to_plot[, paste("tp_pred", state_k, sep = "")],
            "obs" = boot_obs$t0,
            "obs_lower" = lower,
            "obs_upper" = upper)
        }
      }

      ###
      ### Mean calibration
      ###
      if (assess_mean == TRUE){
        boot_mean <- boot::boot(data_raw, calc_mean_blr_boot, R = CI_R_boot,
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
                                w_max_follow = w_max_follow, ...)

        ### Extract confidence bands
        lower <- stats::quantile(boot_mean$t[,1], probs = alpha, na.rm = TRUE)
        upper <- stats::quantile(boot_mean$t[,1], probs = 1 - alpha, na.rm = TRUE)

        ### Put into output object
        output_object_mean[[k]] <- c("mean" = boot_mean$t0, "mean_lower" = as.numeric(lower), "mean_upper" = as.numeric(upper))

      }

      # ###
      # ### PLACE HOLDER FOR Weak calibration (calibration slope)
      # ###
      # boot_weak <- boot::boot(data_raw, calc_weak_blr_boot, R = CI_R_boot,
      #                         data_ms = data_ms,
      #                         state_k = state_k,
      #                         j = j,
      #                         s2 = s,
      #                         t = t,
      #                         weights_provided = weights_provided,
      #                         w_function = w_function,
      #                         w_covs = w_covs,
      #                         w_landmark_type = w_landmark_type,
      #                         w_max = w_max,
      #                         w_stabilised = w_stabilised,
      #                         w_max_follow = w_max_follow, ...)
      #
      # ### Extract confidence bands
      # lower <- stats::quantile(boot_weak$t[,1], probs = alpha, na.rm = TRUE)
      # upper <- stats::quantile(boot_weak$t[,1], probs = 1 - alpha, na.rm = TRUE)
      #
      # ### Put into output object
      # output_object_weak[[k]] <- c("slope" = boot_weak$t0, "slope_lower" = as.numeric(lower), "slope_upper" = as.numeric(upper))

    } else if (CI != FALSE & CI_type == "parametric"){

      ### PLACE HOLDER FOR FUTURE INCLUSION OF PARAMETRIC CONFIDENCE INTERVALS

    } else if (CI == FALSE){

      ### Apply above functions to data_raw (note the chosen set of indices samples every patient once)

      ###
      ### Calibration plots
      ###
      if (assess_moderate == TRUE){
        if (curve_type == "loess"){
          pred_obs <- calc_obs_blr_loess_boot(data_raw = data_raw,
                                              indices = 1:nrow(data_raw),
                                              data_ms = data_ms,
                                              data_to_plot = data_to_plot,
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
                                              loess_span = loess_span,
                                              loess_degree = loess_degree,
                                              loess_surface = loess_surface,
                                              loess_trace_hat = loess_trace_hat,
                                              loess_cell = loess_cell,
                                              loess_iterations = loess_iterations,
                                              loess_iterTrace = loess_iterTrace, ...)
        } else if (curve_type == "rcs"){
          pred_obs <- calc_obs_blr_rcs_boot(data_raw = data_raw,
                                            indices = 1:nrow(data_raw),
                                            data_ms = data_ms,
                                            data_to_plot = data_to_plot,
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
                                            rcs_nk = rcs_nk, ...)
        }

        ### Assign output
        if ("id" %in% colnames(data_to_plot)){
          output_object_plots[[k]] <- data.frame("id" = data_to_plot[, "id"],
                                                 "pred" = data_to_plot[, paste("tp_pred", state_k, sep = "")],
                                                 "obs" = pred_obs)
        } else {
          output_object_plots[[k]] <- data.frame(
            "pred" = data_to_plot[, paste("tp_pred", state_k, sep = "")],
            "obs" = pred_obs)
        }
      }


      ###
      ### Mean calibration
      ###
      if (assess_mean == TRUE){
        mean <- calc_mean_blr_boot(data_raw = data_raw,
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
                                   w_max_follow = w_max_follow, ...)

        ### Put into output object
        output_object_mean[[k]] <- c("mean" = mean)

      }

      # ###
      # ### PLACEHOLDER FOR Weak calibration (calibration slope)
      # ###
      # weak <- calc_weak_blr_boot(data_raw = data_raw,
      #                            indices = 1:nrow(data_raw),
      #                            data_ms = data_ms,
      #                            state_k = state_k,
      #                            j = j,
      #                            s2 = s,
      #                            t = t,
      #                            weights_provided = weights_provided,
      #                            w_function = w_function,
      #                            w_covs = w_covs,
      #                            w_landmark_type = w_landmark_type,
      #                            w_max = w_max,
      #                            w_stabilised = w_stabilised,
      #                            w_max_follow = w_max_follow, ...)
      #
      # ### Put into output object
      # output_object_weak[[k]] <- c("weak" = weak)

    }

  }

  ###
  ### If assess_mean == TRUE, and CI == FALSE, collapse the mean calibrations into a vector (to match output from calib_aj and calib_mlr)
  ###
  if (assess_mean == TRUE){
    if (CI == FALSE){

      ### Put into output object
      output_object_mean <- do.call("c", output_object_mean)
      names(output_object_mean) <- paste("state", transitions_out, sep = "")

    }
  }


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

#' Estimate observed event probabilities for state_k using BLR-IPCW and loess smoothers.
#' @description
#' Estimate observed event probabilities for a specific state using BLR-IPCW when `curve_type = 'loess'` specified.
#' Function is called by calibmsm::calib_blr_ipcw, which is called by calibmsm::calib_msm.
#' This function does the heavy lifting, and is what estimates the observed event probabilities.
#'
#' @details
#' Function written in a format so that it can be used in combination with \code{\link[boot]{boot}}
#' for bootstrapping. Specifying `indices = 1:nrow(data_raw)` will return the calibration
#' of interest.
#'
#' @returns A vector of observed event probabilities for state_k. Observed event probabilities
#' are estimated for data points in data_to_plot.
#'
#' @noRd
calc_obs_blr_loess_boot <- function(data_raw,
                                    indices,
                                    data_ms,
                                    data_to_plot,
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
                                    loess_span,
                                    loess_degree,
                                    loess_surface,
                                    loess_trace_hat,
                                    loess_cell,
                                    loess_iterations,
                                    loess_iterTrace, ...){

  # Create bootstrapped dataset
  data_boot <- data_raw[indices, ]

  ## Create landmarked dataset
  data_boot_lmk_js_uncens <-  apply_landmark(data_raw = data_boot, data_ms = data_ms, j = j, s = s2, t = t, exclude_cens_t = TRUE)

  ## Calculate weights
  ## Note this is done in the entire dataset data_raw, which has its own functionality (w_landmark_type) to landmark on j and s, or just s, before
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

  ## Define equation
  eq_loess <- stats::formula(paste("state", state_k, "_bin ~ tp_pred", state_k, sep = ""))

  ## Fit model
  loess_model <- stats::loess(eq_loess,
                              data = data_boot_lmk_js_uncens,
                              weights = data_boot_lmk_js_uncens[, "ipcw"],
                              span = loess_span,
                              degree = loess_degree,
                              control = stats::loess.control(surface = loess_surface,
                                                             statistics = "none",
                                                             trace.hat = loess_trace_hat,
                                                             cell = loess_cell,
                                                             iterations = loess_iterations,
                                                             iterTrace = loess_iterTrace))

  ## Create predicted observed probabilities.
  loess_pred_obs <- predict(loess_model, newdata = data_to_plot)

  return(loess_pred_obs)

}

#' Estimate observed event probabilities for state_k using BLR-IPCW and restricted cubic splines.
#' @description
#' Estimate observed event probabilities for a specific state using BLR-IPCW when `curve_type = 'rcs'` specified.
#' Function is called by calibmsm::calib_blr_ipcw, which is called by calibmsm::calib_msm.
#' This function does the heavy lifting, and is what estimates the observed event probabilities.
#'
#' @details
#' Function written in a format so that it can be used in combination with \code{\link[boot]{boot}}
#' for bootstrapping. Specifying `indices = 1:nrow(data_raw)` will return the calibration
#' of interest.
#'
#' @returns A vector of observed event probabilities for state_k. Observed event probabilities
#' are estimated for data points in data_to_plot.
#'
#' @noRd
calc_obs_blr_rcs_boot <- function(data_raw,
                                  indices,
                                  data_ms,
                                  data_to_plot,
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
                                  rcs_nk, ...){

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

  ## Create restricted cubic splines for the cloglog of the linear predictor for the state of interest
  rcs_mat <- Hmisc::rcspline.eval(data_boot_lmk_js_uncens[,paste("tp_pred_logit", state_k, sep = "")],nk=rcs_nk,inclx=T)
  colnames(rcs_mat) <- paste("rcs_x", 1:ncol(rcs_mat), sep = "")
  knots <- attr(rcs_mat,"knots")

  ## Add the cubic splines for logit of the predicted probability to  data_boot_lmk_js_uncens
  data_boot_lmk_js_uncens <- data.frame(cbind(data_boot_lmk_js_uncens, rcs_mat))

  ## Define equation
  eq_LHS <- paste("state", state_k, "_bin ~ ", sep = "")
  eq_RHS <- paste("rcs_x", 1:ncol(rcs_mat), sep = "", collapse = "+")
  eq_rcs <- stats::formula(paste(eq_LHS, eq_RHS, sep = ""))

  ## Fit model
  ## NB: Warnings are suppressed beceause rms gives the following warning:
  ## In rms::lrm(eq_rcs, data = data_boot_lmk_js_uncens, weights = data_boot_lmk_js_uncens[,:
  ## currently weights are ignored in model validation and bootstrapping lrm fits
  ## We are not using the model validation or bootstrapping aspects of rms::lrm (we are applying bootstrapping ourselves),
  ## meaning the warning is not necessary.
  suppressWarnings(
    rcs_model <- rms::lrm(eq_rcs,
                          data = data_boot_lmk_js_uncens,
                          weights = data_boot_lmk_js_uncens[, "ipcw"])
  )

  ## Create predicted observed probabilities for data_to_plot.

  ## For this, need to calculate the same knot locations calculated for the fitted model.
  rcs_mat_data_to_plot <- Hmisc::rcspline.eval(data_to_plot[,paste("tp_pred_logit", state_k, sep = "")],knots = knots,inclx=T)
  colnames(rcs_mat_data_to_plot) <- paste("rcs_x", 1:ncol(rcs_mat), sep = "")
  #attr(rcs_mat_data_to_plot,"knots")

  ## Add the cubic splines for logit of the predicted probability to data_to_plot
  data_to_plot <- data.frame(cbind(data_to_plot, rcs_mat_data_to_plot))

  ## Calculate observed event probabilities
  rcs_pred_obs <- predict(rcs_model, newdata = data_to_plot, type = "fitted")

  return(rcs_pred_obs)

}


#' Estimate mean calibration using BLR-IPCW for state k.
#' @description
#' Estimate mean calibration using BLR-IPCW for state k.
#'
#' @details
#' Function written in a format so that it can be used in combination with \code{\link[boot]{boot}}
#' for bootstrapping. Specifying `indices = 1:nrow(data_raw)` will return the calibration
#' of interest.
#'
#' @returns Mean calibration and calibration slope for a given state.
#'
#' @noRd
calc_mean_blr_boot <- function(data_raw,
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
                               w_max_follow, ...){

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

  ###
  ### Estimate difference between predicted-observed (calculated using intercept model with offset) and predicted risks

  ## Define equation
  eq_int <- stats::formula(paste("state", state_k, "_bin ~ offset(tp_pred_logit", state_k, ")", sep = ""))

  ## Fit recalibration models to the uncensored observations at time t to calculate intercept
  lrm_int <- stats::glm(eq_int, data = data_boot_lmk_js_uncens, weights = ipcw, family = stats::quasibinomial(link = "logit"))

  ## Generate predicted-observed risks
  pred_obs <- predict(lrm_int, newdata = data_boot_lmk_js_uncens, type = "response")

  ## Extract difference
  int_diff <- mean(pred_obs - data_boot_lmk_js_uncens[, paste("tp_pred", state_k, sep = "")])

  return(int_diff)

}


#' Estimate calibration slope using BLR-IPCW for state k.
#' @description
#' Estimate calibration slope using BLR-IPCW for state k.
#'
#' @details
#' Function written in a format so that it can be used in combination with \code{\link[boot]{boot}}
#' for bootstrapping. Specifying `indices = 1:nrow(data_raw)` will return the calibration
#' of interest.
#'
#' @returns Mean calibration and calibration slope for a given state.
#'
#' @noRd
calc_weak_blr_boot <- function(data_raw,
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
                               w_max_follow, ...){

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

  ###
  ### Estiamte slope

  ## Define equation
  eq_slope <- stats::formula(paste("state", state_k, "_bin ~ tp_pred_logit", state_k, sep = ""))

  ## Fit recalibration models to the uncensored observations at time t to calculate slope
  lrm_slope <- stats::glm(eq_slope, data = data_boot_lmk_js_uncens, weights = ipcw, family = stats::quasibinomial(link = "logit"))

  ## Extract slope and confidence interval
  slope <- c(as.numeric(stats::coefficients(lrm_slope)[2]))

  return(slope)

}




