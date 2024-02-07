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
calib_blr_ipcw <- function(data.raw,
                           data.mstate,
                           tp.pred.plot,
                           j,
                           s,
                           t,
                           curve.type,
                           rcs.nk,
                           loess.span,
                           loess.degree,
                           weights.provided,
                           w.function,
                           w.covs,
                           w.landmark.type,
                           w.max,
                           w.stabilised,
                           w.max.follow,
                           CI,
                           CI.type,
                           CI.R.boot,
                           CI.seed,
                           transitions.out,
                           assess.moderate,
                           assess.mean, ...){

  ###
  ### Create the object of predicted risks over which the calibration plots will be plotted

  ### For calib_blr_ipcw, this is the landmarked cohort of individuals who are also uncensored at time t,
  ### or tp.pred.plot if specified
  if (is.null(tp.pred.plot)){
    ## Note that tp.pred has been added to data.raw and the predicted transition probabilities (and relevant transformations)
    ## are contained within this dataset.
    data.to.plot <- apply_landmark(data.raw = data.raw, data.mstate = data.mstate, j = j, s = s, t = t, exclude.cens.t = TRUE)
  } else if (!is.null(tp.pred.plot)){
    data.to.plot <- tp.pred.plot
  }

  ### Create object to store output
  output.object.plots <- vector("list", length(transitions.out))
  names(output.object.plots) <- paste("state", transitions.out, sep = "")

  output.object.mean <- vector("list", length(transitions.out))
  names(output.object.mean) <- paste("state", transitions.out, sep = "")

  # output.object.weak <- vector("list", length(transitions.out))
  # names(output.object.weak) <- paste("state", transitions.out, sep = "")

  ### Loop through and fit models if BLR selected
  for (k in 1:length(transitions.out)){

    ### Assign state of interest
    state.k <- transitions.out[k]

    ### Calculate confidence intervals
    if (CI != FALSE & CI.type == "bootstrap"){

      ### Define alpha for CI's
      alpha <- (1-CI/100)/2

      ### Set seed for bootstrapping
      set.seed(CI.seed)

      ###
      ### Calibration plots
      ###
      if (assess.moderate == TRUE){
        if (curve.type == "loess"){
          boot.obs <- boot::boot(data.raw, calc_obs_blr_loess_boot, R = CI.R.boot,
                                 data.mstate = data.mstate,
                                 data.to.plot = data.to.plot,
                                 state.k = state.k,
                                 j = j,
                                 s2 = s,
                                 t = t,
                                 weights.provided = weights.provided,
                                 w.function = w.function,
                                 w.covs = w.covs,
                                 w.landmark.type = w.landmark.type,
                                 w.max = w.max,
                                 w.stabilised = w.stabilised,
                                 w.max.follow = w.max.follow,
                                 loess.span = loess.span,
                                 loess.degree = loess.degree, ...)
        } else if (curve.type == "rcs"){
          boot.obs <- boot::boot(data.raw, calc_obs_blr_rcs_boot, R = CI.R.boot,
                                 data.mstate = data.mstate,
                                 data.to.plot = data.to.plot,
                                 state.k = state.k,
                                 j = j,
                                 s2 = s,
                                 t = t,
                                 weights.provided = weights.provided,
                                 w.function = w.function,
                                 w.covs = w.covs,
                                 w.landmark.type = w.landmark.type,
                                 w.max = w.max,
                                 w.stabilised = w.stabilised,
                                 w.max.follow = w.max.follow,
                                 rcs.nk = rcs.nk, ...)
        }

        ### Extract confidence bands
        lower <- apply(boot.obs$t, 2, stats::quantile, probs = alpha, na.rm = TRUE)
        upper <- apply(boot.obs$t, 2, stats::quantile, probs = 1-alpha, na.rm = TRUE)

        ### Produce a warning if any NA values
        if(sum(is.na(boot.obs$t)) > 0){
          warning(paste("WARNING, SOME BOOTSTRAPPED OBSERVED EVENT PROBABILITIES WERE NA FOR STATE", state.k, "\n",
                        "THERE ARE ", sum(apply(boot.obs$t, 1, function(x) {sum(is.na(x)) > 0})), " ITERATIONS WITH NA's \n",
                        "THE MEAN NUMBER OF NA's IN EACH ITERATION IS", mean(apply(boot.obs$t, 1, function(x) {sum(is.na(x))}))
          ))
        }

        ### Assign output
        if ("id" %in% colnames(data.to.plot)){
          output.object.plots[[k]] <- data.frame(
            "id" = data.to.plot[, "id"],
            "pred" = data.to.plot[, paste("tp.pred", state.k, sep = "")],
            "obs" = boot.obs$t0,
            "obs.lower" = lower,
            "obs.upper" = upper)
        } else {
          output.object.plots[[k]] <- data.frame(
            "pred" = data.to.plot[, paste("tp.pred", state.k, sep = "")],
            "obs" = boot.obs$t0,
            "obs.lower" = lower,
            "obs.upper" = upper)
        }
      }

      ###
      ### Mean calibration
      ###
      if (assess.mean == TRUE){
        boot.mean <- boot::boot(data.raw, calc_mean_blr_boot, R = CI.R.boot,
                                data.mstate = data.mstate,
                                state.k = state.k,
                                j = j,
                                s2 = s,
                                t = t,
                                weights.provided = weights.provided,
                                w.function = w.function,
                                w.covs = w.covs,
                                w.landmark.type = w.landmark.type,
                                w.max = w.max,
                                w.stabilised = w.stabilised,
                                w.max.follow = w.max.follow, ...)

        ### Extract confidence bands
        lower <- stats::quantile(boot.mean$t[,1], probs = alpha, na.rm = TRUE)
        upper <- stats::quantile(boot.mean$t[,1], probs = 1 - alpha, na.rm = TRUE)

        ### Put into output object
        output.object.mean[[k]] <- c("mean" = boot.mean$t0, "mean.lower" = as.numeric(lower), "mean.upper" = as.numeric(upper))

      }

      # ###
      # ### PLACE HOLDER FOR Weak calibration (calibration slope)
      # ###
      # boot.weak <- boot::boot(data.raw, calc_weak_blr_boot, R = CI.R.boot,
      #                         data.mstate = data.mstate,
      #                         state.k = state.k,
      #                         j = j,
      #                         s2 = s,
      #                         t = t,
      #                         weights.provided = weights.provided,
      #                         w.function = w.function,
      #                         w.covs = w.covs,
      #                         w.landmark.type = w.landmark.type,
      #                         w.max = w.max,
      #                         w.stabilised = w.stabilised,
      #                         w.max.follow = w.max.follow, ...)
      #
      # ### Extract confidence bands
      # lower <- stats::quantile(boot.weak$t[,1], probs = alpha, na.rm = TRUE)
      # upper <- stats::quantile(boot.weak$t[,1], probs = 1 - alpha, na.rm = TRUE)
      #
      # ### Put into output object
      # output.object.weak[[k]] <- c("slope" = boot.weak$t0, "slope.lower" = as.numeric(lower), "slope.upper" = as.numeric(upper))

    } else if (CI != FALSE & CI.type == "parametric"){

      ### PLACE HOLDER FOR FUTURE INCLUSION OF PARAMETRIC CONFIDENCE INTERVALS

    } else if (CI == FALSE){

      ### Apply above functions to data.raw (note the chosen set of indices samples every patient once)

      ###
      ### Calibration plots
      ###
      if (assess.moderate == TRUE){
        if (curve.type == "loess"){
          pred.obs <- calc_obs_blr_loess_boot(data.raw = data.raw,
                                              indices = 1:nrow(data.raw),
                                              data.mstate = data.mstate,
                                              data.to.plot = data.to.plot,
                                              state.k = state.k,
                                              j = j,
                                              s2 = s,
                                              t = t,
                                              weights.provided = weights.provided,
                                              w.function = w.function,
                                              w.covs = w.covs,
                                              w.landmark.type = w.landmark.type,
                                              w.max = w.max,
                                              w.stabilised = w.stabilised,
                                              w.max.follow = w.max.follow,
                                              loess.span = loess.span,
                                              loess.degree = loess.degree, ...)
        } else if (curve.type == "rcs"){
          pred.obs <- calc_obs_blr_rcs_boot(data.raw = data.raw,
                                            indices = 1:nrow(data.raw),
                                            data.mstate = data.mstate,
                                            data.to.plot = data.to.plot,
                                            state.k = state.k,
                                            j = j,
                                            s2 = s,
                                            t = t,
                                            weights.provided = weights.provided,
                                            w.function = w.function,
                                            w.covs = w.covs,
                                            w.landmark.type = w.landmark.type,
                                            w.max = w.max,
                                            w.stabilised = w.stabilised,
                                            w.max.follow = w.max.follow,
                                            rcs.nk = rcs.nk, ...)
        }

        ### Assign output
        if ("id" %in% colnames(data.to.plot)){
          output.object.plots[[k]] <- data.frame("id" = data.to.plot[, "id"],
                                                 "pred" = data.to.plot[, paste("tp.pred", state.k, sep = "")],
                                                 "obs" = pred.obs)
        } else {
          output.object.plots[[k]] <- data.frame(
            "pred" = data.to.plot[, paste("tp.pred", state.k, sep = "")],
            "obs" = pred.obs)
        }
      }


      ###
      ### Mean calibration
      ###
      if (assess.mean == TRUE){
        mean <- calc_mean_blr_boot(data.raw = data.raw,
                                   indices = 1:nrow(data.raw),
                                   data.mstate = data.mstate,
                                   state.k = state.k,
                                   j = j,
                                   s2 = s,
                                   t = t,
                                   weights.provided = weights.provided,
                                   w.function = w.function,
                                   w.covs = w.covs,
                                   w.landmark.type = w.landmark.type,
                                   w.max = w.max,
                                   w.stabilised = w.stabilised,
                                   w.max.follow = w.max.follow, ...)

        ### Put into output object
        output.object.mean[[k]] <- c("mean" = mean)

      }

      # ###
      # ### PLACEHOLDER FOR Weak calibration (calibration slope)
      # ###
      # weak <- calc_weak_blr_boot(data.raw = data.raw,
      #                            indices = 1:nrow(data.raw),
      #                            data.mstate = data.mstate,
      #                            state.k = state.k,
      #                            j = j,
      #                            s2 = s,
      #                            t = t,
      #                            weights.provided = weights.provided,
      #                            w.function = w.function,
      #                            w.covs = w.covs,
      #                            w.landmark.type = w.landmark.type,
      #                            w.max = w.max,
      #                            w.stabilised = w.stabilised,
      #                            w.max.follow = w.max.follow, ...)
      #
      # ### Put into output object
      # output.object.weak[[k]] <- c("weak" = weak)

    }

  }

  ###
  ### If assess.mean == TRUE, and CI == FALSE, collapse the mean calibrations into a vector (to match output from calib_aj and calib_mlr)
  ###
  if (assess.mean == TRUE){
    if (CI == FALSE){

      ### Put into output object
      output.object.mean <- do.call("c", output.object.mean)
      names(output.object.mean) <- paste("state", transitions.out, sep = "")

    }
  }


  ### Define combined output object
  output.object.comb <- vector("list")

  if (assess.moderate == TRUE){
    output.object.comb[["plotdata"]] <- output.object.plots
  }
  if (assess.mean == TRUE){
    output.object.comb[["mean"]] <- output.object.mean
  }
  # if (assess.weak == TRUE){
  #   output.object.comb[["weak"]] <- output.object.weak
  # }

  return(output.object.comb)

}

#' Estimate observed event probabilities for state.k using BLR-IPCW and loess smoothers.
#' @description
#' Estimate observed event probabilities for a specific state using BLR-IPCW when `curve.type = 'loess'` specified.
#' Function is called by calibmsm::calib_blr_ipcw, which is called by calibmsm::calib_msm.
#' This function does the heavy lifting, and is what estimates the observed event probabilities.
#'
#' @details
#' Function written in a format so that it can be used in combination with \code{\link[boot]{boot}}
#' for bootstrapping. Specifying `indices = 1:nrow(data.raw)` will return the calibration
#' of interest.
#'
#' @returns A vector of observed event probabilities for state.k. Observed event probabilities
#' are estimated for data points in data.to.plot.
#'
#' @noRd
calc_obs_blr_loess_boot <- function(data.raw,
                                    indices,
                                    data.mstate,
                                    data.to.plot,
                                    state.k,
                                    j,
                                    s2, # can't use 's' because it matches an argument for the boot function
                                    t,
                                    weights.provided,
                                    w.function,
                                    w.covs,
                                    w.landmark.type,
                                    w.max,
                                    w.stabilised,
                                    w.max.follow,
                                    loess.span,
                                    loess.degree, ...){

  # Create bootstrapped dataset
  data.boot <- data.raw[indices, ]

  ## Create landmarked dataset
  data.boot.lmk.js.uncens <-  apply_landmark(data.raw = data.boot, data.mstate = data.mstate, j = j, s = s2, t = t, exclude.cens.t = TRUE)

  ## Calculate weights
  ## Note this is done in the entire dataset data.raw, which has its own functionality (w.landmark.type) to landmark on j and s, or just s, before
  ## calculating the weights
  if (weights.provided == FALSE){

    ## If custom function for estimating weights has been inputted ("w.function"), replace "calc_weights" with this function
    if (!is.null(w.function)){
      calc_weights <- w.function
    }

    ## Calculate the weights
    weights <- calc_weights(data.mstate = data.mstate,
                            data.raw = data.boot,
                            covs = w.covs,
                            t = t,
                            s = s2,
                            landmark.type = w.landmark.type,
                            j = j,
                            max.weight = w.max,
                            stabilised = w.stabilised,
                            max.follow = w.max.follow,
                            ...)

    ## Add weights to data.boot
    data.boot.lmk.js.uncens <- dplyr::left_join(data.boot.lmk.js.uncens, dplyr::distinct(weights), by = dplyr::join_by(id))

  }

  ## Define equation
  eq.loess <- stats::formula(paste("state", state.k, ".bin ~ tp.pred", state.k, sep = ""))

  ## Fit model
  loess.model <- stats::loess(eq.loess,
                              data = data.boot.lmk.js.uncens,
                              weights = data.boot.lmk.js.uncens[, "ipcw"],
                              span = loess.span,
                              degree = loess.degree)

  ## Create predicted observed probabilities.
  loess.pred.obs <- predict(loess.model, newdata = data.to.plot)

  return(loess.pred.obs)

}

#' Estimate observed event probabilities for state.k using BLR-IPCW and restricted cubic splines.
#' @description
#' Estimate observed event probabilities for a specific state using BLR-IPCW when `curve.type = 'rcs'` specified.
#' Function is called by calibmsm::calib_blr_ipcw, which is called by calibmsm::calib_msm.
#' This function does the heavy lifting, and is what estimates the observed event probabilities.
#'
#' @details
#' Function written in a format so that it can be used in combination with \code{\link[boot]{boot}}
#' for bootstrapping. Specifying `indices = 1:nrow(data.raw)` will return the calibration
#' of interest.
#'
#' @returns A vector of observed event probabilities for state.k. Observed event probabilities
#' are estimated for data points in data.to.plot.
#'
#' @noRd
calc_obs_blr_rcs_boot <- function(data.raw,
                                  indices,
                                  data.mstate,
                                  data.to.plot,
                                  state.k,
                                  j,
                                  s2, # can't use 's' because it matches an argument for the boot function
                                  t,
                                  weights.provided,
                                  w.function,
                                  w.covs,
                                  w.landmark.type,
                                  w.max,
                                  w.stabilised,
                                  w.max.follow,
                                  rcs.nk, ...){

  ## Create bootstrapped dataset
  data.boot <- data.raw[indices, ]

  ## Create landmarked dataset
  data.boot.lmk.js.uncens <-  apply_landmark(data.raw = data.boot, data.mstate = data.mstate, j = j, s = s2, t = t, exclude.cens.t = TRUE)

  ## Calculate weights
  ## Note this is done in the entire dataset data.boot, which has its own functionality (w.landmark.type) to landmark on j and s, or just s, before
  ## calculating the weights
  if (weights.provided == FALSE){

    ## If custom function for estimating weights has been inputted ("w.function"), replace "calc_weights" with this function
    if (!is.null(w.function)){
      calc_weights <- w.function
    }

    ## Calculate the weights
    weights <- calc_weights(data.mstate = data.mstate,
                            data.raw = data.boot,
                            covs = w.covs,
                            t = t,
                            s = s2,
                            landmark.type = w.landmark.type,
                            j = j,
                            max.weight = w.max,
                            stabilised = w.stabilised,
                            max.follow = w.max.follow,
                            ...)

    ## Add weights to data.boot
    data.boot.lmk.js.uncens <- dplyr::left_join(data.boot.lmk.js.uncens, dplyr::distinct(weights), by = dplyr::join_by(id))

  }

  ## Create restricted cubic splines for the cloglog of the linear predictor for the state of interest
  rcs.mat <- Hmisc::rcspline.eval(data.boot.lmk.js.uncens[,paste("tp.pred.logit", state.k, sep = "")],nk=rcs.nk,inclx=T)
  colnames(rcs.mat) <- paste("rcs.x", 1:ncol(rcs.mat), sep = "")
  knots <- attr(rcs.mat,"knots")

  ## Add the cubic splines for logit of the predicted probability to  data.boot.lmk.js.uncens
  data.boot.lmk.js.uncens <- data.frame(cbind(data.boot.lmk.js.uncens, rcs.mat))

  ## Define equation
  eq.LHS <- paste("state", state.k, ".bin ~ ", sep = "")
  eq.RHS <- paste("rcs.x", 1:ncol(rcs.mat), sep = "", collapse = "+")
  eq.rcs <- stats::formula(paste(eq.LHS, eq.RHS, sep = ""))

  ## Fit model
  ## NB: Warnings are suppressed beceause rms gives the following warning:
  ## In rms::lrm(eq.rcs, data = data.boot.lmk.js.uncens, weights = data.boot.lmk.js.uncens[,:
  ## currently weights are ignored in model validation and bootstrapping lrm fits
  ## We are not using the model validation or bootstrapping aspects of rms::lrm (we are applying bootstrapping ourselves),
  ## meaning the warning is not necessary.
  suppressWarnings(
    rcs.model <- rms::lrm(eq.rcs,
                          data = data.boot.lmk.js.uncens,
                          weights = data.boot.lmk.js.uncens[, "ipcw"])
  )

  ## Create predicted observed probabilities for data.to.plot.

  ## For this, need to calculate the same knot locations calculated for the fitted model.
  rcs.mat.data.to.plot <- Hmisc::rcspline.eval(data.to.plot[,paste("tp.pred.logit", state.k, sep = "")],knots = knots,inclx=T)
  colnames(rcs.mat.data.to.plot) <- paste("rcs.x", 1:ncol(rcs.mat), sep = "")
  #attr(rcs.mat.data.to.plot,"knots")

  ## Add the cubic splines for logit of the predicted probability to data.to.plot
  data.to.plot <- data.frame(cbind(data.to.plot, rcs.mat.data.to.plot))

  ## Calculate observed event probabilities
  rcs.pred.obs <- predict(rcs.model, newdata = data.to.plot, type = "fitted")

  return(rcs.pred.obs)

}


#' Estimate mean calibration using BLR-IPCW for state k.
#' @description
#' Estimate mean calibration using BLR-IPCW for state k.
#'
#' @details
#' Function written in a format so that it can be used in combination with \code{\link[boot]{boot}}
#' for bootstrapping. Specifying `indices = 1:nrow(data.raw)` will return the calibration
#' of interest.
#'
#' @returns Mean calibration and calibration slope for a given state.
#'
#' @noRd
calc_mean_blr_boot <- function(data.raw,
                               indices,
                               data.mstate,
                               state.k,
                               j,
                               s2, # can't use 's' because it matches an argument for the boot function
                               t,
                               weights.provided,
                               w.function,
                               w.covs,
                               w.landmark.type,
                               w.max,
                               w.stabilised,
                               w.max.follow, ...){

  ## Create bootstrapped dataset
  data.boot <- data.raw[indices, ]

  ## Create landmarked dataset
  data.boot.lmk.js.uncens <-  apply_landmark(data.raw = data.boot, data.mstate = data.mstate, j = j, s = s2, t = t, exclude.cens.t = TRUE)

  ## Calculate weights
  ## Note this is done in the entire dataset data.boot, which has its own functionality (w.landmark.type) to landmark on j and s, or just s, before
  ## calculating the weights
  if (weights.provided == FALSE){

    ## If custom function for estimating weights has been inputted ("w.function"), replace "calc_weights" with this function
    if (!is.null(w.function)){
      calc_weights <- w.function
    }

    ## Calculate the weights
    weights <- calc_weights(data.mstate = data.mstate,
                            data.raw = data.boot,
                            covs = w.covs,
                            t = t,
                            s = s2,
                            landmark.type = w.landmark.type,
                            j = j,
                            max.weight = w.max,
                            stabilised = w.stabilised,
                            max.follow = w.max.follow,
                            ...)

    ## Add weights to data.boot
    data.boot.lmk.js.uncens <- dplyr::left_join(data.boot.lmk.js.uncens, dplyr::distinct(weights), by = dplyr::join_by(id))

  }

  ###
  ### Estimate difference between predicted-observed (calculated using intercept model with offset) and predicted risks

  ## Define equation
  eq.int <- stats::formula(paste("state", state.k, ".bin ~ offset(tp.pred.logit", state.k, ")", sep = ""))

  ## Fit recalibration models to the uncensored observations at time t to calculate intercept
  lrm.int <- stats::glm(eq.int, data = data.boot.lmk.js.uncens, weights = ipcw, family = stats::quasibinomial(link = "logit"))

  ## Generate predicted-observed risks
  pred.obs <- predict(lrm.int, newdata = data.boot.lmk.js.uncens, type = "response")

  ## Extract difference
  int.diff <- mean(pred.obs - data.boot.lmk.js.uncens[, paste("tp.pred", state.k, sep = "")])

  return(int.diff)

}


#' Estimate calibration slope using BLR-IPCW for state k.
#' @description
#' Estimate calibration slope using BLR-IPCW for state k.
#'
#' @details
#' Function written in a format so that it can be used in combination with \code{\link[boot]{boot}}
#' for bootstrapping. Specifying `indices = 1:nrow(data.raw)` will return the calibration
#' of interest.
#'
#' @returns Mean calibration and calibration slope for a given state.
#'
#' @noRd
calc_weak_blr_boot <- function(data.raw,
                               indices,
                               data.mstate,
                               state.k,
                               j,
                               s2, # can't use 's' because it matches an argument for the boot function
                               t,
                               weights.provided,
                               w.function,
                               w.covs,
                               w.landmark.type,
                               w.max,
                               w.stabilised,
                               w.max.follow, ...){

  ## Create bootstrapped dataset
  data.boot <- data.raw[indices, ]

  ## Create landmarked dataset
  data.boot.lmk.js.uncens <-  apply_landmark(data.raw = data.boot, data.mstate = data.mstate, j = j, s = s2, t = t, exclude.cens.t = TRUE)

  ## Calculate weights
  ## Note this is done in the entire dataset data.boot, which has its own functionality (w.landmark.type) to landmark on j and s, or just s, before
  ## calculating the weights
  if (weights.provided == FALSE){

    ## If custom function for estimating weights has been inputted ("w.function"), replace "calc_weights" with this function
    if (!is.null(w.function)){
      calc_weights <- w.function
    }

    ## Calculate the weights
    weights <- calc_weights(data.mstate = data.mstate,
                            data.raw = data.boot,
                            covs = w.covs,
                            t = t,
                            s = s2,
                            landmark.type = w.landmark.type,
                            j = j,
                            max.weight = w.max,
                            stabilised = w.stabilised,
                            max.follow = w.max.follow,
                            ...)

    ## Add weights to data.boot
    data.boot.lmk.js.uncens <- dplyr::left_join(data.boot.lmk.js.uncens, dplyr::distinct(weights), by = dplyr::join_by(id))

  }

  ###
  ### Estiamte slope

  ## Define equation
  eq.slope <- stats::formula(paste("state", state.k, ".bin ~ tp.pred.logit", state.k, sep = ""))

  ## Fit recalibration models to the uncensored observations at time t to calculate slope
  lrm.slope <- stats::glm(eq.slope, data = data.boot.lmk.js.uncens, weights = ipcw, family = stats::quasibinomial(link = "logit"))

  ## Extract slope and confidence interval
  slope <- c(as.numeric(stats::coefficients(lrm.slope)[2]))

  return(slope)

}




