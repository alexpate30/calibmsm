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
calib_mlr_ipcw <- function(data.raw,
                           data.ms,
                           j,
                           s,
                           t,
                           weights.provided,
                           w.covs,
                           w.function,
                           w.landmark.type,
                           w.max,
                           w.stabilised,
                           w.max.follow,
                           mlr.smoother.type,
                           mlr.ps.int,
                           mlr.degree,
                           mlr.s.df,
                           mlr.niknots,
                           CI,
                           CI.R.boot,
                           CI.seed,
                           valid.transitions,
                           assess.moderate,
                           assess.mean, ...){

  ## Create landmarked dataset
  data.raw.lmk.js.uncens <-  apply_landmark(data.raw = data.raw, data.ms = data.ms, j = j, s = s, t = t, exclude.cens.t = TRUE)

  ## Calculate weights
  ## Note this is done in the entire dataset data.boot, which has its own functionality (w.landmark.type) to landmark on j and s, or just s, before
  ## calculating the weights
  if (weights.provided == FALSE){

    ## If custom function for estimating weights has been inputted ("w.function"), replace "calc_weights" with this function
    if (!is.null(w.function)){
      calc_weights <- w.function
    }

    ## Calculate the weights
    weights <- calc_weights(data.ms = data.ms,
                            data.raw = data.raw,
                            covs = w.covs,
                            t = t,
                            s = s,
                            landmark.type = w.landmark.type,
                            j = j,
                            max.weight = w.max,
                            stabilised = w.stabilised,
                            max.follow = w.max.follow,
                            ...)

    ## Add weights to data.boot
    data.raw.lmk.js.uncens <- dplyr::left_join(data.raw.lmk.js.uncens, dplyr::distinct(weights), by = dplyr::join_by(id))

  }

  ### Assign reference category
  ref.cat <- paste(valid.transitions[1])

  ###
  ### Calibration plots/Moderate Calibration
  ###
  if (assess.moderate == TRUE){

    ### Define equation
    eq.LHS <- paste("state.poly.fac ~ ")
    if (mlr.smoother.type == "s"){
      eq.RHS <- paste("s(mlr.lp", 1:(length(valid.transitions) - 1), ", df = mlr.s.df)", sep = "", collapse = "+")
    } else if (mlr.smoother.type == "sm.ps"){
      eq.RHS <- paste("sm.ps(mlr.lp", 1:(length(valid.transitions) - 1), ", ps.int = mlr.ps.int, degree = mlr.degree)", sep = "", collapse = "+")
    } else if (mlr.smoother.type == "sm.os"){
      eq.RHS <- paste("sm.os(mlr.lp", 1:(length(valid.transitions) - 1), ", niknots = mlr.niknots)", sep = "", collapse = "+")
    }
    eq.mlr <- stats::as.formula(paste(eq.LHS, eq.RHS, sep =""))

    ### Apply nominal recalibration framework of van Hoorde et al., (2014)
    calib.model <- VGAM::vgam(eq.mlr, weights = data.raw.lmk.js.uncens[, "ipcw"],
                              data = data.raw.lmk.js.uncens, family = VGAM::multinomial(refLevel = ref.cat))

    ###
    ### Generate predicted-observed risks and add to data.raw.lmk.js.uncens

    ## Create dataframe to store
    dat.mlr.pred.obs <- data.frame(matrix(NA, ncol = length(valid.transitions), nrow = nrow(data.raw.lmk.js.uncens)))
    ## Assign colnames
    colnames(dat.mlr.pred.obs) <- paste("mlr.pred.obs", valid.transitions, sep = "")
    ## Calc pred.obs for those who are uncensored
    mlr.pred.obs <- VGAM::predictvglm(calib.model, newdata = data.raw.lmk.js.uncens, type = "response")
    ## Assign to appropriate individuals
    dat.mlr.pred.obs[, ] <- mlr.pred.obs

    ### Then add it to data.raw.lmk.js.uncens
    data.raw.lmk.js.uncens <- cbind(data.raw.lmk.js.uncens, dat.mlr.pred.obs)

    ### Assign output
    output.object <- dplyr::select(data.raw.lmk.js.uncens, id, paste("tp.pred", valid.transitions, sep = ""), paste("mlr.pred.obs", valid.transitions, sep = ""))

    ### Get plotdata in same format as calib_blr
    ## Start by creating new output object
    output.object.plots <- vector("list", length(valid.transitions))
    names(output.object.plots) <- paste("state", valid.transitions, sep = "")

    ## Loop through and create output object for each valid transition (same format as for BLR-IPCW and PV)
    for (k in 1:length(valid.transitions)){

      ## Assign state of interest
      state.k <- valid.transitions[k]

      ## Create output object
      output.object.plots[[k]] <- data.frame("id" = output.object[, "id"],
                                             "pred" = output.object[, paste("tp.pred", valid.transitions[k], sep = "")],
                                             "obs" = output.object[, paste("mlr.pred.obs", valid.transitions[k], sep = "")])
    }

  }

  ###
  ### Calibration intercept
  ###
  if (assess.mean == TRUE){

    if (CI != FALSE){

      ### Create object to store plot data
      output.object.mean <- vector("list", length(valid.transitions))
      names(output.object.mean) <- paste("state", valid.transitions, sep = "")

      ### Define alpha for CI's
      alpha <- (1-CI/100)/2

      ### Set seed for bootstrapping
      set.seed(CI.seed)

      ### Apply bootstrapping
      boot.mean <- boot::boot(data.raw, calc_mean_mlr_boot, R = CI.R.boot,
                              data.ms = data.ms,
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
                              valid.transitions = valid.transitions, ...)

      ### Extract confidence bands
      lower <- apply(boot.mean$t, 2, stats::quantile, probs = alpha, na.rm = TRUE)
      upper <- apply(boot.mean$t, 2, stats::quantile, probs = 1-alpha, na.rm = TRUE)

      ### Cycle through states and assign to correct part of output object
      for (state in 1:length(valid.transitions)){

        ### Put into output object
        output.object.mean[[state]] <- c("mean" = as.numeric(boot.mean$t0[state]),
                                         "mean.lower" = as.numeric(lower[state]),
                                         "mean.upper" = as.numeric(upper[state]))

      }

    } else if (CI == FALSE){

      ### Assess mean calibration
      output.object.mean <- calc_mean_mlr_boot(data.raw = data.raw,
                                               indices = 1:nrow(data.raw),
                                               data.ms = data.ms,
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
                                               valid.transitions = valid.transitions, ...)

      ### Put into output object
      names(output.object.mean) <- paste("state", valid.transitions, sep = "")

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
  # clist <- list("(Intercept)" = i, "mlr.lp1" = i1, "mlr.lp2" = i2, "mlr.lp3" = i3, "mlr.lp4" = i4)
  # clist
  #
  # ### Apply nominal recalibration framework
  # ###
  #
  # ### Define equation
  # eq.LHS <- paste("state.poly.fac ~ ")
  # eq.RHS <- paste("mlr.lp", 1:(length(valid.transitions) - 1), sep = "", collapse = "+")
  # eq.mlr <- stats::as.formula(paste(eq.LHS, eq.RHS, sep =""))
  #
  # ### Fit model to estimate slopes
  # mlr.slope.model <- VGAM::vgam(state.poly.fac ~ mlr.lp1 + mlr.lp2 + mlr.lp3 + mlr.lp4, constraints = clist, weights = ipcw,
  #                           data = data.raw.lmk.js.uncens, family = VGAM::multinomial(refLevel = ref.cat))
  #
  # ### Extract slopes and standard errors
  # mlr.slopes <- mlr.slope.model@coefficients[paste("mlr.lp", 1:(length(valid.transitions) - 1), sep = "")]
  # mlr.slopes.se <- sqrt(diag(stats::vcov(mlr.slope.model))[paste("mlr.lp", 1:(length(valid.transitions) - 1), sep = "")])
  #
  # ### Define output object
  # output.object.weak <- list("slopes" = mlr.slopes, "slopes.se" = mlr.slopes.se)

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



#' Estimate mean calibration using MLR-IPCW.
#' @description
#' Estimate mean calibration using MLR-IPCW.
#'
#' @details
#' Function written in a format so that it can be used in combination with \code{\link[boot]{boot}}
#' for bootstrapping. Specifying `indices = 1:nrow(data.raw)` will return the calibration
#' of interest.
#'
#' @returns Mean calibration and calibration slope for a given state.
#'
#' @noRd
calc_mean_mlr_boot <- function(data.raw,
                               indices,
                               data.ms,
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
                               valid.transitions, ...){

  ## Create bootstrapped dataset
  data.boot <- data.raw[indices, ]

  ## Create landmarked dataset
  data.boot.lmk.js.uncens <-  apply_landmark(data.raw = data.boot, data.ms = data.ms, j = j, s = s2, t = t, exclude.cens.t = TRUE)

  ## Calculate weights
  ## Note this is done in the entire dataset data.boot, which has its own functionality (w.landmark.type) to landmark on j and s, or just s, before
  ## calculating the weights
  if (weights.provided == FALSE){

    ## If custom function for estimating weights has been inputted ("w.function"), replace "calc_weights" with this function
    if (!is.null(w.function)){
      calc_weights <- w.function
    }

    ## Calculate the weights
    weights <- calc_weights(data.ms = data.ms,
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

  ### Assign reference category
  ref.cat <- paste(valid.transitions[1])

  ### Define equation
  eq.LHS <- paste("state.poly.fac ~ ")
  eq.RHS <- paste("mlr.lp", 1:(length(valid.transitions) - 1), sep = "", collapse = "+")
  eq.mlr <- stats::as.formula(paste(eq.LHS, eq.RHS, sep =""))

  ### Fit model to estimate intercepts
  mlr.int.model <- VGAM::vgam(data.boot.lmk.js.uncens$state.poly.fac ~ 1,
                              offset = as.matrix(data.boot.lmk.js.uncens[, paste("mlr.lp", 1:(length(valid.transitions) - 1), sep = "")]),
                              weights = data.boot.lmk.js.uncens$ipcw,
                              family = VGAM::multinomial(refLevel = ref.cat))

  ###  Now need to calculate predicted-observed values for each state using this model, and calculate mean difference from predicted risks

  ### Create a new dataset with linear predictors from the model (intercept plus offset)
  data.pred.obs.lp <- matrix( nrow = nrow(data.boot.lmk.js.uncens), ncol = length(valid.transitions) - 1)
  for (state in 1:(length(valid.transitions) - 1)){
    data.pred.obs.lp[,state] <-
      exp(stats::coefficients(mlr.int.model)[state] + data.boot.lmk.js.uncens[,paste("mlr.lp", state, sep = "")])
  }

  ### Now calculate the predicted probabilities
  p1 <- 1/(1 + rowSums(data.pred.obs.lp))
  p.rest <- data.pred.obs.lp/(1 + rowSums(data.pred.obs.lp))

  ### Put together into dataset of predicted-observed values from the intercept model
  int.pred.obs <- data.frame(cbind(p1, p.rest))
  colnames(int.pred.obs) <- paste("pred.obs.", valid.transitions, sep = "")

  ### Check rows sum to 1 as expected
  # rowSums(int.pred.obs)

  ### Get difference between predicted-observed and predicted risks
  int.diff <- colMeans(int.pred.obs - data.boot.lmk.js.uncens[,paste("tp.pred", valid.transitions, sep = "")])
  names(int.diff) <- paste("state", valid.transitions, sep = "")

  return(int.diff)

}
