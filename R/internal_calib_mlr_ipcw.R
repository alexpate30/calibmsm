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
                           data.mstate,
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
                           valid.transitions,
                           assess.moderate,
                           assess.mean, ...){

  ## Create landmarked dataset
  data.raw.lmk.js.uncens <-  apply_landmark(data.raw = data.raw, data.mstate = data.mstate, j = j, s = s, t = t, exclude.cens.t = TRUE)

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

    ### Define equation
    eq.LHS <- paste("state.poly.fac ~ ")
    eq.RHS <- paste("mlr.lp", 1:(length(valid.transitions) - 1), sep = "", collapse = "+")
    eq.mlr <- stats::as.formula(paste(eq.LHS, eq.RHS, sep =""))

    ### Fit model to estimate intercepts
    mlr.int.model <- VGAM::vgam(data.raw.lmk.js.uncens$state.poly.fac ~ 1,
                                offset = as.matrix(data.raw.lmk.js.uncens[, paste("mlr.lp", 1:(length(valid.transitions) - 1), sep = "")]),
                                weights = data.raw.lmk.js.uncens$ipcw,
                                family = VGAM::multinomial(refLevel = ref.cat))

    ###  Now need to calculate predicted-observed values for each state using this model, and calculate mean difference from predicted risks

    ### Create a new dataset with linear predictors from the model (intercept plus offset)
    data.pred.obs.lp <- matrix( nrow = nrow(data.raw.lmk.js.uncens), ncol = length(valid.transitions) - 1)
    for (state in 1:(length(valid.transitions) - 1)){
      data.pred.obs.lp[,state] <-
        exp(stats::coefficients(mlr.int.model)[state] + data.raw.lmk.js.uncens[,paste("mlr.lp", state, sep = "")])
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
    output.object.mean <- colMeans(int.pred.obs - data.raw.lmk.js.uncens[,paste("tp.pred", valid.transitions, sep = "")])
    names(output.object.mean) <- paste("state", valid.transitions, sep = "")

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

