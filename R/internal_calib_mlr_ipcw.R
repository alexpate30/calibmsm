### Internal functions to allow estimation of calibration curves using the MLR-IPCW method.

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
                           valid.transitions, ...){

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

  ### Assign reference category
  ref.cat <- paste(valid.transitions[1])

  ### Apply nominal recalibration framework of van Hoorde et al., (2014)
  calib.model <- VGAM::vgam(eq.mlr, weights = data.raw.lmk.js.uncens[, "ipcw"],
                            data = data.raw.lmk.js.uncens, family = VGAM::multinomial(refLevel = ref.cat))

  ###
  ### Generate predicted-observed risks and add to data.raw.lmk.js.uncens
  ###

  ### For all other functions, I just generate predicted risks for all individuals, and then just plot for those who were uncensored
  ### However, some of the censored individuals are causing an error, therefore I must generate NA vectors, then assign the
  ### predicted observed probabilities to the correct individuals

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
  output.object2 <- vector("list", length(valid.transitions))
  names(output.object2) <- paste("state", valid.transitions, sep = "")

  ## Loop through and create output for each valid transition
  for (k in 1:length(valid.transitions)){

    ## Assign state of interest
    state.k <- valid.transitions[k]

    ## Create output object
    output.object2[[k]] <- data.frame("id" = output.object[, "id"],
                                      "pred" = output.object[, paste("tp.pred", valid.transitions[k], sep = "")],
                                      "obs" = output.object[, paste("mlr.pred.obs", valid.transitions[k], sep = "")])
  }

  return(output.object2)

}

