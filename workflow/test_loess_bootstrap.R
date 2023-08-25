### Test bootstrap with loess
load_all()

calib_blr_test <- function(data.mstate,
                      data.raw,
                      j,
                      s,
                      t,
                      tp.pred,
                      curve.type = "rcs",
                      rcs.nk = 3,
                      loess.span = 0.75,
                      loess.degree = 2,
                      weights = NULL,
                      w.function = NULL,
                      w.covs = NULL,
                      w.landmark.type = "state",
                      w.max = 10,
                      w.stabilised = FALSE,
                      w.max.follow = NULL,
                      CI = FALSE,
                      CI.type = "bootstrap",
                      CI.R.boot = NULL,
                      data.pred.plot = NULL,
                      transitions.out = NULL, ...){

  ###
  ### Warnings and errors
  ###

  ### Stop if patients in data.raw are not in data.mstate
  if (!base::all(unique(data.raw$id) %in% unique(data.mstate$id))){
    stop("All patients in data.raw are not contained in data.mstate. Landmarking cannot be applied.")
  }

  ### Warning if weights inputted manually, and confidence interval requested internally
  if ((CI != FALSE) & !is.null(weights)){
    stop("Estimation of confidence interval using internal bootstrapping procedure was requested. This is not possible with user-inputted weights.")
  }

  ### Stop if CI.type = "bootstrap" but CI.R.boot not specified
  if (CI != FALSE){
    if (!(CI.type %in% c("parametric", "bootstrap"))){
      stop("CI.type takes values in 'parametric' and 'bootstrap'")
    } else {
      if (CI.type == "bootstrap" & is.null(CI.R.boot)){
        stop("Must specify number of bootstrap replicates for confidence interval using CI.R.boot.")
      }
    }
  }

  ### Check if transitions.out is only specified for non-zero columns
  if (!is.null(transitions.out)){
    if (sum(c(colSums(tp.pred) == 0)[transitions.out] == TRUE) > 0){
      stop("Calibraiton curves have been requested for transitions into states which have zero probability of occuring.")
    }
  }

  ### If vector of weights and custom function for specifying weights both inputted, give error
  if (!is.null(weights) & !is.null(w.function)){
    stop("Cannot specify weights manually, and specify a custom function for estimating the weights. Choose one or the other.")
  }

  ### If a vector of weights has been provided, add it to the dataset
  if (!is.null(weights)){
    ### First check whether it is the correct length (NA's should be present)
    if (length(weights) != nrow(data.raw)){
      stop("Weights vector not same length as data.raw")
    } else {
      data.raw$ipcw <- weights
    }
  }

  ### If custom function for estimating weights has been inputted ("w.function"), replace "calc_weights" with this function
  if (!is.null(w.function)){
    ### stop if w.function doesn't have correct arguments
    if(!all(names(formals(calc_weights)) %in% names(formals(w.function)))){
      stop("Arguments for w.function does not contain those from calibmsm::calc_weights")
    }
    calc_weights <- w.function
  }

  ### Extract transition matrix from msdata object
  tmat <- attributes(data.mstate)$trans

  ### Assign the maximum state an individual may enter
  max.state <- max(data.mstate$to)

  ### Assign colnames to predicted transition probabilities (and in data.pred.plot)
  colnames(tp.pred) <- paste("tp.pred", 1:ncol(tp.pred), sep = "")
  if (!is.null(data.pred.plot)){
    colnames(data.pred.plot) <- paste("tp.pred", 1:ncol(tp.pred), sep = "")
  }

  ### Extract what states an individual can move into from state j (states with a non-zero predicted risk)
  valid.transitions <- which(colSums(tp.pred) != 0)

  ### Add the predicted risks, and the logit transformation of the predicted risks to data.raw and data.pred.plot
  tp.pred.logit <- log(tp.pred[,valid.transitions]/(1-tp.pred[,valid.transitions]))
  colnames(tp.pred.logit) <- paste("tp.pred.logit", valid.transitions, sep = "")
  data.raw <- data.frame(data.raw, tp.pred[,valid.transitions], tp.pred.logit)

  if (!is.null(data.pred.plot)){
    data.pred.plot.logit <- log(data.pred.plot[,valid.transitions]/(1-data.pred.plot[,valid.transitions]))
    colnames(data.pred.plot.logit) <- paste("tp.pred.logit", valid.transitions, sep = "")
    data.pred.plot <- data.frame(data.pred.plot[,valid.transitions], data.pred.plot.logit)

  }

  ### Extract which state individuals are in at time t
  ids.state.list <- vector("list", max.state)
  for (k in valid.transitions){
    ids.state.list[[k]] <- extract_ids_states(data.mstate, tmat, k, t)
  }

  ### Create a variable to say which state an individual was in at the time of interest
  ## Create list containing the relevant data
  v1 <- data.raw$id
  m1 <- outer(v1, ids.state.list, FUN = Vectorize('%in%'))
  state.poly <- lapply(split(m1, row(m1)), function(x) (1:max.state)[x])

  ## Change integer(0) values to NA's
  idx <- !sapply(state.poly, length)
  state.poly[idx] <- NA

  ## Add to data.raw
  data.raw <- dplyr::mutate(data.raw, state.poly = unlist(state.poly),
                            state.poly.fac = factor(state.poly))

  ### Create binary variables for each possible state that can be transitioned to
  ## Start by creating NA data.frame
  temp.dummy <- data.frame(matrix(NA, ncol = length(valid.transitions), nrow = nrow(data.raw)))

  ## Create dummy variables
  temp.dummy.calc <- stats::model.matrix(~state.poly.fac - 1, data.raw)

  ## Assign to temp.dummy, for rows where data.raw is not NA
  temp.dummy[!is.na(data.raw$state.poly), ] <- temp.dummy.calc

  ## Assign colnames
  colnames(temp.dummy) <- paste("state", valid.transitions, ".bin", sep = "")

  ### Add to dataset
  data.raw <- cbind(data.raw, temp.dummy)
  rm(temp.dummy)

  ### Identify individuals who are in state j at time s (will be used for landmarking)
  ids.state.js <- base::subset(data.mstate, from == j & Tstart <= s & s < Tstop) |>
    dplyr::select(id) |>
    dplyr::distinct(id) |>
    dplyr::pull(id)

  ### Reduce data.raw to landmarked dataset of individuals who are uncensored at time t,
  ### this is the set of predicted risks over which we plot calibration curves
  data.raw.lmk.js.uncens <- data.raw |> base::subset(id %in% ids.state.js) |> base::subset(!is.na(state.poly))


  ############################################################
  ### Write functions to estimate predicted observed risks ###
  ############################################################

  ### Functions are written to allow bootstrapping and have the following arguments:
  ### Arg1: data.in = dataset which we want to assess calibration in
  ### Arg2: indices = vector of indices to sample dataset
  ### Arg3: state we want to assess calibration for
  ### Arg4: data.in.uncens = dataset of uncensored individuals at time t, it is these we will generate the predicted observed risks for

  ###
  ### Function 1 uses loess smoothers
  ### Weights are calculated within the function to allow for proper bootstrapping
  ###
  calc_obs_loess_weights_func <- function(data.in, indices, state.k, data.in.uncens, data.pred.plot = NULL){

    ## Create bootstrapped dataset
    data.boot <- data.in[indices, ]

    ## Create landmarked dataset
    data.boot.lmk.js <-  data.boot |> base::subset(id %in% ids.state.js)

    ## Calculate weights
    ## Note this is done in the entire dataset data.boot, which has its own functionality (w.landmark.type) to landmark on j and s, or just s, before
    ## calculating the weights
    weights <- calc_weights(data.mstate = data.mstate,
                            data.raw = data.boot,
                            covs = w.covs,
                            t = t,
                            s = s,
                            landmark.type = w.landmark.type,
                            j = j,
                            max.weight = w.max,
                            stabilised = w.stabilised,
                            max.follow = w.max.follow,
                            ...)

    ## Add to data.boot
    data.boot.lmk.js <- dplyr::left_join(data.boot.lmk.js, dplyr::distinct(weights), by = dplyr::join_by(id))

    ## Reduce data.boot to individuals who are uncensored at time t
    data.boot.lmk.js.uncens <- base::subset(data.boot.lmk.js, !is.na(state.poly))

    ## Define equation
    eq.loess <- stats::formula(paste("state", state.k, ".bin ~ tp.pred", state.k, sep = ""))

    ## Fit model
    if (w.stabilised == FALSE){
      loess.model <- stats::loess(eq.loess,
                                  data = data.boot.lmk.js.uncens,
                                  weights = data.boot.lmk.js.uncens[, "ipcw"],
                                  span = loess.span,
                                  degree = loess.degree)
    } else if (w.stabilised == TRUE){
      loess.model <- stats::loess(eq.loess,
                                  data = data.boot.lmk.js.uncens,
                                  weights = data.boot.lmk.js.uncens[, "ipcw.stab"],
                                  span = loess.span,
                                  degree = loess.degree)
    }

    ## Create predicted observed probabilities. Create predictions for the vector of predicted probabilities for uncensored individuals from original dataset.
    loess.pred.obs <- predict(loess.model, newdata = data.in.uncens)

    return(loess.pred.obs)
  }

  ###
  ### Function 2 uses restricted cubic splines
  ### Weights are calculated within the function to allow for proper bootstrapping
  ###
  calc_obs_rcs_weights_func <- function(data.in, indices, state.k, data.in.uncens, data.pred.plot = NULL){

    ## Create bootstrapped dataset
    data.boot <- data.in[indices, ]

    # data.boot <- data.raw
    # state.k <- 1
    # data.in.uncens = data.raw.lmk.js.uncens

    ## Create landmarked dataset
    data.boot.lmk.js <-  data.boot |> base::subset(id %in% ids.state.js)

    ## Calculate weights
    ## Note this is done in the entire dataset data.boot, which has its own functionality (w.landmark.type) to landmark on j and s, or just s, before
    ## calculating the weights
    weights <- calc_weights(data.mstate = data.mstate,
                            data.raw = data.boot,
                            covs = w.covs,
                            t = t,
                            s = s,
                            landmark.type = w.landmark.type,
                            j = j,
                            max.weight = w.max,
                            stabilised = w.stabilised,
                            max.follow = w.max.follow,
                            ...)

    ## Add to data.boot
    data.boot.lmk.js <- dplyr::left_join(data.boot.lmk.js, dplyr::distinct(weights), by = dplyr::join_by(id))

    ## Reduce data.boot to individuals who are uncensored at time t
    data.boot.lmk.js.uncens <- base::subset(data.boot.lmk.js, !is.na(state.poly))

    ## Create restricted cubic splines for the cloglog of the linear predictor for the state of interst
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
    ## meaning the warning is not neccesary.
    suppressWarnings(
      if (w.stabilised == FALSE){
        rcs.model <- rms::lrm(eq.rcs,
                              data = data.boot.lmk.js.uncens,
                              weights = data.boot.lmk.js.uncens[, "ipcw"])
      } else if (w.stabilised == TRUE){
        rcs.model <- rms::lrm(eq.rcs,
                              data = data.boot.lmk.js.uncens,
                              weights = data.boot.lmk.js.uncens[, "ipcw.stab"])
      }
    )

    ## Create predicted observed probabilities.
    if (is.null(data.pred.plot)){

      ## Create predictions for the vector of predicted probabilities for uncensored individuals from original dataset.
      ## For this, need to calculate the correct splines for the original dataset.
      rcs.mat.data.in.uncens <- Hmisc::rcspline.eval(data.in.uncens[,paste("tp.pred.logit", state.k, sep = "")],knots = knots,inclx=T)
      colnames(rcs.mat.data.in.uncens) <- paste("rcs.x", 1:ncol(rcs.mat), sep = "")
      #attr(rcs.mat.data.in.uncens,"knots")

      ## Add the cubic splines for logit of the predicted probability to data.in.uncens
      data.in.uncens <- data.frame(cbind(data.in.uncens, rcs.mat.data.in.uncens))

      ## Calculate predicted observed probabilities
      rcs.pred.obs <- predict(rcs.model, newdata = data.in.uncens, type = "fitted")

    } else {
      ## Create spline terms for data.pred.plot
      rcs.mat.data.pred.plot <- Hmisc::rcspline.eval(data.pred.plot[,paste("tp.pred.logit", state.k, sep = "")],knots = knots,inclx=T)
      colnames(rcs.mat.data.pred.plot) <- paste("rcs.x", 1:ncol(rcs.mat), sep = "")

      ## Add to dataset
      data.pred.plot <- data.frame(cbind(data.pred.plot, rcs.mat.data.pred.plot))

      ## Calculate predicted observed probabilities
      rcs.pred.obs <- predict(rcs.model, newdata = data.pred.plot, type = "fitted")
    }

    return(rcs.pred.obs)
  }


  ###
  ### Function 3 uses loess smoothers
  ### Weights are not calculated within the function (assumes they have been inputted by user)
  ###
  calc_obs_loess_func <- function(data.in, indices, state.k, data.in.uncens, data.pred.plot = NULL){

    ## Create bootstrapped dataset
    data.boot <- data.in[indices, ]

    ## Create landmarked dataset
    data.boot.lmk.js <-  data.boot |> base::subset(id %in% ids.state.js)

    ## Reduce data.boot to individuals who are uncensored at time t
    data.boot.lmk.js.uncens <- base::subset(data.boot.lmk.js, !is.na(state.poly))

    ## Define equation
    eq.loess <- stats::formula(paste("state", state.k, ".bin ~ tp.pred", state.k, sep = ""))

    ## Fit model
    loess.model <- stats::loess(eq.loess,
                                data = data.boot.lmk.js.uncens,
                                weights = data.boot.lmk.js.uncens[, "ipcw"],
                                span = loess.span,
                                degree = loess.degree)

    ## Create predictions for the vector of predicted probabilities for people included in original the calibration curve (this is the individuals
    ## who are in the unbootstrapped dataset, and are uncensored at time t. This would be the vector of predicted probabilities for the
    ## calibration curve if no bootstrapping was done)

    ## Create predicted observed risks
    if (is.null(data.pred.plot)){
      loess.pred.obs <- predict(loess.model, newdata = data.in.uncens)
    } else {
      loess.pred.obs <- predict(loess.model, newdata = data.pred.plot)
    }

    return(loess.pred.obs)
  }


  ###
  ### Function 4 uses restricted cubic splines
  ### Weights are not calculated within the function (assumes they have been inputted by user)
  ###
  calc_obs_rcs_func <- function(data.in, indices, state.k, data.in.uncens, data.pred.plot = NULL){

    ## Create bootstrapped dataset
    data.boot <- data.in[indices, ]

    ## Create landmarked dataset
    data.boot.lmk.js <-  data.boot |> base::subset(id %in% ids.state.js)

    ## Reduce data.boot to individuals who are uncensored at time t
    data.boot.lmk.js.uncens <- base::subset(data.boot.lmk.js, !is.na(state.poly))

    ## Create restricted cubic splines for the cloglog of the linear predictor for the state of interst
    rcs.mat <- Hmisc::rcspline.eval(data.boot.lmk.js.uncens[,paste("tp.pred.logit", state.k, sep = "")],nk=rcs.nk,inclx=T)
    colnames(rcs.mat) <- paste("rcs.x", 1:ncol(rcs.mat), sep = "")
    knots <- attr(rcs.mat,"knots")

    ## Add the cubic splines for logit of the predicted probability to data.boot.lmk.js.uncens
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
    ## meaning the warning is not neccesary.
    suppressWarnings(
      rcs.model <- rms::lrm(eq.rcs,
                            data = data.boot.lmk.js.uncens,
                            weights = data.boot.lmk.js.uncens[, "ipcw"])
    )


    ## Create predicted observed probabilities. Create predictions for the vector of predicted probabilities form the original model
    ## So that we always plot over the same range.

    ## Calculate predicted observed probabilities
    if (is.null(data.pred.plot)){

      ## For this, need to calculate the correct splines for the original dataset.
      rcs.mat.data.in.uncens <- Hmisc::rcspline.eval(data.in.uncens[,paste("tp.pred.logit", state.k, sep = "")],knots = knots,inclx=T)
      colnames(rcs.mat.data.in.uncens) <- paste("rcs.x", 1:ncol(rcs.mat), sep = "")
      #attr(rcs.mat.data.in.uncens,"knots")

      ## Add the cubic splines for logit of the predicted probability to data.in.uncens
      data.in.uncens <- data.frame(cbind(data.in.uncens, rcs.mat.data.in.uncens))

      ## Calculate predicted observed probabilities
      rcs.pred.obs <- predict(rcs.model, newdata = data.in.uncens, type = "fitted")
    } else {
      ## Create spline terms for data.pred.plot
      rcs.mat.data.pred.plot <- Hmisc::rcspline.eval(data.pred.plot[,paste("tp.pred.logit", state.k, sep = "")],knots = knots,inclx=T)
      colnames(rcs.mat.data.pred.plot) <- paste("rcs.x", 1:ncol(rcs.mat), sep = "")

      ## Add to dataset
      data.pred.plot <- data.frame(cbind(data.pred.plot, rcs.mat.data.pred.plot))

      ## Calculate predicted observed probabilities
      rcs.pred.obs <- predict(rcs.model, newdata = data.pred.plot, type = "fitted")
    }

    return(rcs.pred.obs)
  }

  ###################################
  ### FUNCTIONS HAVE BEEN DEFINED ###
  ###################################

  ### Assign function to calculate predicted observed risks based on user input
  if (is.null(weights)){
    if (curve.type == "loess"){
      calc_obs_func <- calc_obs_loess_weights_func
    } else if (curve.type == "rcs"){
      calc_obs_func <- calc_obs_rcs_weights_func
    }
  } else {
    if (curve.type == "loess"){
      calc_obs_func <- calc_obs_loess_func
    } else if (curve.type == "rcs"){
      calc_obs_func <- calc_obs_rcs_func
    }
  }

  ### Create object to store output
  if (is.null(transitions.out)){
    transitions.out <- valid.transitions
  }

  output.object <- vector("list", length(transitions.out))
  output.boot.object <- vector("list", length(transitions.out))
  names(output.object) <- paste("state", transitions.out, sep = "")
  names(output.boot.object) <- paste("state", transitions.out, sep = "")

  ### Loop through and fit models
  for (k in 1:length(transitions.out)){

    ### Assign state of interest
    state.k <- transitions.out[k]

    ### Calculate predicted observed probabilities for calibration curve (note the chosen set of indices samples every patient once)
    rcs.pred.obs <- calc_obs_func(data.in = data.raw, indices = 1:nrow(data.raw), state.k = state.k,
                                  data.in.uncens = data.raw.lmk.js.uncens, data.pred.plot = data.pred.plot)

    ### Assign output
    if (is.null(data.pred.plot)){
      output.object[[k]] <- data.frame("id" = data.raw.lmk.js.uncens[, "id"],
                                       "pred" = data.raw.lmk.js.uncens[, paste("tp.pred", transitions.out[k], sep = "")],
                                       "obs" = rcs.pred.obs)
    } else {
      output.object[[k]] <- data.frame(
        "pred" = data.pred.plot[, paste("tp.pred", transitions.out[k], sep = "")],
        "obs" = rcs.pred.obs)
    }

    ### Calculate confidence intervals (only if user-specifies them, and user didn't input their own weights, which would lead to incorrect confidence intervals)
    #     if (CI != FALSE & (is.null(weights))){
    if (CI != FALSE & CI.type == "bootstrap"){

      ### Define alpha for CI's
      alpha <- (1-CI/100)/2

      ### Run bootstrapping
      if (is.null(data.pred.plot)){
        boot.obs <- boot::boot(data.raw, calc_obs_func, R = CI.R.boot, state.k = transitions.out[k], data.in.uncens = data.raw.lmk.js.uncens)$t
      } else if (!is.null(data.pred.plot)){
        boot.obs <- boot::boot(data.raw, calc_obs_func, R = CI.R.boot, state.k = transitions.out[k], data.in.uncens = data.raw.lmk.js.uncens,
                               data.pred.plot = data.pred.plot)$t
      }

      output.boot.object[[k]] <- boot.obs

      ### Extract confidence bands
      lower <- apply(boot.obs, 2, stats::quantile, probs = alpha, na.rm = TRUE)
      upper <- apply(boot.obs, 2, stats::quantile, probs = 1-alpha, na.rm = TRUE)

      ### Produce a warning if any NA values
      if(sum(is.na(boot.obs)) > 0){
        print(paste("WARNING, SOME BOOTSTRAPPED OBSERVED RISKS WERE NA FOR STATE", transitions.out[k]))
        print(paste("THERE ARE ", sum(apply(boot.obs, 1, function(x) {sum(is.na(x)) > 0})), " ITERATIONS WITH NA's FOR OBSERVED RISKS"))
        print(paste("THE MEAN NUMBER OF NA's IN EACH ITERATION IS", mean(apply(boot.obs, 1, function(x) {sum(is.na(x))}))))
      }

      ### Assign output
      if (is.null(data.pred.plot)){
        output.object[[k]] <- data.frame("id" = data.raw.lmk.js.uncens[, "id"],
                                         "pred" = data.raw.lmk.js.uncens[, paste("tp.pred", transitions.out[k], sep = "")],
                                         "obs" = rcs.pred.obs,
                                         "obs.lower" = lower,
                                         "obs.upper" = upper)
      } else if (!is.null(data.pred.plot)){
        output.object[[k]] <- data.frame(
          "pred" = data.raw.lmk.js.uncens[, paste("tp.pred", transitions.out[k], sep = "")],
          "obs" = rcs.pred.obs,
          "obs.lower" = lower,
          "obs.upper" = upper)
      }

    } else if (CI != FALSE & CI.type == "parametric"){

      ### PLACE HOLDER FOR FUTURE INCLUSION OF PARAMETRIC CONFIDENCE INTERVALS

    }
  }

  ### Create metadata object
  metadata <- list("valid.transitions" = as.numeric(valid.transitions),
                   "assessed.transitions" = as.numeric(transitions.out),
                   "CI" = CI,
                   "CI.R.boot" = CI.R.boot,
                   "curve.type" = curve.type,
                   "j" = j,
                   "s" = s,
                   "t" = t)

  ### Crate a combined output object with metadata, as well as plot data
  output.object.comb <- list("plotdata" = output.object, "metadata" = metadata, "bootdata" = output.boot.object)

  ### Assign calib_blr class
  attr(output.object.comb, "class") <- "calib_blr"

  return(output.object.comb)
}

data("ebmtcal")
head(ebmtcal)

data("msebmtcal")
head(msebmtcal)

data("tps0")
head(tps0)

tp.pred.s0 <- tps0 |>
  dplyr::filter(j == 1) |>
  dplyr::select(any_of(paste("pstate", 1:6, sep = "")))
t.eval <- 1826
dat.calib.blr <-
  calib_blr_test(data.mstate = msebmtcal,
            data.raw = ebmtcal,
            j = 1,
            s = 0,
            t = t.eval,
            tp.pred = tp.pred.s0,
            curve.type = "loess",
            w.covs = c("year", "agecl", "proph", "match"),
            CI = 95,
            CI.type = "bootstrap",
            CI.R.boot = 200,
            transitions.out = c(3))

plot(dat.calib.blr)
str(dat.calib.blr[["bootdata"]])
str(dat.calib.blr[["plotdata"]]$state3)
which.min(dat.calib.blr[["plotdata"]]$state3$obs)
state3boot <- dat.calib.blr[["bootdata"]]$state3
state3boot[,1741]
