calib_mlr_SAVE <- function(data.mstate,
                      data.raw,
                      j,
                      s,
                      t,
                      tp.pred,
                      smoother.type = "sm.ps",
                      ps.int = 4,
                      degree = 3,
                      s.df = 4,
                      niknots = 4,
                      weights = NULL,
                      w.function = NULL,
                      w.covs = NULL,
                      w.landmark.type = "state",
                      w.max = 10,
                      w.stabilised = FALSE,
                      w.max.follow = NULL, ...){

  ### Stop if patients in data.raw are not in data.mstate
  if (!base::all(unique(data.raw$id) %in% unique(data.mstate$id))){
    stop("All patients in data.raw are not contained in data.mstate. Landmarking cannot be applied.")
  }

  ### Warning if patients in data.mstate are not in data.raw
  if (!base::all(unique(data.mstate$id) %in% unique(data.raw$id))){
    warning("All patients in data.mstate are not contained in data.raw. Landmarking can still be applied, but potential mismatch in these two datasets?")
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

  ### Extract transition matrix from msdata object
  tmat <- attributes(data.mstate)$trans

  ### Assign the maximum state an individual may enter
  max.state <- max(data.mstate$to)

  ### Assign colnames to tp.pred
  colnames(tp.pred) <- paste("tp.pred", 1:ncol(tp.pred), sep = "")

  ### Extract what states an individual can move into from state j (states with a non-zero predicted risk)
  valid.transitions <- which(colSums(tp.pred) != 0)

  ### Add linear predictors from a multinomial framework
  ### Start by reducing to non-zero columns
  tp.pred.mlr <- tp.pred[,valid.transitions]

  ### Calculate linear predictors
  tp.pred.mlr <- log(tp.pred.mlr[,2:ncol(tp.pred.mlr)]/tp.pred.mlr[,1])
  colnames(tp.pred.mlr) <- paste("mlr.lp", 1:(ncol(tp.pred.mlr)), sep = "")

  ### Add to data frame
  data.raw <- data.frame(data.raw, tp.pred, tp.pred.mlr)

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
  data.raw <- dplyr::mutate(data.raw, state.poly = base::unlist(state.poly),
                            state.poly.fac = base::factor(state.poly))

  ### Identify individuals who are in state j at time s (will be used for landmarking)
  ids.state.js <- base::subset(data.mstate, from == j & Tstart <= s & s < Tstop) |>
    dplyr::select(id) |>
    dplyr::distinct(id) |>
    dplyr::pull(id)

  ### Reduce data.raw to landmarked dataset of individuals who are uncensored at time t,
  ### this is the set of predicted risks over which we plot calibration curves
  data.raw.lmk.js.uncens <- data.raw |> base::subset(id %in% ids.state.js) |> base::subset(!is.na(state.poly))

  ### Calculate weights if not specified manually
  if (is.null(weights)){

    ### Assign custom function for estimating weights, if specified
    if (!is.null(w.function)){
      ### stop if w.function doesn't have correct arguments
      if(!all(names(formals(calc_weights)) %in% names(formals(w.function)))){
        stop("Arguments for w.function does not contain those from calibmsm::calc_weights")
      }
      calc_weights <- w.function
    }

    ### Estimate the weights
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
    ## Add to data.raw
    data.raw.lmk.js.uncens <- dplyr::left_join(data.raw.lmk.js.uncens, dplyr::distinct(weights), by = dplyr::join_by(id))
  }

  ### Define equation
  eq.LHS <- paste("state.poly.fac ~ ")
  if (smoother.type == "s"){
    eq.RHS <- paste("s(mlr.lp", 1:ncol(tp.pred.mlr), ", df = s.df)", sep = "", collapse = "+")
  } else if (smoother.type == "sm.ps"){
    eq.RHS <- paste("sm.ps(mlr.lp", 1:ncol(tp.pred.mlr), ", ps.int = ps.int, degree = degree)", sep = "", collapse = "+")
  } else if (smoother.type == "sm.os"){
    eq.RHS <- paste("sm.os(mlr.lp", 1:ncol(tp.pred.mlr), ", niknots = niknots)", sep = "", collapse = "+")
  }
  eq.mlr <- stats::as.formula(paste(eq.LHS, eq.RHS, sep =""))

  ### Assign reference category
  ref.cat <- paste(valid.transitions[1])

  ### Apply nominal recalibration framework with vector spline smoothers
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

  ### Create metadata object
  metadata <- list("valid.transitions" = valid.transitions,
                   "j" = j,
                   "s" = s,
                   "t" = t)

  ### Crate a combined output object with metadata, as well as plot data
  output.object.comb <- list("plotdata" = output.object2, "metadata" = metadata)

  ### Assign calib_blr class
  attr(output.object.comb, "class") <- "calib_mlr"

  return(output.object.comb)

}


