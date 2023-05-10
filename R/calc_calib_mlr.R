#' Create data for calibration curves using a multinomial logistic regression framework with inverse probability of censoring weights
#'
#' @description
#' Creates the underlying data for the calibration plots. Observed event
#' probabilities at time `t.eval` are estimated for inputted predicted
#' transition probabilities `tp.pred` out of state `j` at time `s`.
#' `calc_calib_mlr` estimates calibration scatter plots using a multinomial logistic
#' framework in combination with landmarking and inverse probability of
#' censoring weights.
#'
#' Two datasets for the same cohort of inidividuals must be provided. A `msdata`
#' format dataset generated using the `mstate` package. A `data.frame` with one
#' row per individual, relevant variables for estimating the weights, and a time
#' until censoring varaible (`dtcens`) and indicator (`dtcens.s`). Weights are
#' estimated using a cox-proportional hazard model and assuming linear
#' functional form of the variables defined in `w.covs`. We urge users to
#' specify their own modwl for estimating the weights. Confidence intervals for
#' the calibration scatter plots cannot be produced as it is currently unclear how
#' to present such data.
#'
#' @param data.mstate Validation data in `msdata` format
#' @param data.raw Validation data in `data.frame` (one row per individual)
#' @param j Landmark state at which predictions were made
#' @param s Landmark time at which predictions were made
#' @param t.eval Follow up time at which calibration is to be assessed
#' @param tp.pred Vector of predicted transition probabilities at time t.eval
#' @param s.df degrees of freedom of vector spline (see \code{\link[VGAM]{s}})
#' @param smoother.type Type of smoothing applied. Takes values `s` (see \code{\link[VGAM]{s}}), `sm.ps` (see \code{\link[VGAM]{sm.ps}}) or `sm.os` (see \code{\link[VGAM]{sm.os}}).
#' @param ps.int the number of equally-spaced B spline intervals in the vector spline smoother (see \code{\link[VGAM]{sm.ps}})
#' @param degree the degree of B-spline basis in the vector spline smoother (see \code{\link[VGAM]{sm.ps}})
#' @param niknots number of interior knots (see \code{\link[VGAM]{sm.os}})
#' @param weights Vector of inverse probability of censoring weights
#' @param w.covs Character vector of variable names to adjust for when calculating inverse probability of censoring weights
#' @param w.landmark.type Whether weights are estimated in all individuals uncensored at time s ('all') or only in individuals uncensored and in state j at time s ('state')
#' @param w.max Maximum bound for inverse probability of censoring weights
#' @param w.stabilised Indicates whether inverse probability of censoring weights should be stabilised or not
#' @param w.max.follow Maximum follow up for model calculating inverse probability of censoring weights. Reducing this to `t.eval` + 1 may aid in the proportional hazards assumption being met in this model.
#'
#' @details
#' Observed event probabilities at time `t.eval` are estimated for predicted
#' transition probabilities `tp.pred` out of state `j` at time `s`.
#' `calc_calib_mlr` estimates calibration scatter plots uses a technique for assessing the calibration of multinomial logistic
#' regression models, namely the nominal calibration framework of #' van Hoorde et al. (2014, 2015).
#' Landmarking is applied to only assess calibration in individuals who are uncensored
#' and in state `j` at time `s`. Censoring is dealt with using inverse probability of
#' censoring weights.
#'
#' Two datasets for the same cohort of inidividuals must be provided. Firstly `data.mstate` must be a dataset of class `msdata`,
#' generated using the \code{[mstate]} package. This dataset is used to apply the landmarking. Secondly, `data.raw` must be
#' a `data.frame` with one row per individual, containing the desired variables for estimating the weights, and variables for the time
#' until censoring (`dtcens`), and an indicator for censoring `dtcens.s`, where (`dtcens.s = 1`) if
#' an individual is censored at time `dtcens`, and `dtcens.s = 0` otherwise. When an individual
#' enters an absorbing state, this prevents censoring from happening (i.e. dtcens.s = 0). Unless the user specifies
#' the weights using `weights`, the weights are
#' estimated using a cox-proportional hazard model, assuming a linear
#' functional form of the variables defined in `w.covs`. We urge users to
#' specify their own model for estimating the weights. The `weights` argument
#' must be a vector with length equal to the number of rows of `data.raw`.
#'
#' Confidence intervals cannot be generated for the calibration scatter plots. While
#' confidence intervals could be estimated for each data point, it is not clear how
#' these would be plotted cohesively.
#'
#' Calibration plots cannot be produced for specific transitions (i.e. `transitions.out`
#' in \code{\link{calc_calib_blr}}) because the nominal calibration framework (van Hoorde et al., 2014, 2015) assesses
#' the calibration of all states simultaneously.
#'
#' The calibration scatter plots can be plotted using \code{\link{plot.calib_mlr}}.
#'
#'#' @returns \code{\link{calc_calib_mlr}} returns a list containing two elements:
#' \code{plotdata} and \code{metadata}. The \code{plotdata} element contains the
#' data for the calibration scatter plots This will itself be a list with each element
#' containing the data for the transition probabilities into each of the possible
#' states. Each list element contains patient ids (\code{id}), the predicted
#' transition probabilities (\code{pred}) and the estimated observed event
#' probabilities (\code{obs}). The \code{metadata} element contains metadata
#' including a vector of the possible transitions and other user specified information.
#'
#' @examples
#'
#' # Estimate calibration scatter plots for the predicted transition
#' # probabilities at time t.eval = 1826, when predictions were made at time
#' # s = 0 in state j = 1. These predicted transition probabilities are stored in tps0.
#'
#' # Extract the predicted transition probabilities out of state j = 1
#' tp.pred <- dplyr::select(dplyr::filter(tps0, j == 1), any_of(paste("pstate", 1:6, sep = "")))
#'
#' # Now estimate the observed event probabilities for each possible transition.
#'
#' dat.calib.mlr <-
#' calc_calib_mlr(data.mstate = msebmtcal,
#'  data.raw = ebmtcal,
#'  j=1,
#'  s=0,
#'  t.eval = 1826,
#'  tp.pred = tp.pred,
#'  w.covs = c("year", "agecl", "proph", "match"))
#'
#' # The data for each calibration scatter plots are stored in the "plotdata"
#' # list element.
#' str(dat.calib.mlr)
#'

#' @export
calc_calib_mlr <- function(data.mstate, data.raw, j, s, t.eval, tp.pred, smoother.type = "sm.ps", ps.int = 4, degree = 3, s.df = 4, niknots = 4,
                           weights = NULL, w.covs, w.landmark.type = "state", w.max = 10, w.stabilised = FALSE, w.max.follow = NULL){

  # data.mstate <- msebmtcal
  # data.raw <- ebmtcal
  # j<-1
  # s<-0
  # t.eval <- 1826
  # tp.pred = tps0 %>% dplyr::filter(j == 1) %>% dplyr::select(any_of(paste("pstate", 1:6, sep = "")))
  # ps.int <- 4
  # degree <- 3
  # weights <- NULL
  # w.covs = c("year", "agecl", "proph", "match")
  # w.landmark.type = "all"
  # w.max = 10
  # w.stabilised = FALSE
  # smoother.type <- "sm.ps"

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

  ### Extract which state individuals are in at time t.eval
  ids.state.list <- vector("list", max.state)
  for (k in valid.transitions){
    ids.state.list[[k]] <- extract_ids_states(data.mstate, tmat, k, t.eval)
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
  ids.state.js <- base::subset(data.mstate, from == j & Tstart <= s & s < Tstop) %>%
    dplyr::select(id) %>%
    dplyr::distinct(id) %>%
    dplyr::pull(id)

  ### Reduce data.raw to landmarked dataset of individuals who are uncensored at time t.eval,
  ### this is the set of predicted risks over which we plot calibration curves
  data.raw.lmk.js.uncens <- data.raw %>% base::subset(id %in% ids.state.js) %>% base::subset(!is.na(state.poly))

  ## Calculate weights
  ## Note this is done in the entire dataset data.raw, which has its own functionality (w.landmark.type) to landmark on j and s, or just s, before
  ## calculating the weights
  weights <- calc_weights(data.mstate = data.mstate,
                          data.raw = data.raw,
                          covs = w.covs,
                          t.eval = t.eval,
                          s = s,
                          landmark.type = w.landmark.type,
                          j = j,
                          max.weight = w.max,
                          stabilised = w.stabilised,
                          max.follow = w.max.follow)

  ## Add to data.boot
  data.raw.lmk.js.uncens <- dplyr::left_join(data.raw.lmk.js.uncens, dplyr::distinct(weights), by = dplyr::join_by(id))

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

  ### Get plotdata in same format as calc_calib_blr
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
                   "t.eval" = t.eval)

  ### Crate a combined output object with metadata, as well as plot data
  output.object.comb <- list("plotdata" = output.object2, "metadata" = metadata)

  ### Assign calib_blr class
  attr(output.object.comb, "class") <- "calib_mlr"

  return(output.object.comb)

}

#' @export
summary.calib_mlr <- function(object, ...) {

  cat("There were non-zero predicted transition probabilities into states ",
      paste(object[["metadata"]]$valid.transitions, collapse = ","),  sep = " ")

  cat("\n\nCalibration was assessed at time ", object[["metadata"]]$t.eval, " and calibration was assessed in a landmarked cohort of individuals in state j = ", object[["metadata"]]$j,
      " at time s = ", object[["metadata"]]$s, sep = "")

  cat("\n\nThe estimated calibration scatter plots are stored in list element `plotdata`:\n\n")

  print(lapply(object[["plotdata"]], "head"))

}

