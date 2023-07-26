#' Calculate inverse probability of censoring weights at time `t`.
#'
#' @description
#' Estimates the inverse probability of censoring weights by fitting a cox-propotinal hazards model
#' in a landmark cohort of individuals. Primarily used internally, this function has
#' been exported to allow users to reproduce results in the vignette when
#' estimating confidence intervals using bootstrapping manually.
#'
#' @param data.mstate Validation data in msdata format
#' @param data.raw Validation data in data.frame (one row per individual)
#' @param t Follow up time at which to calculate weights
#' @param s Landmark time at which predictions were made
#' @param landmark.type Whether weights are estimated in all individuals uncensored at time s ('all') or only in individuals uncensored and in state j at time s ('state')
#' @param j Landmark state at which predictions were made (only required in landmark.type = 'state')
#' @param covs Character vector of variable names to adjust for when calculating inverse probability of censoring weights
#' @param max.weight Maximum bound for weights
#' @param stabilised Indicates whether weights should be stabilised or not
#' @param max.follow Maximum follow up for model calculating inverse probability of censoring weights. Reducing this to `t` + 1 may aid in the proportional hazards assumption being met in this model.
#'
#' @returns A dataframe with three columns. `id` corresponds to the patient ids from `data.raw`. `ipcw` contains the inverse probability
#' of censoring weights (specifically the inverse of the probability of being uncesored). `pcw = 1/ipcw`. If `stabilised = TRUE` was specified,
#' a fourth variable `ipcw.stab` will be returned, which is the stabilised inverse probability of censoring weights.
#'
#' @details
#' Fits a cox proportional hazards model to individuals in a landmark cohort, predicting the probability of being censored
#' at time `t`. This landmark cohort may either be all individuals uncensored at time `s`, or those uncensored
#' and in state `j` at time `s`. All predictors in `w.covs` are assumed to have a linear effect on the hazard.
#' Weights are estimated for all individuals in `data.raw`, even if they will not be used in the analysis as they do not meet the landmarking
#' requirements. If an individual enters an absorbing state prior to `t`, we estimate the probability of being censored
#' before the time of entry into the absorbing state, rather than at `t`. Details on all the above this are provided in
#' (vignette OVERVIEW-AND-THEORY-XXXX)
#'
#'
#' @examples
#' # Estimate inverse probability of censoring weights for individual in cohort ebmtcal.
#' # Specifically the probability of being uncensored at t = 1826 days.
#' # Weights are estimated using a model fitted in all individuals uncensored at time s = 0.
#' weights.manual <-
#' calc_weights(data.mstate = msebmtcal,
#'   data.raw = ebmtcal,
#'   covs = c("year", "agecl", "proph", "match"),
#'   t = 1826,
#'   s = 0,
#'   landmark.type = "state",
#'   j = 1)
#'
#'   str(weights.manual)
#'
#' @export
calc_weights <- function(data.mstate, data.raw, covs = NULL, t, s, landmark.type = "state", j = NULL, max.weight = 10, stabilised = FALSE, max.follow = NULL){

  ### Modify everybody to be censored after time t, if a max.follow has been specified
  if(!is.null(max.follow)){

    ### Stop if max follow is smaller than t
    if (t < max.follow){
      stop("Max follow cannot be smaller than t")
    } else {
      data.raw <- dplyr::mutate(data.raw,
                                dtcens.s = dplyr::case_when(dtcens < max.follow + 2 ~ dtcens.s,
                                                            dtcens >= max.follow + 2 ~ 0),
                                dtcens = dplyr::case_when(dtcens < max.follow + 2 ~ dtcens,
                                                          dtcens >= max.follow + 2 ~ max.follow + 2))
    }
  }

  ### Create a new outcome, which is the time until censored from s
  data.raw$dtcens.modified <- data.raw$dtcens - s

  ### Save a copy of data.raw
  data.raw.save <- data.raw

  ### If landmark.type = "state", calculate weights only in individuals in state j at time s
  ### If landmark.type = "all", calculate weights in all uncensored individuals at time s (note that this excludes individuals
  ### who have reached absorbing states, who have been 'censored' from the survival distribution is censoring)
  if (landmark.type == "state"){
    ### Identify individuals who are uncensored in state j at time s
    ids.uncens <- base::subset(data.mstate, from == j & Tstart <= s & s < Tstop) |>
      dplyr::select(id) |>
      dplyr::distinct(id) |>
      dplyr::pull(id)

  } else if (landmark.type == "all"){
    ### Identify individuals who are uncensored time s
    ids.uncens <- base::subset(data.mstate, Tstart <= s & s < Tstop) |>
      dplyr::select(id) |>
      dplyr::distinct(id) |>
      dplyr::pull(id)

  }

  ### Subset data.mstate and data.raw to these individuals
  data.mstate <- data.mstate |> base::subset(id %in% ids.uncens)
  data.raw <- data.raw |> base::subset(id %in% ids.uncens)

  ###
  ### Create models for censoring in order to calculate the IPCW weights
  ### Seperate models for estimating the weights, and stabilising the weights (intercept only model)
  ###
  if (!is.null(covs)){
    ### A model where we adjust for predictor variables
    cens.model <- survival::coxph(stats::as.formula(paste("survival::Surv(dtcens.modified, dtcens.s) ~ ",
                                                          paste(covs, collapse = "+"),
                                                          sep = "")),
                                  data = data.raw)

    ### Intercept only model (numerator for stabilised weights)
    cens.model.int <- survival::coxph(stats::as.formula(paste("survival::Surv(dtcens.modified, dtcens.s) ~ 1",
                                                              sep = "")),
                                      data = data.raw)
  } else if (is.null(covs)){
    ### If user has not input any predictors for estimating weights, the model for estimating the weights is the intercept only model (i.e. Kaplan Meier estimator)

    ### Intercept only model (numerator for stabilised weights)
    cens.model.int <- survival::coxph(stats::as.formula(paste("survival::Surv(dtcens.modified, dtcens.s) ~ 1",
                                                              sep = "")),
                                      data = data.raw)
    ### Assign cens.model to be the same
    cens.model <- cens.model.int


  }

  ### Calculate a data frame containing probability of censored and uncenosred at each time point
  ### The weights will be the probability of being uncensored, at the time of the event for each individual

  ## Extract baseline hazard
  data.weights <- survival::basehaz(cens.model, centered = FALSE)
  ## Add lp to data.raw.save
  data.raw.save$lp <- stats::predict(cens.model, newdata = data.raw.save, type = "lp", reference = "zero")

  ### Create weights for the cohort at time t - s
  ### Note for individuals who reached an absorbing state, we take the probability of them being uncensored at the time of reached the
  ### abosrbing state. For individuals still alive, we take the probability of being uncensored at time t - s.

  ### Get location of individuals who entered absorbing states or were censored prior to evaluation time
  obs.absorbed.prior <- which(data.raw.save$dtcens <= t & data.raw.save$dtcens.s == 0)
  obs.censored.prior <- which(data.raw.save$dtcens <= t & data.raw.save$dtcens.s == 1)

  ###
  ### Now create unstabilised probability of (un)censoring weights
  ### Note that weights are the probability of being uncensored, so if an individual has low probability of being uncesored,
  ### the inervse of this will be big, weighting them strongly
  ###

  ### First assign all individuals a weight of the probability of being uncensored at time t
  ### This is the linear predictor times the cumulative hazard at time t, and appropriate transformation to get a risk
  data.raw.save$pcw <- as.numeric(exp(-exp(data.raw.save$lp)*data.weights$hazard[max(which(data.weights$time <= t - s))]))

  ## Write a function which will extract the uncensored probability for an individual with linear predictor lp at a given time t
  prob.uncens.func <- function(input){

    ## Assign t and person_id
    t <- input[1]
    lp <- input[2]

    if (t <= 0){
      return(NA)
    } else if (t > 0){
      ## Get hazard at appropriate time
      if (t < min(data.weights$time)){
        bhaz.t <- 0
      } else if (t >= min(data.weights$time)){
        bhaz.t <- data.weights$hazard[max(which(data.weights$time <= t))]
      }

      ## Return risk
      return(exp(-exp(lp)*bhaz.t))
    }
  }

  ### Apply this function to all the times at which individuals have entered an absorbing state prior to censoring
  data.raw.save$pcw[obs.absorbed.prior] <- apply(data.raw.save[obs.absorbed.prior, c("dtcens.modified", "lp")], 1, FUN = prob.uncens.func)

  ### For individuals who were censored prior to t, assign the weight as NA
  data.raw.save$pcw[obs.censored.prior] <- NA

  ### Invert these
  data.raw.save$ipcw <- 1/data.raw.save$pcw

  ###
  ### Stabilise these weights dependent on user-input
  ###
  if (stabilised == TRUE){

    ## Extract baseline hazard
    data.weights.numer <- survival::basehaz(cens.model.int, centered = TRUE)

    ### Assign all individuals a weight of the probability of being uncesored at time t
    data.raw.save$pcw.numer <- as.numeric(exp(-data.weights.numer$hazard[max(which(data.weights.numer$time <= t - s))]))

    ### Create stabilised weight
    data.raw.save$ipcw.stab <- data.raw.save$pcw.numer*data.raw.save$ipcw
  }

  ### Finally cap these at 10 and create output object

  ### Create output object
  if (stabilised == FALSE){
    data.raw.save$ipcw <- pmin(data.raw.save$ipcw, max.weight)
    output.weights <- data.frame("id" = data.raw.save$id, "ipcw" = data.raw.save$ipcw, "pcw" = data.raw.save$pcw)
  } else if (stabilised == TRUE){
    data.raw.save$ipcw <- pmin(data.raw.save$ipcw, max.weight)
    data.raw.save$ipcw.stab <- pmin(data.raw.save$ipcw.stab, max.weight)
    output.weights <- data.frame("id" = data.raw.save$id, "ipcw" = data.raw.save$ipcw, "ipcw.stab" = data.raw.save$ipcw.stab, "pcw" = data.raw.save$pcw)
  }

  return(output.weights)

}
