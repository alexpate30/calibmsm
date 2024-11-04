#' Calculate inverse probability of censoring weights at time `t`.
#'
#' @description
#' Estimates the inverse probability of censoring weights by fitting a cox-propotinal hazards model
#' in a landmark cohort of individuals. Primarily used internally, this function has
#' been exported to allow users to reproduce results in the vignette when
#' estimating confidence intervals using bootstrapping manually.
#'
#' @param data_ms Validation data in msdata format
#' @param data_raw Validation data in data.frame (one row per individual)
#' @param t Follow up time at which to calculate weights
#' @param s Landmark time at which predictions were made
#' @param landmark_type Whether weights are estimated in all individuals uncensored at time s ('all') or only in individuals uncensored and in state j at time s ('state')
#' @param j Landmark state at which predictions were made (only required in landmark_type = 'state')
#' @param covs Character vector of variable names to adjust for when calculating inverse probability of censoring weights
#' @param max_weight Maximum bound for weights
#' @param stabilised Indicates whether weights should be stabilised or not
#' @param max_follow Maximum follow up for model calculating inverse probability of censoring weights. Reducing this to `t` + 1 may aid in the proportional hazards assumption being met in this model.
#'
#' @returns A data frame with two columns. `id` corresponds to the patient ids from `data_raw`. `ipcw` contains the inverse probability
#' of censoring weights (specifically the inverse of the probability of being uncesored). If `stabilised = TRUE` was specified,
#' a third variable `ipcw_stab` will be returned, which is the stabilised inverse probability of censoring weights.
#'
#' @details
#' Estimates inverse probability of censoring weights (Hernan M, Robins J, 2020).
#' Fits a cox proportional hazards model to individuals in a landmark cohort, predicting the probability of being censored
#' at time `t`. This landmark cohort may either be all individuals uncensored at time `s`, or those uncensored
#' and in state `j` at time `s`. All predictors in `w_covs` are assumed to have a linear effect on the hazard.
#' Weights are estimated for all individuals in `data_raw`, even if they will not be used in the analysis as they do not meet the landmarking
#' requirements. If an individual enters an absorbing state prior to `t`, we estimate the probability of being censored
#' before the time of entry into the absorbing state, rather than at `t`. Details on all the above this are provided in
#' vignette \emph{overview}.
#'
#' @references
#' Hernan M, Robins J (2020). “12.2 Estimating IP weights via modeling.” In \emph{Causal Inference:
#' What If}, chapter 12.2. Chapman Hall/CRC, Boca Raton.
#'
#' @examples
#' # Estimate inverse probability of censoring weights for individual in cohort ebmtcal.
#' # Specifically the probability of being uncensored at t = 1826 days.
#' # Weights are estimated using a model fitted in all individuals uncensored at time s = 0.
#' weights_manual <-
#' calc_weights(data_ms = msebmtcal,
#'   data_raw = ebmtcal,
#'   covs = c("year", "agecl", "proph", "match"),
#'   t = 1826,
#'   s = 0,
#'   landmark_type = "state",
#'   j = 1)
#'
#'  str(weights_manual)
#'
#' @export
calc_weights <- function(data_ms, data_raw, covs = NULL, t, s, landmark_type = "state", j = NULL, max_weight = 10, stabilised = FALSE, max_follow = NULL){

  ### Modify everybody to be censored after time t, if a max_follow has been specified
  if(!is.null(max_follow)){

    ### Stop if max follow is smaller than t
    if (max_follow < t){
      stop("Max follow cannot be smaller than t")
    } else {
      data_raw <- dplyr::mutate(data_raw,
                                dtcens_s = dplyr::case_when(dtcens < max_follow + 2 ~ dtcens_s,
                                                            dtcens >= max_follow + 2 ~ 0),
                                dtcens = dplyr::case_when(dtcens < max_follow + 2 ~ dtcens,
                                                          dtcens >= max_follow + 2 ~ max_follow + 2))
    }
  }

  ### Create a new outcome, which is the time until censored from s
  data_raw$dtcens_modified <- data_raw$dtcens - s

  ### Save a copy of data_raw
  data_raw_save <- data_raw

  ### If landmark_type = "state", calculate weights only in individuals in state j at time s
  ### If landmark_type = "all", calculate weights in all uncensored individuals at time s (note that this excludes individuals
  ### who have reached absorbing states, who have been 'censored' from the survival distribution is censoring)
  if (landmark_type == "state"){
    ### Identify individuals who are uncensored in state j at time s
    ids_uncens <- base::subset(data_ms, from == j & Tstart <= s & s < Tstop) |>
      dplyr::select(id) |>
      dplyr::distinct(id) |>
      dplyr::pull(id)

  } else if (landmark_type == "all"){
    ### Identify individuals who are uncensored time s
    ids_uncens <- base::subset(data_ms, Tstart <= s & s < Tstop) |>
      dplyr::select(id) |>
      dplyr::distinct(id) |>
      dplyr::pull(id)

  }

  ### Subset data_ms and data_raw to these individuals
  data_ms <- data_ms |> base::subset(id %in% ids_uncens)
  data_raw <- data_raw |> base::subset(id %in% ids_uncens)

  ###
  ### Create models for censoring in order to calculate the IPCW weights
  ### Seperate models for estimating the weights, and stabilising the weights (intercept only model)
  ###
  if (!is.null(covs)){
    ### A model where we adjust for predictor variables
    cens_model <- survival::coxph(stats::as.formula(paste("survival::Surv(dtcens_modified, dtcens_s) ~ ",
                                                          paste(covs, collapse = "+"),
                                                          sep = "")),
                                  data = data_raw)

    ### Intercept only model (numerator for stabilised weights)
    cens_model_int <- survival::coxph(stats::as.formula(paste("survival::Surv(dtcens_modified, dtcens_s) ~ 1",
                                                              sep = "")),
                                      data = data_raw)
  } else if (is.null(covs)){
    ### If user has not input any predictors for estimating weights, the model for estimating the weights is the intercept only model (i.e. Kaplan Meier estimator)

    ### Intercept only model (numerator for stabilised weights)
    cens_model_int <- survival::coxph(stats::as.formula(paste("survival::Surv(dtcens_modified, dtcens_s) ~ 1",
                                                              sep = "")),
                                      data = data_raw)
    ### Assign cens_model to be the same
    cens_model <- cens_model_int


  }

  ### Calculate a data frame containing probability of censored and uncenosred at each time point
  ### The weights will be the probability of being uncensored, at the time of the event for each individual

  ## Extract baseline hazard
  data_weights <- survival::basehaz(cens_model, centered = FALSE)
  ## Add lp to data_raw_save
  data_raw_save$lp <- stats::predict(cens_model, newdata = data_raw_save, type = "lp", reference = "zero")

  ### Create weights for the cohort at time t - s
  ### Note for individuals who reached an absorbing state, we take the probability of them being uncensored at the time of reached the
  ### abosrbing state. For individuals still alive, we take the probability of being uncensored at time t - s.

  ### Get location of individuals who entered absorbing states or were censored prior to evaluation time
  obs_absorbed_prior <- which(data_raw_save$dtcens <= t & data_raw_save$dtcens_s == 0)
  obs_censored_prior <- which(data_raw_save$dtcens <= t & data_raw_save$dtcens_s == 1)

  ###
  ### Now create unstabilised probability of (un)censoring weights
  ### Note that weights are the probability of being uncensored, so if an individual has low probability of being uncesored,
  ### the inervse of this will be big, weighting them strongly
  ###

  ### First assign all individuals a weight of the probability of being uncensored at time t
  ### This is the linear predictor times the cumulative hazard at time t, and appropriate transformation to get a risk
  data_raw_save$pcw <- as.numeric(exp(-exp(data_raw_save$lp)*data_weights$hazard[max(which(data_weights$time <= t - s))]))

  ## Write a function which will extract the uncensored probability for an individual with linear predictor lp at a given time t
  prob_uncens_func <- function(input){

    ## Assign t and person_id
    t <- input[1]
    lp <- input[2]

    if (t <= 0){
      return(NA)
    } else if (t > 0){
      ## Get hazard at appropriate time
      if (t < min(data_weights$time)){
        bhaz_t <- 0
      } else if (t >= min(data_weights$time)){
        bhaz_t <- data_weights$hazard[max(which(data_weights$time <= t))]
      }

      ## Return risk
      return(exp(-exp(lp)*bhaz_t))
    }
  }

  ### Apply this function to all the times at which individuals have entered an absorbing state prior to censoring
  data_raw_save$pcw[obs_absorbed_prior] <- apply(data_raw_save[obs_absorbed_prior, c("dtcens_modified", "lp")], 1, FUN = prob_uncens_func)

  ### For individuals who were censored prior to t, assign the weight as NA
  data_raw_save$pcw[obs_censored_prior] <- NA

  ### Invert these
  data_raw_save$ipcw <- 1/data_raw_save$pcw

  ###
  ### Stabilise these weights dependent on user-input
  ###
  if (stabilised == TRUE){

    ## Extract baseline hazard
    data_weights_numer <- survival::basehaz(cens_model_int, centered = TRUE)

    ### Assign all individuals a weight of the probability of being uncesored at time t
    data_raw_save$pcw_numer <- as.numeric(exp(-data_weights_numer$hazard[max(which(data_weights_numer$time <= t - s))]))

    ### Create stabilised weight
    data_raw_save$ipcw_stab <- data_raw_save$pcw_numer*data_raw_save$ipcw
  }

  ### Finally cap these at 10 and create output object

  ### Create output object
  if (stabilised == FALSE){
    data_raw_save$ipcw <- pmin(data_raw_save$ipcw, max_weight)
    output_weights <- data.frame("id" = data_raw_save$id, "ipcw" = data_raw_save$ipcw)
  } else if (stabilised == TRUE){
    data_raw_save$ipcw <- pmin(data_raw_save$ipcw, max_weight)
    data_raw_save$ipcw_stab <- pmin(data_raw_save$ipcw_stab, max_weight)
    output_weights <- data.frame("id" = data_raw_save$id, "ipcw" = data_raw_save$ipcw_stab)
  }

  return(output_weights)

}
