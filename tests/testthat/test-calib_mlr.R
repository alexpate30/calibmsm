###
### Tests for calibration curves produced using MLR-IPCW (calib_type = 'mlr')
###

test_that("check calib_msm output", {

  ## Reduce to 300 individuals
  # Extract the predicted transition probabilities out of state j = 1 for first 300 individuals
  tp_pred <- tps0 |>
    dplyr::filter(j == 1) |>
    dplyr::filter(id %in% 1:300) |>
    dplyr::select(any_of(paste("pstate", 1:6, sep = "")))

  ### Reduce ebmtcal
  ebmtcal <- ebmtcal |> dplyr::filter(id %in% 1:300)
  msebmtcal <- msebmtcal |> dplyr::filter(id %in% 1:300)

  ## Expect error if generate with CI
  expect_error(calib_msm(data_ms = msebmtcal,
                         data_raw = ebmtcal,
                         j=1,
                         s=0,
                         t = 1826,
                         tp_pred = tp_pred, calib_type = 'mlr',
                         w_covs = c("year", "agecl", "proph", "match"),
                         mlr_ps_int = 2,
                         mlr_degree = 2,
                         CI = 95,
                         CI_R_boot = 5))

  ## Calculate observed event probabilities
  dat_calib_mlr <-
    suppressWarnings(calib_msm(data_ms = msebmtcal,
                               data_raw = ebmtcal,
                               j=1,
                               s=0,
                               t = 1826,
                               tp_pred = tp_pred, calib_type = 'mlr',
                               w_covs = c("year", "agecl", "proph", "match"),
                               mlr_ps_int = 2,
                               mlr_degree = 2))

  expect_type(dat_calib_mlr, "list")
  expect_equal(class(dat_calib_mlr), c("calib_mlr", "calib_msm"))
  expect_length(dat_calib_mlr[["mean"]], 6)
  expect_length(dat_calib_mlr[["plotdata"]], 6)
  expect_length(dat_calib_mlr[["plotdata"]][["state3"]]$id, 265)
  expect_length(dat_calib_mlr[["plotdata"]][["state6"]]$id, 265)
  expect_no_error(summary(dat_calib_mlr))

})

test_that("check calib_msm output with CI = TRUE (for assess_mean)", {

  ## Reduce to 300 individuals
  # Extract the predicted transition probabilities out of state j = 1 for first 300 individuals
  tp_pred <- tps0 |>
    dplyr::filter(j == 1) |>
    dplyr::filter(id %in% 1:300) |>
    dplyr::select(any_of(paste("pstate", 1:6, sep = "")))

  ### Reduce ebmtcal
  ebmtcal <- ebmtcal |> dplyr::filter(id %in% 1:300)
  msebmtcal <- msebmtcal |> dplyr::filter(id %in% 1:300)

  ## Assess mean calibration
  dat_calib_mlr <- suppressWarnings(calib_msm(data_ms = msebmtcal,
                                              data_raw = ebmtcal,
                                              j=1,
                                              s=0,
                                              t = 1826,
                                              tp_pred = tp_pred,
                                              calib_type = 'mlr',
                                              w_covs = c("year", "agecl", "proph", "match"),
                                              mlr_ps_int = 2,
                                              mlr_degree = 2,
                                              CI = 95,
                                              CI_R_boot = 5,
                                              assess_moderate = FALSE,
                                              assess_mean = TRUE))

  expect_type(dat_calib_mlr, "list")
  expect_equal(class(dat_calib_mlr), c("calib_mlr", "calib_msm"))
  expect_length(dat_calib_mlr[["mean"]], 6)
  expect_no_error(summary(dat_calib_mlr))

})

test_that("check calib_msm output, (j = 1, s = 0),
          Manually define function to estimate weights", {

            ## Reduce to 300 individuals
            # Extract the predicted transition probabilities out of state j = 1 for first 300 individuals
            tp_pred <- tps0 |>
              dplyr::filter(j == 1) |>
              dplyr::filter(id %in% 1:300) |>
              dplyr::select(any_of(paste("pstate", 1:6, sep = "")))

            ### Reduce ebmtcal
            ebmtcal <- ebmtcal |> dplyr::filter(id %in% 1:300)
            msebmtcal <- msebmtcal |> dplyr::filter(id %in% 1:300)

            ###
            ### Calculate observed event probabilities
            dat_calib_mlr <-
              suppressWarnings(calib_msm(data_ms = msebmtcal,
                                         data_raw = ebmtcal,
                                         j = 1,
                                         s = 0,
                                         t = 1826,
                                         tp_pred = tp_pred, calib_type = 'mlr',
                                         w_covs = c("year", "agecl", "proph", "match"),
                                         mlr_ps_int = 2,
                                         mlr_degree = 2))

            ###
            ### Calculate manual weights
            weights_manual <-
              calc_weights(data_ms = msebmtcal,
                           data_raw = ebmtcal,
                           covs =  c("year", "agecl", "proph", "match"),
                           t = 1826,
                           s = 0,
                           landmark_type = "state",
                           j = 1,
                           max_weight = 10,
                           stabilised = FALSE)

            ###
            ### Calculate observed event probabilities same function as internal procedure, and check it agrees with dat_calib_mlr
            dat_calib_mlr_w_manual <-
              suppressWarnings(calib_msm(data_ms = msebmtcal,
                                         data_raw = ebmtcal,
                                         j = 1,
                                         s = 0,
                                         t = 1826,
                                         tp_pred = tp_pred, calib_type = 'mlr',
                                         weights = weights_manual$ipcw,
                                         mlr_ps_int = 2,
                                         mlr_degree = 2))

            expect_equal(dat_calib_mlr[["plotdata"]][[1]], dat_calib_mlr_w_manual[["plotdata"]][[1]])

            ###
            ### Calculate observed event probabilities using an incorrect vector of weights, and see if its different from dat_calib_mlr
            dat_calib_mlr_w_manual <-
              suppressWarnings(calib_msm(data_ms = msebmtcal,
                                         data_raw = ebmtcal,
                                         j = 1,
                                         s = 0,
                                         t = 1826,
                                         tp_pred = tp_pred, calib_type = 'mlr',
                                         weights = rep(1,nrow(weights_manual)),
                                         mlr_ps_int = 2,
                                         mlr_degree = 2))

            expect_false(any(dat_calib_mlr[["plotdata"]][[1]]$obs == dat_calib_mlr_w_manual[["plotdata"]][[1]]$obs))

            ###
            ### Calculate observed event probabilities with w_function, where calc_weights_manual = calc_weights (exactly same as internal procedure)
            calc_weights_manual <- calc_weights

            dat_calib_mlr_w_function <-
              suppressWarnings(calib_msm(data_ms = msebmtcal,
                                         data_raw = ebmtcal,
                                         j = 1,
                                         s = 0,
                                         t = 1826,
                                         tp_pred = tp_pred, calib_type = 'mlr',
                                         w_function = calc_weights_manual,
                                         w_covs = c("year", "agecl", "proph", "match"),
                                         mlr_ps_int = 2,
                                         mlr_degree = 2))

            expect_type(dat_calib_mlr_w_function, "list")
            expect_equal(class(dat_calib_mlr_w_function), c("calib_mlr", "calib_msm"))
            expect_length(dat_calib_mlr_w_function[["mean"]], 6)
            expect_length(dat_calib_mlr_w_function[["plotdata"]], 6)
            expect_length(dat_calib_mlr_w_function[["plotdata"]][[1]]$id, 265)
            expect_length(dat_calib_mlr_w_function[["plotdata"]][[4]]$id, 265)

            ## Check answer is same whether w_function used or not
            expect_equal(dat_calib_mlr[["plotdata"]][[1]], dat_calib_mlr_w_function[["plotdata"]][[1]])
            expect_equal(dat_calib_mlr[["plotdata"]][[4]], dat_calib_mlr_w_function[["plotdata"]][[4]])

            ###
            ### Redefine calc_weights, but change order of all the input arguments (this shouldn't make a difference)
            calc_weights_manual <- function(stabilised = FALSE, max_follow = NULL, data_ms, covs = NULL, landmark_type = "state", j = NULL, t, s, max_weight = 10, data_raw){

              ### Modify everybody to be censored after time t, if a max_follow has been specified
              if(!is.null(max_follow)){
                if (max_follow == "t"){
                  data_raw <- dplyr::mutate(data_raw,
                                            dtcens_s = dplyr::case_when(dtcens < t + 2 ~ dtcens_s,
                                                                        dtcens >= t + 2 ~ 0),
                                            dtcens = dplyr::case_when(dtcens < t + 2 ~ dtcens,
                                                                      dtcens >= t + 2 ~ t + 2))
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
                ### If user has not input any predictors for estimating weights, the model for estimating the weights is the intercept only model (i_e_ Kaplan Meier estimator)

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
              ### abosrbing state_ For individuals still alive, we take the probability of being uncensored at time t - s_

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
                output_weights <- data.frame("id" = data_raw_save$id, "ipcw" = data_raw_save$ipcw, "pcw" = data_raw_save$pcw)
              } else if (stabilised == TRUE){
                data_raw_save$ipcw <- pmin(data_raw_save$ipcw, max_weight)
                data_raw_save$ipcw_stab <- pmin(data_raw_save$ipcw_stab, max_weight)
                output_weights <- data.frame("id" = data_raw_save$id, "ipcw" = data_raw_save$ipcw, "ipcw_stab" = data_raw_save$ipcw_stab, "pcw" = data_raw_save$pcw)
              }

              return(output_weights)

            }

            ### Calculate observed event probabilities with new w_function
            dat_calib_mlr_w_function <-
              suppressWarnings(calib_msm(data_ms = msebmtcal,
                                         data_raw = ebmtcal,
                                         j = 1,
                                         s = 0,
                                         t = 1826,
                                         tp_pred = tp_pred, calib_type = 'mlr',
                                         w_function = calc_weights_manual,
                                         w_covs = c("year", "agecl", "proph", "match"),
                                         mlr_ps_int = 2,
                                         mlr_degree = 2))

            expect_type(dat_calib_mlr_w_function, "list")
            expect_equal(class(dat_calib_mlr_w_function), c("calib_mlr", "calib_msm"))
            expect_length(dat_calib_mlr_w_function[["mean"]], 6)
            expect_length(dat_calib_mlr_w_function[["plotdata"]], 6)
            expect_length(dat_calib_mlr_w_function[["plotdata"]][[1]]$id, 265)
            expect_length(dat_calib_mlr_w_function[["plotdata"]][[4]]$id, 265)

            ## Check answer is same whether w_function used or not
            expect_equal(dat_calib_mlr[["plotdata"]][[1]], dat_calib_mlr_w_function[["plotdata"]][[1]])
            expect_equal(dat_calib_mlr[["plotdata"]][[4]], dat_calib_mlr_w_function[["plotdata"]][[4]])

            ###
            ### Repeat this process (manual definition of calc_weights), again arguments are in different order, but this time an extra argument is added, which adds 'extra_arg' to every weight_
            ### This extra arguments is something that could be inputted by user, and want to check it does actually change the answer_ It should no longer agree with dat_calb_mlr_
            calc_weights_manual_extra <- function(stabilised = FALSE, max_follow = NULL, data_ms, covs = NULL, landmark_type = "state", j = NULL, t, s, max_weight = 10, data_raw, extra_arg = NULL){

              ### Modify everybody to be censored after time t, if a max_follow has been specified
              if(!is.null(max_follow)){
                if (max_follow == "t"){
                  data_raw <- dplyr::mutate(data_raw,
                                            dtcens_s = dplyr::case_when(dtcens < t + 2 ~ dtcens_s,
                                                                        dtcens >= t + 2 ~ 0),
                                            dtcens = dplyr::case_when(dtcens < t + 2 ~ dtcens,
                                                                      dtcens >= t + 2 ~ t + 2))
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
                ### If user has not input any predictors for estimating weights, the model for estimating the weights is the intercept only model (i_e_ Kaplan Meier estimator)

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
              ### abosrbing state_ For individuals still alive, we take the probability of being uncensored at time t - s_

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
                output_weights <- data.frame("id" = data_raw_save$id, "ipcw" = data_raw_save$ipcw, "pcw" = data_raw_save$pcw)
              } else if (stabilised == TRUE){
                data_raw_save$ipcw <- pmin(data_raw_save$ipcw, max_weight)
                data_raw_save$ipcw_stab <- pmin(data_raw_save$ipcw_stab, max_weight)
                output_weights <- data.frame("id" = data_raw_save$id, "ipcw" = data_raw_save$ipcw, "ipcw_stab" = data_raw_save$ipcw_stab, "pcw" = data_raw_save$pcw)
              }

              ### Add this extra argument to the weights, to check it does something
              output_weights$ipcw <- output_weights$ipcw + extra_arg

              return(output_weights)

            }

            ### Calculate observed event probabilities with new w_function
            dat_calib_mlr_w_function <-
              suppressWarnings(calib_msm(data_ms = msebmtcal,
                                         data_raw = ebmtcal,
                                         j = 1,
                                         s = 0,
                                         t = 1826,
                                         tp_pred = tp_pred, calib_type = 'mlr',
                                         w_function = calc_weights_manual_extra,
                                         w_covs = c("year", "agecl", "proph", "match"),
                                         mlr_ps_int = 2,
                                         mlr_degree = 2,
                                         extra_arg = 0.1))

            expect_type(dat_calib_mlr_w_function, "list")
            expect_equal(class(dat_calib_mlr_w_function), c("calib_mlr", "calib_msm"))
            expect_length(dat_calib_mlr_w_function[["mean"]], 6)
            expect_length(dat_calib_mlr_w_function[["plotdata"]], 6)
            expect_length(dat_calib_mlr_w_function[["plotdata"]][[1]]$id, 265)
            expect_length(dat_calib_mlr_w_function[["plotdata"]][[4]]$id, 265)

            ## Check answer is same whether w_function used or not
            expect_false(any(dat_calib_mlr[["plotdata"]][[1]]$obs == dat_calib_mlr_w_function[["plotdata"]][[1]]$obs))

          })

