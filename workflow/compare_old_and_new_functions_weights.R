### Here we will compare the new functions to old functions to ensure results are the same
rm(list=ls())

devtools::load_all()
data("ebmtcal")
data("msebmtcal")
data("tps0")
data("tps100")

data.raw <- ebmtcal
data.mstate <- msebmtcal
tp.pred.s0 <- tps0 |>
  subset(j == 1) |>
  dplyr::select(paste("pstate", 1:6, sep = ""))
tp.pred.s100 <- tps100 |>
  subset(j == 1) |>
  dplyr::select(paste("pstate", 1:6, sep = ""))
j <- 1
s <- 0
t.eval <- 1826

## Extract data for plot manually
ids.uncens <- ebmtcal |>
  subset(dtcens > t.eval | (dtcens < t.eval & dtcens.s == 0)) |>
  dplyr::pull(id)
tp.pred.plot <- tps0 |>
  dplyr::filter(j == 1 & id %in% ids.uncens) |>
  dplyr::select(any_of(paste("pstate", 1:6, sep = "")))


tp.pred.plot2 <- tps0 |>
  dplyr::filter(j == 1) |>
  dplyr::select(any_of(paste("pstate", 1:6, sep = "")))


## Calculate manual weights
weights.manual <-
  calc_weights(data.mstate = msebmtcal,
               data.raw = ebmtcal,
               t = 1826,
               s = 0,
               landmark.type = "state",
               j = 1,
               max.weight = 10,
               stabilised = FALSE)

## Calculate observed event probabilities using weights.manual
dat.calib.blr.w.manual <-
  calibmsm(data.mstate = msebmtcal,
           data.raw = ebmtcal,
           j=1,
           s=0,
           t = 1826,
           tp.pred = tp.pred.s0,
           calib.type = 'blr',
           curve.type = "rcs",
           rcs.nk = 3,
           weights = weights.manual$ipcw,
           transitions.out = c(1,2,3,4,5,6))



#### Define manual functions:

###
### Redefine calc_weights, but change order of all the input arguments (this shouldn't make a difference)
calc_weights_manual <- function(stabilised = FALSE, max.follow = NULL, data.mstate, covs = NULL, landmark.type = "state", j = NULL, t, s, max.weight = 10, data.raw){

  ### Modify everybody to be censored after time t, if a max.follow has been specified
  if(!is.null(max.follow)){
    if (max.follow == "t"){
      data.raw <- dplyr::mutate(data.raw,
                                dtcens.s = dplyr::case_when(dtcens < t + 2 ~ dtcens.s,
                                                            dtcens >= t + 2 ~ 0),
                                dtcens = dplyr::case_when(dtcens < t + 2 ~ dtcens,
                                                          dtcens >= t + 2 ~ t + 2))
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
    output.weights <- data.frame("id" = data.raw.save$id, "ipcw" = data.raw.save$ipcw)
  } else if (stabilised == TRUE){
    data.raw.save$ipcw <- pmin(data.raw.save$ipcw, max.weight)
    data.raw.save$ipcw.stab <- pmin(data.raw.save$ipcw.stab, max.weight)
    output.weights <- data.frame("id" = data.raw.save$id, "ipcw" = data.raw.save$ipcw.stab)
  }

  return(output.weights)

}

###
### Repeat this process (manual definition of calc_weights), again arguments are in different order, but this time an extra argument is added, which adds 'extra.arg' to every weight.
### This extra arguments is something that could be inputted by user, and want to check it does actually change the answer. It should no longer agree with dat.calb.mlr.
calc_weights_manual_extra <- function(stabilised = FALSE, max.follow = NULL, data.mstate, covs = NULL, landmark.type = "state", j = NULL, t, s, max.weight = 10, data.raw, extra.arg = NULL){

  ### Modify everybody to be censored after time t, if a max.follow has been specified
  if(!is.null(max.follow)){
    if (max.follow == "t"){
      data.raw <- dplyr::mutate(data.raw,
                                dtcens.s = dplyr::case_when(dtcens < t + 2 ~ dtcens.s,
                                                            dtcens >= t + 2 ~ 0),
                                dtcens = dplyr::case_when(dtcens < t + 2 ~ dtcens,
                                                          dtcens >= t + 2 ~ t + 2))
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
    output.weights <- data.frame("id" = data.raw.save$id, "ipcw" = data.raw.save$ipcw)
  } else if (stabilised == TRUE){
    data.raw.save$ipcw <- pmin(data.raw.save$ipcw, max.weight)
    data.raw.save$ipcw.stab <- pmin(data.raw.save$ipcw.stab, max.weight)
    output.weights <- data.frame("id" = data.raw.save$id, "ipcw" = data.raw.save$ipcw.stab)
  }

  ### Add this extra argument to the weights, to check it does something
  output.weights$ipcw <- output.weights$ipcw + extra.arg

  return(output.weights)

}


weights.manual1 <-
  calc_weights_manual(data.mstate = msebmtcal,
               data.raw = ebmtcal,
               t = 1826,
               s = 0,
               landmark.type = "state",
               j = 1,
               max.weight = 10,
               stabilised = FALSE)

weights.manual2 <-
  calc_weights_manual_extra(data.mstate = msebmtcal,
               data.raw = ebmtcal,
               t = 1826,
               s = 0,
               landmark.type = "state",
               j = 1,
               max.weight = 10,
               stabilised = FALSE,
               extra.arg = 0.1)
head(weights.manual2)

testthat::test_that('test', {

  expect_false(any(weights.manual1$ipcw[!is.na(weights.manual1$ipcw)] ==
                                       weights.manual2$ipcw[!is.na(weights.manual2$ipcw)]))
})
head(weights.manual1)
head(weights.manual2)


######################
###
######################
devtools::load_all()

dat.calib.blr <-
  calibmsm(data.mstate = msebmtcal,
           data.raw = ebmtcal,
           j=1,
           s=0,
           t = 1826,
           tp.pred = tp.pred.s0, calib.type = 'blr',
           curve.type = "rcs",
           rcs.nk = 3,
           w.covs = c("year", "agecl", "proph", "match"))


dat.calib.blr.w.function1 <-
  calibmsm(data.mstate = msebmtcal,
           data.raw = ebmtcal,
           j=1,
           s=0,
           t = 1826,
           tp.pred = tp.pred.s0, calib.type = 'blr',
           curve.type = "rcs",
           rcs.nk = 3,
           w.function = calc_weights_manual,
           w.covs = c("year", "agecl", "proph", "match"))


dat.calib.blr.w.function2 <-
  calibmsm(data.mstate = msebmtcal,
           data.raw = ebmtcal,
           j=1,
           s=0,
           t = 1826,
           tp.pred = tp.pred.s0, calib.type = 'blr',
           curve.type = "rcs",
           rcs.nk = 3,
           w.function = calc_weights_manual_extra,
           w.covs = c("year", "agecl", "proph", "match"),
           extra.arg = 10)

testthat::test_that('test1', {

  expect_true(all(dat.calib.blr[["plotdata"]][[1]]$obs == dat.calib.blr.w.function1[["plotdata"]][[1]]$obs))
})

testthat::test_that('test2', {

  expect_false(any(dat.calib.blr[["plotdata"]][[1]]$obs == dat.calib.blr.w.function2[["plotdata"]][[1]]$obs))
})


testthat::test_that('test3', {

  expect_false(any(dat.calib.blr.w.function1[["plotdata"]][[1]]$obs == dat.calib.blr.w.function2[["plotdata"]][[1]]$obs))
})


dat.calib.blr <-
  calibmsm(data.mstate = msebmtcal,
           data.raw = ebmtcal,
           j=1,
           s=0,
           t = 1826,
           tp.pred = tp.pred.s0, calib.type = 'blr',
           curve.type = "rcs",
           rcs.nk = 3,
           w.covs = c("year", "agecl", "proph", "match"))


dat.calib.blr.w.function1 <-
  calibmsm(data.mstate = msebmtcal,
           data.raw = ebmtcal,
           j=1,
           s=0,
           t = 1826,
           tp.pred = tp.pred.s0, calib.type = 'blr',
           curve.type = "rcs",
           rcs.nk = 3,
           w.function = calc_weights_manual,
           w.covs = c("year", "agecl", "proph", "match"))


dat.calib.blr.w.function2 <-
  calibmsm(data.mstate = msebmtcal,
           data.raw = ebmtcal,
           j=1,
           s=0,
           t = 1826,
           tp.pred = tp.pred.s0, calib.type = 'blr',
           curve.type = "rcs",
           rcs.nk = 3,
           w.function = calc_weights_manual_extra,
           w.covs = c("year", "agecl", "proph", "match"),
           extra.arg = 10)

testthat::test_that('test1', {

  expect_true(all(dat.calib.blr[["plotdata"]][[1]]$obs == dat.calib.blr.w.function1[["plotdata"]][[1]]$obs))
})

testthat::test_that('test2', {

  expect_false(any(dat.calib.blr[["plotdata"]][[1]]$obs == dat.calib.blr.w.function2[["plotdata"]][[1]]$obs))
})


testthat::test_that('test3', {

  expect_false(any(dat.calib.blr.w.function1[["plotdata"]][[1]]$obs == dat.calib.blr.w.function2[["plotdata"]][[1]]$obs))
})


head(dat.calib.blr.w.function1[[1]][[1]])
head(dat.calib.blr.w.function2[[1]][[1]])
