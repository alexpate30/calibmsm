### WRAPPER FOR CALC_CALB_BLR


### Start by defining predicted risks
  tp.pred <- tps100 %>%
    dplyr::filter(j == 3) %>%
    dplyr::select(any_of(paste("pstate", 1:6, sep = "")))

    ## Extract ids for individuals uncensored at t.eval
    ids.uncens <- ebmtcal %>%
      subset(dtcens > t.eval | (dtcens < t.eval & dtcens.s == 0)) %>%
      dplyr::pull(id)
    ## Extract the predicted risks out of state 1 for these individuals
    data.pred.plot <- tps100 %>%
      dplyr::filter(j == 3 & id %in% ids.uncens) %>%
      dplyr::select(any_of(paste("pstate", 1:6, sep = "")))


      weights.manual <-
      calc_weights(data.mstate = msebmtcal,
                   data.raw = ebmtcal,
                   covs = c("year", "agecl", "proph", "match"),
                   t.eval = t.eval,
                   s = 0,
                   landmark.type = "state",
                   j = 1,
                   max.weight = 10,
                   stabilised = FALSE)$ipcw


      dat.calib.boot.manual <-
      calc_calib_blr(data.mstate = msebmtcal,
                     data.raw = ebmtcal,
                     j=1,
                     s=0,
                     t.eval = t.eval,
                     tp.pred = tp.pred,
                     curve.type = "rcs",
                     rcs.nk = 3,
                     weights = weights.manual)

user_defined_function_for_weights <- function(data.in){

  ### A model where we adjust for predictor variables
  cens.model <- survival::coxph(survival::Surv(dtcens.modified, dtcens.s) ~ agecl + year + match + proph, data = data.in)

  ## Extract baseline hazard
  data.weights <- survival::basehaz(cens.model, centered = FALSE)

  ## Make predictions of the linear predictor from this model
  data.raw.save$lp <- stats::predict(cens.model, newdata = data.raw.save, type = "lp", reference = "zero")

  ### Get location of individuals who entered absorbing states or were censored prior to evaluation time
  obs.absorbed.prior <- which(data.raw.save$dtcens <= t.eval & data.raw.save$dtcens.s == 0)
  obs.censored.prior <- which(data.raw.save$dtcens <= t.eval & data.raw.save$dtcens.s == 1)

  ### First assign all individuals a weight of the probability of being uncensored at time t.eval
  ### This is the linear predictor times the cumulative hazard at time t.eval, and appropriate transformation to get a risk
  data.raw.save$pcw <- as.numeric(exp(-exp(data.raw.save$lp)*data.weights$hazard[max(which(data.weights$time <= t.eval))]))

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

  ### For individuals who were censored prior to t.eval, assign the weight as NA
  data.raw.save$pcw[obs.censored.prior] <- NA

  ### Invert these
  data.raw.save$ipcw <- 1/data.raw.save$pcw

}


### I actually think in my bootstrapping procedure, I want to do:

weights <- user_defined_function(data.raw = data.boot,
                                 ...)

my_wrapper <- function(func_in, data.raw, data.mstate, tp.pred){

  ### Function to calculate confidence interval using bootstrapping
  calc_obs_boot <- function(data, indices, tp.pred, state.k){

    ## Bootstrap dataset and predicted transition probabilities
    data.boot <- data[indices,]
    tp.pred.boot <- tp.pred[indices, ]

    ## Calculate weights
    ## In practice - replace this function with your own
    weights.manual <-
      func_in(data.mstate = msebmtcal,
                   data.raw = data.boot,
                   covs = c("year", "agecl", "proph", "match"),
                   t.eval = t.eval,
                   s = 0,
                   landmark.type = "state",
                   j = 1,
                   max.weight = 10,
                   stabilised = FALSE)$ipcw

    ## Estimate bootstrapped calibration curve
    curve.est <-
      calc_calib_blr(data.mstate = msebmtcal,
                     data.raw = data.boot,
                     j=1,
                     s=0,
                     t.eval = t.eval,
                     tp.pred = tp.pred.boot,
                     curve.type = "rcs",
                     rcs.nk = 3,
                     weights = weights.manual,
                     data.pred.plot = data.pred.plot,
                     transitions.out = state.k)

    ## Extract observed event probabilities
    curve.obs <-
      curve.est[["plotdata"]][[paste("state", state.k, sep = "")]]$obs

    return(curve.obs)

  }

}




      alpha <- (1-95/100)/2
    valid.transitions <- which(colSums(tp.pred) != 0)
    plot.data.list <- vector("list", length(valid.transitions))

      for (k in 1:length(valid.transitions)){

        ## Assign state k
        state.k <- valid.transitions[k]

        ## Run bootstrapping
        boot.obs <- boot::boot(ebmtcal,
                               user_defined_function_for_weights,
                               R = 200,
                               tp.pred = tp.pred,
                               state.k = state.k)$t

        ## Extract confidence bands
        lower <- apply(boot.obs, 2, stats::quantile, probs = alpha, na.rm = TRUE)
        upper <- apply(boot.obs, 2, stats::quantile, probs = 1-alpha, na.rm = TRUE)

        ## Assign output
        plot.data.list[[k]] <- data.frame(
          "pred" = dat.calib.boot.manual[["plotdata"]][[k]]$pred,
          "obs" = dat.calib.boot.manual[["plotdata"]][[k]]$obs,
          "obs.lower" = lower,
          "obs.upper" = upper)

      }

      metadata <- list("valid.transitions"= valid.transitions,
                       "CI" = 95,
                       "curve.type" = "rcs")
    dat.calib.blr.manual <- list("plotdata" = plot.data.list, "metadata" = metadata)
    attr(dat.calib.blr.manual, "class") <- "calib_blr"

