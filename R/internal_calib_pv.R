### Internal functions to allow estimation of calibration curves using the pseudo-value method.

#' Estimate data for calibration plots using pseudo-values.
#' @description
#' Called in calibmsm::calib_msm to apply the pseudo-value method.
#'
#' @details
#' Calls heavily on calibmsm::calc_obs_pv_boot to estimate observed transition probabilities.
#' Bootstrapping may be applied depending on user input.
#'
#' @returns A list of datasets for each calibration plot.
#'
#' @noRd
calib_pv <- function(data.mstate,
                     data.raw,
                     tp.pred.plot,
                     j,
                     s,
                     t,
                     curve.type = "rcs",
                     rcs.nk = 3,
                     loess.span = 0.75,
                     loess.degree = 2,
                     pv.group.vars = NULL,
                     pv.n.pctls = NULL,
                     CI = FALSE,
                     CI.type = 'parametric',
                     CI.R.boot = NULL,
                     CI.seed = 1,
                     transitions.out = NULL){

  ###
  ### Create the object of predicted risks over which the calibration plots will be plotted

  ### For calib_pv, this is the landmarked cohort of individuals, including those censored at time t,
  ### or tp.pred.plot if specified
  if (is.null(tp.pred.plot)){
    ## Note that tp.pred has been added to data.raw and the predicted transition probabilities (and relevant transformations)
    ## are contained within this dataset.
    data.to.plot <- apply_landmark(data.raw = data.raw, data.mstate = data.mstate, j = j, s = s, t = t, exclude.cens.t = FALSE)
  } else if (!is.null(tp.pred.plot)){
    data.to.plot <- tp.pred.plot
  }

  ###
  ### Create plotdata object.

  ### 1) If a confidence interval was requested using bootstrapping, use calc_obs_pv_boot in conjuction with boot::boot.
  ### This must be done separately for each state.

  ### 2) If a confidence interval was requested using parametric form, call calc_obs_pv_boot once, specifying indices to be
  ### the sequence 1:nrow(data.raw), in order to calculate the plot data.

  ### 3) If a confidence interval was not requested, call calc_obs_pv_boot once, specifying indices to be
  ### the sequence 1:nrow(data.raw), in order to calculate the plot data.

  ### Note that 2) and 3) are the same. This is because the function calib_pseudo_func is dependent on CI and CI.type,
  ### which were defined as input into calib_pv. They will there give different output (as they should) when it is run.

  ### If a confidence interval was not requested, run this function once,
  if (CI != FALSE & CI.type == "bootstrap"){

    ### Define alpha for CI's
    alpha <- (1-CI/100)/2

    ### Create object to store plot data
    plotdata <- vector("list", length(transitions.out))

    ### Cycle through states
    for (state in 1:length(transitions.out)){

      ### Assign state.k
      state.k <- transitions.out[state]

      ### Print progress
      print(paste("Beginning bootstrapping for state = ", state.k, Sys.time()))

      ### Set seed for bootstrapping
      set.seed(CI.seed)

      ### Put function through bootstrap
      boot.obs <- boot::boot(data.raw,
                             calc_obs_pv_boot,
                             R = CI.R.boot,
                             data.mstate = data.mstate,
                             data.to.plot = data.to.plot,
                             j = j,
                             s2 = s,
                             t = t,
                             curve.type = curve.type,
                             rcs.nk = rcs.nk,
                             loess.span = loess.span,
                             loess.degree = loess.degree,
                             pv.group.vars = pv.group.vars,
                             pv.n.pctls = pv.n.pctls,
                             CI = FALSE,
                             transitions.out = state.k,
                             boot.format = TRUE)

      ### Extract confidence bands
      lower <- apply(boot.obs$t, 2, stats::quantile, probs = alpha, na.rm = TRUE)
      upper <- apply(boot.obs$t, 2, stats::quantile, probs = 1-alpha, na.rm = TRUE)

      ### Produce a warning if any NA values
      if(sum(is.na(boot.obs$t)) > 0){
        warning(paste("WARNING, SOME BOOTSTRAPPED OBSERVED EVENT PROBABILITIES WERE NA FOR STATE", state.k, "\n",
                      "THERE ARE ", sum(apply(boot.obs$t, 1, function(x) {sum(is.na(x)) > 0})), " ITERATIONS WITH NA's \n",
                      "THE MEAN NUMBER OF NA's IN EACH ITERATION IS", mean(apply(boot.obs$t, 1, function(x) {sum(is.na(x))}))
        ))
      }

      ### Assign output
      if ("id" %in% colnames(data.to.plot)) {
        plotdata[[state]] <- data.frame("id" = data.to.plot$id,
                                        "pred" = data.to.plot[,paste("tp.pred", state.k, sep = "")],
                                        "obs" = boot.obs$t0,
                                        "obs.lower" = lower,
                                        "obs.upper" = upper)
      } else {
        plotdata[[state]] <- data.frame(
          "pred" = data.to.plot[,paste("tp.pred", state.k, sep = "")],
          "obs" = boot.obs$t0,
          "obs.lower" = lower,
          "obs.upper" = upper)
      }
    }
  } else {

    ### Note that calc_obs_pv_boot is dependent on curve.type and CI.type, parameters input to calib_pv
    plotdata <- calc_obs_pv_boot(data.raw = data.raw,
                                 indices = 1:nrow(data.raw),
                                 data.mstate = data.mstate,
                                 data.to.plot = data.to.plot,
                                 j = j,
                                 s2 = s,
                                 t = t,
                                 curve.type = curve.type,
                                 rcs.nk = rcs.nk,
                                 loess.span = loess.span,
                                 loess.degree = loess.degree,
                                 pv.group.vars = pv.group.vars,
                                 pv.n.pctls = pv.n.pctls,
                                 CI = CI,
                                 CI.type = CI.type,
                                 transitions.out = transitions.out,
                                 boot.format = FALSE)
  }

  return(plotdata)

}


#' Calculate pseudo-values and estimate observed event probabilities using pseudo-values.
#' @description
#' Estimate observed event probabilities for all states using pseudo-values.
#' Function is called by calibmsm::calib_pv, which is called by calibmsm::calib_msm.
#' This function does the heavy lifting, and has two major steps. First the pseudo-values
#' are calculated, taking advantage of functions calibmsm::calc_aj and calibmsm::calc_pv_aj.
#' Secondly, the pseudo-values are regressed on the predicted transition
#' probabilities, in order to generate fitted values (the observed event probabilities).
#' This second stage is implemented through the functions calibmsm::calc_obs_pv_rcs_model
#' or calibmsm::calc_obs_pv_loess_model.
#'
#' @details
#' Function written in a format so that it can be used in combination with \code{\link[boot]{boot}}
#' for bootstrapping. Specifying `indices = 1:nrow(data.raw)` will produce calibration
#' curves as normal.
#'
#' @returns If `boot.format = FALSE` a data.frame of predicted and observed event probabilities
#' is returned for each state in `transitions.out`. If boot.format = TRUE, a vector of observed
#' event probabilities is returned. Observed event probabilities are estimated for data points in
#' data.to.plot.
#'
#' @noRd
calc_obs_pv_boot <- function(data.raw,
                             indices,
                             data.mstate,
                             data.to.plot,
                             j,
                             s2, # can't use 's' because it matches an argument for the boot function
                             t,
                             curve.type,
                             rcs.nk,
                             loess.span,
                             loess.degree,
                             pv.group.vars,
                             pv.n.pctls,
                             CI,
                             CI.type,
                             transitions.out,
                             boot.format = FALSE){

  ### Create object 's' from 's2'
  s <- s2

  ### If boot.format = TRUE and requested more than one state, stop
  if (boot.format == TRUE & (length(transitions.out) > 1)){
    stop("CANNOT OUTPUT IN BOOT FORMAT IF REQUESTING OBSERVED EVENT PROBABILITIES FOR MORE
         THAN ONE STATE")
  }

  ### Create bootstrapped dataset
  data.raw.boot <- data.raw[indices, ]

  ### Create a new id for these individuals (calc_pv_aj relies on each individual having a unique identifier),
  ### meaning the duplicate values in the bootstrapped datasets will cause problems
  data.raw.boot$id2 <- 1:nrow(data.raw.boot)

  ### Create bootstrapped data.mstate (we replicate the choice of patients that was chosen in data.raw)
  data.mstate.boot <-
    do.call("rbind",
            lapply(1:nrow(data.raw.boot),
                   function(x) {
                     base::subset(data.mstate, id == data.raw.boot$id[x]) |>
                       dplyr::mutate(id2 = data.raw.boot$id2[x])
                   }
            )
    )

  ###
  ### Apply bootstrapping and landmarking

  ### Extract transition matrix from msdata object
  tmat <- attributes(data.mstate)$trans

  ### Apply attribute tmat to the bootstrapped data.mstate dataset
  attributes(data.mstate.boot)$trans <- tmat

  ### Set 'id' to be same as 'id2' in bootstrapped datasets, as the function calc_pv_aj works by removing individual
  ### with the 'id' variable
  data.mstate.boot$id <- data.mstate.boot$id2
  data.raw.boot$id <- data.raw.boot$id2

  ### Relabel data.mstate.boot and data.raw.boot and remove '.boot' datasets
  data.raw <- data.raw.boot
  data.mstate <- data.mstate.boot
  rm(data.raw.boot, data.mstate.boot)

  ### For calib_pv, we need to apply landmarking to both data.raw and data.mstate
  ### We model the pseudo-values on the predicted transition probabilities in the bootstrapped data.raw dataset
  ### However the calculation of the pseudo-values must be done in the bootstrapped data.mstate dataset

  ### Apply landmarking to data.raw and data.mstate
  data.raw.lmk.js <- apply_landmark(data.raw = data.raw,
                                    data.mstate = data.mstate,
                                    j = j,
                                    s = s,
                                    t = t,
                                    exclude.cens.t = FALSE,
                                    data.return = "data.raw")
  data.mstate.lmk.js <- apply_landmark(data.raw = data.raw,
                                       data.mstate = data.mstate,
                                       j = j,
                                       s = s,
                                       t = t,
                                       exclude.cens.t = FALSE,
                                       data.return = "data.mstate")

  ###
  ### Restructure mstate data so that time s = time 0, and relabel transitions to 1, 2,...
  ### This is required in order to estimate Aalen-Johansene estimator and calculate pseudo-values

  ### Reduce transition times by s and remove observations which now occur entirely prior to start up
  data.mstate.lmk.js <-
    dplyr::mutate(data.mstate.lmk.js,
                  Tstart = pmax(0, Tstart - s),
                  Tstop = pmax(0, Tstop - s),
                  time = Tstop - Tstart) |>
    base::subset(!(Tstart == 0 & Tstop == 0))

  ###
  ### Remove observations for transitions which are not made in the landmarked cohort
  ### Otherwise mstate::msfit will throw out an unneccesary (in this context) warning

  ### Start by identifying which transitions these are
  suppressMessages(zero.transition.table <- data.mstate.lmk.js |>
                     dplyr::group_by(from, to) |>
                     dplyr::summarise(Frequency = sum(status)))

  ### Only edit data.mstate if some transitions have a frequency of zero
  if (any(zero.transition.table$Frequency == 0)){

    ### Extract the transitions
    zero.transition.from <- zero.transition.table$from[zero.transition.table$Frequency == 0]
    zero.transition.to <- zero.transition.table$to[zero.transition.table$Frequency == 0]

    ### Remove them from dataset
    for (i in 1:length(zero.transition.from)){
      data.mstate.lmk.js <- base::subset(data.mstate.lmk.js, !(from == zero.transition.from[i] & to == zero.transition.to[i]))
      rm(i)
    }
  }

  ### Fit csh's with no predictors
  strata <- survival::strata
  csh.aj <- survival::coxph(survival::Surv(Tstart, Tstop, status) ~ strata(trans), data.mstate.lmk.js)

  ### Extract numeric values for transitions that can occur in the landmarked cohort
  landmark.transitions <- as.numeric(sapply(csh.aj[["xlevels"]]$`strata(trans)`, gsub, pattern = ".*=", replacement =  ""))

  ### Create a mapping from the old transition numbers to new transition numbers which are in sequence
  map.transitions <- data.frame("new" = 1:length(landmark.transitions),
                                "old" = landmark.transitions)

  ### Write a function to apply the mapping
  map.func <- function(x){
    if(!is.na(x)){
      if(!(x %in% landmark.transitions)){
        return(NA)
      } else if (x %in% landmark.transitions)
        return(map.transitions$new[map.transitions$old == x])
    } else if (is.na(x))
      return(NA)
  }

  ### Create new tmat for the new transition numbers
  tmat.lmk.js <- apply(tmat, c(1,2), map.func)

  ### Define max.state (note this be will the same as ncol(tmat))
  max.state <- ncol(tmat.lmk.js)

  ######################################
  ### A) CALCULATE THE PSEUDO VALUES ###
  ######################################

  ### Data must now be split up into groups defined by predictor variables (pv.group.vars) and/or predicted risks (pv.n.pctls)

  ### Pseudo-values will be calculated seperately within each of these groups. We will also calculate
  ### the Aalen-Johansen estimate of observed risk within each of these groups to enable quicker
  ### estimation of pseudo-values

  ### To maximise code efficiency, there are some differences depending on whether groups have been defined using predictor
  ### variables or predicted risks.

  ### 1) If no grouping at all, just need to calculate pseudo-values for each individual within the entire group
  ### (don't need to do pseudo-values for each transition seperately, because the grouping is the same)

  ### 2) If grouping is only within variables, again, just need to calculate pseudo-values for each individual within the groups
  ### defined by the variables (don't need to do pseudo-values for each transition seperately, because the grouping is the same)

  ### 3) If grouping is done by predicted risk of each transition (with or without grouping by baseline variables),
  ### need to calculate pseudo-values for each individual seperately for each transition,
  ### as the ordering of individuals, and therefore group, will be different for each transition.

  ### Some references to other functions.
  ### calc_aj: function to calculate to Aalen-Johanser estimator
  ### calc_pv_aj: calculate pseudo-value for an individual based on the Aalen-Johansen estimator


  ### There should reallly just be one function, which calcualtes Aalen-Johansen for a group, then calculates pseudo-values for individuals in that group
  calc_pv_subgroup <- function(subset.ids){

    ### Calcuate Aalen-Johansen
    obs.aj <- calc_aj(data.mstate = base::subset(data.mstate.lmk.js, id %in% subset.ids),
                      tmat = tmat.lmk.js,
                      t = t - s,
                      j = j)[["obs.aj"]]

    ### Now calculate pseudo-values for each individual
    ### Calculate pseudo-values (lapply part of function) and combine into dataset (rbind part of function)
    pv.temp <- do.call("rbind",
                       lapply(subset.ids, calc_pv_aj,
                              data.mstate = base::subset(data.mstate.lmk.js, id %in% subset.ids),
                              obs.aj,
                              tmat = tmat.lmk.js,
                              n.cohort = length(subset.ids),
                              t = t - s,
                              j = j)
    )

    ### Add id and columns names
    pv.temp <- data.frame(subset.ids, pv.temp)
    colnames(pv.temp) <- c("id", paste("pstate", 1:max.state, sep = ""))

    return(pv.temp)

  }


  ###
  ### APPLY calc_pv_subgroup WITHIN SUBGROUPS TO ESTIMATE PSEUDO-VALUES
  ###

  if (is.null(pv.group.vars) & is.null(pv.n.pctls)){

    ###
    ### 1) No grouping
    ###

    ### Calculate psuedo-value for each individual
    pv.out <- calc_pv_subgroup(data.raw.lmk.js$id)

  } else if (!is.null(pv.group.vars) & is.null(pv.n.pctls)) {

    ###
    ### 2) Grouping only by baseline variables
    ###

    ### Split data into groups defined by the variables in pv.group.vars
    ## Create formula to split the dataset by (by pv.group.vars)
    split.formula <- stats::as.formula(paste("~ ", paste(pv.group.vars, collapse = "+"), sep = ""))
    ## Split the dataset into the respective groups
    data.groups <- split(data.raw.lmk.js, split.formula)

    ### Get group ids for subgroups
    group.ids <- sapply(data.groups, function(x) as.numeric(x[,c("id")]))

    ### Calculate pseudo-values in each subgroup
    pv.out <- lapply(group.ids, calc_pv_subgroup)

    ### Combine into single dataset
    pv.out <- do.call("rbind", pv.out)

    ### Sort by "id"
    pv.out <- dplyr::arrange(pv.out, id)

  } else if (is.null(pv.group.vars) & !is.null(pv.n.pctls)) {

    ###
    ### 3) Grouping only by predicted risk
    ###

    ### Write a function to calculate the pseudo-values for the transition probailities into a particular state,
    ### within subgroups defined by percentiles of the predicted transition probabilities into that state.

    ### Note that we now need to apply the function seperately to each state because the subgroups will change depending on the state of interest.
    ### In 1) and 2), we could calculate the pseudo-values for all states simultaneously within each subgroup.
    apply_calc_pv_subgroup_pctls <- function(state.k){

      ### Split data by predicted risk of state k
      data.pctls <- base::split(data.raw.lmk.js,
                                cut(data.raw.lmk.js[,paste("tp.pred", state.k, sep = "")],
                                    breaks =  stats::quantile(data.raw.lmk.js[,paste("tp.pred", state.k, sep = "")],
                                                              seq(0,1,1/pv.n.pctls)),
                                    include.lowest = TRUE))

      ### Get group ids for subgroups
      group.ids <- sapply(data.pctls, function(x) as.numeric(x[,c("id")]))

      ### Calculate pseudo-values in each subgroup
      pv.temp <- lapply(group.ids, calc_pv_subgroup)

      ### Combine into single dataset
      pv.temp <- do.call("rbind", pv.temp)

      ### Add id, and only retain the pseudo-value for the state of interest (that we sorted the data by)
      ### The pseudo-values for each state are calculated seperately
      pv.temp <- pv.temp[, c("id", paste("pstate", state.k, sep = ""))]

      ### Sort by "id"
      pv.temp <- dplyr::arrange(pv.temp, id)

      return(pv.temp)

    }

    ### Calculate pseudo-values in each subgroup
    pv.out <- lapply(transitions.out, apply_calc_pv_subgroup_pctls)

    ### Combine into a single dataset
    pv.out <- Reduce(function(...) merge(..., by = "id", all.x = TRUE), pv.out)

    ### Arrange
    pv.out <- dplyr::arrange(pv.out, id)

  } else if (!is.null(pv.group.vars) & !is.null(pv.n.pctls)) {

    ###
    ### 4) Grouping by baseline variables and predicted risk
    ###

    ### Again, we must go seperate for each state
    apply_calc_pv_subgroup_pctls_vars <- function(state.k){

      ###
      ### Split data into groups defined by the variables in pv.group.vars, and then predicted risk of transition k

      ###
      ### Start by splitting up data by baseline variables

      ### Create formula to split the dataset by (by pv.group.vars)
      split.formula <- stats::as.formula(paste("~ ", paste(pv.group.vars, collapse = "+"), sep = ""))
      ### Split the dataset into the respective groups
      data.groups <- split(data.raw.lmk.js, split.formula)

      ###
      ### Split each dataset of data.groups into groups defined by percentile of predicted risk for state k

      ### Write a function to do this
      split_group_by_pctl <- function(data.in){
        base::split(data.in,
                    cut(data.in[,paste("tp.pred", state.k, sep = "")],
                        breaks =  stats::quantile(data.in[,paste("tp.pred", state.k, sep = "")],
                                                  seq(0,1,1/pv.n.pctls)),
                        include.lowest = TRUE))
      }

      ### Apply to each group in data.groups
      data.groups.pctls <- lapply(data.groups, split_group_by_pctl)

      ### Create a single list containing each of these datasets
      data.groups.pctls <- unlist(data.groups.pctls, recursive = FALSE)

      ### Get group ids for subgroups
      group.ids <- sapply(data.groups.pctls, function(x) as.numeric(x[,c("id")]))

      ### Calculate pseudo-values in each subgroup
      pv.temp <- lapply(group.ids, calc_pv_subgroup)

      ### Combine into single dataset
      pv.temp <- do.call("rbind", pv.temp)

      ### Add id, and only retain the pseudo-value for the state of interest (that we sorted the data by)
      ### The pseudo-values for each state are calculated seperately
      pv.temp <- pv.temp[, c("id", paste("pstate", state.k, sep = ""))]

      ### Sort by "id"
      pv.temp <- dplyr::arrange(pv.temp, id)

      return(pv.temp)

    }

    ### Calculate pseudo-values in each subgroup
    pv.out <- lapply(transitions.out, apply_calc_pv_subgroup_pctls_vars)

    ### Combine into a single dataset
    pv.out <- Reduce(function(...) merge(..., by = "id", all.x = TRUE), pv.out)

    ### Arrange
    pv.out <- dplyr::arrange(pv.out, id)

  }

  ##########################################
  ### PSEUDO-VALUES HAVE BEEN CALCULATED ###
  ### STORED IN PV.OUT                   ###
  ##########################################

  ###
  ### Now generate the observed event probabilities for each state by regressing the calculated pseudo-values
  ### on the predicted transition probabilities

  ###
  ### Create object to store output
  output.object <- vector("list", length(transitions.out))
  names(output.object) <- paste("state", transitions.out, sep = "")

  ###
  ### Loop through and generate observed event probabilities
  for (state in 1:length(transitions.out)){

    ### Assign state.k
    state.k <- transitions.out[state]

    ### Calculate observed event probabilities
    if (curve.type == "loess"){
      obs <- calc_obs_pv_loess_model(pred = data.raw.lmk.js[,paste("tp.pred", state.k, sep = "")],
                                     pv = pv.out[,paste("pstate", state.k, sep = "")],
                                     data.to.plot = data.to.plot[,paste("tp.pred", state.k, sep = "")],
                                     loess.span = loess.span,
                                     loess.degree = loess.degree,
                                     CI = CI,
                                     CI.type = CI.type)
    } else if (curve.type == "rcs"){
      obs <- calc_obs_pv_rcs_model(pred = data.raw.lmk.js[,paste("tp.pred", state.k, sep = "")],
                                   pv = pv.out[,paste("pstate", state.k, sep = "")],
                                   data.to.plot = data.to.plot[,paste("tp.pred", state.k, sep = "")],
                                   rcs.nk = rcs.nk,
                                   CI = CI,
                                   CI.type = CI.type)
    }

    ### Create output object
    if ("id" %in% colnames(data.to.plot)) {
      output.object[[state]] <- data.frame(
        "id" = data.to.plot$id,
        "pred" = data.to.plot[,paste("tp.pred", state.k, sep = "")],
        obs)

    } else {
      output.object[[state]] <- data.frame(
        "pred" = data.to.plot[,paste("tp.pred", state.k, sep = "")],
        obs)
    }

  }

  ###
  ### If boot.format is true, just return the observed event probabilities for the first state
  if(boot.format == TRUE){
    output.object <- output.object[[1]]$obs
  }

  return(output.object)

}


#' Estimate Aalen-Johansen estimator for a cohort of individuals
#'
#' @description
#' Estimates Aalen-Johansen estimator for the transition probabilities in cohort data.mstate.
#' Estimates transition probabilities at time t if in state j at time 0
#' The Aalen-Johansen estimator for the entire cohort (including individual person_id.eval)
#' is inputted manually (obs.aj), to speed up computational time if calculating pseudo-values
#' for multiple individuals from the same cohort.
#'
#' Function is called in calibmsm::calc_obs_pv_boot
#'
#' @param data.mstate Validation data in `msdata` format
#' @param tmat Transition probability matrix
#' @param t Follow up time at which calibration is to be assessed
#' @param j Landmark state at which predictions were made
#'
#' @noRd
calc_aj <- function(data.mstate, tmat, t, j){

  ### Assign max state number
  max.state <- ncol(tmat)

  ### Fit csh's with no predictors
  strata <- survival::strata
  csh.aj <- survival::coxph(survival::Surv(Tstart, Tstop, status) ~ strata(trans), data.mstate)

  ### Calculate cumulative incidence functions using the new transition matrix
  suppressWarnings(
    msfit.aj <- mstate::msfit(csh.aj, trans = tmat)
  )

  ### Calculate Aalen-Johansen estimator
  suppressWarnings(
    pt.aj <- mstate::probtrans(msfit.aj, predt = 0)
  )

  ### Note that warnings are suppressed at both these stages because user will be warned if there are states which can possibly be moved to, but no individual
  ### makes this transition, resulting in zero probabilities. For example in our vignette example, this happens when individuals are in
  ### starting state for 100 days, by definition they can no longer have an adverse event, and mstate gives a warning:
  ### "In max(x[!is.na(x)]) : no non-missing arguments to max; returning -Inf"
  ### There are no problems with this, as it just returns a zero probability of being in that state in the next step (mstate::probtrans), which
  ### A) is correct, and B) we aren't interested in those states anyway

  ### Extract the closest time in the data to the time we want to evaluate at
  t.dat <- pt.aj[[j]]$time[max(which(pt.aj[[j]]$time <= t))]

  ### Extract AJ estimator at this time point
  obs.aj <- pt.aj[[j]][pt.aj[[j]]$time == t.dat, paste("pstate", 1:max.state, sep = "")]

  ### Extract AJ standard error  at this time point
  obs.aj.se <- pt.aj[[j]][pt.aj[[j]]$time == t.dat, paste("se", 1:max.state, sep = "")]

  ### Create output object
  output.object <- list("obs.aj" = obs.aj, "obs.aj.se" = obs.aj.se)

  return(output.object)

}



#' Estimate pseudo-values for the transition probabilities based on the Aalen-Johansen estimator
#'
#' @description
#' Estimates the pseudo-values for an individual (person_id.eval) from cohort data.mstate.
#' Calculates psuedo-values for transition probabilities at time t if in state j at time 0
#' The Aalen-Johansen estimator for the entire cohort (including individual person_id.eval)
#' is inputted manually (obs.aj), to speed up computaitonal time if calculating pseudo-values
#' for multiple individuals from the same cohort.
#'
#' Function is called in calibmsm::calc_obs_pv_boot
#'
#' @param person_id.eval id of individual to calculate the pseudo-value for
#' @param data.mstate Validation data in `msdata` format
#' @param obs.aj Aalen-Johansen estimator of the transition probabilities in the entire cohort (not excluding person_id.eval)
#' @param tmat Transition probability matrix
#' @param n.cohort Size of cohort (number of unique entries in data.mstate)
#' @param t Follow up time at which calibration is to be assessed
#' @param j Landmark state at which predictions were made
#'
#' @noRd
calc_pv_aj <- function(person_id.eval, data.mstate, obs.aj, tmat, n.cohort, t, j){

  ### Calculate AJ estimate without patient in dataset
  est.drop.pat <- calc_aj(subset(data.mstate, id != person_id.eval),
                          tmat = tmat,
                          t = t,
                          j = j)

  ### Retain just the estimate (not the standard error)
  est.drop.pat <- est.drop.pat[["obs.aj"]]

  ### Calculate the pseudo-value
  pv.pat <- n.cohort*obs.aj - (n.cohort-1)*est.drop.pat

  return(pv.pat)

}



#' Estimate observed event probabilities using pseudo-values and loess smoothers.
#' @description
#' Estimate observed event probabilities for a given input vector of pseudo-values and
#' predicted transition probabilities. This function is called in calibmsm::calc_obs_pv_boot.
#'
#' @returns A vector of observed event probabilities.
#'
#' @noRd
calc_obs_pv_loess_model <- function(pred, pv, data.to.plot, loess.span, loess.degree, CI, CI.type){

  ### Fit model
  loess.model <- stats::loess(pv ~ pred,
                              span = loess.span,
                              degree = loess.degree)

  ## Calculate predicted observed probabilities (and confidence intervals if requested using parametric approach)
  ## Note we do not calculate standard errors if confidence interval has been requested using the bootstrap
  if (CI == FALSE){
    ## Predict observed
    obs <- predict(loess.model, newdata = data.to.plot)
    ## Put into dataframe
    obs.data <- data.frame("obs" = obs)
  } else if (CI != FALSE){
    if (CI.type == "bootstrap"){
      ## Predict observed
      obs <- predict(loess.model, newdata = data.to.plot)
      ## Put into dataframe
      obs.data <- data.frame("obs" = obs)
    } else if (CI.type == "parametric"){
      ## Predict observed
      obs <- predict(loess.model, newdata = data.to.plot, se = TRUE)
      ## Define alpha for CIs
      alpha <- (1-CI/100)/2
      ## Put into dataframe
      obs.data <- data.frame("obs" = obs$fit,
                             "obs.lower" = obs$fit - stats::qnorm(1-alpha)*obs$se,
                             "obs.upper" = obs$fit + stats::qnorm(1-alpha)*obs$se)
    }
  }

  ### Return obs.data
  return(obs.data)

}


#' Estimate observed event probabilities using pseudo-values and restricted cubic splines.
#' @description
#' Estimate observed event probabilities for a given input vector of pseudo-values and
#' predicted transition probabilities. This function is called in calibmsm::calc_obs_pv_boot.
#'
#' @returns A vector of observed event probabilities.
#'
#' @noRd
calc_obs_pv_rcs_model <- function(pred, pv, data.to.plot, rcs.nk, CI, CI.type){

  ### Create spline terms based on predicted risks
  rcs.pred <- Hmisc::rcspline.eval(pred, nk=rcs.nk, inclx=T)
  colnames(rcs.pred) <- paste("rcs.x", 1:ncol(rcs.pred), sep = "")
  knots.pred <- attr(rcs.pred,"knots")

  ### Create spline terms in data.to.plot (using same knot locations derived from the predicted risks)
  ### Note that if data.to.plot == pred, these will be the same
  rcs.data.to.plot <- data.frame(Hmisc::rcspline.eval(data.to.plot ,knots = knots.pred, inclx=T))
  colnames(rcs.data.to.plot) <- paste("rcs.x", 1:ncol(rcs.data.to.plot), sep = "")

  ### Create dataset in which to fit the model
  data.rcs <- data.frame("pv" = pv, rcs.pred)

  ### Define equation
  eq.LHS <- paste("pv ~ ", sep = "")
  eq.RHS <- paste("rcs.x", 1:ncol(rcs.data.to.plot), sep = "", collapse = "+")
  eq.rcs <- stats::formula(paste(eq.LHS, eq.RHS, sep = ""))

  ## Fit the model
  rcs.model <- stats::glm(eq.rcs, data = data.rcs, family = stats::gaussian(link = "logit"), start = rep(0, ncol(rcs.pred) + 1))

  ## Calculate predicted observed probabilities (and confidence intervals if requested using parametric approach)
  ## Note we do not calculate standard errors if confidence interval has been requested using the bootstrap
  if (CI == FALSE){
    ## Predict observed
    obs <- predict(rcs.model, newdata = rcs.data.to.plot, type = "link")
    ## Put into dataframe
    obs.data <- data.frame("obs" = 1/(1+exp(-obs)))
  } else if (CI != FALSE){
    if (CI.type == "bootstrap"){
      ## Predict observed
      obs <- predict(rcs.model, newdata = rcs.data.to.plot, type = "link")
      ## Put into dataframe
      obs.data <- data.frame("obs" = 1/(1+exp(-obs)))
    } else if (CI.type == "parametric"){
      ## Predict observed
      obs <- predict(rcs.model, newdata = rcs.data.to.plot, type = "link", se.fit = TRUE)
      ## Define alpha for CIs
      alpha <- (1-CI/100)/2
      ## Put into dataframe
      obs.data <- data.frame("obs" = 1/(1+exp(-obs$fit)),
                             "obs.lower" = 1/(1+exp(-(obs$fit - stats::qnorm(1-alpha)*obs$se.fit))),
                             "obs.upper" = 1/(1+exp(-(obs$fit + stats::qnorm(1-alpha)*obs$se.fit)))
      )
    }

  }

  ### Return obs.data
  return(obs.data)

}


