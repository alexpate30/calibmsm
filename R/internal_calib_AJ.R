### Internal functions to allow estimation of mean calibration using the Landmark Aalen-Johasen estimator.

#' Estimate data for calibration plots using pseudo-values.
#' @description
#' Called in calibmsm::calib_msm to apply the Aalen-Johansen method.
#'
#' @details
#' Calls heavily on calibmsm::calc_obs_pv_boot to estimate observed transition probabilities.
#' Bootstrapping may be applied depending on user input.
#'
#' @returns A list of datasets for each calibration plot.
#'
#' @noRd
calib_aj <- function(data.mstate,
                     data.raw,
                     j,
                     s,
                     t,
                     pv.group.vars = NULL,
                     pv.n.pctls = NULL,
                     CI = FALSE,
                     CI.type = 'bootstrap',
                     CI.R.boot = NULL,
                     CI.seed = 1,
                     transitions.out = NULL,
                     valid.transitions){

  ### 1) If a confidence interval was requested using bootstrapping, use calib_AJ_boot in conjuction with boot::boot.
  ### This must be done separately for each state.

  ### 2) If a confidence interval was requested using parametric form, call calib_AJ_boot once, specifying indices to be
  ### the sequence 1:nrow(data.raw)
  ### PLACEHOLDER - THIS NEEDS TO BE ADDED XXXXX

  ### 3) If a confidence interval was not requested, call calib_AJ_boot once, specifying indices to be
  ### the sequence 1:nrow(data.raw).

  ### Note that 2) and 3) are the same. This is because the function calib_AJ_boot is dependent on CI and CI.type,
  ### which were defined as input into calib_pAJ. They will there give different output (as they should) when it is run.

  ### If a confidence interval was not requested, run this function once,
  if (CI != FALSE & CI.type == "bootstrap"){

    ### Define alpha for CI's
    alpha <- (1-CI/100)/2

    ### Create object to store plot data
    output.object.mean <- vector("list", length(transitions.out))
    names(output.object.mean) <- paste("state", transitions.out, sep = "")

    ### Put function through bootstrap
    boot.mean <- boot::boot(data.raw,
                            calib_AJ_boot,
                            R = CI.R.boot,
                            data.mstate = data.mstate,
                            j = j,
                            s2 = s,
                            t = t,
                            pv.group.vars = pv.group.vars,
                            pv.n.pctls = pv.n.pctls,
                            transitions.out = transitions.out,
                            valid.transitions = valid.transitions)

    ### Extract confidence bands
    lower <- apply(boot.mean$t, 2, stats::quantile, probs = alpha, na.rm = TRUE)
    upper <- apply(boot.mean$t, 2, stats::quantile, probs = 1-alpha, na.rm = TRUE)

    ### Cycle through states and assign to correct part of output object
    for (state in 1:length(transitions.out)){

      ### Put into output object
      output.object.mean[[state]] <- c("mean" = as.numeric(boot.mean$t0[state]),
                                       "mean.lower" = as.numeric(lower[state]),
                                       "mean.upper" = as.numeric(upper[state]))
    }
  } else {

    ### This will just estimate mean calibration without applying bootstrapping
    output.object.mean <- calib_AJ_boot(data.raw = data.raw,
                                        indices = 1:nrow(data.raw),
                                        data.mstate = data.mstate,
                                        j = j,
                                        s2 = s,
                                        t = t,
                                        pv.group.vars = pv.group.vars,
                                        pv.n.pctls = pv.n.pctls,
                                        transitions.out = transitions.out,
                                        valid.transitions = valid.transitions)

    names(output.object.mean) <- paste("state", transitions.out, sep = "")

  }

  ### Define output object
  output.object.mean <- list("mean" = output.object.mean)

  return(output.object.mean)

}

#' Estimate mean calibration in landmarked cohort using Aalen-Johansen estimator, within subgroups if specified.
#'
#' @description
#' Estimate mean calibration in landmarked cohort using Aalen-Johansen estimator, within subgroups if specified.
#' Function is called by calibmsm::calib_AJ, which is called by calibmsm::calib_msm.
#'
#' @details
#' Function written in a format so that it can be used in combination with \code{\link[boot]{boot}}
#' for bootstrapping. Specifying `indices = 1:nrow(data.raw)` will produce mean calibration for
#' data.raw with no bootstrapping applied.
#'
#' @returns A vector of mean calibration.
#'
#' @noRd
calib_AJ_boot <- function(data.raw,
                          indices,
                          data.mstate,
                          j,
                          s2, # can't use 's' because it matches an argument for the boot function
                          t,
                          pv.group.vars,
                          pv.n.pctls,
                          transitions.out,
                          valid.transitions){

  ### Create object 's' from 's2'
  s <- s2

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
  calc_aj_subgroup <- function(subset.ids, state.k = NULL){

    ### Calculate Aalen-Johansen
    obs.aj <- calc_aj(data.mstate = base::subset(data.mstate.lmk.js, id %in% subset.ids),
                      tmat = tmat.lmk.js,
                      t = t - s,
                      j = j)[["obs.aj"]]

    ### Extract predicted risks in subgroup
    pred <- data.raw.lmk.js[data.raw.lmk.js$id %in% subset.ids, ]

    ### Calculate the difference netween observed from AJ and mean predicted risk
    calib.aj <- obs.aj[paste("pstate", valid.transitions, sep = "")] - colMeans(pred[, paste("tp.pred", valid.transitions, sep = "")])
    if (!is.null(state.k)){
      calib.aj <- calib.aj[paste("pstate", state.k, sep = "")]
    }

    return(calib.aj)

  }

  ###
  ### APPLY calc_aj_subgroup WITHIN SUBGROUPS TO ESTIMATE PSEUDO-VALUES
  ###

  if (is.null(pv.group.vars) & is.null(pv.n.pctls)){

    ###
    ### 1) No grouping
    ###

    ### Calculate psuedo-value for each individual
    aj.out <- calc_aj_subgroup(data.raw.lmk.js$id)
    aj.out <- unlist(aj.out)

    ### Reduce to transitions.out
    aj.out <- aj.out[paste("pstate", transitions.out, sep = "")]

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
    group.ids <- lapply(data.groups, function(x) as.numeric(x[,c("id")]))

    ### Calculate pseudo-values in each subgroup
    aj.out <- lapply(group.ids, calc_aj_subgroup)

    ### Take the average of these and get into correct format
    aj.out <- colMeans(do.call("rbind", aj.out))

    ### Reduce to transitions.out
    aj.out <- aj.out[paste("pstate", transitions.out, sep = "")]

  } else if (is.null(pv.group.vars) & !is.null(pv.n.pctls)) {

    ###
    ### 3) Grouping only by predicted risk
    ###

    ### Write a function to calculate the pseudo-values for the transition probailities into a particular state,
    ### within subgroups defined by percentiles of the predicted transition probabilities into that state.

    ### Note that we now need to apply the function seperately to each state because the subgroups will change depending on the state of interest.
    ### In 1) and 2), we could calculate the pseudo-values for all states simultaneously within each subgroup.
    apply_calc_aj_subgroup_pctls <- function(state.k){

      ### Split data by predicted risk of state k
      data.pctls <- base::split(data.raw.lmk.js,
                                cut(data.raw.lmk.js[,paste("tp.pred", state.k, sep = "")],
                                    breaks =  stats::quantile(data.raw.lmk.js[,paste("tp.pred", state.k, sep = "")],
                                                              seq(0,1,1/pv.n.pctls)),
                                    include.lowest = TRUE))

      ### Get group ids for subgroups
      group.ids <- lapply(data.pctls, function(x) as.numeric(x[,c("id")]))

      ### Calculate pseudo-values in each subgroup
      aj.temp <- lapply(group.ids, calc_aj_subgroup, state.k = state.k)

      return(aj.temp)

    }

    ### Calculate pseudo-values in each subgroup
    aj.out <- lapply(transitions.out, apply_calc_aj_subgroup_pctls)

    ### Take the average of these and get into correct format, and also only retain calibration for the state by which
    ### the data was sorted
    aj.out <- lapply(1:length(transitions.out), function(x) {colMeans(do.call("rbind", aj.out[[x]]))})

    ### Concatenate into a vector
    aj.out <- do.call("c", aj.out)

  } else if (!is.null(pv.group.vars) & !is.null(pv.n.pctls)) {

    ###
    ### 4) Grouping by baseline variables and predicted risk
    ###

    ### Again, we must go seperate for each state
    apply_calc_aj_subgroup_pctls_vars <- function(state.k){

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
      group.ids <- lapply(data.groups.pctls, function(x) as.numeric(x[,c("id")]))

      ### Calculate pseudo-values in each subgroup
      aj.temp <- lapply(group.ids, calc_aj_subgroup, state.k = state.k)

      return(aj.temp)

    }

    ### Calculate pseudo-values in each subgroup
    aj.out <- lapply(transitions.out, apply_calc_aj_subgroup_pctls_vars)

    ### Take the average of these and get into correct format
    aj.out <- lapply(1:length(transitions.out), function(x) {colMeans(do.call("rbind", aj.out[[x]]))})

    ### Concatenate into a vector
    aj.out <- do.call("c", aj.out)

  }

  ##########################################
  ### Return output                      ###
  ##########################################

  return(aj.out)

}
