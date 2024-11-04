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
calib_aj <- function(data_ms,
                     data_raw,
                     j,
                     s,
                     t,
                     pv_group_vars = NULL,
                     pv_n_pctls = NULL,
                     CI = FALSE,
                     CI_type = 'bootstrap',
                     CI_R_boot = NULL,
                     CI_seed = 1,
                     transitions_out = NULL,
                     valid_transitions){

  ### 1) If a confidence interval was requested using bootstrapping, use calib_AJ_boot in conjuction with boot::boot.
  ### This must be done separately for each state.

  ### 2) If a confidence interval was requested using parametric form, call calib_AJ_boot once, specifying indices to be
  ### the sequence 1:nrow(data_raw)
  ### PLACEHOLDER - THIS NEEDS TO BE ADDED XXXXX

  ### 3) If a confidence interval was not requested, call calib_AJ_boot once, specifying indices to be
  ### the sequence 1:nrow(data_raw).

  ### Note that 2) and 3) are the same. This is because the function calib_AJ_boot is dependent on CI and CI_type,
  ### which were defined as input into calib_aj. They will there give different output (as they should) when it is run.

  ### If a confidence interval was not requested, run this function once,
  if (CI != FALSE & CI_type == "bootstrap"){

    ### Define alpha for CI's
    alpha <- (1-CI/100)/2

    ### Create object to store plot data
    output_object_mean <- vector("list", length(transitions_out))
    names(output_object_mean) <- paste("state", transitions_out, sep = "")

    ### Put function through bootstrap
    boot_mean <- boot::boot(data_raw,
                            calib_AJ_boot,
                            R = CI_R_boot,
                            data_ms = data_ms,
                            j = j,
                            s2 = s,
                            t = t,
                            pv_group_vars = pv_group_vars,
                            pv_n_pctls = pv_n_pctls,
                            transitions_out = transitions_out,
                            valid_transitions = valid_transitions)

    ### Extract confidence bands
    lower <- apply(boot_mean$t, 2, stats::quantile, probs = alpha, na.rm = TRUE)
    upper <- apply(boot_mean$t, 2, stats::quantile, probs = 1-alpha, na.rm = TRUE)

    ### Cycle through states and assign to correct part of output object
    for (state in 1:length(transitions_out)){

      ### Put into output object
      output_object_mean[[state]] <- c("mean" = as.numeric(boot_mean$t0[state]),
                                       "mean_lower" = as.numeric(lower[state]),
                                       "mean_upper" = as.numeric(upper[state]))
    }
  } else {

    ### This will just estimate mean calibration without applying bootstrapping
    output_object_mean <- calib_AJ_boot(data_raw = data_raw,
                                        indices = 1:nrow(data_raw),
                                        data_ms = data_ms,
                                        j = j,
                                        s2 = s,
                                        t = t,
                                        pv_group_vars = pv_group_vars,
                                        pv_n_pctls = pv_n_pctls,
                                        transitions_out = transitions_out,
                                        valid_transitions = valid_transitions)

    names(output_object_mean) <- paste("state", transitions_out, sep = "")

  }

  ### Define output object
  output_object_mean <- list("mean" = output_object_mean)

  return(output_object_mean)

}

#' Estimate mean calibration in landmarked cohort using Aalen-Johansen estimator, within subgroups if specified.
#'
#' @description
#' Estimate mean calibration in landmarked cohort using Aalen-Johansen estimator, within subgroups if specified.
#' Function is called by calibmsm::calib_AJ, which is called by calibmsm::calib_msm.
#'
#' @details
#' Function written in a format so that it can be used in combination with \code{\link[boot]{boot}}
#' for bootstrapping. Specifying `indices = 1:nrow(data_raw)` will produce mean calibration for
#' data.raw with no bootstrapping applied.
#'
#' @returns A vector of mean calibration.
#'
#' @noRd
calib_AJ_boot <- function(data_raw,
                          indices,
                          data_ms,
                          j,
                          s2, # can't use 's' because it matches an argument for the boot function
                          t,
                          pv_group_vars,
                          pv_n_pctls,
                          transitions_out,
                          valid_transitions){

  ### Create object 's' from 's2'
  s <- s2

  ###
  ### Apply bootstrapping

  ### Create bootstrapped dataset
  data_raw_boot <- data_raw[indices, ]

  ### Create a new id for these individuals (calc_pv_aj relies on each individual having a unique identifier),
  ### meaning the duplicate values in the bootstrapped datasets will cause problems
  data_raw_boot$id2 <- 1:nrow(data_raw_boot)

  ### Create bootstrapped data_ms (we replicate the choice of patients that was chosen in data_raw)
  data_ms_boot <- apply_bootstrap_msdata(data_ms = data_ms, indices = indices)

  ### Extract transition matrix from original msdata object, as this will have been lost when bootstrapping
  tmat <- attributes(data_ms)$trans

  ### Apply attribute tmat to the bootstrapped data_ms dataset
  attributes(data_ms_boot)$trans <- tmat

  ### Set 'id' to be same as 'id2' in bootstrapped datasets, as the function calc_pv_aj works by removing individual
  ### with the 'id' variable
  data_ms_boot$id <- data_ms_boot$id2
  data_raw_boot$id <- data_raw_boot$id2

  ### Relabel data_ms_boot and data_raw_boot and remove '_boot' datasets
  data_raw <- data_raw_boot
  data_ms <- data_ms_boot
  rm(data_raw_boot, data_ms_boot)

  ###
  ### Apply landmarking

  ### For calib_aj, we need to apply landmarking to both data_raw and data_ms
  ### We model the pseudo-values on the predicted transition probabilities in the bootstrapped data_raw dataset
  ### However the calculation of the pseudo-values must be done in the bootstrapped data_ms dataset

  ### Apply landmarking to data_raw and data_ms
  data_raw_lmk_js <- apply_landmark(data_raw = data_raw,
                                    data_ms = data_ms,
                                    j = j,
                                    s = s,
                                    t = t,
                                    exclude_cens_t = FALSE,
                                    data_return = "data_raw")
  data_ms_lmk_js <- apply_landmark(data_raw = data_raw,
                                       data_ms = data_ms,
                                       j = j,
                                       s = s,
                                       t = t,
                                       exclude_cens_t = FALSE,
                                       data_return = "data_ms")

  ###
  ### Restructure mstate data so that time s = time 0, and relabel transitions to 1, 2,...
  ### This is required in order to estimate Aalen-Johansene estimator and calculate pseudo-values

  ### Reduce transition times by s and remove observations which now occur entirely prior to start up
  data_ms_lmk_js <-
    dplyr::mutate(data_ms_lmk_js,
                  Tstart = pmax(0, Tstart - s),
                  Tstop = pmax(0, Tstop - s),
                  time = Tstop - Tstart) |>
    base::subset(!(Tstart == 0 & Tstop == 0))

  ###
  ### Remove observations for transitions which are not made in the landmarked cohort
  ### Otherwise mstate::msfit will throw out an unneccesary (in this context) warning

  ### Start by identifying which transitions these are
  suppressMessages(zero_transition_table <- data_ms_lmk_js |>
                     dplyr::group_by(from, to) |>
                     dplyr::summarise(Frequency = sum(status)))

  ### Only edit data_ms if some transitions have a frequency of zero
  if (any(zero_transition_table$Frequency == 0)){

    ### Extract the transitions
    zero_transition_from <- zero_transition_table$from[zero_transition_table$Frequency == 0]
    zero_transition_to <- zero_transition_table$to[zero_transition_table$Frequency == 0]

    ### Remove them from dataset
    for (i in 1:length(zero_transition_from)){
      data_ms_lmk_js <- base::subset(data_ms_lmk_js, !(from == zero_transition_from[i] & to == zero_transition_to[i]))
      rm(i)
    }
  }

  ### Fit csh's with no predictors
  strata <- survival::strata
  csh_aj <- survival::coxph(survival::Surv(Tstart, Tstop, status) ~ strata(trans), data_ms_lmk_js)

  ### Extract numeric values for transitions that can occur in the landmarked cohort
  landmark_transitions <- as.numeric(sapply(csh_aj[["xlevels"]]$`strata(trans)`, gsub, pattern = ".*=", replacement =  ""))

  ### Create a mapping from the old transition numbers to new transition numbers which are in sequence
  map_transitions <- data.frame("new" = 1:length(landmark_transitions),
                                "old" = landmark_transitions)

  ### Write a function to apply the mapping
  map_func <- function(x){
    if(!is.na(x)){
      if(!(x %in% landmark_transitions)){
        return(NA)
      } else if (x %in% landmark_transitions)
        return(map_transitions$new[map_transitions$old == x])
    } else if (is.na(x))
      return(NA)
  }

  ### Create new tmat for the new transition numbers
  tmat_lmk_js <- apply(tmat, c(1,2), map_func)

  ### Define max_state (note this be will the same as ncol(tmat))
  max_state <- ncol(tmat_lmk_js)

  ######################################
  ### A) CALCULATE THE PSEUDO VALUES ###
  ######################################

  ### Data must now be split up into groups defined by predictor variables (pv_group_vars) and/or predicted risks (pv_n_pctls)

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
  calc_aj_subgroup <- function(subset_ids, state_k = NULL){

    ### Calculate Aalen-Johansen
    obs_aj <- calc_aj(data_ms = base::subset(data_ms_lmk_js, id %in% subset_ids),
                      tmat = tmat_lmk_js,
                      t = t - s,
                      j = j)[["obs_aj"]]

    ### Extract predicted risks in subgroup
    pred <- data_raw_lmk_js[data_raw_lmk_js$id %in% subset_ids, ]

    ### Calculate the difference between observed from AJ and mean predicted risk
    calib_aj <- obs_aj[paste("pstate", valid_transitions, sep = "")] - colMeans(pred[, paste("tp_pred", valid_transitions, sep = "")])
    if (!is.null(state_k)){
      calib_aj <- calib_aj[paste("pstate", state_k, sep = "")]
    }

    return(calib_aj)

  }

  ###
  ### APPLY calc_aj_subgroup WITHIN SUBGROUPS TO ESTIMATE PSEUDO-VALUES
  ###

  if (is.null(pv_group_vars) & is.null(pv_n_pctls)){

    ###
    ### 1) No grouping
    ###

    ### Calculate psuedo-value for each individual
    aj_out <- calc_aj_subgroup(data_raw_lmk_js$id)
    aj_out <- unlist(aj_out)

    ### Reduce to transitions_out
    aj_out <- aj_out[paste("pstate", transitions_out, sep = "")]

  } else if (!is.null(pv_group_vars) & is.null(pv_n_pctls)) {

    ###
    ### 2) Grouping only by baseline variables
    ###

    ### Split data into groups defined by the variables in pv_group_vars
    ## Create formula to split the dataset by (by pv_group_vars)
    split_formula <- stats::as.formula(paste("~ ", paste(pv_group_vars, collapse = "+"), sep = ""))
    ## Split the dataset into the respective groups
    data_groups <- split(data_raw_lmk_js, split_formula)

    ### Get group ids for subgroups
    group_ids <- lapply(data_groups, function(x) as.numeric(x[,c("id")]))

    ### Calculate pseudo-values in each subgroup
    aj_out <- lapply(group_ids, calc_aj_subgroup)

    ### Take the average of these and get into correct format
    aj_out <- colMeans(do.call("rbind", aj_out))

    ### Reduce to transitions_out
    aj_out <- aj_out[paste("pstate", transitions_out, sep = "")]

  } else if (is.null(pv_group_vars) & !is.null(pv_n_pctls)) {

    ###
    ### 3) Grouping only by predicted risk
    ###

    ### Write a function to calculate the pseudo-values for the transition probailities into a particular state,
    ### within subgroups defined by percentiles of the predicted transition probabilities into that state.

    ### Note that we now need to apply the function seperately to each state because the subgroups will change depending on the state of interest.
    ### In 1) and 2), we could calculate the pseudo-values for all states simultaneously within each subgroup.
    apply_calc_aj_subgroup_pctls <- function(state_k){

      ### Split data by predicted risk of state k
      data_pctls <- base::split(data_raw_lmk_js,
                                cut(data_raw_lmk_js[,paste("tp_pred", state_k, sep = "")],
                                    breaks =  stats::quantile(data_raw_lmk_js[,paste("tp_pred", state_k, sep = "")],
                                                              seq(0,1,1/pv_n_pctls)),
                                    include.lowest = TRUE))

      ### Get group ids for subgroups
      group_ids <- lapply(data_pctls, function(x) as.numeric(x[,c("id")]))

      ### Calculate pseudo-values in each subgroup
      aj_temp <- lapply(group_ids, calc_aj_subgroup, state_k = state_k)

      return(aj_temp)

    }

    ### Calculate pseudo-values in each subgroup
    aj_out <- lapply(transitions_out, apply_calc_aj_subgroup_pctls)

    ### Take the average of these and get into correct format, and also only retain calibration for the state by which
    ### the data was sorted
    aj_out <- lapply(1:length(transitions_out), function(x) {colMeans(do.call("rbind", aj_out[[x]]))})

    ### Concatenate into a vector
    aj_out <- do.call("c", aj_out)

  } else if (!is.null(pv_group_vars) & !is.null(pv_n_pctls)) {

    ###
    ### 4) Grouping by baseline variables and predicted risk
    ###

    ### Again, we must go seperate for each state
    apply_calc_aj_subgroup_pctls_vars <- function(state_k){

      ###
      ### Split data into groups defined by the variables in pv_group_vars, and then predicted risk of transition k

      ###
      ### Start by splitting up data by baseline variables

      ### Create formula to split the dataset by (by pv_group_vars)
      split_formula <- stats::as.formula(paste("~ ", paste(pv_group_vars, collapse = "+"), sep = ""))
      ### Split the dataset into the respective groups
      data_groups <- split(data_raw_lmk_js, split_formula)

      ###
      ### Split each dataset of data_groups into groups defined by percentile of predicted risk for state k

      ### Write a function to do this
      split_group_by_pctl <- function(data_in){
        base::split(data_in,
                    cut(data_in[,paste("tp_pred", state_k, sep = "")],
                        breaks =  stats::quantile(data_in[,paste("tp_pred", state_k, sep = "")],
                                                  seq(0,1,1/pv_n_pctls)),
                        include.lowest = TRUE))
      }

      ### Apply to each group in data_groups
      data_groups_pctls <- lapply(data_groups, split_group_by_pctl)

      ### Create a single list containing each of these datasets
      data_groups_pctls <- unlist(data_groups_pctls, recursive = FALSE)

      ### Get group ids for subgroups
      group_ids <- lapply(data_groups_pctls, function(x) as.numeric(x[,c("id")]))

      ### Calculate pseudo-values in each subgroup
      aj_temp <- lapply(group_ids, calc_aj_subgroup, state_k = state_k)

      return(aj_temp)

    }

    ### Calculate pseudo-values in each subgroup
    aj_out <- lapply(transitions_out, apply_calc_aj_subgroup_pctls_vars)

    ### Take the average of these and get into correct format
    aj_out <- lapply(1:length(transitions_out), function(x) {colMeans(do.call("rbind", aj_out[[x]]))})

    ### Concatenate into a vector
    aj_out <- do.call("c", aj_out)

  }

  ##########################################
  ### Return output                      ###
  ##########################################

  return(aj_out)

}
