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
calib_pv <- function(data_ms,
                     data_raw,
                     tp_pred_plot,
                     j,
                     s,
                     t,
                     curve_type,
                     rcs_nk,
                     loess_span,
                     loess_degree,
                     loess_surface,
                     loess_statistics,
                     loess_trace_hat,
                     loess_cell,
                     loess_iterations,
                     loess_iterTrace,
                     pv_group_vars,
                     pv_n_pctls,
                     pv_precalc,
                     pv_ids,
                     CI,
                     CI_type,
                     CI_R_boot,
                     CI_seed,
                     transitions_out){

  ###
  ### Create the object of predicted risks over which the calibration plots will be plotted

  ### For calib_pv, this is the landmarked cohort of individuals, including those censored at time t,
  ### or tp_pred_plot if specified
  if (is.null(tp_pred_plot)){
    ## Note that tp_pred has been added to data_raw and the predicted transition probabilities (and relevant transformations)
    ## are contained within this dataset.
    data_to_plot <- apply_landmark(data_raw = data_raw, data_ms = data_ms, j = j, s = s, t = t, exclude_cens_t = FALSE)
  } else if (!is.null(tp_pred_plot)){
    data_to_plot <- tp_pred_plot
  }

  ###
  ### Create plotdata object.

  ### 1) If a confidence interval was requested using bootstrapping, use calc_obs_pv_boot in conjuction with boot::boot.
  ### This must be done separately for each state.

  ### 2) If a confidence interval was requested using parametric form, call calc_obs_pv_boot once, specifying indices to be
  ### the sequence 1:nrow(data_raw), in order to calculate the plot data.

  ### 3) If a confidence interval was not requested, call calc_obs_pv_boot once, specifying indices to be
  ### the sequence 1:nrow(data_raw), in order to calculate the plot data.

  ### Note that 2) and 3) are the same. This is because the function calib_pseudo_func is dependent on CI and CI_type,
  ### which were defined as input into calib_pv. They will there give different output (as they should) when it is run.
  if (CI != FALSE & CI_type == "bootstrap" & is.null(pv_ids)){

    ### Define alpha for CI's
    alpha <- (1-CI/100)/2

    ### Create object to store plot data
    plotdata <- vector("list", length(transitions_out))

    ### Cycle through states
    for (state in 1:length(transitions_out)){

      ### Assign state_k
      state_k <- transitions_out[state]

      ### Print progress
      print(paste("Beginning bootstrapping for state = ", state_k, Sys.time()))

      ### Set seed for bootstrapping
      if (!is.null(CI_seed)){
        set.seed(CI_seed)
      }

      ### Put function through bootstrap
      boot_obs <- boot::boot(data_raw,
                             calc_obs_pv_boot,
                             R = CI_R_boot,
                             data_ms = data_ms,
                             data_to_plot = data_to_plot,
                             j = j,
                             s2 = s,
                             t = t,
                             curve_type = curve_type,
                             rcs_nk = rcs_nk,
                             loess_span = loess_span,
                             loess_degree = loess_degree,
                             loess_surface = loess_surface,
                             loess_statistics = loess_statistics,
                             loess_trace_hat = loess_trace_hat,
                             loess_cell = loess_cell,
                             loess_iterations = loess_iterations,
                             loess_iterTrace = loess_iterTrace,
                             pv_group_vars = pv_group_vars,
                             pv_n_pctls = pv_n_pctls,
                             pv_precalc = pv_precalc,
                             pv_ids = pv_ids,
                             CI = FALSE,
                             transitions_out = state_k,
                             boot_format = TRUE)

      ### Extract confidence bands
      lower <- apply(boot_obs$t, 2, stats::quantile, probs = alpha, na.rm = TRUE)
      upper <- apply(boot_obs$t, 2, stats::quantile, probs = 1-alpha, na.rm = TRUE)

      ### Produce a warning if any NA values
      if(sum(is.na(boot_obs$t)) > 0){
        warning(paste("WARNING, SOME BOOTSTRAPPED OBSERVED EVENT PROBABILITIES WERE NA FOR STATE", state_k, "\n",
                      "THERE ARE ", sum(apply(boot_obs$t, 1, function(x) {sum(is.na(x)) > 0})), " ITERATIONS WITH NA's \n",
                      "THE MEAN NUMBER OF NA's IN EACH ITERATION IS", mean(apply(boot_obs$t, 1, function(x) {sum(is.na(x))}))
        ))
      }

      ### Assign output
      if ("id" %in% colnames(data_to_plot)) {
        plotdata[[state]] <- data.frame("id" = data_to_plot$id,
                                        "pred" = data_to_plot[,paste("tp_pred", state_k, sep = "")],
                                        "obs" = boot_obs$t0,
                                        "obs_lower" = lower,
                                        "obs_upper" = upper)
      } else {
        plotdata[[state]] <- data.frame(
          "pred" = data_to_plot[,paste("tp_pred", state_k, sep = "")],
          "obs" = boot_obs$t0,
          "obs_lower" = lower,
          "obs_upper" = upper)
      }
    }
  } else {

    ### NB: calc_obs_pv_boot has the ability to output calibration curve with confidence interval estimated parametrically, as well as outputting
    ### data in boot format (a vector), which was utilised when CI_type = "bootstrap". Here, we specify boot_format = FALSE, which will allow the
    ### confidence interval to be calculated parametrically.

    ### NB: If !is.null(pv_ids), i.e. pv_ids was specified, the function will just return a the pseudo-values themselves.
    plotdata <- calc_obs_pv_boot(data_raw = data_raw,
                                 indices = 1:nrow(data_raw),
                                 data_ms = data_ms,
                                 data_to_plot = data_to_plot,
                                 j = j,
                                 s2 = s,
                                 t = t,
                                 curve_type = curve_type,
                                 rcs_nk = rcs_nk,
                                 loess_span = loess_span,
                                 loess_degree = loess_degree,
                                 loess_surface = loess_surface,
                                 loess_statistics = loess_statistics,
                                 loess_trace_hat = loess_trace_hat,
                                 loess_cell = loess_cell,
                                 loess_iterations = loess_iterations,
                                 loess_iterTrace = loess_iterTrace,
                                 pv_group_vars = pv_group_vars,
                                 pv_n_pctls = pv_n_pctls,
                                 pv_precalc = pv_precalc,
                                 pv_ids = pv_ids,
                                 CI = CI,
                                 CI_type = CI_type,
                                 transitions_out = transitions_out,
                                 boot_format = FALSE)
  }

  ### Define combined output object
  output_object_comb = list("plotdata" = plotdata)

  return(output_object_comb)

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
#' for bootstrapping. Specifying `indices = 1:nrow(data_raw)` will produce calibration
#' curves as normal.
#'
#' @returns If `boot_format = FALSE` a data.frame of predicted and observed event probabilities
#' is returned for each state in `transitions_out`. If boot_format = TRUE, a vector of observed
#' event probabilities is returned. Observed event probabilities are estimated for data points in
#' data_to_plot.
#'
#' @noRd
calc_obs_pv_boot <- function(data_raw,
                             indices,
                             data_ms,
                             data_to_plot,
                             j,
                             s2, # can't use 's' because it matches an argument for the boot function
                             t,
                             curve_type,
                             rcs_nk,
                             loess_span,
                             loess_degree,
                             loess_surface,
                             loess_statistics,
                             loess_trace_hat,
                             loess_cell,
                             loess_iterations,
                             loess_iterTrace,
                             pv_group_vars,
                             pv_n_pctls,
                             pv_precalc,
                             pv_ids,
                             CI,
                             CI_type,
                             transitions_out,
                             boot_format = FALSE){

  ### Create object 's' from 's2'
  s <- s2

  ### If boot_format = TRUE and requested more than one state, stop
  if (boot_format == TRUE & (length(transitions_out) > 1)){
    stop("CANNOT OUTPUT IN BOOT FORMAT IF REQUESTING OBSERVED EVENT PROBABILITIES FOR MORE
         THAN ONE STATE")
  }

  ### The following steps will calculate the pseudo-values before fitting the calibration model
  if (is.null(pv_precalc)){

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

    ### For calib_pv, we need to apply landmarking to both data_raw and data_ms
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
    ### calc_aj: function to calculate to Aalen-Johansen estimator
    ### calc_pv_aj: calculate pseudo-value for an individual based on the Aalen-Johansen estimator

    ### Write one function, which calculates Aalen-Johansen for a group (subset_ids), then calculates pseudo-values for individuals in that group
    ### Can also specify specific individuals to calculate pseudo-values for (pv_ids)
    calc_pv_subgroup <- function(subset_ids, pv_ids = NULL){

      ### Calcuate Aalen-Johansen
      obs_aj <- calc_aj(data_ms = base::subset(data_ms_lmk_js, id %in% subset_ids),
                        tmat = tmat_lmk_js,
                        t = t - s,
                        j = j)[["obs_aj"]]

      ### Now calculate pseudo-values for each individual
      if (is.null(pv_ids)){
        ### Calculate pseudo-values (lapply part of function) and combine into dataset (rbind part of function)
        pv_temp <- do.call("rbind",
                           lapply(subset_ids, calc_pv_aj,
                                  data_ms = base::subset(data_ms_lmk_js, id %in% subset_ids),
                                  obs_aj,
                                  tmat = tmat_lmk_js,
                                  n_cohort = length(subset_ids),
                                  t = t - s,
                                  j = j)
        )

        ### Add id and columns names
        pv_temp <- data.frame(subset_ids, pv_temp)
        colnames(pv_temp) <- c("id", paste("pstate", 1:max_state, sep = ""))

      } else if (!is.null(pv_ids)){
        ### Calculate pseudo-values (lapply part of function) and combine into dataset (rbind part of function)
        pv_temp <- do.call("rbind",
                           lapply(pv_ids, calc_pv_aj,
                                  data_ms = base::subset(data_ms_lmk_js, id %in% subset_ids),
                                  obs_aj,
                                  tmat = tmat_lmk_js,
                                  n_cohort = length(subset_ids),
                                  t = t - s,
                                  j = j)
        )

        ### Add id and columns names
        pv_temp <- data.frame(pv_ids, pv_temp)
        colnames(pv_temp) <- c("id", paste("pstate", 1:max_state, sep = ""))
      }

      return(pv_temp)

    }

    ###
    ### APPLY calc_pv_subgroup WITHIN SUBGROUPS TO ESTIMATE PSEUDO-VALUES
    ###

    if (is.null(pv_group_vars) & is.null(pv_n_pctls)){

      ###
      ### 1) No grouping
      ###

      ### For all individuals
      if (is.null(pv_ids)){

        ### Calculate psuedo-value for each individual
        pv_out <- calc_pv_subgroup(data_raw_lmk_js$id)

        ### For just individual specified in pv_ids
      } else if (!is.null(pv_ids)){

        ### Calculate psuedo-value for each individual in pv_ids
        pv_out <- calc_pv_subgroup(data_raw_lmk_js$id, pv_ids)

      }

    } else if (!is.null(pv_group_vars) & is.null(pv_n_pctls)) {

      ###
      ### 2) Grouping only by baseline variables
      ###

      ### Split data into groups defined by the variables in pv_group_vars
      ## Create formula to split the dataset by (by pv_group_vars)
      split_formula <- stats::as.formula(paste("~ ", paste(pv_group_vars, collapse = "+"), sep = ""))
      ## Split the dataset into the respective groups
      data_groups <- split(data_raw_lmk_js, split_formula)

      ### Get group ids for subgroups in which the pseudo-values need to be calculated
      group_ids <- lapply(data_groups, function(x) as.numeric(x[,c("id")]))

      ### For all individuals
      if (is.null(pv_ids)){

        ### Calculate pseudo-values in each subgroup
        pv_out <- lapply(group_ids, calc_pv_subgroup)

        ### For just individual specified in pv_ids
      } else if (!is.null(pv_ids)){

        ### Identify which pv_ids fit into which subgroup
        group_pv_ids <- lapply(group_ids, function(x) x[x %in% pv_ids])

        ### Calculate pseudo-values in each subgroup, just for pv_ids
        pv_out <- lapply(1:length(group_ids), function(x) {
          if (length(group_pv_ids[[x]] > 0)){
            calc_pv_subgroup(subset_ids = group_ids[[x]], pv_ids = group_pv_ids[[x]])
          } else if (length(group_pv_ids[[x]] == 0)){
            c()
          }
        })

      }

      ### Combine into single dataset
      pv_out <- do.call("rbind", pv_out)

      ### Sort by "id"
      pv_out <- dplyr::arrange(pv_out, id)

    } else if (is.null(pv_group_vars) & !is.null(pv_n_pctls)) {

      ###
      ### 3) Grouping only by predicted risk
      ###

      ### Write a function to calculate the pseudo-values for the transition probailities into a particular state,
      ### within subgroups defined by percentiles of the predicted transition probabilities into that state.

      ### Note that we now need to apply the function seperately to each state because the subgroups will change depending on the state of interest.
      ### In 1) and 2), we could calculate the pseudo-values for all states simultaneously within each subgroup.
      apply_calc_pv_subgroup_pctls <- function(state_k){

        ### Split data by predicted risk of state k
        data_pctls <- base::split(data_raw_lmk_js,
                                  cut(data_raw_lmk_js[,paste("tp_pred", state_k, sep = "")],
                                      breaks =  stats::quantile(data_raw_lmk_js[,paste("tp_pred", state_k, sep = "")],
                                                                seq(0,1,1/pv_n_pctls)),
                                      include.lowest = TRUE))

        ### Get group ids for subgroups
        group_ids <- lapply(data_pctls, function(x) as.numeric(x[,c("id")]))

        ### For all individuals
        if (is.null(pv_ids)){

          ### Calculate pseudo-values in each subgroup
          pv_temp <- lapply(group_ids, calc_pv_subgroup)

          ### For just individual specified in pv_ids
        } else if (!is.null(pv_ids)){

          ### Identify which pv_ids fit into which subgroup
          group_pv_ids <- lapply(group_ids, function(x) x[x %in% pv_ids])

          ### Calculate pseudo-values in each subgroup, just for pv_ids
          pv_temp <- lapply(1:length(group_ids), function(x) {
            if (length(group_pv_ids[[x]] > 0)){
              calc_pv_subgroup(subset_ids = group_ids[[x]], pv_ids = group_pv_ids[[x]])
            } else if (length(group_pv_ids[[x]] == 0)){
              c()
            }
          })

        }

        ### Combine into single dataset
        pv_temp <- do.call("rbind", pv_temp)

        ### Add id, and only retain the pseudo-value for the state of interest (that we sorted the data by)
        ### The pseudo-values for each state are calculated seperately
        pv_temp <- pv_temp[, c("id", paste("pstate", state_k, sep = ""))]

        ### Sort by "id"
        pv_temp <- dplyr::arrange(pv_temp, id)

        return(pv_temp)

      }

      ### Calculate pseudo-values in each subgroup
      pv_out <- lapply(transitions_out, apply_calc_pv_subgroup_pctls)

      ### Combine into a single dataset
      pv_out <- Reduce(function(...) merge(..., by = "id", all.x = TRUE), pv_out)

      ### Arrange
      pv_out <- dplyr::arrange(pv_out, id)

    } else if (!is.null(pv_group_vars) & !is.null(pv_n_pctls)) {

      ###
      ### 4) Grouping by baseline variables and predicted risk
      ###

      ### Again, we must go separate for each state
      apply_calc_pv_subgroup_pctls_vars <- function(state_k){

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

        ### For all individuals
        if (is.null(pv_ids)){

          ### Calculate pseudo-values in each subgroup
          pv_temp <- lapply(group_ids, calc_pv_subgroup)

          ### For just individual specified in pv_ids
        } else if (!is.null(pv_ids)){

          ### Identify which pv_ids fit into which subgroup
          group_pv_ids <- lapply(group_ids, function(x) x[x %in% pv_ids])

          ### Calculate pseudo-values in each subgroup, just for pv_ids
          pv_temp <- lapply(1:length(group_ids), function(x) {
            if (length(group_pv_ids[[x]] > 0)){
              calc_pv_subgroup(subset_ids = group_ids[[x]], pv_ids = group_pv_ids[[x]])
            } else if (length(group_pv_ids[[x]] == 0)){
              c()
            }
          })

        }

        ### Combine into single dataset
        pv_temp <- do.call("rbind", pv_temp)

        ### Add id, and only retain the pseudo-value for the state of interest (that we sorted the data by)
        ### The pseudo-values for each state are calculated seperately
        pv_temp <- pv_temp[, c("id", paste("pstate", state_k, sep = ""))]

        ### Sort by "id"
        pv_temp <- dplyr::arrange(pv_temp, id)

        return(pv_temp)

      }

      ### Calculate pseudo-values in each subgroup
      pv_out <- lapply(transitions_out, apply_calc_pv_subgroup_pctls_vars)

      ### Combine into a single dataset
      pv_out <- Reduce(function(...) merge(..., by = "id", all.x = TRUE), pv_out)

      ### Arrange
      pv_out <- dplyr::arrange(pv_out, id)

    }

    ##########################################
    ### PSEUDO-VALUES HAVE BEEN CALCULATED ###
    ### STORED IN PV_OUT                   ###
    ##########################################

    ### If pv_ids was specified, we just want to return
    ### the pseudo-values and nothing else,
    ### as no point modelling on a subset of the dataset
    if (!is.null(pv_ids)){

      ### Assign column names
      colnames(pv_out)[-1] <- paste("pv_state", transitions_out, sep = "")

      ### Assign output_object
      output_object <- pv_out
    }

    ### If pv_ids was not specified, and we have calculated psudo-values for entire cohort,
    ### generate the observed event probabilities for each state by regressing the calculated pseudo-values
    ### on the predicted transition probabilities
    if (is.null(pv_ids)){

      ###
      ### Create object to store output
      output_object <- vector("list", length(transitions_out))
      names(output_object) <- paste("state", transitions_out, sep = "")

      ###
      ### Loop through and generate observed event probabilities
      for (state in 1:length(transitions_out)){

        ### Assign state_k
        state_k <- transitions_out[state]

        ### Calculate observed event probabilities
        if (curve_type == "loess"){
          obs <- calc_obs_pv_loess_model(pred = data_raw_lmk_js[,paste("tp_pred", state_k, sep = "")],
                                         pv = pv_out[,paste("pstate", state_k, sep = "")],
                                         data_to_plot = data_to_plot[,paste("tp_pred", state_k, sep = "")],
                                         loess_span = loess_span,
                                         loess_degree = loess_degree,
                                         loess_surface = loess_surface,
                                         loess_statistics = loess_statistics,
                                         loess_trace_hat = loess_trace_hat,
                                         loess_cell = loess_cell,
                                         loess_iterations = loess_iterations,
                                         loess_iterTrace = loess_iterTrace,
                                         CI = CI,
                                         CI_type = CI_type)
        } else if (curve_type == "rcs"){
          obs <- calc_obs_pv_rcs_model(pred = data_raw_lmk_js[,paste("tp_pred", state_k, sep = "")],
                                       pv = pv_out[,paste("pstate", state_k, sep = "")],
                                       data_to_plot = data_to_plot[,paste("tp_pred", state_k, sep = "")],
                                       rcs_nk = rcs_nk,
                                       CI = CI,
                                       CI_type = CI_type)
        }

        ### Create output object
        if ("id" %in% colnames(data_to_plot)) {
          output_object[[state]] <- data.frame(
            "id" = data_to_plot$id,
            "pred" = data_to_plot[,paste("tp_pred", state_k, sep = "")],
            obs,
            "pv" = pv_out[,paste("pstate", state_k, sep = "")])

        } else {
          output_object[[state]] <- data.frame(
            "pred" = data_to_plot[,paste("tp_pred", state_k, sep = "")],
            obs,
            "pv" = pv_out[,paste("pstate", state_k, sep = "")])
        }

      }

      ###
      ### If boot_format is true, just return the observed event probabilities for the first state
      if(boot_format == TRUE){
        output_object <- output_object[[1]]$obs
      }

      ### If pseudo-values have been user-inputted, skip the majority of steps and just fit the calibration model using pv_precalc and data_raw
    }
  } else if (!is.null(pv_precalc)){

    ###
    ### Create object to store output
    output_object <- vector("list", length(transitions_out))
    names(output_object) <- paste("state", transitions_out, sep = "")

    ###
    ### Loop through and generate observed event probabilities
    for (state in 1:length(transitions_out)){

      ### Assign state_k
      state_k <- transitions_out[state]

      ### Calculate observed event probabilities
      if (curve_type == "loess"){
        obs <- calc_obs_pv_loess_model(pred = data_raw[,paste("tp_pred", state_k, sep = "")],
                                       pv = pv_precalc[,paste("pstate", state_k, sep = "")],
                                       data_to_plot = data_to_plot[,paste("tp_pred", state_k, sep = "")],
                                       loess_span = loess_span,
                                       loess_degree = loess_degree,
                                       loess_surface = loess_surface,
                                       loess_statistics = loess_statistics,
                                       loess_trace_hat = loess_trace_hat,
                                       loess_cell = loess_cell,
                                       loess_iterations = loess_iterations,
                                       loess_iterTrace = loess_iterTrace,
                                       CI = CI,
                                       CI_type = CI_type)
      } else if (curve_type == "rcs"){
        obs <- calc_obs_pv_rcs_model(pred = data_raw[,paste("tp_pred", state_k, sep = "")],
                                     pv = pv_precalc[,paste("pstate", state_k, sep = "")],
                                     data_to_plot = data_to_plot[,paste("tp_pred", state_k, sep = "")],
                                     rcs_nk = rcs_nk,
                                     CI = CI,
                                     CI_type = CI_type)
      }

      ### Create output object
      if ("id" %in% colnames(data_to_plot)){
        output_object[[state]] <- data.frame(
          "id" = data_to_plot$id,
          "pred" = data_to_plot[,paste("tp_pred", state_k, sep = "")],
          obs,
          "pv" = pv_precalc[,paste("pstate", state_k, sep = "")])

      } else {
        output_object[[state]] <- data.frame(
          "pred" = data_to_plot[,paste("tp_pred", state_k, sep = "")],
          obs,
          "pv" = pv_precalc[,paste("pstate", state_k, sep = "")])
      }
    }
  }

  return(output_object)

}


#' Estimate Aalen-Johansen estimator for a cohort of individuals
#'
#' @description
#' Estimates Aalen-Johansen estimator for the transition probabilities in cohort data_ms.
#' Estimates transition probabilities at time t if in state j at time 0
#' The Aalen-Johansen estimator for the entire cohort (including individual person_id_eval)
#' is inputted manually (obs_aj), to speed up computational time if calculating pseudo-values
#' for multiple individuals from the same cohort.
#'
#' Function is called in calibmsm::calc_obs_pv_boot
#'
#' @param data_ms Validation data in `msdata` format
#' @param tmat Transition probability matrix
#' @param t Follow up time at which calibration is to be assessed
#' @param j Landmark state at which predictions were made
#'
#' @noRd
calc_aj <- function(data_ms, tmat, t, j){

  ### Assign max state number
  max_state <- ncol(tmat)

  ### Fit csh's with no predictors
  strata <- survival::strata
  csh_aj <- survival::coxph(survival::Surv(Tstart, Tstop, status) ~ strata(trans), data_ms)

  ### Calculate cumulative incidence functions using the new transition matrix
  suppressWarnings(
    msfit_aj <- mstate::msfit(csh_aj, trans = tmat)
  )

  ### Calculate Aalen-Johansen estimator
  suppressWarnings(
    pt_aj <- mstate::probtrans(msfit_aj, predt = 0)
  )

  ### Note that warnings are suppressed at both these stages because user will be warned if there are states which can possibly be moved to, but no individual
  ### makes this transition, resulting in zero probabilities. For example in our vignette example, this happens when individuals are in
  ### starting state for 100 days, by definition they can no longer have an adverse event, and mstate gives a warning:
  ### "In max(x[!is.na(x)]) : no non-missing arguments to max; returning -Inf"
  ### There are no problems with this, as it just returns a zero probability of being in that state in the next step (mstate::probtrans), which
  ### A) is correct, and B) we aren't interested in those states anyway

  ### Extract the closest time in the data to the time we want to evaluate at
  t_dat <- pt_aj[[j]]$time[max(which(pt_aj[[j]]$time <= t))]

  ### Extract AJ estimator at this time point
  obs_aj <- pt_aj[[j]][pt_aj[[j]]$time == t_dat, paste("pstate", 1:max_state, sep = "")]

  ### Extract AJ standard error  at this time point
  obs_aj_se <- pt_aj[[j]][pt_aj[[j]]$time == t_dat, paste("se", 1:max_state, sep = "")]

  ### Create output object
  output_object <- list("obs_aj" = obs_aj, "obs_aj_se" = obs_aj_se)

  return(output_object)

}



#' Estimate pseudo-values for the transition probabilities based on the Aalen-Johansen estimator
#'
#' @description
#' Estimates the pseudo-values for an individual (person_id_eval) from cohort data_ms.
#' Calculates psuedo-values for transition probabilities at time t if in state j at time 0
#' The Aalen-Johansen estimator for the entire cohort (including individual person_id_eval)
#' is inputted manually (obs_aj), to speed up computaitonal time if calculating pseudo-values
#' for multiple individuals from the same cohort.
#'
#' Function is called in calibmsm::calc_obs_pv_boot
#'
#' @param person_id_eval id of individual to calculate the pseudo-value for
#' @param data_ms Validation data in `msdata` format
#' @param obs_aj Aalen-Johansen estimator of the transition probabilities in the entire cohort (not excluding person_id_eval)
#' @param tmat Transition probability matrix
#' @param n_cohort Size of cohort (number of unique entries in data_ms)
#' @param t Follow up time at which calibration is to be assessed
#' @param j Landmark state at which predictions were made
#'
#' @noRd
calc_pv_aj <- function(person_id_eval, data_ms, obs_aj, tmat, n_cohort, t, j){

  ### Calculate AJ estimate without patient in dataset
  est_drop_pat <- calc_aj(subset(data_ms, id != person_id_eval),
                          tmat = tmat,
                          t = t,
                          j = j)

  ### Retain just the estimate (not the standard error)
  est_drop_pat <- est_drop_pat[["obs_aj"]]

  ### Calculate the pseudo-value
  pv_pat <- n_cohort*obs_aj - (n_cohort-1)*est_drop_pat

  return(pv_pat)

}



#' Estimate observed event probabilities using pseudo-values and loess smoothers.
#' @description
#' Estimate observed event probabilities for a given input vector of pseudo-values and
#' predicted transition probabilities. This function is called in calibmsm::calc_obs_pv_boot.
#'
#' @returns A vector of observed event probabilities.
#'
#' @noRd
calc_obs_pv_loess_model <- function(pred, pv, data_to_plot,
                                    loess_span,
                                    loess_degree,
                                    loess_surface,
                                    loess_statistics,
                                    loess_trace_hat,
                                    loess_cell,
                                    loess_iterations,
                                    loess_iterTrace,
                                    CI,
                                    CI_type){

  ### Fit model
  if (CI != FALSE){
    if (CI_type == "parametric"){
      loess_model <- stats::loess(pv ~ pred,
                                  span = loess_span,
                                  degree = loess_degree,
                                  control = stats::loess.control(surface = loess_surface,
                                                                 statistics = loess_statistics,
                                                                 trace.hat = loess_trace_hat,
                                                                 cell = loess_cell,
                                                                 iterations = loess_iterations,
                                                                 iterTrace = loess_iterTrace))
    }
  } else {
    ### If not requiring standard errors for parametric confidence, let statistics = "none" for computational efficiency
    loess_model <- stats::loess(pv ~ pred,
                                span = loess_span,
                                degree = loess_degree,
                                control = stats::loess.control(surface = loess_surface,
                                                               statistics = "none",
                                                               trace.hat = loess_trace_hat,
                                                               cell = loess_cell,
                                                               iterations = loess_iterations,
                                                               iterTrace = loess_iterTrace))
  }


  ## Calculate predicted observed probabilities (and confidence intervals if requested using parametric approach)
  ## Note we do not calculate standard errors if confidence interval has been requested using the bootstrap (or if no CI requested)
  if (CI != FALSE){
    if (CI_type == "parametric"){
      ## Need to split up individuals into smaller grouped otherwise predict.loess with SE = TRUE will give an
      ## error, as it will create a matrix that is too large.

      ## Split data into groups of size 10000
      data_to_plot_list <- split(data_to_plot,
                                 rep(1:ceiling(length(data_to_plot)/10000), each = 10000, length.out = length(data_to_plot)))

      ## Define alpha for CIs
      alpha <- (1-CI/100)/2

      ## Predict observed and create data frame
      obs_data <- lapply(1:length(data_to_plot_list),
                         function(x) {
                           ## Predict observed
                           obs <- predict(loess_model, newdata = data_to_plot_list[[x]], se = TRUE)

                           ## Put into dataframe
                           obs_df <- data.frame("obs" = obs$fit,
                                                "obs_lower" = obs$fit - stats::qnorm(1-alpha)*obs$se,
                                                "obs_upper" = obs$fit + stats::qnorm(1-alpha)*obs$se)
                           ## Return
                           return(obs_df)
                         })

      ## Combine into one data frame
      obs_data <- do.call("rbind", obs_data)

    }
  } else {
    ## Predict observed
    obs <- predict(loess_model, newdata = data_to_plot)
    ## Put into dataframe
    obs_data <- data.frame("obs" = obs)
  }

  ### Return obs_data
  return(obs_data)

}


#' Estimate observed event probabilities using pseudo-values and restricted cubic splines.
#' @description
#' Estimate observed event probabilities for a given input vector of pseudo-values and
#' predicted transition probabilities. This function is called in calibmsm::calc_obs_pv_boot.
#'
#' @returns A vector of observed event probabilities.
#'
#' @noRd
calc_obs_pv_rcs_model <- function(pred, pv, data_to_plot, rcs_nk, CI, CI_type){

  ### Start by transforming pred onto logit scale
  pred_logit <- log(pred/(1-pred))

  ### Create spline terms based on predicted risks
  rcs_pred <- Hmisc::rcspline.eval(pred_logit, nk=rcs_nk, inclx=T)
  colnames(rcs_pred) <- paste("rcs_x", 1:ncol(rcs_pred), sep = "")
  knots_pred <- attr(rcs_pred,"knots")

  ### Create spline terms in data_to_plot (using same knot locations derived from the predicted risks)
  ### Note that if data_to_plot == pred, these will be the same

  ### First transform onto logit scale
  data_to_plot <- log(data_to_plot/(1 - data_to_plot))

  ### Create spline terms
  rcs_data_to_plot <- data.frame(Hmisc::rcspline.eval(data_to_plot, knots = knots_pred, inclx=T))
  colnames(rcs_data_to_plot) <- paste("rcs_x", 1:ncol(rcs_data_to_plot), sep = "")

  ### Create dataset in which to fit the model
  data_rcs <- data.frame("pv" = pv, rcs_pred)

  ### Define equation
  eq_LHS <- paste("pv ~ ", sep = "")
  eq_RHS <- paste("rcs_x", 1:ncol(rcs_data_to_plot), sep = "", collapse = "+")
  eq_rcs <- stats::formula(paste(eq_LHS, eq_RHS, sep = ""))

  ## Fit the model using logit link function
  rcs_model <- stats::glm(eq_rcs, data = data_rcs, family = stats::gaussian(link = "logit"), start = rep(0, ncol(rcs_pred) + 1))

  ## Calculate predicted observed probabilities (and confidence intervals if requested using parametric approach)
  ## Note we do not calculate standard errors if confidence interval has been requested using the bootstrap
  if (CI == FALSE){
    ## Predict observed
    obs <- predict(rcs_model, newdata = rcs_data_to_plot, type = "link")
    ## Put into dataframe
    obs_data <- data.frame("obs" = 1/(1+exp(-obs)))
  } else if (CI != FALSE){
    if (CI_type == "bootstrap"){
      ## Predict observed
      obs <- predict(rcs_model, newdata = rcs_data_to_plot, type = "link")
      ## Put into data frame
      obs_data <- data.frame("obs" = 1/(1+exp(-obs)))
    } else if (CI_type == "parametric"){
      ## Predict observed
      obs <- predict(rcs_model, newdata = rcs_data_to_plot, type = "link", se.fit = TRUE)
      ## Define alpha for CIs
      alpha <- (1-CI/100)/2
      ## Put into dataframe
      obs_data <- data.frame("obs" = 1/(1+exp(-obs$fit)),
                             "obs_lower" = 1/(1+exp(-(obs$fit - stats::qnorm(1-alpha)*obs$se.fit))),
                             "obs_upper" = 1/(1+exp(-(obs$fit + stats::qnorm(1-alpha)*obs$se.fit)))
      )
    }

  }

  ### Return obs_data
  return(obs_data)

}


