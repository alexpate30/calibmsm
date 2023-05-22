# ###
# ### Preliminary code for assessing calibration using pseudo-values
# ###
#
# ### Clear workspace
# rm(list=ls())
#
# ### Load calibmsm and required data
# library("calibmsm")
# library("dplyr")
# library("mstate")
# data("ebmtcal")
# data("msebmtcal")
# data("tps0")
# data("tps100")

# ###
# ### Define a function to calculate Aalen-Johansen estimator of transition probabilities
# ### for a landmark group of patients in state j at time s
# ###
# calc_aj_js <- function(data.mstate, tmat, t.eval, j, s){
#
#   # data.mstate <- msebmtcal
#   # tmat <- attributes(msebmtcal)$trans
#   # t.eval <- 1826
#   # j <- 3
#   # s <- 100
#
#   ### Assign max state number
#   max.state <- max(data.mstate$to)
#
#   ### Identify which individuals are in state j at time s
#   ids.j.s <- extract_ids_states(data.mstate = data.mstate, tmat = tmat, j = j, t.eval = s)
#
#   ### Create new dataset containing just these individuals
#   data.mstate.j.s <- base::subset(data.mstate, id %in% ids.j.s)
#
#   ### Reduce transition times by s and remove observations which now occur entirely prior to start up
#   data.mstate.j.s <-
#     dplyr::mutate(data.mstate.j.s,
#                   Tstart = pmax(0, Tstart - s),
#                   Tstop = pmax(0, Tstop - s),
#                   time = Tstop - Tstart) %>%
#     base::subset(!(Tstart == 0 & Tstop == 0))
#
#
#   ###
#   ### Remove observations for transitions where no individuals make that transition
#   ### Otherwise mstate::msfit will throw out an unneccesary (in this context) warning
#   ###
#
#   ### Start by identifying which transitions these are
#   suppressMessages(zero.transition.table <- data.mstate.j.s %>%
#                      dplyr::group_by(from, to) %>%
#                      dplyr::summarise(Frequency = sum(status)))
#
#   ### Only edit the dataset if some transitions have a drequency of zero
#   if (any(zero.transition.table$Frequency == 0)){
#
#     ## Extract the transitions
#     zero.transition.from <- zero.transition.table$from[zero.transition.table$Frequency == 0]
#     zero.transition.to <- zero.transition.table$to[zero.transition.table$Frequency == 0]
#
#     ## Remove them from dataset
#     for (i in 1:length(zero.transition.from)){
#       data.mstate.j.s <- base::subset(data.mstate.j.s, !(from == zero.transition.from[i] & to == zero.transition.to[i]))
#       rm(i)
#     }
#   }
#
#   ### Fit csh's with no predictors
#   strata <- survival::strata
#   csh.aj <- survival::coxph(survival::Surv(Tstart, Tstop, status) ~ strata(trans), data.mstate.j.s)
#
#   ###
#   ### Extract transitions that can occur after landmarking
#   landmark.transitions <- as.numeric(sapply(csh.aj[["xlevels"]]$`strata(trans)`, gsub, pattern = ".*=", replacement =  ""))
#
#   ###
#   ### Create a mapping
#   map.transitions <- data.frame("new" = 1:length(landmark.transitions),
#                                 "old" = landmark.transitions)
#
#   ###
#   ### Write a function for the mapping
#   map.func <- function(x){
#     if(!is.na(x)){
#       if(!(x %in% landmark.transitions)){
#         return(NA)
#       } else if (x %in% landmark.transitions)
#         return(map.transitions$new[map.transitions$old == x])
#     } else if (is.na(x))
#       return(NA)
#   }
#
#   ###
#   ### Create new tmat
#   tmat.new <- apply(tmat, c(1,2), map.func)
#
#   ### Calculate cumulative incidence functions using the new transition matrix
#   msfit.aj <- mstate::msfit(csh.aj, trans = tmat.new)
#
#   ### Calculate Aalen-Johansen estimator
#   pt.aj <- mstate::probtrans(msfit.aj, predt = 0)
#
#   ### Extract the closest time in the data to the time we want to evaluate at
#   t.eval.dat <- pt.aj[[j]]$time[max(which(pt.aj[[j]]$time <= t.eval - s))]
#
#   ### Extract AJ estimator at this time point
#   obs.aj <- pt.aj[[j]][pt.aj[[j]]$time == t.eval.dat, paste("pstate", 1:max.state, sep = "")]
#
#   ### Extract AJ standard error  at this time point
#   obs.aj.se <- pt.aj[[j]][pt.aj[[j]]$time == t.eval.dat, paste("se", 1:max.state, sep = "")]
#
#   ### Create output object
#   output.object <- list("obs.aj" = obs.aj, "obs.aj.se" = obs.aj.se)
#
#   return(output.object)
# }

calc_aj <- function(data.mstate, tmat, t.eval, j){

  ### Assign max state number
  max.state <- ncol(tmat)

  ### Fit csh's with no predictors
  strata <- survival::strata
  csh.aj <- survival::coxph(survival::Surv(Tstart, Tstop, status) ~ strata(trans), data.mstate)

  ### Calculate cumulative incidence functions using the new transition matrix
  msfit.aj <- mstate::msfit(csh.aj, trans = tmat)

  ### Calculate Aalen-Johansen estimator
  pt.aj <- mstate::probtrans(msfit.aj, predt = 0)

  ### Extract the closest time in the data to the time we want to evaluate at
  t.eval.dat <- pt.aj[[j]]$time[max(which(pt.aj[[j]]$time <= t.eval))]

  ### Extract AJ estimator at this time point
  obs.aj <- pt.aj[[j]][pt.aj[[j]]$time == t.eval.dat, paste("pstate", 1:max.state, sep = "")]

  ### Extract AJ standard error  at this time point
  obs.aj.se <- pt.aj[[j]][pt.aj[[j]]$time == t.eval.dat, paste("se", 1:max.state, sep = "")]

  ### Create output object
  output.object <- list("obs.aj" = obs.aj, "obs.aj.se" = obs.aj.se)

  return(output.object)

}



###
### Define a function to calculate pseudo-value for an individual, using the Aalen-Johansen estimator
calc_pv_aj <- function(person_id.eval, data.mstate, obs.aj, tmat, n.cohort, t.eval, j){

  # person_id.eval <- 3
  # data.mstate <- data.mstate.lmk.js
  # tmat <- tmat.lmk.js
  # obs.aj <- obs.aj.save
  # t.eval <- 1826 - s
  # n.cohort <- nrow(data.raw.lmk.js)
  # j <- 3
  # s <- 100

  ### Calculate AJ estimate without patient in dataset
  est.drop.pat <- calc_aj(subset(data.mstate, id != person_id.eval),
                          tmat = tmat,
                          t.eval = t.eval,
                          j = j)

  ### Retain just the estimate (not the standard error)
  est.drop.pat <- est.drop.pat[["obs.aj"]]

  ### Calculate the pseudo-value
  pv.pat <- n.cohort*obs.aj - (n.cohort-1)*est.drop.pat

  return(pv.pat)

}

###
### Define a function to calculate the calibration plot data for a cohort of individuals
calc_calib_pv <- function(data.mstate, data.raw, j, s, t.eval, tp.pred,
                          group.vars = NULL,
                          n.pctls = NULL,
                          loess.span = 0.75, loess.degree = 2,
                          data.pred.plot = NULL, transitions.out = NULL){

  # data.mstate <- msebmtcal
  # data.raw <- ebmtcal
  #
  # indices <- sample(1:nrow(data.raw), nrow(data.raw), replace = TRUE)
  # data.raw.boot <- data.raw[indices, ]
  # data.mstate.boot <-
  #   do.call("rbind", lapply(
  #     data.raw.boot$id, function(x) {base::subset(data.mstate, id == x)})
  #   )
  # attributes(data.mstate.boot)$trans <- attributes(data.mstate)$trans
  #
  # j <- 3
  # s <- 100
  # t.eval <- 1826
  # tp.pred <- tps100 %>% dplyr::filter(j == 3) %>% dplyr::select(any_of(paste("pstate", 1:6, sep = "")))
  # id.lmk <- extract_ids_states(data.mstate = data.mstate,
  #                              tmat = attributes(data.mstate)$trans,
  #                              j = j,
  #                              t.eval = s)
  # data.pred.plot <- tps100 %>%
  #   dplyr::filter(id %in% id.lmk) %>%
  #   dplyr::filter(j == 3) %>%
  #   dplyr::select(any_of(paste("pstate", 1:6, sep = "")))
  # data.pred.plot <- NULL
  # transitions.out <- 3
  # group.vars <- "year"
  # n.pctls <- 2
  # loess.span <- 0.75
  # loess.degree <- 2
  #
  # data.mstate <- data.mstate.boot
  # data.raw <- data.raw.boot
  #
  #   ###
  #   ### apply a bootstrap temporarily
  #   data.raw <- data.raw[sample(1:nrow(data.raw), nrow(data.raw), replace = TRUE),]

  ###
  ### Warnings and errors

  ### Check if transitions.out is only specified for non-zero columns
  if (!is.null(transitions.out)){
    if (sum(c(colSums(tp.pred) == 0)[transitions.out] == TRUE) > 0){
      stop("Calibraiton curves have been requested for transitions into states which have zero probability of occuring.")
    }
  }

  ###
  ### Combine data and apply landmarking

  ### Assign colnames to predicted transition probabilities (and in data.pred.plot)
  colnames(tp.pred) <- paste("tp.pred", 1:ncol(tp.pred), sep = "")
  if (!is.null(data.pred.plot)){
    colnames(data.pred.plot) <- paste("tp.pred", 1:ncol(tp.pred), sep = "")
  }

  ### Extract what states an individual can move into from state j (states with a non-zero predicted risk)
  valid.transitions <- which(colSums(tp.pred) != 0)

  ### Merge data.raw with predicted risks
  data.raw <- data.frame(data.raw, tp.pred[,valid.transitions])

  ### Extract transition matrix from msdata object
  tmat <- attributes(data.mstate)$trans

  ### Identify which individuals are in state j at time s
  ids.state.js <- extract_ids_states(data.mstate = data.mstate, tmat = tmat, j = j, t.eval = s)

  ### Apply landmarking to data.raw and data.mstate
  data.raw.lmk.js <- data.raw %>% base::subset(id %in% ids.state.js)
  data.mstate.lmk.js <- base::subset(data.mstate, id %in% ids.state.js)

  ###
  ### Restructure mstate data so that time s = time 0, and relabel transitions to 1, 2,...

  ### Reduce transition times by s and remove observations which now occur entirely prior to start up
  data.mstate.lmk.js <-
    dplyr::mutate(data.mstate.lmk.js,
                  Tstart = pmax(0, Tstart - s),
                  Tstop = pmax(0, Tstop - s),
                  time = Tstop - Tstart) %>%
    base::subset(!(Tstart == 0 & Tstop == 0))

  ###
  ### Remove observations for transitions where no individuals make that transition
  ### Otherwise mstate::msfit will throw out an unneccesary (in this context) warning
  ### This does happen, for example as no patients in state 1 after 100 days move into state 3, despite this being
  ### a possible transition.

  ### Start by identifying which transitions these are
  suppressMessages(zero.transition.table <- data.mstate.lmk.js %>%
                     dplyr::group_by(from, to) %>%
                     dplyr::summarise(Frequency = sum(status)))

  ### Only edit the dataset if some transitions have a frequency of zero
  if (any(zero.transition.table$Frequency == 0)){

    ## Extract the transitions
    zero.transition.from <- zero.transition.table$from[zero.transition.table$Frequency == 0]
    zero.transition.to <- zero.transition.table$to[zero.transition.table$Frequency == 0]

    ## Remove them from dataset
    for (i in 1:length(zero.transition.from)){
      data.mstate.lmk.js <- base::subset(data.mstate.lmk.js, !(from == zero.transition.from[i] & to == zero.transition.to[i]))
      rm(i)
    }
  }

  ### Fit csh's with no predictors
  strata <- survival::strata
  csh.aj <- survival::coxph(survival::Surv(Tstart, Tstop, status) ~ strata(trans), data.mstate.lmk.js)

  ### Extract transitions that can occur after landmarking
  landmark.transitions <- as.numeric(sapply(csh.aj[["xlevels"]]$`strata(trans)`, gsub, pattern = ".*=", replacement =  ""))

  ### Create a mapping
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

  ### Create new tmat
  tmat.lmk.js <- apply(tmat, c(1,2), map.func)

  ### Define max.state (note this be will the same as ncol(tmat))
  max.state <- ncol(tmat.lmk.js)

  ###################################
  ### CALCULATE THE PSEUDO VALUES ###
  ###################################

  ### Data must now be split up into groups defined by predictor variables and/or predicted risks
  ### Pseudo-values will be calculated seperately within each of these groups. We will also calculate
  ### the Aalen-Johansen estimate of observed risk within each of these groups to enable quicker
  ### estimation of pseudo-values

  ### Note to self when coding:
  ### 1) If no grouping at all, just need to calculate pseudo-values for each individual within the entire group
  ### (don't need to do pseudo-values for each transition seperately, because the grouping is the same)

  ### 2) If grouping is only within variables, just need to calculate pseudo-values for each individual within the groups
  ### defined by the variables (don't need to do pseudo-values for each transition seperately, because the grouping is the same)

  ### 3) If grouping is done by pctl (with or without grouping by baseline variables),
  ### need to calculate pseudo-values for each individual seperately for each transition,
  ### as the ordering of individuals, and therefore group, will be different for each transition.

  if (is.null(group.vars) & is.null(n.pctls)){

    ### 1) No grouping

    if (!is.null(transitions.out)){
      message("Pseudo-values for each state are calculated simultaenously for each individual, because individuals
              are not being grouped by predicted risk. No efficiency gains to only estimate pseudo-values for a subset
              of states, transitions.out has been ignored")
    }

    ### Calculate the observed Aalen-Johansen once to enable quicker calculation for the pseudo-values
    obs.aj.save <- calc_aj(data.mstate = data.mstate.lmk.js,
                           tmat = tmat.lmk.js,
                           t.eval = t.eval - s,
                           j = j)[["obs.aj"]]


    ### Calculate psuedo-value for each individual
    pv.out <- lapply(data.raw.lmk.js$id, calc_pv_aj,
                     data.mstate = data.mstate.lmk.js,
                     obs.aj = obs.aj.save,
                     tmat = tmat.lmk.js,
                     n.cohort = nrow(data.raw.lmk.js),
                     t.eval = t.eval - s,
                     j = j)

    ### Combine into dataset
    pv.out <- data.frame("id" = data.raw.lmk.js$id, do.call("rbind", pv.out))

  } else if (!is.null(group.vars) & is.null(n.pctls)) {

    ### 2) Grouping only by baseline variables

    if (!is.null(transitions.out)){
      message("Pseudo-values for each state are calculated simultaenously for each individual, because individuals
              are not being grouped by predicted risk. No efficiency gains to only estimate pseudo-values for a subset
              of states, transitions.out has been ignored")
    }

    ###
    ### Split data into groups defined by the variables in group.vars

    ### Create formula to split the dataset by (by group.vars)
    split.formula <- as.formula(paste("~ ", paste(group.vars, collapse = "+"), sep = ""))
    ### Split the dataset into the respective groups
    data.groups <- split(data.raw.lmk.js, split.formula)

    ###
    ### Calculate the Aalen-Johansen estimator within each group

    ### Write a function to calculate Aalen-Johansen for patients in a subgrouped dataset,
    ### defined by the baseline variables
    calc_aj_group <- function(group){
      calc_aj(data.mstate = base::subset(data.mstate.lmk.js, id %in% data.groups[[group]]$id),
              tmat = tmat.lmk.js,
              t.eval = t.eval - s,
              j = j)[["obs.aj"]]
    }

    ### Apply this function to each of the datasets in data.groups
    obs.aj.groups <- lapply(1:length(data.groups), calc_aj_group)

    ###
    ### Now calculate pseudo-values within these groups seperately

    ### Create a function which will extract pseudo-values for all individuals in a
    ### subgrouped dataset, defined by the baseline variables
    calc_pv_group <- function(group){

      ### Calculate pseudo-values (lapply part of function) and combine into dataset (rbind part of function)
      pv.temp <- do.call("rbind",
                         lapply(data.groups[[group]]$id, calc_pv_aj,
                                data.mstate = base::subset(data.mstate.lmk.js, id %in% data.groups[[group]]$id),
                                obs.aj = obs.aj.groups[[group]],
                                tmat = tmat.lmk.js,
                                n.cohort = nrow(data.groups[[group]]),
                                t.eval = t.eval - s,
                                j = j)
      )

      ### Add id and columns names
      pv.temp <- data.frame(data.groups[[group]]$id, pv.temp)
      colnames(pv.temp) <- c("id", paste("pstate", 1:max.state, sep = ""))

      return(pv.temp)

    }

    ### Apply this function to each of the datasets in data.groups
    pv.out <- lapply(1:length(data.groups), calc_pv_group)

    ### Combine into single dataset
    pv.out <- do.call("rbind", pv.out)

    ### Sort by "id"
    pv.out <- arrange(pv.out, id)

  } else if (is.null(group.vars) & !is.null(n.pctls)) {

    ### 3) Grouping only by predicted risk

    ### Because we have to do this seperately for each state [in 1) and 2) we
    ### just fitted a single AJ within each group for all states], this process will only loop through
    ### the transitions specified by transitions.out (or valid.transitions if unspecified)

    ### Assign transitions.out if required
    if (is.null(transitions.out)){
      transitions.out <- valid.transitions
    }

    ### Create object to store data grouped by percentiles of predicted risk
    data.pctls <- vector("list", length(transitions.out))

    ### Create object to store Aalen-Johansen estimator of predicted risk for each group
    obs.aj.pctls <- vector("list", length(transitions.out))

    ### Create object to store pseudo-values for each individual for each transition probability
    pv.out.pctls <- vector("list", length(transitions.out))

    ### Loop trough, and for each state k, split into groups defined by predicted risk of state k
    ### Using a for loop as unlikely to be looping through a high number of states, and effeciency gains should not be an issue
    for (state in 1:length(transitions.out)){

      ### Assign the state of interest
      state.k <- as.numeric(transitions.out[state])
      print(paste("state = ", state.k, Sys.time()))

      ### Split data by predicted risk of state k
      data.pctls[[state]] <- base::split(data.raw.lmk.js,
                                         cut(data.raw.lmk.js[,paste("tp.pred", state.k, sep = "")],
                                             breaks =  quantile(data.raw.lmk.js[,paste("tp.pred", state.k, sep = "")],
                                                                seq(0,1,1/n.pctls)),
                                             include.lowest = TRUE))

      ### Write a function to calculate Aalen-Johansen for patients in a subgrouped dataset,
      ### defined by the percentile/group of predicted risk
      calc_aj_pctl <- function(pctl){
        calc_aj(data.mstate = base::subset(data.mstate.lmk.js, id %in% data.pctls[[state]][[pctl]]$id),
                tmat = tmat.lmk.js,
                t.eval = t.eval - s,
                j = j)[["obs.aj"]]
      }

      ### Apply this function to each of the datasets in data.pctls[[state]]
      obs.aj.pctls[[state]] <- lapply(1:n.pctls, calc_aj_pctl)

      ### Create a function which will extract pseudo-values for all individuals in a
      ### subgrouped dataset, defined by the percentile/group of predicted risk
      calc_pv_pctl <- function(pctl){

        ### Calculate pseudo-values (lapply part of function) and combine into dataset (rbind part of function)
        pv.temp <- do.call("rbind",
                           lapply(data.pctls[[state]][[pctl]]$id, calc_pv_aj,
                                  data.mstate = base::subset(data.mstate.lmk.js, id %in% data.pctls[[state]][[pctl]]$id),
                                  obs.aj = obs.aj.pctls[[state]][[pctl]],
                                  tmat = tmat.lmk.js,
                                  n.cohort = nrow(data.pctls[[state]][[pctl]]),
                                  t.eval = t.eval - s,
                                  j = j)
        )

        ### Add id, and only retain the pseudo-value for the state of interest (that we sorted the data by)
        ### The pseudo-values for each state are calculated seperately
        pv.temp <- data.frame(data.pctls[[state]][[pctl]]$id, pv.temp[,state.k])
        colnames(pv.temp) <- c("id", paste("pstate", state.k, sep = ""))

        return(pv.temp)

      }

      ### Calculate the pseudo-values for individuals in each risk group, combine into a single dataset,
      ### and assign to output object pv.out.pctls
      pv.out.temp <- lapply(1:n.pctls, calc_pv_pctl)
      pv.out.pctls[[state]] <- do.call("rbind", pv.out.temp)

      ### END LOOP FOR VARIABLE 'state'

    }

    ### Combine into a single dataset
    pv.out <- Reduce(function(...) merge(..., by = "id", all.x = TRUE), pv.out.pctls)

  } else if (!is.null(group.vars) & !is.null(n.pctls)) {

    ### 4) Grouping by baseline variables and predicted risk

    ### Assign transitions.out if required
    if (is.null(transitions.out)){
      transitions.out <- valid.transitions
    }

    ### Create object to store data grouped by percentiles of predicted risk
    data.groups.pctls <- vector("list", length(transitions.out))

    ### Create object to store Aalen-Johansen estimator of predicted risk for each group
    obs.aj.groups.pctls <- vector("list", length(transitions.out))

    ### Create object to store pseudo-values for each individual for each transition probability
    pv.out.groups.pctls <- vector("list", length(transitions.out))

    ### Loop trough, and for each state k, split into groups defined by predicted risk of state k
    ### Using a for loop as unlikely to be looping through a high number of states, and effeciency gains should not be an issue
    for (state in 1:length(transitions.out)){

      ### Assign the state of interest
      state.k <- as.numeric(transitions.out[state])
      print(paste("state = ", state.k, Sys.time()))

      ###
      ### Split data into groups defined by the variables in group.vars

      ### Create formula to split the dataset by (by group.vars)
      split.formula <- as.formula(paste("~ ", paste(group.vars, collapse = "+"), sep = ""))
      ### Split the dataset into the respective groups
      data.groups.pctls[[state]] <- split(data.raw.lmk.js, split.formula)

      ### Define n.groups
      n.groups <- length(data.groups.pctls[[state]])

      ###
      ### Split each dataset of data.groups.pctls[[state]] into groups defined by percentile of predicted risk for state k

      ### Write a function to do this
      split_group_by_pctl <- function(data.in){
        base::split(data.in,
                    cut(data.in[,paste("tp.pred", state.k, sep = "")],
                        breaks =  quantile(data.in[,paste("tp.pred", state.k, sep = "")],
                                           seq(0,1,1/n.pctls)),
                        include.lowest = TRUE))
      }

      ### Apply to each group in data.groups.pctls[[state]]
      data.groups.pctls[[state]] <- lapply(data.groups.pctls[[state]], split_group_by_pctl)

      ### NB: Each dataset can be accessed by data.groups.pctls[[state]][[group]][[pctl]]

      ###
      ### Calculate Aalen-Johansen estimator within each group

      ### Create storage object
      obs.aj.groups.pctls[[state]] <- vector("list", n.groups)
      pv.out.groups.pctls[[state]] <- vector("list", n.groups)

      ### Iterate over group
      for (group in 1:n.groups){

        ### Create sub-storage object
        obs.aj.groups.pctls[[state]][[group]] <- vector("list", n.pctls)
        pv.out.groups.pctls[[state]][[group]] <- vector("list", n.pctls)

        ### Iterate over pctl
        for (pctl in 1:n.pctls){

          ### Calculate Aalen-Johansen estimator
          obs.aj.groups.pctls[[state]][[group]][[pctl]] <-
            calc_aj(data.mstate = base::subset(data.mstate.lmk.js, id %in% data.groups.pctls[[state]][[group]][[pctl]]$id),
                    tmat = tmat.lmk.js,
                    t.eval = t.eval - s,
                    j = j)[["obs.aj"]]

          ### Calculate pseudo-values
          pv.temp <- do.call("rbind",
                             lapply(data.groups.pctls[[state]][[group]][[pctl]]$id, calc_pv_aj,
                                    data.mstate = base::subset(data.mstate.lmk.js, id %in% data.groups.pctls[[state]][[group]][[pctl]]$id),
                                    obs.aj = obs.aj.groups.pctls[[state]][[group]][[pctl]],
                                    tmat = tmat.lmk.js,
                                    n.cohort = nrow(data.groups.pctls[[state]][[group]][[pctl]]),
                                    t.eval = t.eval - s,
                                    j = j)
          )

          ### Add id, and only retain the pseudo-value for the state of interest (that we sorted the data by)
          ### The pseudo-values for each state are calculated seperately
          pv.out.groups.pctls[[state]][[group]][[pctl]] <-
            data.frame(data.groups.pctls[[state]][[group]][[pctl]]$id, pv.temp[,state.k])
          colnames(pv.out.groups.pctls[[state]][[group]][[pctl]]) <-
            c("id", paste("pstate", state.k, sep = ""))

          ### END LOOP FOR vARIABLE 'pctl'
        }

        ### Combine the psuedo-values in the groups defined by predicted risk
        pv.out.groups.pctls[[state]][[group]] <- do.call("rbind", pv.out.groups.pctls[[state]][[group]])

        ### END LOOP FOR VARIABLE 'group'
      }

      ### Combine the psuedo-values in the groups defined by baseline variables
      pv.out.groups.pctls[[state]] <- do.call("rbind", pv.out.groups.pctls[[state]])

      ### END LOOP FOR VARIABLE 'state'
    }

    ### Combine into a single dataset
    pv.out <- Reduce(function(...) merge(..., by = "id", all.x = TRUE), pv.out.groups.pctls)

  }

  ##############################################
  ### PSEUDO VALUES HAVE NOW BEEN CALCULATED ###
  ##############################################

  ###
  ### Function to calculate observed event probabilities/calibration plot data
  calc_obs_loess_func <- function(pred, pv, data.pred.plot = NULL){

    ### Fit model
    loess.model <- loess(pv ~ pred,
                         span = loess.span,
                         degree = loess.degree)

    ### Created observed event probabilities for each individual
    if (is.null(data.pred.plot)){
      obs <- predict(loess.model, newdata = pred)
    } else {
      obs <- predict(loess.model, newdata = data.pred.plot)
    }

    ### Only return observed
    return(obs)

  }

  ###
  ### Create object to store output
  output.object <- vector("list", length(transitions.out))
  names(output.object) <- paste("state", transitions.out, sep = "")

  ###
  ### Loop through and generate observed event probabilities
  for (k in 1:length(transitions.out)){

    ### Assign state.k
    state.k <- transitions.out[k]

    ### Calculate observed event probabilities
    obs <- calc_obs_loess_func(pred = data.raw.lmk.js[,paste("tp.pred", state.k, sep = "")],
                               pv = pv.out[,paste("pstate", state.k, sep = "")],
                               data.pred.plot = data.pred.plot[,paste("tp.pred", state.k, sep = "")])


    ### Create output object
    if (is.null(data.pred.plot)){
      output.object[[k]] <- data.frame("id" = data.raw.lmk.js[, "id"],
                                       "pred" = data.raw.lmk.js[,paste("tp.pred", state.k, sep = "")],
                                       "obs" = obs)
    } else {
      output.object[[k]] <- data.frame(
        "pred" = data.pred.plot[,paste("tp.pred", state.k, sep = "")],
        "obs" = obs)
    }

  }

  ### Create metadata object
  metadata <- list("valid.transitions" = as.numeric(valid.transitions),
                   "assessed.transitions" = as.numeric(transitions.out),
                   "CI" = FALSE,
                   "j" = j,
                   "s" = s,
                   "t.eval" = t.eval)

  ### Crate a combined output object with metadata, as well as plot data
  output.object.comb <- list("plotdata" = output.object, "metadata" = metadata)

  ### Assign calib_blr class
  attr(output.object.comb, "class") <- "calib_psuedo"

  return(output.object.comb)

}

### TO DO

### Matts suggestion about helping with the bootstrap
### Apply bootstrap for pseudo values?
