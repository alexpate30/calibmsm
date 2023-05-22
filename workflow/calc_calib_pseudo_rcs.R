
load_all()

calc_calib_pv <- function(data.mstate,
                          data.raw,
                          j,
                          s,
                          t.eval,
                          tp.pred,
                          group.vars = NULL,
                          n.pctls = NULL,
                          curve.type = "rcs",
                          loess.span = 0.75,
                          loess.degree = 2,
                          rcs.nk = 3,
                          CI = FALSE,
                          CI.type = NULL,
                          CI.R.boot = NULL,
                          data.pred.plot = NULL,
                          transitions.out = NULL){

  data.mstate <- msebmtcal
  data.raw <- ebmtcal
  #
  # indices <- sample(1:nrow(data.raw), nrow(data.raw), replace = TRUE)
  # data.raw.boot <- data.raw[indices, ]
  # data.mstate.boot <-
  #   do.call("rbind", lapply(
  #     data.raw.boot$id, function(x) {base::subset(data.mstate, id == x)})
  #   )
  # attributes(data.mstate.boot)$trans <- attributes(data.mstate)$trans
  #
  j <- 3
  s <- 100
  t.eval <- 1826
  tp.pred <- tps100 %>% dplyr::filter(j == 3) %>% dplyr::select(any_of(paste("pstate", 1:6, sep = "")))
  # id.lmk <- extract_ids_states(data.mstate = data.mstate,
  #                              tmat = attributes(data.mstate)$trans,
  #                              j = j,
  #                              t.eval = s)
  # data.pred.plot <- tps100 %>%
  #   dplyr::filter(id %in% id.lmk) %>%
  #   dplyr::filter(j == 3) %>%
  #   dplyr::select(any_of(paste("pstate", 1:6, sep = "")))
  data.pred.plot <- NULL
  transitions.out <- NULL
  group.vars <- "year"
  n.pctls <- 2
  loess.span <- 0.75
  loess.degree <- 2
  CI <- FALSE
  CI.type <- parametric
  CI.R.boot <- 2
  rcs.nk <- 3
  curve.type <- "rcs"

  #
  # data.mstate <- data.mstate.boot
  # data.raw <- data.raw.boot
  #
  #   ###
  #   ### apply a bootstrap temporarily
  #   data.raw <- data.raw[sample(1:nrow(data.raw), nrow(data.raw), replace = TRUE),]

  ###
  ### Warnings and errors

  ### Error if CI requested by CI.type ignored
  if (CI != FALSE & is.null(CI.type)){
    stop("Confidence interval requested but CI.type not specified. Choose either 'parametric' or 'bootstrap'.
         If curve.type = 'loess', CI.type cannot be equal to 'parametric'.")
  }

  if (!is.null(CI.type) & !(CI.type %in% c("parametric", "bootstrap"))){
    stop("CI.type must be 'parametric' or 'bootstrap'.")
  }

  if (!is.null(transitions.out)){
    if (sum(c(colSums(tp.pred) == 0)[transitions.out] == TRUE) > 0){
      stop("Calibraiton curves have been requested for transitions into states which have zero probability of occuring.")
    }
  }

  ### Check if transitions.out is only specified for non-zero columns
  if (!is.null(transitions.out)){
    if (sum(c(colSums(tp.pred) == 0)[transitions.out] == TRUE) > 0){
      stop("Calibraiton curves have been requested for transitions into states which have zero probability of occuring.")
    }
  }

  ### Check if CI defined between 0 and 100
  if (CI != FALSE){
    if (CI >= 100 | CI <= 0){
      stop("CI should be a number taking values in (0,100)")
    }
  }

  ###
  ### Initialise inputted data

  ### Assign colnames to predicted transition probabilities (and in data.pred.plot if specified)
  colnames(tp.pred) <- paste("tp.pred", 1:ncol(tp.pred), sep = "")
  if (!is.null(data.pred.plot)){
    colnames(data.pred.plot) <- paste("tp.pred", 1:ncol(tp.pred), sep = "")
  }

  ### Extract what states an individual can move into from state j (states with a non-zero predicted risk)
  valid.transitions <- which(colSums(tp.pred) != 0)

  ### Assign transitions.out if required to do so
  if (is.null(transitions.out)){
    transitions.out <- valid.transitions
  }

  ### Merge data.raw with predicted risks
  data.raw <- data.frame(data.raw, tp.pred[,valid.transitions])

  ### Extract transition matrix from msdata object
  tmat <- attributes(data.mstate)$trans

  ###
  ### Define data.pred.plot if not specified
  ### This is what the "predicted transition probabilities" will be on the x-axis of the calibration plots
  if (is.null(data.pred.plot)){

    ### Identify landmarked individuals
    id.lmk.js <- extract_ids_states(data.mstate = data.mstate,
                                    tmat = attributes(data.mstate)$trans,
                                    j = j,
                                    t.eval = s)

    ### Extract predicted risks of these individuals
    data.pred.plot <- data.raw %>%
      dplyr::filter(id %in% id.lmk.js) %>%
      dplyr::select(dplyr::any_of(paste("tp.pred", 1:6, sep = "")))

    ### Indicator to remember data.pred.plot was defined manually
    manual.data.pred.plot <- FALSE

  } else if (!is.null(data.pred.plot)){

    ### Indicator to remember data.pred.plot was defined manually
    manual.data.pred.plot <- TRUE

  }

  ###
  ### calib_pseudo_func calculates the calibration curves by estimating observed event probabilities
  ### through the pseudo-value approach. It is written in a way to allow bootstrapping to be applied.
  ### If the argument indices is specified to be 1:nrow(data.raw), this will result in the calibration curve
  ### for data.raw with no bootstrapping applied.

  ### The argument "boot.format" specifies whether only the observed event probabilities are outputted (boot.format = TRUE),
  ### or if data is outputted in a dataframe with patient ids and predicted transition probabilities (boot.format = FALSE).

  ### The function can only be used within the boot::boot functionality if boot.format = TRUE, and transitions.out has been specified
  ### to be one state. If a confidence interval is not requested, the data will be outputted using boot.format = FALSE, and > 1 state
  ### can be specified in transitions.out

  ### This allows us to use the same function to generate observed event probabilities, whether
  ### a confidence interval has been requested via bootstrapping or not.

  calib_pseudo_func <- function(data.raw, indices, data.mstate, transitions.out, boot.format = FALSE){

    ### If boot.format = TRUE and requested more than one state, stop
    if (boot.format == TRUE & (length(transitions.out) > 1)){
      stop("CANNOT OUTPUT IN BOOT FORMAT IF REQUESTING OBSERVED EVENT PROBABILITIES FOR MORE
         THAN ONE STATE")
    }

    print("START BOOT SAMPLE")
    # data.raw.in <- data.raw
    # data.mstate.in <- data.mstate
    # indices <- sample(1:nrow(data.raw), nrow(data.raw), replace = TRUE)
    # indices <- 1:nrow(data.raw)
    # boot.format = FALSE
    # transitions.out <- NULL

    ### Create bootstrapped dataset
    data.raw.boot <- data.raw[indices, ]

    ### Create a new id for these individuals (calc_pv_aj relies on each individual having a unique identifier),
    ### meaning the duplicate values in the bootstrapped datasets will cause problems
    data.raw.boot$id2 <- 1:nrow(data.raw.boot)

    ### Create bootstrapped data.mstate
    data.mstate.boot <-
      do.call("rbind",
              lapply(1:nrow(data.raw.boot),
                     function(x) {
                       base::subset(data.mstate, id == data.raw.boot$id[x]) %>%
                         dplyr::mutate(id2 = data.raw.boot$id2[x])
                     }
              )
      )

    ### Apply attribute tmat
    attributes(data.mstate.boot)$trans <- attributes(data.mstate)$trans

    ### Set 'id' to be same as 'id2' in both datasets, as the function calc_pv_aj works by removing individual
    ### with the 'id' variable
    data.mstate.boot$id <- data.mstate.boot$id2
    data.raw.boot$id <- data.raw.boot$id2

    ### Relabel data.mstate.boot and data.raw.boot and remove '.boot' datasets
    data.raw <- data.raw.boot
    data.mstate <- data.mstate.boot
    rm(data.raw.boot, data.mstate.boot)

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

      if (length(transitions.out) < length(valid.transitions)){
        message("Pseudo-values for each state are calculated simultaenously for each individual, because individuals
              are not being grouped by predicted risk. No efficiency gains to only estimate pseudo-values for a subset
              of states, argument 'transitions.out' has been ignored")
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

      if (length(transitions.out) < length(valid.transitions)){
        message("Pseudo-values for each state are calculated simultaenously for each individual, because individuals
              are not being grouped by predicted risk. No efficiency gains to only estimate pseudo-values for a subset
              of states, argument 'transitions.out' has been ignored")
      }

      ###
      ### Split data into groups defined by the variables in group.vars

      ### Create formula to split the dataset by (by group.vars)
      split.formula <- stats::as.formula(paste("~ ", paste(group.vars, collapse = "+"), sep = ""))
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
      pv.out <- dplyr::arrange(pv.out, id)

    } else if (is.null(group.vars) & !is.null(n.pctls)) {

      ### 3) Grouping only by predicted risk

      ### Because we have to do this seperately for each state [in 1) and 2) we
      ### just fitted a single AJ within each group for all states], this process will only loop through
      ### the transitions specified by transitions.out (or valid.transitions if unspecified)

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
                                               breaks =  stats::quantile(data.raw.lmk.js[,paste("tp.pred", state.k, sep = "")],
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
        split.formula <- stats::as.formula(paste("~ ", paste(group.vars, collapse = "+"), sep = ""))
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
                          breaks =  stats::quantile(data.in[,paste("tp.pred", state.k, sep = "")],
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

        ### Combine the pseudo-values in the groups defined by baseline variables
        pv.out.groups.pctls[[state]] <- do.call("rbind", pv.out.groups.pctls[[state]])

        ### END LOOP FOR VARIABLE 'state'
      }

      ### Combine into a single dataset
      pv.out <- Reduce(function(...) merge(..., by = "id", all.x = TRUE), pv.out.groups.pctls)
      ### Arrange manually, because if only one set of p-values to merge (length(transitions.out == 1)),
      ### the merge command will not also arrange the values
      pv.out <- dplyr::arrange(pv.out, id)

    }

    ##############################################
    ### PSEUDO VALUES HAVE NOW BEEN CALCULATED ###
    ##############################################

    ###
    ### Define function to calculate observed event probabilities/calibration plot data,
    ### for given set of pseudo-values (pv) and predicted risks (pred) using loess smoothers
    calc_obs_loess_func <- function(pred, pv, plotdat){

      ### Fit model
      loess.model <- stats::loess(pv ~ pred,
                                  span = loess.span,
                                  degree = loess.degree)

      ### Created observed event probabilities for each individual
      obs <- predict(loess.model, newdata = plotdat)
      obs.data <- data.frame("obs" = obs)

      ### Return obs.data
      return(obs.data)

    }

    ###
    ### Define function to calculate observed event probabilities/calibration plot data,
    ### for given set of pseudo-values (pv) and predicted risks (pred), using a logit link
    ### function and restricted cubic splines
    calc_obs_rcs_func <- function(pred, pv, plotdat, rcs.nk){
      # rcs.nk <- 3
      # pred <- data.raw.lmk.js[,paste("tp.pred", 3, sep = "")]
      # pv <- pv.out[,paste("pstate", state.k, sep = "")]
      # plotdat <- data.pred.plot[,paste("tp.pred", 3, sep = "")]

      ### Create spline terms based on predicted risks
      rcs.pred <- Hmisc::rcspline.eval(pred, nk=rcs.nk, inclx=T)
      colnames(rcs.pred) <- paste("rcs.x", 1:ncol(rcs.pred), sep = "")
      knots.pred <- attr(rcs.pred,"knots")

      ### Create spline terms in plotdat (using same knot locations derived from the predicted risks)
      ### Note that if plotdat == pred, these will be the same
      rcs.plotdat <- data.frame(Hmisc::rcspline.eval(plotdat ,knots = knots.pred, inclx=T))
      colnames(rcs.plotdat) <- paste("rcs.x", 1:ncol(rcs.plotdat), sep = "")

      ### Create dataset in which to fit the model
      data.rcs <- data.frame("pv" = pv, rcs.pred)

      ### Define equation
      eq.LHS <- paste("pv ~ ", sep = "")
      eq.RHS <- paste("rcs.x", 1:ncol(rcs.logit.pred), sep = "", collapse = "+")
      eq.rcs <- stats::formula(paste(eq.LHS, eq.RHS, sep = ""))

      ## Fit the model
      rcs.model <- glm(eq.rcs, data = data.rcs, family = gaussian(link = "logit"), start = rep(0, ncol(rcs.pred) + 1))

      ## Calculate predicted observed probabilities (and confidence intervals if requested using parametric approach)
      ## Note we do not calculate standard errors if confidence interval has been requested using the bootstrap
      if (CI == FALSE){
        ## Predict observed
        obs <- predict(rcs.model, newdata = rcs.plotdat, type = "link")
        ## Put into dataframe
        obs.data <- data.frame("obs" = 1/(1+exp(-obs)))
      } else if (CI != FALSE){
        if (CI.type == "bootstrap"){
          ## Predict observed
          obs <- predict(rcs.model, newdata = rcs.plotdat, type = "link")
          ## Put into dataframe
          obs.data <- data.frame("obs" = 1/(1+exp(-obs)))
        } else if (CI.type == "parametric"){
          ## Predict observed
          obs <- predict(rcs.model, newdata = rcs.plotdat, type = "link", se.fit = TRUE)
          ## Define alpha for CIs
          alpha <- (1-CI/100)/2
          ## Put into dataframe
          obs.data <- data.frame("obs" = 1/(1+exp(-obs$fit)),
                                 "lower" = 1/(1+exp(-(obs$fit - qnorm(1-alpha)*obs$se.fit))),
                                 "upper" = 1/(1+exp(-(obs$fit + qnorm(1-alpha)*obs$se.fit)))
          )
        }

      }

      ### Return obs.data
      return(obs.data)

    }

    ###
    ### Generate observed event probabilities for each state using the calculated pseudo-values
    ### stored in pv.out, and the function for fitting the calibration model: calc_obs_loess_func

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
        obs <- calc_obs_loess_func(pred = data.raw.lmk.js[,paste("tp.pred", state.k, sep = "")],
                                   pv = pv.out[,paste("pstate", state.k, sep = "")],
                                   plotdat = data.pred.plot[,paste("tp.pred", state.k, sep = "")])
      } else if (curve.type == "rcs"){
        obs <- calc_obs_rcs_func(pred = data.raw.lmk.js[,paste("tp.pred", state.k, sep = "")],
                                 pv = pv.out[,paste("pstate", state.k, sep = "")],
                                 plotdat = data.pred.plot[,paste("tp.pred", state.k, sep = "")],
                                 rcs.nk = rcs.nk)
      }

      ### Create output object
      if (manual.data.pred.plot == FALSE) {
        output.object[[state]] <- data.frame(
          "id" = id.lmk.js,
          "pred" = data.pred.plot[,paste("tp.pred", state.k, sep = "")],
          obs)

      } else if (manual.data.pred.plot == TRUE) {
        output.object[[state]] <- data.frame(
          "pred" = data.pred.plot[,paste("tp.pred", state.k, sep = "")],
          obs)
      }

    }

    ###
    ### If boot.format is true, just return the observed event probabilities
    if(boot.format == TRUE){
      output.object <- output.object[[1]]$obs
    }

    return(output.object)

  }

  ##############################################
  ##############################################
  ### END OF FUNCTION calib_pseudo_func ###
  ##############################################
  ##############################################

  ###
  ### Create plotdata object (this contains the observed event probabilities/calibration plot data through
  ### application of calib_pseudo_func).

  ### 1) If a confidence interval was requested using bootstrapping, use calib_pseudo_func in conjuction with boot::boot.
  ### This must be done separately for each state.

  ### 2) If a confidence interval was requested using parametric form, call calib_pseudo_func once, specifying indices to be
  ### the sequence 1:nrow(data.raw), in order to calculate the plot data.

  ### 3) If a confidence interval was not requested, call calib_pseudo_func once, specifying indices to be
  ### the sequence 1:nrow(data.raw), in order to calculate the plot data.

  ### Note that 2) and 3) are the same. This is because the function calib_pseudo_func is dependent on CI and CI.type,
  ### which were defined as input into calc_calib_pv. They will there give different output (as they should) when it is run.

  ### Note that if curve.type = "loess", a request for CI.type = "parametric" will be ignored and no CI is reported.
  ### A stop error is given if a user tries to request this.

  ### If a confidence interval was not requested, run this function once,
  if (CI != FALSE) {

    if (CI.type == "bootstrap"){

      ### Define alpha for CI's
      alpha <- (1-CI/100)/2

      ### Create object to store plot data
      plotdata <- vector("list", length(transitions.out))

      ### Cycle through states
      for (state in 1:length(transitions.out)){

        ### Assign state.k
        state.k <- transitions.out[state]

        ### Put function through bootstrap
        boot.obs <- boot::boot(data.raw,
                               calib_pseudo_func,
                               R = 2,
                               data.mstate = data.mstate,
                               transitions.out = state.k,
                               boot.format = TRUE)

        ### Extract confidence bands
        lower <- apply(boot.obs$t, 2, stats::quantile, probs = alpha, na.rm = TRUE)
        upper <- apply(boot.obs$t, 2, stats::quantile, probs = 1-alpha, na.rm = TRUE)

        ### Produce a warning if any NA values
        if(sum(is.na(boot.obs$t)) > 0){
          print(paste("WARNING, SOME BOOTSTRAPPED OBSERVED RISKS WERE NA FOR STATE", transitions.out[state]))
          print(paste("THERE ARE ", sum(apply(boot.obs$t, 1, function(x) {sum(is.na(x)) > 0})), " ITERATIONS WITH NA's FOR OBSERVED RISKS"))
          print(paste("THE MEAN NUMBER OF NA's IN EACH ITERATION IS", mean(apply(boot.obs$t, 1, function(x) {sum(is.na(x))}))))
        }

        ### Assign output
        if (manual.data.pred.plot == FALSE) {
          plotdata[[state]] <- data.frame("id" = id.lmk.js,
                                          "pred" = data.pred.plot[,paste("tp.pred", state.k, sep = "")],
                                          "obs" = boot.obs$t0,
                                          "obs.lower" = lower,
                                          "obs.upper" = upper)
        } else if (manual.data.pred.plot == TRUE) {
          plotdata[[state]] <- data.frame(
            "pred" = data.pred.plot[,paste("tp.pred", state.k, sep = "")],
            "obs" = boot.obs$t0,
            "obs.lower" = lower,
            "obs.upper" = upper)
        }
      }
    } else if (CI.type == "parametric"){

      ### Note that calib_pseudo_func is dependent on curve.type and CI.type, parameters input to calc_calib_pv
      plotdata <- calib_pseudo_func(data.raw = data.raw,
                                    indices = 1:nrow(data.raw),
                                    data.mstate = data.mstate,
                                    transitions.out = transitions.out)
    }

  } else if (CI == FALSE) {

    ### Calculate calibration plot data
    plotdata <- calib_pseudo_func(data.raw = data.raw,
                                  indices = 1:nrow(data.raw),
                                  data.mstate = data.mstate,
                                  transitions.out = transitions.out)

  }

  ###
  ### Create metadata object
  metadata <- list("valid.transitions" = as.numeric(valid.transitions),
                   "assessed.transitions" = as.numeric(transitions.out),
                   "CI" = CI,
                   "CI.type" = CI.type,
                   "CI.R.boot" = CI.R.boot,
                   "j" = j,
                   "s" = s,
                   "t.eval" = t.eval)

  ###
  ### Crate a combined output object with metadata, as well as plot data
  output.object.comb <- list("plotdata" = plotdata, "metadata" = metadata)

  ### Assign calib_pseudo class
  attr(output.object.comb, "class") <- "calib_psuedo"

  return(output.object.comb)

}


