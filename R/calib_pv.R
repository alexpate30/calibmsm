#' Estimate Aalen-Johansen estimator for a cohort of individuals
#'
#' @description
#' Estimates Aalen-Johansen estimator for the transition probabilities in cohort data.mstate.
#' Estimates transition probabilities at time t if in state j at time 0
#' The Aalen-Johansen estimator for the entire cohort (including individual person_id.eval)
#' is inputted manually (obs.aj), to speed up computational time if calculating pseudo-values
#' for multiple individuals from the same cohort.
#'
#' @param data.mstate Validation data in `msdata` format
#' @param tmat Transition probability matrix
#' @param t Follow up time at which calibration is to be assessed
#' @param j Landmark state at which predictions were made
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
#' @param person_id.eval id of individual to calculate the pseudo-value for
#' @param data.mstate Validation data in `msdata` format
#' @param obs.aj Aalen-Johansen estimator of the transition probabilities in the entire cohort (not excluding person_id.eval)
#' @param tmat Transition probability matrix
#' @param n.cohort Size of cohort (number of unique entries in data.mstate)
#' @param t Follow up time at which calibration is to be assessed
#' @param j Landmark state at which predictions were made
calc_pv_aj <- function(person_id.eval, data.mstate, obs.aj, tmat, n.cohort, t, j){

  # person_id.eval <- 3
  # data.mstate <- data.mstate.lmk.js
  # tmat <- tmat.lmk.js
  # obs.aj <- obs.aj.save
  # t <- 1826 - s
  # n.cohort <- nrow(data.raw.lmk.js)
  # j <- 3
  # s <- 100

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

#' Estimate calibration curves for a multistate model using pseudo-values.
#'
#' @description
#' Creates the underlying data for the calibration curves. `calib_pv`
#' estimates the
#' observed event probabilities for a given set of predicted transition probabilities
#' in a cohort of interest. This is done using techniques for assessing calibration of binary logistic regression models,
#' in combination with inverse probability of censoring weights and landmarking.
#'
#' @param data.mstate Validation data in `msdata` format
#' @param data.raw Validation data in `data.frame` (one row per individual)
#' @param j Landmark state at which predictions were made
#' @param s Landmark time at which predictions were made
#' @param t Follow up time at which calibration is to be assessed
#' @param tp.pred Matrix of predicted transition probabilities at time t, if in state j at time s. There must be a seperate column for the predicted transition probabilities into every state, even if these predicted transition probabilities are 0.
#' @param curve.type Whether calibration curves are estimated using restricted cubic splines ('rcs') or loess smoothers ('loess')
#' @param rcs.nk Number of knots when curves are estimated using restricted cubic splines
#' @param loess.span Span when curves are estimated using loess smoothers
#' @param loess.degree Degree when curves are estimated. using loess smoothers
#' @param group.vars Baseline variables to define groups within which to estimate pseudo-values
#' @param n.pctls Number of percentiles to group individuals by with respect to predicted transition probabilities when estimating pseudo-values
#' @param CI Size of confidence intervals as a %
#' @param CI.type Method for estimating confidence interval (`bootstrap` or `parametric`)
#' @param CI.R.boot Number of bootstrap replicates when estimating the confidence interval for the calibration curve using bootstrapping
#' @param data.pred.plot Data frame or matrix of predicted risks for each possible transition over which to plot the calibration curves. Must have one column for every possible transition.
#' @param transitions.out Transitions for which to calculate calibration curves. Will do all possible transitions if left as NULL.
#'
#' @details
#' Observed event probabilities at time `t` are estimated for predicted
#' transition probabilities `tp.pred` out of state `j` at time `s`.
#' `calib_pv` estimates the observed event probabilities using pseudo-values.
#' Calibraiton curves are generatd by regression the pseudo-values for the transition
#' probabilities on the predicted transition probabilities. REF XXXX. Currently calibration
#' curves can only be produced using loess smoothers. This will be updated to include
#' restricted cubic splines. XXXX Landmarking is applied to only assess calibration
#' in individuals who are uncensored and in state `j` at time `s`.
#'
#' Two datasets for the same cohort of inidividuals must be provided. Firstly `data.mstate` must be a dataset of class `msdata`,
#' generated using the \code{[mstate]} package. This dataset is used to apply the landmarking. Secondly, `data.raw` must be
#' a `data.frame` with one row per individual, containing the desired variables for
#' calculating pseudo-values within (no baseline variables required if `group.vars = NULL`).
#' Confidence intervals for the calibration curves can be estimated using bootstrapping.
#'
#' The calibration curves can be plotted using \code{\link{plot.calib_pv}}.
#'
#' @returns \code{\link{calib_pv}} returns a list containing two elements:
#' \code{plotdata} and \code{metadata}. The \code{plotdata} element contains the
#' data for the calibration curves. This will itself be a list with each element
#' containing calibration plot data for the transition probabilities into each of the possible
#' states. Each list element contains patient ids (\code{id}) from `data.raw`, the predicted
#' transition probabilities (\code{pred}) and the estimated observed event
#' probabilities (\code{obs}). If a confidence interval is requested, upper (`obs.upper`)
#' and lower (`obs.lower`) bounds for the observed event probabilities are also returned.
#' If data.pred.plot is defined manually, column (\code{id}) is not returned.
#' The \code{metadata} element contains metadata including: a vector of the possible transitions,
#' a vector of which transitions calibration curves have been estimated for, the
#' size of the confidence interval, the method for estimating the calibration curve
#' and other user specified information.
#'
#' @examples
#' # Estimate calibration curves for the predicted transition
#' # probabilities at time t = 1826, when predictions were made at time
#' # s = 100 in state j = 3. These predicted transition probabilities are stored in tps100.
#'
#' # Extract the predicted transition probabilities out of state j = 1
#' tp.pred <- dplyr::select(dplyr::filter(tps100, j == 3), any_of(paste("pstate", 1:6, sep = "")))
#'
#' # Now estimate the observed event probabilities for each possible transition.
#'
#' dat.calib.pseudo <- calib_pv(data.mstate = msebmtcal,
#'   data.raw = ebmtcal,
#'   j = 3,
#'   s = 100,
#'   t = 1826,
#'   tp.pred = tp.pred,
#'   group.vars = c("year"),
#'   n.pctls = 2)
#'
#' # The data for each calibration curve are stored in the "plotdata" list
#' # element.
#' str(dat.calib.pseudo)
#'
#' #'
#' @export
calib_pv <- function(data.mstate,
                     data.raw,
                     j,
                     s,
                     t,
                     tp.pred,
                     curve.type = "rcs",
                     rcs.nk = 3,
                     loess.span = 0.75,
                     loess.degree = 2,
                     group.vars = NULL,
                     n.pctls = NULL,
                     CI = FALSE,
                     CI.type = 'parametric',
                     CI.R.boot = NULL,
                     data.pred.plot = NULL,
                     transitions.out = NULL){

  ###########################
  ### Warnings and errors ###
  ###########################

  ### Stop if patients in data.raw are not in data.mstate
  if (!base::all(unique(data.raw$id) %in% unique(data.mstate$id))){
    stop("All patients in data.raw are not contained in data.mstate. Landmarking cannot be applied.")
  }

  ### Warning if patients in data.mstate are not in data.raw
  if (!base::all(unique(data.mstate$id) %in% unique(data.raw$id))){
    warning("All patients in data.mstate are not contained in data.raw. Landmarking can still be applied, but potential mismatch in these two datasets?")
  }

  ### Make sure CI.type is specified to one of the required values
  if (!is.null(CI.type)){
    if (!(CI.type %in% c("parametric", "bootstrap"))){
      stop("CI.type must be 'parametric' or 'bootstrap'.")
    }
  }

  ### Stop if curve.type = "loess" and CI requested but CI.type not specified to be 'bootstrap'.
  if (curve.type == "loess" & CI != FALSE & !is.null(CI.type)){
    if (CI.type != "bootstrap"){
      stop("For curve.type = 'loess', CI.type must be 'bootstrap'.")
    }
  }

  ### Stop if CI.type = "bootstrap" but CI.R.boot not specified
  if (!(CI.type %in% c("parametric", "bootstrap"))){
    stop("CI.type takes values in 'parametric' and 'bootstrap'")
  } else {
    if (CI.type == "bootstrap" & is.null(CI.R.boot)){
      stop("Must specify number of bootstrap replicates for confidence interval using CI.R.boot.")
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

  ################################
  ### Initialise inputted data ###
  ################################

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
                                    t = s)

    ### Extract predicted risks of these individuals
    data.pred.plot <- data.raw |>
      dplyr::filter(id %in% id.lmk.js) |>
      dplyr::select(dplyr::any_of(paste("tp.pred", 1:6, sep = "")))

    ### Indicator to remember data.pred.plot was defined manually
    manual.data.pred.plot <- FALSE

  } else if (!is.null(data.pred.plot)){

    ### Indicator to remember data.pred.plot was defined manually
    manual.data.pred.plot <- TRUE

  }

  ########################################################################
  ### Create function to calculate pseudo-values for a cohort data.raw ###
  ########################################################################

  ### This functions has two steps:
  ###   A) Calculate the pseudo values
  ###   B) Calculate the observed event probabilities based off the pseudo-values

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
                       base::subset(data.mstate, id == data.raw.boot$id[x]) |>
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
    ids.state.js <- extract_ids_states(data.mstate = data.mstate, tmat = tmat, j = j, t = s)

    ### Apply landmarking to data.raw and data.mstate
    data.raw.lmk.js <- data.raw |> base::subset(id %in% ids.state.js)
    data.mstate.lmk.js <- base::subset(data.mstate, id %in% ids.state.js)

    ###
    ### Restructure mstate data so that time s = time 0, and relabel transitions to 1, 2,...

    ### Reduce transition times by s and remove observations which now occur entirely prior to start up
    data.mstate.lmk.js <-
      dplyr::mutate(data.mstate.lmk.js,
                    Tstart = pmax(0, Tstart - s),
                    Tstop = pmax(0, Tstop - s),
                    time = Tstop - Tstart) |>
      base::subset(!(Tstart == 0 & Tstop == 0))

    ###
    ### Remove observations for transitions where no individuals make that transition
    ### Otherwise mstate::msfit will throw out an unneccesary (in this context) warning
    ### This does happen, for example as no patients in state 1 after 100 days move into state 3, despite this being
    ### a possible transition.

    ### Start by identifying which transitions these are
    suppressMessages(zero.transition.table <- data.mstate.lmk.js |>
                       dplyr::group_by(from, to) |>
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

    ######################################
    ### A) CALCULATE THE PSEUDO VALUES ###
    ######################################

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

      ### Calculate the observed Aalen-Johansen once to enable quicker calculation for the pseudo-values
      obs.aj.save <- calc_aj(data.mstate = data.mstate.lmk.js,
                             tmat = tmat.lmk.js,
                             t = t - s,
                             j = j)[["obs.aj"]]


      ### Calculate psuedo-value for each individual
      pv.out <- lapply(data.raw.lmk.js$id, calc_pv_aj,
                       data.mstate = data.mstate.lmk.js,
                       obs.aj = obs.aj.save,
                       tmat = tmat.lmk.js,
                       n.cohort = nrow(data.raw.lmk.js),
                       t = t - s,
                       j = j)

      ### Combine into dataset
      pv.out <- data.frame("id" = data.raw.lmk.js$id, do.call("rbind", pv.out))

    } else if (!is.null(group.vars) & is.null(n.pctls)) {

      ### 2) Grouping only by baseline variables

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
                t = t - s,
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
                                  t = t - s,
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

      ### Because we have to do this seperately for each state, [in 1) and 2) we
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
        print(paste("Calculate pseudo values for state = ", state.k, Sys.time()))

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
                  t = t - s,
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
                                    t = t - s,
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
                      t = t - s,
                      j = j)[["obs.aj"]]

            ### Calculate pseudo-values
            pv.temp <- do.call("rbind",
                               lapply(data.groups.pctls[[state]][[group]][[pctl]]$id, calc_pv_aj,
                                      data.mstate = base::subset(data.mstate.lmk.js, id %in% data.groups.pctls[[state]][[group]][[pctl]]$id),
                                      obs.aj = obs.aj.groups.pctls[[state]][[group]][[pctl]],
                                      tmat = tmat.lmk.js,
                                      n.cohort = nrow(data.groups.pctls[[state]][[group]][[pctl]]),
                                      t = t - s,
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


    #################################################
    ### B) Calculate observed event probabilities ###
    #################################################

    ###
    ### Define a function to do this depending on whether loess or rcs was requested

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
      eq.RHS <- paste("rcs.x", 1:ncol(rcs.plotdat), sep = "", collapse = "+")
      eq.rcs <- stats::formula(paste(eq.LHS, eq.RHS, sep = ""))

      ## Fit the model
      rcs.model <- stats::glm(eq.rcs, data = data.rcs, family = stats::gaussian(link = "logit"), start = rep(0, ncol(rcs.pred) + 1))

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
                                 "obs.lower" = 1/(1+exp(-(obs$fit - stats::qnorm(1-alpha)*obs$se.fit))),
                                 "obs.upper" = 1/(1+exp(-(obs$fit + stats::qnorm(1-alpha)*obs$se.fit)))
          )
        }

      }

      ### Return obs.data
      return(obs.data)

    }

    ###
    ### Now generate the observed event probabilities for each state using the calculated pseudo-values
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

  #########################################
  #########################################
  ### END OF FUNCTION calib_pseudo_func ###
  #########################################
  #########################################

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
  ### which were defined as input into calib_pv. They will there give different output (as they should) when it is run.

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

        ### Print progress
        print(paste("Beginning bootstrapping for state = ", state.k, Sys.time()))

        ### Put function through bootstrap
        boot.obs <- boot::boot(data.raw,
                               calib_pseudo_func,
                               R = CI.R.boot,
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

      ### Note that calib_pseudo_func is dependent on curve.type and CI.type, parameters input to calib_pv
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
                   "curve.type" = curve.type,
                   "CI" = CI,
                   "CI.type" = CI.type,
                   "CI.R.boot" = CI.R.boot,
                   "j" = j,
                   "s" = s,
                   "t" = t,
                   "group.vars" = group.vars,
                   "n.pctls" = n.pctls)

  ###
  ### Crate a combined output object with metadata, as well as plot data
  output.object.comb <- list("plotdata" = plotdata, "metadata" = metadata)

  ### Assign calib_pseudo class
  attr(output.object.comb, "class") <- "calib_pv"

  return(output.object.comb)

}

#' @export
summary.calib_pv <- function(object, ...) {

  cat("There were non-zero predicted transition probabilities into states ",
      paste(object[["metadata"]]$valid.transitions, collapse = ","),  sep = " ")

  cat("\n\nCalibration curves have been estimated for transitions into states ",
      paste(object[["metadata"]]$assessed.transitions, collapse = ","), sep = " ")

  cat("\n\nCalibration was assessed at time ", object[["metadata"]]$t, " and calibration was assessed in a landmarked cohort of individuals in state j = ", object[["metadata"]]$j,
      " at time s = ", object[["metadata"]]$s, sep = "")

  if (isFALSE(object[["metadata"]]$CI)){
    cat("\n\nA confidence interval was not estimated")
  } else {
    cat("\n\nA ",  object[["metadata"]]$CI, "% confidence interval was estimated with", object[["metadata"]]$CI.R.boot, " bootstrap replicates", sep = "")
  }

  if (!is.null(object[["metadata"]]$group.vars)){
    cat("\n\nPseudo-values were calculated within groups specified by covariates", paste(object[["metadata"]]$group.vars, collapse = ","), sep = "")
  }

  if (!is.null(object[["metadata"]]$n.pctls)){
    cat("\n\nPseudo-values were calculated within groups defined by predicted risk of each transition probability. the numbre of groups was", object[["metadata"]]$n.pctls, sep = "")
  }

  cat("\n\nThe estimated calibration curves are stored in list element `plotdata`:\n\n")

  print(lapply(object[["plotdata"]], "head"))

}
