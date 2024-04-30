#' Calibration plots for a multistate model'
#' @description
#' Calculates the underlying data for calibration plots of the predicted transition
#' probabilities from a multistate model using three methods.
#' 1) BLR-IPCW: Binary logistic regression with inverse probability of censoring weights.
#' 2) MLR-IPCW: Multinomial logistic regression with inverse probability of censoring
#' weights, based on the nominal calibration  framework of van Hoorde et al. (2014, 2015)
#' 3) Pseudo-values: Pseudo-values estimated using the Aalen-Johansen estimator (Aalen OO, Johansen S, 1978).
#'
#' @param data.raw Validation data in `data.frame` (one row per individual)
#' @param data.mstate Validation data in `msdata` format
#' @param j Landmark state at which predictions were made
#' @param s Landmark time at which predictions were made
#' @param t Follow up time at which calibration is to be assessed
#' @param tp.pred Data frame or matrix of predicted transition probabilities at time t, if in state j at time s. There must be a separate column for the predicted transition probabilities into every state, even if these predicted transition probabilities are 0.
#' @param tp.pred.plot Data frame or matrix of predicted risks for each possible transition over which to plot the calibration curves. Argument provided to enable user to apply bootstrapping manually.
#' @param calib.type Whether calibration plots are estimated using BLR-IPCW ('blr'), MLR-IPCW ('mlr') or pseudo-values ('pv')
#' @param curve.type Whether calibration curves are estimated using restricted cubic splines ('rcs') or loess smoothers ('loess')
#' @param rcs.nk Number of knots when curves are estimated using restricted cubic splines
#' @param loess.span Span when curves are estimated using loess smoothers
#' @param loess.degree Degree when curves are estimated. using loess smoothers
#' @param loess.surface see \code{\link[stats]{loess.control}}
#' @param loess.statistics see \code{\link[stats]{loess.control}}
#' @param loess.trace.hat see \code{\link[stats]{loess.control}}
#' @param loess.cell see \code{\link[stats]{loess.control}}
#' @param loess.iterations see \code{\link[stats]{loess.control}}
#' @param loess.iterTrace see \code{\link[stats]{loess.control}}
#' @param mlr.s.df degrees of freedom of vector spline (see \code{\link[VGAM]{s}})
#' @param mlr.smoother.type Type of smoothing applied. Takes values `s` (see \code{\link[VGAM]{s}}), `sm.ps` (see \code{\link[VGAM]{sm.ps}}) or `sm.os` (see \code{\link[VGAM]{sm.os}}).
#' @param mlr.ps.int the number of equally-spaced B spline intervals in the vector spline smoother (see \code{\link[VGAM]{sm.ps}})
#' @param mlr.degree the degree of B-spline basis in the vector spline smoother (see \code{\link[VGAM]{sm.ps}})
#' @param mlr.niknots number of interior knots (see \code{\link[VGAM]{sm.os}})
#' @param weights Vector of inverse probability of censoring weights
#' @param w.function Custom function for estimating the inverse probability of censoring weights
#' @param w.covs Character vector of variable names to adjust for when calculating inverse probability of censoring weights
#' @param w.landmark.type Whether weights are estimated in all individuals uncensored at time s ('all') or only in individuals uncensored and in state j at time s ('state')
#' @param w.max Maximum bound for inverse probability of censoring weights
#' @param w.stabilised Indicates whether inverse probability of censoring weights should be stabilised or not
#' @param w.max.follow Maximum follow up for model calculating inverse probability of censoring weights. Reducing this to `t` + 1 may aid in the proportional hazards assumption being met in this model.
#' @param pv.group.vars Variables to group by before calculating pseudo-values
#' @param pv.n.pctls Number of percentiles of predicted risk to group by before calculating pseudo-values
#' @param pv.precalc Pre-calculated pseudo-values
#' @param pv.ids Id's of individuals to calculate pseudo-values for
#' @param CI Size of confidence intervals as a %
#' @param CI.type Method for estimating confidence interval (currently restricted to `bootstrap`)
#' @param CI.R.boot Number of bootstrap replicates when estimating the confidence interval for the calibration curve
#' @param CI.seed Seed for bootstrapping procedure
#' @param transitions.out Transitions for which to calculate calibration curves. Will do all possible transitions if left as NULL.
#' @param assess.moderate TRUE/FALSE whether to estimate data for calibration plots
#' @param assess.mean TRUE/FALSE whether to estimate mean calibration
#' @param ... Extra arguments to be passed to w.function (custom function for estimating weights)
#'
#' @details
#' Observed event probabilities at time `t` are estimated for predicted
#' transition probabilities `tp.pred` out of state `j` at time `s`.
#'
#' `calib.type = 'blr'` estimates calibration curves using techniques for assessing
#' the calibration of a binary logistic regression model (Van Calster et al., 2016).
#' A choice between restricted cubic splines and loess smoothers for estimating the
#' calibration curve can be made using `curve.type`. Landmarking (van Houwelingen HC, 2007)
#' is applied to only assess calibration in individuals who are uncensored and in state `j` at time `s`.
#' Calibration can only be assessed in individuals who are also uncensored at time `t`,
#' which is accounted for using inverse probability of censoring weights (Hernan M, Robins J, 2020).
#' See method BLR-IPCW from Pate et al XXXX for a full explanation of the approach.
#'
#' `calib.type = 'mlr'` estimates calibration scatter plots using a technique for assessing
#' the calibration of multinomial logistic regression models, namely the nominal
#' calibration framework of van Hoorde et al. (2014, 2015). Landmarking (van Houwelingen HC, 2007)
#' is applied to only assess calibration in individuals who are uncensored
#' and in state `j` at time `s`. Calibration can only be assessed in individuals
#' who are also uncensored at time `t`, which is accounted for using inverse probability
#' of censoring weights (Hernan M, Robins J, 2020). See method BLR-IPCW from Pate et al XXXX
#' for a full explanation of the approach.
#'
#' `calib.type = 'pv'` estimates calibration curves using using pseudo-values (Andersen PK, Pohar Perme M, 2010)
#' calculated using the Aalen-Johansen estimator (Aalen OO, Johansen S, 1978).
#' Calibration curves are generated by regressing the pseudo-values on the predicted transition probabilities.
#' A choice between restricted cubic splines and loess smoothers for estimating the
#' calibration curve can be made using `curve.type`. Landmarking (van Houwelingen HC, 2007)
#' is applied to only assess calibration in individuals who are uncensored and in state `j` at time `s`.
#' The nature of pseudo-values means calibration can be assessed in all landmarked individuals,
#' regardless of their censoring time. See method PV from Pate et al XXXX
#' for a full explanation of the approach.
#'
#' Two datasets for the same cohort of inidividuals must be provided. Firstly, `data.raw`
#' must be a `data.frame` with one row per individual containing the variables for
#' the time until censoring (`dtcens`), and an indicator for censoring `dtcens.s`,
#' where (`dtcens.s = 1`) if an individual is censored at time `dtcens`, and `dtcens.s = 0`
#' otherwise. When an individual enters an absorbing state, this prevents censoring
#' from happening (i.e. dtcens.s = 0). `data.raw` must also contain the desired variables
#' for estimating the weights. Secondly, `data.mstate` must be a dataset of class `msdata`,
#' generated using the \code{[mstate]} package. This dataset is used to apply the landmarking
#' and identify which state individuals are in at time `t`. While `data.mstate` can be
#' derived from `data.raw`, it would be inefficient to do this within `calibmsm::calib_msm`
#' due to the bootstrapping procedure, and therefore they must be inputted seperately.
#'
#' Unless the user specifies the weights using `weights`, the weights are
#' estimated using a cox-proportional hazard model, assuming a linear
#' functional form of the variables defined in `w.covs`. We urge users to
#' specify their own model for estimating the weights. The `weights` argument
#' must be a vector with length equal to the number of rows of `data.raw`.
#'
#' Confidence intervals cannot be produced for the calibration scatter plots (`calib.type = 'mlr'`).
#' For calibration curves estimated using `calib.type = 'blr'`, confidence intervals
#' can only be estimated using bootstrapping (`CI.type = 'bootstrap`). This procedure uses the internal method for
#' estimating weights, we therefore encourage users to specify their own bootstrapping
#' procedure, which incorporates their own model for estimating the weights. Details
#' on how to do this are provided in the vignette \emph{BLR-IPCW-manual-bootstrap}.
#' For calibration curves estimated using `calib.type = 'pv'`, confidence intervals
#' can be estimated using bootstrapping (`CI.type = 'bootstrap`) or parametric formulae (`CI.type = 'parametric`).
#' For computational reasons we recommend using the parametric approach.
#'
#' The calibration plots can be plotted using \code{\link{plot.calib_msm}} and \code{\link{plot.calib_mlr}}.
#'
#' @returns \code{\link{calib_msm}} returns a list containing two elements:
#' \code{plotdata} and \code{metadata}. The \code{plotdata} element contains the
#' data for the calibration plots. This will itself be a list with each element
#' containing calibration plot data for the transition probabilities into each of the possible
#' states. Each list element contains patient ids (\code{id}) from `data.raw`, the predicted
#' transition probabilities (\code{pred}) and the estimated observed event
#' probabilities (\code{obs}). If a confidence interval is requested, upper (`obs.upper`)
#' and lower (`obs.lower`) bounds for the observed event probabilities are also returned.
#' If tp.pred.plot is specified, column (\code{id}) is not returned.
#' The \code{metadata} element contains metadata including: a vector of the possible transitions,
#' a vector of which transitions calibration curves have been estimated for, the
#' size of the confidence interval, the method for estimating the calibration curve
#' and other user specified information.
#'
#' @references
#' Aalen OO, Johansen S. An Empirical Transition Matrix for Non-Homogeneous Markov Chains Based on Censored Observations.
#' \emph{Scand J Stat}. 1978;5(3):141-150.
#'
#' Andersen PK, Pohar Perme M. Pseudo-observations in survival analysis.
#' \emph{Stat Methods Med Res}. 2010;19(1):71-99. doi:10.1177/0962280209105020
#'
#' Hernan M, Robins J (2020). “12.2 Estimating IP weights via modeling.” In \emph{Causal Inference:
#' What If}, chapter 12.2. Chapman Hall/CRC, Boca Raton.
#'
#' Van Calster B, Nieboer D, Vergouwe Y, De Cock B, Pencina MJ, Steyerberg EW (2016). “A
#' calibration hierarchy for risk models was defined: From utopia to empirical data.” \emph{Journal
#' of Clinical Epidemiology}, 74, 167–176. ISSN 18785921. doi:10.1016/j.jclinepi.2015.
#' 12.005. URL http://dx.doi.org/10.1016/j.jclinepi.2015.12.005
#'
#' Van Hoorde K, Vergouwe Y, Timmerman D, Van Huffel S, Steyerberg W, Van Calster B
#' (2014). “Assessing calibration of multinomial risk prediction models.” \emph{Statistics in Medicine},
#' 33(15), 2585–2596. doi:10.1002/sim.6114.
#'
#' Van Hoorde K, Van Huffel S, Timmerman D, Bourne T, Van Calster B (2015).
#' “A spline-based tool to assess and visualize the calibration of multiclass risk predictions.”
#' \emph{Journal of Biomedical Informatics}, 54, 283–293. ISSN 15320464. doi:10.1016/j.jbi.2014.12.016.
#' URL http://dx.doi.org/10.1016/j.jbi.2014.12.016.
#'
#' van Houwelingen HC (2007). “Dynamic Prediction by Landmarking in Event History Analysis.”
#' \emph{Scandinavian Journal of Statistics}, 34(1), 70–85.
#'
#' Yee TW (2015). \emph{Vector Generalized Linear and Additive Models}. 1 edition.
#' Springer New, NY. ISBN 978-1-4939-4198-8. doi:10.1007/978-1-4939-2818-7.
#' URL https://link.springer.com/book/10.1007/978-1-4939-2818-7.
#'
#' @examples
#' # Estimate BLR-IPCW calibration curves for the predicted transition
#' # probabilities at time t = 1826, when predictions were made at time
#' # s = 0 in state j = 1. These predicted transition probabilities are stored in tps0.
#'
#' # Extract the predicted transition probabilities out of state j = 1
#' tp.pred <- dplyr::select(dplyr::filter(tps0, j == 1), any_of(paste("pstate", 1:6, sep = "")))
#'
#' # Now estimate the observed event probabilities for each possible transition.
#' dat.calib <-
#' calib_msm(data.mstate = msebmtcal,
#'  data.raw = ebmtcal,
#'  j=1,
#'  s=0,
#'  t = 1826,
#'  tp.pred = tp.pred,
#'  w.covs = c("year", "agecl", "proph", "match"))
#'
#' # Summarise the output
#' summary(dat.calib)
#'
#' @export
calib_msm <- function(data.mstate,
                      data.raw,
                      j,
                      s,
                      t,
                      tp.pred,
                      tp.pred.plot = NULL,
                      calib.type = "blr",
                      curve.type = "rcs",
                      rcs.nk = 3,
                      loess.span = 0.75,
                      loess.degree = 2,
                      loess.surface = c("interpolate", "direct"),
                      loess.statistics = c("approximate", "exact", "none"),
                      loess.trace.hat = c("exact", "approximate"),
                      loess.cell = 0.2,
                      loess.iterations = 4,
                      loess.iterTrace = FALSE,
                      mlr.smoother.type = "sm.ps",
                      mlr.ps.int = 4,
                      mlr.degree = 3,
                      mlr.s.df = 4,
                      mlr.niknots = 4,
                      weights = NULL,
                      w.function = NULL,
                      w.covs = NULL,
                      w.landmark.type = "state",
                      w.max = 10,
                      w.stabilised = FALSE,
                      w.max.follow = NULL,
                      pv.group.vars = NULL,
                      pv.n.pctls = NULL,
                      pv.precalc = NULL,
                      pv.ids = NULL,
                      CI = FALSE,
                      CI.type = "bootstrap",
                      CI.R.boot = NULL,
                      CI.seed = 1,
                      transitions.out = NULL,
                      assess.moderate = TRUE,
                      assess.mean = TRUE,
                      ...){

  # rm(list=ls())
  #
  # devtools::load_all()
  #   data("ebmtcal")
  #   data("msebmtcal")
  #   data("tps0")
  #   data("tps100")
  # calib.type <- "pv"
  # data.raw <- ebmtcal
  # data.mstate <- msebmtcal
  # tp.pred <- tps0 |>
  #   subset(j == 1) |>
  #   dplyr::select(paste("pstate", 1:6, sep = ""))
  #
  # data.raw <- ebmtcal[ebmtcal$id %in% 1:50, ]
  # data.mstate <- msebmtcal[msebmtcal$id %in% 1:50, ]
  # tp.pred <- dplyr::select(dplyr::filter(tps0, j == 1), any_of(paste("pstate", 1:6, sep = "")))
  # tp.pred <- tp.pred[1:50, ]
  #
  #
  # j <- 1
  # s <- 0
  # t <- 1826
  #
  # curve.type = "rcs"
  # tp.pred.plot = NULL
  # transitions.out = NULL
  # weights = NULL
  #
  # w.covs = NULL
  # w.landmark.type = "state"
  # w.max = 10
  # w.stabilised = FALSE
  # w.max.follow = NULL
  # w.function = NULL
  #
  # CI = FALSE
  # # CI = 95
  # CI.R.boot = 2
  # rcs.nk = 3
  # CI.type = "bootstrap"
  #
  # CI.seed = 1
  # loess.span = 1
  # loess.degree = 1
  #
  # pv.group.vars = c("year")
  # pv.n.pctls = 2
  #
  # mlr.smoother.type = "sm.ps"
  # mlr.ps.int = 4
  # mlr.degree = 3
  # mlr.s.df = 4
  # mlr.niknots = 4
  # assess.moderate = TRUE
  # assess.mean = TRUE
  #
  # pv.precalc = NULL
  # str(data.raw)
  #
  # pv.group.vars = NULL
  # tp.pred <- readRDS("P:/Documents/aaa_incline/DEBUG.tp.pred.rds")
  # data.raw <- readRDS("P:/Documents/aaa_incline/DEBUG.data.raw.rds")
  # data.mstate <- readRDS("P:/Documents/aaa_incline/DEBUG.data.mstate.rds")
  # pv.n.pctls = 10
  # t <- 2557

  # str(ebmtcal)
  # str(readRDS("P:/Documents/aaa_incline/data.raw.reduc.rds"))
  # calib.type <- "blr"
  #
  # # ## Calculate manual weights
  # # weights.manual <-
  # #   calc_weights(data.mstate = msebmtcal,
  # #                data.raw = ebmtcal,
  # #                t = 1826,
  # #                s = 0,
  # #                landmark.type = "state",
  # #                j = 1,
  # #                max.weight = 10,
  # #                stabilised = FALSE)
  # # weights <- weights.manual$ipcw

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

  ### Stop if variables dtcens and dtcens.s do not exist, and if any NA values for dtcens
  if (!("dtcens" %in% colnames(data.raw)) | !("dtcens.s" %in% colnames(data.raw))){
    stop("data.raw should contains variables dtcens and dtcens.s")
  } else if (!(sum(is.na(data.raw$dtcens)) == 0)){
    stop("NA values found in dtcens. Censoring time must be known for all individuals.")
  }

  ### Stop if weights inputted manually, and confidence interval requested internally
  if ((CI != FALSE) & !is.null(weights)){
    stop("Estimation of confidence interval using internal bootstrapping procedure was requested.
         This is not possible with fixed user-inputted weights as the calibration curves will be incorrect.
         Weights must be calculated seperately within each bootstrapped dataset, this can be done using the internal procedure,
         or with a user-specified function (w.function)")
  }

  ### Ensure appropriate confidence type has been specified
  if (!isFALSE(CI)){
    if (CI >= 100 | CI <= 0){
      stop("CI should be a number taking values in (0,100)")
    } else if (!(CI.type %in% c("parametric", "bootstrap"))){
      stop("CI.type must takes values in 'parametric' or 'bootstrap'")
    } else if (CI.type == "bootstrap" & is.null(CI.R.boot)){
      stop("Must specify number of bootstrap replicates for confidence interval using CI.R.boot.")
    } else if (calib.type %in% c("blr", "mlr") & CI.type == "parametric"){
      stop("Cannot produce a parametric confidence interva for calib.type = 'blr' or 'mlr'")
    } else if (calib.type == "mlr" & assess.moderate == TRUE){
      stop("Cannot produce a confidence interval for moderate calibration plots using method calib.type = 'mlr'")
    }
  }

  ### Check if transitions.out is only specified for non-zero columns
  if (!is.null(transitions.out)){
    if (sum(c(colSums(tp.pred) == 0)[transitions.out] == TRUE) > 0){
      stop("Calibraiton curves have been requested for transitions into states which have zero probability of occuring.")
    }
  }

  ### If vector of weights and custom function for specifying weights both inputted, give error
  if (!is.null(weights) & !is.null(w.function)){
    stop("Cannot specify weights manually and specify a custom function for estimating the weights. Choose one or the other.")
  }

  ### If pseudo-values does not have same number of columns as tp.pred give error
  ### If pseudo-values does not have same number of rows as data.raw give error
  ### If pseudo-values pre-calculated and bootstrapping requested give error
  if (!is.null(pv.precalc)){
    if (nrow(pv.precalc) != nrow(data.raw)){
      stop("pv.precalc must have same number of rows as data.raw. calib_msm assumes landmarking has already been applied to data.raw as part of estimating the pseudo-values")
    } else if (ncol(pv.precalc) != ncol(tp.pred)){
      stop("pv.precalc must have same number of columns as tp.pred")
    } else if (!isFALSE(CI) & CI.type == "bootstrap"){
      stop("Cannot estimate a bootstrapped confidence interval if inputting pre-calculating pseudo-values.")
    }
  }

  ### Stop if calib.type = "AJ" and assess.moderate = TRUE, or parametric confidence interval requested
  if (calib.type == "AJ" & assess.moderate == TRUE){
    stop("Cannot assess moderate calibration for calib.type = 'AJ'")
  } else if (calib.type == "AJ" & CI != FALSE & CI.type == "parametric"){
    stop("Cannot produce parametric confidence intervals for mean calibration assessd using calib.type = 'AJ'")
  }

  ##########################################################
  ### Data preparation and further warnings/error checks ###
  ##########################################################

  ## If a vector of weights has been provided, add it to the dataset
  if (!is.null(weights)){
    ### First check whether it is the correct length (NA's should be present)
    if (length(weights) != nrow(data.raw)){
      stop("Weights vector not same length as data.raw")
    } else {
      data.raw$ipcw <- weights
      weights.provided <- TRUE
    }
  } else if (is.null(weights)){
    weights.provided <- FALSE
  }

  ### If custom function for estimating weights has been inputted ("w.function"),
  ### stop if it does not contain all the arguments from calc_weights
  if (!is.null(w.function)){
    ### stop if w.function doesn't have correct arguments
    if(!all(names(formals(calc_weights)) %in% names(formals(w.function)))){
      stop("Arguments for w.function does not contain those from calibmsm::calc_weights")
    }
    # calc_weights <- w.function
    # print(calc_weights)
  }

  ### If tp.pred.plot is user specified, ensure it has correct number of columns
  if (!is.null(tp.pred.plot)){
    if (ncol(tp.pred.plot) != ncol(tp.pred)){
      stop("Data pred plot must have same number of columns as tp.pred")
    }
  }

  ### Identify valid transitions
  valid.transitions <- identify_valid_transitions(data.raw = data.raw, data.mstate = data.mstate, j = j, s = s, t = t)

  ### Check there are individuals in state.k at time t for the transitions with non-zero predicted probability
  for (state.k in 1:max(data.mstate$to)){
    if (sum(tp.pred[,state.k]) > 0 & !(state.k %in% valid.transitions)){
      stop(paste("There are no individuals in state ", state.k, " at time point ", t,
                 " but there are non-zero predicted probabilities of being in this state according to tp.pred. ",
                 "This issue must be resolved before assessing calibration. ",
                 "Is there a difference in possible transitions between the cohort the model was developed on, and the validation cohort?", sep = "")
      )
    }
  }

  ### Check if there are any valid transitions which have zero predicted probability
  for (state.k in valid.transitions){
    if (sum(tp.pred[,state.k]) == 0){
      stop(paste("There are individuals in state ", state.k, " at time point ", t,
                 " but there is zero predicted probability of being in this state according to tp.pred. ",
                 "This issue must be resolved before assessing calibration. ",
                 "Is there a difference in possible transitions between the cohort the model was developed on, and the validation cohort?", sep = "")
      )
    }
  }

  ### Assign column names to pv.precalc
  if (!is.null(pv.precalc)){
    if(ncol(pv.precalc) != ncol (tp.pred)){
      stop("pv.precalc must have same number of columns as tp.pred, even if the probability (and therefore pseudo-values) of entering these states is zero")
    } else {
      colnames(pv.precalc) <- paste("pstate", 1:ncol(pv.precalc), sep = "")
    }
  }

  ########################
  ### DATA PREPERATION ###
  ########################

  ### Extract transition matrix from msdata object
  tmat <- attributes(data.mstate)$trans

  ### Assign the maximum state an individual may enter
  max.state <- max(data.mstate$to)

  ### Assign colnames to predicted transition probabilities (and in tp.pred.plot)
  colnames(tp.pred) <- paste("tp.pred", 1:ncol(tp.pred), sep = "")
  if (!is.null(tp.pred.plot)){
    colnames(tp.pred.plot) <- paste("tp.pred", 1:ncol(tp.pred), sep = "")
  }

  ### Estimate logit transformation of the predicted risks, which are the linear predictors in the binary logistic regression models
  tp.pred.logit <- log(tp.pred[,valid.transitions]/(1-tp.pred[,valid.transitions]))
  colnames(tp.pred.logit) <- paste("tp.pred.logit", valid.transitions, sep = "")

  ### Estimate log-ratios of the predicted risks, which are the linear predictors in the multinomial logistic model
  tp.pred.mlr <- tp.pred[,valid.transitions]
  tp.pred.mlr <- log(tp.pred.mlr[,2:ncol(tp.pred.mlr)]/tp.pred.mlr[,1])
  colnames(tp.pred.mlr) <- paste("mlr.lp", 1:(ncol(tp.pred.mlr)), sep = "")

  ### Add these, along with the predicted risks, to data.raw
  data.raw <- data.frame(data.raw, tp.pred[,valid.transitions], tp.pred.logit, tp.pred.mlr)


  ### If specified, add the predicted risks and logit transformation to tp.pred.plot
  ### Note we do not add the log-ratios, because the tp.pred.plot argument does not work for calib_type = "mlr"
  if (!is.null(tp.pred.plot)){
    tp.pred.plot.logit <- log(tp.pred.plot[,valid.transitions]/(1-tp.pred.plot[,valid.transitions]))
    colnames(tp.pred.plot.logit) <- paste("tp.pred.logit", valid.transitions, sep = "")
    tp.pred.plot <- data.frame(tp.pred.plot[,valid.transitions], tp.pred.plot.logit)
  }

  ### Extract which state individuals are in at time t
  ids.state.list <- vector("list", max.state)
  for (k in valid.transitions){
    ids.state.list[[k]] <- extract_ids_states(data.mstate, tmat, k, t)
  }

  ### Create a variable to say which state an individual was in at the time of interest
  ## Create list containing the relevant data
  v1 <- data.raw$id
  m1 <- outer(v1, ids.state.list, FUN = Vectorize('%in%'))
  state.poly <- lapply(split(m1, row(m1)), function(x) (1:max.state)[x])

  ## Change integer(0) values to NA's
  idx <- !sapply(state.poly, length)
  state.poly[idx] <- NA

  ## Add to data.raw
  data.raw <- dplyr::mutate(data.raw, state.poly = unlist(state.poly),
                            state.poly.fac = factor(state.poly))

  ### Create binary variables for each possible state that can be transitioned to
  ## Start by creating NA data.frame
  temp.dummy <- data.frame(matrix(NA, ncol = length(valid.transitions), nrow = nrow(data.raw)))

  ## Create dummy variables
  temp.dummy.calc <- stats::model.matrix(~state.poly.fac - 1, data.raw)

  ## Assign to temp.dummy, for rows where data.raw is not NA
  temp.dummy[!is.na(data.raw$state.poly), ] <- temp.dummy.calc

  ## Assign colnames
  colnames(temp.dummy) <- paste("state", valid.transitions, ".bin", sep = "")

  ### Add to dataset
  data.raw <- cbind(data.raw, temp.dummy)
  rm(temp.dummy)

  ### Define the transitions for which we will be making plots for
  if (is.null(transitions.out)){
    transitions.out <- valid.transitions
  }

  ##########################
  ### Assess calibration ###
  ##########################
  if (calib.type == "blr"){
    output.object <- calib_blr_ipcw(data.raw = data.raw,
                                    data.mstate = data.mstate,
                                    tp.pred.plot = tp.pred.plot,
                                    j = j,
                                    s = s,
                                    t = t,
                                    curve.type = curve.type,
                                    rcs.nk = rcs.nk,
                                    loess.span = loess.span,
                                    loess.degree = loess.degree,
                                    loess.surface = loess.surface, # no need for loess.statistics argument as never producing parametric confidence interval for calib_blr
                                    loess.trace.hat = loess.trace.hat,
                                    loess.cell = loess.cell,
                                    loess.iterations = loess.iterations,
                                    loess.iterTrace = loess.iterTrace,
                                    weights.provided = weights.provided,
                                    w.function = w.function,
                                    w.covs = w.covs,
                                    w.landmark.type,
                                    w.max = w.max,
                                    w.stabilised = w.stabilised,
                                    w.max.follow = w.max.follow,
                                    CI = CI,
                                    CI.type = CI.type,
                                    CI.R.boot = CI.R.boot,
                                    CI.seed = CI.seed,
                                    transitions.out = transitions.out,
                                    assess.moderate = assess.moderate,
                                    assess.mean = assess.mean, ...)
  } else if (calib.type == "mlr"){
    ### Estimate predicted-obsserved probabilities using the MLR-IPCW method
    output.object <- calib_mlr_ipcw(data.raw = data.raw,
                                    data.mstate = data.mstate,
                                    j = j,
                                    s = s,
                                    t = t,
                                    weights.provided = weights.provided,
                                    w.function = w.function,
                                    w.covs = w.covs,
                                    w.landmark.type = w.landmark.type,
                                    w.max = w.max,
                                    w.stabilised = w.stabilised,
                                    w.max.follow = w.max.follow,
                                    mlr.smoother.type = mlr.smoother.type,
                                    mlr.ps.int = mlr.ps.int,
                                    mlr.degree = mlr.degree,
                                    mlr.s.df = mlr.s.df,
                                    mlr.niknots = mlr.niknots,
                                    CI = CI,
                                    CI.R.boot = CI.R.boot,
                                    CI.seed = CI.seed,
                                    valid.transitions = valid.transitions,
                                    assess.moderate = assess.moderate,
                                    assess.mean = assess.mean, ...)
  } else if (calib.type == "pv"){
    output.object <- calib_pv(data.raw = data.raw,
                              data.mstate = data.mstate,
                              tp.pred.plot = tp.pred.plot,
                              j = j,
                              s = s,
                              t = t,
                              curve.type = curve.type,
                              rcs.nk = rcs.nk,
                              loess.span = loess.span,
                              loess.degree = loess.degree,
                              loess.surface = loess.surface,
                              loess.statistics = loess.statistics,
                              loess.trace.hat = loess.trace.hat,
                              loess.cell = loess.cell,
                              loess.iterations = loess.iterations,
                              loess.iterTrace = loess.iterTrace,
                              pv.group.vars = pv.group.vars,
                              pv.n.pctls = pv.n.pctls,
                              pv.precalc = pv.precalc,
                              pv.ids = pv.ids,
                              CI = CI,
                              CI.type = CI.type,
                              CI.R.boot = CI.R.boot,
                              CI.seed = CI.seed,
                              transitions.out = transitions.out)
  } else if (calib.type == "aj"){
    output.object <- calib_aj(data.raw = data.raw,
                              data.mstate = data.mstate,
                              j = j,
                              s = s,
                              t = t,
                              pv.group.vars = pv.group.vars,
                              pv.n.pctls = pv.n.pctls,
                              CI = CI,
                              CI.type = CI.type,
                              CI.R.boot = CI.R.boot,
                              CI.seed = CI.seed,
                              transitions.out = transitions.out,
                              valid.transitions = valid.transitions)
  }

  ### Create metadata object
  metadata <- list("valid.transitions" = as.numeric(valid.transitions),
                   "assessed.transitions" = as.numeric(transitions.out),
                   "CI" = CI,
                   "CI.type" = CI.type,
                   "CI.R.boot" = CI.R.boot,
                   "j" = j,
                   "s" = s,
                   "t" = t,
                   "calib.type" = calib.type)
  if (calib.type %in% c("blr", "mlr")){
    metadata[["curve.type"]] <- curve.type
  } else if (calib.type %in% c("pv", "aj")){
    metadata[["pv.group.vars"]] <- pv.group.vars
    metadata[["pv.n.pctls"]] <- pv.n.pctls
    if (calib.type == "pv" & is.null(pv.ids)){
      metadata[["curve.type"]] <- curve.type
    }
  }
  if (CI != FALSE){
    metadata[["CI"]] <- CI
    metadata[["CI.type"]] <- CI.type
    if (CI.type == "bootstrap"){
      metadata[["CI.R.boot"]] <- CI.R.boot
    }
  }

  ### Crate a combined output object with metadata, as well as plot data
  output.object[["metadata"]] <- metadata

  ### Assign classes
  if (calib.type == "blr"){
    class(output.object) <- c("calib_blr", "calib_msm")
  } else if (calib.type == "mlr"){
    class(output.object) <- c("calib_mlr", "calib_msm")
  } else if (calib.type == "pv"){
    class(output.object) <- c("calib_pv", "calib_msm")
  } else if (calib.type == "aj"){
    class(output.object) <- c("calib_aj", "calib_msm")
  }

  ### Return output object
  return(output.object)

}


#' @export
summary.calib_msm <- function(object, ...) {

  cat("The method used to assess calibration was", ifelse(object[["metadata"]]$calib.type == "blr", "BLR-IPCW",
                                                          ifelse(object[["metadata"]]$calib.type == "mlr", "MLR-IPCW",
                                                                 ifelse(object[["metadata"]]$calib.type == "pv", "Pseudo-values with Aalen-Johansen estimator"))),  sep = " ")

  cat("\n\nThere were non-zero predicted transition probabilities into states ",
      paste(object[["metadata"]]$valid.transitions, collapse = ","),  sep = " ")

  cat("\n\nCalibration curves have been estimated for transitions into states ",
      paste(object[["metadata"]]$assessed.transitions, collapse = ","), sep = " ")

  cat("\n\nCalibration was assessed at time ", object[["metadata"]]$t, " and calibration was assessed in a landmarked cohort of individuals in state j = ", object[["metadata"]]$j,
      " at time s = ", object[["metadata"]]$s, sep = "")

  if (isFALSE(object[["metadata"]]$CI)){
    cat("\n\nA confidence interval was not estimated")
  } else {
    cat("\n\nA ", object[["metadata"]]$CI, "% confidence interval was estimated with ",
        ifelse(object[["metadata"]]$CI.type == "bootstrap",
               paste("bootstrapping with ", object[["metadata"]]$CI.R.boot, " bootstrap replicates", sep = ""),
               ifelse(object[["metadata"]]$CI.type == "parametric",
                      "a parametric approach")), sep = "")
  }

  if (object[["metadata"]]$calib.type == "pv"){
    if (!is.null(object[["metadata"]]$pv.group.vars)){
      cat("\n\nPseudo-values were calculated within groups specified by covariates", paste(object[["metadata"]]$pv.group.vars, collapse = ","), sep = "")
    }

    if (!is.null(object[["metadata"]]$pv.n.pctls)){
      cat("\n\nPseudo-values were calculated within groups defined by predicted risk of each transition probability. the numbre of groups was", object[["metadata"]]$pv.n.pctls, sep = "")
    }
  }

  cat("\n\nThe estimated calibration plots are stored in list element `plotdata`:\n\n")

  print(lapply(object[["plotdata"]], utils::head, 3))

}

#' @export
print.calib_msm <- function(x, ...) {

  print(lapply(x[["plotdata"]], utils::head, 10))

}


#' @export
print.calib_aj <- function(x, ...) {

  print(x[["mean"]])

}
