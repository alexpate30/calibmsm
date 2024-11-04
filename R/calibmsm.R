#' Assess the calibration of a multistate model
#' @description
#' Calculates the underlying data for calibration plots of the predicted transition
#' probabilities from a multistate model using three methods.
#' 1) BLR-IPCW: Binary logistic regression with inverse probability of censoring weights.
#' 2) MLR-IPCW: Multinomial logistic regression with inverse probability of censoring
#' weights, based on the nominal calibration  framework of van Hoorde et al. (2014, 2015)
#' 3) Pseudo-values: Pseudo-values estimated using the Aalen-Johansen estimator (Aalen OO, Johansen S, 1978).
#'
#' @param data_raw Validation data in `data.frame` (one row per individual)
#' @param data_ms Validation data in `msdata` format
#' @param j Landmark state at which predictions were made
#' @param s Landmark time at which predictions were made
#' @param t Follow up time at which calibration is to be assessed
#' @param tp_pred Data frame or matrix of predicted transition probabilities at time t, if in state j at time s. There must be a separate column for the predicted transition probabilities into every state, even if these predicted transition probabilities are 0.
#' @param tp_pred_plot Data frame or matrix of predicted risks for each possible transition over which to plot the calibration curves. Argument provided to enable user to apply bootstrapping manually.
#' @param calib_type Whether calibration plots are estimated using BLR-IPCW ('blr'), MLR-IPCW ('mlr') or pseudo-values ('pv')
#' @param curve_type Whether calibration curves are estimated using restricted cubic splines ('rcs') or loess smoothers ('loess')
#' @param rcs_nk Number of knots when curves are estimated using restricted cubic splines
#' @param loess_span Span when curves are estimated using loess smoothers
#' @param loess_degree Degree when curves are estimated_ using loess smoothers
#' @param loess_surface see \code{\link[stats]{loess.control}}
#' @param loess_statistics see \code{\link[stats]{loess.control}}
#' @param loess_trace_hat see \code{\link[stats]{loess.control}}
#' @param loess_cell see \code{\link[stats]{loess.control}}
#' @param loess_iterations see \code{\link[stats]{loess.control}}
#' @param loess_iterTrace see \code{\link[stats]{loess.control}}
#' @param mlr_s_df degrees of freedom of vector spline (see \code{\link[VGAM]{s}})
#' @param mlr_smoother_type Type of smoothing applied. Takes values `s` (see \code{\link[VGAM]{s}}), `sm.ps` (see \code{\link[VGAM]{sm.ps}}) or `sm.os` (see \code{\link[VGAM]{sm.os}}).
#' @param mlr_ps_int the number of equally-spaced B spline intervals in the vector spline smoother (see \code{\link[VGAM]{sm.ps}})
#' @param mlr_degree the degree of B-spline basis in the vector spline smoother (see \code{\link[VGAM]{sm.ps}})
#' @param mlr_niknots number of interior knots (see \code{\link[VGAM]{sm.os}})
#' @param weights Vector of inverse probability of censoring weights
#' @param w_function Custom function for estimating the inverse probability of censoring weights
#' @param w_covs Character vector of variable names to adjust for when calculating inverse probability of censoring weights
#' @param w_landmark_type Whether weights are estimated in all individuals uncensored at time s ('all') or only in individuals uncensored and in state j at time s ('state')
#' @param w_max Maximum bound for inverse probability of censoring weights
#' @param w_stabilised Indicates whether inverse probability of censoring weights should be stabilised or not
#' @param w_max_follow Maximum follow up for model calculating inverse probability of censoring weights. Reducing this to `t` + 1 may aid in the proportional hazards assumption being met in this model.
#' @param pv_group_vars Variables to group by before calculating pseudo-values
#' @param pv_n_pctls Number of percentiles of predicted risk to group by before calculating pseudo-values
#' @param pv_precalc Pre-calculated pseudo-values
#' @param pv_ids Id's of individuals to calculate pseudo-values for
#' @param CI Size of confidence intervals as a %
#' @param CI_type Method for estimating confidence interval (currently restricted to `bootstrap`)
#' @param CI_R_boot Number of bootstrap replicates when estimating the confidence interval for the calibration curve
#' @param CI_seed Seed for bootstrapping procedure
#' @param transitions_out Transitions for which to calculate calibration curves. Will do all possible transitions if left as NULL.
#' @param assess_moderate TRUE/FALSE whether to estimate data for calibration plots
#' @param assess_mean TRUE/FALSE whether to estimate mean calibration
#' @param ... Extra arguments to be passed to w_function (custom function for estimating weights)
#'
#' @details
#' Observed event probabilities at time `t` are estimated for predicted
#' transition probabilities `tp_pred` out of state `j` at time `s`.
#'
#' `calib_type = 'blr'` estimates calibration curves using techniques for assessing
#' the calibration of a binary logistic regression model (Van Calster et al., 2016).
#' A choice between restricted cubic splines and loess smoothers for estimating the
#' calibration curve can be made using `curve_type`. Landmarking (van Houwelingen HC, 2007)
#' is applied to only assess calibration in individuals who are uncensored and in state `j` at time `s`.
#' Calibration can only be assessed in individuals who are also uncensored at time `t`,
#' which is accounted for using inverse probability of censoring weights (Hernan M, Robins J, 2020).
#' See method BLR-IPCW from Pate et al., (2024) for a full explanation of the approach.
#'
#' `calib_type = 'mlr'` estimates calibration scatter plots using a technique for assessing
#' the calibration of multinomial logistic regression models, namely the nominal
#' calibration framework of van Hoorde et al. (2014, 2015). Landmarking (van Houwelingen HC, 2007)
#' is applied to only assess calibration in individuals who are uncensored
#' and in state `j` at time `s`. Calibration can only be assessed in individuals
#' who are also uncensored at time `t`, which is accounted for using inverse probability
#' of censoring weights (Hernan M, Robins J, 2020). See method BLR-IPCW from Pate et al., (2024)
#' for a full explanation of the approach.
#'
#' `calib_type = 'pv'` estimates calibration curves using using pseudo-values (Andersen PK, Pohar Perme M, 2010)
#' calculated using the Aalen-Johansen estimator (Aalen OO, Johansen S, 1978).
#' Calibration curves are generated by regressing the pseudo-values on the predicted transition probabilities.
#' A choice between restricted cubic splines and loess smoothers for estimating the
#' calibration curve can be made using `curve_type`. Landmarking (van Houwelingen HC, 2007)
#' is applied to only assess calibration in individuals who are uncensored and in state `j` at time `s`.
#' The nature of pseudo-values means calibration can be assessed in all landmarked individuals,
#' regardless of their censoring time. See method Pseudo-value approach from Pate et al., (2024)
#' for a full explanation of the approach.
#'
#' Two datasets for the same cohort of inidividuals must be provided. Firstly, `data_raw`
#' must be a `data.frame` with one row per individual containing the variables for
#' the time until censoring (`dtcens`), and an indicator for censoring `dtcens_s`,
#' where (`dtcens_s = 1`) if an individual is censored at time `dtcens`, and `dtcens_s = 0`
#' otherwise. When an individual enters an absorbing state, this prevents censoring
#' from happening (i.e. dtcens_s = 0). `data_raw` must also contain the desired variables
#' for estimating the weights. Secondly, `data_ms` must be a dataset of class `msdata`,
#' generated using the \code{[mstate]} package. This dataset is used to apply the landmarking
#' and identify which state individuals are in at time `t`. While `data_ms` can be
#' derived from `data_raw`, it would be inefficient to do this within `calibmsm::calib_msm`
#' due to the bootstrapping procedure, and therefore they must be inputted seperately.
#'
#' Unless the user specifies the weights using `weights`, the weights are
#' estimated using a cox-proportional hazard model, assuming a linear
#' functional form of the variables defined in `w_covs`. We urge users to
#' specify their own model for estimating the weights. The `weights` argument
#' must be a vector with length equal to the number of rows of `data_raw`.
#'
#' Confidence intervals cannot be produced for the calibration scatter plots (`calib_type = 'mlr'`).
#' For calibration curves estimated using `calib_type = 'blr'`, confidence intervals
#' can only be estimated using bootstrapping (`CI_type = 'bootstrap`). This procedure uses the internal method for
#' estimating weights, we therefore encourage users to specify their own bootstrapping
#' procedure, which incorporates their own model for estimating the weights. Details
#' on how to do this are provided in the vignette \emph{BLR-IPCW-manual-bootstrap}.
#' For calibration curves estimated using `calib_type = 'pv'`, confidence intervals
#' can be estimated using bootstrapping (`CI_type = 'bootstrap`) or parametric formulae (`CI_type = 'parametric`).
#' For computational reasons we recommend using the parametric approach.
#'
#' The calibration plots can be plotted using \code{\link{plot.calib_msm}} and \code{\link{plot.calib_mlr}}.
#'
#' @returns \code{\link{calib_msm}} returns a list containing two elements:
#' \code{plotdata} and \code{metadata}. The \code{plotdata} element contains the
#' data for the calibration plots. This will itself be a list with each element
#' containing calibration plot data for the transition probabilities into each of the possible
#' states. Each list element contains patient ids (\code{id}) from `data_raw`, the predicted
#' transition probabilities (\code{pred}) and the estimated observed event
#' probabilities (\code{obs}). If a confidence interval is requested, upper (`obs_upper`)
#' and lower (`obs_lower`) bounds for the observed event probabilities are also returned.
#' If tp_pred_plot is specified, column (\code{id}) is not returned.
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
#' Pate, A., Sperrin, M., Riley, R. D., Peek, N., Van Staa, T., Sergeant, J. C., Mamas, M. A.,
#' Lip, G. Y. H., Flaherty, M. O., Barrowman, M., Buchan, I., & Martin, G. P.
#' Calibration plots for multistate risk predictions models.
#' \emph{Statistics in Medicine}. 2024;April:1–23. doi: 10.1002/sim.10094.
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
#' tp_pred <- dplyr::select(dplyr::filter(tps0, j == 1), any_of(paste("pstate", 1:6, sep = "")))
#'
#' # Now estimate the observed event probabilities for each possible transition.
#' dat_calib <-
#' calib_msm(data_ms = msebmtcal,
#'  data_raw = ebmtcal,
#'  j=1,
#'  s=0,
#'  t = 1826,
#'  tp_pred = tp_pred,
#'  w_covs = c("year", "agecl", "proph", "match"))
#'
#' # Summarise the output
#' summary(dat_calib)
#'
#' @export
calib_msm <- function(data_ms,
                      data_raw,
                      j,
                      s,
                      t,
                      tp_pred,
                      tp_pred_plot = NULL,
                      calib_type = "blr",
                      curve_type = "rcs",
                      rcs_nk = 3,
                      loess_span = 0.75,
                      loess_degree = 2,
                      loess_surface = c("interpolate", "direct"),
                      loess_statistics = c("approximate", "exact", "none"),
                      loess_trace_hat = c("exact", "approximate"),
                      loess_cell = 0.2,
                      loess_iterations = 4,
                      loess_iterTrace = FALSE,
                      mlr_smoother_type = c("sm.ps", "sm.os", "s"),
                      mlr_ps_int = 4,
                      mlr_degree = 3,
                      mlr_s_df = 4,
                      mlr_niknots = 4,
                      weights = NULL,
                      w_function = NULL,
                      w_covs = NULL,
                      w_landmark_type = "state",
                      w_max = 10,
                      w_stabilised = FALSE,
                      w_max_follow = NULL,
                      pv_group_vars = NULL,
                      pv_n_pctls = NULL,
                      pv_precalc = NULL,
                      pv_ids = NULL,
                      CI = FALSE,
                      CI_type = "bootstrap",
                      CI_R_boot = NULL,
                      CI_seed = NULL,
                      transitions_out = NULL,
                      assess_moderate = TRUE,
                      assess_mean = TRUE,
                      ...){

  # rm(list=ls())
  #
  # devtools::load_all()
  #   data("ebmtcal")
  #   data("msebmtcal")
  #   data("tps0")
  #   data("tps100")
  # calib_type <- "pv"
  # data_raw <- ebmtcal
  # data_ms <- msebmtcal
  # tp_pred <- tps0 |>
  #   subset(j == 1) |>
  #   dplyr::select(paste("pstate", 1:6, sep = ""))
  #
  # data_raw <- ebmtcal[ebmtcal$id %in% 1:50, ]
  # data_ms <- msebmtcal[msebmtcal$id %in% 1:50, ]
  # tp_pred <- dplyr::select(dplyr::filter(tps0, j == 1), any_of(paste("pstate", 1:6, sep = "")))
  # tp_pred <- tp_pred[1:50, ]
  #
  #
  # j <- 1
  # s <- 0
  # t <- 1826
  #
  # curve_type = "rcs"
  # tp_pred_plot = NULL
  # transitions_out = NULL
  # weights = NULL
  #
  # w_covs = NULL
  # w_landmark_type = "state"
  # w_max = 10
  # w_stabilised = FALSE
  # w_max_follow = NULL
  # w_function = NULL
  #
  # CI = FALSE
  # # CI = 95
  # CI_R_boot = 2
  # rcs_nk = 3
  # CI_type = "bootstrap"
  #
  # CI_seed = 1
  # loess_span = 1
  # loess_degree = 1
  #
  # pv_group_vars = c("year")
  # pv_n_pctls = 2
  #
  # mlr_smoother_type = "sm.ps"
  # mlr_ps_int = 4
  # mlr_degree = 3
  # mlr_s_df = 4
  # mlr_niknots = 4
  # assess_moderate = TRUE
  # assess_mean = TRUE
  #
  # pv_precalc = NULL
  # str(data_raw)
  #
  # pv_group_vars = NULL
  # tp_pred <- readRDS("P:/Documents/aaa_incline/DEBUG.tp.pred.rds")
  # data_raw <- readRDS("P:/Documents/aaa_incline/DEBUG.data.raw.rds")
  # data_ms <- readRDS("P:/Documents/aaa_incline/DEBUG.data.ms.rds")
  # pv_n_pctls = 10
  # t <- 2557

  # str(ebmtcal)
  # str(readRDS("P:/Documents/aaa_incline/data.raw.reduc.rds"))
  # calib_type <- "blr"
  #
  # # ## Calculate manual weights
  # # weights_manual <-
  # #   calc_weights(data_ms = msebmtcal,
  # #                data_raw = ebmtcal,
  # #                t = 1826,
  # #                s = 0,
  # #                landmark_type = "state",
  # #                j = 1,
  # #                max_weight = 10,
  # #                stabilised = FALSE)
  # # weights <- weights_manual$ipcw

  ###########################
  ### Warnings and errors ###
  ###########################

  ### Assign arg
  mlr_smoother_type <- match.arg(mlr_smoother_type)

  ### Stop if data_ms is missing the transition matrix, this can happen when using the subset function on data_ms
  if (!("trans" %in% names(attributes(data_ms)))){
    stop("The is no transition matrix (trans) attribute in data_ms, this may have happened when using the subset function to subset an 'msdata' data frame,
         which should have this attribute")
  }
  ### Stop if patients in data_raw are not in data_ms
  if (!("id" %in% colnames(data_raw) & "id" %in% colnames(data_ms))){
    stop("Variable 'id' must be in both data_raw and data_ms. Individuals are identified across the two datasets using this variable.")
  } else if (!base::all(unique(data_raw$id) %in% unique(data_ms$id))){
    stop("All patients in data_raw are not contained in data_ms. Landmarking cannot be applied.")
  }

  ### Stop if not same number of rows in data_raw and tp_pred
  if (nrow(tp_pred) != nrow(data_raw)){
    stop("Number of rows in tp_pred does not match number of rows in data_raw")
  }

  ### Warning if patients in data_ms are not in data_raw
  if (!base::all(unique(data_ms$id) %in% unique(data_raw$id))){
    warning("All patients in data_ms are not contained in data_raw. Landmarking can still be applied, but potential mismatch in these two datasets?")
  }

  ### Stop if variables dtcens and dtcens_s do not exist, and if any NA values for dtcens
  if (!("dtcens" %in% colnames(data_raw)) | !("dtcens_s" %in% colnames(data_raw))){
    stop("data_raw should contains variables dtcens and dtcens_s")
  } else if (!(sum(is.na(data_raw$dtcens)) == 0)){
    stop("NA values found in dtcens. Censoring time must be known for all individuals.")
  }

  ### Stop if weights inputted manually, and confidence interval requested internally
  if ((CI != FALSE) & !is.null(weights)){
    stop("Estimation of confidence interval using internal bootstrapping procedure was requested.
         This is not possible with fixed user-inputted weights as the calibration curves will be incorrect.
         Weights must be calculated seperately within each bootstrapped dataset, this can be done using the internal procedure,
         or with a user-specified function (w_function)")
  }

  ### Ensure appropriate confidence type has been specified
  if (!isFALSE(CI)){
    if (CI >= 100 | CI <= 0){
      stop("CI should be a number taking values in (0,100)")
    } else if (!(CI_type %in% c("parametric", "bootstrap"))){
      stop("CI_type must takes values in 'parametric' or 'bootstrap'")
    } else if (CI_type == "bootstrap" & is.null(CI_R_boot)){
      stop("Must specify number of bootstrap replicates for confidence interval using CI_R_boot.")
    } else if (calib_type %in% c("blr", "mlr") & CI_type == "parametric"){
      stop("Cannot produce a parametric confidence interva for calib_type = 'blr' or 'mlr'")
    } else if (calib_type == "mlr" & assess_moderate == TRUE){
      stop("Cannot produce a confidence interval for moderate calibration plots using method calib_type = 'mlr'")
    }
  }

  ### Check if transitions_out is only specified for non-zero columns
  if (!is.null(transitions_out)){
    if (sum(c(colSums(tp_pred) == 0)[transitions_out] == TRUE) > 0){
      stop("Calibraiton curves have been requested for transitions into states which have zero probability of occuring.")
    }
  }

  ### If vector of weights and custom function for specifying weights both inputted, give error
  if (!is.null(weights) & !is.null(w_function)){
    stop("Cannot specify weights manually and specify a custom function for estimating the weights. Choose one or the other.")
  }

  ### If pseudo-values does not have same number of columns as tp_pred give error
  ### If pseudo-values does not have same number of rows as data_raw give error
  ### If pseudo-values pre-calculated and bootstrapping requested give error
  if (!is.null(pv_precalc)){
    if (nrow(pv_precalc) != nrow(data_raw)){
      stop("pv_precalc must have same number of rows as data_raw. calib_msm assumes landmarking has already been applied to data_raw as part of estimating the pseudo-values")
    } else if (ncol(pv_precalc) != ncol(tp_pred)){
      stop("pv_precalc must have same number of columns as tp_pred")
    } else if (!isFALSE(CI) & CI_type == "bootstrap"){
      stop("Cannot estimate a bootstrapped confidence interval if inputting pre-calculating pseudo-values.")
    }
  }

  ### Stop if calib_type = "AJ" and assess_moderate = TRUE, or parametric confidence interval requested
  if (calib_type == "AJ" & assess_moderate == TRUE){
    stop("Cannot assess moderate calibration for calib_type = 'AJ'")
  } else if (calib_type == "AJ" & CI != FALSE & CI_type == "parametric"){
    stop("Cannot produce parametric confidence intervals for mean calibration assessd using calib_type = 'AJ'")
  }

  ##########################################################
  ### Data preparation and further warnings/error checks ###
  ##########################################################

  ## If a vector of weights has been provided, add it to the dataset
  if (!is.null(weights)){
    ### First check whether it is the correct length (NA's should be present)
    if (length(weights) != nrow(data_raw)){
      stop("Weights vector not same length as data_raw")
    } else {
      data_raw$ipcw <- weights
      weights_provided <- TRUE
    }
  } else if (is.null(weights)){
    weights_provided <- FALSE
  }

  ### If custom function for estimating weights has been inputted ("w_function"),
  ### stop if it does not contain all the arguments from calc_weights
  if (!is.null(w_function)){
    ### stop if w_function doesn't have correct arguments
    if(!all(names(formals(calc_weights)) %in% names(formals(w_function)))){
      stop("Arguments for w_function does not contain those from calibmsm::calc_weights")
    }
    # calc_weights <- w_function
    # print(calc_weights)
  }

  ### If tp_pred_plot is user specified, ensure it has correct number of columns
  if (!is.null(tp_pred_plot)){
    if (ncol(tp_pred_plot) != ncol(tp_pred)){
      stop("Data pred plot must have same number of columns as tp_pred")
    }
  }

  ### Identify valid transitions
  valid_transitions <- identify_valid_transitions(data_raw = data_raw, data_ms = data_ms, j = j, s = s, t = t)

  ### Check there are individuals in state_k at time t for the transitions with non-zero predicted probability
  for (state_k in 1:max(data_ms$to)){
    if (sum(tp_pred[,state_k]) > 0 & !(state_k %in% valid_transitions)){
      stop(paste("There are no individuals in state ", state_k, " at time point ", t,
                 " but there are non-zero predicted probabilities of being in this state according to tp_pred. ",
                 "This issue must be resolved before assessing calibration. ",
                 "Is there a difference in possible transitions between the cohort the model was developed on, and the validation cohort?", sep = "")
      )
    }
  }

  ### Check if there are any valid transitions which have zero predicted probability
  for (state_k in valid_transitions){
    if (sum(tp_pred[,state_k]) == 0){
      stop(paste("There are individuals in state ", state_k, " at time point ", t,
                 " but there is zero predicted probability of being in this state according to tp_pred. ",
                 "This issue must be resolved before assessing calibration. ",
                 "Is there a difference in possible transitions between the cohort the model was developed on, and the validation cohort?", sep = "")
      )
    }
  }

  ### Check for sufficient numbers in each state at time when calibration is being assessed
  ## Create landmarked dataset
  temp_landmark <-  apply_landmark(data_raw = data_raw, data_ms = data_ms, j = j, s = s, t = t, exclude_cens_t = TRUE, data_return = "data_ms")

  ## Identify individuals in state j at time s
  temp_ids_lmk <- lapply(valid_transitions, extract_ids_states, data_ms = temp_landmark, tmat = attributes(data_ms)$trans, t = t)
  if (any(unlist(lapply(temp_ids_lmk, length)) < 50)){
    warning("In the landmark cohort of individuals uncensored and in state j at time s,
    there are some states have less than 50 people at the time at which calibration is being assessed (t).
    Warnings and errors may occur when the models are fitted to estimate the calibration curves due to small sample size.
    This warning has been written to try and intercept some uninformative error messages when the underlying statistical models fail.
    The number to flag this warning (50) has been chosen arbitrarily, and does not constitute a sufficient sample size from a statistical point of view.")
  }
  rm(temp_landmark, temp_ids_lmk)

  ### Assign column names to pv_precalc
  if (!is.null(pv_precalc)){
    if(ncol(pv_precalc) != ncol (tp_pred)){
      stop("pv_precalc must have same number of columns as tp_pred, even if the probability (and therefore pseudo-values) of entering these states is zero")
    } else {
      colnames(pv_precalc) <- paste("pstate", 1:ncol(pv_precalc), sep = "")
    }
  }

  ########################
  ### DATA PREPERATION ###
  ########################

  ### Extract transition matrix from msdata object
  tmat <- attributes(data_ms)$trans

  ### Assign the maximum state an individual may enter
  max_state <- max(data_ms$to)

  ### Assign colnames to predicted transition probabilities (and in tp_pred_plot)
  colnames(tp_pred) <- paste("tp_pred", 1:ncol(tp_pred), sep = "")
  if (!is.null(tp_pred_plot)){
    colnames(tp_pred_plot) <- paste("tp_pred", 1:ncol(tp_pred), sep = "")
  }

  ### Estimate logit transformation of the predicted risks, which are the linear predictors in the binary logistic regression models
  tp_pred_logit <- log(tp_pred[,valid_transitions]/(1-tp_pred[,valid_transitions]))
  colnames(tp_pred_logit) <- paste("tp_pred_logit", valid_transitions, sep = "")

  ### Estimate log-ratios of the predicted risks, which are the linear predictors in the multinomial logistic model
  tp_pred_mlr <- tp_pred[,valid_transitions]
  tp_pred_mlr <- log(tp_pred_mlr[,2:ncol(tp_pred_mlr)]/tp_pred_mlr[,1])
  colnames(tp_pred_mlr) <- paste("mlr_lp", 1:(ncol(tp_pred_mlr)), sep = "")

  ### Add these, along with the predicted risks, to data_raw
  data_raw <- data.frame(data_raw, tp_pred[,valid_transitions], tp_pred_logit, tp_pred_mlr)


  ### If specified, add the predicted risks and logit transformation to tp_pred_plot
  ### Note we do not add the log-ratios, because the tp_pred_plot argument does not work for calib_type = "mlr"
  if (!is.null(tp_pred_plot)){
    tp_pred_plot_logit <- log(tp_pred_plot[,valid_transitions]/(1-tp_pred_plot[,valid_transitions]))
    colnames(tp_pred_plot_logit) <- paste("tp_pred_logit", valid_transitions, sep = "")
    tp_pred_plot <- data.frame(tp_pred_plot[,valid_transitions], tp_pred_plot_logit)
  }

  ### Extract which state individuals are in at time t
  ids_state_list <- vector("list", max_state)
  for (k in valid_transitions){
    ids_state_list[[k]] <- extract_ids_states(data_ms, tmat, k, t)
  }

  ### Create a variable to say which state an individual was in at the time of interest
  ## Create list containing the relevant data
  v1 <- data_raw$id
  m1 <- outer(v1, ids_state_list, FUN = Vectorize('%in%'))
  state_poly <- lapply(split(m1, row(m1)), function(x) (1:max_state)[x])

  ## Change integer(0) values to NA's
  idx <- !sapply(state_poly, length)
  state_poly[idx] <- NA

  ## Add to data_raw
  data_raw <- dplyr::mutate(data_raw, state_poly = unlist(state_poly),
                            state_poly_fac = factor(state_poly))

  ### Create binary variables for each possible state that can be transitioned to
  ## Start by creating NA data.frame
  temp_dummy <- data.frame(matrix(NA, ncol = length(valid_transitions), nrow = nrow(data_raw)))

  ## Create dummy variables
  temp_dummy_calc <- stats::model.matrix(~state_poly_fac - 1, data_raw)

  ## Assign to temp_dummy, for rows where data_raw is not NA
  temp_dummy[!is.na(data_raw$state_poly), ] <- temp_dummy_calc

  ## Assign colnames
  colnames(temp_dummy) <- paste("state", valid_transitions, "_bin", sep = "")

  ### Add to dataset
  data_raw <- cbind(data_raw, temp_dummy)
  rm(temp_dummy)

  ### Define the transitions for which we will be making plots for
  if (is.null(transitions_out)){
    transitions_out <- valid_transitions
  }

  ##########################
  ### Assess calibration ###
  ##########################
  if (calib_type == "blr"){
    output_object <- calib_blr_ipcw(data_raw = data_raw,
                                    data_ms = data_ms,
                                    tp_pred_plot = tp_pred_plot,
                                    j = j,
                                    s = s,
                                    t = t,
                                    curve_type = curve_type,
                                    rcs_nk = rcs_nk,
                                    loess_span = loess_span,
                                    loess_degree = loess_degree,
                                    loess_surface = loess_surface, # no need for loess_statistics argument as never producing parametric confidence interval for calib_blr
                                    loess_trace_hat = loess_trace_hat,
                                    loess_cell = loess_cell,
                                    loess_iterations = loess_iterations,
                                    loess_iterTrace = loess_iterTrace,
                                    weights_provided = weights_provided,
                                    w_function = w_function,
                                    w_covs = w_covs,
                                    w_landmark_type,
                                    w_max = w_max,
                                    w_stabilised = w_stabilised,
                                    w_max_follow = w_max_follow,
                                    CI = CI,
                                    CI_type = CI_type,
                                    CI_R_boot = CI_R_boot,
                                    CI_seed = CI_seed,
                                    transitions_out = transitions_out,
                                    assess_moderate = assess_moderate,
                                    assess_mean = assess_mean, ...)
  } else if (calib_type == "mlr"){
    ### Estimate predicted-obsserved probabilities using the MLR-IPCW method
    output_object <- calib_mlr_ipcw(data_raw = data_raw,
                                    data_ms = data_ms,
                                    j = j,
                                    s = s,
                                    t = t,
                                    weights_provided = weights_provided,
                                    w_function = w_function,
                                    w_covs = w_covs,
                                    w_landmark_type = w_landmark_type,
                                    w_max = w_max,
                                    w_stabilised = w_stabilised,
                                    w_max_follow = w_max_follow,
                                    mlr_smoother_type = mlr_smoother_type,
                                    mlr_ps_int = mlr_ps_int,
                                    mlr_degree = mlr_degree,
                                    mlr_s_df = mlr_s_df,
                                    mlr_niknots = mlr_niknots,
                                    CI = CI,
                                    CI_R_boot = CI_R_boot,
                                    CI_seed = CI_seed,
                                    valid_transitions = valid_transitions,
                                    assess_moderate = assess_moderate,
                                    assess_mean = assess_mean, ...)
  } else if (calib_type == "pv"){
    output_object <- calib_pv(data_raw = data_raw,
                              data_ms = data_ms,
                              tp_pred_plot = tp_pred_plot,
                              j = j,
                              s = s,
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
                              CI_R_boot = CI_R_boot,
                              CI_seed = CI_seed,
                              transitions_out = transitions_out)
  } else if (calib_type == "aj"){
    output_object <- calib_aj(data_raw = data_raw,
                              data_ms = data_ms,
                              j = j,
                              s = s,
                              t = t,
                              pv_group_vars = pv_group_vars,
                              pv_n_pctls = pv_n_pctls,
                              CI = CI,
                              CI_type = CI_type,
                              CI_R_boot = CI_R_boot,
                              CI_seed = CI_seed,
                              transitions_out = transitions_out,
                              valid_transitions = valid_transitions)
  }

  ### Create metadata object
  metadata <- list("valid_transitions" = as.numeric(valid_transitions),
                   "assessed_transitions" = as.numeric(transitions_out),
                   "CI" = CI,
                   "CI_type" = CI_type,
                   "CI_R_boot" = CI_R_boot,
                   "j" = j,
                   "s" = s,
                   "t" = t,
                   "calib_type" = calib_type)
  if (calib_type %in% c("blr", "mlr")){
    metadata[["curve_type"]] <- curve_type
  } else if (calib_type %in% c("pv", "aj")){
    metadata[["pv_group_vars"]] <- pv_group_vars
    metadata[["pv_n_pctls"]] <- pv_n_pctls
    if (calib_type == "pv" & is.null(pv_ids)){
      metadata[["curve_type"]] <- curve_type
    }
  }
  if (CI != FALSE){
    metadata[["CI"]] <- CI
    metadata[["CI_type"]] <- CI_type
    if (CI_type == "bootstrap"){
      metadata[["CI_R_boot"]] <- CI_R_boot
    }
  }

  ### Crate a combined output object with metadata, as well as plot data
  output_object[["metadata"]] <- metadata

  ### Assign classes
  if (calib_type == "blr"){
    class(output_object) <- c("calib_blr", "calib_msm")
  } else if (calib_type == "mlr"){
    class(output_object) <- c("calib_mlr", "calib_msm")
  } else if (calib_type == "pv"){
    class(output_object) <- c("calib_pv", "calib_msm")
  } else if (calib_type == "aj"){
    class(output_object) <- c("calib_aj", "calib_msm")
  }

  ### Return output object
  return(output_object)

}


#' @export
summary.calib_msm <- function(object, ...) {

  cat("The method used to assess calibration was", ifelse(object[["metadata"]]$calib_type == "blr", "BLR-IPCW",
                                                          ifelse(object[["metadata"]]$calib_type == "mlr", "MLR-IPCW",
                                                                 ifelse(object[["metadata"]]$calib_type == "pv", "Pseudo-values with Aalen-Johansen estimator",
                                                                        ifelse(object[["metadata"]]$calib_type == "aj", "Aalen-Johansen estimator")))),  sep = " ")

  cat("\n\nThere were non-zero predicted transition probabilities into states ",
      paste(object[["metadata"]]$valid_transitions, collapse = ","),  sep = " ")

  cat("\n\nCalibration curves have been estimated for transitions into states ",
      paste(object[["metadata"]]$assessed_transitions, collapse = ","), sep = " ")

  cat("\n\nCalibration was assessed at time ", object[["metadata"]]$t, " and calibration was assessed in a landmarked cohort of individuals in state j = ", object[["metadata"]]$j,
      " at time s = ", object[["metadata"]]$s, sep = "")

  if (isFALSE(object[["metadata"]]$CI)){
    cat("\n\nA confidence interval was not estimated")
  } else {
    cat("\n\nA ", object[["metadata"]]$CI, "% confidence interval was estimated with ",
        ifelse(object[["metadata"]]$CI_type == "bootstrap",
               paste("bootstrapping with ", object[["metadata"]]$CI_R_boot, " bootstrap replicates", sep = ""),
               ifelse(object[["metadata"]]$CI_type == "parametric",
                      "a parametric approach")), sep = "")
  }

  if (object[["metadata"]]$calib_type == "pv"){
    if (!is.null(object[["metadata"]]$pv_group_vars)){
      cat("\n\nPseudo-values were calculated within groups specified by covariates", paste(object[["metadata"]]$pv_group_vars, collapse = ","), sep = "")
    }

    if (!is.null(object[["metadata"]]$pv_n_pctls)){
      cat("\n\nPseudo-values were calculated within groups defined by predicted risk of each transition probability. the numbre of groups was", object[["metadata"]]$pv_n_pctls, sep = "")
    }
  }

  if ("plotdata" %in% names(object)){
    cat("\n\nThe estimated data for calibration plots are stored in list element `plotdata`:\n\n")

    print(lapply(object[["plotdata"]], utils::head, 2))
  }

  if ("mean" %in% names(object)){
    cat("\n\nThe estimated mean calibration are stored in list element `mean`:\n\n")

    print(object[["mean"]])
  }


}

#' @export
print.calib_msm <- function(x, ...) {

  print(lapply(x[["plotdata"]], utils::head, 3))

}


#' @export
print.calib_aj <- function(x, ...) {

  print(x[["mean"]])

}

#' @export
plot.calib_aj <- function(x, ...) {

  print("Calibration plots are not available for calib_type = 'aj'")

}

#' Create S3 generic for printing metadata
#'
#' @param x Object generated from \code{\link{calib_msm}}.
#' @param ... Extra arguments
#'
#' @export
metadata <- function(x, ...) {
  UseMethod("metadata")
}

#' @export
metadata.calib_msm <- function(x, ...) {

  print(x[["metadata"]])

}
