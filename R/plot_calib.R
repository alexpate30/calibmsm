#' Plots calibration curves estimated using \code{\link{calib_blr}}.
#'
#' @description
#' Plots calibration curves for the transition probabilities of a multistate model
#' estimated using \code{\link{calib_blr}}.
#'
#' @param x Object of class 'calib_blr' generated from \code{\link{calib_blr}}.
#' @param ... Other
#' @param combine Whether to combine into one plot using ggarrange, or return as a list of individual plots
#' @param ncol Number of columns for combined calibration plot
#' @param nrow Number of rows for combined calibration plot
#' @param transparency.rug Degree of transparency for the density rug plot along each axis
#'
#' @returns If `combine = TRUE`, returns an object of classes `gg`, `ggplot`, and `ggarrange`,
#' as all ggplots have been combined into one object. If `combine = FALSE`, returns an object of
#' class `list`, each element containing an object of class `gg` and `ggplot`.
#'
#' @examples
#' # Estimate and plot BLR-IPCW calibration curves for the predicted transition
#' # probabilities at time t = 1826, when predictions were made at time
#' # s = 0 in state j = 1. These predicted transition probabilities are stored in tps0.
#'
#' # Extract the predicted transition probabilities out of state j = 1
#' tp.pred <- dplyr::select(dplyr::filter(tps0, j == 1), any_of(paste("pstate", 1:6, sep = "")))
#'
#' # Now estimate the observed event probabilities for each possible transition.
#' dat.calib.blr <-
#' calib_blr(data.mstate = msebmtcal,
#'  data.raw = ebmtcal,
#'  j=1,
#'  s=0,
#'  t = 1826,
#'  tp.pred = tp.pred,
#'  w.covs = c("year", "agecl", "proph", "match"))
#'
#'  # These are then plotted
#'  plot(dat.calib.blr, combine = TRUE, nrow = 2, ncol = 3)
#'
#' @importFrom graphics plot
#' @export
plot.calib_blr <- function(x, ..., combine = TRUE, ncol = NULL, nrow = NULL, transparency.rug = 0.1){

  ### Extract plot data and relevant metadata
  object.in <- x
  plot.data <- object.in[["plotdata"]]
  assessed.transitions <- object.in[["metadata"]][["assessed.transitions"]]
  CI <- object.in[["metadata"]][["CI"]]

  if (CI != FALSE){
    ### Create list to store plots
    plots.list <- vector("list", length(assessed.transitions))

    for (k in 1:length(assessed.transitions)){
      ### Assign plot data
      plot.data.k <- plot.data[[k]]

      ### Assign state of interest
      state.k <- assessed.transitions[k]

      ### Pivot longer to create data for ggplot and assign appropriate labels
      plot.data.k.longer <- tidyr::pivot_longer(plot.data.k, cols = c(obs, obs.upper, obs.lower), names_to = "line.group")
      plot.data.k.longer <- dplyr::mutate(plot.data.k.longer,
                                   line.group = base::factor(line.group),
                                   mapping = dplyr::case_when(line.group == "obs" ~ 1,
                                                       line.group %in% c("obs.upper", "obs.lower") ~ 2),
                                   mapping = base::factor(mapping))

      levels(plot.data.k.longer$line.group) <- c("Calibration", "Upper", "Lower")
      levels(plot.data.k.longer$mapping) <- c("Calibration", "95% CI")

      ### Create the plots
      plots.list[[k]] <- ggplot2::ggplot(data = plot.data.k.longer |> dplyr::arrange(pred) |> dplyr::select(pred, line.group, value, mapping)) +
        ggplot2::geom_line(ggplot2::aes(x = pred, y = value, group = line.group, color = mapping)) +
        ggplot2::geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
        ggplot2::xlab("Predicted risk") + ggplot2::ylab("Observed risk") +
        ggplot2::xlim(c(0, max(plot.data.k.longer$pred))) +
        ggplot2::ylim(c(min(plot.data.k.longer$value), max(plot.data.k.longer$value))) +
        ggplot2::geom_rug(data = plot.data.k.longer |> dplyr::arrange(pred) |> dplyr::select(pred, line.group, value, mapping) |> subset(line.group == "Calibration"),
                          ggplot2::aes(x = pred, y = value), col = grDevices::rgb(1, 0, 0, alpha = transparency.rug)) +
        ggplot2::ggtitle(paste("State ", state.k, sep = ""))
    }
  } else if (CI == FALSE){
    ### Create list to store plots
    plots.list <- vector("list", length(assessed.transitions))
    for (k in 1:length(assessed.transitions)){

      ### Assign plot data
      plot.data.k <- plot.data[[k]]

      ### Assign state of interest
      state.k <- assessed.transitions[k]

      ### Create the plots
      plots.list[[k]] <- ggplot2::ggplot(data = plot.data.k |> dplyr::arrange(pred) |> dplyr::select(id, pred, obs)) +
        ggplot2::geom_line(ggplot2::aes(x = pred, y = obs)) +
        ggplot2::geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
        ggplot2::xlab("Predicted risk") + ggplot2::ylab("Observed risk") +
        ggplot2::xlim(c(0, max(plot.data.k$pred))) +
        ggplot2::ylim(c(min(plot.data.k$obs), max(plot.data.k$obs))) +
        ggplot2::geom_rug(ggplot2::aes(x = pred, y = obs), col = grDevices::rgb(1, 0, 0, alpha = transparency.rug)) +
        ggplot2::theme(legend.position = "none") +
        ggplot2::ggtitle(paste("State ", state.k, sep = ""))
    }
  }

  ### Assign nrow and ncol if not provided by user
  if (is.null(nrow)){
    nrow <- 2
  }
  if (is.null(ncol)){
    ncol <- base::ceiling(length(plots.list)/2)
  }

  ### Combine plots into single ggplot
  if (combine == TRUE){
    plots.list <- ggpubr::ggarrange(plotlist = plots.list, nrow = nrow, ncol = ncol, common.legend = TRUE)
  }

  ### Return output object
  return(plots.list)

}


#' Plots calibration scatter plots estimated using \code{\link{calib_mlr}}.
#'
#' @description
#' Plots calibration scatter plots for the transition probabilities of a multistate model
#' estimated using \code{\link{calib_mlr}}.
#'
#' @param x Object of class 'calib_mlr' generated from \code{\link{calib_mlr}}
#' @param ... Other
#' @param combine Whether to combine into one plot using ggarrange, or return as a list of individual plots
#' @param ncol Number of columns for combined calibration plot
#' @param nrow Number of rows for combined calibration plot
#' @param transparency.plot Degree of transparency for the calibration scatter plot
#' @param transparency.rug Degree of transparency for the density rug plot along each axis
#'
#' @returns If `combine = TRUE`, returns an object of classes `gg`, `ggplot`, and `ggarrange`,
#' as all ggplots have been combined into one object. If `combine = FALSE`, returns an object of
#' class `list`, each element containing an object of class `gg` and `ggplot`.
#'
#' @examples
#' # Using competing risks data out of initial state (see vignette: ... -in-competing-risk-setting).
#' # Estimate and plot MLR-IPCW calibration scatter plots for the predicted transition
#' # probabilities at time t = 1826, when predictions were made at time
#' # s = 0 in state j = 1. These predicted transition probabilities are stored in tp.cmprsk.j0.
#'
#' # To minimise example time we reduce the datasets to 150 individuals.
#' # Extract the predicted transition probabilities out of state j = 1 for first 150 individuals
#' tp.pred <- tp.cmprsk.j0 |>
#'  dplyr::filter(id %in% 1:150) |>
#'  dplyr::select(any_of(paste("pstate", 1:6, sep = "")))
#' # Reduce ebmtcal to first 150 individuals
#' ebmtcal <- ebmtcal |> dplyr::filter(id %in% 1:150)
#' # Reduce msebmtcal.cmprsk to first 150 individuals
#' msebmtcal.cmprsk <- msebmtcal.cmprsk |> dplyr::filter(id %in% 1:150)
#'
#' # Now estimate the observed event probabilities for each possible transition.
#' dat.calib.mlr <-
#' calib_mlr(data.mstate = msebmtcal.cmprsk,
#'  data.raw = ebmtcal,
#'  j=1,
#'  s=0,
#'  t = 1826,
#'  tp.pred = tp.pred,
#'  w.covs = c("year", "agecl", "proph", "match"),
#'  ps.int = 2,
#'  degree = 2)
#'
#'  # These are then plotted
#'  plot(dat.calib.mlr, combine = TRUE, nrow = 2, ncol = 3)
#'
#' @importFrom graphics plot
#' @export
plot.calib_mlr <- function(x, ..., combine = TRUE, ncol = NULL, nrow = NULL, transparency.plot = 0.25, transparency.rug = 0.1){

  ### Extract plot data and relevant metadata
  object.in <- x
  plot.data <- object.in[["plotdata"]]
  valid.transitions <- object.in[["metadata"]][["valid.transitions"]]

  ### Create list to store plots
  plots.list <- vector("list", length(valid.transitions))
  for (k in 1:length(valid.transitions)){

    ### Assign state of interest
    state.k <- valid.transitions[k]

    ### Assign plot data
    plot.data.k <- plot.data[[k]]

    ### Create the plots
    plots.list[[k]] <- ggplot2::ggplot(data = plot.data.k |> dplyr::arrange(pred) |>  dplyr::select(id, pred, obs)) +
      ggplot2::geom_point(ggplot2::aes(x = pred, y = obs), color = "red", alpha = transparency.plot, size = 0.5) +
      ggplot2::geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
      ggplot2::xlab("Predicted risk") + ggplot2::ylab("Observed risk") +
      ggplot2::xlim(c(0, max(plot.data.k$pred))) +
      ggplot2::ylim(c(min(plot.data.k$obs), max(plot.data.k$obs))) +
      ggplot2::geom_rug(ggplot2::aes(x = pred, y = obs), col = grDevices::rgb(1, 0, 0, alpha = .3), alpha = transparency.rug) +
      ggplot2::theme(legend.position = "none") +
      ggplot2::ggtitle(paste("State ", state.k, sep = ""))

  }

  ### Assign nrow and ncol if not provided by user
  if (is.null(nrow)){
    nrow <- 2
  }
  if (is.null(ncol)){
    ncol <- base::ceiling(length(plots.list)/2)
  }

  ### Combine plots into single ggplot
  if (combine == TRUE){
    plots.list <- ggpubr::ggarrange(plotlist = plots.list, nrow = nrow, ncol = ncol)
  }

  ### Return output object
  return(plots.list)

}


#' Plots calibration curves estimated using \code{\link{calib_pv}}.
#'
#' @description
#' Plots calibration curves for the transition probabilities of a multistate model
#' estimated using \code{\link{calib_pv}}.
#'
#' @param x Object of class 'calib_pseudo' generated from \code{\link{calib_pv}}.
#' @param ... Other
#' @param combine Whether to combine into one plot using ggarrange, or return as a list of individual plots
#' @param ncol Number of columns for combined calibration plot
#' @param nrow Number of rows for combined calibration plot
#' @param transparency.rug Degree of transparency for the density rug plot along each axis
#'
#' @returns If `combine = TRUE`, returns an object of classes `gg`, `ggplot`, and `ggarrange`,
#' as all ggplots have been combined into one object. If `combine = FALSE`, returns an object of
#' class `list`, each element containing an object of class `gg` and `ggplot`.
#'
#' @examples
#' # Using competing risks data out of initial state (see vignette: ... -in-competing-risk-setting).
#' # Estimate and plot pseudo-value calibration curves for the predicted transition
#' # probabilities at time t = 1826, when predictions were made at time
#' # s = 0 in state j = 1. These predicted transition probabilities are stored in tp.cmprsk.j0.
#'
#' # To minimise example time we reduce the datasets to 50 individuals.
#' # Extract the predicted transition probabilities out of state j = 1 for first 50 individuals
#' tp.pred <- tp.cmprsk.j0 |>
#'  dplyr::filter(id %in% 1:50) |>
#'  dplyr::select(any_of(paste("pstate", 1:6, sep = "")))
#' # Reduce ebmtcal to first 50 individuals
#' ebmtcal <- ebmtcal |> dplyr::filter(id %in% 1:50)
#' # Reduce msebmtcal.cmprsk to first 50 individuals
#' msebmtcal.cmprsk <- msebmtcal.cmprsk |> dplyr::filter(id %in% 1:50)
#'
#' # Now estimate the observed event probabilities for each possible transition.
#' dat.calib.pv <- calib_pv(data.mstate = msebmtcal.cmprsk,
#'   data.raw = ebmtcal,
#'   j = 1,
#'   s = 0,
#'   t = 1826,
#'   tp.pred = tp.pred,
#'   curve.type = "loess",
#'   loess.span = 1,
#'   loess.degree = 1)
#'
#'  # These are then plotted
#'  plot(dat.calib.pv, combine = TRUE, nrow = 2, ncol = 3)
#'
#' @importFrom graphics plot
#' @export
plot.calib_pv <- function(x, ..., combine = TRUE, ncol = NULL, nrow = NULL, transparency.rug = 0.1){

  ### Extract plot data and relevant metadata
  object.in <- x
  plot.data <- object.in[["plotdata"]]
  assessed.transitions <- object.in[["metadata"]][["assessed.transitions"]]
  CI <- object.in[["metadata"]][["CI"]]

  if (CI != FALSE){
    ### Create list to store plots
    plots.list <- vector("list", length(assessed.transitions))

    for (k in 1:length(assessed.transitions)){
      ### Assign plot data
      plot.data.k <- plot.data[[k]]

      ### Assign state of interest
      state.k <- assessed.transitions[k]

      ### Pivot longer to create data for ggplot and assign appropriate labels
      plot.data.k.longer <- tidyr::pivot_longer(plot.data.k, cols = c(obs, obs.upper, obs.lower), names_to = "line.group")
      plot.data.k.longer <- dplyr::mutate(plot.data.k.longer,
                                          line.group = base::factor(line.group),
                                          mapping = dplyr::case_when(line.group == "obs" ~ 1,
                                                                     line.group %in% c("obs.upper", "obs.lower") ~ 2),
                                          mapping = base::factor(mapping))

      levels(plot.data.k.longer$line.group) <- c("Calibration", "Upper", "Lower")
      levels(plot.data.k.longer$mapping) <- c("Calibration", "95% CI")

      ### Create the plots
      plots.list[[k]] <- ggplot2::ggplot(data = plot.data.k.longer |> dplyr::arrange(pred) |> dplyr::select(pred, line.group, value, mapping)) +
        ggplot2::geom_line(ggplot2::aes(x = pred, y = value, group = line.group, color = mapping)) +
        ggplot2::geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
        ggplot2::xlab("Predicted risk") + ggplot2::ylab("Observed risk") +
        ggplot2::xlim(c(0, max(plot.data.k.longer$pred))) +
        ggplot2::ylim(c(min(plot.data.k.longer$value), max(plot.data.k.longer$value))) +
        ggplot2::geom_rug(data = plot.data.k.longer |> dplyr::arrange(pred) |> dplyr::select(pred, line.group, value, mapping) |> subset(line.group == "Calibration"),
                          ggplot2::aes(x = pred, y = value), col = grDevices::rgb(1, 0, 0, alpha = transparency.rug)) +
        ggplot2::ggtitle(paste("State ", state.k, sep = ""))
    }
  } else if (CI == FALSE){
    ### Create list to store plots
    plots.list <- vector("list", length(assessed.transitions))
    for (k in 1:length(assessed.transitions)){

      ### Assign plot data
      plot.data.k <- plot.data[[k]]

      ### Assign state of interest
      state.k <- assessed.transitions[k]

      ### Create the plots
      plots.list[[k]] <- ggplot2::ggplot(data = plot.data.k |> dplyr::arrange(pred) |> dplyr::select(id, pred, obs)) +
        ggplot2::geom_line(ggplot2::aes(x = pred, y = obs)) +
        ggplot2::geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
        ggplot2::xlab("Predicted risk") + ggplot2::ylab("Observed risk") +
        ggplot2::xlim(c(0, max(plot.data.k$pred))) +
        ggplot2::ylim(c(min(plot.data.k$obs), max(plot.data.k$obs))) +
        ggplot2::geom_rug(ggplot2::aes(x = pred, y = obs), col = grDevices::rgb(1, 0, 0, alpha = transparency.rug)) +
        ggplot2::theme(legend.position = "none") +
        ggplot2::ggtitle(paste("State ", state.k, sep = ""))
    }
  }

  ### Assign nrow and ncol if not provided by user
  if (is.null(nrow)){
    nrow <- 2
  }
  if (is.null(ncol)){
    ncol <- base::ceiling(length(plots.list)/2)
  }

  ### Combine plots into single ggplot
  if (combine == TRUE){
    plots.list <- ggpubr::ggarrange(plotlist = plots.list, nrow = nrow, ncol = ncol, common.legend = TRUE)
  }

  ### Return output object
  return(plots.list)

}



