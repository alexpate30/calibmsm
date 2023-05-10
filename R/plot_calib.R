#' Plots calibration curves estimated using \code{\link{calc_calib_blr}}.
#'
#' @description
#' Plots calibration curves for the transition probabilities of a multistate model
#' using data estimated from \code{\link{calc_calib_blr}}.
#'
#' @param x Object of class 'calib_blr' generated from \code{\link{calc_calib_blr}}.
#' @param ... Other
#' @param combine Whether to combine into one plot using ggarrange, or return as a list of individual plots
#' @param ncol Number of columns for combined calibration plot
#' @param nrow Number of rows for combined calibration plot
#'
#' @returns If `combine = TRUE`, returns an object of classes `gg`, `ggplot`, and `ggarrange`,
#' as all ggplots have been combined into one object. If `combine = FALSE`, returns an object of
#' class `list`, each element containing an object of class `gg` and `ggplot`.
#'
#' @importFrom graphics plot
#' @export
plot.calib_blr <- function(x, ..., combine = TRUE, ncol = NULL, nrow = NULL){

  ### Extract plot data and relevant metadata
  object.in <- x
  plot.data <- object.in[["plotdata"]]
  valid.transitions <- object.in[["metadata"]][["valid.transitions"]]
  CI <- object.in[["metadata"]][["CI"]]

  if (CI != FALSE){
    ### Create list to store plots
    plots.list <- vector("list", length(valid.transitions))

    for (k in 1:length(valid.transitions)){
      ### Assign plot data
      plot.data.k <- plot.data[[k]]

      ### Assign state of interest
      state.k <- valid.transitions[k]

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
      plots.list[[k]] <- ggplot2::ggplot(data = plot.data.k.longer %>% dplyr::arrange(pred) %>% dplyr::select(pred, line.group, value, mapping)) +
        ggplot2::geom_line(ggplot2::aes(x = pred, y = value, group = line.group, color = mapping)) +
        ggplot2::geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
        ggplot2::xlab("Predicted risk") + ggplot2::ylab("Observed risk") +
        ggplot2::xlim(c(0, max(plot.data.k.longer$pred,
                                plot.data.k.longer$value))) +
        ggplot2::ylim(c(0, max(plot.data.k.longer$pred,
                                plot.data.k.longer$value))) +
        ggplot2::geom_rug(data = plot.data.k.longer %>% dplyr::arrange(pred) %>% dplyr::select(pred, line.group, value, mapping) %>% subset(line.group == "Calibration"),
                          ggplot2::aes(x = pred, y = value), col = grDevices::rgb(1, 0, 0, alpha = 0.3)) +
        ggplot2::ggtitle(paste("State ", state.k, sep = ""))
    }
  } else if (CI == FALSE){
    ### Create list to store plots
    plots.list <- vector("list", length(valid.transitions))
    for (k in 1:length(valid.transitions)){

      ### Assign plot data
      plot.data.k <- plot.data[[k]]

      ### Assign state of interest
      state.k <- valid.transitions[k]

      ### Create the plots
      plots.list[[k]] <- ggplot2::ggplot(data = plot.data.k %>% dplyr::arrange(pred) %>% dplyr::select(id, pred, obs)) +
        ggplot2::geom_line(ggplot2::aes(x = pred, y = obs)) +
        ggplot2::geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
        ggplot2::xlab("Predicted risk") + ggplot2::ylab("Observed risk") +
        ggplot2::xlim(c(0, max(plot.data.k$pred,
                                plot.data.k$obs))) +
        ggplot2::ylim(c(0, max(plot.data.k$pred,
                                plot.data.k$obs))) +
        ggplot2::geom_rug(ggplot2::aes(x = pred, y = obs), col = grDevices::rgb(1, 0, 0, alpha = 0.3)) +
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


#' Plots calibration scatter plots estimated using \code{\link{calc_calib_mlr}}.
#'
#' @description
#' Plots calibration scatter plots for the transition probabilities of a multistate model
#' using data estimated from \code{\link{calc_calib_mlr}}.
#'
#' @param x Object of class 'calib_mlr' generated from \code{\link{calc_calib_mlr}}
#' @param ... Other
#' @param combine Whether to combine into one plot using ggarrange, or return as a list of individual plots
#' @param ncol Number of columns for combined calibration plot
#' @param nrow Number of rows for combined calibration plot
#'
#' @returns If `combine = TRUE`, returns an object of classes `gg`, `ggplot`, and `ggarrange`,
#' as all ggplots have been combined into one object. If `combine = FALSE`, returns an object of
#' class `list`, each element containing an object of class `gg` and `ggplot`.
#'
#' @importFrom graphics plot
#' @export
plot.calib_mlr <- function(x, ..., combine = TRUE, ncol = NULL, nrow = NULL){

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
    plots.list[[k]] <- ggplot2::ggplot(data = plot.data.k %>% dplyr::arrange(pred) %>%  dplyr::select(id, pred, obs)) +
      ggplot2::geom_point(ggplot2::aes(x = pred, y = obs), color = "red", size = 0.5) +
      ggplot2::geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
      ggplot2::xlab("Predicted risk") + ggplot2::ylab("Observed risk") +
      ggplot2::xlim(c(0, max(plot.data.k$pred,
                    plot.data.k$obs))) +
      ggplot2::ylim(c(0, max(plot.data.k$pred,
                    plot.data.k$obs))) +
      ggplot2::geom_rug(ggplot2::aes(x = pred, y = obs), col = grDevices::rgb(1, 0, 0, alpha = .3)) +
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



