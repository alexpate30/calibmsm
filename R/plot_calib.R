#' Plots calibration curves estimated using \code{\link{calib_msm}}.
#'
#' @description
#' Plots calibration curves for the transition probabilities of a multistate model
#' estimated using BLR-IPCW and pseudo-value approaches.
#'
#' @param x Object of class 'calib_msm' generated from \code{\link{calib_msm}}.
#' @param ... Other
#' @param combine Whether to combine into one plot using ggarrange, or return as a list of individual plots
#' @param ncol Number of columns for combined calibration plot
#' @param nrow Number of rows for combined calibration plot
#' @param size_text Size of text in plot
#' @param size_line Size of line plots
#' @param marg_density Whether to produce marginal density plots TRUE/FALSE
#' @param marg_density_size Size of the main plot relative to the density plots (see \code{\link[ggExtra]{ggMarginal}})
#' @param marg_density_type What type of marginal plot to show (see \code{\link[ggExtra]{ggMarginal}})
#' @param marg_rug Whether to produce marginal rug plots TRUE/FALSE
#' @param marg_rug_transparency Degree of transparency for the density rug plot along each axis
#' @param titles_include Whether to include titles for each individual calibration plots
#' @param titles Vector of titles for the calibration plots. Defaults to "State k" for each plot.
#' @param axis_titles_x Position of plots for which to include title on x-axis
#' @param axis_titles_text_x x-axis title
#' @param axis_titles_y Position of plots for which to include title on y-axis
#' @param axis_titles_text_y y-axis title
#' @param legend_include Whether to produce a legend
#' @param legend_seperate = Whether to include legend in plot (FALSE) or as a seperate object (TRUE)
#' @param legend_title Title of legend
#' @param legend_position Position of legend
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
#'  # These are then plotted
#'  plot(dat_calib, combine = TRUE, nrow = 2, ncol = 3)
#'
#' @importFrom graphics plot
#' @export
plot.calib_msm <- function(x, ..., combine = TRUE, ncol = NULL, nrow = NULL, size_line = 0.5, size_text = 12,
                           marg_density = TRUE, marg_density_size = 5, marg_density_type = "density",
                           marg_rug = FALSE, marg_rug_transparency = 0.1,
                           titles_include = TRUE, titles = NULL,
                           axis_titles_x = NULL, axis_titles_text_x = "Predicted risk",
                           axis_titles_y = NULL, axis_titles_text_y = "Predicted-observed risk",
                           legend_include = TRUE, legend_seperate = FALSE, legend_title = NULL, legend_position = "bottom"){

  # x <- dat_calib_blr
  # str(x)
  #
  # ncol = 5
  # nrow = 1
  # marg_density = FALSE
  # marg_density_size = 5
  # marg_density_type = "density"
  # marg_rug = FALSE
  # marg_rug_transparency = 0.1
  # titles_include = TRUE
  # legend_seperate = TRUE
  # legend_title = NULL
  # axis_titles_x = NULL
  # axis_titles_text_x = "Predicted risk"
  # axis_titles_y = NULL
  # axis_titles_text_y = "Observed risk"
  # size = 12

  ### Extract plot data and relevant metadata
  object_in <- x
  plot_data <- object_in[["plotdata"]]
  assessed_transitions <- object_in[["metadata"]][["assessed_transitions"]]
  CI <- object_in[["metadata"]][["CI"]]

  ### Create list to store plots
  plots_list <- vector("list", length(assessed_transitions))

  for (k in 1:length(assessed_transitions)){

    if (CI != FALSE){

      ### Assign plot data
      plot_data_k <- plot_data[[k]]

      ### Assign state of interest
      state_k <- assessed_transitions[k]

      ### Pivot longer to create data for ggplot and assign appropriate labels
      plot_data_k_longer <- tidyr::pivot_longer(plot_data_k, cols = c(obs, obs_upper, obs_lower), names_to = "line_group")
      plot_data_k_longer <- dplyr::mutate(plot_data_k_longer,
                                          line_group = base::factor(line_group),
                                          mapping = dplyr::case_when(line_group == "obs" ~ 1,
                                                                     line_group %in% c("obs_upper", "obs_lower") ~ 2),
                                          mapping = base::factor(mapping))

      levels(plot_data_k_longer$line_group) <- c("Calibration", "Upper", "Lower")
      levels(plot_data_k_longer$mapping) <- c("Calibration", "95% CI")

      ### Create the plots
      plots_list[[k]] <- ggplot2::ggplot(data = plot_data_k_longer |> dplyr::arrange(pred) |> dplyr::select(pred, line_group, value, mapping)) +
        ggplot2::geom_line(ggplot2::aes(x = pred, y = value, group = line_group, color = mapping), linewidth = size_line) +
        ggplot2::geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
        ggplot2::xlim(c(min(min(plot_data_k_longer$value), min(plot_data_k_longer$pred)),
                        max(max(plot_data_k_longer$value), max(plot_data_k_longer$pred)))) +
        ggplot2::ylim(c(min(min(plot_data_k_longer$value), min(plot_data_k_longer$pred)),
                        max(max(plot_data_k_longer$value), max(plot_data_k_longer$pred)))) +
        ggplot2::labs(x = NULL, y = NULL)

    } else if (CI == FALSE){

      ### Assign plot data
      plot_data_k <- plot_data[[k]]

      ### Assign state of interest
      state_k <- assessed_transitions[k]

      ### Create the plots
      plots_list[[k]] <- ggplot2::ggplot(data = plot_data_k |> dplyr::arrange(pred) |> dplyr::select(id, pred, obs) |> dplyr::rename(value = obs)) +
        ggplot2::geom_line(ggplot2::aes(x = pred, y = value), colour = "red", linewidth = size_line) +
        ggplot2::geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
        ggplot2::xlim(c(min(min(plot_data_k$obs), min(plot_data_k$pred)),
                        max(max(plot_data_k$obs), max(plot_data_k$pred)))) +
        ggplot2::ylim(c(min(min(plot_data_k$obs), min(plot_data_k$pred)),
                        max(max(plot_data_k$obs), max(plot_data_k$pred)))) +
        ggplot2::labs(x = NULL, y = NULL)


    }

    ### Add text size
    plots_list[[k]] <- plots_list[[k]] +
      ggplot2::theme(text = ggplot2::element_text(size = size_text),
                     legend.text = ggplot2::element_text(size = size_text))

    ### Add ggtitles if specified
    if (titles_include == TRUE){
      if (is.null(titles)){
        plots_list[[k]] <- plots_list[[k]] + ggplot2::ggtitle(paste("State ", state_k, sep = ""))
      } else {
        plots_list[[k]] <- plots_list[[k]] + ggplot2::ggtitle(titles[k])
      }
    }

    ### Add axis titles
    if (is.null(axis_titles_x)){
      plots_list[[k]] <- plots_list[[k]] + ggplot2::xlab(axis_titles_text_x)
    } else if (!is.null(axis_titles_x)){
      if (k %in% axis_titles_x){
        plots_list[[k]] <- plots_list[[k]] + ggplot2::xlab(axis_titles_text_x)
      }
    }

    if (is.null(axis_titles_y)){
      if ("calib_pv" %in% class(object_in)){axis_titles_text_y <- "Pseudo-observed risk"}
      plots_list[[k]] <- plots_list[[k]] + ggplot2::ylab(axis_titles_text_y)
    } else if (!is.null(axis_titles_y)){
      if ("calib_pv" %in% class(object_in)){axis_titles_text_y <- "Pseudo-observed risk"}
      if (k %in% axis_titles_y){
        plots_list[[k]] <- plots_list[[k]] + ggplot2::ylab(axis_titles_text_y)
      }
    }

    ### Add legend title if specified
    if (is.null(legend_title)){
      plots_list[[k]] <- plots_list[[k]] + ggplot2::theme(legend.title = ggplot2::element_blank())
    } else if (!is.null(legend_title)){
      plots_list[[k]] <- plots_list[[k]] + ggplot2::theme(legend.title = ggplot2::element_text(size = size_text, face = "bold")) +
        ggplot2::guides(color = ggplot2::guide_legend(title = legend_title), lty = ggplot2::guide_legend(title = legend_title))
    }

    ## Save legend
    if (k == 1){
      legend_save <- ggpubr::get_legend(plots_list[[k]], position = legend_position)
    }

    ## Remove legend if requested to do so.
    if (legend_include == FALSE){
      plots_list[[k]] <- plots_list[[k]] +
        ggplot2::theme(legend.position = "none")
    }

    ### If marginal density plot has been requested add an invisible scatter plot on which to base this
    ### We also remove legend forcefully. When adding marginal density plots, the legend must be added using an arrangeGrob.
    if (marg_density == TRUE){
      ### IF CI != FALSE want to only extract calibration line for the density plot
      if (CI != FALSE){
        ## Add scatter
        plots_list[[k]] <- plots_list[[k]] +
          ggplot2::geom_point(data = plot_data_k_longer |> subset(mapping == "Calibration"),
                              ggplot2::aes(x = pred, y = value),
                              col = grDevices::rgb(0, 0, 0, alpha = 0)) +
          ## Remove legend
          ggplot2::theme(legend.position = "none")
      } else {
        ## Add scatter
        plots_list[[k]] <- plots_list[[k]] +
          ggplot2::geom_point(ggplot2::aes(x = pred, y = value),
                              col = grDevices::rgb(0, 0, 0, alpha = 0)) +
          ## Remove legend
          ggplot2::theme(legend.position = "none")
      }

      ## Add ggMarginal
      plots_list[[k]] <- ggExtra::ggMarginal(plots_list[[k]], margins = "x", size = marg_density_size, type = marg_density_type, colour = "red")

    } else if (marg_rug == TRUE){

      ## Add the marginal rug plot
      if (CI != FALSE){
        plots_list[[k]] <- plots_list[[k]] +
          ggplot2::geom_rug(data = plot_data_k_longer |> subset(mapping == "Calibration"),
                            ggplot2::aes(x = pred, y = value), col = grDevices::rgb(1, 0, 0, alpha = marg_rug_transparency))
      } else {
        plots_list[[k]] <- plots_list[[k]] +
          ggplot2::geom_rug(ggplot2::aes(x = pred, y = value), col = grDevices::rgb(1, 0, 0, alpha = marg_rug_transparency))
      }

    }

  }

  ### Assign nrow and ncol if not provided by user
  if (is.null(nrow)){
    nrow <- 2
  }
  if (is.null(ncol)){
    ncol <- base::ceiling(length(plots_list)/2)
  }

  ### Combine plots into single ggplot
  if (combine == TRUE){
    if (marg_density == FALSE){
      if (legend_include == TRUE){
        if (legend_seperate == FALSE){
          ## Combine with common legend
          plots_list <- ggpubr::ggarrange(plotlist = plots_list, nrow = nrow, ncol = ncol, common.legend = TRUE, legend = legend_position)
        } else {
          ## Combine without legend
          plots_list <- ggpubr::ggarrange(plotlist = plots_list, nrow = nrow, ncol = ncol, legend = "none")
          ## Add legend as seperate list element
          plots_list <- list("plots" = plots_list, "legend" = legend_save)
        }
      } else if (legend_include == FALSE){
        plots_list <- ggpubr::ggarrange(plotlist = plots_list, nrow = nrow, ncol = ncol, legend = "none")
      }
    } else if (marg_density == TRUE){
      plots_list <- gridExtra::arrangeGrob(grobs = plots_list,
                                           layout_matrix = base::matrix(base::seq_len(nrow*ncol),
                                                                        nrow = nrow,
                                                                        ncol = ncol,
                                                                        byrow = TRUE),
                                           top = NULL)
      ### Marginal density plots require legend to be added manually, because otherwise you get the scatter plot which ggMarginal relies on.
      ### We only add legend if CI != FALSE, as there is no legend when CI == FALSE
      if (legend_include == TRUE & !isFALSE(CI)){
        if (legend_seperate == FALSE){
          plots_list <- gridExtra::arrangeGrob(plots_list, legend_save, nrow = 2, heights = c(15, 1))
        } else {
          plots_list <- list("plots" = plots_list, "legend" = legend_save)
        }
      }
    }
  } else if (combine == FALSE){
    ### If combine == FALSE, but marginal_density == TRUE, need to manually legend for each calibration plot legends are requested
    if (marg_density == TRUE & legend_include == TRUE & !isFALSE(CI)){
      if (legend_seperate == FALSE){
        for (k in 1:length(assessed_transitions)){
          plots_list[[k]] <- gridExtra::arrangeGrob(plots_list[[k]], legend_save, nrow = 2, heights = c(15, 1))
        }
      } else if (legend_seperate == TRUE){
        plots_list <- list("plots" = plots_list, "legend" = legend_save)
      }
    }
  }

  ### Return output object
  return(plots_list)

}


#' Plots calibration scatter plots for objects of class `calib_mlr` estimated using
#' using \code{\link{calib_msm}}.
#'
#' @description
#' Plots calibration scatter plots for the transition probabilities of a multistate model
#' estimated using the MLR-IPCW approach.
#'
#' @param x Object of class `calib_mlr` generated from \code{\link{calib_msm}}
#' @param ... Other
#' @param combine Whether to combine into one plot using ggarrange, or return as a list of individual plots
#' @param ncol Number of columns for combined calibration plot
#' @param nrow Number of rows for combined calibration plot
#' @param size_point Size of points in scatter plot
#' @param size_text Size of text in plot
#' @param transparency_plot Degree of transparency for points in the calibration scatter plot
#' @param marg_density Whether to produce marginal density plots TRUE/FALSE
#' @param marg_density_size Size of the main plot relative to the density plots (see \code{\link[ggExtra]{ggMarginal}})
#' @param marg_density_type What type of marginal plot to show (see \code{\link[ggExtra]{ggMarginal}})
#' @param marg_rug Whether to produce marginal rug plots TRUE/FALSE
#' @param marg_rug_transparency Degree of transparency for the density rug plot along each axis
#' @param titles_include Whether to include titles for each individual calibration plots
#' @param titles Vector of titles for the calibration plots_ Defaults to "State k" for each plot_
#' @param axis_titles_x Position of plots for which to include title on x-axis
#' @param axis_titles_text_x x-axis title
#' @param axis_titles_y Position of plots for which to include title on y-axis
#' @param axis_titles_text_y y-axis title
#'
#' @returns If `combine = TRUE`, returns an object of classes `gg`, `ggplot`, and `ggarrange`,
#' as all ggplots have been combined into one object. If `combine = FALSE`, returns an object of
#' class `list`, each element containing an object of class `gg` and `ggplot`.
#'
#' @examples
#' # Using competing risks data out of initial state (see vignette: ... -in-competing-risk-setting).
#' # Estimate and plot MLR-IPCW calibration scatter plots for the predicted transition
#' # probabilities at time t = 1826, when predictions were made at time
#' # s = 0 in state j = 1. These predicted transition probabilities are stored in tp_cmprsk_j0.
#'
#' # To minimise example time we reduce the datasets to 150 individuals.
#' # Extract the predicted transition probabilities out of state j = 1 for first 150 individuals
#' tp_pred <- tp_cmprsk_j0 |>
#'  dplyr::filter(id %in% 1:150) |>
#'  dplyr::select(any_of(paste("pstate", 1:6, sep = "")))
#' # Reduce ebmtcal to first 150 individuals
#' ebmtcal <- ebmtcal |> dplyr::filter(id %in% 1:150)
#' # Reduce msebmtcal_cmprsk to first 150 individuals
#' msebmtcal_cmprsk <- msebmtcal_cmprsk |> dplyr::filter(id %in% 1:150)
#'
#' # Now estimate the observed event probabilities for each possible transition.
#' dat_calib <-
#' calib_msm(data_ms = msebmtcal_cmprsk,
#'  data_raw = ebmtcal,
#'  j=1,
#'  s=0,
#'  t = 1826,
#'  tp_pred = tp_pred,
#'  calib_type = "mlr",
#'  w_covs = c("year", "agecl", "proph", "match"),
#'  mlr_ps_int = 2,
#'  mlr_degree = 2)
#'
#'  # These are then plotted
#'  plot(dat_calib, combine = TRUE, nrow = 2, ncol = 3)
#'
#' @importFrom graphics plot
#' @export
plot.calib_mlr <- function(x, ..., combine = TRUE, ncol = NULL, nrow = NULL, size_point = 0.5, size_text = 12, transparency_plot = 0.25,
                           marg_density = FALSE, marg_density_size = 5, marg_density_type = "density",
                           marg_rug = FALSE, marg_rug_transparency = 0.1,
                           titles_include = TRUE, titles = NULL,
                           axis_titles_x = NULL, axis_titles_text_x = "Predicted risk",
                           axis_titles_y = NULL, axis_titles_text_y = "Predicted-observed risk"){

  ### Extract plot data and relevant metadata
  object_in <- x
  plot_data <- object_in[["plotdata"]]
  valid_transitions <- object_in[["metadata"]][["valid_transitions"]]

  ### Create list to store plots
  plots_list <- vector("list", length(valid_transitions))
  for (k in 1:length(valid_transitions)){

    ### Assign state of interest
    state_k <- valid_transitions[k]

    ### Assign plot data
    plot_data_k <- plot_data[[k]]

    ### Create the plots
    plots_list[[k]] <- ggplot2::ggplot(data = plot_data_k |> dplyr::arrange(pred) |>  dplyr::select(id, pred, obs)) +
      ggplot2::geom_point(ggplot2::aes(x = pred, y = obs), color = "red", alpha = transparency_plot, size = size_point) +
      ggplot2::geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
      ggplot2::xlim(c(0, max(plot_data_k$pred))) +
      ggplot2::ylim(c(min(plot_data_k$obs), max(plot_data_k$obs))) +
      ggplot2::theme(legend.position = "none") +
      ggplot2::labs(x = NULL, y = NULL) +
      ggplot2::theme(text = ggplot2::element_text(size = size_text),
                     legend.text = ggplot2::element_text(size = size_text))

    ### Add ggtitles if specified
    if (titles_include == TRUE){
      if (is.null(titles)){
        plots_list[[k]] <- plots_list[[k]] + ggplot2::ggtitle(paste("State ", state_k, sep = ""))
      } else {
        plots_list[[k]] <- plots_list[[k]] + ggplot2::ggtitle(titles[k])
      }
    }

    ### Add axis titles
    if (is.null(axis_titles_x)){
      plots_list[[k]] <- plots_list[[k]] + ggplot2::xlab(axis_titles_text_x)
    } else if (!is.null(axis_titles_x)){
      if (k %in% axis_titles_x){
        plots_list[[k]] <- plots_list[[k]] + ggplot2::xlab(axis_titles_text_x)
      }
    }

    if (is.null(axis_titles_y)){
      plots_list[[k]] <- plots_list[[k]] + ggplot2::ylab(axis_titles_text_y)
    } else if (!is.null(axis_titles_y)){
      if (k %in% axis_titles_y){
        plots_list[[k]] <- plots_list[[k]] + ggplot2::ylab(axis_titles_text_y)
      }
    }

    ### If marginal density plot has been requested add density plot
    if (marg_density == TRUE){
      plots_list[[k]] <- plots_list[[k]] +
        ## Add a geom_point object of the line and set to invisible (scatter plot required for marginal density using ggMarginal)
        ## Subset to ignore the confidence intervals when doing the density plots
        ggplot2::geom_point(ggplot2::aes(x = pred, y = obs), col = grDevices::rgb(0, 0, 0, alpha = 0))

      ## Add ggMarginal
      plots_list[[k]] <- ggExtra::ggMarginal(plots_list[[k]], margins = "x", size = marg_density_size, type = marg_density_type, colour = "red")
      ### If marginal rug plot has been requested
    } else if (marg_rug == TRUE){
      plots_list[[k]] <- plots_list[[k]] +
        ggplot2::geom_rug(ggplot2::aes(x = pred, y = obs), col = grDevices::rgb(1, 0, 0, alpha = marg_rug_transparency))
    }

  }

  ### Assign nrow and ncol if not provided by user
  if (is.null(nrow)){
    nrow <- 2
  }
  if (is.null(ncol)){
    ncol <- base::ceiling(length(plots_list)/2)
  }

  ### Combine plots into single ggplot
  if (combine == TRUE){
    if (marg_density == FALSE){
      plots_list <- ggpubr::ggarrange(plotlist = plots_list, nrow = nrow, ncol = ncol, common.legend = TRUE)
    } else if (marg_density == TRUE){
      plots_list <- gridExtra::arrangeGrob(grobs = plots_list,
                                           layout_matrix = base::matrix(base::seq_len(nrow*ncol),
                                                                        nrow = nrow,
                                                                        ncol = ncol,
                                                                        byrow = TRUE),
                                           top = NULL)
    }
  }

  ### Return output object
  return(plots_list)

}
