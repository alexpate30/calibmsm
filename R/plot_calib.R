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
#' @param size.text Size of text in plot
#' @param size.line Size of line plots
#' @param marg.density Whether to produce marginal density plots TRUE/FALSE
#' @param marg.density.size Size of the main plot relative to the density plots (see \code{\link[ggExtra]{ggMarginal}})
#' @param marg.density.type What type of marginal plot to show (see \code{\link[ggExtra]{ggMarginal}})
#' @param marg.rug Whether to produce marginal rug plots TRUE/FALSE
#' @param marg.rug.transparency Degree of transparency for the density rug plot along each axis
#' @param titles.include Whether to include titles for each individual calibration plots
#' @param titles Vector of titles for the calibration plots. Defaults to "State k" for each plot.
#' @param axis.titles.x Position of plots for which to include title on x-axis
#' @param axis.titles.text.x x-axis title
#' @param axis.titles.y Position of plots for which to include title on y-axis
#' @param axis.titles.text.y y-axis title
#' @param legend.include Whether to produce a legend
#' @param legend.seperate = Whether to include legend in plot (FALSE) or as a seperate object (TRUE)
#' @param legend.title Title of legend
#' @param legend.position Position of legend
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
#' dat.calib <-
#' calib_msm(data.mstate = msebmtcal,
#'  data.raw = ebmtcal,
#'  j=1,
#'  s=0,
#'  t = 1826,
#'  tp.pred = tp.pred,
#'  w.covs = c("year", "agecl", "proph", "match"))
#'
#'  # These are then plotted
#'  plot(dat.calib, combine = TRUE, nrow = 2, ncol = 3)
#'
#' @importFrom graphics plot
#' @export
plot.calib_msm <- function(x, ..., combine = TRUE, ncol = NULL, nrow = NULL, size.line = 0.5, size.text = 12,
                           marg.density = FALSE, marg.density.size = 5, marg.density.type = "density",
                           marg.rug = FALSE, marg.rug.transparency = 0.1,
                           titles.include = TRUE, titles = NULL,
                           axis.titles.x = NULL, axis.titles.text.x = "Predicted risk",
                           axis.titles.y = NULL, axis.titles.text.y = "Observed risk",
                           legend.include = TRUE, legend.seperate = FALSE, legend.title = NULL, legend.position = "bottom"){

  # x <- dat.calib.blr
  # str(x)
  #
  # ncol = 5
  # nrow = 1
  # marg.density = FALSE
  # marg.density.size = 5
  # marg.density.type = "density"
  # marg.rug = FALSE
  # marg.rug.transparency = 0.1
  # titles.include = TRUE
  # legend.seperate = TRUE
  # legend.title = NULL
  # axis.titles.x = NULL
  # axis.titles.text.x = "Predicted risk"
  # axis.titles.y = NULL
  # axis.titles.text.y = "Observed risk"
  # size = 12


  ### Extract plot data and relevant metadata
  object.in <- x
  plot.data <- object.in[["plotdata"]]
  assessed.transitions <- object.in[["metadata"]][["assessed.transitions"]]
  CI <- object.in[["metadata"]][["CI"]]

  ### Create list to store plots
  plots.list <- vector("list", length(assessed.transitions))

  for (k in 1:length(assessed.transitions)){

    if (CI != FALSE){

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
        ggplot2::geom_line(ggplot2::aes(x = pred, y = value, group = line.group, color = mapping), linewidth = size.line) +
        ggplot2::geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
        ggplot2::xlim(c(min(min(plot.data.k.longer$value), min(plot.data.k.longer$pred)),
                        max(max(plot.data.k.longer$value), max(plot.data.k.longer$pred)))) +
        ggplot2::ylim(c(min(min(plot.data.k.longer$value), min(plot.data.k.longer$pred)),
                        max(max(plot.data.k.longer$value), max(plot.data.k.longer$pred)))) +
        ggplot2::labs(x = NULL, y = NULL)

    } else if (CI == FALSE){

      ### Assign plot data
      plot.data.k <- plot.data[[k]]

      ### Assign state of interest
      state.k <- assessed.transitions[k]

      ### Create the plots
      plots.list[[k]] <- ggplot2::ggplot(data = plot.data.k |> dplyr::arrange(pred) |> dplyr::select(id, pred, obs) |> dplyr::rename(value = obs)) +
        ggplot2::geom_line(ggplot2::aes(x = pred, y = value), colour = "red", linewidth = size.line) +
        ggplot2::geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
        ggplot2::xlim(c(min(min(plot.data.k$obs), min(plot.data.k$pred)),
                        max(max(plot.data.k$obs), max(plot.data.k$pred)))) +
        ggplot2::ylim(c(min(min(plot.data.k$obs), min(plot.data.k$pred)),
                        max(max(plot.data.k$obs), max(plot.data.k$pred)))) +
        ggplot2::labs(x = NULL, y = NULL)


    }

    ### Add text size
    plots.list[[k]] <- plots.list[[k]] +
      ggplot2::theme(text = ggplot2::element_text(size = size.text),
                     legend.text = ggplot2::element_text(size = size.text))

    ### Add ggtitles if specified
    if (titles.include == TRUE){
      if (is.null(titles)){
        plots.list[[k]] <- plots.list[[k]] + ggplot2::ggtitle(paste("State ", state.k, sep = ""))
      } else {
        plots.list[[k]] <- plots.list[[k]] + ggplot2::ggtitle(titles[k])
      }
    }

    ### Add axis titles
    if (is.null(axis.titles.x)){
      plots.list[[k]] <- plots.list[[k]] + ggplot2::xlab(axis.titles.text.x)
    } else if (!is.null(axis.titles.x)){
      if (k %in% axis.titles.x){
        plots.list[[k]] <- plots.list[[k]] + ggplot2::xlab(axis.titles.text.x)
      }
    }

    if (is.null(axis.titles.y)){
      plots.list[[k]] <- plots.list[[k]] + ggplot2::ylab(axis.titles.text.y)
    } else if (!is.null(axis.titles.y)){
      if (k %in% axis.titles.y){
        plots.list[[k]] <- plots.list[[k]] + ggplot2::ylab(axis.titles.text.y)
      }
    }

    ### Add legend title if specified
    if (is.null(legend.title)){
      plots.list[[k]] <- plots.list[[k]] + ggplot2::theme(legend.title = ggplot2::element_blank())
    } else if (!is.null(legend.title)){
      plots.list[[k]] <- plots.list[[k]] + ggplot2::theme(legend.title = ggplot2::element_text(size = size.text, face = "bold")) +
        ggplot2::guides(color = ggplot2::guide_legend(title = legend.title), lty = ggplot2::guide_legend(title = legend.title))
    }

    ## Save legend
    if (k == 1){
      legend.save <- ggpubr::get_legend(plots.list[[k]], position = legend.position)
    }

    ## Remove legend if requested to do so.
    if (legend.include == FALSE){
      plots.list[[k]] <- plots.list[[k]] +
        ggplot2::theme(legend.position = "none")
    }

    ### If marginal density plot has been requested add an invisible scatter plot on which to base this
    ### We also remove legend forcefully. When adding marginal density plots, the legend must be added using an arrangeGrob.
    if (marg.density == TRUE){
      ### IF CI != FALSE want to only extract calibration line for the density plot
      if (CI != FALSE){
        ## Add scatter
        plots.list[[k]] <- plots.list[[k]] +
          ggplot2::geom_point(data = plot.data.k.longer |> subset(mapping == "Calibration"),
                              ggplot2::aes(x = pred, y = value),
                              col = grDevices::rgb(0, 0, 0, alpha = 0)) +
          ## Remove legend
          ggplot2::theme(legend.position = "none")
      } else {
        ## Add scatter
        plots.list[[k]] <- plots.list[[k]] +
          ggplot2::geom_point(ggplot2::aes(x = pred, y = value),
                              col = grDevices::rgb(0, 0, 0, alpha = 0)) +
          ## Remove legend
          ggplot2::theme(legend.position = "none")
      }

      ## Add ggMarginal
      plots.list[[k]] <- ggExtra::ggMarginal(plots.list[[k]], margins = "x", size = marg.density.size, type = marg.density.type, colour = "red")

    } else if (marg.rug == TRUE){

      ## Add the marginal rug plot
      if (CI != FALSE){
        plots.list[[k]] <- plots.list[[k]] +
          ggplot2::geom_rug(data = plot.data.k.longer |> subset(mapping == "Calibration"),
                            ggplot2::aes(x = pred, y = value), col = grDevices::rgb(1, 0, 0, alpha = marg.rug.transparency))
      } else {
        plots.list[[k]] <- plots.list[[k]] +
          ggplot2::geom_rug(ggplot2::aes(x = pred, y = value), col = grDevices::rgb(1, 0, 0, alpha = marg.rug.transparency))
      }

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
    if (marg.density == FALSE){
      if (legend.include == TRUE){
        if (legend.seperate == FALSE){
          ## Combine with common legend
          plots.list <- ggpubr::ggarrange(plotlist = plots.list, nrow = nrow, ncol = ncol, common.legend = TRUE, legend = legend.position)
        } else {
          ## Combine without legend
          plots.list <- ggpubr::ggarrange(plotlist = plots.list, nrow = nrow, ncol = ncol, legend = "none")
          ## Add legend as seperate list element
          plots.list <- list("plots" = plots.list, "legend" = legend.save)
        }
      } else if (legend.include == FALSE){
        plots.list <- ggpubr::ggarrange(plotlist = plots.list, nrow = nrow, ncol = ncol, legend = "none")
      }
    } else if (marg.density == TRUE){
      plots.list <- gridExtra::arrangeGrob(grobs = plots.list,
                                           layout_matrix = base::matrix(base::seq_len(nrow*ncol),
                                                                        nrow = nrow,
                                                                        ncol = ncol,
                                                                        byrow = TRUE),
                                           top = NULL)
      ### Marginal density plots require legend to be added manually, because otherwise you get the scatter plot which ggMarginal relies on.
      ### We only add legend if CI != FALSE, as there is no legend when CI == FALSE
      if (legend.include == TRUE & !isFALSE(CI)){
        if (legend.seperate == FALSE){
          plots.list <- gridExtra::arrangeGrob(plots.list, legend.save, nrow = 2, heights = c(15, 1))
        } else {
          plots.list <- list("plots" = plots.list, "legend" = legend.save)
        }
      }
    }
  } else if (combine == FALSE){
    ### If combine == FALSE, but marginal.density == TRUE, need to manually legend for each calibration plot legends are requested
    if (marg.density == TRUE & legend.include == TRUE & !isFALSE(CI)){
      if (legend.seperate == FALSE){
        for (k in 1:length(assessed.transitions)){
          plots.list[[k]] <- gridExtra::arrangeGrob(plots.list[[k]], legend.save, nrow = 2, heights = c(15, 1))
        }
      } else if (legend.seperate == TRUE){
        plots.list <- list("plots" = plots.list, "legend" = legend.save)
      }
    }
  }

  ### Return output object
  return(plots.list)

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
#' @param size.point Size of points in scatter plot
#' @param size.text Size of text in plot
#' @param transparency.plot Degree of transparency for points in the calibration scatter plot
#' @param marg.density Whether to produce marginal density plots TRUE/FALSE
#' @param marg.density.size Size of the main plot relative to the density plots (see \code{\link[ggExtra]{ggMarginal}})
#' @param marg.density.type What type of marginal plot to show (see \code{\link[ggExtra]{ggMarginal}})
#' @param marg.rug Whether to produce marginal rug plots TRUE/FALSE
#' @param marg.rug.transparency Degree of transparency for the density rug plot along each axis
#' @param titles.include Whether to include titles for each individual calibration plots
#' @param titles Vector of titles for the calibration plots. Defaults to "State k" for each plot.
#' @param axis.titles.x Position of plots for which to include title on x-axis
#' @param axis.titles.text.x x-axis title
#' @param axis.titles.y Position of plots for which to include title on y-axis
#' @param axis.titles.text.y y-axis title
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
#' dat.calib <-
#' calib_msm(data.mstate = msebmtcal.cmprsk,
#'  data.raw = ebmtcal,
#'  j=1,
#'  s=0,
#'  t = 1826,
#'  tp.pred = tp.pred,
#'  calib.type = "mlr",
#'  w.covs = c("year", "agecl", "proph", "match"),
#'  mlr.ps.int = 2,
#'  mlr.degree = 2)
#'
#'  # These are then plotted
#'  plot(dat.calib, combine = TRUE, nrow = 2, ncol = 3)
#'
#' @importFrom graphics plot
#' @export
plot.calib_mlr <- function(x, ..., combine = TRUE, ncol = NULL, nrow = NULL, size.point = 0.5, size.text = 12, transparency.plot = 0.25,
                           marg.density = FALSE, marg.density.size = 5, marg.density.type = "density",
                           marg.rug = FALSE, marg.rug.transparency = 0.1,
                           titles.include = TRUE, titles = NULL,
                           axis.titles.x = NULL, axis.titles.text.x = "Predicted risk",
                           axis.titles.y = NULL, axis.titles.text.y = "Observed risk"){

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
      ggplot2::geom_point(ggplot2::aes(x = pred, y = obs), color = "red", alpha = transparency.plot, size = size.point) +
      ggplot2::geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
      ggplot2::xlim(c(0, max(plot.data.k$pred))) +
      ggplot2::ylim(c(min(plot.data.k$obs), max(plot.data.k$obs))) +
      ggplot2::theme(legend.position = "none") +
      ggplot2::labs(x = NULL, y = NULL) +
      ggplot2::theme(text = ggplot2::element_text(size = size.text),
                     legend.text = ggplot2::element_text(size = size.text))

    ### Add ggtitles if specified
    if (titles.include == TRUE){
      if (is.null(titles)){
        plots.list[[k]] <- plots.list[[k]] + ggplot2::ggtitle(paste("State ", state.k, sep = ""))
      } else {
        plots.list[[k]] <- plots.list[[k]] + ggplot2::ggtitle(titles[k])
      }
    }

    ### Add axis titles
    if (is.null(axis.titles.x)){
      plots.list[[k]] <- plots.list[[k]] + ggplot2::xlab(axis.titles.text.x)
    } else if (!is.null(axis.titles.x)){
      if (k %in% axis.titles.x){
        plots.list[[k]] <- plots.list[[k]] + ggplot2::xlab(axis.titles.text.x)
      }
    }

    if (is.null(axis.titles.y)){
      plots.list[[k]] <- plots.list[[k]] + ggplot2::ylab(axis.titles.text.y)
    } else if (!is.null(axis.titles.y)){
      if (k %in% axis.titles.y){
        plots.list[[k]] <- plots.list[[k]] + ggplot2::ylab(axis.titles.text.y)
      }
    }

    ### If marginal density plot has been requested add density plot
    if (marg.density == TRUE){
      plots.list[[k]] <- plots.list[[k]] +
        ## Add a geom_point object of the line and set to invisible (scatter plot required for marginal density using ggMarginal)
        ## Subset to ignore the confidence intervals when doing the density plots
        ggplot2::geom_point(ggplot2::aes(x = pred, y = obs), col = grDevices::rgb(0, 0, 0, alpha = 0))

      ## Add ggMarginal
      plots.list[[k]] <- ggExtra::ggMarginal(plots.list[[k]], margins = "x", size = marg.density.size, type = marg.density.type, colour = "red")
      ### If marginal rug plot has been requested
    } else if (marg.rug == TRUE){
      plots.list[[k]] <- plots.list[[k]] +
        ggplot2::geom_rug(ggplot2::aes(x = pred, y = obs), col = grDevices::rgb(1, 0, 0, alpha = marg.rug.transparency))
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
    if (marg.density == FALSE){
      plots.list <- ggpubr::ggarrange(plotlist = plots.list, nrow = nrow, ncol = ncol, common.legend = TRUE)
    } else if (marg.density == TRUE){
      plots.list <- gridExtra::arrangeGrob(grobs = plots.list,
                                           layout_matrix = base::matrix(base::seq_len(nrow*ncol),
                                                                        nrow = nrow,
                                                                        ncol = ncol,
                                                                        byrow = TRUE),
                                           top = NULL)
    }
  }

  ### Return output object
  return(plots.list)

}
