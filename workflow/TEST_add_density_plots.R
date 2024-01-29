library(calibmsm)

# s = 0 in state j = 1. These predicted transition probabilities are stored in tps0.
devtools::load_all()
devtools::document()
 # Extract the predicted transition probabilities out of state j = 1
 tp.pred <- dplyr::select(dplyr::filter(tps0, j == 1), any_of(paste("pstate", 1:6, sep = "")))

dat.calib <-
calib_msm(data.mstate = msebmtcal,
data.raw = ebmtcal,
j=1,
s=0,
 t = 1826,
tp.pred = tp.pred,
w.covs = c("year", "agecl", "proph", "match"))
library(ggExtra)
library(gridExtra)
png("workflow/plottest.png")
test <- plot(dat.calib, marg.density = TRUE)

class(test)
dev.off()

png("workflow/plottest2.png")
test <- plot(dat.calib)
test
dev.off()


class(test)
test[[1]]
test[[2]]
str(test[[1]])
test.comb <- ggpubr::ggarrange(test[[1]], test[[2]], nrow = 1, ncol = 2, common.legend = TRUE)
test.comb

#extract legend
#https://github.com/hadley/ggplot2/wiki/Share-a-legend-between-two-ggplot2-graphs
g_legend<-function(a.gplot){
  tmp <- ggplot2::ggplot_gtable(ggplot2::ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}

mylegend<-g_legend(test[[1]])

p3 <- grid.arrange(arrangeGrob(p1 + theme(legend.position="none"),
                               p2 + theme(legend.position="none"),
                               nrow=1),
                   mylegend, nrow=2,heights=c(10, 1))

test.comb <- gridExtra::marrangeGrob(grobs = list(test[[1]], test[[2]]), nrow = 2, ncol = 1)
test.comb

plot(dat.calib)
plot(dat.calib, marg.density = TRUE)
png("workflow/plottest.png")
plot(dat.calib, marg.rug = TRUE)
dev.off()
plot.calib_msm <- function(x, ..., combine = TRUE, ncol = NULL, nrow = NULL, transparency.rug = 0.1){

  x <- dat.calib

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
    k <- 1
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




plots.out <- ggplot2::ggplot(data = plot.data.k |> dplyr::arrange(pred) |> dplyr::select(id, pred, obs)) +
  ggplot2::geom_line(ggplot2::aes(x = pred, y = obs)) +
  ggplot2::geom_point(ggplot2::aes(x = pred, y = obs), col = grDevices::rgb(0, 0, 0, alpha = 0)) +
  ggplot2::geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  ggplot2::xlab("Predicted risk") + ggplot2::ylab("Observed risk") +
  ggplot2::xlim(c(0, max(plot.data.k$pred))) +
  ggplot2::ylim(c(min(plot.data.k$obs), max(plot.data.k$obs))) +
  ggplot2::theme(legend.position = "none") +
  ggplot2::ggtitle(paste("State ", state.k, sep = ""))

plots.out.marg <- ggExtra::ggMarginal(plots.out, margins = "x")

