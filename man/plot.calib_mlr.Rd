% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_calib.R
\name{plot.calib_mlr}
\alias{plot.calib_mlr}
\title{Plots calibration scatter plots estimated using \code{\link{calib_mlr}}.}
\usage{
\method{plot}{calib_mlr}(
  x,
  ...,
  combine = TRUE,
  ncol = NULL,
  nrow = NULL,
  transparency.plot = 0.25,
  transparency.rug = 0.1
)
}
\arguments{
\item{x}{Object of class 'calib_mlr' generated from \code{\link{calib_mlr}}}

\item{...}{Other}

\item{combine}{Whether to combine into one plot using ggarrange, or return as a list of individual plots}

\item{ncol}{Number of columns for combined calibration plot}

\item{nrow}{Number of rows for combined calibration plot}

\item{transparency.plot}{Degree of transparency for the calibration scatter plot}

\item{transparency.rug}{Degree of transparency for the density rug plot along each axis}
}
\value{
If \code{combine = TRUE}, returns an object of classes \code{gg}, \code{ggplot}, and \code{ggarrange},
as all ggplots have been combined into one object. If \code{combine = FALSE}, returns an object of
class \code{list}, each element containing an object of class \code{gg} and \code{ggplot}.
}
\description{
Plots calibration scatter plots for the transition probabilities of a multistate model
using data estimated from \code{\link{calib_mlr}}.
}
