# Plots calibration curves estimated using [`calib_msm`](https://alexpate30.github.io/calibmsm/reference/calib_msm.md).

Plots calibration curves for the transition probabilities of a
multistate model estimated using BLR-IPCW and pseudo-value approaches.

## Usage

``` r
# S3 method for class 'calib_msm'
plot(
  x,
  ...,
  combine = TRUE,
  ncol = NULL,
  nrow = NULL,
  size_line = 0.5,
  size_text = 12,
  marg_density = TRUE,
  marg_density_size = 5,
  marg_density_type = "density",
  marg_rug = FALSE,
  marg_rug_transparency = 0.1,
  titles_include = TRUE,
  titles = NULL,
  axis_titles_x = NULL,
  axis_titles_text_x = "Predicted risk",
  axis_titles_y = NULL,
  axis_titles_text_y = "Predicted-observed risk",
  legend_include = TRUE,
  legend_seperate = FALSE,
  legend_title = NULL,
  legend_position = "bottom"
)
```

## Arguments

- x:

  Object of class 'calib_msm' generated from
  [`calib_msm`](https://alexpate30.github.io/calibmsm/reference/calib_msm.md).

- ...:

  Other

- combine:

  Whether to combine into one plot using ggarrange, or return as a list
  of individual plots

- ncol:

  Number of columns for combined calibration plot

- nrow:

  Number of rows for combined calibration plot

- size_line:

  Size of line plots

- size_text:

  Size of text in plot

- marg_density:

  Whether to produce marginal density plots TRUE/FALSE

- marg_density_size:

  Size of the main plot relative to the density plots (see
  [`ggMarginal`](https://rdrr.io/pkg/ggExtra/man/ggMarginal.html))

- marg_density_type:

  What type of marginal plot to show (see
  [`ggMarginal`](https://rdrr.io/pkg/ggExtra/man/ggMarginal.html))

- marg_rug:

  Whether to produce marginal rug plots TRUE/FALSE

- marg_rug_transparency:

  Degree of transparency for the density rug plot along each axis

- titles_include:

  Whether to include titles for each individual calibration plots

- titles:

  Vector of titles for the calibration plots. Defaults to "State k" for
  each plot.

- axis_titles_x:

  Position of plots for which to include title on x-axis

- axis_titles_text_x:

  x-axis title

- axis_titles_y:

  Position of plots for which to include title on y-axis

- axis_titles_text_y:

  y-axis title

- legend_include:

  Whether to produce a legend

- legend_seperate:

  = Whether to include legend in plot (FALSE) or as a seperate object
  (TRUE)

- legend_title:

  Title of legend

- legend_position:

  Position of legend

## Value

If `combine = TRUE`, returns an object of classes `gg`, `ggplot`, and
`ggarrange`, as all ggplots have been combined into one object. If
`combine = FALSE`, returns an object of class `list`, each element
containing an object of class `gg` and `ggplot`.

## Examples

``` r
# Estimate and plot BLR-IPCW calibration curves for the predicted transition
# probabilities at time t = 1826, when predictions were made at time
# s = 0 in state j = 1. These predicted transition probabilities are stored in tps0.

# Extract the predicted transition probabilities out of state j = 1
tp_pred <- dplyr::select(dplyr::filter(tps0, j == 1), any_of(paste("pstate", 1:6, sep = "")))

# Now estimate the observed event probabilities for each possible transition.
dat_calib <-
calib_msm(data_ms = msebmtcal,
 data_raw = ebmtcal,
 j=1,
 s=0,
 t = 1826,
 tp_pred = tp_pred,
 w_covs = c("year", "agecl", "proph", "match"))

 # These are then plotted
 plot(dat_calib, combine = TRUE, nrow = 2, ncol = 3)
#> TableGrob (2 x 3) "arrange": 6 grobs
#>   z     cells    name                grob
#> 1 1 (1-1,1-1) arrange ggExtraPlot[layout]
#> 2 2 (1-1,2-2) arrange ggExtraPlot[layout]
#> 3 3 (1-1,3-3) arrange ggExtraPlot[layout]
#> 4 4 (2-2,1-1) arrange ggExtraPlot[layout]
#> 5 5 (2-2,2-2) arrange ggExtraPlot[layout]
#> 6 6 (2-2,3-3) arrange ggExtraPlot[layout]
```
