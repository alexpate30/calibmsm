
<!-- README.md is generated from README.Rmd. Please edit that file -->

# calibmsm

<!-- badges: start -->

[![R-CMD-check](https://github.com/alexpate30/calibmsm/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/alexpate30/calibmsm/actions/workflows/R-CMD-check.yaml)
[![Codecov test
coverage](https://codecov.io/gh/alexpate30/calibmsm/branch/master/graph/badge.svg)](https://app.codecov.io/gh/alexpate30/calibmsm/?branch=master)
<!-- badges: end -->

The goal of calibmsm is to provide a set of tools for producing
calibration plots for validating an existing (i.e. previously developed)
multistate model.

## Installation

You can install the development version of calibmsm from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("alexpate30/calibmsm")
```

## Example

This is a basic example which shows you how to assess the calibration of
the transition probabilities at 5 years follow up for individuals out of
state `j` at time `s`. The predicted transition probabilities are stored
in `tps0`, the individuals data are stored in `ebmtcal`, and the data is
`msdata` format are stored in `msebmtcal`. Inverse probability of
censoring weights are calculated based on variables year, age,
prophylaxis and donor gender match. The calibration curves are produced
using restricted cubic splines. A 95% confidence interval is caulcated
using bootstrapping with 200 bootstrap replicates.

``` r
## Load calibmsm
library(calibmsm)

## Extract relevant predicted risks from tps0
tp.pred <- dplyr::select(dplyr::filter(tps0, j == 1), any_of(paste("pstate", 1:6, sep = "")))

## Calculate observed event probabilities
dat.calib.blr <-
  calc_calib_blr(data.mstate = msebmtcal,
                 data.raw = ebmtcal,
                 j=1,
                 s=0,
                 t.eval = 1826,
                 tp.pred = tp.pred,
                 curve.type = "rcs",
                 rcs.nk = 3,
                 w.covs = c("year", "agecl", "proph", "match"),
                 CI = 95,
                 CI.R.boot = 200)

## Plot calibration plots
plot(dat.calib.blr, combine = TRUE, nrow = 2, ncol = 3)
```

<img src="man/figures/README-example-1.png" width="100%" />

## Getting help

If you encounter a bug, please file an issue with a minimal reproducible
example on [GitHub](https://github.com/alexpate30/calibmsm).

NB TO DELETE - A USEFUL NOTE FOR NOW:: You’ll still need to render
`README.Rmd` regularly, to keep `README.md` up-to-date.
`devtools::build_readme()` is handy for this. You could also use GitHub
Actions to re-render `README.Rmd` every time you push. An example
workflow can be found here:
<https://github.com/r-lib/actions/tree/v1/examples>.
