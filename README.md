
<!-- README.md is generated from README.Rmd. Please edit that file -->

# calibmsm

<!-- badges: start -->

[![R-CMD-check](https://github.com/alexpate30/calibmsm/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/alexpate30/calibmsm/actions/workflows/R-CMD-check.yaml)
[![Codecov test
coverage](https://codecov.io/gh/alexpate30/calibmsm/branch/master/graph/badge.svg)](https://app.codecov.io/gh/alexpate30/calibmsm/?branch=master)
<!-- badges: end -->

The goal of **calibmsm** is to provide a set of tools for estimating
calibration plots when validating an existing (i.e. previously
developed) multistate model.

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
state `j = 1` at time `s = 0` using the BLR-IPCW approach. Please see
the *overview* vignette for examples of how to assess calibration using
the pseudo-value and MLR-IPCW approaches, as well as details of the
methodology.

The predicted transition probabilities are stored in `tps0`, the
individuals data are stored in `ebmtcal`, and the data in `msdata`
format are stored in `msebmtcal`. Calibration curves are estimated using
`calib_blr`. Inverse probability of censoring weights are calculated
based on variables `year`, `age`, `prophylaxis` and donor gender
`match`. The calibration curves are estimated using restricted cubic
splines with 3 knots. A 95% confidence interval is calculated using
bootstrapping with 200 bootstrap replicates.

``` r
## Load calibmsm
library(calibmsm)

## Extract relevant predicted risks from tps0
tp.pred <- dplyr::select(dplyr::filter(tps0, j == 1), any_of(paste("pstate", 1:6, sep = "")))

## Calculate observed event probabilities
dat.calib.blr <-
  calib_blr(data.mstate = msebmtcal,
                 data.raw = ebmtcal,
                 j = 1,
                 s = 0,
                 t = 1826,
                 tp.pred = tp.pred,
                 curve.type = "rcs",
                 rcs.nk = 3,
                 w.covs = c("year", "agecl", "proph", "match"),
                 CI = 95,
                 CI.R.boot = 200)

## Produce summary
summary(dat.calib.blr)
#> There were non-zero predicted transition probabilities into states  1,2,3,4,5,6
#> 
#> Calibration curves have been estimated for transitions into states  1,2,3,4,5,6
#> 
#> Calibration was assessed at time 1826 and calibration was assessed in a landmarked cohort of individuals in state j = 1 at time s = 0
#> 
#> A 95% confidence interval was estimated withb200 bootstrap replicates
#> 
#> The estimated calibration curves are stored in list element `plotdata`:
#> 
#> $state1
#>    id       pred       obs  obs.lower obs.upper
#> 2   2 0.11401890 0.1095897 0.08772698 0.1372000
#> 4   4 0.13838778 0.1036308 0.07914540 0.1271335
#> 5   5 0.12332255 0.1051035 0.08375497 0.1293558
#> 7   7 0.09737975 0.1236322 0.08791651 0.1699968
#> 10 10 0.11371889 0.1097779 0.08789552 0.1373795
#> 13 13 0.11385388 0.1096929 0.08781939 0.1372984
#> 
#> $state2
#>    id      pred       obs obs.lower obs.upper
#> 2   2 0.2316569 0.1698031 0.1225450 0.2234683
#> 4   4 0.1836189 0.1855591 0.1552150 0.2131623
#> 5   5 0.1609740 0.1759804 0.1409099 0.2063167
#> 7   7 0.2121470 0.1785688 0.1460522 0.2130925
#> 10 10 0.2315632 0.1698443 0.1227148 0.2234114
#> 13 13 0.2316571 0.1698030 0.1225445 0.2234685
#> 
#> $state3
#>    id       pred        obs  obs.lower obs.upper
#> 2   2 0.08442692 0.12485834 0.09574042 0.1551600
#> 4   4 0.07579429 0.11666056 0.08513277 0.1475401
#> 5   5 0.05508100 0.09189341 0.04939423 0.1381486
#> 7   7 0.06154308 0.10011560 0.06387762 0.1394670
#> 10 10 0.08440940 0.12484341 0.09571208 0.1551223
#> 13 13 0.08257284 0.12323792 0.09251243 0.1526126
#> 
#> $state4
#>    id      pred       obs obs.lower obs.upper
#> 2   2 0.2328398 0.2427580 0.2008174 0.2819309
#> 4   4 0.2179331 0.2243106 0.1904264 0.2554722
#> 5   5 0.1828176 0.1851051 0.1564639 0.2140929
#> 7   7 0.2206335 0.2275985 0.1931240 0.2599470
#> 10 10 0.2326989 0.2425807 0.2008259 0.2816719
#> 13 13 0.2326047 0.2424622 0.2008316 0.2814986
#> 
#> $state5
#>    id      pred       obs obs.lower obs.upper
#> 2   2 0.1481977 0.1909795 0.1619616 0.2204901
#> 4   4 0.1538475 0.1654523 0.1458289 0.1843488
#> 5   5 0.1425950 0.2215190 0.1770159 0.2703490
#> 7   7 0.1441960 0.2123460 0.1721405 0.2552151
#> 10 10 0.1488068 0.1879278 0.1608139 0.2164994
#> 13 13 0.1505092 0.1797461 0.1574633 0.2043647
#> 
#> $state6
#>    id      pred       obs obs.lower obs.upper
#> 2   2 0.1888598 0.2069354 0.1871640 0.2301865
#> 4   4 0.2304185 0.2542212 0.2325763 0.2838723
#> 5   5 0.3352099 0.3163102 0.2789292 0.3542136
#> 7   7 0.2641006 0.2800368 0.2607361 0.3056078
#> 10 10 0.1888028 0.2068586 0.1870630 0.2300992
#> 13 13 0.1888022 0.2068578 0.1870620 0.2300984

## Plot calibration plots
plot(dat.calib.blr, combine = TRUE, nrow = 2, ncol = 3)
```

<img src="man/figures/README-example-1.png" width="100%" />

## Getting help

If you encounter a bug, please file an issue with a minimal reproducible
example on [GitHub](https://github.com/alexpate30/calibmsm).
