# Calibration-curves-estimated-with-loess-smoothers

## Data preperation

This vignette showcases how to estimate BLR-IPCW and pseudo-value
calibration curves when using loess smoothers for the model to estimate
the observed event probabilities. This is in contrast to using
restricted cubic splines, which was done so in the [Overview
vignette](https://alexpate30.github.io/calibmsm/articles/overview_vignette.pdf).
Both approaches are viable and we refer to the literature for a
discussion of each method (Austin, Harrell, and Klaveren 2020; Harrell
2015). We mimic exactly the analysis from the [Overview
vignette](https://alexpate30.github.io/calibmsm/articles/overview_vignette.pdf),
specifying loess smoother rather than restricted cubic splines, and
assess whether the conclusions would be any different. This is primarily
to showcase the utility of the package.

We use data from the European Society for Blood and Marrow
Transplantation (EBMT 2023), which contains multistate survival data
after a transplant for patients with blood cancer. The start of follow
up is the day of the transplant and the initial state is alive and in
remission. There are three intermediate events ($2$: recovery, $3$:
adverse event, or $4$: recovery + adverse event), and two absorbing
states ($5$: relapse and $6$: death). This data was originally made
available from the `mstate` package (Wreede, Fiocco, and Putter 2011).

We start by reminding ourselves of the EBMT (Wreede, Fiocco, and Putter
2011; EBMT 2023) validation datasets `ebmtcal` and `msebmtcal`, and the
predicted transition probabilities out of state `j = 1` made at time
`s = 0` from the multistate model (`tps0`). We also set the time at
which we want to evaluate the transition probabilities `t`. Please refer
to the [Overview
vignette](https://alexpate30.github.io/calibmsm/articles/overview_vignette.pdf)
for a more detailed description of this data.

``` r
set.seed(101)

library(calibmsm)

data("ebmtcal")
head(ebmtcal)
#>   id  rec rec.s   ae ae.s recae recae.s  rel rel.s  srv srv.s      year agecl
#> 1  1   22     1  995    0   995       0  995     0  995     0 1995-1998 20-40
#> 2  2   29     1   12    1    29       1  422     1  579     1 1995-1998 20-40
#> 3  3 1264     0   27    1  1264       0 1264     0 1264     0 1995-1998 20-40
#> 4  4   50     1   42    1    50       1   84     1  117     1 1995-1998 20-40
#> 5  5   22     1 1133    0  1133       0  114     1 1133     0 1995-1998   >40
#> 6  6   33     1   27    1    33       1 1427     0 1427     0 1995-1998 20-40
#>   proph              match dtcens dtcens_s
#> 1    no no gender mismatch    995        1
#> 2    no no gender mismatch    422        0
#> 3    no no gender mismatch   1264        1
#> 4    no    gender mismatch     84        0
#> 5    no    gender mismatch    114        0
#> 6    no no gender mismatch   1427        1

data("msebmtcal")
head(msebmtcal)
#>   id from to trans Tstart Tstop time status
#> 1  1    1  2     1      0    22   22      1
#> 2  1    1  3     2      0    22   22      0
#> 3  1    1  5     3      0    22   22      0
#> 4  1    1  6     4      0    22   22      0
#> 5  1    2  4     5     22   995  973      0
#> 6  1    2  5     6     22   995  973      0

data("tps0")
head(tps0)
#>   id   pstate1   pstate2    pstate3   pstate4   pstate5   pstate6        se1
#> 1  1 0.1139726 0.2295006 0.08450376 0.2326861 0.1504855 0.1888514 0.01291133
#> 2  2 0.1140189 0.2316569 0.08442692 0.2328398 0.1481977 0.1888598 0.01291552
#> 3  3 0.1136646 0.2317636 0.08274331 0.2325663 0.1504787 0.1887834 0.01289444
#> 4  4 0.1383878 0.1836189 0.07579429 0.2179331 0.1538475 0.2304185 0.01857439
#> 5  5 0.1233226 0.1609740 0.05508100 0.1828176 0.1425950 0.3352099 0.01944967
#> 6  6 0.1136646 0.2317636 0.08462424 0.2305854 0.1505534 0.1888087 0.01289444
#>          se2        se3        se4        se5        se6 j
#> 1 0.02369584 0.01257251 0.02323376 0.01648630 0.01601795 1
#> 2 0.02374329 0.01256056 0.02324869 0.01632797 0.01603703 1
#> 3 0.02375770 0.01245752 0.02322375 0.01647890 0.01601525 1
#> 4 0.03004447 0.01462570 0.03018673 0.02124071 0.02416121 1
#> 5 0.03419721 0.01367768 0.03423941 0.02329644 0.03688586 1
#> 6 0.02375770 0.01257276 0.02317348 0.01649531 0.01602438 1

t_eval <- 1826
```

## Estimation of calibration curves out of state j = 1 at time s = 0 using loess smoothers

We now produce calibration curves for the predicted transition
probabilities out of state $j = 1$ at time $s = 0$. Given all
individuals start in state $1$, there is no need to consider the
transition probabilities out of states $j \neq 1$ at $s = 0$.
Calibration is assessed at follow up time ($t = 1826$ days). We start by
extracting the predicted transition probabilities from state $j = 1$ at
time $s = 0$ from the object . These are the transition probabilities we
aim to assess the calibration of.

``` r
tp_pred_s0 <- tps0 |>
  dplyr::filter(j == 1) |>
  dplyr::select(any_of(paste("pstate", 1:6, sep = "")))
```

We now produce calibration curves using loess smoothers for the
regression model where the observed event probabilities are estimated.
This is specified by . Behaviour of the curve can be controlled using
and , which we leave as default in this example. Please refere to
documentation for the function from the package for details of these
parameters. A smaller span and/or degree will create a more smoother
line. This is done for both the BLR-IPCW and pseudo-value calibration
curves.

``` r
dat_calib_blr <-
  calib_msm(data_ms = msebmtcal,
           data_raw = ebmtcal,
           j = 1,
           s = 0,
           t = t_eval,
           tp_pred = tp_pred_s0,
           calib_type = 'blr',
           curve_type = "loess",
           w_covs = c("year", "agecl", "proph", "match"),
           CI = 95,
           CI_type = "bootstrap",
           CI_R_boot = 200)
#> Warning in calib_blr_ipcw(data_raw = data_raw, data_ms = data_ms, tp_pred_plot = tp_pred_plot, : WARNING, SOME BOOTSTRAPPED OBSERVED EVENT PROBABILITIES WERE NA FOR STATE 1 
#>  THERE ARE  85  ITERATIONS WITH NA's 
#>  THE MEAN NUMBER OF NA's IN EACH ITERATION IS 0.835
#> Warning in calib_blr_ipcw(data_raw = data_raw, data_ms = data_ms, tp_pred_plot = tp_pred_plot, : WARNING, SOME BOOTSTRAPPED OBSERVED EVENT PROBABILITIES WERE NA FOR STATE 2 
#>  THERE ARE  125  ITERATIONS WITH NA's 
#>  THE MEAN NUMBER OF NA's IN EACH ITERATION IS 1.205
#> Warning in calib_blr_ipcw(data_raw = data_raw, data_ms = data_ms, tp_pred_plot = tp_pred_plot, : WARNING, SOME BOOTSTRAPPED OBSERVED EVENT PROBABILITIES WERE NA FOR STATE 3 
#>  THERE ARE  133  ITERATIONS WITH NA's 
#>  THE MEAN NUMBER OF NA's IN EACH ITERATION IS 1.29
#> Warning in calib_blr_ipcw(data_raw = data_raw, data_ms = data_ms, tp_pred_plot = tp_pred_plot, : WARNING, SOME BOOTSTRAPPED OBSERVED EVENT PROBABILITIES WERE NA FOR STATE 4 
#>  THERE ARE  115  ITERATIONS WITH NA's 
#>  THE MEAN NUMBER OF NA's IN EACH ITERATION IS 1.025
#> Warning in calib_blr_ipcw(data_raw = data_raw, data_ms = data_ms, tp_pred_plot = tp_pred_plot, : WARNING, SOME BOOTSTRAPPED OBSERVED EVENT PROBABILITIES WERE NA FOR STATE 5 
#>  THERE ARE  117  ITERATIONS WITH NA's 
#>  THE MEAN NUMBER OF NA's IN EACH ITERATION IS 1.16
#> Warning in calib_blr_ipcw(data_raw = data_raw, data_ms = data_ms, tp_pred_plot = tp_pred_plot, : WARNING, SOME BOOTSTRAPPED OBSERVED EVENT PROBABILITIES WERE NA FOR STATE 6 
#>  THERE ARE  111  ITERATIONS WITH NA's 
#>  THE MEAN NUMBER OF NA's IN EACH ITERATION IS 1.115

dat_calib_pv <-
  calib_msm(data_ms = msebmtcal,
           data_raw = ebmtcal,
           j = 1,
           s = 0,
           t = t_eval,
           tp_pred = tp_pred_s0,
           calib_type = 'pv',
           curve_type = "loess",
           pv_group_vars = c("year"),
           pv_n_pctls = 3,
           CI = 95,
           CI_type = 'parametric')
```

We now produce calibration plots using the function.

``` r
plot(dat_calib_blr)
#> TableGrob (2 x 1) "arrange": 2 grobs
#>   z     cells    name              grob
#> 1 1 (1-1,1-1) arrange   gtable[arrange]
#> 2 2 (2-2,1-1) arrange gtable[guide-box]
```

``` r
plot(dat_calib_pv)
#> TableGrob (2 x 1) "arrange": 2 grobs
#>   z     cells    name              grob
#> 1 1 (1-1,1-1) arrange   gtable[arrange]
#> 2 2 (2-2,1-1) arrange gtable[guide-box]
```

The equivalent Figures from the [Overview
vignette](https://alexpate30.github.io/calibmsm/articles/overview_vignette.pdf),
produced using restricted cubic splines, are Figures 2 and 3.

Write a short comparison… XXXX

## Estimation of calibration curves out of state j = 1 at time s = 100 using loess smoothers

We now produce calibration curves for the predicted transition
probabilities out of state $j = 1$ at time $s = 100$. Calibration is
assessed at follow up time ($t = 1826$ days). We start by extracting the
predicted transition probabilities from state $j = 1$ at time $s = 100$
from the object . These are the transition probabilities we aim to
assess the calibration of.

``` r
tp_pred_s100 <- tps100 |>
  dplyr::filter(j == 1) |>
  dplyr::select(any_of(paste("pstate", 1:6, sep = "")))
```

We now produce calibration curves using loess smoothers for the
regression model where the observed event probabilities are estimated.
This is specified by . Behaviour of the curve can be controlled using
and , which we leave as default in this example. Please refere to
documentation for the function from the package for details of these
parameters. A smaller span and/or degree will create a more smoother
line. This is done for both the BLR-IPCW and pseudo-value calibration
curves.

``` r
dat_calib_blr <-
  calib_msm(data_ms = msebmtcal,
           data_raw = ebmtcal,
           j = 1,
           s = 100,
           t = t_eval,
           tp_pred = tp_pred_s100,
           calib_type = 'blr',
           curve_type = "loess",
           w_covs = c("year", "agecl", "proph", "match"),
           CI = 95,
           CI_type = "bootstrap",
           CI_R_boot = 200)
#> In the landmark cohort of individuals uncensored and in state j at time s, states {2} have less than 30 people in them at the time calibration is being assessed (t).
#>               Warnings and errors may occur when the models to estimate the calibration curves are fitted, due to small sample size.
#>               The number to flag this warning (30) has been chosen arbitrarily, and does not constitute a sufficient sample size from
#>               a statistical point of view.
#> Warning in calib_msm(data_ms = msebmtcal, data_raw = ebmtcal, j = 1, s = 100, :
#> Warning in calib_blr_ipcw(data_raw = data_raw, data_ms = data_ms, tp_pred_plot = tp_pred_plot, : WARNING, SOME BOOTSTRAPPED OBSERVED EVENT PROBABILITIES WERE NA FOR STATE 1 
#>  THERE ARE  124  ITERATIONS WITH NA's 
#>  THE MEAN NUMBER OF NA's IN EACH ITERATION IS 1.105
#> Warning in calib_blr_ipcw(data_raw = data_raw, data_ms = data_ms, tp_pred_plot = tp_pred_plot, : WARNING, SOME BOOTSTRAPPED OBSERVED EVENT PROBABILITIES WERE NA FOR STATE 2 
#>  THERE ARE  120  ITERATIONS WITH NA's 
#>  THE MEAN NUMBER OF NA's IN EACH ITERATION IS 1.205
#> Warning in calib_blr_ipcw(data_raw = data_raw, data_ms = data_ms, tp_pred_plot = tp_pred_plot, : WARNING, SOME BOOTSTRAPPED OBSERVED EVENT PROBABILITIES WERE NA FOR STATE 5 
#>  THERE ARE  99  ITERATIONS WITH NA's 
#>  THE MEAN NUMBER OF NA's IN EACH ITERATION IS 0.955
#> Warning in calib_blr_ipcw(data_raw = data_raw, data_ms = data_ms, tp_pred_plot = tp_pred_plot, : WARNING, SOME BOOTSTRAPPED OBSERVED EVENT PROBABILITIES WERE NA FOR STATE 6 
#>  THERE ARE  140  ITERATIONS WITH NA's 
#>  THE MEAN NUMBER OF NA's IN EACH ITERATION IS 1.375

dat_calib_pv <-
  calib_msm(data_ms = msebmtcal,
           data_raw = ebmtcal,
           j = 1,
           s = 100,
           t = t_eval,
           tp_pred = tp_pred_s100,
           calib_type = 'pv',
           curve_type = "loess",
           pv_group_vars = c("year"),
           CI = 95,
           CI_type = 'parametric')
#> In the landmark cohort of individuals uncensored and in state j at time s, states {2} have less than 30 people in them at the time calibration is being assessed (t).
#>               Warnings and errors may occur when the models to estimate the calibration curves are fitted, due to small sample size.
#>               The number to flag this warning (30) has been chosen arbitrarily, and does not constitute a sufficient sample size from
#>               a statistical point of view.
#> Warning in calib_msm(data_ms = msebmtcal, data_raw = ebmtcal, j = 1, s = 100, :
```

We now produce calibration plots using the function.

``` r
plot(dat_calib_blr)
#> TableGrob (2 x 1) "arrange": 2 grobs
#>   z     cells    name              grob
#> 1 1 (1-1,1-1) arrange   gtable[arrange]
#> 2 2 (2-2,1-1) arrange gtable[guide-box]
```

``` r
plot(dat_calib_pv)
#> TableGrob (2 x 1) "arrange": 2 grobs
#>   z     cells    name              grob
#> 1 1 (1-1,1-1) arrange   gtable[arrange]
#> 2 2 (2-2,1-1) arrange gtable[guide-box]
```

The equivalent Figures from the [Overview
vignette](https://alexpate30.github.io/calibmsm/articles/overview_vignette.pdf),
produced using restricted cubic splines, are Figures 5 and 6.

Write a short comparison…

## Estimation of calibration curves out of state j = 3 at time s = 100 using loess smoothers

We now produce calibration curves for the predicted transition
probabilities out of state $j = 3$ at time $s = 100$. Calibration is
assessed at follow up time ($t = 1826$ days). We start by extracting the
predicted transition probabilities from state $j = 3$ at time $s = 100$
from the object . These are the transition probabilities we aim to
assess the calibration of.

``` r
tp_pred_s100 <- tps100 |>
  dplyr::filter(j == 3) |>
  dplyr::select(any_of(paste("pstate", 1:6, sep = "")))
```

We now produce calibration curves using loess smoothers for the
regression model where the observed event probabilities are estimated.
This is specified by . Behaviour of the curve can be controlled using
and , which we leave as default in this example. Please refere to
documentation for the function from the package for details of these
parameters. A smaller span and/or degree will create a more smoother
line. This is done for both the BLR-IPCW and pseudo-value calibration
curves.

``` r
dat_calib_blr <-
  calib_msm(data_ms = msebmtcal,
           data_raw = ebmtcal,
           j = 3,
           s = 100,
           t = t_eval,
           tp_pred = tp_pred_s100,
           calib_type = 'blr',
           curve_type = "loess",
           w_covs = c("year", "agecl", "proph", "match"),
           CI = 95,
           CI_type = "bootstrap",
           CI_R_boot = 200)
#> In the landmark cohort of individuals uncensored and in state j at time s, states {4,5} have less than 30 people in them at the time calibration is being assessed (t).
#>               Warnings and errors may occur when the models to estimate the calibration curves are fitted, due to small sample size.
#>               The number to flag this warning (30) has been chosen arbitrarily, and does not constitute a sufficient sample size from
#>               a statistical point of view.
#> Warning in calib_msm(data_ms = msebmtcal, data_raw = ebmtcal, j = 3, s = 100, :
#> Warning in calib_blr_ipcw(data_raw = data_raw, data_ms = data_ms, tp_pred_plot = tp_pred_plot, : WARNING, SOME BOOTSTRAPPED OBSERVED EVENT PROBABILITIES WERE NA FOR STATE 3 
#>  THERE ARE  124  ITERATIONS WITH NA's 
#>  THE MEAN NUMBER OF NA's IN EACH ITERATION IS 1.165
#> Warning in calib_blr_ipcw(data_raw = data_raw, data_ms = data_ms, tp_pred_plot = tp_pred_plot, : WARNING, SOME BOOTSTRAPPED OBSERVED EVENT PROBABILITIES WERE NA FOR STATE 4 
#>  THERE ARE  116  ITERATIONS WITH NA's 
#>  THE MEAN NUMBER OF NA's IN EACH ITERATION IS 1.15
#> Warning in calib_blr_ipcw(data_raw = data_raw, data_ms = data_ms, tp_pred_plot = tp_pred_plot, : WARNING, SOME BOOTSTRAPPED OBSERVED EVENT PROBABILITIES WERE NA FOR STATE 5 
#>  THERE ARE  121  ITERATIONS WITH NA's 
#>  THE MEAN NUMBER OF NA's IN EACH ITERATION IS 1.23
#> Warning in calib_blr_ipcw(data_raw = data_raw, data_ms = data_ms, tp_pred_plot = tp_pred_plot, : WARNING, SOME BOOTSTRAPPED OBSERVED EVENT PROBABILITIES WERE NA FOR STATE 6 
#>  THERE ARE  112  ITERATIONS WITH NA's 
#>  THE MEAN NUMBER OF NA's IN EACH ITERATION IS 1.07

dat_calib_pv <-
  calib_msm(data_ms = msebmtcal,
           data_raw = ebmtcal,
           j = 3,
           s = 100,
           t = t_eval,
           tp_pred = tp_pred_s100,
           calib_type = 'pv',
           curve_type = "loess",
           pv_group_vars = c("year"),
           CI = 95,
           CI_type = 'parametric')
#> In the landmark cohort of individuals uncensored and in state j at time s, states {4,5} have less than 30 people in them at the time calibration is being assessed (t).
#>               Warnings and errors may occur when the models to estimate the calibration curves are fitted, due to small sample size.
#>               The number to flag this warning (30) has been chosen arbitrarily, and does not constitute a sufficient sample size from
#>               a statistical point of view.
#> Warning in calib_msm(data_ms = msebmtcal, data_raw = ebmtcal, j = 3, s = 100, :
```

We now produce calibration plots using the function.

``` r
plot(dat_calib_blr)
#> TableGrob (2 x 1) "arrange": 2 grobs
#>   z     cells    name              grob
#> 1 1 (1-1,1-1) arrange   gtable[arrange]
#> 2 2 (2-2,1-1) arrange gtable[guide-box]
```

``` r
plot(dat_calib_pv)
#> TableGrob (2 x 1) "arrange": 2 grobs
#>   z     cells    name              grob
#> 1 1 (1-1,1-1) arrange   gtable[arrange]
#> 2 2 (2-2,1-1) arrange gtable[guide-box]
```

The equivalent Figures from the [Overview
vignette](https://alexpate30.github.io/calibmsm/articles/overview_vignette.pdf),
produced using restricted cubic splines, are Figures 7 and 8.

Write a short comparison… XXXX

## A note on the fact that confidence intervals may be lower than 0 or bigger than 1

Given the standard error is calculated on the probability scale here, as
opposed to the logit scale, the confidence intervals may exceed 0 and 1.
I have checked and the confidence intervals are calculated in the same
way in the CalibrationCurves package:
<https://cran.r-project.org/web/packages/CalibrationCurves/CalibrationCurves.pdf>
. In particular see code for val.prob.ci.2:
<https://github.com/BavoDC/CalibrationCurves/blob/master/R/val.prob.ci.2.R>
. There parametric confidence intervals are calculated in the same way
we calculate parametric confidence intervals for the pseudo value
approach. I am therefore confident this is correct, but is something i’d
like to come back to and consider the implications of. It’s also
intriguing that even for the bootstrapped BLR-IPCW confidence intervals,
we are seeing confidence intervals well outside 0 and 1. I thought this
would be an issue for the pseudo-value parametric confidence intervals
as we are just adding on 1.96\*se on the probability scale. However, I
didnt actually think that would be possible for the bootstrapped
confidence intervals, as in any given bootstrapped dataset,
predicted-observed values would still lie between 0 and 1, so this is
something to explore in more detail.

Having done a bit more research (see test_loess_bootstrap.R in the test
folder), the values below 0 are driven by a small number of individuals
with the smallest risks. Lots of their predicted-observed values are
NAs, and they have some negative values too. It is evident these small
number of individuals are driving all the NA values we are finding in
each bootstrap iteration. Going to try increasing the span to see if
that helps. Now i know the mechanism, although not something I’m
particularly worried about, as seems to be driven by small sample size
when assessing calibration at time s = 100.

Compare the following graphs to figures 1 and 2 above. First do .

``` r
dat_calib_blr <-
  calib_msm(data_ms = msebmtcal,
           data_raw = ebmtcal,
           j = 1,
           s = 0,
           t = t_eval,
           tp_pred = tp_pred_s0,
           calib_type = 'blr',
           curve_type = "loess",
           loess_span = 0.9,
           w_covs = c("year", "agecl", "proph", "match"),
           CI = 95,
           CI_type = "bootstrap",
           CI_R_boot = 200)
#> Warning in calib_blr_ipcw(data_raw = data_raw, data_ms = data_ms, tp_pred_plot = tp_pred_plot, : WARNING, SOME BOOTSTRAPPED OBSERVED EVENT PROBABILITIES WERE NA FOR STATE 1 
#>  THERE ARE  94  ITERATIONS WITH NA's 
#>  THE MEAN NUMBER OF NA's IN EACH ITERATION IS 0.91
#> Warning in calib_blr_ipcw(data_raw = data_raw, data_ms = data_ms, tp_pred_plot = tp_pred_plot, : WARNING, SOME BOOTSTRAPPED OBSERVED EVENT PROBABILITIES WERE NA FOR STATE 2 
#>  THERE ARE  118  ITERATIONS WITH NA's 
#>  THE MEAN NUMBER OF NA's IN EACH ITERATION IS 1.16
#> Warning in calib_blr_ipcw(data_raw = data_raw, data_ms = data_ms, tp_pred_plot = tp_pred_plot, : WARNING, SOME BOOTSTRAPPED OBSERVED EVENT PROBABILITIES WERE NA FOR STATE 3 
#>  THERE ARE  116  ITERATIONS WITH NA's 
#>  THE MEAN NUMBER OF NA's IN EACH ITERATION IS 0.95
#> Warning in calib_blr_ipcw(data_raw = data_raw, data_ms = data_ms, tp_pred_plot = tp_pred_plot, : WARNING, SOME BOOTSTRAPPED OBSERVED EVENT PROBABILITIES WERE NA FOR STATE 4 
#>  THERE ARE  119  ITERATIONS WITH NA's 
#>  THE MEAN NUMBER OF NA's IN EACH ITERATION IS 1.135
#> Warning in calib_blr_ipcw(data_raw = data_raw, data_ms = data_ms, tp_pred_plot = tp_pred_plot, : WARNING, SOME BOOTSTRAPPED OBSERVED EVENT PROBABILITIES WERE NA FOR STATE 5 
#>  THERE ARE  131  ITERATIONS WITH NA's 
#>  THE MEAN NUMBER OF NA's IN EACH ITERATION IS 1.315
#> Warning in calib_blr_ipcw(data_raw = data_raw, data_ms = data_ms, tp_pred_plot = tp_pred_plot, : WARNING, SOME BOOTSTRAPPED OBSERVED EVENT PROBABILITIES WERE NA FOR STATE 6 
#>  THERE ARE  115  ITERATIONS WITH NA's 
#>  THE MEAN NUMBER OF NA's IN EACH ITERATION IS 1.145

dat_calib_pv <-
  calib_msm(data_ms = msebmtcal,
           data_raw = ebmtcal,
           j = 1,
           s = 0,
           t = t_eval,
           tp_pred = tp_pred_s0,
           calib_type = 'pv',
           curve_type = "loess",
           loess_span = 0.9,
           pv_group_vars = c("year"),
           pv_n_pctls = 3,
           CI = 95,
           CI_type = 'parametric')
```

We now produce calibration plots using the function.

``` r
plot(dat_calib_blr)
#> TableGrob (2 x 1) "arrange": 2 grobs
#>   z     cells    name              grob
#> 1 1 (1-1,1-1) arrange   gtable[arrange]
#> 2 2 (2-2,1-1) arrange gtable[guide-box]
```

``` r
plot(dat_calib_pv)
#> TableGrob (2 x 1) "arrange": 2 grobs
#>   z     cells    name              grob
#> 1 1 (1-1,1-1) arrange   gtable[arrange]
#> 2 2 (2-2,1-1) arrange gtable[guide-box]
```

Then do .

``` r
dat_calib_blr <-
  calib_msm(data_ms = msebmtcal,
           data_raw = ebmtcal,
           j = 1,
           s = 0,
           t = t_eval,
           tp_pred = tp_pred_s0,
           calib_type = 'blr',
           curve_type = "loess",
           loess_span = 1,
           w_covs = c("year", "agecl", "proph", "match"),
           CI = 95,
           CI_type = "bootstrap",
           CI_R_boot = 200)
#> Warning in calib_blr_ipcw(data_raw = data_raw, data_ms = data_ms, tp_pred_plot = tp_pred_plot, : WARNING, SOME BOOTSTRAPPED OBSERVED EVENT PROBABILITIES WERE NA FOR STATE 1 
#>  THERE ARE  86  ITERATIONS WITH NA's 
#>  THE MEAN NUMBER OF NA's IN EACH ITERATION IS 0.805
#> Warning in calib_blr_ipcw(data_raw = data_raw, data_ms = data_ms, tp_pred_plot = tp_pred_plot, : WARNING, SOME BOOTSTRAPPED OBSERVED EVENT PROBABILITIES WERE NA FOR STATE 2 
#>  THERE ARE  130  ITERATIONS WITH NA's 
#>  THE MEAN NUMBER OF NA's IN EACH ITERATION IS 1.27
#> Warning in calib_blr_ipcw(data_raw = data_raw, data_ms = data_ms, tp_pred_plot = tp_pred_plot, : WARNING, SOME BOOTSTRAPPED OBSERVED EVENT PROBABILITIES WERE NA FOR STATE 3 
#>  THERE ARE  128  ITERATIONS WITH NA's 
#>  THE MEAN NUMBER OF NA's IN EACH ITERATION IS 1.255
#> Warning in calib_blr_ipcw(data_raw = data_raw, data_ms = data_ms, tp_pred_plot = tp_pred_plot, : WARNING, SOME BOOTSTRAPPED OBSERVED EVENT PROBABILITIES WERE NA FOR STATE 4 
#>  THERE ARE  123  ITERATIONS WITH NA's 
#>  THE MEAN NUMBER OF NA's IN EACH ITERATION IS 1.125
#> Warning in calib_blr_ipcw(data_raw = data_raw, data_ms = data_ms, tp_pred_plot = tp_pred_plot, : WARNING, SOME BOOTSTRAPPED OBSERVED EVENT PROBABILITIES WERE NA FOR STATE 5 
#>  THERE ARE  132  ITERATIONS WITH NA's 
#>  THE MEAN NUMBER OF NA's IN EACH ITERATION IS 1.285
#> Warning in calib_blr_ipcw(data_raw = data_raw, data_ms = data_ms, tp_pred_plot = tp_pred_plot, : WARNING, SOME BOOTSTRAPPED OBSERVED EVENT PROBABILITIES WERE NA FOR STATE 6 
#>  THERE ARE  112  ITERATIONS WITH NA's 
#>  THE MEAN NUMBER OF NA's IN EACH ITERATION IS 1.02

dat_calib_pv <-
  calib_msm(data_ms = msebmtcal,
           data_raw = ebmtcal,
           j = 1,
           s = 0,
           t = t_eval,
           tp_pred = tp_pred_s0,
           calib_type = 'pv',
           curve_type = "loess",
           loess_span = 1,
           pv_group_vars = c("year"),
           pv_n_pctls = 3,
           CI = 95,
           CI_type = 'parametric')
```

We now produce calibration plots using the function.

``` r
plot(dat_calib_blr)
#> TableGrob (2 x 1) "arrange": 2 grobs
#>   z     cells    name              grob
#> 1 1 (1-1,1-1) arrange   gtable[arrange]
#> 2 2 (2-2,1-1) arrange gtable[guide-box]
```

``` r
plot(dat_calib_pv)
#> TableGrob (2 x 1) "arrange": 2 grobs
#>   z     cells    name              grob
#> 1 1 (1-1,1-1) arrange   gtable[arrange]
#> 2 2 (2-2,1-1) arrange gtable[guide-box]
```

Austin, Peter C., Frank E. Harrell, and David van Klaveren. 2020.
“Graphical calibration curves and the integrated calibration index (ICI)
for survival models.” *Statistics in Medicine* 39 (21): 2714–42.
<https://doi.org/10.1002/sim.8570>.

EBMT. 2023. “Data from the European Society for Blood and Marrow
Transplantation.”
<https://search.r-project.org/CRAN/refmans/mstate/html/EBMT-data.html>.

Harrell, Frank E. 2015. *Regression Modeling Strategies*. Springer S.
Cham: Springer.

Wreede, Liesbeth C de, Marta Fiocco, and Hein Putter. 2011. “mstate: An
R Package for the Analysis of Competing Risks and Multi-State Models.”
*Journal of Statistical Software* 38 (7).
<https://cran.r-project.org/package=mstate>.
