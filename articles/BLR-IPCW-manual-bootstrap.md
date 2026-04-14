# BLR-IPCW-manual-bootstrap

## Data preperation

This vignette showcases how to manually apply a bootstrap procedure to
estimate a confidence interval for the BLR-IPCW calibration curves when
using a custom function to estimate the inverse probability of censoring
weights. We use data from the European Society for Blood and Marrow
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

t <- 1826
```

## Confidence intervals for BLR-IPCW using internal bootstrapping procedure

We now remind ourselves of the procedure for generating a confidence
interval, and estimate confidence intervals for the BLR-IPCW calibration
curves using the internal bootstrapping procedure, as is done in the
[Overview
vignette](https://alexpate30.github.io/calibmsm/articles/overview_vignette.pdf).

1.  Resample validation dataset with replacement
2.  Landmark the dataset for assessment of calibration
3.  Calculate inverse probability of censoring weights
4.  Fit the preferred calibration model in the landmarked dataset
    (restricted cubic splines or loess smoother)
5.  Generate observed event probabilities for a fixed vector of
    predicted transition probabilities (specifically the predicted
    transition probabilities from the non-bootstrapped landmark
    validation dataset)

The code to produce the calibration curves with confidence intervals is
as follows:

``` r
dat_calib_blr <-
  calib_msm(data_ms = msebmtcal,
           data_raw = ebmtcal,
           j=1,
           s=0,
           t = t,
           tp_pred = tps0 |>
             dplyr::filter(j == 1) |>
             dplyr::select(any_of(paste("pstate", 1:6, sep = ""))),
           calib_type = 'blr',
           curve_type = "rcs",
           rcs_nk = 3,
           w_covs = c("year", "agecl", "proph", "match"),
           CI = 95,
           CI_R_boot = 200)
```

``` r
plot(dat_calib_blr, combine = TRUE)
#> TableGrob (2 x 1) "arrange": 2 grobs
#>   z     cells    name              grob
#> 1 1 (1-1,1-1) arrange   gtable[arrange]
#> 2 2 (2-2,1-1) arrange gtable[guide-box]
```

## Confidence intervals for BLR-IPCW by inputting custom weights funciton into calibmsm

While the above approach is accessible, it may lead to a misspecified
model for the weights given the user has no control over the choice of
model (Cox model) or the functional form of the predictor variables
$\textbf{𝐙}$, where all continuous variables are assumed to have linear
effects on the hazard, and there are no interaction terms. Potential
issues with using the internal process for estimating the weights is
exemplified in the [Evaluation-of-estimation-of-IPCWs
vignette](https://alexpate30.github.io/calibmsm/articles/Evaluation-of-estimation-of-IPCWs.md).
We therefore encourage users to explore their data and produce a more
suitable model to estimate the weights. Please see the [Overview
vignette](https://alexpate30.github.io/calibmsm/articles/overview_vignette.pdf)
for a definition of the estimand of the weights and how to estimate
them.

Once a function for estimating the weights has been defined, the
internal bootstrapping procedure can then be applied with this function
to estimate the weights. This function must contain the parameters that
are in the `calc_weights` function, even if you do not plan to use these
parameters. These parameters are then specified through arguments with
the `w.` prefix in `calib_msm`. The user specified function may also
contain extra parameters not in `calc_weights`, and these are specified
through the `...` argument in `calib_msm`. We showcase this by defining
a new function `calc_weights_manual`. For this example, we define it to
be the same as `calc_weights` function. In practice, we do not recommend
this, as this option is to allow the user to specify their own function
for estimating the weights. Given our choice of function to estimate the
weights, the calibration curve is the same as that from Figure 1, and
the confidence interval is a similar size.

``` r
calc_weights_manual <- calibmsm::calc_weights

dat_calib_blr_w_function <-
  calib_msm(data_ms = msebmtcal,
           data_raw = ebmtcal,
           j=1,
           s=0,
           t = t,
           tp_pred = tps0 |>
             dplyr::filter(j == 1) |>
             dplyr::select(any_of(paste("pstate", 1:6, sep = ""))),
           calib_type = 'blr',
           curve_type = "rcs",
           rcs_nk = 3,
           w_function = calc_weights_manual,
           w_covs = c("year", "agecl", "proph", "match"),
           CI = 95,
           CI_R_boot = 200)
```

``` r
plot(dat_calib_blr_w_function, combine = TRUE, nrow = 2, ncol = 3)
#> TableGrob (2 x 1) "arrange": 2 grobs
#>   z     cells    name              grob
#> 1 1 (1-1,1-1) arrange   gtable[arrange]
#> 2 2 (2-2,1-1) arrange gtable[guide-box]
```

## Confidence intervals for BLR-IPCW by manually running bootstrapping procedure with custom function for estimating the weights

Finally, we showcase how to run the bootstrapping procedure manually.
This allows users to parallelise their code to reduce computational
time. This is done by using the package *boot* (Canty and Ripley 2022)
in conjunction with `calib_msm`.

We start by creating a permanent object with the predicted risks of each
individual. When the validation cohort `ebmtcal` is resampled, the
appropriate predicted transition probabilities must also be resampled
from this object.

``` r
tp_pred <- tps0 |>
  dplyr::filter(j == 1) |>
  dplyr::select(any_of(paste("pstate", 1:6, sep = "")))
```

Next define the predicted transition probabilities over which the
observed event probabilities will be plotted, which must be the same for
every bootstrapped calibration curve. These are the predicted transition
probabilities of the individuals who are uncensored at time `t`, which
are the predicted transition probabilities for the original calibration
curves.

``` r
## Extract ids for individuals uncensored at t
ids_uncens <- ebmtcal |>
  subset(dtcens > t | (dtcens < t & dtcens_s == 0)) |>
  dplyr::pull(id)
## Extract the predicted risks out of state 1 for these individuals
tp_pred_plot <- tps0 |>
  dplyr::filter(j == 1 & id %in% ids_uncens) |>
  dplyr::select(any_of(paste("pstate", 1:6, sep = "")))
```

The calibration curves (with no confidence interval) are now calculated.
We utilise the `calc_weights` function from **calibmsm** to estimate the
weights, which is the function used internally. However, in practice the
point of this exercise is to create a better performing model to
estimate the weights than the one from the internal procedure. The
vector of weights should have an entry for every row the validation
dataset. For individuals who will not be included in the landmark
cohort, or who are censored prior to time $t$, these weights can take
any value as they will not be included when the calibration model is
fitted.

``` r
weights_manual <-
  calc_weights(data_ms = msebmtcal,
               data_raw = ebmtcal,
               covs = c("year", "agecl", "proph", "match"),
               t = t,
               s = 0,
               landmark_type = "state",
               j = 1,
               max_weight = 10,
               stabilised = FALSE)$ipcw
str(weights_manual)
#>  num [1:2279] NA 1.14 NA 1.01 1.03 ...
```

The observed event probabilities can then be estimated using `calib_msm`
and a call to the `weights` argument.

``` r
dat_calib_boot_manual <-
  calib_msm(data_ms = msebmtcal,
           data_raw = ebmtcal,
           j = 1,
           s = 0,
           t = t,
           tp_pred = tp_pred,
           calib_type = 'blr',
           curve_type = "rcs",
           rcs_nk = 3,
           weights = weights_manual)
```

It is now time to estimate the confidence interval for this curve. We
define a function `calc_obs_boot` to generate a bootstrapped calibration
curve which is compatible with the `boot` function from the *boot*
(Canty and Ripley 2022) package. Within each bootstrapped dataset,
weights are calculated and then a calibration curve is estimated using
`calib_msm`. The use of the argument `tp.pred.plot` is essential here,
to ensure the bootstrapped observed event probabilities are generated
for the same vectors of predicted transition probabilities. The function
is also written to estimate a curve for the transition probabilities
into a specific state $k$, as the output for `boot` must be a vector,
rather than a matrix. This is done by utilising the `transitions.out`
argument in `calib_msm`.

``` r
calc_obs_boot <- function(data, indices, tp_pred, state_k){
  
  ## Bootstrap dataset and predicted transition probabilities
  data_boot <- data[indices,]
  tp_pred_boot <- tp_pred[indices, ]
  
  ## Calculate weights
  ## In practice - replace this function with your own
  weights_manual <-
    calc_weights(data_ms = msebmtcal[msebmtcal$id %in% data_boot$id, ],
                 data_raw = data_boot,
                 covs = c("year", "agecl", "proph", "match"),
                 t = t,
                 s = 0,
                 landmark_type = "state",
                 j = 1,
                 max_weight = 10,
                 stabilised = FALSE)$ipcw
  
  ## Estimate bootstrapped calibration curve
  curve_est <-
    calib_msm(data_ms = msebmtcal[msebmtcal$id %in% data_boot$id, ],
             data_raw = data_boot,
             j=1,
             s=0,
             t = t,
             tp_pred = tp_pred_boot,
             calib_type = 'blr',
             curve_type = "rcs",
             rcs_nk = 3,
             weights = weights_manual,
             tp_pred_plot = tp_pred_plot,
             transitions_out = state_k)
  
  ## Extract observed event probabilities
  curve_obs <-
    curve_est[["plotdata"]][[paste("state", state_k, sep = "")]]$obs
  
  return(curve_obs)
  
}
```

The size of the confidence interval, the states $k$ which individuals
may transition into from state $j$, and a list to store the data for the
plots are then defined:

``` r
alpha <- (1-95/100)/2
valid_transitions <- which(colSums(tp_pred) != 0)
plot_data_list <- vector("list", length(valid_transitions))
```

The bootstrap procedure is now applied to the validation dataset
`ebmtcal` for each state $k$ in `valid_transitions`, and the upper and
lower confidence bands are stored in a `data.frame` along with the
calibration curve calculated earlier.

``` r
for (k in 1:length(valid_transitions)){
  
  ## Assign state k
  state_k <- valid_transitions[k]
  
  ## Run bootstrapping
  boot_obs <- boot::boot(ebmtcal,
                         calc_obs_boot,
                         R = 200,
                         tp_pred = tp_pred,
                         state_k = state_k)$t
  
  ## Extract confidence bands
  lower <- apply(boot_obs, 2, stats::quantile, probs = alpha, na.rm = TRUE)
  upper <- apply(boot_obs, 2, stats::quantile, probs = 1-alpha, na.rm = TRUE)
  
  ## Assign output
  plot_data_list[[k]] <- data.frame(
    "pred" = dat_calib_boot_manual[["plotdata"]][[k]]$pred,
    "obs" = dat_calib_boot_manual[["plotdata"]][[k]]$obs,
    "obs_lower" = lower,
    "obs_upper" = upper)
  
}
```

Finally, metadata is added and the appropriate `class` is defined so
that the data is in the same format as the output from `calib_msm`,
meaning it can be used with the S3 generic `plot`. Note that in this
example `assessed.transitions` = `valid_transitions`, but this may not
be the case if you only estimated calibration curves for a subset of the
possible transitions.

``` r
metadata <- list("valid_transitions"= valid_transitions,
                 "assessed_transitions" = valid_transitions,
                 "CI" = 95,
                 "curve_type" = "rcs")
dat_calib_blr_manual <- list("plotdata" = plot_data_list, "metadata" = metadata)
attr(dat_calib_blr_manual, "class") <- c("calib_blr", "calib_msm")
```

The calibration curves with manually estimated confidence intervals are
plotted in Figure 3. This Figure is very similar to Figures 1 and 2.
This is to be expected given we have used the same model to estimate the
weights that is used in the internal procedure, and verifies that the
manual procedure for calculating the confidence interval has been
successful.

``` r
plot(dat_calib_blr_manual)
#> TableGrob (2 x 1) "arrange": 2 grobs
#>   z     cells    name              grob
#> 1 1 (1-1,1-1) arrange   gtable[arrange]
#> 2 2 (2-2,1-1) arrange gtable[guide-box]
```

Code from these examples can be re-utilised, all that is required is to
specify a function to estimate the inverse probability of censoring
weights.

Canty, Angelo, and Brian Ripley. 2022. “boot: Bootstrap R (S-Plus)
Functions.” <https://cran.r-project.org/package=boot>.

EBMT. 2023. “Data from the European Society for Blood and Marrow
Transplantation.”
<https://search.r-project.org/CRAN/refmans/mstate/html/EBMT-data.html>.

Wreede, Liesbeth C de, Marta Fiocco, and Hein Putter. 2011. “mstate: An
R Package for the Analysis of Competing Risks and Multi-State Models.”
*Journal of Statistical Software* 38 (7).
<https://cran.r-project.org/package=mstate>.
