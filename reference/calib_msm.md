# Assess the calibration of a multistate model

Calculates the underlying data for calibration plots of the predicted
transition probabilities from a multistate model using three methods.

1.  BLR-IPCW: Binary logistic regression with inverse probability of
    censoring weights.

2.  MLR-IPCW: Multinomial logistic regression with inverse probability
    of censoring weights, based on the nominal calibration framework of
    van Hoorde et al. (2014, 2015)

3.  Pseudo-values: Pseudo-values estimated using the Aalen-Johansen
    estimator (Aalen OO, Johansen S, 1978).

## Usage

``` r
calib_msm(
  data_ms,
  data_raw,
  j,
  s,
  t,
  tp_pred,
  tp_pred_plot = NULL,
  calib_type = "blr",
  curve_type = "rcs",
  rcs_nk = 3,
  loess_span = 0.75,
  loess_degree = 2,
  loess_surface = c("interpolate", "direct"),
  loess_statistics = c("approximate", "exact", "none"),
  loess_trace_hat = c("exact", "approximate"),
  loess_cell = 0.2,
  loess_iterations = 4,
  loess_iterTrace = FALSE,
  mlr_smoother_type = c("sm.ps", "sm.os", "s"),
  mlr_ps_int = 4,
  mlr_degree = 3,
  mlr_s_df = 4,
  mlr_niknots = 4,
  weights = NULL,
  w_function = NULL,
  w_covs = NULL,
  w_landmark_type = "state",
  w_max = 10,
  w_stabilised = FALSE,
  w_max_follow = NULL,
  pv_group_vars = NULL,
  pv_n_pctls = NULL,
  pv_precalc = NULL,
  pv_ids = NULL,
  CI = FALSE,
  CI_type = "bootstrap",
  CI_R_boot = NULL,
  CI_seed = NULL,
  transitions_out = NULL,
  assess_moderate = TRUE,
  assess_mean = TRUE,
  ...
)
```

## Arguments

- data_ms:

  Validation data in `msdata` format

- data_raw:

  Validation data in `data.frame` (one row per individual)

- j:

  Landmark state at which predictions were made

- s:

  Landmark time at which predictions were made

- t:

  Follow up time at which calibration is to be assessed

- tp_pred:

  Data frame or matrix of predicted transition probabilities at time t,
  if in state j at time s. There must be a separate column for the
  predicted transition probabilities into every state, even if these
  predicted transition probabilities are 0.

- tp_pred_plot:

  Data frame or matrix of predicted risks for each possible transition
  over which to plot the calibration curves. Argument provided to enable
  user to apply bootstrapping manually.

- calib_type:

  Whether calibration plots are estimated using BLR-IPCW ('blr'),
  MLR-IPCW ('mlr') or pseudo-values ('pv')

- curve_type:

  Whether calibration curves are estimated using restricted cubic
  splines ('rcs') or loess smoothers ('loess')

- rcs_nk:

  Number of knots when curves are estimated using restricted cubic
  splines

- loess_span:

  Span when curves are estimated using loess smoothers

- loess_degree:

  Degree when curves are estimated\_ using loess smoothers

- loess_surface:

  see [`loess.control`](https://rdrr.io/r/stats/loess.control.html)

- loess_statistics:

  see [`loess.control`](https://rdrr.io/r/stats/loess.control.html)

- loess_trace_hat:

  see [`loess.control`](https://rdrr.io/r/stats/loess.control.html)

- loess_cell:

  see [`loess.control`](https://rdrr.io/r/stats/loess.control.html)

- loess_iterations:

  see [`loess.control`](https://rdrr.io/r/stats/loess.control.html)

- loess_iterTrace:

  see [`loess.control`](https://rdrr.io/r/stats/loess.control.html)

- mlr_smoother_type:

  Type of smoothing applied. Takes values `s` (see
  [`s`](https://rdrr.io/pkg/VGAM/man/s.html)), `sm.ps` (see
  [`sm.ps`](https://rdrr.io/pkg/VGAM/man/sm.ps.html)) or `sm.os` (see
  [`sm.os`](https://rdrr.io/pkg/VGAM/man/sm.os.html)).

- mlr_ps_int:

  the number of equally-spaced B spline intervals in the vector spline
  smoother (see [`sm.ps`](https://rdrr.io/pkg/VGAM/man/sm.ps.html))

- mlr_degree:

  the degree of B-spline basis in the vector spline smoother (see
  [`sm.ps`](https://rdrr.io/pkg/VGAM/man/sm.ps.html))

- mlr_s_df:

  degrees of freedom of vector spline (see
  [`s`](https://rdrr.io/pkg/VGAM/man/s.html))

- mlr_niknots:

  number of interior knots (see
  [`sm.os`](https://rdrr.io/pkg/VGAM/man/sm.os.html))

- weights:

  Vector of inverse probability of censoring weights

- w_function:

  Custom function for estimating the inverse probability of censoring
  weights

- w_covs:

  Character vector of variable names to adjust for when calculating
  inverse probability of censoring weights

- w_landmark_type:

  Whether weights are estimated in all individuals uncensored at time s
  ('all') or only in individuals uncensored and in state j at time s
  ('state')

- w_max:

  Maximum bound for inverse probability of censoring weights

- w_stabilised:

  Indicates whether inverse probability of censoring weights should be
  stabilised or not

- w_max_follow:

  Maximum follow up for model calculating inverse probability of
  censoring weights. Reducing this to `t` + 1 may aid in the
  proportional hazards assumption being met in this model.

- pv_group_vars:

  Variables to group by before calculating pseudo-values

- pv_n_pctls:

  Number of percentiles of predicted risk to group by before calculating
  pseudo-values

- pv_precalc:

  Pre-calculated pseudo-values

- pv_ids:

  Id's of individuals to calculate pseudo-values for

- CI:

  Size of confidence intervals as a %

- CI_type:

  Method for estimating confidence interval (currently restricted to
  `bootstrap`)

- CI_R_boot:

  Number of bootstrap replicates when estimating the confidence interval
  for the calibration curve

- CI_seed:

  Seed for bootstrapping procedure

- transitions_out:

  Transitions for which to calculate calibration curves. Will do all
  possible transitions if left as NULL.

- assess_moderate:

  TRUE/FALSE whether to estimate data for calibration plots

- assess_mean:

  TRUE/FALSE whether to estimate mean calibration

- ...:

  Extra arguments to be passed to w_function (custom function for
  estimating weights)

## Value

`calib_msm` returns a list containing two elements: `plotdata` and
`metadata`. The `plotdata` element contains the data for the calibration
plots. This will itself be a list with each element containing
calibration plot data for the transition probabilities into each of the
possible states. Each list element contains patient ids (`id`) from
`data_raw`, the predicted transition probabilities (`pred`) and the
estimated observed event probabilities (`obs`). If a confidence interval
is requested, upper (`obs_upper`) and lower (`obs_lower`) bounds for the
observed event probabilities are also returned. If tp_pred_plot is
specified, column (`id`) is not returned. The `metadata` element
contains metadata including: a vector of the possible transitions, a
vector of which transitions calibration curves have been estimated for,
the size of the confidence interval, the method for estimating the
calibration curve and other user specified information.

## Details

Observed event probabilities at time `t` are estimated for predicted
transition probabilities `tp_pred` out of state `j` at time `s`.

`calib_type = 'blr'` estimates calibration curves using techniques for
assessing the calibration of a binary logistic regression model (Van
Calster et al., 2016). A choice between restricted cubic splines and
loess smoothers for estimating the calibration curve can be made using
`curve_type`. Landmarking (van Houwelingen HC, 2007) is applied to only
assess calibration in individuals who are uncensored and in state `j` at
time `s`. Calibration can only be assessed in individuals who are also
uncensored at time `t`, which is accounted for using inverse probability
of censoring weights (Hernan M, Robins J, 2020). See method BLR-IPCW
from Pate et al., (2024) for a full explanation of the approach.

`calib_type = 'mlr'` estimates calibration scatter plots using a
technique for assessing the calibration of multinomial logistic
regression models, namely the nominal calibration framework of van
Hoorde et al. (2014, 2015). Landmarking (van Houwelingen HC, 2007) is
applied to only assess calibration in individuals who are uncensored and
in state `j` at time `s`. Calibration can only be assessed in
individuals who are also uncensored at time `t`, which is accounted for
using inverse probability of censoring weights (Hernan M, Robins J,
2020). See method BLR-IPCW from Pate et al., (2024) for a full
explanation of the approach.

`calib_type = 'pv'` estimates calibration curves using using
pseudo-values (Andersen PK, Pohar Perme M, 2010) calculated using the
Aalen-Johansen estimator (Aalen OO, Johansen S, 1978). Calibration
curves are generated by regressing the pseudo-values on the predicted
transition probabilities. A choice between restricted cubic splines and
loess smoothers for estimating the calibration curve can be made using
`curve_type`. Landmarking (van Houwelingen HC, 2007) is applied to only
assess calibration in individuals who are uncensored and in state `j` at
time `s`. The nature of pseudo-values means calibration can be assessed
in all landmarked individuals, regardless of their censoring time. See
method Pseudo-value approach from Pate et al., (2024) for a full
explanation of the approach.

Two datasets for the same cohort of inidividuals must be provided.
Firstly, `data_raw` must be a `data.frame` with one row per individual
containing the variables for the time until censoring (`dtcens`), and an
indicator for censoring `dtcens_s`, where (`dtcens_s = 1`) if an
individual is censored at time `dtcens`, and `dtcens_s = 0` otherwise.
When an individual enters an absorbing state, this prevents censoring
from happening (i.e. dtcens_s = 0). `data_raw` must also contain the
desired variables for estimating the weights. Secondly, `data_ms` must
be a dataset of class `msdata`, generated using the `[mstate]` package.
This dataset is used to apply the landmarking and identify which state
individuals are in at time `t`. While `data_ms` can be derived from
`data_raw`, it would be inefficient to do this within
`calibmsm::calib_msm` due to the bootstrapping procedure, and therefore
they must be inputted seperately.

Unless the user specifies the weights using `weights`, the weights are
estimated using a cox-proportional hazard model, assuming a linear
functional form of the variables defined in `w_covs`. We urge users to
specify their own model for estimating the weights. The `weights`
argument must be a vector with length equal to the number of rows of
`data_raw`.

Confidence intervals cannot be produced for the calibration scatter
plots (`calib_type = 'mlr'`). For calibration curves estimated using
`calib_type = 'blr'`, confidence intervals can only be estimated using
bootstrapping (`CI_type = 'bootstrap`). This procedure uses the internal
method for estimating weights, we therefore encourage users to specify
their own bootstrapping procedure, which incorporates their own model
for estimating the weights. Details on how to do this are provided in
the vignette *BLR-IPCW-manual-bootstrap*. For calibration curves
estimated using `calib_type = 'pv'`, confidence intervals can be
estimated using bootstrapping (`CI_type = 'bootstrap`) or parametric
formulae (`CI_type = 'parametric`). For computational reasons we
recommend using the parametric approach.

The calibration plots can be plotted using
[`plot.calib_msm`](https://alexpate30.github.io/calibmsm/reference/plot.calib_msm.md)
and
[`plot.calib_mlr`](https://alexpate30.github.io/calibmsm/reference/plot.calib_mlr.md).

## References

Aalen OO, Johansen S. An Empirical Transition Matrix for Non-Homogeneous
Markov Chains Based on Censored Observations. *Scand J Stat*.
1978;5(3):141-150.

Andersen PK, Pohar Perme M. Pseudo-observations in survival analysis.
*Stat Methods Med Res*. 2010;19(1):71-99. doi:10.1177/0962280209105020

Hernan M, Robins J (2020). “12.2 Estimating IP weights via modeling.” In
*Causal Inference: What If*, chapter 12.2. Chapman Hall/CRC, Boca Raton.

Pate, A., Sperrin, M., Riley, R. D., Peek, N., Van Staa, T., Sergeant,
J. C., Mamas, M. A., Lip, G. Y. H., Flaherty, M. O., Barrowman, M.,
Buchan, I., & Martin, G. P. Calibration plots for multistate risk
predictions models. *Statistics in Medicine*. 2024;April:1–23. doi:
10.1002/sim.10094.

Van Calster B, Nieboer D, Vergouwe Y, De Cock B, Pencina MJ, Steyerberg
EW (2016). “A calibration hierarchy for risk models was defined: From
utopia to empirical data.” *Journal of Clinical Epidemiology*, 74,
167–176. ISSN 18785921. doi:10.1016/j.jclinepi.2015. 12.005. URL
http://dx.doi.org/10.1016/j.jclinepi.2015.12.005

Van Hoorde K, Vergouwe Y, Timmerman D, Van Huffel S, Steyerberg W, Van
Calster B (2014). “Assessing calibration of multinomial risk prediction
models.” *Statistics in Medicine*, 33(15), 2585–2596.
doi:10.1002/sim.6114.

Van Hoorde K, Van Huffel S, Timmerman D, Bourne T, Van Calster B (2015).
“A spline-based tool to assess and visualize the calibration of
multiclass risk predictions.” *Journal of Biomedical Informatics*, 54,
283–293. ISSN 15320464. doi:10.1016/j.jbi.2014.12.016. URL
http://dx.doi.org/10.1016/j.jbi.2014.12.016.

van Houwelingen HC (2007). “Dynamic Prediction by Landmarking in Event
History Analysis.” *Scandinavian Journal of Statistics*, 34(1), 70–85.

Yee TW (2015). *Vector Generalized Linear and Additive Models*. 1
edition. Springer New, NY. ISBN 978-1-4939-4198-8.
doi:10.1007/978-1-4939-2818-7. URL
https://link.springer.com/book/10.1007/978-1-4939-2818-7.

## Examples

``` r
# Estimate BLR-IPCW calibration curves for the predicted transition
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

# Summarise the output
summary(dat_calib)
#> The method used to assess calibration was BLR-IPCW
#> 
#> There were non-zero predicted transition probabilities into states  1,2,3,4,5,6
#> 
#> Calibration curves have been estimated for transitions into states  1,2,3,4,5,6
#> 
#> Calibration was assessed at time 1826 and calibration was assessed in a landmarked cohort of individuals in state j = 1 at time s = 0
#> 
#> A confidence interval was not estimated
#> 
#> The estimated data for calibration plots are stored in list element `plotdata`:
#> 
#> $state1
#>   id      pred       obs
#> 2  2 0.1140189 0.1095897
#> 4  4 0.1383878 0.1036308
#> 
#> $state2
#>   id      pred       obs
#> 2  2 0.2316569 0.1698031
#> 4  4 0.1836189 0.1855591
#> 
#> $state3
#>   id       pred       obs
#> 2  2 0.08442692 0.1248583
#> 4  4 0.07579429 0.1166606
#> 
#> $state4
#>   id      pred       obs
#> 2  2 0.2328398 0.2427580
#> 4  4 0.2179331 0.2243106
#> 
#> $state5
#>   id      pred       obs
#> 2  2 0.1481977 0.1909795
#> 4  4 0.1538475 0.1654523
#> 
#> $state6
#>   id      pred       obs
#> 2  2 0.1888598 0.2069354
#> 4  4 0.2304185 0.2542212
#> 
#> 
#> 
#> The estimated mean calibration are stored in list element `mean`:
#> 
#>        state1        state2        state3        state4        state5 
#> -0.0216273416 -0.0152282576  0.0254839288  0.0097158314 -0.0003011927 
#>        state6 
#>  0.0032309988 
```
