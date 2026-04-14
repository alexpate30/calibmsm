# Calculate inverse probability of censoring weights at time `t`.

Estimates the inverse probability of censoring weights by fitting a
cox-propotinal hazards model in a landmark cohort of individuals.
Primarily used internally, this function has been exported to allow
users to reproduce results in the vignette when estimating confidence
intervals using bootstrapping manually.

## Usage

``` r
calc_weights(
  data_ms,
  data_raw,
  covs = NULL,
  t,
  s,
  landmark_type = "state",
  j = NULL,
  max_weight = 10,
  stabilised = FALSE,
  max_follow = NULL
)
```

## Arguments

- data_ms:

  Validation data in msdata format

- data_raw:

  Validation data in data.frame (one row per individual)

- covs:

  Character vector of variable names to adjust for when calculating
  inverse probability of censoring weights

- t:

  Follow up time at which to calculate weights

- s:

  Landmark time at which predictions were made

- landmark_type:

  Whether weights are estimated in all individuals uncensored at time s
  ('all') or only in individuals uncensored and in state j at time s
  ('state')

- j:

  Landmark state at which predictions were made (only required in
  landmark_type = 'state')

- max_weight:

  Maximum bound for weights

- stabilised:

  Indicates whether weights should be stabilised or not

- max_follow:

  Maximum follow up for model calculating inverse probability of
  censoring weights. Reducing this to `t` + 1 may aid in the
  proportional hazards assumption being met in this model.

## Value

A data frame with two columns. `id` corresponds to the patient ids from
`data_raw`. `ipcw` contains the inverse probability of censoring weights
(specifically the inverse of the probability of being uncesored). If
`stabilised = TRUE` was specified, a third variable `ipcw_stab` will be
returned, which is the stabilised inverse probability of censoring
weights.

## Details

Estimates inverse probability of censoring weights (Hernan M, Robins J,
2020). Fits a cox proportional hazards model to individuals in a
landmark cohort, predicting the probability of being censored at time
`t`. This landmark cohort may either be all individuals uncensored at
time `s`, or those uncensored and in state `j` at time `s`. All
predictors in `w_covs` are assumed to have a linear effect on the
hazard. Weights are estimated for all individuals in `data_raw`, even if
they will not be used in the analysis as they do not meet the
landmarking requirements. If an individual enters an absorbing state
prior to `t`, we estimate the probability of being censored before the
time of entry into the absorbing state, rather than at `t`. Details on
all the above this are provided in vignette *overview*.

## References

Hernan M, Robins J (2020). “12.2 Estimating IP weights via modeling.” In
*Causal Inference: What If*, chapter 12.2. Chapman Hall/CRC, Boca Raton.

## Examples

``` r
# Estimate inverse probability of censoring weights for individual in cohort ebmtcal.
# Specifically the probability of being uncensored at t = 1826 days.
# Weights are estimated using a model fitted in all individuals uncensored at time s = 0.
weights_manual <-
calc_weights(data_ms = msebmtcal,
  data_raw = ebmtcal,
  covs = c("year", "agecl", "proph", "match"),
  t = 1826,
  s = 0,
  landmark_type = "state",
  j = 1)

 str(weights_manual)
#> 'data.frame':    2279 obs. of  2 variables:
#>  $ id  : int  1 2 3 4 5 6 7 8 9 10 ...
#>  $ ipcw: num  NA 1.14 NA 1.01 1.03 ...
```
