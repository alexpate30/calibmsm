###
### Tests for calibration curves produced using pseudo-values (calib.type = 'AJ')
###

### Run tests for pv_n_pctls = NULL and pv_group_vars = NULL
test_that("check calib_aj, pv_n_pctls = NULL and pv_group_vars = NULL", {

  ## Reduce to 50 individuals
  # Extract the predicted transition probabilities out of state j = 1 for first 50 individuals
  tp_pred <- tps0 |>
    dplyr::filter(id %in% 1:50) |>
    dplyr::filter(j == 1) |>
    dplyr::select(any_of(paste("pstate", 1:6, sep = "")))
  # Reduce ebmtcal to first 50 individuals
  ebmtcal <- ebmtcal |> dplyr::filter(id %in% 1:50)
  # Reduce msebmtcal_cmprsk to first 100 individuals
  msebmtcal <- msebmtcal |> dplyr::filter(id %in% 1:50)

  ## Calculate observed event probabilities using transitions_out = NULL
  dat_calib_aj_1 <- suppressWarnings(calib_msm(data_ms = msebmtcal,
                                               data_raw = ebmtcal,
                                               j = 1,
                                               s = 0,
                                               t = 1826,
                                               tp_pred = tp_pred,
                                               calib_type = 'aj',
                                               tp_pred_plot = NULL, transitions_out = NULL))

  expect_equal(class(dat_calib_aj_1), c("calib_aj", "calib_msm"))
  expect_equal(length(dat_calib_aj_1[["mean"]]), 6)

  ## Calculate observed event probabilities using transitions_out = NULL
  dat_calib_aj_CI_1 <- suppressWarnings(calib_msm(data_ms = msebmtcal,
                                                  data_raw = ebmtcal,
                                                  j = 1,
                                                  s = 0,
                                                  t = 1826,
                                                  tp_pred = tp_pred,
                                                  calib_type = 'aj',
                                                  CI = 95,
                                                  CI_R_boot = 10,
                                                  tp_pred_plot = NULL, transitions_out = NULL))

  expect_equal(class(dat_calib_aj_CI_1), c("calib_aj", "calib_msm"))
  expect_equal(length(dat_calib_aj_CI_1[["mean"]]), 6)
  expect_equal(length(dat_calib_aj_CI_1[["mean"]][[1]]), 3)
  expect_equal(as.numeric(dat_calib_aj_1[["mean"]][1]), as.numeric(dat_calib_aj_CI_1[["mean"]][[1]][1]))
  expect_equal(as.numeric(dat_calib_aj_1[["mean"]][6]), as.numeric(dat_calib_aj_CI_1[["mean"]][[6]][1]))

})


### Run tets pv_n_pctls specified
test_that("check calib_pv output, pv_n_pctls specified", {

  ## Reduce to 50 individuals
  # Extract the predicted transition probabilities out of state j = 1 for first 50 individuals
  tp_pred <- tps0 |>
    dplyr::filter(id %in% 1:50) |>
    dplyr::filter(j == 1) |>
    dplyr::select(any_of(paste("pstate", 1:6, sep = "")))
  # Reduce ebmtcal to first 50 individuals
  ebmtcal <- ebmtcal |> dplyr::filter(id %in% 1:50)
  # Reduce msebmtcal_cmprsk to first 100 individuals
  msebmtcal <- msebmtcal |> dplyr::filter(id %in% 1:50)

  ## Calculate observed event probabilities using transitions_out = NULL
  dat_calib_aj_1 <- suppressWarnings(calib_msm(data_ms = msebmtcal,
                                               data_raw = ebmtcal,
                                               j = 1,
                                               s = 0,
                                               t = 1826,
                                               tp_pred = tp_pred,
                                               calib_type = 'aj',
                                               pv_n_pctls = 2,
                                               tp_pred_plot = NULL, transitions_out = NULL))

  expect_equal(class(dat_calib_aj_1), c("calib_aj", "calib_msm"))
  expect_equal(length(dat_calib_aj_1[["mean"]]), 6)

  ## Calculate observed event probabilities using transitions_out = NULL
  dat_calib_aj_CI_1 <- suppressWarnings(calib_msm(data_ms = msebmtcal,
                                                  data_raw = ebmtcal,
                                                  j = 1,
                                                  s = 0,
                                                  t = 1826,
                                                  tp_pred = tp_pred,
                                                  calib_type = 'aj',
                                                  pv_n_pctls = 2,
                                                  CI = 95,
                                                  CI_R_boot = 10,
                                                  tp_pred_plot = NULL, transitions_out = NULL))

  expect_equal(class(dat_calib_aj_CI_1), c("calib_aj", "calib_msm"))
  expect_equal(length(dat_calib_aj_CI_1[["mean"]]), 6)
  expect_equal(length(dat_calib_aj_CI_1[["mean"]][[1]]), 3)
  expect_equal(as.numeric(dat_calib_aj_1[["mean"]][1]), as.numeric(dat_calib_aj_CI_1[["mean"]][[1]][1]))
  expect_equal(as.numeric(dat_calib_aj_1[["mean"]][6]), as.numeric(dat_calib_aj_CI_1[["mean"]][[6]][1]))

})


### Run tests pv_group_vars specified
test_that("check calib_pv output, pv_group_vars specified", {

  skip_on_cran()

  ## Reduce to 50 individuals
  # Extract the predicted transition probabilities out of state j = 1 for first 50 individuals
  tp_pred <- tps0 |>
    dplyr::filter(id %in% 1:50) |>
    dplyr::filter(j == 1) |>
    dplyr::select(any_of(paste("pstate", 1:6, sep = "")))
  # Reduce ebmtcal to first 50 individuals
  ebmtcal <- ebmtcal |> dplyr::filter(id %in% 1:50)
  # Reduce msebmtcal_cmprsk to first 100 individuals
  msebmtcal <- msebmtcal |> dplyr::filter(id %in% 1:50)

  ## Calculate observed event probabilities using transitions_out = NULL
  dat_calib_aj_1 <- suppressWarnings(calib_msm(data_ms = msebmtcal,
                                               data_raw = ebmtcal,
                                               j = 1,
                                               s = 0,
                                               t = 1826,
                                               tp_pred = tp_pred,
                                               calib_type = 'aj',
                                               pv_group_vars = c("year"),
                                               tp_pred_plot = NULL, transitions_out = NULL))

  expect_equal(class(dat_calib_aj_1), c("calib_aj", "calib_msm"))
  expect_equal(length(dat_calib_aj_1[["mean"]]), 6)

  ## Calculate observed event probabilities using transitions_out = NULL
  dat_calib_aj_CI_1 <- suppressWarnings(calib_msm(data_ms = msebmtcal,
                                                  data_raw = ebmtcal,
                                                  j = 1,
                                                  s = 0,
                                                  t = 1826,
                                                  tp_pred = tp_pred,
                                                  calib_type = 'aj',
                                                  pv_group_vars = c("year"),
                                                  CI = 95,
                                                  CI_R_boot = 10,
                                                  tp_pred_plot = NULL, transitions_out = NULL))

  expect_equal(class(dat_calib_aj_CI_1), c("calib_aj", "calib_msm"))
  expect_equal(length(dat_calib_aj_CI_1[["mean"]]), 6)
  expect_equal(length(dat_calib_aj_CI_1[["mean"]][[1]]), 3)
  expect_equal(as.numeric(dat_calib_aj_1[["mean"]][1]), as.numeric(dat_calib_aj_CI_1[["mean"]][[1]][1]))
  expect_equal(as.numeric(dat_calib_aj_1[["mean"]][6]), as.numeric(dat_calib_aj_CI_1[["mean"]][[6]][1]))

})


### Run tests pv_group_vars and pv_n_pctls specified
test_that("check calib_pv output, pv_group_vars and pv_n_pctls specified ", {

  skip_on_cran()

  ## Set seed
  set.seed(101)

  ## Reduce to 50 individuals
  # Extract the predicted transition probabilities out of state j = 1 for first 50 individuals
  tp_pred <- tps0 |>
    dplyr::filter(id %in% 1:100) |>
    dplyr::filter(j == 1) |>
    dplyr::select(any_of(paste("pstate", 1:6, sep = "")))
  # Reduce ebmtcal to first 50 individuals
  ebmtcal <- ebmtcal |> dplyr::filter(id %in% 1:100)
  # Reduce msebmtcal_cmprsk to first 100 individuals
  msebmtcal <- msebmtcal |> dplyr::filter(id %in% 1:100)

  ## Calculate observed event probabilities using transitions_out = NULL
  dat_calib_aj_1 <- suppressWarnings(calib_msm(data_ms = msebmtcal,
                                               data_raw = ebmtcal,
                                               j = 1,
                                               s = 0,
                                               t = 1826,
                                               tp_pred = tp_pred,
                                               calib_type = 'aj',
                                               pv_n_pctls = 2,
                                               pv_group_vars = c("year"),
                                               tp_pred_plot = NULL, transitions_out = NULL))

  expect_equal(class(dat_calib_aj_1), c("calib_aj", "calib_msm"))
  expect_equal(length(dat_calib_aj_1[["mean"]]), 6)

  ## Calculate observed event probabilities using transitions_out = NULL
  dat_calib_aj_CI_1 <- suppressWarnings(calib_msm(data_ms = msebmtcal,
                                                  data_raw = ebmtcal,
                                                  j = 1,
                                                  s = 0,
                                                  t = 1826,
                                                  tp_pred = tp_pred,
                                                  calib_type = 'aj',
                                                  pv_n_pctls = 2,
                                                  pv_group_vars = c("year"),
                                                  CI = 95,
                                                  CI_R_boot = 10,
                                                  tp_pred_plot = NULL, transitions_out = NULL))

  expect_equal(class(dat_calib_aj_CI_1), c("calib_aj", "calib_msm"))
  expect_equal(length(dat_calib_aj_CI_1[["mean"]]), 6)
  expect_equal(length(dat_calib_aj_CI_1[["mean"]][[1]]), 3)
  expect_equal(as.numeric(dat_calib_aj_1[["mean"]][1]), as.numeric(dat_calib_aj_CI_1[["mean"]][[1]][1]))
  expect_equal(as.numeric(dat_calib_aj_1[["mean"]][6]), as.numeric(dat_calib_aj_CI_1[["mean"]][[6]][1]))

})

### Finish tests
