###
### Tests for calibration curves produced using pseudo-values (calib_type = 'pv')
###

### Run tests for when curve_type = "loess" and CI_type = "bootstrap"_
test_that("check calib_pv output, (j = 1, s = 0), curve_type = loess, CI_type = bootstrap", {

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
  dat_calib_pv_1 <- calib_msm(data_ms = msebmtcal,
                              data_raw = ebmtcal,
                              j = 1,
                              s = 0,
                              t = 1826,
                              tp_pred = tp_pred,
                              calib_type = 'pv',
                              curve_type = "loess",
                              tp_pred_plot = NULL, transitions_out = NULL)

  expect_equal(class(dat_calib_pv_1), c("calib_pv", "calib_msm"))
  expect_equal(dat_calib_pv_1[["metadata"]][["curve_type"]], "loess")
  expect_equal(ncol(dat_calib_pv_1[["plotdata"]][[1]]), 4)
  expect_no_error(summary(dat_calib_pv_1))

  ## Check same results when just calculating pseudo-values for first three individuals
  dat_calib_pv_ids_1 <- calib_msm(data_ms = msebmtcal,
                                  data_raw = ebmtcal,
                                  j = 1,
                                  s = 0,
                                  t = 1826,
                                  tp_pred = tp_pred,
                                  calib_type = 'pv',
                                  pv_ids = 1:3,
                                  tp_pred_plot = NULL, transitions_out = NULL)

  expect_equal(dat_calib_pv_1[["plotdata"]][[1]][1:3, "pv"], dat_calib_pv_ids_1[[1]][,2])
  expect_equal(dat_calib_pv_1[["plotdata"]][[6]][1:3, "pv"], dat_calib_pv_ids_1[[1]][,7])

  ## Calculate observed event probabilities with a confidence interval using bootstrapping and transitions_out = NULL
  expect_warning(calib_msm(data_ms = msebmtcal,
                           data_raw = ebmtcal,
                           j = 1,
                           s = 0,
                           t = 1826,
                           tp_pred = tp_pred,
                           calib_type = 'pv',
                           curve_type = "loess",
                           CI = 95, CI_type = "bootstrap", CI_R_boot = 3,
                           tp_pred_plot = NULL, transitions_out = c(1)))

  dat_calib_pv_4 <- suppressWarnings(calib_msm(data_ms = msebmtcal,
                                               data_raw = ebmtcal,
                                               j = 1,
                                               s = 0,
                                               t = 1826,
                                               tp_pred = tp_pred,
                                               calib_type = 'pv',
                                               curve_type = "loess",
                                               CI = 95, CI_type = "bootstrap", CI_R_boot = 3,
                                               tp_pred_plot = NULL, transitions_out = c(1,2)))

  expect_equal(class(dat_calib_pv_4), c("calib_pv", "calib_msm"))
  expect_equal(ncol(dat_calib_pv_4[["plotdata"]][[1]]), 5)

  expect_equal(dat_calib_pv_1[["plotdata"]][[1]]$obs, dat_calib_pv_4[["plotdata"]][[1]]$obs)
  expect_equal(dat_calib_pv_1[["plotdata"]][[1]]$pred, dat_calib_pv_4[["plotdata"]][[1]]$pred)

  expect_equal(dat_calib_pv_1[["plotdata"]][[2]]$obs, dat_calib_pv_4[["plotdata"]][[2]]$obs)
  expect_equal(dat_calib_pv_1[["plotdata"]][[2]]$pred, dat_calib_pv_4[["plotdata"]][[2]]$pred)

  expect_no_error(summary(dat_calib_pv_4))

  ## Calculate observed event probabilities with a confidence interval using bootstrapping, transitions_out = NULL and defining tp_pred_plot manually

  ### Create landmark ids and extract tp_pred_plot correct
  id_lmk <- 1:50
  tp_pred_plot <- tps0 |>
    dplyr::filter(id %in% id_lmk) |>
    dplyr::filter(j == 1) |>
    dplyr::select(any_of(paste("pstate", 1:6, sep = "")))

  ## No confidence interval
  dat_calib_pv_9 <- suppressWarnings(calib_msm(data_ms = msebmtcal,
                                               data_raw = ebmtcal,
                                               j = 1,
                                               s = 0,
                                               t = 1826,
                                               tp_pred = tp_pred,
                                               calib_type = 'pv',
                                               curve_type = "loess",
                                               tp_pred_plot = tp_pred_plot, transitions_out = NULL))

  ## Should be one less column in plotdata (no patient ids)
  expect_equal(class(dat_calib_pv_9), c("calib_pv", "calib_msm"))
  expect_equal(ncol(dat_calib_pv_9[["plotdata"]][[1]]), 3)
  expect_equal(nrow(dat_calib_pv_9[["plotdata"]][[1]]), 50)
  expect_equal(dat_calib_pv_1[["plotdata"]][[1]]$obs, dat_calib_pv_9[["plotdata"]][[1]]$obs)
  expect_equal(dat_calib_pv_1[["plotdata"]][[1]]$pred, dat_calib_pv_9[["plotdata"]][[1]]$pred)

  expect_no_error(summary(dat_calib_pv_9))

  ## With confidence interval
  dat_calib_pv_10 <- suppressWarnings(calib_msm(data_ms = msebmtcal,
                                                data_raw = ebmtcal,
                                                j = 1,
                                                s = 0,
                                                t = 1826,
                                                tp_pred = tp_pred,
                                                calib_type = 'pv',
                                                curve_type = "loess",
                                                CI = 95, CI_type = "bootstrap", CI_R_boot = 3,
                                                tp_pred_plot = tp_pred_plot, transitions_out = NULL))

  expect_equal(class(dat_calib_pv_10), c("calib_pv", "calib_msm"))
  expect_equal(ncol(dat_calib_pv_10[["plotdata"]][[1]]), 4)
  expect_equal(nrow(dat_calib_pv_10[["plotdata"]][[1]]), 50)
  expect_equal(dat_calib_pv_1[["plotdata"]][[1]]$obs, dat_calib_pv_10[["plotdata"]][[1]]$obs)
  expect_equal(dat_calib_pv_1[["plotdata"]][[1]]$pred, dat_calib_pv_10[["plotdata"]][[1]]$pred)

  expect_no_error(summary(dat_calib_pv_10))

})

### Run tests for when curve_type = "loess" and CI_type = "bootstrap"_
test_that("check calib_pv output, (j = 1, s = 0), curve_type = loess, CI_type = parametric", {

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
  dat_calib_pv_1 <- calib_msm(data_ms = msebmtcal,
                              data_raw = ebmtcal,
                              j = 1,
                              s = 0,
                              t = 1826,
                              tp_pred = tp_pred,
                              calib_type = 'pv',
                              curve_type = "loess",
                              tp_pred_plot = NULL, transitions_out = NULL)

  expect_equal(dat_calib_pv_1[["metadata"]][["curve_type"]], "loess")
  expect_equal(ncol(dat_calib_pv_1[["plotdata"]][[1]]), 4)
  expect_no_error(summary(dat_calib_pv_1))

  ## Calculate observed event probabilities with a confidence interval using parametric approach
  dat_calib_pv_5 <- calib_msm(data_ms = msebmtcal,
                              data_raw = ebmtcal,
                              j = 1,
                              s = 0,
                              t = 1826,
                              tp_pred = tp_pred,
                              calib_type = 'pv',
                              curve_type = "loess",
                              CI = 95, CI_type = "parametric",
                              tp_pred_plot = NULL, transitions_out = c(1,2))

  expect_equal(ncol(dat_calib_pv_5[["plotdata"]][[1]]), 6)

  expect_equal(dat_calib_pv_1[["plotdata"]][[1]]$obs, dat_calib_pv_5[["plotdata"]][[1]]$obs)
  expect_equal(dat_calib_pv_1[["plotdata"]][[1]]$pred, dat_calib_pv_5[["plotdata"]][[1]]$pred)

  expect_equal(dat_calib_pv_1[["plotdata"]][[2]]$obs, dat_calib_pv_5[["plotdata"]][[2]]$obs)
  expect_equal(dat_calib_pv_1[["plotdata"]][[2]]$pred, dat_calib_pv_5[["plotdata"]][[2]]$pred)

  expect_no_error(summary(dat_calib_pv_5))

  ## Calculate observed event probabilities with a confidence interval using bootstrapping, transitions_out = NULL and defining tp_pred_plot manually

  ### Create landmark ids and extract tp_pred_plot correct
  id_lmk <- 1:50
  tp_pred_plot <- tps0 |>
    dplyr::filter(id %in% id_lmk) |>
    dplyr::filter(j == 1) |>
    dplyr::select(any_of(paste("pstate", 1:6, sep = "")))

  ## With confidence interval
  dat_calib_pv_10 <- calib_msm(data_ms = msebmtcal,
                               data_raw = ebmtcal,
                               j = 1,
                               s = 0,
                               t = 1826,
                               tp_pred = tp_pred,
                               calib_type = 'pv',
                               curve_type = "loess",
                               CI = 95, CI_type = "parametric",
                               tp_pred_plot = tp_pred_plot, transitions_out = NULL)

  str(dat_calib_pv_10)
  expect_equal(ncol(dat_calib_pv_10[["plotdata"]][[1]]), 5)
  expect_equal(dat_calib_pv_1[["plotdata"]][[1]]$obs, dat_calib_pv_10[["plotdata"]][[1]]$obs)
  expect_equal(dat_calib_pv_1[["plotdata"]][[1]]$pred, dat_calib_pv_10[["plotdata"]][[1]]$pred)

  expect_no_error(summary(dat_calib_pv_10))

})


### Run tests for when curve_type = "rcs" and CI_type = "bootstrap" (not rerunning all of them for curve_type = rcs)
test_that("check calib_pv output, (j = 1, s = 0), curve_type = rcs, CI_type = bootstrap_", {

  skip_on_cran()

  ## Reduce to 150 individuals
  # Extract the predicted transition probabilities out of state j = 1 for first 150 individuals
  tp_pred <- tps0 |>
    dplyr::filter(id %in% 1:150) |>
    dplyr::filter(j == 1) |>
    dplyr::select(any_of(paste("pstate", 1:6, sep = "")))
  # Reduce ebmtcal to first 150 individuals
  ebmtcal <- ebmtcal |> dplyr::filter(id %in% 1:150)
  # Reduce msebmtcal_cmprsk to first 150 individuals
  msebmtcal <- msebmtcal |> dplyr::filter(id %in% 1:150)

  ## Calculate observed event probabilities using transitions_out = NULL
  dat_calib_pv_1 <- calib_msm(data_ms = msebmtcal,
                              data_raw = ebmtcal,
                              j = 1,
                              s = 0,
                              t = 1826,
                              tp_pred = tp_pred,
                              calib_type = 'pv',
                              curve_type = "rcs",
                              tp_pred_plot = NULL, transitions_out = c(1))

  expect_equal(dat_calib_pv_1[["metadata"]][["curve_type"]], "rcs")
  expect_equal(ncol(dat_calib_pv_1[["plotdata"]][[1]]), 4)
  expect_no_error(summary(dat_calib_pv_1))

  ## Calculate observed event probabilities with a confidence interval using bootstrapping
  dat_calib_pv_4 <- suppressWarnings(calib_msm(data_ms = msebmtcal,
                                               data_raw = ebmtcal,
                                               j = 1,
                                               s = 0,
                                               t = 1826,
                                               tp_pred = tp_pred,
                                               calib_type = 'pv',
                                               curve_type = "rcs",
                                               CI = 95, CI_type = "bootstrap", CI_R_boot = 3,
                                               tp_pred_plot = NULL, transitions_out = c(1)))

  expect_equal(ncol(dat_calib_pv_4[["plotdata"]][[1]]), 5)

  expect_equal(dat_calib_pv_1[["plotdata"]][[1]]$obs, dat_calib_pv_4[["plotdata"]][[1]]$obs)
  expect_equal(dat_calib_pv_1[["plotdata"]][[1]]$pred, dat_calib_pv_4[["plotdata"]][[1]]$pred)

  expect_no_error(summary(dat_calib_pv_4))

})

### Run tests for when curve_type = "rcs" and CI_type = "parametric" (not rerunning all of them for curve_type = rcs)
test_that("check calib_pv output, (j = 1, s = 0), curve_type = rcs, CI_type = bootstrap_", {

  skip_on_cran()

  ## Reduce to 150 individuals
  # Extract the predicted transition probabilities out of state j = 1 for first 150 individuals
  tp_pred <- tps0 |>
    dplyr::filter(id %in% 1:150) |>
    dplyr::filter(j == 1) |>
    dplyr::select(any_of(paste("pstate", 1:6, sep = "")))
  # Reduce ebmtcal to first 150 individuals
  ebmtcal <- ebmtcal |> dplyr::filter(id %in% 1:150)
  # Reduce msebmtcal_cmprsk to first 150 individuals
  msebmtcal <- msebmtcal |> dplyr::filter(id %in% 1:150)

  ## Calculate observed event probabilities using transitions_out = NULL
  dat_calib_pv_1 <- calib_msm(data_ms = msebmtcal,
                              data_raw = ebmtcal,
                              j = 1,
                              s = 0,
                              t = 1826,
                              tp_pred = tp_pred,
                              calib_type = 'pv',
                              curve_type = "rcs",
                              tp_pred_plot = NULL, transitions_out = c(1))

  expect_equal(dat_calib_pv_1[["metadata"]][["curve_type"]], "rcs")
  expect_equal(ncol(dat_calib_pv_1[["plotdata"]][[1]]), 4)
  expect_no_error(summary(dat_calib_pv_1))

  ## Calculate observed event probabilities with a confidence interval using parametric approach
  dat_calib_pv_4 <- calib_msm(data_ms = msebmtcal,
                              data_raw = ebmtcal,
                              j = 1,
                              s = 0,
                              t = 1826,
                              tp_pred = tp_pred,
                              calib_type = 'pv',
                              curve_type = "rcs",
                              CI = 95, CI_type = "parametric",
                              tp_pred_plot = NULL, transitions_out = c(1))

  expect_equal(ncol(dat_calib_pv_4[["plotdata"]][[1]]), 6)

  expect_equal(dat_calib_pv_1[["plotdata"]][[1]]$obs, dat_calib_pv_4[["plotdata"]][[1]]$obs)
  expect_equal(dat_calib_pv_1[["plotdata"]][[1]]$pred, dat_calib_pv_4[["plotdata"]][[1]]$pred)

  expect_no_error(summary(dat_calib_pv_4))

})


### Add some tests for when each of group_vars and pv_n_pctls are specified
test_that("check calib_pv output, (j = 1, s = 0), groups_vars and pv_n_pctls specified", {

  skip_on_cran()

  ## Reduce to 50 individuals
  # Extract the predicted transition probabilities out of state j = 1 for first 100 individuals
  tp_pred <- tps0 |>
    dplyr::filter(id %in% 1:50) |>
    dplyr::filter(j == 1) |>
    dplyr::select(any_of(paste("pstate", 1:6, sep = "")))
  # Reduce ebmtcal to first 50 individuals
  ebmtcal <- ebmtcal |> dplyr::filter(id %in% 1:50)
  # Reduce msebmtcal_cmprsk to first 100 individuals
  msebmtcal <- msebmtcal |> dplyr::filter(id %in% 1:50)

  ## Calculate observed event probabilities when both pv_group_vars and pv_n_pctls are specified
  dat_calib_pv_1 <- calib_msm(data_ms = msebmtcal,
                              data_raw = ebmtcal,
                              j = 1,
                              s = 0,
                              t = 1826,
                              tp_pred = tp_pred,
                              calib_type = 'pv',
                              curve_type = "loess",
                              loess_span = 1,
                              loess_degree = 1,
                              pv_group_vars = c("year"),
                              pv_n_pctls = 2,
                              tp_pred_plot = NULL, transitions_out = NULL)

  expect_equal(ncol(dat_calib_pv_1[["plotdata"]][[1]]), 4)
  expect_equal(length(dat_calib_pv_1[["plotdata"]]), 6)

  ## Check same results when just calculating pseudo-values for first three individuals
  dat_calib_pv_ids_1 <- calib_msm(data_ms = msebmtcal,
                                  data_raw = ebmtcal,
                                  j = 1,
                                  s = 0,
                                  t = 1826,
                                  tp_pred = tp_pred,
                                  calib_type = 'pv',
                                  pv_group_vars = c("year"),
                                  pv_n_pctls = 2,
                                  pv_ids = 1:3,
                                  tp_pred_plot = NULL, transitions_out = NULL)

  expect_equal(dat_calib_pv_1[["plotdata"]][[1]][1:3, "pv"], dat_calib_pv_ids_1[[1]][,2])
  expect_equal(dat_calib_pv_1[["plotdata"]][[6]][1:3, "pv"], dat_calib_pv_ids_1[[1]][,7])

  ## Check same results when just calculating pseudo-values for first three individuals, but specify transitions 1 and 6
  dat_calib_pv_ids_1_tout <- calib_msm(data_ms = msebmtcal,
                                       data_raw = ebmtcal,
                                       j = 1,
                                       s = 0,
                                       t = 1826,
                                       tp_pred = tp_pred,
                                       calib_type = 'pv',
                                       pv_group_vars = c("year"),
                                       pv_n_pctls = 2,
                                       pv_ids = 1:3,
                                       tp_pred_plot = NULL, transitions_out = c(1,6))

  expect_equal(dat_calib_pv_ids_1_tout[[1]][,2], dat_calib_pv_ids_1[[1]][,2])
  expect_equal(dat_calib_pv_ids_1_tout[[1]][,3], dat_calib_pv_ids_1[[1]][,7])
  expect_equal(ncol(dat_calib_pv_ids_1_tout[["plotdata"]]), 3)

  ## Calculate observed event probabilities for pv_group_vars
  dat_calib_pv_2 <- calib_msm(data_ms = msebmtcal,
                              data_raw = ebmtcal,
                              j = 1,
                              s = 0,
                              t = 1826,
                              tp_pred = tp_pred,
                              calib_type = 'pv',
                              curve_type = "loess",
                              loess_span = 1,
                              loess_degree = 1,
                              pv_group_vars = c("year"),
                              tp_pred_plot = NULL, transitions_out = NULL)

  expect_equal(ncol(dat_calib_pv_2[["plotdata"]][[1]]), 4)
  expect_equal(length(dat_calib_pv_2[["plotdata"]]), 6)

  ## Check same results when just calculating pseudo-values for first three individuals
  dat_calib_pv_ids_2 <- calib_msm(data_ms = msebmtcal,
                                  data_raw = ebmtcal,
                                  j = 1,
                                  s = 0,
                                  t = 1826,
                                  tp_pred = tp_pred,
                                  calib_type = 'pv',
                                  pv_group_vars = c("year"),
                                  pv_ids = 1:3,
                                  tp_pred_plot = NULL, transitions_out = NULL)

  expect_equal(dat_calib_pv_2[["plotdata"]][[1]][1:3, "pv"], dat_calib_pv_ids_2[[1]][,2])
  expect_equal(dat_calib_pv_2[["plotdata"]][[6]][1:3, "pv"], dat_calib_pv_ids_2[[1]][,7])

  ## No need to test for transitions_out when pv_n_pctls not specified, because there are no computational gains and
  ## pseudo-values are just calculated for all states anyway_

  ## Calculate observed event probabilities for pv_n_pctls
  dat_calib_pv_3 <- calib_msm(data_ms = msebmtcal,
                              data_raw = ebmtcal,
                              j = 1,
                              s = 0,
                              t = 1826,
                              tp_pred = tp_pred,
                              calib_type = 'pv',
                              curve_type = "loess",
                              loess_span = 1,
                              loess_degree = 1,
                              pv_n_pctls = 2,
                              tp_pred_plot = NULL, transitions_out = NULL)

  expect_equal(ncol(dat_calib_pv_3[["plotdata"]][[1]]), 4)
  expect_equal(length(dat_calib_pv_3[["plotdata"]]), 6)

  ## Check same results when just calculating pseudo-values for first three individuals
  dat_calib_pv_ids_3 <- calib_msm(data_ms = msebmtcal,
                                  data_raw = ebmtcal,
                                  j = 1,
                                  s = 0,
                                  t = 1826,
                                  tp_pred = tp_pred,
                                  calib_type = 'pv',
                                  pv_n_pctls = 2,
                                  pv_ids = 1:3,
                                  tp_pred_plot = NULL, transitions_out = NULL)

  expect_equal(dat_calib_pv_3[["plotdata"]][[1]][1:3, "pv"], dat_calib_pv_ids_3[[1]][,2])
  expect_equal(dat_calib_pv_3[["plotdata"]][[6]][1:3, "pv"], dat_calib_pv_ids_3[[1]][,7])

  ## Check same results when just calculating pseudo-values for first three individuals, but specify transitions 1 and 6
  dat_calib_pv_ids_3_tout <- calib_msm(data_ms = msebmtcal,
                                       data_raw = ebmtcal,
                                       j = 1,
                                       s = 0,
                                       t = 1826,
                                       tp_pred = tp_pred,
                                       calib_type = 'pv',
                                       pv_n_pctls = 2,
                                       pv_ids = 1:3,
                                       tp_pred_plot = NULL, transitions_out = c(1,6))

  expect_equal(dat_calib_pv_ids_3_tout[[1]][,2], dat_calib_pv_ids_3[[1]][,2])
  expect_equal(dat_calib_pv_ids_3_tout[[1]][,3], dat_calib_pv_ids_3[[1]][,7])
  expect_equal(ncol(dat_calib_pv_ids_3_tout[["plotdata"]]), 3)

})



### Add some tests where we expect errors, if requesting things that aren't possible
test_that("check calib_pv output, (j = 1, s = 0), cause errors", {

  skip_on_cran()

  ## Reduce to 50 individuals
  # Extract the predicted transition probabilities out of state j = 1 for first 100 individuals
  tp_pred <- tps0 |>
    dplyr::filter(id %in% 1:50) |>
    dplyr::filter(j == 1) |>
    dplyr::select(any_of(paste("pstate", 1:6, sep = "")))
  # Reduce ebmtcal to first 50 individuals
  ebmtcal <- ebmtcal |> dplyr::filter(id %in% 1:50)
  # Reduce msebmtcal_cmprsk to first 100 individuals
  msebmtcal <- msebmtcal |> dplyr::filter(id %in% 1:50)

  ## Request bootstrap confidence interval and don't give number of bootstrap replicates (for either rcs or parametric)
  expect_error(calib_msm(data_ms = msebmtcal,
                         data_raw = ebmtcal,
                         j = 1,
                         s = 0,
                         t = 1826,
                         tp_pred = tp_pred,
                         calib_type = 'pv',
                         curve_type = "loess",
                         CI = 95,
                         CI_type = "bootstrap",
                         tp_pred_plot = NULL, transitions_out = NULL))

  expect_error(calib_msm(data_ms = msebmtcal,
                         data_raw = ebmtcal,
                         j = 1,
                         s = 0,
                         t = 1826,
                         tp_pred = tp_pred,
                         calib_type = 'pv',
                         curve_type = "rcs",
                         CI = 95,
                         CI_type = "bootstrap",
                         tp_pred_plot = NULL, transitions_out = NULL))

})


test_that("check calib_pv output, (j = 3, s = 100), pv_group_vars defined", {

  skip_on_cran()

  ## Extract relevant predicted risks from tps100
  tp_pred <- dplyr::select(dplyr::filter(tps100, j == 3), any_of(paste("pstate", 1:6, sep = "")))

  ## Calculate observed event probabilities
  dat_calib_pv <-
    calib_msm(data_ms = msebmtcal,
              data_raw = ebmtcal,
              j=3,
              s=100,
              t = 1826,
              tp_pred = tp_pred, calib_type = 'pv',
              curve_type = "rcs",
              rcs_nk = 3,
              pv_group_vars = c("year"))

  expect_type(dat_calib_pv, "list")
  expect_equal(class(dat_calib_pv), c("calib_pv", "calib_msm"))
  expect_length(dat_calib_pv[["plotdata"]], 4)
  expect_length(dat_calib_pv[["plotdata"]][["state3"]]$id, 413)
  expect_length(dat_calib_pv[["plotdata"]][["state6"]]$id, 413)
  expect_error(dat_calib_pv[["plotdata"]][[6]])
  expect_false(dat_calib_pv[["metadata"]]$CI)


})


test_that("check calib_pv output, (j = 3, s = 100), pv_n_pctls defined", {

  skip_on_cran()

  ## Extract relevant predicted risks from tps100
  tp_pred <- dplyr::select(dplyr::filter(tps100, j == 3), any_of(paste("pstate", 1:6, sep = "")))

  ## Calculate observed event probabilities
  dat_calib_pv <-
    calib_msm(data_ms = msebmtcal,
              data_raw = ebmtcal,
              j=3,
              s=100,
              t = 1826,
              tp_pred = tp_pred, calib_type = 'pv',
              curve_type = "rcs",
              rcs_nk = 3,
              pv_n_pctls = 2)

  expect_type(dat_calib_pv, "list")
  expect_equal(class(dat_calib_pv), c("calib_pv", "calib_msm"))
  expect_length(dat_calib_pv[["plotdata"]], 4)
  expect_length(dat_calib_pv[["plotdata"]][["state3"]]$id, 413)
  expect_length(dat_calib_pv[["plotdata"]][["state6"]]$id, 413)
  expect_error(dat_calib_pv[["plotdata"]][[6]])
  expect_false(dat_calib_pv[["metadata"]]$CI)

})


test_that("check calib_pv output, (j = 3, s = 100), pv_group_vars and pv_n_pctls defined", {

  skip_on_cran()

  ## Extract relevant predicted risks from tps100
  tp_pred <- dplyr::select(dplyr::filter(tps100, j == 3), any_of(paste("pstate", 1:6, sep = "")))

  ## Calculate observed event probabilities
  dat_calib_pv <-
    calib_msm(data_ms = msebmtcal,
              data_raw = ebmtcal,
              j=3,
              s=100,
              t = 1826,
              tp_pred = tp_pred, calib_type = 'pv',
              curve_type = "rcs",
              rcs_nk = 3,
              pv_group_vars = c("year"),
              pv_n_pctls = 2)

  expect_type(dat_calib_pv, "list")
  expect_equal(class(dat_calib_pv), c("calib_pv", "calib_msm"))
  expect_length(dat_calib_pv[["plotdata"]], 4)
  expect_length(dat_calib_pv[["plotdata"]][["state3"]]$id, 413)
  expect_length(dat_calib_pv[["plotdata"]][["state6"]]$id, 413)
  expect_error(dat_calib_pv[["plotdata"]][[6]])
  expect_false(dat_calib_pv[["metadata"]]$CI)

})


test_that("check calib_pv output, (j = 1, s = 0), pv_precalc", {

  skip_on_cran()

  ## Extract relevant predicted risks from tps100
  tp_pred <- dplyr::select(dplyr::filter(tps0, j == 1), any_of(paste("pstate", 1:6, sep = "")))

  ## Define pv_precalc to be the estimated predicted probabilities
  pv_precalc <- tp_pred

  ## Calculate observed event probabilities
  dat_calib_pv <-
    calib_msm(data_ms = msebmtcal,
              data_raw = ebmtcal,
              j = 1,
              s = 0,
              t = 1826,
              tp_pred = tp_pred,
              calib_type = 'pv',
              pv_precalc = tp_pred,
              curve_type = "rcs",
              rcs_nk = 3)

  expect_type(dat_calib_pv, "list")
  expect_equal(class(dat_calib_pv), c("calib_pv", "calib_msm"))
  expect_length(dat_calib_pv[["plotdata"]], 6)
  expect_length(dat_calib_pv[["plotdata"]][["state3"]]$id, 2279)
  expect_length(dat_calib_pv[["plotdata"]][["state6"]]$id, 2279)
  expect_false(dat_calib_pv[["metadata"]]$CI)

})

