### Run tests for when curve.type = "loess".
test_that("check calib_pv output, (j = 1, s = 0), curve.type = loess", {

  ## Reduce to 50 individuals
  # Extract the predicted transition probabilities out of state j = 1 for first 50 individuals
  tp.pred <- tps0 |>
    dplyr::filter(id %in% 1:50) |>
    dplyr::filter(j == 1) |>
    dplyr::select(any_of(paste("pstate", 1:6, sep = "")))
  # Reduce ebmtcal to first 50 individuals
  ebmtcal <- ebmtcal |> dplyr::filter(id %in% 1:50)
  # Reduce msebmtcal.cmprsk to first 100 individuals
  msebmtcal <- msebmtcal |> dplyr::filter(id %in% 1:50)

  ## Calculate observed event probabilities using transitions.out = NULL
  dat.calib.pv.1 <- calib_pv(data.mstate = msebmtcal,
                                  data.raw = ebmtcal,
                                  j = 1,
                                  s = 0,
                                  t = 1826,
                                  tp.pred = tp.pred,
                                  curve.type = "loess",
                                  data.pred.plot = NULL, transitions.out = NULL)

  # ## Calculate observed event probabilities using transitions.out = c(3,4,5,6)
  # dat.calib.pv.2 <- calib_pv(data.mstate = msebmtcal,
  #                                  data.raw = ebmtcal,
  #                                  j = 3,
  #                                  s = 100,
  #                                  t = 1826,
  #                                  tp.pred = tp.pred,
  #                                  curve.type = "loess",
  #                                  group.vars = c("year"),
  #                                  n.pctls = 2,
  #                                  data.pred.plot = NULL, transitions.out = c(3,4,5,6))

  expect_equal(dat.calib.pv.1[["metadata"]][["curve.type"]], "loess")
  expect_equal(ncol(dat.calib.pv.1[["plotdata"]][[1]]), 3)
  expect_no_error(summary(dat.calib.pv.1))
  # expect_equal(ncol(dat.calib.pv.2[["plotdata"]][[1]]), 3)
  # expect_equal(dat.calib.pv.1[["plotdata"]][[1]], dat.calib.pv.2[["plotdata"]][[1]])

  # ## Calculate observed event probabilities using transitions.out = 3
  # dat.calib.pv.3 <- calib_pv(data.mstate = msebmtcal,
  #                                 data.raw = ebmtcal,
  #                                 j = 3,
  #                                 s = 100,
  #                                 t = 1826,
  #                                 tp.pred = tp.pred,
  #                                 curve.type = "loess",
  #                                 group.vars = c("year"),
  #                                 n.pctls = 2,
  #                                 data.pred.plot = NULL, transitions.out = 3)
  #
  # expect_equal(length(dat.calib.pv.3[["plotdata"]]), 1)
  # expect_equal(dat.calib.pv.1[["plotdata"]][[1]], dat.calib.pv.3[["plotdata"]][[1]])

  ## Calculate observed event probabilities with a confidence interval using bootstrapping and transitions.out = NULL
  dat.calib.pv.4 <- calib_pv(data.mstate = msebmtcal,
                                  data.raw = ebmtcal,
                                  j = 1,
                                  s = 0,
                                  t = 1826,
                                  tp.pred = tp.pred,
                                  curve.type = "loess",
                                  CI = 95, CI.type = "bootstrap", CI.R.boot = 3,
                                  data.pred.plot = NULL, transitions.out = c(1,2))

  expect_equal(ncol(dat.calib.pv.4[["plotdata"]][[1]]), 5)

  expect_equal(dat.calib.pv.1[["plotdata"]][[1]]$obs, dat.calib.pv.4[["plotdata"]][[1]]$obs)
  expect_equal(dat.calib.pv.1[["plotdata"]][[1]]$pred, dat.calib.pv.4[["plotdata"]][[1]]$pred)

  expect_equal(dat.calib.pv.1[["plotdata"]][[2]]$obs, dat.calib.pv.4[["plotdata"]][[2]]$obs)
  expect_equal(dat.calib.pv.1[["plotdata"]][[2]]$pred, dat.calib.pv.4[["plotdata"]][[2]]$pred)

  expect_no_error(summary(dat.calib.pv.4))

  ## Calculate observed event probabilities with a confidence interval using parametric approach
  dat.calib.pv.5 <- calib_pv(data.mstate = msebmtcal,
                             data.raw = ebmtcal,
                             j = 1,
                             s = 0,
                             t = 1826,
                             tp.pred = tp.pred,
                             curve.type = "loess",
                             CI = 95, CI.type = "parametric",
                             data.pred.plot = NULL, transitions.out = c(1,2))

  expect_equal(ncol(dat.calib.pv.5[["plotdata"]][[1]]), 5)

  expect_equal(dat.calib.pv.1[["plotdata"]][[1]]$obs, dat.calib.pv.5[["plotdata"]][[1]]$obs)
  expect_equal(dat.calib.pv.1[["plotdata"]][[1]]$pred, dat.calib.pv.5[["plotdata"]][[1]]$pred)

  expect_equal(dat.calib.pv.1[["plotdata"]][[2]]$obs, dat.calib.pv.5[["plotdata"]][[2]]$obs)
  expect_equal(dat.calib.pv.1[["plotdata"]][[2]]$pred, dat.calib.pv.5[["plotdata"]][[2]]$pred)

  expect_no_error(summary(dat.calib.pv.5))

  ## Calculate observed event probabilities with a confidence interval using bootstrapping, transitions.out = NULL and defining data.pred.plot manually

  ### Create landmark ids and extract data.pred.plot correct
  id.lmk <- 1:50
  data.pred.plot.in <- tps0 |>
    dplyr::filter(id %in% id.lmk) |>
    dplyr::filter(j == 1) |>
    dplyr::select(any_of(paste("pstate", 1:6, sep = "")))

  ## No confidence interval
  dat.calib.pv.9 <- calib_pv(data.mstate = msebmtcal,
                                  data.raw = ebmtcal,
                                  j = 1,
                                  s = 0,
                                  t = 1826,
                                  tp.pred = tp.pred,
                                  curve.type = "loess",
                                  data.pred.plot = data.pred.plot.in, transitions.out = NULL)

  ## Should be one less column in plotdata (no patient ids)
  expect_equal(ncol(dat.calib.pv.9[["plotdata"]][[1]]), 2)
  expect_equal(dat.calib.pv.1[["plotdata"]][[1]]$obs, dat.calib.pv.9[["plotdata"]][[1]]$obs)
  expect_equal(dat.calib.pv.1[["plotdata"]][[1]]$pred, dat.calib.pv.9[["plotdata"]][[1]]$pred)

  expect_no_error(summary(dat.calib.pv.9))

  ## With confidence interval
  dat.calib.pv.10 <- calib_pv(data.mstate = msebmtcal,
                                  data.raw = ebmtcal,
                                  j = 1,
                                  s = 0,
                                  t = 1826,
                                  tp.pred = tp.pred,
                                  curve.type = "loess",
                                  CI = 95, CI.type = "bootstrap", CI.R.boot = 3,
                                  data.pred.plot = data.pred.plot.in, transitions.out = NULL)

  expect_equal(ncol(dat.calib.pv.10[["plotdata"]][[1]]), 4)
  expect_equal(dat.calib.pv.1[["plotdata"]][[1]]$obs, dat.calib.pv.10[["plotdata"]][[1]]$obs)
  expect_equal(dat.calib.pv.1[["plotdata"]][[1]]$pred, dat.calib.pv.10[["plotdata"]][[1]]$pred)

  expect_no_error(summary(dat.calib.pv.10))

})

### Run some of these tests for when curve.type = "rcs".

### Run tests for when curve.type = "rcs" and CI.type = "bootstrap".
test_that("check calib_pv output, (j = 1, s = 0), curve.type = rcs, CI.type = bootstrap.", {

  ## Reduce to 150 individuals
  # Extract the predicted transition probabilities out of state j = 1 for first 150 individuals
  tp.pred <- tps0 |>
    dplyr::filter(id %in% 1:150) |>
    dplyr::filter(j == 1) |>
    dplyr::select(any_of(paste("pstate", 1:6, sep = "")))
  # Reduce ebmtcal to first 150 individuals
  ebmtcal <- ebmtcal |> dplyr::filter(id %in% 1:150)
  # Reduce msebmtcal.cmprsk to first 150 individuals
  msebmtcal <- msebmtcal |> dplyr::filter(id %in% 1:150)

  ## Calculate observed event probabilities using transitions.out = NULL
  dat.calib.pv.1 <- calib_pv(data.mstate = msebmtcal,
                                  data.raw = ebmtcal,
                                  j = 1,
                                  s = 0,
                                  t = 1826,
                                  tp.pred = tp.pred,
                                  curve.type = "rcs",
                                  data.pred.plot = NULL, transitions.out = c(1))

  expect_equal(dat.calib.pv.1[["metadata"]][["curve.type"]], "rcs")
  expect_equal(ncol(dat.calib.pv.1[["plotdata"]][[1]]), 3)
  expect_no_error(summary(dat.calib.pv.1))

  ## Calculate observed event probabilities with a confidence interval using bootstrapping
  dat.calib.pv.4 <- calib_pv(data.mstate = msebmtcal,
                                  data.raw = ebmtcal,
                                  j = 1,
                                  s = 0,
                                  t = 1826,
                                  tp.pred = tp.pred,
                                  curve.type = "rcs",
                                  CI = 95, CI.type = "bootstrap", CI.R.boot = 3,
                                  data.pred.plot = NULL, transitions.out = c(1))

  expect_equal(ncol(dat.calib.pv.4[["plotdata"]][[1]]), 5)

  expect_equal(dat.calib.pv.1[["plotdata"]][[1]]$obs, dat.calib.pv.4[["plotdata"]][[1]]$obs)
  expect_equal(dat.calib.pv.1[["plotdata"]][[1]]$pred, dat.calib.pv.4[["plotdata"]][[1]]$pred)

  expect_no_error(summary(dat.calib.pv.4))

  ## Calculate observed event probabilities with a confidence interval using parametric approach
  dat.calib.pv.4 <- calib_pv(data.mstate = msebmtcal,
                                  data.raw = ebmtcal,
                                  j = 1,
                                  s = 0,
                                  t = 1826,
                                  tp.pred = tp.pred,
                                  curve.type = "rcs",
                                  CI = 95, CI.type = "parametric",
                                  data.pred.plot = NULL, transitions.out = c(1))

  expect_equal(ncol(dat.calib.pv.4[["plotdata"]][[1]]), 5)

  expect_equal(dat.calib.pv.1[["plotdata"]][[1]]$obs, dat.calib.pv.4[["plotdata"]][[1]]$obs)
  expect_equal(dat.calib.pv.1[["plotdata"]][[1]]$pred, dat.calib.pv.4[["plotdata"]][[1]]$pred)

  expect_no_error(summary(dat.calib.pv.4))

})


### Add some tests for when each of group.vars and n.pctls are specified
test_that("check calib_pv output, (j = 1, s = 0), groups.vars and n.pctls = NULL", {

  ## Reduce to 50 individuals
  # Extract the predicted transition probabilities out of state j = 1 for first 100 individuals
  tp.pred <- tps0 |>
    dplyr::filter(id %in% 1:50) |>
    dplyr::filter(j == 1) |>
    dplyr::select(any_of(paste("pstate", 1:6, sep = "")))
  # Reduce ebmtcal to first 50 individuals
  ebmtcal <- ebmtcal |> dplyr::filter(id %in% 1:50)
  # Reduce msebmtcal.cmprsk to first 100 individuals
  msebmtcal <- msebmtcal |> dplyr::filter(id %in% 1:50)

  ## Calculate observed event probabilities using transitions.out = NULL
  dat.calib.pv.1 <- calib_pv(data.mstate = msebmtcal,
                                  data.raw = ebmtcal,
                                  j = 1,
                                  s = 0,
                                  t = 1826,
                                  tp.pred = tp.pred,
                                  curve.type = "loess",
                                  loess.span = 1,
                                  loess.degree = 1,
                                  group.vars = c("year"),
                                  n.pctls = 2,
                                  data.pred.plot = NULL, transitions.out = NULL)

  expect_equal(ncol(dat.calib.pv.1[["plotdata"]][[1]]), 3)
  expect_equal(length(dat.calib.pv.1[["plotdata"]]), 6)

  ## Calculate observed event probabilities using transitions.out = NULL
  dat.calib.pv.2 <- calib_pv(data.mstate = msebmtcal,
                                  data.raw = ebmtcal,
                                  j = 1,
                                  s = 0,
                                  t = 1826,
                                  tp.pred = tp.pred,
                                  curve.type = "loess",
                                  loess.span = 1,
                                  loess.degree = 1,
                                  group.vars = c("year"),
                                  data.pred.plot = NULL, transitions.out = NULL)

  expect_equal(ncol(dat.calib.pv.2[["plotdata"]][[1]]), 3)
  expect_equal(length(dat.calib.pv.2[["plotdata"]]), 6)

  ## Calculate observed event probabilities using transitions.out = NULL
  dat.calib.pv.3 <- calib_pv(data.mstate = msebmtcal,
                                  data.raw = ebmtcal,
                                  j = 1,
                                  s = 0,
                                  t = 1826,
                                  tp.pred = tp.pred,
                                  curve.type = "loess",
                                  loess.span = 1,
                                  loess.degree = 1,
                                  n.pctls = 2,
                                  data.pred.plot = NULL, transitions.out = NULL)

  expect_equal(ncol(dat.calib.pv.3[["plotdata"]][[1]]), 3)
  expect_equal(length(dat.calib.pv.3[["plotdata"]]), 6)

})


### Add some tests where we expect errors, if requesting things that aren't possible
test_that("check calib_pv output, (j = 1, s = 0), cause errors", {

  ## Reduce to 50 individuals
  # Extract the predicted transition probabilities out of state j = 1 for first 100 individuals
  tp.pred <- tps0 |>
    dplyr::filter(id %in% 1:50) |>
    dplyr::filter(j == 1) |>
    dplyr::select(any_of(paste("pstate", 1:6, sep = "")))
  # Reduce ebmtcal to first 50 individuals
  ebmtcal <- ebmtcal |> dplyr::filter(id %in% 1:50)
  # Reduce msebmtcal.cmprsk to first 100 individuals
  msebmtcal <- msebmtcal |> dplyr::filter(id %in% 1:50)

  ## Request bootstrap confidence interval and don't give number of bootstrap replicates (for either rcs or parametric)
  expect_error(calib_pv(data.mstate = msebmtcal,
                             data.raw = ebmtcal,
                             j = 1,
                             s = 0,
                             t = 1826,
                             tp.pred = tp.pred,
                             curve.type = "loess",
                             CI = 95,
                             CI.type = "bootstrap",
                             data.pred.plot = NULL, transitions.out = NULL))

  expect_error(calib_pv(data.mstate = msebmtcal,
                             data.raw = ebmtcal,
                             j = 1,
                             s = 0,
                             t = 1826,
                             tp.pred = tp.pred,
                             curve.type = "rcs",
                             CI = 95,
                             CI.type = "bootstrap",
                             data.pred.plot = NULL, transitions.out = NULL))

})
