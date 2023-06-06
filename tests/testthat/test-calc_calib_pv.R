### Run tests for when curve.type = "loess".
test_that("check calc_calib_pv output, (j = 3, s = 100), curve.type = loess", {

  ## Extract relevant predicted risks from tps0
  tp.pred <- dplyr::select(dplyr::filter(tps100, j == 3), dplyr::any_of(paste("pstate", 1:6, sep = "")))

  ## Calculate observed event probabilities using transitions.out = NULL
  dat.calib.pv.1 <- calc_calib_pv(data.mstate = msebmtcal,
                                  data.raw = ebmtcal,
                                  j = 3,
                                  s = 100,
                                  t.eval = 1826,
                                  tp.pred = tp.pred,
                                  curve.type = "loess",
                                  group.vars = c("year"),
                                  n.pctls = 2,
                                  data.pred.plot = NULL, transitions.out = NULL)

  # ## Calculate observed event probabilities using transitions.out = c(3,4,5,6)
  # dat.calib.pv.2 <- calc_calib_pv(data.mstate = msebmtcal,
  #                                  data.raw = ebmtcal,
  #                                  j = 3,
  #                                  s = 100,
  #                                  t.eval = 1826,
  #                                  tp.pred = tp.pred,
  #                                  curve.type = "loess",
  #                                  group.vars = c("year"),
  #                                  n.pctls = 2,
  #                                  data.pred.plot = NULL, transitions.out = c(3,4,5,6))

  expect_equal(dat.calib.pv.1[["metadata"]][["curve.type"]], "loess")
  expect_equal(ncol(dat.calib.pv.1[["plotdata"]][[1]]), 3)
  # expect_equal(ncol(dat.calib.pv.2[["plotdata"]][[1]]), 3)
  # expect_equal(dat.calib.pv.1[["plotdata"]][[1]], dat.calib.pv.2[["plotdata"]][[1]])

  # ## Calculate observed event probabilities using transitions.out = 3
  # dat.calib.pv.3 <- calc_calib_pv(data.mstate = msebmtcal,
  #                                 data.raw = ebmtcal,
  #                                 j = 3,
  #                                 s = 100,
  #                                 t.eval = 1826,
  #                                 tp.pred = tp.pred,
  #                                 curve.type = "loess",
  #                                 group.vars = c("year"),
  #                                 n.pctls = 2,
  #                                 data.pred.plot = NULL, transitions.out = 3)
  #
  # expect_equal(length(dat.calib.pv.3[["plotdata"]]), 1)
  # expect_equal(dat.calib.pv.1[["plotdata"]][[1]], dat.calib.pv.3[["plotdata"]][[1]])

  ## Calculate observed event probabilities with a confidence interval using bootstrapping and transitions.out = NULL
  dat.calib.pv.4 <- calc_calib_pv(data.mstate = msebmtcal,
                                  data.raw = ebmtcal,
                                  j = 3,
                                  s = 100,
                                  t.eval = 1826,
                                  tp.pred = tp.pred,
                                  curve.type = "loess",
                                  group.vars = c("year"),
                                  n.pctls = 2,
                                  CI = 95, CI.type = "bootstrap", CI.R.boot = 3,
                                  data.pred.plot = NULL, transitions.out = c(3,4,5,6))

  expect_equal(ncol(dat.calib.pv.4[["plotdata"]][[1]]), 5)

  expect_equal(dat.calib.pv.1[["plotdata"]][[1]]$obs, dat.calib.pv.4[["plotdata"]][[1]]$obs)
  expect_equal(dat.calib.pv.1[["plotdata"]][[1]]$pred, dat.calib.pv.4[["plotdata"]][[1]]$pred)

  expect_equal(dat.calib.pv.1[["plotdata"]][[4]]$obs, dat.calib.pv.4[["plotdata"]][[4]]$obs)
  expect_equal(dat.calib.pv.1[["plotdata"]][[4]]$pred, dat.calib.pv.4[["plotdata"]][[4]]$pred)

  # ## Calculate observed event probabilities with a confidence interval using bootstrapping and transitions.out = c(3,4,5,6)
  # dat.calib.pv.5 <- calc_calib_pv(data.mstate = msebmtcal,
  #                                 data.raw = ebmtcal,
  #                                 j = 3,
  #                                 s = 100,
  #                                 t.eval = 1826,
  #                                 tp.pred = tp.pred,
  #                                 curve.type = "loess",
  #                                 group.vars = c("year"),
  #                                 n.pctls = 2,
  #                                 CI = 95, CI.type = "bootstrap", CI.R.boot = 3,
  #                                 data.pred.plot = NULL, transitions.out = c(3,4,5,6))
  #
  # expect_equal(dat.calib.pv.4[["plotdata"]][[1]]$obs, dat.calib.pv.5[["plotdata"]][[1]]$obs)
  # expect_equal(dat.calib.pv.4[["plotdata"]][[1]]$pred, dat.calib.pv.5[["plotdata"]][[1]]$pred)
  #
  # expect_equal(dat.calib.pv.4[["plotdata"]][[4]]$obs, dat.calib.pv.5[["plotdata"]][[4]]$obs)
  # expect_equal(dat.calib.pv.4[["plotdata"]][[4]]$pred, dat.calib.pv.5[["plotdata"]][[4]]$pred)

  # ## Calculate observed event probabilities with a confidence interval using bootstrapping and transitions.out = 3
  # dat.calib.pv.6 <- calc_calib_pv(data.mstate = msebmtcal,
  #                                 data.raw = ebmtcal,
  #                                 j = 3,
  #                                 s = 100,
  #                                 t.eval = 1826,
  #                                 tp.pred = tp.pred,
  #                                 curve.type = "loess",
  #                                 group.vars = c("year"),
  #                                 n.pctls = 2,
  #                                 CI = 95, CI.type = "bootstrap", CI.R.boot = 3,
  #                                 data.pred.plot = NULL, transitions.out = 3)
  #
  # expect_equal(length(dat.calib.pv.6[["plotdata"]]), 1)
  # expect_equal(ncol(dat.calib.pv.6[["plotdata"]][[1]]), 5)
  # expect_equal(dat.calib.pv.4[["plotdata"]][[1]]$obs, dat.calib.pv.6[["plotdata"]][[1]]$obs)
  # expect_equal(dat.calib.pv.4[["plotdata"]][[1]]$pred, dat.calib.pv.6[["plotdata"]][[1]]$pred)

  # ## Calculate observed event probabilities with a confidence interval using bootstrapping and transitions.out = 6
  # dat.calib.pv.7 <- calc_calib_pv(data.mstate = msebmtcal,
  #                                 data.raw = ebmtcal,
  #                                 j = 3,
  #                                 s = 100,
  #                                 t.eval = 1826,
  #                                 tp.pred = tp.pred,
  #                                 curve.type = "loess",
  #                                 group.vars = c("year"),
  #                                 n.pctls = 2,
  #                                 CI = 95, CI.type = "bootstrap", CI.R.boot = 10,
  #                                 data.pred.plot = NULL, transitions.out = 6)
  #
  # expect_equal(length(dat.calib.pv.7[["plotdata"]]), 1)
  # expect_equal(ncol(dat.calib.pv.7[["plotdata"]][[1]]), 5)
  # expect_false(isTRUE(all.equal(dat.calib.pv.4[["plotdata"]][[1]]$obs, dat.calib.pv.7[["plotdata"]][[1]]$obs)))
  #
  # ## Calculate observed event probabilities with a confidence interval using bootstrapping, transitions.out = 6 and a smaller confidence interval
  # dat.calib.pv.8 <- calc_calib_pv(data.mstate = msebmtcal,
  #                                 data.raw = ebmtcal,
  #                                 j = 3,
  #                                 s = 100,
  #                                 t.eval = 1826,
  #                                 tp.pred = tp.pred,
  #                                 curve.type = "loess",
  #                                 group.vars = c("year"),
  #                                 n.pctls = 2,
  #                                 CI = 50, CI.type = "bootstrap", CI.R.boot = 10,
  #                                 data.pred.plot = NULL, transitions.out = 6)
  #
  # expect_equal(ncol(dat.calib.pv.8[["plotdata"]][[1]]), 5)
  # expect_true(mean(dat.calib.pv.7[["plotdata"]][[1]]$obs.upper - dat.calib.pv.7[["plotdata"]][[1]]$obs.lower, na.rm = TRUE) >
  #   mean(dat.calib.pv.8[["plotdata"]][[1]]$obs.upper - dat.calib.pv.8[["plotdata"]][[1]]$obs.lower, na.rm = TRUE))

  ## Calculate observed event probabilities with a confidence interval using bootstrapping, transitions.out = NULL and defining data.pred.plot manually

  ### Create landmark ids and extract data.pred.plot correct
  id.lmk <- dat.calib.pv.1[["plotdata"]][["state3"]]$id
  data.pred.plot.in <- tps100 |>
    dplyr::filter(id %in% id.lmk) |>
    dplyr::filter(j == 3) |>
    dplyr::select(any_of(paste("pstate", 1:6, sep = "")))

  ## No confidence interval
  dat.calib.pv.9 <- calc_calib_pv(data.mstate = msebmtcal,
                                  data.raw = ebmtcal,
                                  j = 3,
                                  s = 100,
                                  t.eval = 1826,
                                  tp.pred = tp.pred,
                                  curve.type = "loess",
                                  group.vars = c("year"),
                                  n.pctls = 2,
                                  data.pred.plot = data.pred.plot.in, transitions.out = NULL)

  ## Should be one less column in plotdata (no patient ids)
  expect_equal(ncol(dat.calib.pv.9[["plotdata"]][[1]]), 2)
  expect_equal(dat.calib.pv.1[["plotdata"]][[1]]$obs, dat.calib.pv.9[["plotdata"]][[1]]$obs)
  expect_equal(dat.calib.pv.1[["plotdata"]][[1]]$pred, dat.calib.pv.9[["plotdata"]][[1]]$pred)

  ## With confidence interval
  dat.calib.pv.10 <- calc_calib_pv(data.mstate = msebmtcal,
                                  data.raw = ebmtcal,
                                  j = 3,
                                  s = 100,
                                  t.eval = 1826,
                                  tp.pred = tp.pred,
                                  curve.type = "loess",
                                  group.vars = c("year"),
                                  n.pctls = 2,
                                  CI = 95, CI.type = "bootstrap", CI.R.boot = 3,
                                  data.pred.plot = data.pred.plot.in, transitions.out = NULL)

  expect_equal(ncol(dat.calib.pv.10[["plotdata"]][[1]]), 4)
  expect_equal(dat.calib.pv.1[["plotdata"]][[1]]$obs, dat.calib.pv.10[["plotdata"]][[1]]$obs)
  expect_equal(dat.calib.pv.1[["plotdata"]][[1]]$pred, dat.calib.pv.10[["plotdata"]][[1]]$pred)

})

### Run some of these tests for when curve.type = "rcs".
### Note for this, the CI.type can either be 'bootstrap' or 'parametric'

### Run tests for when curve.type = "rcs" and CI.type = "bootstrap".
test_that("check calc_calib_pv output, (j = 3, s = 100), curve.type = rcs, CI.type = bootstrap.", {

  ## Extract relevant predicted risks from tps0
  tp.pred <- dplyr::select(dplyr::filter(tps100, j == 3), dplyr::any_of(paste("pstate", 1:6, sep = "")))

  ## Calculate observed event probabilities using transitions.out = NULL
  dat.calib.pv.1 <- calc_calib_pv(data.mstate = msebmtcal,
                                  data.raw = ebmtcal,
                                  j = 3,
                                  s = 100,
                                  t.eval = 1826,
                                  tp.pred = tp.pred,
                                  curve.type = "rcs",
                                  group.vars = c("year"),
                                  n.pctls = 2,
                                  data.pred.plot = NULL, transitions.out = NULL)


  expect_equal(dat.calib.pv.1[["metadata"]][["curve.type"]], "rcs")
  expect_equal(ncol(dat.calib.pv.1[["plotdata"]][[1]]), 3)

  ## Calculate observed event probabilities with a confidence interval using bootstrapping
  dat.calib.pv.4 <- calc_calib_pv(data.mstate = msebmtcal,
                                  data.raw = ebmtcal,
                                  j = 3,
                                  s = 100,
                                  t.eval = 1826,
                                  tp.pred = tp.pred,
                                  curve.type = "rcs",
                                  group.vars = c("year"),
                                  n.pctls = 2,
                                  CI = 95, CI.type = "bootstrap", CI.R.boot = 3,
                                  data.pred.plot = NULL, transitions.out = NULL)

  expect_equal(ncol(dat.calib.pv.4[["plotdata"]][[1]]), 5)

  expect_equal(dat.calib.pv.1[["plotdata"]][[1]]$obs, dat.calib.pv.4[["plotdata"]][[1]]$obs)
  expect_equal(dat.calib.pv.1[["plotdata"]][[1]]$pred, dat.calib.pv.4[["plotdata"]][[1]]$pred)

  expect_equal(dat.calib.pv.1[["plotdata"]][[4]]$obs, dat.calib.pv.4[["plotdata"]][[4]]$obs)
  expect_equal(dat.calib.pv.1[["plotdata"]][[4]]$pred, dat.calib.pv.4[["plotdata"]][[4]]$pred)

  ## Calculate observed event probabilities with a confidence interval using parametric approach
  dat.calib.pv.4 <- calc_calib_pv(data.mstate = msebmtcal,
                                  data.raw = ebmtcal,
                                  j = 3,
                                  s = 100,
                                  t.eval = 1826,
                                  tp.pred = tp.pred,
                                  curve.type = "rcs",
                                  group.vars = c("year"),
                                  n.pctls = 2,
                                  CI = 95, CI.type = "parametric",
                                  data.pred.plot = NULL, transitions.out = NULL)

  expect_equal(ncol(dat.calib.pv.4[["plotdata"]][[1]]), 5)

  expect_equal(dat.calib.pv.1[["plotdata"]][[1]]$obs, dat.calib.pv.4[["plotdata"]][[1]]$obs)
  expect_equal(dat.calib.pv.1[["plotdata"]][[1]]$pred, dat.calib.pv.4[["plotdata"]][[1]]$pred)

  expect_equal(dat.calib.pv.1[["plotdata"]][[4]]$obs, dat.calib.pv.4[["plotdata"]][[4]]$obs)
  expect_equal(dat.calib.pv.1[["plotdata"]][[4]]$pred, dat.calib.pv.4[["plotdata"]][[4]]$pred)

})


### Add some tests for when each of group.vars and n.pctls are left as NULL
test_that("check calc_calib_pv output, (j = 3, s = 100), groups.vars and n.pctls = NULL", {

  ## Extract relevant predicted risks from tps0
  tp.pred <- dplyr::select(dplyr::filter(tps100, j == 3), dplyr::any_of(paste("pstate", 1:6, sep = "")))

  ## Calculate observed event probabilities using transitions.out = NULL
  dat.calib.pv.1 <- calc_calib_pv(data.mstate = msebmtcal,
                                  data.raw = ebmtcal,
                                  j = 3,
                                  s = 100,
                                  t.eval = 1826,
                                  tp.pred = tp.pred,
                                  curve.type = "rcs",
                                  data.pred.plot = NULL, transitions.out = NULL)

  expect_equal(ncol(dat.calib.pv.1[["plotdata"]][[1]]), 3)
  expect_equal(length(dat.calib.pv.1[["plotdata"]]), 4)

  ## Calculate observed event probabilities using transitions.out = NULL
  dat.calib.pv.2 <- calc_calib_pv(data.mstate = msebmtcal,
                                  data.raw = ebmtcal,
                                  j = 3,
                                  s = 100,
                                  t.eval = 1826,
                                  tp.pred = tp.pred,
                                  curve.type = "rcs",
                                  group.vars = c("year"),
                                  data.pred.plot = NULL, transitions.out = NULL)

  expect_equal(ncol(dat.calib.pv.2[["plotdata"]][[1]]), 3)
  expect_equal(length(dat.calib.pv.2[["plotdata"]]), 4)

  ## Calculate observed event probabilities using transitions.out = NULL
  dat.calib.pv.3 <- calc_calib_pv(data.mstate = msebmtcal,
                                  data.raw = ebmtcal,
                                  j = 3,
                                  s = 100,
                                  t.eval = 1826,
                                  tp.pred = tp.pred,
                                  curve.type = "rcs",
                                  n.pctls = 2,
                                  data.pred.plot = NULL, transitions.out = NULL)

  expect_equal(ncol(dat.calib.pv.3[["plotdata"]][[1]]), 3)
  expect_equal(length(dat.calib.pv.3[["plotdata"]]), 4)

})


### Add some tests where we expect errors, if requesting things that aren't possible
test_that("check calc_calib_pv output, (j = 3, s = 100), cause errors", {

  ## Extract relevant predicted risks from tps0
  tp.pred <- dplyr::select(dplyr::filter(tps100, j == 3), dplyr::any_of(paste("pstate", 1:6, sep = "")))

  ## Request confidence interval but don't state what type
  expect_error(calc_calib_pv(data.mstate = msebmtcal,
                             data.raw = ebmtcal,
                             j = 3,
                             s = 100,
                             t.eval = 1826,
                             tp.pred = tp.pred,
                             curve.type = "rcs",
                             CI = 95,
                             data.pred.plot = NULL, transitions.out = NULL))

  ## Request confidence interval for loess, and specify parametric
  expect_error(calc_calib_pv(data.mstate = msebmtcal,
                             data.raw = ebmtcal,
                             j = 3,
                             s = 100,
                             t.eval = 1826,
                             tp.pred = tp.pred,
                             curve.type = "loess",
                             CI = 95,
                             CI.type = "parametric",
                             data.pred.plot = NULL, transitions.out = NULL))

  ## Request bootstrap confidence interval and don't give number of bootstrap replicates (for either rcs or parametric)
  expect_error(calc_calib_pv(data.mstate = msebmtcal,
                             data.raw = ebmtcal,
                             j = 3,
                             s = 100,
                             t.eval = 1826,
                             tp.pred = tp.pred,
                             curve.type = "loess",
                             CI = 95,
                             CI.type = "bootstrap",
                             data.pred.plot = NULL, transitions.out = NULL))

  expect_error(calc_calib_pv(data.mstate = msebmtcal,
                             data.raw = ebmtcal,
                             j = 3,
                             s = 100,
                             t.eval = 1826,
                             tp.pred = tp.pred,
                             curve.type = "rcs",
                             CI = 95,
                             CI.type = "bootstrap",
                             data.pred.plot = NULL, transitions.out = NULL))

})
