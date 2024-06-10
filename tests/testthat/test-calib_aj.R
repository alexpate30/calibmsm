###
### Tests for calibration curves produced using pseudo-values (calib.type = 'AJ')
###

### Run tests for pv.n.pctls = NULL and pv.group.vars = NULL
test_that("check calib_aj, pv.n.pctls = NULL and pv.group.vars = NULL", {

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
  dat.calib.aj.1 <- calib_msm(data.ms = msebmtcal,
                              data.raw = ebmtcal,
                              j = 1,
                              s = 0,
                              t = 1826,
                              tp.pred = tp.pred,
                              calib.type = 'aj',
                              tp.pred.plot = NULL, transitions.out = NULL)

  expect_equal(class(dat.calib.aj.1), c("calib_aj", "calib_msm"))
  expect_equal(length(dat.calib.aj.1[["mean"]]), 6)

  ## Calculate observed event probabilities using transitions.out = NULL
  dat.calib.aj.CI.1 <- calib_msm(data.ms = msebmtcal,
                              data.raw = ebmtcal,
                              j = 1,
                              s = 0,
                              t = 1826,
                              tp.pred = tp.pred,
                              calib.type = 'aj',
                              CI = 95,
                              CI.R.boot = 10,
                              tp.pred.plot = NULL, transitions.out = NULL)

  expect_equal(class(dat.calib.aj.CI.1), c("calib_aj", "calib_msm"))
  expect_equal(length(dat.calib.aj.CI.1[["mean"]]), 6)
  expect_equal(length(dat.calib.aj.CI.1[["mean"]][[1]]), 3)
  expect_equal(as.numeric(dat.calib.aj.1[["mean"]][1]), as.numeric(dat.calib.aj.CI.1[["mean"]][[1]][1]))
  expect_equal(as.numeric(dat.calib.aj.1[["mean"]][6]), as.numeric(dat.calib.aj.CI.1[["mean"]][[6]][1]))

})


### Run tets pv.n.pctls specified
test_that("check calib_pv output, pv.n.pctls specified", {

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
  dat.calib.aj.1 <- calib_msm(data.ms = msebmtcal,
                              data.raw = ebmtcal,
                              j = 1,
                              s = 0,
                              t = 1826,
                              tp.pred = tp.pred,
                              calib.type = 'aj',
                              pv.n.pctls = 2,
                              tp.pred.plot = NULL, transitions.out = NULL)

  expect_equal(class(dat.calib.aj.1), c("calib_aj", "calib_msm"))
  expect_equal(length(dat.calib.aj.1[["mean"]]), 6)

  ## Calculate observed event probabilities using transitions.out = NULL
  dat.calib.aj.CI.1 <- calib_msm(data.ms = msebmtcal,
                                 data.raw = ebmtcal,
                                 j = 1,
                                 s = 0,
                                 t = 1826,
                                 tp.pred = tp.pred,
                                 calib.type = 'aj',
                                 pv.n.pctls = 2,
                                 CI = 95,
                                 CI.R.boot = 10,
                                 tp.pred.plot = NULL, transitions.out = NULL)

  expect_equal(class(dat.calib.aj.CI.1), c("calib_aj", "calib_msm"))
  expect_equal(length(dat.calib.aj.CI.1[["mean"]]), 6)
  expect_equal(length(dat.calib.aj.CI.1[["mean"]][[1]]), 3)
  expect_equal(as.numeric(dat.calib.aj.1[["mean"]][1]), as.numeric(dat.calib.aj.CI.1[["mean"]][[1]][1]))
  expect_equal(as.numeric(dat.calib.aj.1[["mean"]][6]), as.numeric(dat.calib.aj.CI.1[["mean"]][[6]][1]))

})


### Run tests pv.group.vars specified
test_that("check calib_pv output, pv.group.vars specified", {

  skip_on_cran()

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
  dat.calib.aj.1 <- calib_msm(data.ms = msebmtcal,
                              data.raw = ebmtcal,
                              j = 1,
                              s = 0,
                              t = 1826,
                              tp.pred = tp.pred,
                              calib.type = 'aj',
                              pv.group.vars = c("year"),
                              tp.pred.plot = NULL, transitions.out = NULL)

  expect_equal(class(dat.calib.aj.1), c("calib_aj", "calib_msm"))
  expect_equal(length(dat.calib.aj.1[["mean"]]), 6)

  ## Calculate observed event probabilities using transitions.out = NULL
  dat.calib.aj.CI.1 <- calib_msm(data.ms = msebmtcal,
                                 data.raw = ebmtcal,
                                 j = 1,
                                 s = 0,
                                 t = 1826,
                                 tp.pred = tp.pred,
                                 calib.type = 'aj',
                                 pv.group.vars = c("year"),
                                 CI = 95,
                                 CI.R.boot = 10,
                                 tp.pred.plot = NULL, transitions.out = NULL)

  expect_equal(class(dat.calib.aj.CI.1), c("calib_aj", "calib_msm"))
  expect_equal(length(dat.calib.aj.CI.1[["mean"]]), 6)
  expect_equal(length(dat.calib.aj.CI.1[["mean"]][[1]]), 3)
  expect_equal(as.numeric(dat.calib.aj.1[["mean"]][1]), as.numeric(dat.calib.aj.CI.1[["mean"]][[1]][1]))
  expect_equal(as.numeric(dat.calib.aj.1[["mean"]][6]), as.numeric(dat.calib.aj.CI.1[["mean"]][[6]][1]))

})


### Run tests pv.group.vars and pv.n.pctls specified
test_that("check calib_pv output, pv.group.vars and pv.n.pctls specified", {

  skip_on_cran()

  ## Reduce to 50 individuals
  # Extract the predicted transition probabilities out of state j = 1 for first 50 individuals
  tp.pred <- tps0 |>
    dplyr::filter(id %in% 1:100) |>
    dplyr::filter(j == 1) |>
    dplyr::select(any_of(paste("pstate", 1:6, sep = "")))
  # Reduce ebmtcal to first 50 individuals
  ebmtcal <- ebmtcal |> dplyr::filter(id %in% 1:100)
  # Reduce msebmtcal.cmprsk to first 100 individuals
  msebmtcal <- msebmtcal |> dplyr::filter(id %in% 1:100)

  ## Calculate observed event probabilities using transitions.out = NULL
  dat.calib.aj.1 <- calib_msm(data.ms = msebmtcal,
                              data.raw = ebmtcal,
                              j = 1,
                              s = 0,
                              t = 1826,
                              tp.pred = tp.pred,
                              calib.type = 'aj',
                              pv.n.pctls = 2,
                              pv.group.vars = c("year"),
                              tp.pred.plot = NULL, transitions.out = NULL)

  expect_equal(class(dat.calib.aj.1), c("calib_aj", "calib_msm"))
  expect_equal(length(dat.calib.aj.1[["mean"]]), 6)

  ## Calculate observed event probabilities using transitions.out = NULL
  dat.calib.aj.CI.1 <- calib_msm(data.ms = msebmtcal,
                                 data.raw = ebmtcal,
                                 j = 1,
                                 s = 0,
                                 t = 1826,
                                 tp.pred = tp.pred,
                                 calib.type = 'aj',
                                 pv.n.pctls = 2,
                                 pv.group.vars = c("year"),
                                 CI = 95,
                                 CI.R.boot = 10,
                                 tp.pred.plot = NULL, transitions.out = NULL)

  expect_equal(class(dat.calib.aj.CI.1), c("calib_aj", "calib_msm"))
  expect_equal(length(dat.calib.aj.CI.1[["mean"]]), 6)
  expect_equal(length(dat.calib.aj.CI.1[["mean"]][[1]]), 3)
  expect_equal(as.numeric(dat.calib.aj.1[["mean"]][1]), as.numeric(dat.calib.aj.CI.1[["mean"]][[1]][1]))
  expect_equal(as.numeric(dat.calib.aj.1[["mean"]][6]), as.numeric(dat.calib.aj.CI.1[["mean"]][[6]][1]))

})
