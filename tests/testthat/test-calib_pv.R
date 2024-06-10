###
### Tests for calibration curves produced using pseudo-values (calib.type = 'pv')
###

### Run tests for when curve.type = "loess" and CI.type = "bootstrap".
test_that("check calib_pv output, (j = 1, s = 0), curve.type = loess, CI.type = bootstrap", {

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
  dat.calib.pv.1 <- calib_msm(data.ms = msebmtcal,
                              data.raw = ebmtcal,
                              j = 1,
                              s = 0,
                              t = 1826,
                              tp.pred = tp.pred,
                              calib.type = 'pv',
                              curve.type = "loess",
                              tp.pred.plot = NULL, transitions.out = NULL)

  expect_equal(class(dat.calib.pv.1), c("calib_pv", "calib_msm"))
  expect_equal(dat.calib.pv.1[["metadata"]][["curve.type"]], "loess")
  expect_equal(ncol(dat.calib.pv.1[["plotdata"]][[1]]), 4)
  expect_no_error(summary(dat.calib.pv.1))

  ## Check same results when just calculating pseudo-values for first three individuals
  dat.calib.pv.ids.1 <- calib_msm(data.ms = msebmtcal,
                                  data.raw = ebmtcal,
                                  j = 1,
                                  s = 0,
                                  t = 1826,
                                  tp.pred = tp.pred,
                                  calib.type = 'pv',
                                  pv.ids = 1:3,
                                  tp.pred.plot = NULL, transitions.out = NULL)

  expect_equal(dat.calib.pv.1[["plotdata"]][[1]][1:3, "pv"], dat.calib.pv.ids.1[[1]][,2])
  expect_equal(dat.calib.pv.1[["plotdata"]][[6]][1:3, "pv"], dat.calib.pv.ids.1[[1]][,7])

  ## Calculate observed event probabilities with a confidence interval using bootstrapping and transitions.out = NULL
  expect_warning(calib_msm(data.ms = msebmtcal,
                           data.raw = ebmtcal,
                           j = 1,
                           s = 0,
                           t = 1826,
                           tp.pred = tp.pred,
                           calib.type = 'pv',
                           curve.type = "loess",
                           CI = 95, CI.type = "bootstrap", CI.R.boot = 3,
                           tp.pred.plot = NULL, transitions.out = c(1)))

  dat.calib.pv.4 <- suppressWarnings(calib_msm(data.ms = msebmtcal,
                                               data.raw = ebmtcal,
                                               j = 1,
                                               s = 0,
                                               t = 1826,
                                               tp.pred = tp.pred,
                                               calib.type = 'pv',
                                               curve.type = "loess",
                                               CI = 95, CI.type = "bootstrap", CI.R.boot = 3,
                                               tp.pred.plot = NULL, transitions.out = c(1,2)))

  expect_equal(class(dat.calib.pv.4), c("calib_pv", "calib_msm"))
  expect_equal(ncol(dat.calib.pv.4[["plotdata"]][[1]]), 5)

  expect_equal(dat.calib.pv.1[["plotdata"]][[1]]$obs, dat.calib.pv.4[["plotdata"]][[1]]$obs)
  expect_equal(dat.calib.pv.1[["plotdata"]][[1]]$pred, dat.calib.pv.4[["plotdata"]][[1]]$pred)

  expect_equal(dat.calib.pv.1[["plotdata"]][[2]]$obs, dat.calib.pv.4[["plotdata"]][[2]]$obs)
  expect_equal(dat.calib.pv.1[["plotdata"]][[2]]$pred, dat.calib.pv.4[["plotdata"]][[2]]$pred)

  expect_no_error(summary(dat.calib.pv.4))

  ## Calculate observed event probabilities with a confidence interval using bootstrapping, transitions.out = NULL and defining tp.pred.plot manually

  ### Create landmark ids and extract tp.pred.plot correct
  id.lmk <- 1:50
  tp.pred.plot <- tps0 |>
    dplyr::filter(id %in% id.lmk) |>
    dplyr::filter(j == 1) |>
    dplyr::select(any_of(paste("pstate", 1:6, sep = "")))

  ## No confidence interval
  dat.calib.pv.9 <- suppressWarnings(calib_msm(data.ms = msebmtcal,
                                               data.raw = ebmtcal,
                                               j = 1,
                                               s = 0,
                                               t = 1826,
                                               tp.pred = tp.pred,
                                               calib.type = 'pv',
                                               curve.type = "loess",
                                               tp.pred.plot = tp.pred.plot, transitions.out = NULL))

  ## Should be one less column in plotdata (no patient ids)
  expect_equal(class(dat.calib.pv.9), c("calib_pv", "calib_msm"))
  expect_equal(ncol(dat.calib.pv.9[["plotdata"]][[1]]), 3)
  expect_equal(nrow(dat.calib.pv.9[["plotdata"]][[1]]), 50)
  expect_equal(dat.calib.pv.1[["plotdata"]][[1]]$obs, dat.calib.pv.9[["plotdata"]][[1]]$obs)
  expect_equal(dat.calib.pv.1[["plotdata"]][[1]]$pred, dat.calib.pv.9[["plotdata"]][[1]]$pred)

  expect_no_error(summary(dat.calib.pv.9))

  ## With confidence interval
  dat.calib.pv.10 <- suppressWarnings(calib_msm(data.ms = msebmtcal,
                                                data.raw = ebmtcal,
                                                j = 1,
                                                s = 0,
                                                t = 1826,
                                                tp.pred = tp.pred,
                                                calib.type = 'pv',
                                                curve.type = "loess",
                                                CI = 95, CI.type = "bootstrap", CI.R.boot = 3,
                                                tp.pred.plot = tp.pred.plot, transitions.out = NULL))

  expect_equal(class(dat.calib.pv.10), c("calib_pv", "calib_msm"))
  expect_equal(ncol(dat.calib.pv.10[["plotdata"]][[1]]), 4)
  expect_equal(nrow(dat.calib.pv.10[["plotdata"]][[1]]), 50)
  expect_equal(dat.calib.pv.1[["plotdata"]][[1]]$obs, dat.calib.pv.10[["plotdata"]][[1]]$obs)
  expect_equal(dat.calib.pv.1[["plotdata"]][[1]]$pred, dat.calib.pv.10[["plotdata"]][[1]]$pred)

  expect_no_error(summary(dat.calib.pv.10))

})

### Run tests for when curve.type = "loess" and CI.type = "bootstrap".
test_that("check calib_pv output, (j = 1, s = 0), curve.type = loess, CI.type = parametric", {

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
  dat.calib.pv.1 <- calib_msm(data.ms = msebmtcal,
                              data.raw = ebmtcal,
                              j = 1,
                              s = 0,
                              t = 1826,
                              tp.pred = tp.pred,
                              calib.type = 'pv',
                              curve.type = "loess",
                              tp.pred.plot = NULL, transitions.out = NULL)

  expect_equal(dat.calib.pv.1[["metadata"]][["curve.type"]], "loess")
  expect_equal(ncol(dat.calib.pv.1[["plotdata"]][[1]]), 4)
  expect_no_error(summary(dat.calib.pv.1))

  ## Calculate observed event probabilities with a confidence interval using parametric approach
  dat.calib.pv.5 <- calib_msm(data.ms = msebmtcal,
                              data.raw = ebmtcal,
                              j = 1,
                              s = 0,
                              t = 1826,
                              tp.pred = tp.pred,
                              calib.type = 'pv',
                              curve.type = "loess",
                              CI = 95, CI.type = "parametric",
                              tp.pred.plot = NULL, transitions.out = c(1,2))

  expect_equal(ncol(dat.calib.pv.5[["plotdata"]][[1]]), 6)

  expect_equal(dat.calib.pv.1[["plotdata"]][[1]]$obs, dat.calib.pv.5[["plotdata"]][[1]]$obs)
  expect_equal(dat.calib.pv.1[["plotdata"]][[1]]$pred, dat.calib.pv.5[["plotdata"]][[1]]$pred)

  expect_equal(dat.calib.pv.1[["plotdata"]][[2]]$obs, dat.calib.pv.5[["plotdata"]][[2]]$obs)
  expect_equal(dat.calib.pv.1[["plotdata"]][[2]]$pred, dat.calib.pv.5[["plotdata"]][[2]]$pred)

  expect_no_error(summary(dat.calib.pv.5))

  ## Calculate observed event probabilities with a confidence interval using bootstrapping, transitions.out = NULL and defining tp.pred.plot manually

  ### Create landmark ids and extract tp.pred.plot correct
  id.lmk <- 1:50
  tp.pred.plot <- tps0 |>
    dplyr::filter(id %in% id.lmk) |>
    dplyr::filter(j == 1) |>
    dplyr::select(any_of(paste("pstate", 1:6, sep = "")))

  ## With confidence interval
  dat.calib.pv.10 <- calib_msm(data.ms = msebmtcal,
                               data.raw = ebmtcal,
                               j = 1,
                               s = 0,
                               t = 1826,
                               tp.pred = tp.pred,
                               calib.type = 'pv',
                               curve.type = "loess",
                               CI = 95, CI.type = "parametric",
                               tp.pred.plot = tp.pred.plot, transitions.out = NULL)

  str(dat.calib.pv.10)
  expect_equal(ncol(dat.calib.pv.10[["plotdata"]][[1]]), 5)
  expect_equal(dat.calib.pv.1[["plotdata"]][[1]]$obs, dat.calib.pv.10[["plotdata"]][[1]]$obs)
  expect_equal(dat.calib.pv.1[["plotdata"]][[1]]$pred, dat.calib.pv.10[["plotdata"]][[1]]$pred)

  expect_no_error(summary(dat.calib.pv.10))

})


### Run tests for when curve.type = "rcs" and CI.type = "bootstrap" (not rerunning all of them for curve.type = rcs)
test_that("check calib_pv output, (j = 1, s = 0), curve.type = rcs, CI.type = bootstrap.", {

  skip_on_cran()

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
  dat.calib.pv.1 <- calib_msm(data.ms = msebmtcal,
                              data.raw = ebmtcal,
                              j = 1,
                              s = 0,
                              t = 1826,
                              tp.pred = tp.pred,
                              calib.type = 'pv',
                              curve.type = "rcs",
                              tp.pred.plot = NULL, transitions.out = c(1))

  expect_equal(dat.calib.pv.1[["metadata"]][["curve.type"]], "rcs")
  expect_equal(ncol(dat.calib.pv.1[["plotdata"]][[1]]), 4)
  expect_no_error(summary(dat.calib.pv.1))

  ## Calculate observed event probabilities with a confidence interval using bootstrapping
  dat.calib.pv.4 <- suppressWarnings(calib_msm(data.ms = msebmtcal,
                                               data.raw = ebmtcal,
                                               j = 1,
                                               s = 0,
                                               t = 1826,
                                               tp.pred = tp.pred,
                                               calib.type = 'pv',
                                               curve.type = "rcs",
                                               CI = 95, CI.type = "bootstrap", CI.R.boot = 3,
                                               tp.pred.plot = NULL, transitions.out = c(1)))

  expect_equal(ncol(dat.calib.pv.4[["plotdata"]][[1]]), 5)

  expect_equal(dat.calib.pv.1[["plotdata"]][[1]]$obs, dat.calib.pv.4[["plotdata"]][[1]]$obs)
  expect_equal(dat.calib.pv.1[["plotdata"]][[1]]$pred, dat.calib.pv.4[["plotdata"]][[1]]$pred)

  expect_no_error(summary(dat.calib.pv.4))

})

### Run tests for when curve.type = "rcs" and CI.type = "parametric" (not rerunning all of them for curve.type = rcs)
test_that("check calib_pv output, (j = 1, s = 0), curve.type = rcs, CI.type = bootstrap.", {

  skip_on_cran()

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
  dat.calib.pv.1 <- calib_msm(data.ms = msebmtcal,
                              data.raw = ebmtcal,
                              j = 1,
                              s = 0,
                              t = 1826,
                              tp.pred = tp.pred,
                              calib.type = 'pv',
                              curve.type = "rcs",
                              tp.pred.plot = NULL, transitions.out = c(1))

  expect_equal(dat.calib.pv.1[["metadata"]][["curve.type"]], "rcs")
  expect_equal(ncol(dat.calib.pv.1[["plotdata"]][[1]]), 4)
  expect_no_error(summary(dat.calib.pv.1))

  ## Calculate observed event probabilities with a confidence interval using parametric approach
  dat.calib.pv.4 <- calib_msm(data.ms = msebmtcal,
                              data.raw = ebmtcal,
                              j = 1,
                              s = 0,
                              t = 1826,
                              tp.pred = tp.pred,
                              calib.type = 'pv',
                              curve.type = "rcs",
                              CI = 95, CI.type = "parametric",
                              tp.pred.plot = NULL, transitions.out = c(1))

  expect_equal(ncol(dat.calib.pv.4[["plotdata"]][[1]]), 6)

  expect_equal(dat.calib.pv.1[["plotdata"]][[1]]$obs, dat.calib.pv.4[["plotdata"]][[1]]$obs)
  expect_equal(dat.calib.pv.1[["plotdata"]][[1]]$pred, dat.calib.pv.4[["plotdata"]][[1]]$pred)

  expect_no_error(summary(dat.calib.pv.4))

})


### Add some tests for when each of group.vars and pv.n.pctls are specified
test_that("check calib_pv output, (j = 1, s = 0), groups.vars and pv.n.pctls specified", {

  skip_on_cran()

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

  ## Calculate observed event probabilities when both pv.group.vars and pv.n.pctls are specified
  dat.calib.pv.1 <- calib_msm(data.ms = msebmtcal,
                              data.raw = ebmtcal,
                              j = 1,
                              s = 0,
                              t = 1826,
                              tp.pred = tp.pred,
                              calib.type = 'pv',
                              curve.type = "loess",
                              loess.span = 1,
                              loess.degree = 1,
                              pv.group.vars = c("year"),
                              pv.n.pctls = 2,
                              tp.pred.plot = NULL, transitions.out = NULL)

  expect_equal(ncol(dat.calib.pv.1[["plotdata"]][[1]]), 4)
  expect_equal(length(dat.calib.pv.1[["plotdata"]]), 6)

  ## Check same results when just calculating pseudo-values for first three individuals
  dat.calib.pv.ids.1 <- calib_msm(data.ms = msebmtcal,
                                  data.raw = ebmtcal,
                                  j = 1,
                                  s = 0,
                                  t = 1826,
                                  tp.pred = tp.pred,
                                  calib.type = 'pv',
                                  pv.group.vars = c("year"),
                                  pv.n.pctls = 2,
                                  pv.ids = 1:3,
                                  tp.pred.plot = NULL, transitions.out = NULL)

  expect_equal(dat.calib.pv.1[["plotdata"]][[1]][1:3, "pv"], dat.calib.pv.ids.1[[1]][,2])
  expect_equal(dat.calib.pv.1[["plotdata"]][[6]][1:3, "pv"], dat.calib.pv.ids.1[[1]][,7])

  ## Check same results when just calculating pseudo-values for first three individuals, but specify transitions 1 and 6
  dat.calib.pv.ids.1.tout <- calib_msm(data.ms = msebmtcal,
                                       data.raw = ebmtcal,
                                       j = 1,
                                       s = 0,
                                       t = 1826,
                                       tp.pred = tp.pred,
                                       calib.type = 'pv',
                                       pv.group.vars = c("year"),
                                       pv.n.pctls = 2,
                                       pv.ids = 1:3,
                                       tp.pred.plot = NULL, transitions.out = c(1,6))

  expect_equal(dat.calib.pv.ids.1.tout[[1]][,2], dat.calib.pv.ids.1[[1]][,2])
  expect_equal(dat.calib.pv.ids.1.tout[[1]][,3], dat.calib.pv.ids.1[[1]][,7])
  expect_equal(ncol(dat.calib.pv.ids.1.tout[["plotdata"]]), 3)

  ## Calculate observed event probabilities for pv.group.vars
  dat.calib.pv.2 <- calib_msm(data.ms = msebmtcal,
                              data.raw = ebmtcal,
                              j = 1,
                              s = 0,
                              t = 1826,
                              tp.pred = tp.pred,
                              calib.type = 'pv',
                              curve.type = "loess",
                              loess.span = 1,
                              loess.degree = 1,
                              pv.group.vars = c("year"),
                              tp.pred.plot = NULL, transitions.out = NULL)

  expect_equal(ncol(dat.calib.pv.2[["plotdata"]][[1]]), 4)
  expect_equal(length(dat.calib.pv.2[["plotdata"]]), 6)

  ## Check same results when just calculating pseudo-values for first three individuals
  dat.calib.pv.ids.2 <- calib_msm(data.ms = msebmtcal,
                                  data.raw = ebmtcal,
                                  j = 1,
                                  s = 0,
                                  t = 1826,
                                  tp.pred = tp.pred,
                                  calib.type = 'pv',
                                  pv.group.vars = c("year"),
                                  pv.ids = 1:3,
                                  tp.pred.plot = NULL, transitions.out = NULL)

  expect_equal(dat.calib.pv.2[["plotdata"]][[1]][1:3, "pv"], dat.calib.pv.ids.2[[1]][,2])
  expect_equal(dat.calib.pv.2[["plotdata"]][[6]][1:3, "pv"], dat.calib.pv.ids.2[[1]][,7])

  ## No need to test for transitions.out when pv.n.pctls not specified, because there are no computational gains and
  ## pseudo-values are just calculated for all states anyway.

  ## Calculate observed event probabilities for pv.n.pctls
  dat.calib.pv.3 <- calib_msm(data.ms = msebmtcal,
                              data.raw = ebmtcal,
                              j = 1,
                              s = 0,
                              t = 1826,
                              tp.pred = tp.pred,
                              calib.type = 'pv',
                              curve.type = "loess",
                              loess.span = 1,
                              loess.degree = 1,
                              pv.n.pctls = 2,
                              tp.pred.plot = NULL, transitions.out = NULL)

  expect_equal(ncol(dat.calib.pv.3[["plotdata"]][[1]]), 4)
  expect_equal(length(dat.calib.pv.3[["plotdata"]]), 6)

  ## Check same results when just calculating pseudo-values for first three individuals
  dat.calib.pv.ids.3 <- calib_msm(data.ms = msebmtcal,
                                  data.raw = ebmtcal,
                                  j = 1,
                                  s = 0,
                                  t = 1826,
                                  tp.pred = tp.pred,
                                  calib.type = 'pv',
                                  pv.n.pctls = 2,
                                  pv.ids = 1:3,
                                  tp.pred.plot = NULL, transitions.out = NULL)

  expect_equal(dat.calib.pv.3[["plotdata"]][[1]][1:3, "pv"], dat.calib.pv.ids.3[[1]][,2])
  expect_equal(dat.calib.pv.3[["plotdata"]][[6]][1:3, "pv"], dat.calib.pv.ids.3[[1]][,7])

  ## Check same results when just calculating pseudo-values for first three individuals, but specify transitions 1 and 6
  dat.calib.pv.ids.3.tout <- calib_msm(data.ms = msebmtcal,
                                       data.raw = ebmtcal,
                                       j = 1,
                                       s = 0,
                                       t = 1826,
                                       tp.pred = tp.pred,
                                       calib.type = 'pv',
                                       pv.n.pctls = 2,
                                       pv.ids = 1:3,
                                       tp.pred.plot = NULL, transitions.out = c(1,6))

  expect_equal(dat.calib.pv.ids.3.tout[[1]][,2], dat.calib.pv.ids.3[[1]][,2])
  expect_equal(dat.calib.pv.ids.3.tout[[1]][,3], dat.calib.pv.ids.3[[1]][,7])
  expect_equal(ncol(dat.calib.pv.ids.3.tout[["plotdata"]]), 3)

})



### Add some tests where we expect errors, if requesting things that aren't possible
test_that("check calib_pv output, (j = 1, s = 0), cause errors", {

  skip_on_cran()

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
  expect_error(calib_msm(data.ms = msebmtcal,
                         data.raw = ebmtcal,
                         j = 1,
                         s = 0,
                         t = 1826,
                         tp.pred = tp.pred,
                         calib.type = 'pv',
                         curve.type = "loess",
                         CI = 95,
                         CI.type = "bootstrap",
                         tp.pred.plot = NULL, transitions.out = NULL))

  expect_error(calib_msm(data.ms = msebmtcal,
                         data.raw = ebmtcal,
                         j = 1,
                         s = 0,
                         t = 1826,
                         tp.pred = tp.pred,
                         calib.type = 'pv',
                         curve.type = "rcs",
                         CI = 95,
                         CI.type = "bootstrap",
                         tp.pred.plot = NULL, transitions.out = NULL))

})


test_that("check calib_pv output, (j = 3, s = 100), pv.group.vars defined", {

  skip_on_cran()

  ## Extract relevant predicted risks from tps100
  tp.pred <- dplyr::select(dplyr::filter(tps100, j == 3), any_of(paste("pstate", 1:6, sep = "")))

  ## Calculate observed event probabilities
  dat.calib.pv <-
    calib_msm(data.ms = msebmtcal,
              data.raw = ebmtcal,
              j=3,
              s=100,
              t = 1826,
              tp.pred = tp.pred, calib.type = 'pv',
              curve.type = "rcs",
              rcs.nk = 3,
              pv.group.vars = c("year"))

  expect_type(dat.calib.pv, "list")
  expect_equal(class(dat.calib.pv), c("calib_pv", "calib_msm"))
  expect_length(dat.calib.pv[["plotdata"]], 4)
  expect_length(dat.calib.pv[["plotdata"]][["state3"]]$id, 413)
  expect_length(dat.calib.pv[["plotdata"]][["state6"]]$id, 413)
  expect_error(dat.calib.pv[["plotdata"]][[6]])
  expect_false(dat.calib.pv[["metadata"]]$CI)


})


test_that("check calib_pv output, (j = 3, s = 100), pv.n.pctls defined", {

  skip_on_cran()

  ## Extract relevant predicted risks from tps100
  tp.pred <- dplyr::select(dplyr::filter(tps100, j == 3), any_of(paste("pstate", 1:6, sep = "")))

  ## Calculate observed event probabilities
  dat.calib.pv <-
    calib_msm(data.ms = msebmtcal,
              data.raw = ebmtcal,
              j=3,
              s=100,
              t = 1826,
              tp.pred = tp.pred, calib.type = 'pv',
              curve.type = "rcs",
              rcs.nk = 3,
              pv.n.pctls = 2)

  expect_type(dat.calib.pv, "list")
  expect_equal(class(dat.calib.pv), c("calib_pv", "calib_msm"))
  expect_length(dat.calib.pv[["plotdata"]], 4)
  expect_length(dat.calib.pv[["plotdata"]][["state3"]]$id, 413)
  expect_length(dat.calib.pv[["plotdata"]][["state6"]]$id, 413)
  expect_error(dat.calib.pv[["plotdata"]][[6]])
  expect_false(dat.calib.pv[["metadata"]]$CI)

})


test_that("check calib_pv output, (j = 3, s = 100), pv.group.vars and pv.n.pctls defined", {

  skip_on_cran()

  ## Extract relevant predicted risks from tps100
  tp.pred <- dplyr::select(dplyr::filter(tps100, j == 3), any_of(paste("pstate", 1:6, sep = "")))

  ## Calculate observed event probabilities
  dat.calib.pv <-
    calib_msm(data.ms = msebmtcal,
              data.raw = ebmtcal,
              j=3,
              s=100,
              t = 1826,
              tp.pred = tp.pred, calib.type = 'pv',
              curve.type = "rcs",
              rcs.nk = 3,
              pv.group.vars = c("year"),
              pv.n.pctls = 2)

  expect_type(dat.calib.pv, "list")
  expect_equal(class(dat.calib.pv), c("calib_pv", "calib_msm"))
  expect_length(dat.calib.pv[["plotdata"]], 4)
  expect_length(dat.calib.pv[["plotdata"]][["state3"]]$id, 413)
  expect_length(dat.calib.pv[["plotdata"]][["state6"]]$id, 413)
  expect_error(dat.calib.pv[["plotdata"]][[6]])
  expect_false(dat.calib.pv[["metadata"]]$CI)

})


test_that("check calib_pv output, (j = 1, s = 0), pv.precalc", {

  skip_on_cran()

  ## Extract relevant predicted risks from tps100
  tp.pred <- dplyr::select(dplyr::filter(tps0, j == 1), any_of(paste("pstate", 1:6, sep = "")))

  ## Define pv.precalc to be the estimated predicted probabilities
  pv.precalc <- tp.pred

  ## Calculate observed event probabilities
  dat.calib.pv <-
    calib_msm(data.ms = msebmtcal,
              data.raw = ebmtcal,
              j = 1,
              s = 0,
              t = 1826,
              tp.pred = tp.pred,
              calib.type = 'pv',
              pv.precalc = tp.pred,
              curve.type = "rcs",
              rcs.nk = 3)

  expect_type(dat.calib.pv, "list")
  expect_equal(class(dat.calib.pv), c("calib_pv", "calib_msm"))
  expect_length(dat.calib.pv[["plotdata"]], 6)
  expect_length(dat.calib.pv[["plotdata"]][["state3"]]$id, 2279)
  expect_length(dat.calib.pv[["plotdata"]][["state6"]]$id, 2279)
  expect_false(dat.calib.pv[["metadata"]]$CI)

})

