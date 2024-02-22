test_that("check plot.calib_msm output (j = 1, s = 0)", {

  ## Extract relevant predicted risks from tps0
  tp.pred <- dplyr::select(dplyr::filter(tps0, j == 1), any_of(paste("pstate", 1:6, sep = "")))

  ## Calculate observed event probabilities
  dat.calib.blr <-
    calib_msm(data.mstate = msebmtcal,
             data.raw = ebmtcal,
             j=1,
             s=0,
             t = 1826,
             tp.pred = tp.pred,
             calib.type = "blr",
             curve.type = "rcs",
             rcs.nk = 3,
             w.covs = c("year", "agecl", "proph", "match"))

  ## Plot calibration plots and run tests
  plot.object <- plot(dat.calib.blr, combine = TRUE, nrow = 2, ncol = 3)
  expect_equal(class(plot.object), c("gg", "ggplot", "ggarrange"))
  plot.object <- plot(dat.calib.blr, combine = FALSE, nrow = 2, ncol = 3)
  expect_length(plot.object, 6)
  expect_type(plot.object, "list")

  ## Plot calibration plots and run tests with marginal density plots
  plot.object <- plot(dat.calib.blr, combine = TRUE, nrow = 2, ncol = 3, marg.density = TRUE, marg.density.size = 1)
  expect_length(plot.object, 6)
  expect_equal(class(plot.object), c("gtable", "gTree", "grob", "gDesc"))

  ## Plot calibration plots and run tests with marginal rug plots
  plot.object <- plot(dat.calib.blr, combine = TRUE, nrow = 2, ncol = 3, marg.rug = TRUE)
  expect_equal(class(plot.object), c("gg", "ggplot", "ggarrange"))

})

test_that("check plot.calib_msm output (j = 1, s = 0) with CI", {

  ## Extract relevant predicted risks from tps0
  tp.pred <- dplyr::select(dplyr::filter(tps0, j == 1), any_of(paste("pstate", 1:6, sep = "")))

  ## Calculate observed event probabilities
  dat.calib.blr <-
    calib_msm(data.mstate = msebmtcal,
             data.raw = ebmtcal,
             j=1,
             s=0,
             t = 1826,
             tp.pred = tp.pred,
             calib.type = "blr",
             curve.type = "rcs",
             rcs.nk = 3,
             w.covs = c("year", "agecl", "proph", "match"),
             CI = 95,
             CI.R.boot = 5)

  ## Plot calibration plots and run tests without marginal density plots
  plot.object <- plot(dat.calib.blr, combine = TRUE, nrow = 2, ncol = 3)
  expect_equal(class(plot.object), c("gg", "ggplot", "ggarrange"))
  plot.object <- plot(dat.calib.blr, combine = FALSE, nrow = 2, ncol = 3)
  expect_length(plot.object, 6)
  expect_type(plot.object, "list")

  ## Plot calibration plots and run tests with marginal density plots
  plot.object <- plot(dat.calib.blr, combine = TRUE, nrow = 2, ncol = 3, marg.density = TRUE, marg.density.size = 1)
  expect_equal(class(plot.object), c("gtable", "gTree", "grob", "gDesc"))

  ## Plot calibration plots and run tests with marginal rug plots
  plot.object <- plot(dat.calib.blr, combine = TRUE, nrow = 2, ncol = 3, marg.rug = TRUE)
  expect_equal(class(plot.object), c("gg", "ggplot", "ggarrange"))

})


test_that("check plot.calib_msm output (j = 3, s = 100)", {

  ## Extract relevant predicted risks from tps0
  tp.pred <- dplyr::select(dplyr::filter(tps100, j == 3), any_of(paste("pstate", 1:6, sep = "")))

  ## Calculate observed event probabilities
  dat.calib.blr <-
    calib_msm(data.mstate = msebmtcal,
             data.raw = ebmtcal,
             j=3,
             s=100,
             t = 1826,
             tp.pred = tp.pred,
             calib.type = "blr",
             curve.type = "rcs",
             rcs.nk = 3,
             w.covs = c("year", "agecl", "proph", "match"))

  ## Plot calibration plots and run tests
  plot.object <- plot(dat.calib.blr, combine = TRUE, nrow = 2, ncol = 3)
  expect_equal(class(plot.object), c("gg", "ggplot", "ggarrange"))
  plot.object <- plot(dat.calib.blr, combine = FALSE, nrow = 2, ncol = 3)
  expect_length(plot.object, 4)
  expect_type(plot.object, "list")

})


test_that("check plot.calib_pv output (j = 1, s = 0)", {

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

  ## Calculate observed event probabilities
  dat.calib.pv <-
    suppressWarnings(calib_msm(data.mstate = msebmtcal,
                              data.raw = ebmtcal,
                              j=1,
                              s=0,
                              t = 1826,
                              tp.pred = tp.pred,
                              calib.type = "pv",
                              curve.type = "rcs",
                              rcs.nk = 3))

  ## Plot calibration plots and run tests
  plot.object <- plot(dat.calib.pv, combine = TRUE)
  expect_equal(class(plot.object), c("gg", "ggplot", "ggarrange"))
  plot.object <- plot(dat.calib.pv, combine = FALSE)
  expect_length(plot.object, 6)
  expect_type(plot.object, "list")

})

test_that("check plot.calib_pv output (j = 1, s = 0) with CI", {

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

  ## Calculate observed event probabilities
  dat.calib.pv <-
    suppressWarnings(calib_msm(data.mstate = msebmtcal,
                              data.raw = ebmtcal,
                              j=1,
                              s=0,
                              t = 1826,
                              tp.pred = tp.pred,
                              calib.type = "pv",
                              curve.type = "rcs",
                              rcs.nk = 3,
                              CI = 95,
                              CI.type = "parametric"))

  ## Plot calibration plots and run tests
  plot.object <- plot(dat.calib.pv, combine = TRUE)
  expect_equal(class(plot.object), c("gg", "ggplot", "ggarrange"))
  plot.object <- plot(dat.calib.pv, combine = FALSE)
  expect_length(plot.object, 6)
  expect_type(plot.object, "list")

})


test_that("check plot.calib_pv output (j = 3, s = 100) with CI", {

  ## Reduce to 500 individuals
  # Extract the predicted transition probabilities out of state j = 1 for first 500 individuals
  tp.pred <- tps0 |>
    dplyr::filter(id %in% 1:500) |>
    dplyr::filter(j == 3) |>
    dplyr::select(any_of(paste("pstate", 1:6, sep = "")))
  # Reduce ebmtcal to first 500 individuals
  ebmtcal <- ebmtcal |> dplyr::filter(id %in% 1:500)
  # Reduce msebmtcal.cmprsk to first 100 individuals
  msebmtcal <- msebmtcal |> dplyr::filter(id %in% 1:500)

  ## Calculate observed event probabilities
  dat.calib.pv <-
    calib_msm(data.mstate = msebmtcal,
             data.raw = ebmtcal,
             j=3,
             s=100,
             t = 1826,
             tp.pred = tp.pred,
             calib.type = "pv",
             curve.type = "rcs",
             rcs.nk = 3,
             CI = 95,
             CI.type = "parametric")

  ## Plot calibration plots and run tests
  plot.object <- plot(dat.calib.pv, combine = TRUE)
  expect_equal(class(plot.object), c("gg", "ggplot", "ggarrange"))
  plot.object <- plot(dat.calib.pv, combine = FALSE)
  expect_length(plot.object, 4)
  expect_type(plot.object, "list")

})



test_that("check plot.calib_mlr output (j = 1, s = 0)", {

  ## Reduce to 500 individuals
  # Extract the predicted transition probabilities out of state j = 1 for first 500 individuals
  tp.pred <- tps0 |>
    dplyr::filter(id %in% 1:500) |>
    dplyr::filter(j == 1) |>
    dplyr::select(any_of(paste("pstate", 1:6, sep = "")))
  # Reduce ebmtcal to first 500 individuals
  ebmtcal <- ebmtcal |> dplyr::filter(id %in% 1:500)
  # Reduce msebmtcal.cmprsk to first 100 individuals
  msebmtcal <- msebmtcal |> dplyr::filter(id %in% 1:500)

  # ## Extract relevant predicted risks from tps0
  # tp.pred <- dplyr::select(dplyr::filter(tps100, j == 3), dplyr::any_of(paste("pstate", 1:6, sep = "")))

  ## Calculate observed event probabilities
  suppressWarnings(
    dat.calib.mlr <-
      calib_msm(data.mstate = msebmtcal,
               data.raw = ebmtcal,
               j=1,
               s=0,
               t = 1826,
               tp.pred = tp.pred,
               calib.type = "mlr",
               w.covs = c("year", "agecl", "proph", "match"))
  )

  ## Plot calibration plots and run tests
  plot.object <- plot(dat.calib.mlr, combine = TRUE, nrow = 2, ncol = 3)
  expect_equal(class(plot.object), c("gg", "ggplot", "ggarrange"))
  plot.object <- plot(dat.calib.mlr, combine = FALSE, nrow = 2, ncol = 3)
  expect_length(plot.object, 6)
  expect_type(plot.object, "list")

})


test_that("check plot.calib_mlr output (j = 3, s = 100)", {

  ## Extract relevant predicted risks from tps0
  tp.pred <- dplyr::select(dplyr::filter(tps100, j == 3), dplyr::any_of(paste("pstate", 1:6, sep = "")))

  ## Calculate observed event probabilities
  suppressWarnings(
    dat.calib.mlr <-
      calib_msm(data.mstate = msebmtcal,
               data.raw = ebmtcal,
               j=3,
               s=100,
               t = 1826,
               tp.pred = tp.pred,
               calib.type = "mlr",
               w.covs = c("year", "agecl", "proph", "match"))
  )

  ## Plot calibration plots and run tests
  plot.object <- plot(dat.calib.mlr, combine = TRUE, nrow = 2, ncol = 3)
  expect_equal(class(plot.object), c("gg", "ggplot", "ggarrange"))
  plot.object <- plot(dat.calib.mlr, combine = FALSE, nrow = 2, ncol = 3)
  expect_length(plot.object, 4)
  expect_type(plot.object, "list")

})
