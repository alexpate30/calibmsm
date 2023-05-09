test_that("check plot.calib_blr output (j = 1, s = 0)", {

  ## Extract relevant predicted risks from tps0
  tp.pred <- dplyr::select(dplyr::filter(tps0, j == 1), any_of(paste("pstate", 1:6, sep = "")))

  ## Calculate observed event probabilities
  dat.calib.blr <-
    calc_calib_blr(data.mstate = msebmtcal,
                   data.raw = ebmtcal,
                   j=1,
                   s=0,
                   t.eval = 1826,
                   tp.pred = tp.pred,
                   curve.type = "rcs",
                   rcs.nk = 3,
                   w.covs = c("year", "agecl", "proph", "match"))

  ## Plot calibration plots and run tests
  plot.object <- plot(dat.calib.blr, combine = TRUE, nrow = 2, ncol = 3)
  expect_equal(class(plot.object), c("gg", "ggplot", "ggarrange"))
  plot.object <- plot(dat.calib.blr, combine = FALSE, nrow = 2, ncol = 3)
  length(plot.object)
  expect_length(plot.object, 6)
  expect_type(plot.object, "list")

})

test_that("check plot.calib_blr output (j = 1, s = 0) with CI", {

  ## Extract relevant predicted risks from tps0
  tp.pred <- dplyr::select(dplyr::filter(tps0, j == 1), any_of(paste("pstate", 1:6, sep = "")))

  ## Calculate observed event probabilities
  dat.calib.blr <-
    calc_calib_blr(data.mstate = msebmtcal,
                   data.raw = ebmtcal,
                   j=1,
                   s=0,
                   t.eval = 1826,
                   tp.pred = tp.pred,
                   curve.type = "rcs",
                   rcs.nk = 3,
                   w.covs = c("year", "agecl", "proph", "match"),
                   CI = 95,
                   CI.R.boot = 5)

  ## Plot calibration plots and run tests
  plot.object <- plot(dat.calib.blr, combine = TRUE, nrow = 2, ncol = 3)
  expect_equal(class(plot.object), c("gg", "ggplot", "ggarrange"))
  plot.object <- plot(dat.calib.blr, combine = FALSE, nrow = 2, ncol = 3)
  length(plot.object)
  expect_length(plot.object, 6)
  expect_type(plot.object, "list")

})


test_that("check plot.calib_blr output (j = 3, s = 100)", {

  ## Extract relevant predicted risks from tps0
  tp.pred <- dplyr::select(dplyr::filter(tps100, j == 3), any_of(paste("pstate", 1:6, sep = "")))

  ## Calculate observed event probabilities
  dat.calib.blr <-
    calc_calib_blr(data.mstate = msebmtcal,
                   data.raw = ebmtcal,
                   j=3,
                   s=100,
                   t.eval = 1826,
                   tp.pred = tp.pred,
                   curve.type = "rcs",
                   rcs.nk = 3,
                   w.covs = c("year", "agecl", "proph", "match"))

  ## Plot calibration plots and run tests
  plot.object <- plot(dat.calib.blr, combine = TRUE, nrow = 2, ncol = 3)
  expect_equal(class(plot.object), c("gg", "ggplot", "ggarrange"))
  plot.object <- plot(dat.calib.blr, combine = FALSE, nrow = 2, ncol = 3)
  length(plot.object)
  expect_length(plot.object, 4)
  expect_type(plot.object, "list")

})

test_that("check plot.calib_mlr output (j = 3, s = 100)", {

  ## Extract relevant predicted risks from tps0
  tp.pred <- dplyr::select(dplyr::filter(tps100, j == 3), any_of(paste("pstate", 1:6, sep = "")))

  ## Calculate observed event probabilities
  dat.calib.mlr <-
    calc_calib_mlr(data.mstate = msebmtcal,
                   data.raw = ebmtcal,
                   j=3,
                   s=100,
                   t.eval = 1826,
                   tp.pred = tp.pred,
                   w.covs = c("year", "agecl", "proph", "match"))

  ## Plot calibration plots and run tests
  plot.object <- plot(dat.calib.mlr, combine = TRUE, nrow = 2, ncol = 3)
  expect_equal(class(plot.object), c("gg", "ggplot", "ggarrange"))
  plot.object <- plot(dat.calib.mlr, combine = FALSE, nrow = 2, ncol = 3)
  length(plot.object)
  expect_length(plot.object, 4)
  expect_type(plot.object, "list")
})
