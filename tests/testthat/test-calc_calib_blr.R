test_that("check calc_calib_blr output, (j = 1, s = 0)", {

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

  expect_type(dat.calib.blr, "list")
  expect_equal(class(dat.calib.blr), "calib_blr")
  expect_length(dat.calib.blr, 2)
  expect_length(dat.calib.blr[["plotdata"]], 6)
  expect_length(dat.calib.blr[["plotdata"]][[1]]$id, 1778)
  expect_length(dat.calib.blr[["plotdata"]][[6]]$id, 1778)
  expect_length(dat.calib.blr[["metadata"]], 3)
  expect_false(dat.calib.blr[["metadata"]]$CI)

})

test_that("check calc_calib_blr output, (j = 1, s = 0) CI", {

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

  expect_type(dat.calib.blr, "list")
  expect_equal(class(dat.calib.blr), "calib_blr")
  expect_length(dat.calib.blr, 2)
  expect_length(dat.calib.blr[["plotdata"]], 6)
  expect_equal(ncol(dat.calib.blr[["plotdata"]][[1]]), 5)
  expect_equal(ncol(dat.calib.blr[["plotdata"]][[6]]), 5)
  expect_length(dat.calib.blr[["plotdata"]][[1]]$id, 1778)
  expect_length(dat.calib.blr[["plotdata"]][[6]]$id, 1778)
  expect_length(dat.calib.blr[["metadata"]], 3)
  expect_equal(dat.calib.blr[["metadata"]]$CI, 95)

})

test_that("check calc_calib_blr output, (j = 3, s = 100)", {

  ## Extract relevant predicted risks from tps0
  tp.pred <- dplyr::select(dplyr::filter(tps0, j == 3), any_of(paste("pstate", 1:6, sep = "")))

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

  expect_type(dat.calib.blr, "list")
  expect_equal(class(dat.calib.blr), "calib_blr")
  expect_length(dat.calib.blr, 2)
  expect_length(dat.calib.blr[["plotdata"]], 4)
  expect_length(dat.calib.blr[["plotdata"]][["state3"]]$id, 359)
  expect_length(dat.calib.blr[["plotdata"]][["state6"]]$id, 359)
  expect_error(dat.calib.blr[["plotdata"]][[6]])
  expect_length(dat.calib.blr[["metadata"]], 3)
  expect_false(dat.calib.blr[["metadata"]]$CI)
  names(dat.calib.blr[["plotdata"]])

})
