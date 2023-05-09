test_that("check calc_calib_mlr output", {

  ## Extract relevant predicted risks from tps0
  tp.pred <- dplyr::select(dplyr::filter(tps100, j == 3), any_of(paste("pstate", 1:6, sep = "")))

  ## Expect error if generate with CI
  expect_error(calc_calib_mlr(data.mstate = msebmtcal,
                              data.raw = ebmtcal,
                              j=3,
                              s=100,
                              t.eval = 1826,
                              tp.pred = tp.pred,
                              w.covs = c("year", "agecl", "proph", "match"),
                              CI = 95,
                              CI.R.boot = 5))

  ## Calculate observed event probabilities (run it on j = 3 and s = 100 so its a bit quicker, as smaller number of individuals)
  dat.calib.mlr <-
    calc_calib_mlr(data.mstate = msebmtcal,
                   data.raw = ebmtcal,
                   j=3,
                   s=100,
                   t.eval = 1826,
                   tp.pred = tp.pred,
                   w.covs = c("year", "agecl", "proph", "match"))

  str(dat.calib.mlr)
  expect_type(dat.calib.mlr, "list")
  expect_equal(class(dat.calib.mlr), "calib_mlr")
  expect_length(dat.calib.mlr, 2)
  expect_length(dat.calib.mlr[["plotdata"]], 4)
  expect_length(dat.calib.mlr[["plotdata"]][["state3"]]$id, 359)
  expect_length(dat.calib.mlr[["plotdata"]][["state6"]]$id, 359)
  expect_error(dat.calib.mlr[["plotdata"]][[6]])
  expect_length(dat.calib.mlr[["metadata"]], 1)

})
