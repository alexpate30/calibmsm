test_that("check calc_calib_blr output, (j = 1, s = 0), curve.type = rcs", {

  ## Extract relevant predicted risks from tps0
  tp.pred <- dplyr::select(dplyr::filter(tps0, j == 1), dplyr::any_of(paste("pstate", 1:6, sep = "")))

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
  expect_length(dat.calib.blr[["metadata"]], 8)
  expect_false(dat.calib.blr[["metadata"]]$CI)

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
                   w.stabilised = TRUE)

  expect_type(dat.calib.blr, "list")
  expect_equal(class(dat.calib.blr), "calib_blr")
  expect_length(dat.calib.blr, 2)
  expect_length(dat.calib.blr[["plotdata"]], 6)
  expect_length(dat.calib.blr[["plotdata"]][[1]]$id, 1778)
  expect_length(dat.calib.blr[["plotdata"]][[6]]$id, 1778)
  expect_length(dat.calib.blr[["metadata"]], 8)
  expect_false(dat.calib.blr[["metadata"]]$CI)

})


test_that("check calc_calib_blr output, (j = 1, s = 0), curve.type = loess", {

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
                   curve.type = "loess",
                   w.covs = c("year", "agecl", "proph", "match"))

  expect_type(dat.calib.blr, "list")
  expect_equal(class(dat.calib.blr), "calib_blr")
  expect_length(dat.calib.blr, 2)
  expect_length(dat.calib.blr[["plotdata"]], 6)
  expect_length(dat.calib.blr[["plotdata"]][[1]]$id, 1778)
  expect_length(dat.calib.blr[["plotdata"]][[6]]$id, 1778)
  expect_length(dat.calib.blr[["metadata"]], 8)
  expect_false(dat.calib.blr[["metadata"]]$CI)

  ## Calculate observed event probabilities
  dat.calib.blr <-
    calc_calib_blr(data.mstate = msebmtcal,
                   data.raw = ebmtcal,
                   j=1,
                   s=0,
                   t.eval = 1826,
                   tp.pred = tp.pred,
                   curve.type = "loess",
                   w.covs = c("year", "agecl", "proph", "match"),
                   w.stabilised = TRUE)

  expect_type(dat.calib.blr, "list")
  expect_equal(class(dat.calib.blr), "calib_blr")
  expect_length(dat.calib.blr, 2)
  expect_length(dat.calib.blr[["plotdata"]], 6)
  expect_length(dat.calib.blr[["plotdata"]][[1]]$id, 1778)
  expect_length(dat.calib.blr[["plotdata"]][[6]]$id, 1778)
  expect_length(dat.calib.blr[["metadata"]], 8)
  expect_false(dat.calib.blr[["metadata"]]$CI)

})


test_that("check calc_calib_blr output, (j = 1, s = 0), with CI", {

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
  expect_length(dat.calib.blr[["metadata"]], 8)
  expect_equal(dat.calib.blr[["metadata"]]$CI, 95)

})


test_that("check calc_calib_blr output, (j = 3, s = 100)", {

  ## Extract relevant predicted risks from tps100
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

  expect_type(dat.calib.blr, "list")
  expect_equal(class(dat.calib.blr), "calib_blr")
  expect_length(dat.calib.blr, 2)
  expect_length(dat.calib.blr[["plotdata"]], 4)
  expect_length(dat.calib.blr[["plotdata"]][["state3"]]$id, 359)
  expect_length(dat.calib.blr[["plotdata"]][["state6"]]$id, 359)
  expect_error(dat.calib.blr[["plotdata"]][[6]])
  expect_length(dat.calib.blr[["metadata"]], 8)
  expect_false(dat.calib.blr[["metadata"]]$CI)
  names(dat.calib.blr[["plotdata"]])

})


test_that("check calc_calib_blr output, (j = 1, s = 0), null covs", {

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
                   rcs.nk = 3)

  expect_type(dat.calib.blr, "list")
  expect_equal(class(dat.calib.blr), "calib_blr")
  expect_length(dat.calib.blr, 2)
  expect_length(dat.calib.blr[["plotdata"]], 6)
  expect_length(dat.calib.blr[["plotdata"]][[1]]$id, 1778)
  expect_length(dat.calib.blr[["plotdata"]][[6]]$id, 1778)
  expect_length(dat.calib.blr[["metadata"]], 8)
  expect_false(dat.calib.blr[["metadata"]]$CI)

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
                   w.stabilised = TRUE)

  expect_type(dat.calib.blr, "list")
  expect_equal(class(dat.calib.blr), "calib_blr")
  expect_length(dat.calib.blr, 2)
  expect_length(dat.calib.blr[["plotdata"]], 6)
  expect_length(dat.calib.blr[["plotdata"]][[1]]$id, 1778)
  expect_length(dat.calib.blr[["plotdata"]][[6]]$id, 1778)
  expect_length(dat.calib.blr[["metadata"]], 8)
  expect_false(dat.calib.blr[["metadata"]]$CI)

})


test_that("check calc_calib_blr output, (j = 1, s = 0),
          manual weights,
          manually define vector of predicted probabilities,
          manually define transition out", {

  ## Extract relevant predicted risks from tps0
  tp.pred <- dplyr::select(dplyr::filter(tps0, j == 1), any_of(paste("pstate", 1:6, sep = "")))

  ## Define t.eval
  t.eval <- 1826

  ## Extract data for plot manually
  ids.uncens <- ebmtcal %>%
    subset(dtcens > t.eval | (dtcens < t.eval & dtcens.s == 0)) %>%
    dplyr::pull(id)
  data.pred.plot <- tps0 %>%
    dplyr::filter(j == 1 & id %in% ids.uncens) %>%
    dplyr::select(any_of(paste("pstate", 1:6, sep = "")))

  ## Calculate manual weights
  weights.manual <-
    calc_weights(data.mstate = msebmtcal,
                 data.raw = ebmtcal,
                 t.eval = t.eval,
                 s = 0,
                 landmark.type = "all",
                 j = 1,
                 max.weight = 10,
                 stabilised = FALSE)

  ## Calculate observed event probabilities
  dat.calib.blr <-
    calc_calib_blr(data.mstate = msebmtcal,
                   data.raw = ebmtcal,
                   j=1,
                   s=0,
                   t.eval = t.eval,
                   tp.pred = tp.pred,
                   curve.type = "rcs",
                   rcs.nk = 3,
                   weights = weights.manual$ipcw,
                   data.pred.plot = data.pred.plot,
                   transitions.out = c(1,2,3,4,5,6))

  expect_type(dat.calib.blr, "list")
  expect_equal(class(dat.calib.blr), "calib_blr")
  expect_length(dat.calib.blr, 2)
  expect_length(dat.calib.blr[["plotdata"]], 6)
  expect_length(dat.calib.blr[["plotdata"]][[1]]$pred, 1778)
  expect_length(dat.calib.blr[["plotdata"]][[6]]$pred, 1778)
  expect_length(dat.calib.blr[["metadata"]], 8)
  expect_false(dat.calib.blr[["metadata"]]$CI)

})


test_that("check calc_calib_blr output, (j = 1, s = 0),
          with CI,
          manually define vector of predicted probabilities,
          manually define transition out", {

            ## Extract relevant predicted risks from tps0
            tp.pred <- dplyr::select(dplyr::filter(tps0, j == 1), any_of(paste("pstate", 1:6, sep = "")))

            ## Define t.eval
            t.eval <- 1826

            ## Extract data for plot manually
            ids.uncens <- ebmtcal %>%
              subset(dtcens > t.eval | (dtcens < t.eval & dtcens.s == 0)) %>%
              dplyr::pull(id)
            data.pred.plot <- tps0 %>%
              dplyr::filter(j == 1 & id %in% ids.uncens) %>%
              dplyr::select(any_of(paste("pstate", 1:6, sep = "")))

            ## Calculate observed event probabilities
            dat.calib.blr <-
              calc_calib_blr(data.mstate = msebmtcal,
                             data.raw = ebmtcal,
                             j=1,
                             s=0,
                             t.eval = t.eval,
                             tp.pred = tp.pred,
                             curve.type = "rcs",
                             rcs.nk = 3,
                             w.covs = c("year", "agecl", "proph", "match"),
                             CI = 95,
                             CI.R.boot = 5,
                             data.pred.plot = data.pred.plot,
                             transitions.out = c(1,2,3,4,5,6))

            expect_type(dat.calib.blr, "list")
            expect_equal(class(dat.calib.blr), "calib_blr")
            expect_length(dat.calib.blr, 2)
            expect_length(dat.calib.blr[["plotdata"]], 6)
            expect_length(dat.calib.blr[["plotdata"]][[1]]$pred, 1778)
            expect_length(dat.calib.blr[["plotdata"]][[6]]$pred, 1778)
            expect_length(dat.calib.blr[["metadata"]], 8)
            expect_equal(dat.calib.blr[["metadata"]]$CI, 95)

          })


test_that("test warnings and errors", {

  ## Extract relevant predicted risks from tps0
  tp.pred <- dplyr::select(dplyr::filter(tps0, j == 3), any_of(paste("pstate", 1:6, sep = "")))

  ## Calculate observed event probabilities
  expect_error(
    calc_calib_blr(data.mstate = msebmtcal,
                   data.raw = ebmtcal,
                   j=1,
                   s=0,
                   t.eval = 1826,
                   tp.pred = tp.pred,
                   curve.type = "rcs",
                   rcs.nk = 3,
                   w.covs = c("year", "agecl", "proph", "match"),
                   transitions.out = c(1,2,3,4))
  )

  ## Calculate observed event probabilities
  weights.manual <-
    calc_weights(data.mstate = msebmtcal,
                 data.raw = ebmtcal,
                 t.eval = 1826,
                 s = 0,
                 landmark.type = "all",
                 j = 1,
                 max.weight = 10,
                 stabilised = FALSE)$ipcw[-1]

  expect_error(
    calc_calib_blr(data.mstate = msebmtcal,
                   data.raw = ebmtcal,
                   j=1,
                   s=0,
                   t.eval = 1826,
                   tp.pred = tp.pred,
                   curve.type = "rcs",
                   rcs.nk = 3,
                   weights = weights.manual)
  )

})

