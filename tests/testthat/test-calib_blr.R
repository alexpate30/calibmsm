###
### Tests for calibration curves produced using BLR-IPCW (calib_type = 'blr')
###

test_that("check calib_msm output, (j = 1, s = 0), curve_type = rcs, stabilised vs unstabilised", {

  ## Extract relevant predicted risks from tps0
  tp_pred <- dplyr::select(dplyr::filter(tps0, j == 1), dplyr::any_of(paste("pstate", 1:6, sep = "")))

  ###
  ### Calculate observed event probabilities
  dat_calib_blr <-
    calib_msm(data_ms = msebmtcal,
             data_raw = ebmtcal,
             j=1,
             s=0,
             t = 1826,
             tp_pred = tp_pred, calib_type = 'blr',
             curve_type = "rcs",
             rcs_nk = 3,
             w_covs = c("year", "agecl", "proph", "match"))

  expect_type(dat_calib_blr, "list")
  expect_equal(class(dat_calib_blr), c("calib_blr", "calib_msm"))
  expect_length(dat_calib_blr[["mean"]], 6)
  expect_length(dat_calib_blr[["plotdata"]], 6)
  expect_length(dat_calib_blr[["plotdata"]][[1]]$id, 1778)
  expect_length(dat_calib_blr[["plotdata"]][[6]]$id, 1778)
  expect_false(dat_calib_blr[["metadata"]]$CI)
  expect_no_error(summary(dat_calib_blr))

  ###
  ### Calculate observed event probabilities with stabilised weights
  dat_calib_blr_stab <-
    calib_msm(data_ms = msebmtcal,
             data_raw = ebmtcal,
             j=1,
             s=0,
             t = 1826,
             tp_pred = tp_pred, calib_type = 'blr',
             curve_type = "rcs",
             rcs_nk = 3,
             w_covs = c("year", "agecl", "proph", "match"),
             w_stabilised = TRUE)

  expect_type(dat_calib_blr_stab, "list")
  expect_equal(class(dat_calib_blr_stab), c("calib_blr", "calib_msm"))
  expect_length(dat_calib_blr_stab[["plotdata"]], 6)
  expect_length(dat_calib_blr_stab[["plotdata"]][[1]]$id, 1778)
  expect_length(dat_calib_blr_stab[["plotdata"]][[6]]$id, 1778)
  expect_false(dat_calib_blr_stab[["metadata"]]$CI)

  ### Check answer is same whether stabilisation used or not
  expect_equal(dat_calib_blr[["plotdata"]][[1]], dat_calib_blr_stab[["plotdata"]][[1]])

})


test_that("check calib_msm output, (j = 1, s = 0), curve_type = loess", {

  ## Extract relevant predicted risks from tps0
  tp_pred <- dplyr::select(dplyr::filter(tps0, j == 1), any_of(paste("pstate", 1:6, sep = "")))

  ## Calculate observed event probabilities
  dat_calib_blr <-
    calib_msm(data_ms = msebmtcal,
             data_raw = ebmtcal,
             j=1,
             s=0,
             t = 1826,
             tp_pred = tp_pred, calib_type = 'blr',
             curve_type = "loess",
             w_covs = c("year", "agecl", "proph", "match"))

  expect_type(dat_calib_blr, "list")
  expect_equal(class(dat_calib_blr), c("calib_blr", "calib_msm"))
  expect_length(dat_calib_blr[["mean"]], 6)
  expect_length(dat_calib_blr[["plotdata"]], 6)
  expect_length(dat_calib_blr[["plotdata"]][[1]]$id, 1778)
  expect_length(dat_calib_blr[["plotdata"]][[6]]$id, 1778)
  expect_false(dat_calib_blr[["metadata"]]$CI)
  expect_no_error(summary(dat_calib_blr))

  ## Calculate observed event probabilities
  dat_calib_blr_stab <-
    calib_msm(data_ms = msebmtcal,
             data_raw = ebmtcal,
             j=1,
             s=0,
             t = 1826,
             tp_pred = tp_pred, calib_type = 'blr',
             curve_type = "loess",
             w_covs = c("year", "agecl", "proph", "match"),
             w_stabilised = TRUE)

  expect_type(dat_calib_blr_stab, "list")
  expect_equal(class(dat_calib_blr_stab), c("calib_blr", "calib_msm"))
  expect_length(dat_calib_blr_stab[["plotdata"]], 6)
  expect_length(dat_calib_blr_stab[["plotdata"]][[1]]$id, 1778)
  expect_length(dat_calib_blr_stab[["plotdata"]][[6]]$id, 1778)
  expect_false(dat_calib_blr_stab[["metadata"]]$CI)

  ### Check answer is same whether stabilisation used or not
  expect_equal(dat_calib_blr[["plotdata"]][[1]], dat_calib_blr_stab[["plotdata"]][[1]])

  ## Calculate observed event probabilities
  dat_calib_blr_w_function <-
    calib_msm(data_ms = msebmtcal,
             data_raw = ebmtcal,
             j=1,
             s=0,
             t = 1826,
             tp_pred = tp_pred, calib_type = 'blr',
             curve_type = "loess",
             w_function = calc_weights,
             w_covs = c("year", "agecl", "proph", "match"))

  expect_type(dat_calib_blr_w_function, "list")
  expect_equal(class(dat_calib_blr_w_function), c("calib_blr", "calib_msm"))
  expect_length(dat_calib_blr_w_function[["mean"]], 6)
  expect_length(dat_calib_blr_w_function[["plotdata"]], 6)
  expect_length(dat_calib_blr_w_function[["plotdata"]][[1]]$id, 1778)
  expect_length(dat_calib_blr_w_function[["plotdata"]][[6]]$id, 1778)
  expect_false(dat_calib_blr_w_function[["metadata"]]$CI)

  ### Check answer is same whether stabilisation used or not
  expect_equal(dat_calib_blr[["plotdata"]][[1]], dat_calib_blr_w_function[["plotdata"]][[1]])
  expect_equal(dat_calib_blr[["plotdata"]][[6]], dat_calib_blr_w_function[["plotdata"]][[6]])

})


test_that("check calib_msm output, (j = 1, s = 0), with CI", {

  skip_on_cran()

  ## Extract relevant predicted risks from tps0
  tp_pred <- dplyr::select(dplyr::filter(tps0, j == 1), any_of(paste("pstate", 1:6, sep = "")))

  ## Calculate observed event probabilities no CI
  dat_calib_blr_noCI <-
    calib_msm(data_ms = msebmtcal,
             data_raw = ebmtcal,
             j=1,
             s=0,
             t = 1826,
             tp_pred = tp_pred, calib_type = 'blr',
             curve_type = "rcs",
             rcs_nk = 3,
             w_covs = c("year", "agecl", "proph", "match"))

  ## Calculate observed event probabilities with CI
  dat_calib_blr <-
    calib_msm(data_ms = msebmtcal,
             data_raw = ebmtcal,
             j=1,
             s=0,
             t = 1826,
             tp_pred = tp_pred, calib_type = 'blr',
             curve_type = "rcs",
             rcs_nk = 3,
             w_covs = c("year", "agecl", "proph", "match"),
             CI = 95,
             CI_R_boot = 5)

  expect_type(dat_calib_blr, "list")
  expect_equal(class(dat_calib_blr), c("calib_blr", "calib_msm"))
  expect_length(dat_calib_blr[["mean"]], 6)
  expect_length(dat_calib_blr[["plotdata"]], 6)
  expect_equal(ncol(dat_calib_blr[["plotdata"]][[1]]), 5)
  expect_equal(ncol(dat_calib_blr[["plotdata"]][[6]]), 5)
  expect_length(dat_calib_blr[["plotdata"]][[1]]$id, 1778)
  expect_length(dat_calib_blr[["plotdata"]][[6]]$id, 1778)
  expect_equal(dat_calib_blr[["metadata"]]$CI, 95)
  expect_no_error(summary(dat_calib_blr))

  expect_equal(dat_calib_blr_noCI[["plotdata"]][[1]]$obs, dat_calib_blr[["plotdata"]][[1]]$obs)
  expect_equal(dat_calib_blr_noCI[["plotdata"]][[6]]$obs, dat_calib_blr[["plotdata"]][[6]]$obs)

})


test_that("check calib_msm output, (j = 3, s = 100)", {

  ## Extract relevant predicted risks from tps100
  tp_pred <- dplyr::select(dplyr::filter(tps100, j == 3), any_of(paste("pstate", 1:6, sep = "")))

  ## Calculate observed event probabilities
  dat_calib_blr <-
    calib_msm(data_ms = msebmtcal,
             data_raw = ebmtcal,
             j=3,
             s=100,
             t = 1826,
             tp_pred = tp_pred, calib_type = 'blr',
             curve_type = "rcs",
             rcs_nk = 3,
             w_covs = c("year", "agecl", "proph", "match"))

  expect_type(dat_calib_blr, "list")
  expect_equal(class(dat_calib_blr), c("calib_blr", "calib_msm"))
  expect_length(dat_calib_blr[["mean"]], 4)
  expect_length(dat_calib_blr[["plotdata"]], 4)
  expect_length(dat_calib_blr[["plotdata"]][["state3"]]$id, 359)
  expect_length(dat_calib_blr[["plotdata"]][["state6"]]$id, 359)
  expect_error(dat_calib_blr[["plotdata"]][[6]])
  expect_false(dat_calib_blr[["metadata"]]$CI)
  names(dat_calib_blr[["plotdata"]])

})


test_that("check calib_msm output, (j = 1, s = 0), null covs", {

  ## Extract relevant predicted risks from tps0
  tp_pred <- dplyr::select(dplyr::filter(tps0, j == 1), any_of(paste("pstate", 1:6, sep = "")))

  ## Calculate observed event probabilities
  dat_calib_blr <-
    calib_msm(data_ms = msebmtcal,
             data_raw = ebmtcal,
             j=1,
             s=0,
             t = 1826,
             tp_pred = tp_pred, calib_type = 'blr',
             curve_type = "rcs",
             rcs_nk = 3)

  expect_type(dat_calib_blr, "list")
  expect_equal(class(dat_calib_blr), c("calib_blr", "calib_msm"))
  expect_length(dat_calib_blr[["mean"]], 6)
  expect_length(dat_calib_blr[["plotdata"]], 6)
  expect_length(dat_calib_blr[["plotdata"]][[1]]$id, 1778)
  expect_length(dat_calib_blr[["plotdata"]][[6]]$id, 1778)
  expect_false(dat_calib_blr[["metadata"]]$CI)

  ## Calculate observed event probabilities
  dat_calib_blr <-
    calib_msm(data_ms = msebmtcal,
             data_raw = ebmtcal,
             j=1,
             s=0,
             t = 1826,
             tp_pred = tp_pred, calib_type = 'blr',
             curve_type = "rcs",
             rcs_nk = 3,
             w_covs = c("year", "agecl", "proph", "match"),
             w_stabilised = TRUE)

  expect_type(dat_calib_blr, "list")
  expect_equal(class(dat_calib_blr), c("calib_blr", "calib_msm"))
  expect_length(dat_calib_blr[["mean"]], 6)
  expect_length(dat_calib_blr[["plotdata"]], 6)
  expect_length(dat_calib_blr[["plotdata"]][[1]]$id, 1778)
  expect_length(dat_calib_blr[["plotdata"]][[6]]$id, 1778)
  expect_false(dat_calib_blr[["metadata"]]$CI)

})


### Run tests for manually inputted predicted probailities
test_that("check calib_msm output, (j = 1, s = 0), curve_type = loess, CI_type = bootstrap", {

  skip_on_cran()

  ## Extract relevant predicted risks from tps0 for creating plots
  tp_pred <- dplyr::select(dplyr::filter(tps0, j == 1), any_of(paste("pstate", 1:6, sep = "")))

  ### Create an object of only 50 observations over which to plot, which we specify manually
  id_lmk <- 1:50
  tp_pred_plot <- tps0 |>
    dplyr::filter(id %in% id_lmk) |>
    dplyr::filter(j == 1) |>
    dplyr::select(any_of(paste("pstate", 1:6, sep = "")))

  ## No confidence interval
  dat_calib_blr <- calib_msm(data_ms = msebmtcal,
                            data_raw = ebmtcal,
                            j = 1,
                            s = 0,
                            t = 1826,
                            tp_pred = tp_pred,
                            calib_type = 'blr',
                            curve_type = "loess",
                            tp_pred_plot = tp_pred_plot, transitions_out = NULL)

  ## Should be one less column in plotdata (no patient ids)
  expect_equal(class(dat_calib_blr), c("calib_blr", "calib_msm"))
  expect_equal(ncol(dat_calib_blr[["plotdata"]][[1]]), 2)
  expect_equal(nrow(dat_calib_blr[["plotdata"]][[1]]), 50)
  expect_no_error(summary(dat_calib_blr))

  ## With confidence interval
  dat_calib_blr <- calib_msm(data_ms = msebmtcal,
                            data_raw = ebmtcal,
                            j = 1,
                            s = 0,
                            t = 1826,
                            tp_pred = tp_pred,
                            calib_type = 'blr',
                            curve_type = "loess",
                            CI = 95,
                            CI_R_boot = 3,
                            tp_pred_plot = tp_pred_plot, transitions_out = NULL)

  ## Should be one less column in plotdata (no patient ids)
  expect_equal(class(dat_calib_blr), c("calib_blr", "calib_msm"))
  expect_equal(ncol(dat_calib_blr[["plotdata"]][[1]]), 4)
  expect_equal(nrow(dat_calib_blr[["plotdata"]][[1]]), 50)
  expect_no_error(summary(dat_calib_blr))

})


### Run tests for warnings with small cohort
test_that("Test warnings for bootstrapping with small cohort", {

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

  ## No confidence interval
  expect_warning(calib_msm(data_ms = msebmtcal,
                          data_raw = ebmtcal,
                          j = 1,
                          s = 0,
                          t = 1826,
                          tp_pred = tp_pred,
                          calib_type = 'blr',
                          curve_type = "loess",
                          CI = 95,
                          CI_R_boot = 3,
                          transitions_out = c(1)))

})


test_that("check calib_msm output, (j = 1, s = 0),
          manual weights,
          manually define vector of predicted probabilities,
          manually define transition out,
          estimate curves using rcs", {

            ## Extract relevant predicted risks from tps0
            tp_pred <- dplyr::select(dplyr::filter(tps0, j == 1), any_of(paste("pstate", 1:6, sep = "")))

            ## Define t
            t <- 1826

            ## Extract data for plot manually
            ids_uncens <- ebmtcal |>
              subset(dtcens > t | (dtcens < t & dtcens_s == 0)) |>
              dplyr::pull(id)
            tp_pred_plot <- tps0 |>
              dplyr::filter(j == 1 & id %in% ids_uncens) |>
              dplyr::select(any_of(paste("pstate", 1:6, sep = "")))

            ## Calculate manual weights
            weights_manual <-
              calc_weights(data_ms = msebmtcal,
                           data_raw = ebmtcal,
                           t = 1826,
                           s = 0,
                           landmark_type = "state",
                           j = 1,
                           max_weight = 10,
                           stabilised = FALSE)

            ## Calculate observed event probabilities using weights_manual
            dat_calib_blr_w_manual <-
              calib_msm(data_ms = msebmtcal,
                       data_raw = ebmtcal,
                       j=1,
                       s=0,
                       t = 1826,
                       tp_pred = tp_pred, calib_type = 'blr',
                       curve_type = "rcs",
                       rcs_nk = 3,
                       weights = weights_manual$ipcw,
                       tp_pred_plot = tp_pred_plot,
                       transitions_out = c(1,2,3,4,5,6))

            expect_type(dat_calib_blr_w_manual, "list")
            expect_equal(class(dat_calib_blr_w_manual), c("calib_blr", "calib_msm"))
            expect_length(dat_calib_blr_w_manual[["mean"]], 6)
            expect_length(dat_calib_blr_w_manual[["plotdata"]], 6)
            expect_length(dat_calib_blr_w_manual[["plotdata"]][[1]]$pred, 1778)
            expect_length(dat_calib_blr_w_manual[["plotdata"]][[6]]$pred, 1778)
            expect_false(dat_calib_blr_w_manual[["metadata"]]$CI)

          })

test_that("check calib_msm output, (j = 1, s = 0),
          manual weights,
          manually define vector of predicted probabilities,
          manually define transition out,
          estimate curves using loess", {

            ## Extract relevant predicted risks from tps0
            tp_pred <- dplyr::select(dplyr::filter(tps0, j == 1), any_of(paste("pstate", 1:6, sep = "")))

            ## Define t
            t <- 1826

            ## Extract data for plot manually
            ids_uncens <- ebmtcal |>
              subset(dtcens > t | (dtcens < t & dtcens_s == 0)) |>
              dplyr::pull(id)
            tp_pred_plot <- tps0 |>
              dplyr::filter(j == 1 & id %in% ids_uncens) |>
              dplyr::select(any_of(paste("pstate", 1:6, sep = "")))

            ## Calculate manual weights
            weights_manual <-
              calc_weights(data_ms = msebmtcal,
                           data_raw = ebmtcal,
                           t = 1826,
                           s = 0,
                           landmark_type = "state",
                           j = 1,
                           max_weight = 10,
                           stabilised = FALSE)

            ## Calculate observed event probabilities using weights_manual
            dat_calib_blr_w_manual <-
              calib_msm(data_ms = msebmtcal,
                       data_raw = ebmtcal,
                       j=1,
                       s=0,
                       t = 1826,
                       tp_pred = tp_pred, calib_type = 'blr',
                       curve_type = "loess",
                       rcs_nk = 3,
                       weights = weights_manual$ipcw,
                       tp_pred_plot = tp_pred_plot,
                       transitions_out = c(1,2,3,4,5,6))

            expect_type(dat_calib_blr_w_manual, "list")
            expect_equal(class(dat_calib_blr_w_manual), c("calib_blr", "calib_msm"))
            expect_length(dat_calib_blr_w_manual[["mean"]], 6)
            expect_length(dat_calib_blr_w_manual[["plotdata"]], 6)
            expect_length(dat_calib_blr_w_manual[["plotdata"]][[1]]$pred, 1778)
            expect_length(dat_calib_blr_w_manual[["plotdata"]][[6]]$pred, 1778)
            expect_false(dat_calib_blr_w_manual[["metadata"]]$CI)

          })

test_that("check calib_msm output, (j = 1, s = 0),
          with CI,
          manually define vector of predicted probabilities,
          manually define transition out", {

            skip_on_cran()

            ## Extract relevant predicted risks from tps0
            tp_pred <- dplyr::select(dplyr::filter(tps0, j == 1), any_of(paste("pstate", 1:6, sep = "")))

            ## Define t
            t <- 1826

            ## Extract data for plot manually
            ids_uncens <- ebmtcal |>
              subset(dtcens > t | (dtcens < t & dtcens_s == 0)) |>
              dplyr::pull(id)
            tp_pred_plot <- tps0 |>
              dplyr::filter(j == 1 & id %in% ids_uncens) |>
              dplyr::select(any_of(paste("pstate", 1:6, sep = "")))

            ## Calculate observed event probabilities
            dat_calib_blr <-
              calib_msm(data_ms = msebmtcal,
                       data_raw = ebmtcal,
                       j=1,
                       s=0,
                       t = 1826,
                       tp_pred = tp_pred, calib_type = 'blr',
                       curve_type = "rcs",
                       rcs_nk = 3,
                       w_covs = c("year", "agecl", "proph", "match"),
                       CI = 95,
                       CI_R_boot = 5,
                       tp_pred_plot = tp_pred_plot,
                       transitions_out = c(1,2,3,4,5,6))

            expect_type(dat_calib_blr, "list")
            expect_equal(class(dat_calib_blr), c("calib_blr", "calib_msm"))
            expect_length(dat_calib_blr[["mean"]], 6)
            expect_length(dat_calib_blr[["plotdata"]], 6)
            expect_length(dat_calib_blr[["plotdata"]][[1]]$pred, 1778)
            expect_length(dat_calib_blr[["plotdata"]][[6]]$pred, 1778)
            expect_equal(dat_calib_blr[["metadata"]]$CI, 95)

          })

test_that("check calib_msm output, (j = 1, s = 0),
          Manually define function to estimate weights", {

            skip_on_cran()

            ## Extract relevant predicted risks from tps0
            tp_pred <- dplyr::select(dplyr::filter(tps0, j == 1), dplyr::any_of(paste("pstate", 1:6, sep = "")))

            ###
            ### Calculate observed event probabilities
            dat_calib_blr <-
              calib_msm(data_ms = msebmtcal,
                       data_raw = ebmtcal,
                       j=1,
                       s=0,
                       t = 1826,
                       tp_pred = tp_pred, calib_type = 'blr',
                       curve_type = "rcs",
                       rcs_nk = 3,
                       w_covs = c("year", "agecl", "proph", "match"))

            ###
            ### Calculate manual weights
            weights_manual <-
              calc_weights(data_ms = msebmtcal,
                           data_raw = ebmtcal,
                           covs =  c("year", "agecl", "proph", "match"),
                           t = 1826,
                           s = 0,
                           landmark_type = "state",
                           j = 1,
                           max_weight = 10,
                           stabilised = FALSE)

            ###
            ### Calculate observed event probabilities same function as internal procedure, and check it agrees with dat_calib_blr
            dat_calib_blr_w_manual <-
              calib_msm(data_ms = msebmtcal,
                       data_raw = ebmtcal,
                       j=1,
                       s=0,
                       t = 1826,
                       tp_pred = tp_pred, calib_type = 'blr',
                       curve_type = "rcs",
                       rcs_nk = 3,
                       weights = weights_manual$ipcw)

            expect_equal(dat_calib_blr[["plotdata"]][[1]], dat_calib_blr_w_manual[["plotdata"]][[1]])

            ###
            ### Calculate observed event probabilities using an incorrect vector of weights, and see if its different from dat_calib_blr
            dat_calib_blr_w_manual <-
              calib_msm(data_ms = msebmtcal,
                       data_raw = ebmtcal,
                       j=1,
                       s=0,
                       t = 1826,
                       tp_pred = tp_pred, calib_type = 'blr',
                       curve_type = "rcs",
                       rcs_nk = 3,
                       weights = rep(1,nrow(weights_manual)))

            expect_false(any(dat_calib_blr[["plotdata"]][[1]]$obs == dat_calib_blr_w_manual[["plotdata"]][[1]]$obs))

            ###
            ### Calculate observed event probabilities with w_function, where calc_weights_manual = calc_weights (exactly same as internal procedure)
            calc_weights_manual <- calc_weights

            dat_calib_blr_w_function <-
              calib_msm(data_ms = msebmtcal,
                       data_raw = ebmtcal,
                       j=1,
                       s=0,
                       t = 1826,
                       tp_pred = tp_pred, calib_type = 'blr',
                       curve_type = "rcs",
                       rcs_nk = 3,
                       w_function = calc_weights_manual,
                       w_covs = c("year", "agecl", "proph", "match"))

            expect_type(dat_calib_blr_w_function, "list")
            expect_equal(class(dat_calib_blr_w_function), c("calib_blr", "calib_msm"))
            expect_length(dat_calib_blr_w_function[["mean"]], 6)
            expect_length(dat_calib_blr_w_function[["plotdata"]], 6)
            expect_length(dat_calib_blr_w_function[["plotdata"]][[1]]$id, 1778)
            expect_length(dat_calib_blr_w_function[["plotdata"]][[6]]$id, 1778)
            expect_false(dat_calib_blr_w_function[["metadata"]]$CI)

            ## Check answer is same whether w_function used or not
            expect_equal(dat_calib_blr[["plotdata"]][[1]], dat_calib_blr_w_function[["plotdata"]][[1]])
            expect_equal(dat_calib_blr[["plotdata"]][[6]], dat_calib_blr_w_function[["plotdata"]][[6]])

            ###
            ### Redefine calc_weights, but change order of all the input arguments (this shouldn't make a difference)
            calc_weights_manual <- function(stabilised = FALSE, max_follow = NULL, data_ms, covs = NULL, landmark_type = "state", j = NULL, t, s, max_weight = 10, data_raw){

              ### Modify everybody to be censored after time t, if a max_follow has been specified
              if(!is.null(max_follow)){
                if (max_follow == "t"){
                  data_raw <- dplyr::mutate(data_raw,
                                            dtcens_s = dplyr::case_when(dtcens < t + 2 ~ dtcens_s,
                                                                        dtcens >= t + 2 ~ 0),
                                            dtcens = dplyr::case_when(dtcens < t + 2 ~ dtcens,
                                                                      dtcens >= t + 2 ~ t + 2))
                } else {
                  data_raw <- dplyr::mutate(data_raw,
                                            dtcens_s = dplyr::case_when(dtcens < max_follow + 2 ~ dtcens_s,
                                                                        dtcens >= max_follow + 2 ~ 0),
                                            dtcens = dplyr::case_when(dtcens < max_follow + 2 ~ dtcens,
                                                                      dtcens >= max_follow + 2 ~ max_follow + 2))
                }
              }

              ### Create a new outcome, which is the time until censored from s
              data_raw$dtcens_modified <- data_raw$dtcens - s

              ### Save a copy of data_raw
              data_raw_save <- data_raw

              ### If landmark_type = "state", calculate weights only in individuals in state j at time s
              ### If landmark_type = "all", calculate weights in all uncensored individuals at time s (note that this excludes individuals
              ### who have reached absorbing states, who have been 'censored' from the survival distribution is censoring)
              if (landmark_type == "state"){
                ### Identify individuals who are uncensored in state j at time s
                ids_uncens <- base::subset(data_ms, from == j & Tstart <= s & s < Tstop) |>
                  dplyr::select(id) |>
                  dplyr::distinct(id) |>
                  dplyr::pull(id)

              } else if (landmark_type == "all"){
                ### Identify individuals who are uncensored time s
                ids_uncens <- base::subset(data_ms, Tstart <= s & s < Tstop) |>
                  dplyr::select(id) |>
                  dplyr::distinct(id) |>
                  dplyr::pull(id)

              }

              ### Subset data_ms and data_raw to these individuals
              data_ms <- data_ms |> base::subset(id %in% ids_uncens)
              data_raw <- data_raw |> base::subset(id %in% ids_uncens)

              ###
              ### Create models for censoring in order to calculate the IPCW weights
              ### Seperate models for estimating the weights, and stabilising the weights (intercept only model)
              ###
              if (!is.null(covs)){
                ### A model where we adjust for predictor variables
                cens_model <- survival::coxph(stats::as.formula(paste("survival::Surv(dtcens_modified, dtcens_s) ~ ",
                                                                      paste(covs, collapse = "+"),
                                                                      sep = "")),
                                              data = data_raw)

                ### Intercept only model (numerator for stabilised weights)
                cens_model_int <- survival::coxph(stats::as.formula(paste("survival::Surv(dtcens_modified, dtcens_s) ~ 1",
                                                                          sep = "")),
                                                  data = data_raw)
              } else if (is.null(covs)){
                ### If user has not input any predictors for estimating weights, the model for estimating the weights is the intercept only model (i_e_ Kaplan Meier estimator)

                ### Intercept only model (numerator for stabilised weights)
                cens_model_int <- survival::coxph(stats::as.formula(paste("survival::Surv(dtcens_modified, dtcens_s) ~ 1",
                                                                          sep = "")),
                                                  data = data_raw)
                ### Assign cens_model to be the same
                cens_model <- cens_model_int


              }

              ### Calculate a data frame containing probability of censored and uncenosred at each time point
              ### The weights will be the probability of being uncensored, at the time of the event for each individual

              ## Extract baseline hazard
              data_weights <- survival::basehaz(cens_model, centered = FALSE)
              ## Add lp to data_raw_save
              data_raw_save$lp <- stats::predict(cens_model, newdata = data_raw_save, type = "lp", reference = "zero")

              ### Create weights for the cohort at time t - s
              ### Note for individuals who reached an absorbing state, we take the probability of them being uncensored at the time of reached the
              ### abosrbing state_ For individuals still alive, we take the probability of being uncensored at time t - s_

              ### Get location of individuals who entered absorbing states or were censored prior to evaluation time
              obs_absorbed_prior <- which(data_raw_save$dtcens <= t & data_raw_save$dtcens_s == 0)
              obs_censored_prior <- which(data_raw_save$dtcens <= t & data_raw_save$dtcens_s == 1)

              ###
              ### Now create unstabilised probability of (un)censoring weights
              ### Note that weights are the probability of being uncensored, so if an individual has low probability of being uncesored,
              ### the inervse of this will be big, weighting them strongly
              ###

              ### First assign all individuals a weight of the probability of being uncensored at time t
              ### This is the linear predictor times the cumulative hazard at time t, and appropriate transformation to get a risk
              data_raw_save$pcw <- as.numeric(exp(-exp(data_raw_save$lp)*data_weights$hazard[max(which(data_weights$time <= t - s))]))

              ## Write a function which will extract the uncensored probability for an individual with linear predictor lp at a given time t
              prob_uncens_func <- function(input){

                ## Assign t and person_id
                t <- input[1]
                lp <- input[2]

                if (t <= 0){
                  return(NA)
                } else if (t > 0){
                  ## Get hazard at appropriate time
                  if (t < min(data_weights$time)){
                    bhaz_t <- 0
                  } else if (t >= min(data_weights$time)){
                    bhaz_t <- data_weights$hazard[max(which(data_weights$time <= t))]
                  }

                  ## Return risk
                  return(exp(-exp(lp)*bhaz_t))
                }
              }

              ### Apply this function to all the times at which individuals have entered an absorbing state prior to censoring
              data_raw_save$pcw[obs_absorbed_prior] <- apply(data_raw_save[obs_absorbed_prior, c("dtcens_modified", "lp")], 1, FUN = prob_uncens_func)

              ### For individuals who were censored prior to t, assign the weight as NA
              data_raw_save$pcw[obs_censored_prior] <- NA

              ### Invert these
              data_raw_save$ipcw <- 1/data_raw_save$pcw

              ###
              ### Stabilise these weights dependent on user-input
              ###
              if (stabilised == TRUE){

                ## Extract baseline hazard
                data_weights_numer <- survival::basehaz(cens_model_int, centered = TRUE)

                ### Assign all individuals a weight of the probability of being uncesored at time t
                data_raw_save$pcw_numer <- as.numeric(exp(-data_weights_numer$hazard[max(which(data_weights_numer$time <= t - s))]))

                ### Create stabilised weight
                data_raw_save$ipcw_stab <- data_raw_save$pcw_numer*data_raw_save$ipcw
              }

              ### Finally cap these at 10 and create output object

              ### Create output object
              if (stabilised == FALSE){
                data_raw_save$ipcw <- pmin(data_raw_save$ipcw, max_weight)
                output_weights <- data.frame("id" = data_raw_save$id, "ipcw" = data_raw_save$ipcw, "pcw" = data_raw_save$pcw)
              } else if (stabilised == TRUE){
                data_raw_save$ipcw <- pmin(data_raw_save$ipcw, max_weight)
                data_raw_save$ipcw_stab <- pmin(data_raw_save$ipcw_stab, max_weight)
                output_weights <- data.frame("id" = data_raw_save$id, "ipcw" = data_raw_save$ipcw, "ipcw_stab" = data_raw_save$ipcw_stab, "pcw" = data_raw_save$pcw)
              }

              return(output_weights)

            }

            ### Calculate observed event probabilities with new w_function
            dat_calib_blr_w_function <-
              calib_msm(data_ms = msebmtcal,
                       data_raw = ebmtcal,
                       j=1,
                       s=0,
                       t = 1826,
                       tp_pred = tp_pred, calib_type = 'blr',
                       curve_type = "rcs",
                       rcs_nk = 3,
                       w_function = calc_weights_manual,
                       w_covs = c("year", "agecl", "proph", "match"))

            expect_type(dat_calib_blr_w_function, "list")
            expect_equal(class(dat_calib_blr_w_function), c("calib_blr", "calib_msm"))
            expect_length(dat_calib_blr_w_function[["mean"]], 6)
            expect_length(dat_calib_blr_w_function[["plotdata"]], 6)
            expect_length(dat_calib_blr_w_function[["plotdata"]][[1]]$id, 1778)
            expect_length(dat_calib_blr_w_function[["plotdata"]][[6]]$id, 1778)
            expect_false(dat_calib_blr_w_function[["metadata"]]$CI)

            ## Check answer is same whether w_function used or not
            expect_equal(dat_calib_blr[["plotdata"]][[1]], dat_calib_blr_w_function[["plotdata"]][[1]])
            expect_equal(dat_calib_blr[["plotdata"]][[6]], dat_calib_blr_w_function[["plotdata"]][[6]])


            ###
            ### Repeat this process (manual definition of calc_weights), again arguments are in different order, but this time an extra argument is added, which adds 10 to every weight_
            ### This extra arguments is something that could be inputted by user, and want to check it does actually change the answer_ It should no longer agree with dat_calb_blr_
            calc_weights_manual <- function(stabilised = FALSE, max_follow = NULL, data_ms, covs = NULL, landmark_type = "state", j = NULL, t, s, max_weight = 10, data_raw, extra_arg = NULL){

              ### Modify everybody to be censored after time t, if a max_follow has been specified
              if(!is.null(max_follow)){
                if (max_follow == "t"){
                  data_raw <- dplyr::mutate(data_raw,
                                            dtcens_s = dplyr::case_when(dtcens < t + 2 ~ dtcens_s,
                                                                        dtcens >= t + 2 ~ 0),
                                            dtcens = dplyr::case_when(dtcens < t + 2 ~ dtcens,
                                                                      dtcens >= t + 2 ~ t + 2))
                } else {
                  data_raw <- dplyr::mutate(data_raw,
                                            dtcens_s = dplyr::case_when(dtcens < max_follow + 2 ~ dtcens_s,
                                                                        dtcens >= max_follow + 2 ~ 0),
                                            dtcens = dplyr::case_when(dtcens < max_follow + 2 ~ dtcens,
                                                                      dtcens >= max_follow + 2 ~ max_follow + 2))
                }
              }

              ### Create a new outcome, which is the time until censored from s
              data_raw$dtcens_modified <- data_raw$dtcens - s

              ### Save a copy of data_raw
              data_raw_save <- data_raw

              ### If landmark_type = "state", calculate weights only in individuals in state j at time s
              ### If landmark_type = "all", calculate weights in all uncensored individuals at time s (note that this excludes individuals
              ### who have reached absorbing states, who have been 'censored' from the survival distribution is censoring)
              if (landmark_type == "state"){
                ### Identify individuals who are uncensored in state j at time s
                ids_uncens <- base::subset(data_ms, from == j & Tstart <= s & s < Tstop) |>
                  dplyr::select(id) |>
                  dplyr::distinct(id) |>
                  dplyr::pull(id)

              } else if (landmark_type == "all"){
                ### Identify individuals who are uncensored time s
                ids_uncens <- base::subset(data_ms, Tstart <= s & s < Tstop) |>
                  dplyr::select(id) |>
                  dplyr::distinct(id) |>
                  dplyr::pull(id)

              }

              ### Subset data_ms and data_raw to these individuals
              data_ms <- data_ms |> base::subset(id %in% ids_uncens)
              data_raw <- data_raw |> base::subset(id %in% ids_uncens)

              ###
              ### Create models for censoring in order to calculate the IPCW weights
              ### Seperate models for estimating the weights, and stabilising the weights (intercept only model)
              ###
              if (!is.null(covs)){
                ### A model where we adjust for predictor variables
                cens_model <- survival::coxph(stats::as.formula(paste("survival::Surv(dtcens_modified, dtcens_s) ~ ",
                                                                      paste(covs, collapse = "+"),
                                                                      sep = "")),
                                              data = data_raw)

                ### Intercept only model (numerator for stabilised weights)
                cens_model_int <- survival::coxph(stats::as.formula(paste("survival::Surv(dtcens_modified, dtcens_s) ~ 1",
                                                                          sep = "")),
                                                  data = data_raw)
              } else if (is.null(covs)){
                ### If user has not input any predictors for estimating weights, the model for estimating the weights is the intercept only model (i_e_ Kaplan Meier estimator)

                ### Intercept only model (numerator for stabilised weights)
                cens_model_int <- survival::coxph(stats::as.formula(paste("survival::Surv(dtcens_modified, dtcens_s) ~ 1",
                                                                          sep = "")),
                                                  data = data_raw)
                ### Assign cens_model to be the same
                cens_model <- cens_model_int


              }

              ### Calculate a data frame containing probability of censored and uncenosred at each time point
              ### The weights will be the probability of being uncensored, at the time of the event for each individual

              ## Extract baseline hazard
              data_weights <- survival::basehaz(cens_model, centered = FALSE)
              ## Add lp to data_raw_save
              data_raw_save$lp <- stats::predict(cens_model, newdata = data_raw_save, type = "lp", reference = "zero")

              ### Create weights for the cohort at time t - s
              ### Note for individuals who reached an absorbing state, we take the probability of them being uncensored at the time of reached the
              ### abosrbing state_ For individuals still alive, we take the probability of being uncensored at time t - s_

              ### Get location of individuals who entered absorbing states or were censored prior to evaluation time
              obs_absorbed_prior <- which(data_raw_save$dtcens <= t & data_raw_save$dtcens_s == 0)
              obs_censored_prior <- which(data_raw_save$dtcens <= t & data_raw_save$dtcens_s == 1)

              ###
              ### Now create unstabilised probability of (un)censoring weights
              ### Note that weights are the probability of being uncensored, so if an individual has low probability of being uncesored,
              ### the inervse of this will be big, weighting them strongly
              ###

              ### First assign all individuals a weight of the probability of being uncensored at time t
              ### This is the linear predictor times the cumulative hazard at time t, and appropriate transformation to get a risk
              data_raw_save$pcw <- as.numeric(exp(-exp(data_raw_save$lp)*data_weights$hazard[max(which(data_weights$time <= t - s))]))

              ## Write a function which will extract the uncensored probability for an individual with linear predictor lp at a given time t
              prob_uncens_func <- function(input){

                ## Assign t and person_id
                t <- input[1]
                lp <- input[2]

                if (t <= 0){
                  return(NA)
                } else if (t > 0){
                  ## Get hazard at appropriate time
                  if (t < min(data_weights$time)){
                    bhaz_t <- 0
                  } else if (t >= min(data_weights$time)){
                    bhaz_t <- data_weights$hazard[max(which(data_weights$time <= t))]
                  }

                  ## Return risk
                  return(exp(-exp(lp)*bhaz_t))
                }
              }

              ### Apply this function to all the times at which individuals have entered an absorbing state prior to censoring
              data_raw_save$pcw[obs_absorbed_prior] <- apply(data_raw_save[obs_absorbed_prior, c("dtcens_modified", "lp")], 1, FUN = prob_uncens_func)

              ### For individuals who were censored prior to t, assign the weight as NA
              data_raw_save$pcw[obs_censored_prior] <- NA

              ### Invert these
              data_raw_save$ipcw <- 1/data_raw_save$pcw

              ###
              ### Stabilise these weights dependent on user-input
              ###
              if (stabilised == TRUE){

                ## Extract baseline hazard
                data_weights_numer <- survival::basehaz(cens_model_int, centered = TRUE)

                ### Assign all individuals a weight of the probability of being uncesored at time t
                data_raw_save$pcw_numer <- as.numeric(exp(-data_weights_numer$hazard[max(which(data_weights_numer$time <= t - s))]))

                ### Create stabilised weight
                data_raw_save$ipcw_stab <- data_raw_save$pcw_numer*data_raw_save$ipcw
              }

              ### Finally cap these at 10 and create output object

              ### Create output object
              if (stabilised == FALSE){
                data_raw_save$ipcw <- pmin(data_raw_save$ipcw, max_weight)
                output_weights <- data.frame("id" = data_raw_save$id, "ipcw" = data_raw_save$ipcw, "pcw" = data_raw_save$pcw)
              } else if (stabilised == TRUE){
                data_raw_save$ipcw <- pmin(data_raw_save$ipcw, max_weight)
                data_raw_save$ipcw_stab <- pmin(data_raw_save$ipcw_stab, max_weight)
                output_weights <- data.frame("id" = data_raw_save$id, "ipcw" = data_raw_save$ipcw, "ipcw_stab" = data_raw_save$ipcw_stab, "pcw" = data_raw_save$pcw)
              }

              ### Add this extra argument to the weights, to check it does something
              output_weights$ipcw <- output_weights$ipcw + extra_arg

              return(output_weights)

            }

            ### Calculate observed event probabilities with new w_function
            dat_calib_blr_w_function <-
              calib_msm(data_ms = msebmtcal,
                       data_raw = ebmtcal,
                       j=1,
                       s=0,
                       t = 1826,
                       tp_pred = tp_pred, calib_type = 'blr',
                       curve_type = "rcs",
                       rcs_nk = 3,
                       w_function = calc_weights_manual,
                       w_covs = c("year", "agecl", "proph", "match"),
                       extra_arg = 10)

            expect_type(dat_calib_blr_w_function, "list")
            expect_equal(class(dat_calib_blr_w_function), c("calib_blr", "calib_msm"))
            expect_length(dat_calib_blr_w_function[["mean"]], 6)
            expect_length(dat_calib_blr_w_function[["plotdata"]], 6)
            expect_length(dat_calib_blr_w_function[["plotdata"]][[1]]$id, 1778)
            expect_length(dat_calib_blr_w_function[["plotdata"]][[6]]$id, 1778)
            expect_false(dat_calib_blr_w_function[["metadata"]]$CI)

            ## Check answer is same whether w_function used or not
            expect_false(any(dat_calib_blr[["plotdata"]][[1]]$obs == dat_calib_blr_w_function[["plotdata"]][[1]]$obs))

          })

test_that("test warnings and errors", {

  ## Extract relevant predicted risks from tps0
  tp_pred <- dplyr::select(dplyr::filter(tps0, j == 3), any_of(paste("pstate", 1:6, sep = "")))

  ## Calculate observed event probabilities
  expect_error(
    calib_msm(data_ms = msebmtcal,
             data_raw = ebmtcal,
             j=1,
             s=0,
             t = 1826,
             tp_pred = tp_pred, calib_type = 'blr',
             curve_type = "rcs",
             rcs_nk = 3,
             w_covs = c("year", "agecl", "proph", "match"),
             transitions_out = c(1,2,3,4))
  )

  ## Calculate observed event probabilities
  weights_manual <-
    calc_weights(data_ms = msebmtcal,
                 data_raw = ebmtcal,
                 t = 1826,
                 s = 0,
                 landmark_type = "state",
                 j = 1,
                 max_weight = 10,
                 stabilised = FALSE)$ipcw[-1]

  expect_error(
    calib_msm(data_ms = msebmtcal,
             data_raw = ebmtcal,
             j=1,
             s=0,
             t = 1826,
             tp_pred = tp_pred, calib_type = 'blr',
             curve_type = "rcs",
             rcs_nk = 3,
             weights = weights_manual)
  )

  ## Write a weights function with the wrong variable names
  calc_weights_error <- function(data_ms, data_raw, covs = NULL, t, s, j = NULL, max_weight = 10, stabilised = FALSE, max_follow = NULL){
    return(data_ms)
  }
  expect_error(
    calib_msm(data_ms = msebmtcal,
             data_raw = ebmtcal,
             j=1,
             s=0,
             t = 1826,
             tp_pred = tp_pred, calib_type = 'blr',
             curve_type = "rcs",
             rcs_nk = 3,
             w_function = calc_weights_error,
             w_covs = c("year", "agecl", "proph", "match"))
  )

  ### check  warnings when there are zero predicted probabilities for valid transitions

  ## Extract relevant predicted risks from tps0
  tp_pred <- dplyr::select(dplyr::filter(tps100, j == 1), dplyr::any_of(paste("pstate", 1:6, sep = "")))
  tp_pred[,1] <- rep(0, nrow(tp_pred))

  ###
  ### Calculate observed event probabilities
  expect_error(calib_msm(data_ms = msebmtcal,
                        data_raw = ebmtcal,
                        j=1,
                        s=100,
                        t = 1826,
                        tp_pred = tp_pred, calib_type = 'blr',
                        curve_type = "loess",
                        rcs_nk = 3,
                        w_covs = c("year", "agecl", "proph", "match")))


  ### check error when there are non-zero predicted probabilities for transitions which do not happen

  ## Extract relevant predicted risks from tps0
  tp_pred <- dplyr::select(dplyr::filter(tps100, j == 1), dplyr::any_of(paste("pstate", 1:6, sep = "")))
  tp_pred[,3] <- runif(nrow(tp_pred), 0, 1)

  ###
  ### Calculate observed event probabilities
  expect_error(calib_msm(data_ms = msebmtcal,
                        data_raw = ebmtcal,
                        j=1,
                        s=100,
                        t = 1826,
                        tp_pred = tp_pred, calib_type = 'blr',
                        curve_type = "loess",
                        rcs_nk = 3,
                        w_covs = c("year", "agecl", "proph", "match")))

  ## Extract relevant predicted risks from tps0
  tp_pred <- dplyr::select(dplyr::filter(tps100, j == 3), dplyr::any_of(paste("pstate", 1:6, sep = "")))
  tp_pred[,1] <- runif(nrow(tp_pred), 0, 1)

  ###
  ### Calculate observed event probabilities
  expect_error(calib_msm(data_ms = msebmtcal,
                        data_raw = ebmtcal,
                        j=3,
                        s=100,
                        t = 1826,
                        tp_pred = tp_pred, calib_type = 'blr',
                        curve_type = "loess",
                        rcs_nk = 3,
                        w_covs = c("year", "agecl", "proph", "match")))

  ### check error when there are zero predicted probabilities for transitions which do happen in dataset

  ## Extract relevant predicted risks from tps0
  tp_pred <- dplyr::select(dplyr::filter(tps100, j == 3), dplyr::any_of(paste("pstate", 1:6, sep = "")))
  tp_pred[,3] <- rep(0, nrow(tp_pred))

  ###
  ### Calculate observed event probabilities
  expect_error(calib_msm(data_ms = msebmtcal,
                        data_raw = ebmtcal,
                        j=3,
                        s=100,
                        t = 1826,
                        tp_pred = tp_pred, calib_type = 'blr',
                        curve_type = "loess",
                        rcs_nk = 3,
                        w_covs = c("year", "agecl", "proph", "match")))


  ## Extract relevant predicted risks from tps0
  tp_pred <- dplyr::select(dplyr::filter(tps100, j == 3), dplyr::any_of(paste("pstate", 1:6, sep = "")))
  tp_pred[,1] <- rep(0, nrow(tp_pred))

})
