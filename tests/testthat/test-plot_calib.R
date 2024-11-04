test_that("check plot_calib_msm output (j = 1, s = 0)", {

  ## Extract relevant predicted risks from tps0
  tp_pred <- dplyr::select(dplyr::filter(tps0, j == 1), any_of(paste("pstate", 1:6, sep = "")))

  ## Calculate observed event probabilities
  dat_calib_blr <-
    calib_msm(data_ms = msebmtcal,
             data_raw = ebmtcal,
             j=1,
             s=0,
             t = 1826,
             tp_pred = tp_pred,
             calib_type = "blr",
             curve_type = "rcs",
             rcs_nk = 3,
             w_covs = c("year", "agecl", "proph", "match"))

  ## Plot calibration plots and run tests
  plot_object <- plot(dat_calib_blr, combine = TRUE, nrow = 2, ncol = 3)
  expect_equal(class(plot_object), c("gtable", "gTree", "grob", "gDesc"))
  plot_object <- plot(dat_calib_blr, combine = FALSE, nrow = 2, ncol = 3)
  expect_length(plot_object, 6)
  expect_type(plot_object, "list")

  ## Plot calibration plots and run tests with marginal density plots
  plot_object <- plot(dat_calib_blr, combine = TRUE, nrow = 2, ncol = 3, marg_density = FALSE, marg_density_size = 1)
  expect_length(plot_object, 11)
  expect_equal(class(plot_object), c("gg", "ggplot", "ggarrange"))

  ## Plot calibration plots and run tests with marginal rug plots
  plot_object <- plot(dat_calib_blr, combine = TRUE, nrow = 2, ncol = 3, marg_rug = TRUE)
  expect_equal(class(plot_object), c("gtable", "gTree", "grob", "gDesc"))

  ## Add titles
  plot_object <- plot(dat_calib_blr, combine = TRUE, nrow = 2, ncol = 3, marg_rug = TRUE,
                      titles = paste("eggs", 1:6),
                      axis_titles_text_x = paste("eggs_x", 1:6),
                      axis_titles_text_y = paste("eggs_y", 1:6))
  expect_equal(class(plot_object), c("gtable", "gTree", "grob", "gDesc"))

})

test_that("check plot_calib_msm output (j = 1, s = 0) with CI", {

  ## Extract relevant predicted risks from tps0
  tp_pred <- dplyr::select(dplyr::filter(tps0, j == 1), any_of(paste("pstate", 1:6, sep = "")))

  ## Calculate observed event probabilities
  dat_calib_blr <-
    calib_msm(data_ms = msebmtcal,
             data_raw = ebmtcal,
             j=1,
             s=0,
             t = 1826,
             tp_pred = tp_pred,
             calib_type = "blr",
             curve_type = "rcs",
             rcs_nk = 3,
             w_covs = c("year", "agecl", "proph", "match"),
             CI = 95,
             CI_R_boot = 5)

  ## Plot calibration plots and run tests without marginal density plots
  plot_object <- plot(dat_calib_blr, combine = TRUE, nrow = 2, ncol = 3)
  expect_equal(class(plot_object), c("gtable", "gTree", "grob", "gDesc"))
  plot_object <- plot(dat_calib_blr, combine = FALSE, nrow = 2, ncol = 3)
  expect_length(plot_object, 6)
  expect_type(plot_object, "list")

  ## Plot calibration plots and run tests with marginal density plots
  plot_object <- plot(dat_calib_blr, combine = TRUE, nrow = 2, ncol = 3, marg_density = FALSE, marg_density_size = 1)
  expect_equal(class(plot_object), c("gg", "ggplot", "ggarrange"))

  ## Plot calibration plots and run tests with marginal rug plots
  plot_object <- plot(dat_calib_blr, combine = TRUE, nrow = 2, ncol = 3, marg_rug = TRUE)
  expect_equal(class(plot_object), c("gtable", "gTree", "grob", "gDesc"))

})


test_that("check plot_calib_msm output (j = 3, s = 100)", {

  ## Extract relevant predicted risks from tps0
  tp_pred <- dplyr::select(dplyr::filter(tps100, j == 3), any_of(paste("pstate", 1:6, sep = "")))

  ## Calculate observed event probabilities
  dat_calib_blr <-
    calib_msm(data_ms = msebmtcal,
             data_raw = ebmtcal,
             j=3,
             s=100,
             t = 1826,
             tp_pred = tp_pred,
             calib_type = "blr",
             curve_type = "rcs",
             rcs_nk = 3,
             w_covs = c("year", "agecl", "proph", "match"))

  ## Plot calibration plots and run tests
  plot_object <- plot(dat_calib_blr, combine = TRUE, nrow = 2, ncol = 3)
  expect_equal(class(plot_object), c("gtable", "gTree", "grob", "gDesc"))
  plot_object <- plot(dat_calib_blr, combine = FALSE, nrow = 2, ncol = 3)
  expect_length(plot_object, 4)
  expect_type(plot_object, "list")

})


test_that("check plot_calib_pv output (j = 1, s = 0)", {

  ## Reduce to 50 individuals
  # Extract the predicted transition probabilities out of state j = 1 for first 50 individuals
  tp_pred <- tps0 |>
    dplyr::filter(id %in% 1:50) |>
    dplyr::filter(j == 1) |>
    dplyr::select(any_of(paste("pstate", 1:6, sep = "")))
  # Reduce ebmtcal to first 50 individuals
  ebmtcal <- ebmtcal |> dplyr::filter(id %in% 1:50)
  # Reduce msebmtcal_cmprsk to first 100 individuals
  msebmtcal <- msebmtcal |> dplyr::filter(id %in% 1:50)

  ## Calculate observed event probabilities
  dat_calib_pv <-
    suppressWarnings(calib_msm(data_ms = msebmtcal,
                              data_raw = ebmtcal,
                              j=1,
                              s=0,
                              t = 1826,
                              tp_pred = tp_pred,
                              calib_type = "pv",
                              curve_type = "rcs",
                              rcs_nk = 3))

  ## Plot calibration plots and run tests
  plot_object <- plot(dat_calib_pv, combine = TRUE)
  expect_equal(class(plot_object), c("gtable", "gTree", "grob", "gDesc"))
  plot_object <- plot(dat_calib_pv, combine = FALSE)
  expect_length(plot_object, 6)
  expect_type(plot_object, "list")

})

test_that("check plot_calib_pv output (j = 1, s = 0) with CI", {

  ## Reduce to 50 individuals
  # Extract the predicted transition probabilities out of state j = 1 for first 50 individuals
  tp_pred <- tps0 |>
    dplyr::filter(id %in% 1:50) |>
    dplyr::filter(j == 1) |>
    dplyr::select(any_of(paste("pstate", 1:6, sep = "")))
  # Reduce ebmtcal to first 50 individuals
  ebmtcal <- ebmtcal |> dplyr::filter(id %in% 1:50)
  # Reduce msebmtcal_cmprsk to first 100 individuals
  msebmtcal <- msebmtcal |> dplyr::filter(id %in% 1:50)

  ## Calculate observed event probabilities
  dat_calib_pv <-
    suppressWarnings(calib_msm(data_ms = msebmtcal,
                              data_raw = ebmtcal,
                              j=1,
                              s=0,
                              t = 1826,
                              tp_pred = tp_pred,
                              calib_type = "pv",
                              curve_type = "rcs",
                              rcs_nk = 3,
                              CI = 95,
                              CI_type = "parametric"))

  ## Plot calibration plots and run tests
  plot_object <- plot(dat_calib_pv, combine = TRUE)
  expect_equal(class(plot_object), c("gtable", "gTree", "grob", "gDesc"))
  plot_object <- plot(dat_calib_pv, combine = FALSE)
  expect_length(plot_object, 6)
  expect_type(plot_object, "list")

})


test_that("check plot_calib_pv output (j = 3, s = 100) with CI", {

  ## Reduce to 500 individuals
  # Extract the predicted transition probabilities out of state j = 1 for first 500 individuals
  tp_pred <- tps100 |>
    dplyr::filter(id %in% 1:500) |>
    dplyr::filter(j == 3) |>
    dplyr::select(any_of(paste("pstate", 1:6, sep = "")))
  # Reduce ebmtcal to first 500 individuals
  ebmtcal <- ebmtcal |> dplyr::filter(id %in% 1:500)
  # Reduce msebmtcal_cmprsk to first 100 individuals
  msebmtcal <- msebmtcal |> dplyr::filter(id %in% 1:500)

  ## Calculate observed event probabilities
  dat_calib_pv <-
    calib_msm(data_ms = msebmtcal,
             data_raw = ebmtcal,
             j=3,
             s=100,
             t = 1826,
             tp_pred = tp_pred,
             calib_type = "pv",
             curve_type = "rcs",
             rcs_nk = 3,
             CI = 95,
             CI_type = "parametric")

  ## Plot calibration plots and run tests
  plot_object <- plot(dat_calib_pv, combine = TRUE)
  expect_equal(class(plot_object), c("gtable", "gTree", "grob", "gDesc"))
  plot_object <- plot(dat_calib_pv, combine = FALSE)
  expect_length(plot_object, 4)
  expect_type(plot_object, "list")

})



test_that("check plot_calib_mlr output (j = 1, s = 0)", {

  ## Reduce to 500 individuals
  # Extract the predicted transition probabilities out of state j = 1 for first 500 individuals
  tp_pred <- tps0 |>
    dplyr::filter(id %in% 1:500) |>
    dplyr::filter(j == 1) |>
    dplyr::select(any_of(paste("pstate", 1:6, sep = "")))
  # Reduce ebmtcal to first 500 individuals
  ebmtcal <- ebmtcal |> dplyr::filter(id %in% 1:500)
  # Reduce msebmtcal_cmprsk to first 100 individuals
  msebmtcal <- msebmtcal |> dplyr::filter(id %in% 1:500)

  # ## Extract relevant predicted risks from tps0
  # tp_pred <- dplyr::select(dplyr::filter(tps100, j == 3), dplyr::any_of(paste("pstate", 1:6, sep = "")))

  ## Calculate observed event probabilities
  suppressWarnings(
    dat_calib_mlr <-
      calib_msm(data_ms = msebmtcal,
               data_raw = ebmtcal,
               j=1,
               s=0,
               t = 1826,
               tp_pred = tp_pred,
               calib_type = "mlr",
               w_covs = c("year", "agecl", "proph", "match"))
  )

  ## Plot calibration plots and run tests
  plot_object <- plot(dat_calib_mlr, combine = TRUE, nrow = 2, ncol = 3)
  expect_equal(class(plot_object),  c("gg", "ggplot", "ggarrange"))
  plot_object <- plot(dat_calib_mlr, combine = FALSE, nrow = 2, ncol = 3)
  expect_length(plot_object, 6)
  expect_type(plot_object, "list")

})


test_that("check plot_calib_mlr output (j = 3, s = 100)", {

  ## Extract relevant predicted risks from tps0
  tp_pred <- dplyr::select(dplyr::filter(tps100, j == 3), dplyr::any_of(paste("pstate", 1:6, sep = "")))

  ## Calculate observed event probabilities
  suppressWarnings(
    dat_calib_mlr <-
      calib_msm(data_ms = msebmtcal,
               data_raw = ebmtcal,
               j=3,
               s=100,
               t = 1826,
               tp_pred = tp_pred,
               calib_type = "mlr",
               w_covs = c("year", "agecl", "proph", "match"))
  )

  ## Plot calibration plots and run tests
  plot_object <- plot(dat_calib_mlr, combine = TRUE, nrow = 2, ncol = 3)
  expect_equal(class(plot_object),  c("gg", "ggplot", "ggarrange"))
  plot_object <- plot(dat_calib_mlr, combine = FALSE, nrow = 2, ncol = 3)
  expect_length(plot_object, 4)
  expect_type(plot_object, "list")

})
