### Section 1

# install.packages("calibmsm")
# devtools::install_github("alexpate30/calibmsm")

### Section 4.1
library(calibmsm)

data("ebmtcal")
head(ebmtcal, n = 3)

data("msebmtcal")
head(msebmtcal, n = 10)

data("tps0")
head(tps0, n = 3)

data("tps100")
head(tps100, n = 3)

### Section 4.2
tp_pred_s0 <- tps0 |>
  dplyr::filter(j == 1) |>
  dplyr::select(any_of(paste("pstate", 1:6, sep = "")))

t_eval <- 1826

dat_calib_blr <-
  calib_msm(data_ms = msebmtcal,
            data_raw = ebmtcal,
            j=1,
            s=0,
            t = t_eval,
            tp_pred = tp_pred_s0,
            calib_type = "blr",
            curve_type = "rcs",
            rcs_nk = 3,
            w_covs = c("year", "agecl", "proph", "match"),
            CI = 95,
            CI_R_boot = 200)

print(dat_calib_blr)

metadata(dat_calib_blr)

plot_blr <- plot(dat_calib_blr, combine = TRUE, nrow = 2, ncol = 3)
grid::grid.draw(plot_blr)

dat_calib_pv <-
  calib_msm(data_ms = msebmtcal,
            data_raw = ebmtcal,
            j=1,
            s=0,
            t = t_eval,
            tp_pred = tp_pred_s0,
            calib_type = "pv",
            curve_type = "rcs",
            rcs_nk = 3,
            pv_group_vars = c("year"),
            pv_n_pctls = 3,
            CI = 95,
            CI_type = "parametric")

plot_pv <- plot(dat_calib_pv, combine = TRUE, nrow = 2, ncol = 3)
grid::grid.draw(plot_pv)

dat_calib_mlr <-
  calib_msm(data_ms = msebmtcal,
            data_raw = ebmtcal,
            j=1,
            s=0,
            t = 1826,
            tp_pred = tp_pred_s0,
            calib_type = "mlr",
            w_covs = c("year", "agecl", "proph", "match"))

plot_mlr <- plot(dat_calib_mlr, combine = TRUE, nrow = 2, ncol = 3)
grid::grid.draw(plot_mlr)

### Section 4.3
tp_pred_j1s100 <- tps100 |>
  dplyr::filter(j == 1) |>
  dplyr::select(any_of(paste("pstate", 1:6, sep = "")))

tp_pred_j3s100 <- tps100 |>
  dplyr::filter(j == 3) |>
  dplyr::select(any_of(paste("pstate", 1:6, sep = "")))

dat_calib_pv_j1_s100 <-
  calib_msm(data_ms = msebmtcal,
            data_raw = ebmtcal,
            j=1,
            s=100,
            t = t_eval,
            tp_pred = tp_pred_j1s100,
            calib_type = "pv",
            curve_type = "rcs",
            rcs_nk = 3,
            pv_group_vars = c("year"),
            CI = 95,
            CI_type = "parametric")

plot_pv_j1_s100 <- plot(dat_calib_pv_j1_s100, combine = TRUE, nrow = 2, ncol = 3)
grid::grid.draw(plot_pv_j1_s100)

dat_calib_pv_j3_s100 <-
  calib_msm(data_ms = msebmtcal,
            data_raw = ebmtcal,
            j=3,
            s=100,
            t = t_eval,
            tp_pred = tp_pred_j3s100,
            calib_type = "pv",
            curve_type = "rcs",
            rcs_nk = 3,
            pv_group_vars = c("year"),
            CI = 95,
            CI_type = "parametric")

plot_pv_j3_s100 <- plot(dat_calib_pv_j3_s100, combine = TRUE, nrow = 2, ncol = 3)
grid::grid.draw(plot_pv_j3_s100)

