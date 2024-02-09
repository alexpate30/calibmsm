# rm(list=ls())
#   rm(list = setdiff(ls(), list("data.raw.lmk.js.normal", "data.pctls.normal")))
# devtools::load_all()
#   data("ebmtcal")
#   data("msebmtcal")
#   data("tps0")
#   data("tps100")
# calib.type <- "aj"
# data.raw <- ebmtcal
# data.mstate <- msebmtcal
# tp.pred <- tps0 |>
#   subset(j == 1) |>
#   dplyr::select(paste("pstate", 1:6, sep = ""))
#
#
# j <- 1
# s <- 0
# t <- 1826
#
# curve.type = "rcs"
# tp.pred.plot = NULL
# transitions.out = NULL
# weights = NULL
#
# w.covs = NULL
# w.landmark.type = "state"
# w.max = 10
# w.stabilised = FALSE
# w.max.follow = NULL
# w.function = NULL
#
# CI = FALSE
# # CI = 95
# CI.R.boot = 2
# rcs.nk = 3
# CI.type = "bootstrap"
#
# CI.seed = 1
# loess.span = 1
# loess.degree = 1
#
# pv.group.vars = c("year")
# pv.n.pctls = 2
#
# mlr.smoother.type = "sm.ps"
# mlr.ps.int = 4
# mlr.degree = 3
# mlr.s.df = 4
# mlr.niknots = 4
# assess.moderate = TRUE
# assess.mean = TRUE
# #
# pv.precalc = NULL
# str(data.raw)
#
# pv.group.vars = c("x12")
# tp.pred <- readRDS("P:/Documents/aaa_incline/tp.pred.reduc.rds")
# data.raw <- readRDS("P:/Documents/aaa_incline/data.raw.reduc.rds")
# data.mstate <- readRDS("P:/Documents/aaa_incline/data.mstate.reduc.rds")
# t <- 2557

# str(ebmtcal)
# str(readRDS("P:/Documents/aaa_incline/data.raw.reduc.rds"))
# calib.type <- "blr"
#
# # ## Calculate manual weights
# # weights.manual <-
# #   calc_weights(data.mstate = msebmtcal,
# #                data.raw = ebmtcal,
# #                t = 1826,
# #                s = 0,
# #                landmark.type = "state",
# #                j = 1,
# #                max.weight = 10,
# #                stabilised = FALSE)
# # weights <- weights.manual$ipcw
