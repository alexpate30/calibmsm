### This contains a quick way to manually define all the input the a function, so I can more easily debug them line by line

###
### calib_blr
###
# load_all()
# data.mstate <- msebmtcal
# data.raw <- ebmtcal
# j <- 1
# j.in <- 1
# s<-0
# t <- 1826
# tp.pred = tps100 |> dplyr::filter(j == 1) |> dplyr::select(any_of(paste("pstate", 1:6, sep = "")))
# curve.type = "rcs"
# rcs.nk = 3
# weights <- NULL
# # weights <- weights.manual
# w.covs = c("year", "agecl", "proph", "match")
# w.landmark.type = "state"
# w.max = 10
# w.stabilised = FALSE
# w.max.follow = NULL
# CI = 95
# CI.type = "parametric"
# str(data.pred.plot)
# data.pred.plot <- tps0 |> dplyr::filter(j == j.in) |> dplyr::select(any_of(paste("pstate", 1:6, sep = "")))
# data.pred.plot$id <- 1:nrow(data.pred.plot)

###
### calib_mlr
###
# data.mstate <- msebmtcal
# data.raw <- ebmtcal
# j<-1
# s<-0
# t <- 1826
# tp.pred = tps0 |> dplyr::filter(j == 1) |> dplyr::select(any_of(paste("pstate", 1:6, sep = "")))
# ps.int <- 4
# degree <- 3
# weights <- NULL
# w.covs = c("year", "agecl", "proph", "match")
# w.landmark.type = "all"
# w.max = 10
# w.stabilised = FALSE
# smoother.type <- "sm.ps"


###
### calib_pv
###
# data.mstate <- msebmtcal
# data.raw <- ebmtcal
# #
# # indices <- sample(1:nrow(data.raw), nrow(data.raw), replace = TRUE)
# # data.raw.boot <- data.raw[indices, ]
# # data.mstate.boot <-
# #   do.call("rbind", lapply(
# #     data.raw.boot$id, function(x) {base::subset(data.mstate, id == x)})
# #   )
# # attributes(data.mstate.boot)$trans <- attributes(data.mstate)$trans
# #
# j <- 3
# s <- 100
# t <- 1826
# tp.pred <- tps100 |> dplyr::filter(j == 3) |> dplyr::select(any_of(paste("pstate", 1:6, sep = "")))
# # id.lmk <- extract_ids_states(data.mstate = data.mstate,
# #                              tmat = attributes(data.mstate)$trans,
# #                              j = j,
# #                              t = s)
# # data.pred.plot <- tps100 |>
# #   dplyr::filter(id %in% id.lmk) |>
# #   dplyr::filter(j == 3) |>
# #   dplyr::select(any_of(paste("pstate", 1:6, sep = "")))
# data.pred.plot <- NULL
# transitions.out <- NULL
# group.vars <- "year"
# n.pctls <- 2
# loess.span <- 0.75
# loess.degree <- 2
# CI <- FALSE
# CI.type <- parametric
# CI.R.boot <- 2
# rcs.nk <- 3
# curve.type <- "rcs"

#
# data.mstate <- data.mstate.boot
# data.raw <- data.raw.boot

###
###
###
