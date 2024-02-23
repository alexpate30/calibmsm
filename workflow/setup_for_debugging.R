rm(list=ls())

devtools::load_all()
  data("ebmtcal")
  data("msebmtcal")
  data("tps0")
  data("tps100")

  head(data.raw)
data.raw <- ebmtcal[rep(1:nrow(ebmtcal), each = 100), ]
str(data.raw)
data.mstate <- msebmtcal
tp.pred <- tps0 |>
  subset(j == 1) |>
  dplyr::select(paste("pstate", 1:6, sep = ""))

pv.comb <- tp.pred

dat.calib.pv.noCI <-
  calib_msm(data.mstate = msebmtcal,
            data.raw = ebmtcal,
            j=1,
            s=0,
            t = 1826,
            tp.pred = tp.pred,
            calib.type = "pv",
            curve.type = "loess",
            assess.moderate = TRUE,
            assess.mean = FALSE,
            pv.precalc = pv.comb)
str(dat.calib.pv.noCI)


dat.calib.pv.CI <-
  calib_msm(data.mstate = msebmtcal,
            data.raw = ebmtcal,
            j=1,
            s=0,
            t = 1826,
            tp.pred = tp.pred,
            calib.type = "pv",
            curve.type = "loess",
            CI = 95,
            CI.type = "parametric",
            assess.moderate = TRUE,
            assess.mean = FALSE,
            pv.precalc = pv.comb)
str(dat.calib.pv.CI)

str(split(1:10, c(1,1,1,1,1,2,2,2,2,2)))
rep(1:ceiling(nrow(data.raw)/10000), each = 10000, length.out = nrow(data.raw))

j <- 1
s <- 0
t <- 1826

curve.type = "rcs"
tp.pred.plot = NULL
transitions.out = NULL
weights = NULL

w.covs = NULL
w.landmark.type = "state"
w.max = 10
w.stabilised = FALSE
w.max.follow = NULL
w.function = NULL

CI = FALSE
# CI = 95
CI.R.boot = 2
rcs.nk = 3
CI.type = "bootstrap"

CI.seed = 1
loess.span = 1
loess.degree = 1

pv.group.vars = c("year")
pv.n.pctls = 2

mlr.smoother.type = "sm.ps"
mlr.ps.int = 4
mlr.degree = 3
mlr.s.df = 4
mlr.niknots = 4
assess.moderate = TRUE
assess.mean = TRUE

pv.precalc = NULL
str(data.raw)

pv.group.vars = c("x12")
tp.pred <- readRDS("P:/Documents/aaa_incline/tp.pred.reduc.rds")
data.raw <- readRDS("P:/Documents/aaa_incline/data.raw.reduc.rds")
data.mstate <- readRDS("P:/Documents/aaa_incline/data.mstate.reduc.rds")
t <- 2557

test <- calib_msm(data.mstate = data.mstate)
