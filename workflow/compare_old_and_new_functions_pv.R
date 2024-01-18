### Here we will compare the new functions to old functions to ensure results are the same
rm(list=ls())

devtools::load_all()
data("ebmtcal")
data("msebmtcal")
data("tps0")
data("tps100")

data.raw <- ebmtcal
data.mstate <- msebmtcal
tp.pred.s0 <- tps0 |>
  subset(j == 1) |>
  dplyr::select(paste("pstate", 1:6, sep = ""))
tp.pred.s100 <- tps100 |>
  subset(j == 1) |>
  dplyr::select(paste("pstate", 1:6, sep = ""))
j <- 1
s <- 0
t.eval <- 1826


################
### PV NO CI ###
### No group.vars ###
### No n.pctls ######
#####################

###
#### RCS
###
dat.calib.pv.old <-
  calib_pv_SAVE(data.mstate = msebmtcal,
                data.raw = ebmtcal,
                j=1,
                s=100,
                t = t.eval,
                tp.pred = tp.pred.s100)


dat.calib.pv.new <-
  calibmsm(data.mstate = msebmtcal,
            data.raw = ebmtcal,
            j=1,
            s=100,
            t = t.eval,
            tp.pred = tp.pred.s100,
            calib.type = "pv",
            curve.type = "rcs")

testthat::expect_equal(dat.calib.pv.old[["plotdata"]][[2]], dat.calib.pv.new[["plotdata"]][[2]])

head(dat.calib.pv.old[["plotdata"]][[2]])
head(dat.calib.pv.new[["plotdata"]][[2]])


###
#### LOESS
###
dat.calib.pv.old <-
  calib_pv_SAVE(data.mstate = msebmtcal,
                data.raw = ebmtcal,
                j=1,
                s=100,
                t = t.eval,
                tp.pred = tp.pred.s100,
                curve.type = "loess")


dat.calib.pv.new <-
  calibmsm(data.mstate = msebmtcal,
            data.raw = ebmtcal,
            j=1,
            s=100,
            t = t.eval,
            tp.pred = tp.pred.s100,
            calib.type = "pv",
            curve.type = "loess")

head(dat.calib.pv.old[["plotdata"]][[2]])
head(dat.calib.pv.new[["plotdata"]][[2]])

testthat::expect_equal(dat.calib.pv.old[["plotdata"]][[2]], dat.calib.pv.new[["plotdata"]][[2]])


################
### PV NO CI ###
### With group.vars ###
### No n.pctls ######
#####################

###
#### RCS
###
dat.calib.pv.old <-
  calib_pv_SAVE(data.mstate = msebmtcal,
                data.raw = ebmtcal,
                j=1,
                s=100,
                t = t.eval,
                tp.pred = tp.pred.s100,
                group.vars = c("year"))


dat.calib.pv.new <-
  calibmsm(data.mstate = msebmtcal,
            data.raw = ebmtcal,
            j=1,
            s=100,
            t = t.eval,
            tp.pred = tp.pred.s100,
            calib.type = "pv",
            curve.type = "rcs",
            pv.group.vars = c("year"))

head(dat.calib.pv.old[["plotdata"]][[2]])
head(dat.calib.pv.new[["plotdata"]][[2]])

testthat::expect_equal(dat.calib.pv.old[["plotdata"]][[2]], dat.calib.pv.new[["plotdata"]][[2]])

################
### PV NO CI ###
### No group.vars ###
### With n.pctls ######
#####################

###
#### RCS
###
dat.calib.pv.old <-
  calib_pv_SAVE(data.mstate = msebmtcal,
                data.raw = ebmtcal,
                j=1,
                s=100,
                t = t.eval,
                tp.pred = tp.pred.s100,
                n.pctls = 2)


dat.calib.pv.new <-
  calibmsm(data.mstate = msebmtcal,
            data.raw = ebmtcal,
            j=1,
            s=100,
            t = t.eval,
            tp.pred = tp.pred.s100,
            calib.type = "pv",
            curve.type = "rcs",
            pv.n.pctls = 2)

head(dat.calib.pv.old[["plotdata"]][[2]])
head(dat.calib.pv.new[["plotdata"]][[2]])

testthat::expect_equal(dat.calib.pv.old[["plotdata"]][[2]], dat.calib.pv.new[["plotdata"]][[2]])


################
### PV NO CI ###
### With group.vars ###
### With n.pctls ######
#####################

###
#### RCS
###
dat.calib.pv.old <-
  calib_pv_SAVE(data.mstate = msebmtcal,
                data.raw = ebmtcal,
                j=1,
                s=100,
                t = t.eval,
                tp.pred = tp.pred.s100,
                group.vars = c("year"),
                n.pctls = 2)


dat.calib.pv.new <-
  calibmsm(data.mstate = msebmtcal,
            data.raw = ebmtcal,
            j=1,
            s=100,
            t = t.eval,
            tp.pred = tp.pred.s100,
            calib.type = "pv",
            curve.type = "rcs",
            pv.group.vars = c("year"),
            pv.n.pctls = 2)

head(dat.calib.pv.old[["plotdata"]][[2]])
head(dat.calib.pv.new[["plotdata"]][[2]])

testthat::expect_equal(dat.calib.pv.old[["plotdata"]][[2]], dat.calib.pv.new[["plotdata"]][[2]])

