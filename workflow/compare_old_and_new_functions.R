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


######################
### BLR-IPCW NO CI ###
######################

###
#### RCS
###
dat.calib.blr.old <-
  calib_blr_SAVE(data.mstate = msebmtcal,
            data.raw = ebmtcal,
            j=1,
            s=0,
            t = t.eval,
            tp.pred = tp.pred.s0,
            curve.type = "rcs",
            w.covs = c("year", "agecl", "proph", "match"))


dat.calib.blr.new <-
  calibmsm(data.mstate = msebmtcal,
                data.raw = ebmtcal,
                j=1,
                s=0,
                t = t.eval,
                tp.pred = tp.pred.s0,
                calib.type = "blr",
                curve.type = "rcs",
                w.covs = c("year", "agecl", "proph", "match"))

head(dat.calib.blr.old[["plotdata"]][[2]])
head(dat.calib.blr.new[["plotdata"]][[2]])

testthat::expect_equal(dat.calib.blr.old[["plotdata"]][[2]], dat.calib.blr.new[["plotdata"]][[2]])

###
#### LOESS
###
dat.calib.blr.old <-
  calib_blr_SAVE(data.mstate = msebmtcal,
            data.raw = ebmtcal,
            j=1,
            s=0,
            t = t.eval,
            tp.pred = tp.pred.s0,
            curve.type = "loess",
            w.covs = c("year", "agecl", "proph", "match"))


dat.calib.blr.new <-
  calibmsm(data.mstate = msebmtcal,
                data.raw = ebmtcal,
                j=1,
                s=0,
                t = t.eval,
                tp.pred = tp.pred.s0,
                calib.type = "blr",
                curve.type = "loess",
                w.covs = c("year", "agecl", "proph", "match"))

head(dat.calib.blr.old[["plotdata"]][[2]])
head(dat.calib.blr.new[["plotdata"]][[2]])

testthat::expect_equal(dat.calib.blr.old[["plotdata"]][[2]], dat.calib.blr.new[["plotdata"]][[2]])

########################
### BLR-IPCW WITH CI ###
########################

###
#### RCS
###
dat.calib.blr.old <-
  calib_blr_SAVE(data.mstate = msebmtcal,
            data.raw = ebmtcal,
            j=1,
            s=0,
            t = t.eval,
            tp.pred = tp.pred.s0,
            curve.type = "rcs",
            w.covs = c("year", "agecl", "proph", "match"),
            CI = 95,
            CI.R.boot = 10)

dat.calib.blr.new <-
  calibmsm(data.mstate = msebmtcal,
                data.raw = ebmtcal,
                j=1,
                s=0,
                t = t.eval,
                tp.pred = tp.pred.s0,
                calib.type = "blr",
                curve.type = "rcs",
                w.covs = c("year", "agecl", "proph", "match"),
                CI = 95,
                CI.R.boot = 10,
                CI.seed = 1)


testthat::expect_equal(dat.calib.blr.old[["plotdata"]][[2]]$obs, dat.calib.blr.new[["plotdata"]][[2]]$obs)

###
#### LOESS
###
dat.calib.blr.old <-
  calib_blr_SAVE(data.mstate = msebmtcal,
            data.raw = ebmtcal,
            j=1,
            s=0,
            t = t.eval,
            tp.pred = tp.pred.s0,
            curve.type = "loess",
            w.covs = c("year", "agecl", "proph", "match"),
            CI = 95,
            CI.R.boot = 10)


dat.calib.blr.new <-
  calibmsm(data.mstate = msebmtcal,
                data.raw = ebmtcal,
                j=1,
                s=0,
                t = t.eval,
                tp.pred = tp.pred.s0,
                curve.type = "loess",
                calib.type = "blr",
                w.covs = c("year", "agecl", "proph", "match"),
                CI = 95,
                CI.R.boot = 10)

head(dat.calib.blr.old[["plotdata"]][[2]])
head(dat.calib.blr.new[["plotdata"]][[2]])

testthat::expect_equal(dat.calib.blr.old[["plotdata"]][[2]]$obs, dat.calib.blr.new[["plotdata"]][[2]]$obs)

#############################
### DATA PRED PLOT, CI    ###
### DATA.PRED.PLOT=tps0   ###
#############################

### Now if we specify data.pred.plot I think we might get a difference becuase of an error in the original code...
## Extract ids for individuals uncensored at t
ids.uncens <- ebmtcal |>
  subset(dtcens > t.eval | (dtcens < t.eval & dtcens.s == 0)) |>
  dplyr::pull(id)
## Extract the predicted risks out of state 1 for these individuals
data.pred.plot <- tps0 |>
  dplyr::filter(j == 1 & id %in% ids.uncens) |>
  dplyr::select(any_of(paste("pstate", 1:6, sep = "")))
## Now pick at random some rows to make it different from the normal cohort
data.pred.plot2 <- dplyr::sample_n(data.pred.plot,100)

###
#### RCS
###
dat.calib.blr.new1 <-
  calibmsm(data.mstate = msebmtcal,
                data.raw = ebmtcal,
                j=1,
                s=0,
                t = t.eval,
                tp.pred = tp.pred.s0,
                calb.type = "blr",
                curve.type = "rcs",
                w.covs = c("year", "agecl", "proph", "match"),
                CI = 95,
                CI.R.boot = 10)


dat.calib.blr.new2 <-
  calibmsm(data.mstate = msebmtcal,
                data.raw = ebmtcal,
                j=1,
                s=0,
                t = t.eval,
                tp.pred = tp.pred.s0,
                calb.type = "blr",
                curve.type = "rcs",
                w.covs = c("year", "agecl", "proph", "match"),
                CI = 95,
                CI.R.boot = 10,
                data.pred.plot = data.pred.plot)

testthat::expect_equal(dat.calib.blr.new1[["plotdata"]][[2]]$obs, dat.calib.blr.new2[["plotdata"]][[2]]$obs)

head(dat.calib.blr.new1[["plotdata"]][[2]])
head(dat.calib.blr.new2[["plotdata"]][[2]])

###
#### LOESS
###
dat.calib.blr.new1 <-
  calibmsm(data.mstate = msebmtcal,
                data.raw = ebmtcal,
                j=1,
                s=0,
                t = t.eval,
                tp.pred = tp.pred.s0,
                calb.type = "blr",
                curve.type = "loess",
                w.covs = c("year", "agecl", "proph", "match"),
                CI = 95,
                CI.R.boot = 10)


dat.calib.blr.new2 <-
  calibmsm(data.mstate = msebmtcal,
                data.raw = ebmtcal,
                j=1,
                s=0,
                t = t.eval,
                tp.pred = tp.pred.s0,
                calb.type = "blr",
                curve.type = "loess",
                w.covs = c("year", "agecl", "proph", "match"),
                CI = 95,
                CI.R.boot = 10,
                data.pred.plot = data.pred.plot)

str(dat.calib.blr.new1)
str(dat.calib.blr.new2)

testthat::expect_equal(dat.calib.blr.new1[["plotdata"]][[2]]$obs, dat.calib.blr.new2[["plotdata"]][[2]]$obs)
testthat::expect_equal(dat.calib.blr.new1[["plotdata"]][[2]]$upper, dat.calib.blr.new2[["plotdata"]][[2]]$upper)
testthat::expect_equal(dat.calib.blr.new1[["plotdata"]][[2]]$lower, dat.calib.blr.new2[["plotdata"]][[2]]$lower)

#############################
### DATA PRED PLOT, CI    ###
### DATA.PRED.PLOT=new    ###
#############################

## Now pick at random some rows to make it different from the normal cohort
data.pred.plot2 <- dplyr::sample_n(data.pred.plot,100)

###
#### RCS NO CI
###
dat.calib.blr.new1 <-
  calibmsm(data.mstate = msebmtcal,
                data.raw = ebmtcal,
                j=1,
                s=0,
                t = t.eval,
                tp.pred = tp.pred.s0,
                calb.type = "blr",
                curve.type = "rcs",
                w.covs = c("year", "agecl", "proph", "match"),
                data.pred.plot = data.pred.plot2)

###
#### RCS CI
###
dat.calib.blr.new2 <-
  calibmsm(data.mstate = msebmtcal,
                data.raw = ebmtcal,
                j=1,
                s=0,
                t = t.eval,
                tp.pred = tp.pred.s0,
                calb.type = "blr",
                curve.type = "rcs",
                w.covs = c("year", "agecl", "proph", "match"),
                CI = 95,
                CI.R.boot = 10,
                data.pred.plot = data.pred.plot2)

str(dat.calib.blr.new1)
str(dat.calib.blr.new2)

testthat::expect_equal(dat.calib.blr.new1[["plotdata"]][[2]], dat.calib.blr.new2[["plotdata"]][[2]][,c(1,2)])

###
#### LOESS NO CI
###
dat.calib.blr.new1 <-
  calibmsm(data.mstate = msebmtcal,
                data.raw = ebmtcal,
                j=1,
                s=0,
                t = t.eval,
                tp.pred = tp.pred.s0,
                calb.type = "blr",
                curve.type = "loess",
                w.covs = c("year", "agecl", "proph", "match"),
                data.pred.plot = data.pred.plot2)

###
#### LOESS CI
###
dat.calib.blr.new2 <-
  calibmsm(data.mstate = msebmtcal,
                data.raw = ebmtcal,
                j=1,
                s=0,
                t = t.eval,
                tp.pred = tp.pred.s0,
                calb.type = "blr",
                curve.type = "loess",
                w.covs = c("year", "agecl", "proph", "match"),
                CI = 95,
                CI.R.boot = 10,
                data.pred.plot = data.pred.plot2)

str(dat.calib.blr.new1)
str(dat.calib.blr.new2)

testthat::expect_equal(dat.calib.blr.new1[["plotdata"]][[2]], dat.calib.blr.new2[["plotdata"]][[2]][,c(1,2)])


###
### SO FAR, EVERYTHING WORKING AS IT SHOULD
###


###
### Lets try MLR
###
dat.calib.mlr.old <-
  calib_mlr_SAVE(data.mstate = msebmtcal,
            data.raw = ebmtcal,
            j=1,
            s=100,
            t = 1826,
            tp.pred = tp.pred.s100,
            w.covs = c("year", "agecl", "proph", "match"))

dat.calib.mlr.new <-
  calibmsm(data.mstate = msebmtcal,
            data.raw = ebmtcal,
            j=1,
            s=100,
            t = 1826,
            tp.pred = tp.pred.s100,
            calib.type = "mlr",
            w.covs = c("year", "agecl", "proph", "match"))


head(dat.calib.mlr.old[["plotdata"]][[2]])
head(dat.calib.mlr.new[["plotdata"]][[2]])

testthat::expect_equal(dat.calib.mlr.old[["plotdata"]][[2]], dat.calib.mlr.new[["plotdata"]][[2]])



###
### Testing the S3 generics
###

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

###
### BLR-IPCW
###
dat.calib.blr <-
  calibmsm(data.mstate = msebmtcal,
           data.raw = ebmtcal,
           j=1,
           s=0,
           t = t.eval,
           tp.pred = tp.pred.s0,
           calib.type = "blr",
           curve.type = "rcs",
           w.covs = c("year", "agecl", "proph", "match"))

str(dat.calib.blr)
class(dat.calib.blr)
inherits(dat.calib.blr, "calibmsm")
summary(dat.calib.blr)
plot(dat.calib.blr)



dat.calib.mlr <-
  calibmsm(data.mstate = msebmtcal,
           data.raw = ebmtcal,
           j=1,
           s=100,
           t = 1826,
           tp.pred = tp.pred.s100,
           calib.type = "mlr",
           w.covs = c("year", "agecl", "proph", "match"))

class(dat.calib.mlr)
summary(dat.calib.mlr)
plot(dat.calib.mlr)


###
### TO DO
###
help(package = "calibmsm")
