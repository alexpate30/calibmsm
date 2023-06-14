###
### Testing package functionality
###
testchange <- 1
a <- 3
b <- 4
rm(list=ls())
load_all()
devtools::document()
devtools::check()
devtools::install()
small.change <- 1
devtools::check(vignettes = FALSE, args = c("--no-tests"))
?calibmsm::calc_calib_mlr
testhaschanged <- 1
whyisntworkflowbeingrecognised <- 3
library(devtools)
data("ebmtcal")
data("msebmtcal")
str(msebmtcal)
data("tps0")
data("tps100")
str(class(msebmtcal))
typeof(msebmtcal)
rownames(ebmtcal) <- NULL
rownames(tps0) <- NULL
rownames(tps100) <- NULL
head(ebmtcal)
head(tps0)
?dplyr::mutate
devtools::document()


### Want to submit to GitHub to test building of package


Sys.getenv("R_ENVIRON")
### Test trying to add transparency
load_all()

tp.pred <- tps100 |>
  dplyr::filter(j == 3) |>
  dplyr::select(any_of(paste("pstate", 1:6, sep = "")))

dat.calib.blr <-
  calib_blr(data.mstate = msebmtcal,
            data.raw = ebmtcal,
            j=3,
            s=100,
            t.eval = 1826,
            tp.pred = tp.pred,
            curve.type = "rcs",
            rcs.nk = 3)

dat.calib.mlr <-
  calib_mlr(data.mstate = msebmtcal,
            data.raw = ebmtcal,
            j=3,
            s=100,
            t.eval = 1826,
            tp.pred = tp.pred)

png("workflow/figures/plotmlralpha01.png", width = 7.5, height = 5, units = "in", res = 300)
plot.calib_mlr(dat.calib.mlr)
dev.off()

png("workflow/figures/plotmlralpha09.png", width = 7.5, height = 5, units = "in", res = 300)
plot.calib_mlr(dat.calib.mlr, transparency.plot = 0.9)
dev.off()

logit.func <- function(x){
  return(log(x/(1-x)))}
inv.logit.func <- function(x){
  return(1/(1+exp(-x)))
}

logit.func(3)

inv.logit.func(-100)
logit.func(0.7)
inv.logit.func(0.84)



####################
### Test weights ###
####################

weights.manual <- calc_weights(data.mstate = msebmtcal,
                               data.raw = ebmtcal,
                               covs = c("year", "agecl", "proph", "match"),
                               t.eval = 1826,
                               s = 0,
                               landmark.type = "state",
                               j = 1,
                               max.weight = 10,
                               stabilised = FALSE)$ipcw

weights.manual.2 <- calc_weights(data.mstate = msebmtcal,
                                 data.raw = ebmtcal,
                                 covs = c("year", "agecl", "proph", "match"),
                                 t.eval = 1826,
                                 s = 0,
                                 landmark.type = "state",
                                 j = 1,
                                 max.weight = 10,
                                 stabilised = FALSE,
                                 max.follow = "t.eval")$ipcw

temp <- cbind(ebmtcal, weights.manual)
temp2 <- cbind(ebmtcal, weights.manual.2)

head(weights.manual)
head(weights.manual.2)
hist(weights.manual, breaks = 50)
hist(weights.manual.2, breaks = 50)
sum(is.na(weights.manual))
sum(is.na(weights.manual.2))
mean()

max(weights.manual, na.rm = TRUE)

##########################
### STAB WEIGHTS STUFF ###
##########################

dat.calib.blr.unstab <- calc_calib_blr(data.mstate = msebmtcal,
                                data.raw = ebmtcal[1:500, ],
                                j=1,
                                s=0,
                                t.eval = t.eval,
                                tp.pred = tps0 %>% dplyr::filter(j == 1) %>% dplyr::select(any_of(paste("pstate", 1:6, sep = ""))) %>% slice(1:500),
                                curve.type = "rcs",
                                rcs.nk = 3,
                                w.covs = c("year", "agecl", "proph", "match"),
                                CI = 95,
                                CI.R.boot = 200)

dat.calib.blr.stab <- calc_calib_blr(data.mstate = msebmtcal,
                                 data.raw = ebmtcal[1:500, ],
                                 j=1,
                                 s=0,
                                 t.eval = t.eval,
                                 tp.pred = tps0 %>% dplyr::filter(j == 1) %>% dplyr::select(any_of(paste("pstate", 1:6, sep = ""))) %>% slice(1:500),
                                 curve.type = "rcs",
                                 rcs.nk = 3,
                                 w.covs = c("year", "agecl", "proph", "match"),
                                 w.stabilised= TRUE,
                                 CI = 95,
                                 CI.R.boot = 200)



plot(dat.calib.blr.unstab[["plotdata"]][[6]]$obs.upper - dat.calib.blr.unstab[["plotdata"]][[6]]$obs.lower, dat.calib.blr.stab[["plotdata"]][[6]]$obs.upper - dat.calib.blr.stab[["plotdata"]][[6]]$obs.lower,
     xlab = "UNSTABILISED", ylab = "STABILISED")
abline(0,1)


testtt3 <- calc_weights(data.mstate = msebmtcal,
                       data.raw = ebmtcal,
                       covs = c("agecl", "year"),
                       j = 1,
                       landmark.type = "all",
                       s = 0,
                       t.eval = 1826,
                       max.weight = 10,
                       stabilised = FALSE)

weights1 <- calc_weights(data.mstate = msebmtcal,
                        data.raw = ebmtcal,
                        covs = c("agecl", "year", "proph", "match"),
                        t.eval = 1826,
                        s = 0,
                        landmark.type = "all",
                        j = 1,
                        max.weight = 10,
                        stabilised = FALSE)
sum(is.na(weights1))

weights2 <- calc_weights(data.mstate = msebmtcal,
                         data.raw = ebmtcal,
                         covs = c("agecl", "year", "proph", "match"),
                         t.eval = t.eval,
                         s = 0,
                         landmark.type = "all",
                         j = 1,
                         max.weight = 10,
                         stabilised = FALSE)


weights3 <- calc_weights(data.mstate = data.mstate,
                        data.raw = data.boot,
                        covs = w.covs,
                        t.eval = t.eval,
                        s = s,
                        landmark.type = w.landmark.type,
                        j = j,
                        max.weight = w.max,
                        stabilised = w.stabilised)


getwd()
t.eval <- 1826
str(ebmtcal)

a <- seq(0.1, 0.9, 0.1)
log(a/(1-a))

load_all()
dat.calib.blr <- calc_calib_blr(data.mstate = msebmtcal,
                                data.raw = ebmtcal,
                                j=1,
                                s=0,
                                t.eval = 1826,
                                tp.pred = tps0 %>% dplyr::filter(j == 1) %>% dplyr::select(any_of(paste("pstate", 1:6, sep = ""))),
                                curve.type = "rcs",
                                rcs.nk = 3,
                                w.covs = c("year", "agecl", "proph", "match"))
plot(dat.calib.blr, combine = TRUE, nrow = 2, ncol = 3)
plot.calib_blr(dat.calib.blr, combine = TRUE, nrow = 2, ncol = 3)
ftype(plot)
class(dat.calib.blr)
library(devtools)
install.packages("devtools")
install.packages("usethis")
install.packages("rlang")
install.packages("Rtools")
install.packages("stringr")
attributes(msebmtcal)
rm(list=ls())
devtools::install()
sessionInfo()
mypaths <- .libPaths()
mypaths
.libPaths()
.libPaths(mypaths[2])
install.packages("rlang")
library(calibmsm)

install.packages("rlang")
devtools::document()
getwd()
rm(list=ls())
load_all()
calc_calib_mlr


devtools::test()
.libPaths()
devtools::install()
testthat::test_file("tests/testthat/test-calib_blr.R")
testthat::test_file("tests/testthat/test-calib_mlr.R")
testthat::test_file("tests/testthat/test-weights.R")
testthat::test_file("tests/testthat/test-calib_pv.R")
testthat::test_file("tests/testthat/test-plot_calib.R")


### Testing manually function for weights with internal CI
load_all()

tp.pred <- dplyr::select(dplyr::filter(tps0, j == 1), dplyr::any_of(paste("pstate", 1:6, sep = "")))

dat.calib.blr <-
  calc_calib_blr(data.mstate = msebmtcal,
                 data.raw = ebmtcal,
                 j=1,
                 s=0,
                 t.eval = 1826,
                 tp.pred = tp.pred,
                 curve.type = "rcs",
                 rcs.nk = 3)

## Calculate observed event probabilities
dat.calib.blr.w.function <-
  calc_calib_blr(data.mstate = msebmtcal,
                 data.raw = ebmtcal,
                 j=1,
                 s=0,
                 t.eval = 1826,
                 tp.pred = tp.pred,
                 curve.type = "rcs",
                 rcs.nk = 3,
                 w.function = calc_weights)

str(dat.calib.blr[["plotdata"]])
str(dat.calib.blr.w.function[["plotdata"]])


### Test summar yfunction for pv
load_all()
## Extract relevant predicted risks from tps0
tp.pred <- dplyr::select(dplyr::filter(tps100, j == 3), dplyr::any_of(paste("pstate", 1:6, sep = "")))

## Calculate observed event probabilities using transitions.out = NULL
dat.calib.pv.1 <- calc_calib_pv(data.mstate = msebmtcal,
                                data.raw = ebmtcal,
                                j = 3,
                                s = 100,
                                t.eval = 1826,
                                tp.pred = tp.pred,
                                curve.type = "loess",
                                group.vars = c("year"),
                                n.pctls = 2,
                                data.pred.plot = NULL, transitions.out = NULL)
str(dat.calib.pv.1[["metadata"]])
summary.calib_pv(dat.calib.pv.1)

setequal(c(1,2,3,4,5), c(3,4,2,5,1))
library(calibmsm)

load_all()
dat.calib.mlr.j1.s0.s <- calc_calib_mlr(data.mstate = msebmtcal,
                                      data.raw = ebmtcal,
                                      j=1,
                                      s=0,
                                      t.eval = 1826,
                                      tp.pred = tps0 %>% dplyr::filter(j == 1) %>% dplyr::select(any_of(paste("pstate", 1:6, sep = ""))),
                                      smoother.type = "s",
                                      ps.int = 3, degree = 2,
                                      w.covs = c("year", "agecl", "proph", "match"),
                                      w.landmark.type = "all")

dat.calib.mlr.j1.s0.sm.ps <- calc_calib_mlr(data.mstate = msebmtcal,
                                        data.raw = ebmtcal,
                                        j=1,
                                        s=0,
                                        t.eval = 1826,
                                        tp.pred = tps0 %>% dplyr::filter(j == 1) %>% dplyr::select(any_of(paste("pstate", 1:6, sep = ""))),
                                        smoother.type = "sm.ps",
                                        ps.int = 3, degree = 2,
                                        w.covs = c("year", "agecl", "proph", "match"),
                                        w.landmark.type = "all")

dat.calib.mlr.j1.s0.sm.os <- calc_calib_mlr(data.mstate = msebmtcal,
                                        data.raw = ebmtcal,
                                        j=1,
                                        s=0,
                                        t.eval = 1826,
                                        tp.pred = tps0 %>% dplyr::filter(j == 1) %>% dplyr::select(any_of(paste("pstate", 1:6, sep = ""))),
                                        smoother.type = "sm.os",
                                        w.covs = c("year", "agecl", "proph", "match"),
                                        w.landmark.type = "all")


plot.calib_mlr(dat.calib.mlr.j1.s0.s)
plot.calib_mlr(dat.calib.mlr.j1.s0.sm.ps)
plot.calib_mlr(dat.calib.mlr.j1.s0.sm.os)


dat.calib.mlr.j1.s0 <- calc_calib_mlr(data.mstate = msebmtcal,
                                      data.raw = ebmtcal,
                                      j=1,
                                      s=0,
                                      t.eval = 1826,
                                      tp.pred = tps0 %>% dplyr::filter(j == 1) %>% dplyr::select(any_of(paste("pstate", 1:6, sep = ""))),
                                      ps.int = 3, degree = 2,
                                      w.covs = c("year", "agecl", "proph", "match"),
                                      w.landmark.type = "all")
plot.calib_mlr(dat.calib.mlr.j1.s0)


dat.calib.mlr.j1.s0 <- calc_calib_mlr(data.mstate = msebmtcal,
                                      data.raw = ebmtcal,
                                      j=1,
                                      s=0,
                                      t.eval = 1826,
                                      tp.pred = tps0 %>% dplyr::filter(j == 1) %>% dplyr::select(any_of(paste("pstate", 1:6, sep = ""))),
                                      ps.int = 3, degree = 2,
                                      w.covs = c("year", "agecl", "proph", "match"),
                                      w.landmark.type = "all")
plot.calib_mlr(dat.calib.mlr.j1.s0)




dat.calib.blr.j1.s100 <- calc_calib_blr(data.mstate = msebmtcal,
                                        data.raw = ebmtcal,
                                        j=1,
                                        s=100,
                                        t.eval = 1826,
                                        tp.pred = tps100 %>% dplyr::filter(j == 1) %>% dplyr::select(any_of(paste("pstate", 1:6, sep = ""))),
                                        curve.type = "rcs",
                                        rcs.nk = 3,
                                        w.covs = c("year", "agecl", "proph", "match"),
                                        w.landmark.type = "all")

dat.calib.blr.j2.s100 <- calc_calib_blr(data.mstate = msebmtcal,
                                        data.raw = ebmtcal,
                                        j=2,
                                        s=100,
                                        t.eval = t.eval,
                                        tp.pred = tps100 %>% filter(j == 2) %>% select(any_of(paste("pstate", 1:6, sep = ""))),
                                        curve.type = "rcs",
                                        rcs.nk = 3,
                                        w.covs = c("year", "agecl", "proph", "match"))

load_all()
dat.calib.blr.j3.s100 <- calc_calib_pv(data.mstate = msebmtcal,
                                        data.raw = ebmtcal,
                                        j=3,
                                        s=100,
                                        t.eval = 1826,
                                        tp.pred = tps100 %>% dplyr::filter(j == 3) %>% dplyr::select(any_of(paste("pstate", 1:6, sep = ""))),
                                        curve.type = "rcs",
                                        rcs.nk = 3,
                                       group.vars = c("year"),
                                       n.pctls = 2)

temp.func <- function(x){
  if (x> 3){
    warning("abcdefg")
  }
}
grepl("abc", "abcdefg")

temp.func(1)
temp.func(4)
warn.func <- function
suppressWarnings(temp.func(4), .f = function(x){grepl("abc", x)})
plot(dat.calib.blr.j1.s100, combine = TRUE, nrow = 2, ncol = 3)
plot(dat.calib.blr.j2.s100, combine = TRUE, nrow = 2, ncol = 3)
plot(dat.calib.blr.j3.s100, combine = TRUE, nrow = 2, ncol = 3)
data("ebmt")
data("ebmtcal")

expect_error(1 / "a")

#   data.mstate <- msebmt
#   data.raw <- ebmt
#   covs <- NULL
#   j
#   landmark.type <- "state"
#   s
#   t.eval <- t.eval
#   max.weight <- 10
#
