###
### Testing the pseudo-value programs
###
rm(list=ls())
load_all()
msebmtcal

####################
### NORMAL
#########################

###
### First must calculate data.pred.plot
###
id.lmk <- extract_ids_states(data.mstate = msebmtcal,
                             tmat = attributes(msebmtcal)$trans,
                             j = 3,
                             t.eval = 100)

data.pred.plot <- tps100 %>%
  dplyr::filter(id %in% id.lmk) %>%
  dplyr::filter(j == 3) %>%
  dplyr::select(any_of(paste("pstate", 1:6, sep = "")))

test1 <-calc_calib_pv(data.mstate = msebmtcal,
                     data.raw = ebmtcal,
                     j = 3,
                     s = 100,
                     t.eval = 1826,
                     tp.pred = tps100 %>%
                       dplyr::filter(j == 3) %>%
                       dplyr::select(any_of(paste("pstate", 1:6, sep = ""))),
                     group.vars = c("year"),
                     n.pctls = 2,
                     data.pred.plot = NULL,
                     transitions.out = 3)

test2 <-calc_calib_pv(data.mstate = msebmtcal,
                     data.raw = ebmtcal,
                     j = 3,
                     s = 100,
                     t.eval = 1826,
                     tp.pred = tps100 %>%
                       dplyr::filter(j == 3) %>%
                       dplyr::select(any_of(paste("pstate", 1:6, sep = ""))),
                     group.vars = c("year"),
                     n.pctls = 2,
                     data.pred.plot = data.pred.plot,
                     transitions.out = 3)

str(test1[["plotdata"]])
str(test2[["plotdata"]])
temp <- plot.calib_pseudo(test)

### Create plots
Cairo::CairoPNG(paste("workflow/figures/supp_pseudo_TEST_j3s100.png", sep = ""),
                dpi = 300, width = 15, height = 10, unit = "in")
print(temp)
dev.off()



##########################################
### BOOTSTRAP ###
####################################

###
### Write a function which will do a bootstrap
###
calc_calib_pseudo_boot <- function(data.in, indices){

  ### Create bootstrapped dataset
  data.raw.boot <- data.in[indices, ]

  ### Create bootstrapped data.mstate
  data.mstate.boot <-
    do.call("rbind", lapply(
      data.raw.boot$id, function(x) {base::subset(msebmtcal, id == x)})
      )

  ### Apply attribute tmat
  attributes(data.mstate.boot)$trans <- attributes(msebmtcal)$trans

  ### Calc calibration for data.boot
  calib.boot <- calc_calib_pv(data.mstate = data.mstate.boot,
                              data.raw = data.raw.boot,
                              j = 3,
                              s = 100,
                              t.eval = 1826,
                              tp.pred = tps100 %>%
                                dplyr::filter(j == 3) %>%
                                dplyr::select(any_of(paste("pstate", 1:6, sep = ""))),
                              group.vars = c("year"),
                              n.pctls = 2,
                              data.pred.plot = NULL,
                              transitions.out = 3)

  ### Return the calibration of interest
  return(calib.boot[["plotdata"]][["state3"]]$obs)

}

library(boot)
test <- boot::boot(ebmtcal, calc_calib_pseudo_boot, R = 3)

temp <- test$t




# data.mstate <- msebmtcal
# data.raw <- ebmtcal
# j <- 3
# s <- 100
# t.evl <- 1826
# tp.pred <-
# data.pred.plot <- NULL
# transitions.out <- NULL


calc.calib.aj.ce.new2(msebmtcal, attributes(msebmtcal)$trans, t.eval = 1826, j = 3, s = 100)

calc.calib.aj.ce.new2(msebmtcal, attributes(msebmtcal)$trans, t.eval = 1826, j = 3, s = 100)
calc.calib.aj.ce.new2(msebmtcal, attributes(msebmtcal)$trans, t.eval = 1826, j = 2, s = 100)
calc.calib.aj.ce.new2(msebmtcal, attributes(msebmtcal)$trans, t.eval = 1826, j = 1, s = 100)


### WARNING MESSAGES TO BE AWARE OF:
# Warning messages:
# 1: In min(diff(time)) : no non-missing arguments to min; returning Inf
# 2: In min(diff(time)) : no non-missing arguments to min; returning Inf
# 3: In max(x[!is.na(x)]) : no non-missing arguments to max; returning -Inf
# 4: In max(x[!is.na(x)]) : no non-missing arguments to max; returning -Inf

### The latter is when there is a transition of which no individuals make. Not a problem.
### Need to identify the former.

### The former occurs when there is a state where no individuals move out of it
