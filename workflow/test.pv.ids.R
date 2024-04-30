devtools::load_all()
# Estimate BLR-IPCW calibration curves for the predicted transition
# probabilities at time t = 1826, when predictions were made at time
# s = 0 in state j = 1. These predicted transition probabilities are stored in tps0.
#'
# Extract the predicted transition probabilities out of state j = 1
tp.pred <- dplyr::select(dplyr::filter(tps0, j == 1), any_of(paste("pstate", 1:6, sep = "")))


test <- msebmtcal[msebmtcal$id %in% 1:100, ]
str(test)
#'
# Now estimate the observed event probabilities for each possible transition.
dat.calib <- calib_msm(data.mstate = msebmtcal[msebmtcal$id %in% 1:100, ],
 data.raw = ebmtcal[ebmtcal$id %in% 1:100, ],
 j=1,
 s=0,
 t = 1826,
 tp.pred = tp.pred[1:100, ],
 calib.type = "pv")

dat.calib <- calib_msm(data.mstate = msebmtcal[msebmtcal$id %in% 1:50, ],
                       data.raw = ebmtcal[ebmtcal$id %in% 1:50, ],
                       j=1,
                       s=0,
                       t = 1826,
                       tp.pred = tp.pred[1:50, ],
                       calib.type = "pv",
                       pv.group.vars = c("match"),
                       pv.ids = NULL)
str(dat.calib[[1]])


dat.calib2 <- calib_msm(data.mstate = msebmtcal[msebmtcal$id %in% 1:50, ],
                       data.raw = ebmtcal[ebmtcal$id %in% 1:50, ],
                       j=1,
                       s=0,
                       t = 1826,
                       tp.pred = tp.pred[1:50, ],
                       calib.type = "pv",
                       pv.group.vars = c("match"),
                       pv.ids = 1:3)
str(dat.calib2[[1]])


dat.calib <- calib_msm(data.mstate = msebmtcal[msebmtcal$id %in% 1:50, ],
                       data.raw = ebmtcal[ebmtcal$id %in% 1:50, ],
                       j=1,
                       s=0,
                       t = 1826,
                       tp.pred = tp.pred[1:50, ],
                       calib.type = "pv",
                       pv.group.vars = c("match"),
                       pv.n.pctls = 2,
                       pv.ids = NULL)
str(dat.calib[[1]])


dat.calib2 <- calib_msm(data.mstate = msebmtcal[msebmtcal$id %in% 1:50, ],
                        data.raw = ebmtcal[ebmtcal$id %in% 1:50, ],
                        j=1,
                        s=0,
                        t = 1826,
                        tp.pred = tp.pred[1:50, ],
                        calib.type = "pv",
                        pv.group.vars = c("match"),
                        pv.n.pctls = 2,
                        pv.ids = 1:3)
str(dat.calib2[[1]])


dat.calib <- calib_msm(data.mstate = msebmtcal[msebmtcal$id %in% 1:50, ],
                       data.raw = ebmtcal[ebmtcal$id %in% 1:50, ],
                       j=1,
                       s=0,
                       t = 1826,
                       tp.pred = tp.pred[1:50, ],
                       calib.type = "pv",
                       pv.n.pctls = 2,
                       pv.ids = NULL)
str(dat.calib[[1]])


dat.calib2 <- calib_msm(data.mstate = msebmtcal[msebmtcal$id %in% 1:50, ],
                        data.raw = ebmtcal[ebmtcal$id %in% 1:50, ],
                        j=1,
                        s=0,
                        t = 1826,
                        tp.pred = tp.pred[1:50, ],
                        calib.type = "pv",
                        pv.n.pctls = 2,
                        pv.ids = 1:3)
str(dat.calib2[[1]])

#'
# Summarise the output
summary(dat.calib)
#'
