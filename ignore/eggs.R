# Extract the predicted transition probabilities out of state j = 1
tp.pred <- dplyr::select(dplyr::filter(tps0, j == 1), any_of(paste("pstate", 1:6, sep = "")))

# Now estimate the observed event probabilities for each possible transition.
dat.calib <-
calib_msm(data.mstate = msebmtcal,
 data.raw = ebmtcal,
 j=1,
 s=0,
 t = 1826,
 tp.pred = tp.pred,
 w.covs = c("year", "agecl", "proph", "match"))

# Summarise the output
plot(dat.calib)
