###################################################################################################################################
### Code to prepare datasets required for the vignette: Comparison-with-graphical-calibration-curves-in-competing-risks-setting ###
###################################################################################################################################

### This code involves generating predicted transition probabilities for a competing risks sub-model out of the initial state.
### All other states are considered absorbing states.

### Clear workspace
rm(list=ls())

### Load mstate
library(mstate)
data("ebmt4")
ebmt <- ebmt4

### Define state which sub-model is coming out of and landmark time
j <- 1
s <- 0

### Define t.eval
t.eval <- 1826

### Define tmat
tmat <- mstate::transMat(x = list(c(2, 3, 5, 6), c(), c(), c(),
                                  c(), c()), names = c("Tx", "Rec", "AE", "Rec+AE", "Rel", "Death"))

### Create data in msdata format
msebmtcal_cmprsk <- mstate::msprep(data = ebmt, trans = tmat, time = c(NA, "rec", "ae","recae", "rel", "srv"),
                                     status = c(NA, "rec.s", "ae.s", "recae.s", "rel.s", "srv.s"),
                                     keep = c("match", "proph", "year", "agecl"))

### Define covariates for model
covs <- c("match", "proph", "year", "agecl")
msebmtcal_cmprsk <- mstate::expand.covs(msebmtcal_cmprsk, covs, longnames = FALSE)

### Assign variables for model we will be fitting
eq.RHS <- paste(do.call(paste0, expand.grid(c("match", "proph", "year1", "year2", "agecl1", "agecl2"), paste(".", 1:sum(!is.na(tmat)), sep = ""))), collapse="+")
strata <- survival::strata
eq <- paste("survival::Surv(Tstart, Tstop, status) ~ ", eq.RHS,  "+ strata(trans)", sep = "")
eq <- as.formula(eq)

### Create dataframe to store predicted risks
tp.all <- data.frame(matrix(NA, ncol = 13, nrow = nrow(ebmt)))
colnames(tp.all) <- c("id", paste("pstate", 1:6, sep = ""), paste("se", 1:6, sep = ""))

### Loop through id.iter
for (id.iter in 1:nrow(ebmt)){

  print(paste("id.iter = ", id.iter, Sys.time()))

  ### Develop a model on entire dataset except individual of interest
  cfull <- survival::coxph(eq, data = subset(msebmtcal_cmprsk, id != id.iter), method = "breslow")

  ### Get location of individual in msebmtcal_cmprsk
  pat.loc <- which(msebmtcal_cmprsk$id == id.iter)

  ### Create a miniture dataset, on which to generate predictions in (must be in mstate format and have a row for every transition)
  pat.dat <- msebmtcal_cmprsk[rep(pat.loc[1], sum(!is.na(tmat))), 9:12]
  pat.dat$trans <- 1:sum(!is.na(tmat))
  attr(pat.dat, "trans") <- tmat
  pat.dat <- mstate::expand.covs(pat.dat, covs, longnames = FALSE)
  pat.dat$strata <- pat.dat$trans

  ### Fit cause-specific hazards
  msf.pat <- mstate::msfit(cfull, pat.dat, trans = tmat)

  ### Generate 5 year transition probabilities for this patient from times s = 0
  pt <- mstate::probtrans(msf.pat, predt = 0)

  ### Write a function to extract the transition probabilities from state j into each state, after followup time f.time
  extract.tp <- function(tp.object, state, f.time){
    ### Create output object
    output.object <- return(base::subset(tp.object[[state]], time > f.time) |> dplyr::slice(1) |> dplyr::select(-c(time)))
  }

  ### Calculate required transition probabilities and store in output dataset
  tp.all[id.iter, ] <- c(id.iter, extract.tp(tp.object = pt, state = j, f.time = t.eval - s))

}

### Rename datasets and get in correct formats
tp_cmprsk_j0 <- tp.all
msebmtcal_cmprsk <- dplyr::select(msebmtcal_cmprsk, c("id", "from", "to", "trans", "Tstart", "Tstop", "time", "status"))
attributes(msebmtcal_cmprsk)$trans <- tmat

### Finally create a new dtcens variable for the competing risk data
### Having an event of interest stops censoring from being observed, therefore the
### competing risks data requires a new censoring variable
ebmtcal_cmprsk <- ebmt
ebmtcal_cmprsk$dtcens <- pmin(ebmtcal_cmprsk$rec, ebmtcal_cmprsk$ae, ebmtcal_cmprsk$rel, ebmtcal_cmprsk$srv)
ebmtcal_cmprsk$dtcens_s <- 1 - pmax(ebmtcal_cmprsk$rec.s, ebmtcal_cmprsk$ae.s, ebmtcal_cmprsk$rel.s, ebmtcal_cmprsk$srv.s)

### Use in package
usethis::use_data(tp_cmprsk_j0, overwrite = TRUE)
usethis::use_data(msebmtcal_cmprsk, overwrite = TRUE)
usethis::use_data(ebmtcal_cmprsk, overwrite = TRUE)
