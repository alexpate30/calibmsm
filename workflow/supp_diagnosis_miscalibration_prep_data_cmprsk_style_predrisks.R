### This preliminary piece of work will look at how to assess calibration across a range of sub-times, r.

### Need to start by estimating the predicted transition probabilities at each time r1 and r2.

###

###

######################################################
### Code to prepare datasets provided with package ###
######################################################
rm(list=ls())

### Load the ebmt4 dataset from mstate package and rename it
### Note that mstate is not an import for calibmsm, so this package may need to be installed.
# install.packages("mstate")
library(mstate)
library(calibmsm)
# data("ebmt4")
# ebmt <- ebmt4
#
# ### Define tmat
# tmat <- transMat(x = list(c(2, 3, 5, 6), c(4, 5, 6), c(4, 5, 6), c(5, 6),
#                           c(), c()), names = c("Tx", "Rec", "AE", "Rec+AE", "Rel", "Death"))
# tmat
#
# ### Create mstate format
# msebmt <- msprep(data = ebmt, trans = tmat, time = c(NA, "rec", "ae","recae", "rel", "srv"),
#                  status = c(NA, "rec.s", "ae.s", "recae.s", "rel.s", "srv.s"),
#                  keep = c("match", "proph", "year", "agecl"))
# msebmt[msebmt$id == 1, c(1:8, 10:12)]
#
# ### Define covariates for model
# covs <- c("match", "proph", "year", "agecl")
# msebmt <- expand.covs(msebmt, covs, longnames = FALSE)
# msebmt[msebmt$id == 1, -c(9, 10, 12:48, 61:84)]
#
# ### Change time scale of model into years
# # msebmt[, c("Tstart", "Tstop", "time")] <- msebmt[, c("Tstart","Tstop", "time")]/365.25
#
# ### Create a variable which is maximum observed follow up time for all individuals, this is when they were either censored, relapsed or died
# ebmt$dtcens <- pmin(ebmt$rel, ebmt$srv)
# ebmt$dtcens.s <- 1 - pmax(ebmt$rel.s, ebmt$srv.s)
#
# ### Assign variables for model we will be fitting
# eq.RHS <- paste(do.call(paste0, expand.grid(c("match", "proph", "year1", "year2", "agecl1", "agecl2"), paste(".", 1:12, sep = ""))), collapse="+")
# eq <- paste("Surv(Tstart, Tstop, status) ~ ", eq.RHS,  "+ strata(trans)", sep = "")
# eq <- as.formula(eq)

### Define ranges for r1 and r2
r1.range <- round(365.25*c(0,1,2,3,4))
r2.range <- round(365.25*c(1,2,3,4,5))


### Create dataframe to store predicted risks
tp_pred <- vector("list", 4)
names(tp_pred) <- paste("j =", 1:4)

for (j in 1:4){
  tp_pred[[j]] <- vector("list", length(r1.range))
  names(tp_pred[[j]]) <- paste("r1 =", r1.range)
  for (r1.count in 1:length(r1.range)){
    tp_pred[[j]][[r1.count]] <- vector("list", length(r2.range))
    names(tp_pred[[j]][[r1.count]]) <- paste("r2 =", r2.range)
    for (r2.count in 1:length(r2.range)){
      tp_pred[[j]][[r1.count]][[r2.count]] <- data.frame(matrix(NA, ncol = 16, nrow = nrow(ebmtcal)))
      colnames(tp_pred[[j]][[r1.count]][[r2.count]]) <- c("id", paste("pstate", 1:6, sep = ""), paste("se", 1:6, sep = ""), "j", "r1", "r2")
    }
  }
}


### Definitely need to do each state seperately now
for (j in 1:4){

  ## Extract transition matrix for multistate data
  tmat.msm <- attributes(msebmtcal)$trans

  ## Create list of correct length for input into mstate::transMat
  list.obj <- vector("list", nrow(tmat.msm))

  ## For row j, add the possible transitions
  list.obj[[j]] <- as.numeric(which(!is.na(tmat.msm[j,])))

  ## Define tmat
  tmat <- mstate::transMat(x = list.obj, names = attributes(tmat.msm)$dimnames$from)

  ## Create data in msdata format
  str(ebmtcal)
  data(ebmtcal)
  data(msebmtcal)
  msebmtcal.cmprsk <- mstate::msprep(data = ebmtcal, trans = tmat.msm, time = c(NA, "rec", "ae","recae", "rel", "srv"),
                                     status = c(NA, "rec.s", "ae.s", "recae.s", "rel.s", "srv.s"),
                                     keep = c("match", "proph", "year", "agecl"))

  tmat.cmprsk <- tmat
 str(tmat.cmprsk)
  tmat.cmprsk[2, c(4,5,6)] <- c(5,6,7)
  msebmtcal.cmprsk <- subset(msebmtcal.cmprsk, from == 2)
  attributes(msebmtcal.cmprsk)$trans <- tmat.cmprsk
  str(msebmtcal.cmprsk)
  names(attributes(msebmtcal.cmprsk))

  ### Define covariates for model
  covs <- c("match", "proph", "year", "agecl")
  msebmtcal.cmprsk <- mstate::expand.covs(msebmtcal.cmprsk, covs, longnames = FALSE)

  ### Assign variables for model we will be fitting
  eq.RHS <- paste(do.call(paste0, expand.grid(c("match", "proph", "year1", "year2", "agecl1", "agecl2"), paste(".", 1:sum(!is.na(tmat)), sep = ""))), collapse="+")
  strata <- survival::strata
  eq <- paste("survival::Surv(Tstart, Tstop, status) ~ ", eq.RHS,  "+ strata(trans)", sep = "")
  eq <- as.formula(eq)

  ### Create dataframe to store predicted risks
  tp.all <- data.frame(matrix(NA, ncol = 13, nrow = nrow(ebmtcal)))
  colnames(tp.all) <- c("id", paste("pstate", 1:6, sep = ""), paste("se", 1:6, sep = ""))

  ### Loop through id.iter
  for (id.iter in 1:5){

    print(paste("j = ", j, "id.iter = ", id.iter, Sys.time()))

    ### Develop a model on entire dataset except individual of interest
    cfull <- survival::coxph(eq, data = subset(msebmtcal.cmprsk, id != id.iter), method = "breslow")

    ### Get location of individual in msebmtcal.cmprsk
    pat.loc <- which(msebmtcal.cmprsk$id == id.iter)

    ### Create a miniture dataset, on which to generate predictions in (must be in mstate format and have a row for every transition)
    pat.dat <- msebmtcal.cmprsk[rep(pat.loc[1], sum(!is.na(tmat))), 9:12]
    pat.dat$trans <- 1:sum(!is.na(tmat))
    attr(pat.dat, "trans") <- tmat
    pat.dat <- mstate::expand.covs(pat.dat, covs, longnames = FALSE)
    pat.dat$strata <- pat.dat$trans

    ### Fit cause-specific hazards
    msf.pat <- mstate::msfit(cfull, pat.dat, trans = tmat)

    ### Write a function to extract the transition probabilities from state j into each state, after followup time f.time
    extract.tp <- function(tp.object, state, f.time){
      ### Create output object
      output.object <- return(base::subset(tp.object[[state]], time > f.time) |> dplyr::slice(1) |> dplyr::select(-c(time)))
    }

    for (r1.count in 1:length(r1.range)){
      r1 <- r1.range[r1.count]
      pt <- mstate::probtrans(msf.pat, predt = r1)
      for (r2.count in 1:length(r2.range)){
        r2 <- r2.range[r2.count]
        if (r2 > r1){
            tp_pred[[j]][[r1.count]][[r2.count]][id.iter, ] <- c(id.iter, extract.tp(tp.object = pt, state = j, f.time = r2 - r1), j, r1, r2)
        }
      }
    }

  }

}


str(pt[[1]])





### Run through each patient id and estimate transition probabilities
for (id.iter in 1:nrow(ebmt)){

  print(paste("id.iter = ", id.iter, Sys.time()))

  ### Develop a model on entire dataset except individual of interest (calculate the cause-specific hazards)
  cfull <- coxph(eq, data = subset(msebmt, id != id.iter), method = "breslow")

  ### Get location of individual in msebmt
  pat.loc <- which(msebmt$id == id.iter)

  ### Create a miniture dataset, on which to generate predictions in (must be in mstate format and have a row for every transition)
  pat.dat <- msebmt[rep(pat.loc[1], 12), 9:12]
  pat.dat$trans <- 1:12
  attr(pat.dat, "trans") <- tmat
  pat.dat <- expand.covs(pat.dat, covs, longnames = FALSE)
  pat.dat$strata <- pat.dat$trans

  ### Obtain cumulative incidence functions for the individual of interest
  msf.pat <- msfit(cfull, pat.dat, trans = tmat)

  ### Write a function to extract the transition probabilities from state j into each state, after followup time f.time
  extract.tp <- function(tp.object, state, f.time){
    ### Create output object
    output.object <- as.numeric(base::subset(tp.object[[state]], time > f.time) |> dplyr::slice(1) |> dplyr::select(-c(time)))
    return(output.object)
  }

  ### Calculate required transition probabilities and store in output dataset
  ### Will generate risks out of every state j and store in tp.id
  for (r1.count in 1:length(r1.range)){
    r1 <- r1.range[r1.count]
    pt <- probtrans(msf.pat, predt = r1)
      for (r2.count in 1:length(r2.range)){
        r2 <- r2.range[r2.count]
        if (r2 > r1){
          for (j in 1:4){
          tp_pred[[j]][[r1.count]][[r2.count]][id.iter, ] <- c(id.iter, extract.tp(tp.object = pt, state = j, f.time = r2 - r1), j, r1, r2)
        }
      }
    }
  }

}

length(tp_pred)
str(tp_pred[[4]])
### Combine into one data frame
for (j in 1:4){
  for (r1.count in 1:length(r1.range)){
    tp_pred[[j]][[r1.count]] <- do.call("rbind", tp_pred[[j]][[r1.count]])
  }
  tp_pred[[j]] <- do.call("rbind", tp_pred[[j]])
}
tp.r1.r2 <- do.call("rbind", tp_pred)

### Remove NA rows
tp.r1.r2 <- tp.r1.r2[complete.cases(tp.r1.r2),]

#################
### Create tp ###
#################
save.image("workflow/supp_diagnosis_miscalibration_prep_data_cmprsk_style_predrisks.RData")

# ### Use this data
# usethis::use_data(tp.r1.r2, overwrite = TRUE)

