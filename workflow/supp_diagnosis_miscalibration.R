#####################
### Preliminaries ###
#####################

### Clear workspace
rm(list=ls())

### Load data
library(mstate)
library(calibmsm)

### Load workspace
load("workflow/supp_diagnosis_miscalibration.RData")
tp.r1.r2 <- tp.r1.r2[complete.cases(tp.r1.r2), ]
table(tp.r1.r2$j)
### Change to following link when ready
#  load("workflow/supp_diagnosis_miscalibration_prep_data.RData")


#########################################################################
### Step 1: Define focus? Maybe this is a goal/aim rather than a step ###
#########################################################################

### We see miscalibration into state 5, both from starting state, and definitely at subsequent states from t = 100, so
### let's focus on calibration of probability of relapse (a key clinical outcome)

##########################################
### Step 2: Look at events frequencies ###
##########################################

events(msebmt)

### Can see lots of people go into recovery or AE states, so assessing calibration out of interest states will be very important
### Lots of people also go into the REC+AE state, so assessing calibration out of state 4 will be important too
### Only 25% of people enter relapse straight from the starting state.

### Lets get cross sections of how many people in each state at each time point
get_states_t <- function(t){
  sapply(lapply(c(1,2,3,4,5,6), extract_ids_states, data.mstate = msebmt, tmat = tmat, t = t), length)
}
lapply(c(0,100,round(c(1,2,3,4,5)*365.25)), get_states_t)

### Not much movement after yeras 1 and 2, so can focus on these times too


####################################################
### Step 3: Calibration out of the initial state ###
####################################################

### This is actually what I already did in supplemtary material comparison with graphical calibration curves
calc.calib.gcc.mod <- function(data.mstate, data.raw, j, s, t, p.est, nk = 3){

  data.mstate = msebmtcal
  data.raw = ebmtcal
  j = 1
  s = 0
  t = round(365.25*5)
  p.est = tp.r1.r2 |> subset(j == 1 & r1 == 0 & r2 == round(365.25*5)) |> dplyr::select(any_of(paste("pstate", 1:6, sep = "")))
  nk = 3

  ###
  ### Need to start by creating a new transition matrix, relevant to the cmprsk submodel of interest
  ### Will use mstate::transMat to do this
  ###

  ## Extract transition matrix for multistate data
  tmat.msm <- attributes(msebmtcal)$trans

  ## Create list of correct length for input into mstate::transMat
  list.obj <- vector("list", nrow(tmat.msm))

  ## For row j, add the possible transitions
  list.obj[[j]] <- as.numeric(which(!is.na(tmat.msm[j,])))

  ## Define tmat
  tmat <- mstate::transMat(x = list.obj, names = attributes(tmat.msm)$dimnames$from)

  ## Create data in msdata format
  msebmtcal.cmprsk <- mstate::msprep(data = ebmtcal, trans = tmat, time = c(NA, "rec", "ae","recae", "rel", "srv"),
                                     status = c(NA, "rec.s", "ae.s", "recae.s", "rel.s", "srv.s"),
                                     keep = c("match", "proph", "year", "agecl"))


  ###
  ### Initialisaiton of other input variables and datasets
  ###

  ### Assign colnames to p.est
  colnames(p.est) <- paste("p.est", 1:ncol(p.est), sep = "")

  ### Extract what states an individual can move into from state j (states with a non-zero predicted risk)
  ### Also drop the staet that an individual is already in, because there is no cmprsk model for staying in the same state
  valid.transitions <- which(colSums(p.est) != 0)
  valid.transitions <- valid.transitions[-(valid.transitions == j)]

  ### Add the predicted risks, and the complementary log log transformation of the predicted risks to data.raw
  p.est.cll <- log(-log(1 - p.est[,valid.transitions]))
  colnames(p.est.cll) <- paste("p.est.cll", valid.transitions, sep = "")

  ###
  ### Apply landmarking
  ###

  ### Identify individuals who are in state j at time s
  ids.state.j <- base::subset(data.mstate, from == j & Tstart <= s & s < Tstop) |>
    dplyr::select(id) |>
    dplyr::distinct(id) |>
    dplyr::pull(id)

  ### Subset data.mstate and data.raw to these individuals.
  data.mstate <- data.mstate |> base::subset(id %in% ids.state.j)
  data.raw <- data.raw |> base::subset(id %in% ids.state.j)

  ### Add the cloglog risks and predicted risks to the landmark dataset
  data.raw <- cbind(data.raw, p.est[,valid.transitions], p.est.cll)

  ### Finally, identify individuals which are censored before experiencing any events (used to maniuplate data for Fine-Gray regression later)
  ids.cens  <- data.mstate |> base::subset(from == j) |> dplyr::group_by(id) |> dplyr::summarize(sum = sum(status)) |> base::subset(sum == 0) |> dplyr::pull(id)

  ###
  ### Produce calibration plots for each possible transition
  ###

  ### Start by creating a list to store the plots
  plots.list <- vector("list", length(valid.transitions))

  for (k in 1:length(valid.transitions)){

    ### Assign state.k
    state.k <- as.numeric(valid.transitions[k])

    ### Create restricted cubic splines for the cloglog of the linear predictor for the state of interst
    rcs.mat <- Hmisc::rcspline.eval(data.raw[,paste("p.est.cll", state.k, sep = "")],nk=nk,inclx=T)
    colnames(rcs.mat) <- paste("rcs.x", 1:ncol(rcs.mat), sep = "")
    knots <- attr(rcs.mat,"knots")

    ### Create a new dataframe for the validation, to avoid recurison with data.raw
    ### Add the cubic splines for thecomplementary loglog of the predicted probability, and the predicted probability itself
    valid.df <- data.frame(data.raw$id, data.raw[,paste("p.est", state.k, sep = "")], rcs.mat)
    colnames(valid.df) <- c("id", "pred", colnames(rcs.mat))

    ### Want to validate the competing risks model out of state j at time s, into state k, so remove individuals not in state k at time s,
    ### and only retain transitions into state k. Also deduct immortal time from time variable
    data.mstate.j.k.s <- base::subset(data.mstate, from == j & to == state.k & Tstart <= s & s < Tstop) |>
      dplyr::mutate(time = Tstop - s) |>
      dplyr::select(c(time, status))

    ### Add to valid.df
    valid.df <- cbind(valid.df, data.mstate.j.k.s)

    ### For individuals who do not have the event of interest, and also are not censored (i.e. they have a different competing event),
    ### set the follow up time to the maximum
    valid.df <- dplyr::mutate(valid.df, time = dplyr::case_when(status == 0 & !(id %in% ids.cens) ~ max(time),
                                                                TRUE ~ time))

    ### Create dataset to fit the recalibration model of Austin et al (Graphical calibration curves, BMC Diagnostic and Prognostic, DOI10.1186/s41512-021-00114-6)
    valid.df.crprep <- mstate::crprep(Tstop='time',status='status',trans=1,
                                      keep=colnames(rcs.mat),valid.df)

    ### Create formula and fit the Fine-Gray recalibration model
    eq.LHS <- paste("survival::Surv(Tstart,Tstop,status==1)~")
    eq.RHS <- paste("rcs.x", 1:ncol(rcs.mat), sep = "", collapse = "+")
    eq <- formula(paste(eq.LHS, eq.RHS, sep = ""))
    model.calibrate.fg <- rms::cph(eq,weights=weight.cens,x=T,y=T,surv=T,data=valid.df.crprep)

    ### Generate predicted probabilities and standard errors
    valid.df$obs.fg <- 1-rms::survest(model.calibrate.fg,newdata=valid.df.crprep,time=t-s)$surv
    valid.df$obs.fg.upper<-1-rms::survest(model.calibrate.fg,newdata=valid.df.crprep,time=t-s)$lower
    valid.df$obs.fg.lower<-1-rms::survest(model.calibrate.fg,newdata=valid.df.crprep,time=t-s)$upper

    ### Produce plots for each and store in a list

    ### Pivot longer to create data for ggplot and assign appropriate labels
    valid.df.longer <- tidyr::pivot_longer(valid.df, cols = c(obs.fg, obs.fg.upper, obs.fg.lower), names_to = "line.group")
    valid.df.longer <- dplyr::mutate(valid.df.longer,
                                     line.group = factor(line.group),
                                     mapping = dplyr::case_when(line.group == "obs.fg" ~ 1,
                                                                line.group %in% c("obs.fg.upper", "obs.fg.lower") ~ 2),
                                     mapping = factor(mapping))

    levels(valid.df.longer$line.group) <- c("Calibration", "Upper", "Lower")
    levels(valid.df.longer$mapping) <- c("Calibration", "95% CI")

    ### Create the plot
    plots.list[[k]] <- ggplot2::ggplot(data = valid.df.longer |> dplyr::arrange(pred) |> dplyr::select(id, pred, line.group, value, mapping)) +
      ggplot2::geom_line(ggplot2::aes(x = pred, y = value, group = line.group, color = mapping)) +
      ggplot2::geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
      ggplot2::xlab("Predicted risk") + ggplot2::ylab("Observed risk") +
      ggplot2::xlim(c(0, max(valid.df.longer$pred,
                             valid.df.longer$value))) +
      ggplot2::ylim(c(0, max(valid.df.longer$pred,
                             valid.df.longer$value))) +
      ggplot2::geom_rug(data = valid.df.longer |> dplyr::arrange(pred) |> dplyr::select(id, pred, line.group, value, mapping) |> base::subset(line.group == "Calibration"),
                        ggplot2::aes(x = pred, y = value), col = grDevices::rgb(1, 0, 0, alpha = .1)) +
      ggplot2::theme(legend.position = "none") +
      ggplot2::ggtitle(paste("State ", state.k, sep = ""))

  }

  ### Return plots
  return(plots.list)
}


calib.j1.r1.0.r2.5 <- calc.calib.gcc.mod(data.mstate = msebmtcal,
                   data.raw = ebmtcal,
                   j = 1,
                   s = 0,
                   t = round(365.25*5),
                   p.est = tp.r1.r2 |> subset(j == 1 & r1 == 0 & r2 == round(365.25*5)) |> dplyr::select(any_of(paste("pstate", 1:6, sep = ""))),
                   nk = 3)

head(msebmt)

max(msebmt$time[msebmt$to == 3 & msebmt$status == 1])


