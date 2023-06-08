###
### Supplementary material 2 gcc
### The aim of this analysis is assess calibration of the sub-models of the multistate model using our framework and graphical calibration curves
### We will focus on the sub-model out of state j = 1 at time s = 0.
###

### Clear workspace
rm(list=ls())

### Load calibmsm and required data
library("calibmsm")

####################################################################
### Create transition probabilities using leave one out approach ###
####################################################################

### Define state which sub-model is coming out of and landmark time
j <- 1
s <- 0

### Define t.eval
t.eval <- 1826

### Define tmat
tmat <- mstate::transMat(x = list(c(2, 3, 5, 6), c(), c(), c(),
                c(), c()), names = c("Tx", "Rec", "AE", "Rec+AE", "Rel", "Death"))

### Load the data
data.raw.cmprsk <- ebmtcal

### Create mstate format
data.mstate.cmprsk <- mstate::msprep(data = data.raw.cmprsk, trans = tmat, time = c(NA, "rec", "ae","recae", "rel", "srv"),
                 status = c(NA, "rec.s", "ae.s", "recae.s", "rel.s", "srv.s"),
                 keep = c("match", "proph", "year", "agecl"))

### Define covariates for model
covs <- c("match", "proph", "year", "agecl")
data.mstate.cmprsk <- mstate::expand.covs(data.mstate.cmprsk, covs, longnames = FALSE)

### Create a variable which is maximum observed follow up time for all individuals, this is when they were either censored, relapsed or died
data.raw.cmprsk$dtcens <- pmin(data.raw.cmprsk$rel, data.raw.cmprsk$srv)
data.raw.cmprsk$dtcens.s <- 1 - pmax(data.raw.cmprsk$rel.s, data.raw.cmprsk$srv.s)

### Assign variables for model we will be fitting
eq.RHS <- paste(do.call(paste0, expand.grid(c("match", "proph", "year1", "year2", "agecl1", "agecl2"), paste(".", 1:sum(!is.na(tmat)), sep = ""))), collapse="+")
strata <- survival::strata
eq <- paste("survival::Surv(Tstart, Tstop, status) ~ ", eq.RHS,  "+ strata(trans)", sep = "")
eq <- as.formula(eq)

### Create dataframe to store predicted risks
tp.all <- data.frame(matrix(NA, ncol = 13, nrow = nrow(data.raw.cmprsk)))
colnames(tp.all) <- c("id", paste("pstate", 1:6, sep = ""), paste("se", 1:6, sep = ""))

### Loop through id.iter
for (id.iter in 1:nrow(data.raw.cmprsk)){

  print(paste("id.iter = ", id.iter, Sys.time()))

  ### Develop a model on entire dataset except individual of interest
  cfull <- survival::coxph(eq, data = subset(data.mstate.cmprsk, id != id.iter), method = "breslow")

  ### Get location of individual in data.mstate.cmprsk
  pat.loc <- which(data.mstate.cmprsk$id == id.iter)

  ### Create a miniture dataset, on which to generate predictions in (must be in mstate format and have a row for every transition)
  pat.dat <- data.mstate.cmprsk[rep(pat.loc[1], sum(!is.na(tmat))), 9:12]
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


################################################################
### Create a calibration plot using calibmsm::calc_calib_blr ###
################################################################

dat.calib.blr <-
  calc_calib_blr(data.mstate = data.mstate.cmprsk,
                 data.raw = data.raw.cmprsk,
                 j=1,
                 s=0,
                 t.eval = t.eval,
                 tp.pred = tp.all |>
                   dplyr::select(any_of(paste("pstate", 1:6, sep = ""))),
                 curve.type = "rcs",
                 rcs.nk = 3,
                 w.covs = c("year", "agecl", "proph", "match"),
                 CI = 95,
                 CI.R.boot = 200)


plot.calibmsm.blr.rcs <- plot(dat.calib.blr, combine = TRUE, nrow = 2, ncol = 3)

Cairo::CairoPNG(paste("workflow/figures/supp_gcc_j1s0_calibmsm_blr_rcs.png", sep = ""),
                dpi = 300, width = 15, height = 10, unit = "in")
print(plot.calibmsm.blr.rcs)
dev.off()

print(paste("BLR-IPCW DONE", Sys.time()))

################################################################
### Create a calibration plot using calibmsm::calc_calib_pv ###
################################################################

dat.calib.pv <-
  calc_calib_pv(data.mstate = data.mstate.cmprsk,
                 data.raw = data.raw.cmprsk,
                 j=1,
                 s=0,
                 t.eval = t.eval,
                 tp.pred = tp.all |>
                   dplyr::select(any_of(paste("pstate", 1:6, sep = ""))),
                 curve.type = "rcs",
                 rcs.nk = 3,
                 n.pctls = 9,
                 CI = 95,
                 CI.type = "parametric")


plot.calibmsm.pv.rcs <- plot(dat.calib.pv, combine = TRUE, nrow = 2, ncol = 3)

Cairo::CairoPNG(paste("workflow/figures/supp_gcc_j1s0_calibmsm_pv_rcs.png", sep = ""),
                dpi = 300, width = 15, height = 10, unit = "in")
print(plot.calibmsm.pv.rcs)
dev.off()

print(paste("PV DONE", Sys.time()))

###########################################
### Create a calibration plot using gcc ###
###########################################
calc.calib.gcc.mod <- function(data.mstate, data.raw, j, s, t.eval, p.est, nk = 3){

        # data.mstate <- data.mstate.cmprsk
        # data.raw <- data.raw.cmprsk
        # p.est <- tp.all[,paste("pstate", 1:6, sep = "")]
        # j <- 1
        # s <- 0
        # t.eval <- 1826
        # nk <- 3

  ### Assign colnames to p.est
  colnames(p.est) <- paste("p.est", 1:ncol(p.est), sep = "")

  ### Extract what states an individual can move into from state j (states with a non-zero predicted risk)
  ### Also drop the staet that an individual is already in, because there is no cmprsk model for staying in the same state
  valid.transitions <- which(colSums(p.est) != 0)
  valid.transitions <- valid.transitions[-(valid.transitions == j)]

  ### Add the predicted risks, and the complementary log log transormation of the predicted risks to data.raw
  p.est.cll <- log(-log(1 - p.est[,valid.transitions]))
  colnames(p.est.cll) <- paste("p.est.cll", valid.transitions, sep = "")

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
    valid.df$obs.fg <- 1-rms::survest(model.calibrate.fg,newdata=valid.df.crprep,time=t.eval-s)$surv
    valid.df$obs.fg.upper<-1-rms::survest(model.calibrate.fg,newdata=valid.df.crprep,time=t.eval-s)$lower
    valid.df$obs.fg.lower<-1-rms::survest(model.calibrate.fg,newdata=valid.df.crprep,time=t.eval-s)$upper

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

  ### DON'T THINK THIS WORKS
  # ### Finally a plot for starting state
  # ### The probability of staying in the starting state is 1 minus the sum of the probability of all other states
  # ### The observed event rate of staying in the starting state is 1 minus the sum of the observed event rates of all other states
  # value.init <- 1-rowSums(do.call("cbind", lapply(plots.list, function(x) {x$data |> dplyr::arrange(id, line.group) |> dplyr::pull(value)})))
  # pred.init <- 1-rowSums(do.call("cbind", lapply(plots.list, function(x) {x$data |> dplyr::arrange(id, line.group) |> dplyr::pull(pred)})))
  #
  # valid.df.init <- data.frame("id" = plots.list[[1]]$data$id,
  #                             "pred" = pred.init,
  #                             "value" = value.init,
  #                             "line.group" = plots.list[[1]]$data$line.group,
  #                             "mapping" = plots.list[[1]]$data$mapping)
  #
  # plot.init <- ggplot2::ggplot(data = valid.df.init |> dplyr::arrange(pred) |> base::subset(line.group == "Calibration") |> dplyr::select(id, pred, line.group, value, mapping)) +
  #   ggplot2::geom_line(ggplot2::aes(x = pred, y = value, group = line.group, color = mapping)) +
  #   ggplot2::geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  #   ggplot2::xlab("Predicted risk") + ggplot2::ylab("Observed risk") +
  #   ggplot2::xlim(c(0, max(valid.df.init$pred,
  #                          valid.df.init$value))) +
  #   ggplot2::ylim(c(0, max(valid.df.init$pred,
  #                          valid.df.init$value))) +
  #   ggplot2::geom_rug(data = valid.df.init |> dplyr::arrange(pred) |> dplyr::select(id, pred, line.group, value, mapping) |> base::subset(line.group == "Calibration"),
  #                     ggplot2::aes(x = pred, y = value), col = grDevices::rgb(1, 0, 0, alpha = .1)) +
  #   ggplot2::theme(legend.position = "none") +
  #   ggplot2::ggtitle(paste("State ", state.k, sep = ""))
  #
  # ### Add to plot list
  # plot.init
  # plots.list <- c(list(plot.init), plots.list)

  ### Return plots
  return(plots.list)
}

### Create plots
plot.gcc.rcs.list <- calc.calib.gcc.mod(data.mstate = data.mstate.cmprsk,
                                    data.raw = data.raw.cmprsk,
                                    j = 1,
                                    s = 0,
                                    t.eval = t.eval,
                                    p.est = tp.all[,paste("pstate", 1:6, sep = "")],
                                    nk = 3)

### Combine into one plot
plot.gcc.rcs <- ggpubr::ggarrange(plotlist = plot.gcc.rcs.list)

### Save image
Cairo::CairoPNG(paste("workflow/figures/supp_gcc_j1s0_gcc.png", sep = ""),
                dpi = 300, width = 15, height = 10, unit = "in")
print(plot.gcc.rcs)
dev.off()

save.image("workflow/supp_gcc.RData")

