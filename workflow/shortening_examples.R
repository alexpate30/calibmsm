library(calibmsm)

data("tps0")
data("ebmtcal")
data("msebmtcal")

# Extract the predicted transition probabilities out of state j = 1
tp.pred <- dplyr::filter(tps0, id %in% 1:50) |> dplyr::filter(j == 1) |>
  dplyr::select(any_of(paste("pstate", 1:6, sep = "")))
ebmtcal <- ebmtcal |> dplyr::filter(id %in% 1:50)
msebmtcal <- msebmtcal |> dplyr::filter(id %in% 1:50)

# Now estimate the observed event probabilities for each possible transition.
time.in <- Sys.time()
dat.calib.mlr <-
calib_mlr(data.mstate = msebmtcal,
 data.raw = ebmtcal,
 j=1,
 s=0,
 t = 1826,
 tp.pred = tp.pred,
 w.covs = c("year", "agecl", "proph", "match"))
time.out <- Sys.time()
time.out - time.in

# The data for each calibration scatter plots are stored in the "plotdata"
# list element.
str(dat.calib.mlr)



# Extract the predicted transition probabilities out of state j = 1
tp.pred <- dplyr::filter(tps0, id %in% 1:50) |> dplyr::filter(j == 1) |>
  dplyr::select(any_of(paste("pstate", 1:6, sep = "")))
ebmtcal <- ebmtcal |> dplyr::filter(id %in% 1:50)
msebmtcal <- msebmtcal |> dplyr::filter(id %in% 1:50)

# Now estimate the observed event probabilities for each possible transition.
dat.calib.pseudo <- calib_pv(data.mstate = msebmtcal,
  data.raw = ebmtcal,
  j = 1,
  s = 0,
  t = 1826,
  tp.pred = tp.pred)
time.out <- Sys.time()
time.out - time.in

# The data for each calibration curve are stored in the "plotdata" list
# element.
str(dat.calib.pseudo)



###
### cmprsk style
###
time.in <- Sys.time()

data("tp.cmprsk.j0")
data("ebmtcal")
data("msebmtcal.cmprsk")

# Extract the predicted transition probabilities out of state j = 1 for first 150 individuals
tp.pred <- tp.cmprsk.j0 |>
  dplyr::filter(id %in% 1:150) |>
  dplyr::select(any_of(paste("pstate", 1:6, sep = "")))
# Reduce ebmtcal to first 150 individuals
ebmtcal <- ebmtcal |> dplyr::filter(id %in% 1:150)
# Reduce msebmtcal.cmprsk to first 150 individuals
msebmtcal.cmprsk <- msebmtcal.cmprsk |> dplyr::filter(id %in% 1:150)

# Now estimate the observed event probabilities for each possible transition.

dat.calib.mlr <-
  calib_mlr(data.mstate = msebmtcal.cmprsk,
            data.raw = ebmtcal,
            j=1,
            s=0,
            t = 1826,
            tp.pred = tp.pred,
            w.covs = c("year", "agecl", "proph", "match"),
            ps.int = 2,
            degree = 2)


# The data for each calibration scatter plots are stored in the "plotdata"
# list element.
str(dat.calib.mlr)

time.out <- Sys.time()
time.out - time.in

###
### cmprsk example
###
time.in <- Sys.time()

data("tp.cmprsk.j0")
data("ebmtcal")
data("msebmtcal.cmprsk")


# Extract the predicted transition probabilities out of state j = 1 for first 75 individuals
tp.pred <- tp.cmprsk.j0 |>
  dplyr::filter(id %in% 1:50) |>
  dplyr::select(any_of(paste("pstate", 1:6, sep = "")))
# Reduce ebmtcal to first 150 individuals
ebmtcal <- ebmtcal |> dplyr::filter(id %in% 1:50)
# Reduce msebmtcal.cmprsk to first 150 individuals
msebmtcal.cmprsk <- msebmtcal.cmprsk |> dplyr::filter(id %in% 1:50)

# Now estimate the observed event probabilities for each possible transition.
dat.calib.pv <-
  calib_pv(data.mstate = msebmtcal.cmprsk,
            data.raw = ebmtcal,
            j=1,
            s=0,
            t = 1826,
            tp.pred = tp.pred,
            curve.type = "loess",
            loess.span = 1,
            loess.degree = 1)

dat.calib.pv[[1]][[1]]
devtools::load_all()
plot(dat.calib.pv, combine = FALSE)[[1]]
plot(dat.calib.pv, combine = FALSE)[[1]]
# The data for each calibration scatter plots are stored in the "plotdata"
# list element.

time.out <- Sys.time()
time.out - time.in
warnings()

#' dat.calib.pseudo <- calib_pv(data.mstate = msebmtcal,
#'   data.raw = ebmtcal,
#'   j = 3,
#'   s = 100,
#'   t = 1826,
#'   tp.pred = tp.pred,
#'   group.vars = c("year"),
#'   n.pctls = 2)
#'
#' # The data for each calibration curve are stored in the "plotdata" list
#' # element.
#' str(dat.calib.pseudo)
