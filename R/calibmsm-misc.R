### Deal with global variables note
utils::globalVariables(c("id", "from", "to", "Tstart", "Tstop", "status", "state.poly.fac", "state.poly", #variables from data inputted into calc_calib_mlr
                         "ipcw", "ipcw.stab", #variable produced by weights function
                         "dtcens", "dtcens.s", #variables from data.raw used for applying landmarking
                         "value", "pred", "obs", "obs.upper", "obs.lower", "line.group", "mapping")) #variables from data inputted into plot_calib_blr
