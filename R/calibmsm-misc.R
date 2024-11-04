### Deal with global variables note
utils::globalVariables(c("id", "from", "to", "Tstart", "Tstop", "status", "state_poly_fac", "state_poly", #variables from data inputted into calc_calib_mlr
                         "ipcw", "ipcw_stab", #variable produced by weights function
                         "dtcens", "dtcens_s", #variables from data_raw used for applying landmarking
                         "value", "pred", "obs", "obs_upper", "obs_lower", "line_group", "mapping")) #variables from data inputted into plot_calib_blr
