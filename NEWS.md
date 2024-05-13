# calibmsm 1.1.0

* Change to package structure.

* Introduced one master function `calibmsm` which assess calibration using BLR-IPCW, MLR-IPCW and pseudo-values methodology.

* `calib_blr` changed to `calib_blr_ipcw`, `calib_mlr` changed to `calib_mlr_ipcw`, and `calib_blr_ipcw`, `calib_mlr_ipcw` and `calib_pv` are all changed into internal functions.

* Return the pseudo-values used to estimate calibration curves in the output for `calib_pv`, this allows users to parallelise code when estimating pseudo-values.

* Some minor changes to code of internal functions, and introduction of some new internal functions for processes which occur a number of times within the package (for example application of landmarking).

# calibmsm 1.0.0

* Submission of first version of package to CRAN and Journal of Statistical Software.

# calibmsm 0.0.0.9000

* Added a `NEWS.md` file to track changes to the package.

# calibmsm 1.0.0.0000

* Submission of first version of package to CRAN. Please see *Overview* vignette for detailed explanation of the package.
