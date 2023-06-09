### Package to do

### DONE 1) Add restricted cubic splines for when using pseudo-values (allows for easier calculation of confidence intervals)
### DONE 2) Add tests for pseudo-values
### NEED TO WRITE THIS VIGNETTE 3) Add a wrapper to help with confidence intervals for calc_calib_blr and calc_calib_mlr. Think I'm actually not going to do this, and allow people to specify their own weights function.
### DONE 4) Anything to do around sub-models? i.e. re-write the supp_gcc program
### DONE 5) Rewrite supplementary material (in particular think pseudo values ones should be removed?)
### DONE 6) Add more tests for pseudo-values when rcs type curve is included
### DONE 7) Add more warnings. A) Is data.mstate and data.raw got the same id's? B) Remove warning about weights and landmark.all = "all", otherwise it comes up when we specify land.mark.type = "all"
### 8) Change some terminology, A) remove the calc_ bits, B) change t.eval to t
### DONE 9) Add the supplementary material gcc to the paper
### DONE 10) Get rid of warnings in pv for when there are zero probability possible transitions
### DONE 11) Write vignette for producing confidence interval using your own method for estimating the weights
### DONE 12) Add a function to allow users to prespecify their weights function within calc_calib_pv. This can be part of the vignette.
### 13) Add the robust sandwich-type estimator
### DONE 14) Why is package failing CMD check now?
### DONE 15) Add something to text manuscript about robust sandwich being conservative
### DONE 16) Get r-cmd-check working, only possibility left is more troubleshooting of the error, or try removing the pdf vignette
### DONE 17) Check all vignettes refer to correct vignettes names of other vignettes
### 18) Update the example on the GitHub website (the Rmd). Also get auto update on the .Rmd file.
### 19) Change manuscript to remove the parametric bootstrap from the BLR-IPCW method, until it's actually programming in. Can't quite figure it out, and still a bit unsure about the theory.

### Notes to add to document
### 1) Verification of assumptions for each method
### 2) Whether to delve into which method I believe is appropriate
### 3) Need to verify Whether its true that LMAJ = Landmark then AJ

### It's basically just to add the robust sandwich...
