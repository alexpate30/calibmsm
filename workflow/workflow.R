###
### Update .Rd files
###
devtools::document()

###
### Add license
###
usethis::use_mit_license()

###
### Use git
###
use_git()

###
### Add dependencies
###
usethis::use_package("survival")
#usethis::use_package("mstate")
usethis::use_package("dplyr")
usethis::use_package("tidyr")
usethis::use_package("stats")
usethis::use_package("ggplot2")
usethis::use_package("ggpubr")
#usethis::use_package("utils")
usethis::use_package("Hmisc")
usethis::use_package("rms")
usethis::use_package("boot")
usethis::use_package("VGAM")
#usethis::use_testthat(3)

###
### Add functions to refer to without ::
###
usethis::use_import_from("magrittr", "%>%")
usethis::use_import_from("stats", "predict")
usethis::use_import_from("VGAM", "sm.ps")
usethis::use_import_from("VGAM", "sm.os")
usethis::use_import_from("VGAM", "s")

###
### Create first vignette
###
usethis::use_vignette("overview")
usethis::use_vignette("vigtest")
devtools::install(build_vignettes = TRUE)
browseVignettes()

###
### Add data
###
usethis::use_data(ebmtcal, overwrite = TRUE)
usethis::use_data(msebmtcal, overwrite = TRUE)
usethis::use_data(tps0, overwrite = TRUE)
usethis::use_data(tps100, overwrite = TRUE)

###
### Add a testing suite
###
usethis::use_testthat(3)

###
### GitHub actions
###

### R-CMD-check
usethis::use_readme_rmd()
usethis::use_github_action_check_standard()

### Test coverage
## For test coverage, used following sources:
## 1) https://www.r-bloggers.com/2017/06/how-to-add-code-coverage-codecov-to-your-r-package/
## 2) The instructions on codecov.io
usethis::use_coverage(type = "codecov")
# covr::codecov(token = "INSERT TOKEN HERE")
use_github_action("test-coverage")
devtools::build_readme()

### Automatically render .Rmd
usethis::use_github_action("render-rmarkdown")

###
### Create website
###
library(dplyr)
usethis::use_pkgdown()
pkgdown::build_site()
usethis::use_pkgdown_github_pages()
gh_token_help()
usethis::create_github_token()
gitcreds::gitcreds_set()
getwd()
hello <- 2
hello <- 3
library(gitcreds)
.libPaths()
whyisntworkflowbeingrecognised <- 2
whyisntworkflowbeingrecognised <- 2

###
### Install package
###
getwd()
load_all()
calc_calib_blr
devtools::document()
devtools::check(vignettes = FALSE)
calibmsm::calc_calib_blr
devtools::install()
?calibmsm::calc_weights()
?calibmsm::calc_weights()
?calibmsm::calc_calib_blr()
devtools::check(build_args = "--compact-vignettes=gs+qpdf")
devtools::check()
install.packages("Formula")
library(dplyr)
getwd()
testchange
testchange2
.libPaths("C:/Program Files/R/R-4.3.0/library")
.libPaths()
install.packages("pkgdown")
install.packages("devtools")
install.packages("gitcreds")
