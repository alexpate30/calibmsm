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

###
### Add data
###
usethis::use_data(ebmtcal, overwrite = TRUE)
usethis::use_data(msebmtcal, overwrite = TRUE)
usethis::use_data(tps0, overwrite = TRUE)
usethis::use_data(tps100, overwrite = TRUE)

###
### Add R-CMD-CHECK
###
usethis::use_readme_rmd()
usethis::use_github_action_check_standard()
devtools::build_readme()

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
devtools::check()
library(dplyr)
getwd()
testchange
testchange2