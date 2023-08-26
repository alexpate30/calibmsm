###
### Required packages for system setup
###
install.packages("devtools")
install.packages("roxygen2")
install.packages("knitr")
install.packages("testthat")

###
### Install
###
devtools::install()

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
usethis::use_package("mstate")
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
install.packages("R.rsp")
usethis::use_vignette("overview")
usethis::use_article("Comparison-with-graphical-calibration-curves-in-competing-risks-setting")
usethis::use_article("BLR-IPCW-calibration-curves-estimated-with-loess-smoothers")
devtools::install(build_vignettes = TRUE)
browseVignettes("calibmsm")

###
### Add file for creation of data
###
usethis::use_data_raw()
usethis::use_data(ebmtcal, overwrite = TRUE)
usethis::use_data(msebmtcal, overwrite = TRUE)
usethis::use_data(tps0, overwrite = TRUE)
usethis::use_data(tps100, overwrite = TRUE)

###
### Add a testing suite
###
usethis::use_testthat(3)

###
### Create readme and news
###
usethis::use_readme_rmd()
usethis::use_news_md()

###
### Add GitHub actions
###

### First link to GitHub
gh_token_help()
usethis::create_github_token()
gitcreds::gitcreds_set()

### R-CMD-check action
usethis::use_github_action_check_standard()

### Test coverage
## For test coverage, used following sources:
## 1) https://www.r-bloggers.com/2017/06/how-to-add-code-coverage-codecov-to-your-r-package/
## 2) The instructions on codecov.io
usethis::use_coverage(type = "codecov")
# covr::codecov(token = "INSERT TOKEN HERE")
use_github_action("test-coverage")

### Build readme
devtools::build_readme()

### Automatically render .Rmd
## This didn't work for me as the "renv" package wouldn't install when running action
## r-lib/actions/setup-renv@v2. I have removed for now.
## If reinstating, you need to remove the git pre-commit hooks which require .Rmd
## and .md files to be the same.
# usethis::use_github_action("render-rmarkdown")


###
### Create website
###

## Now create website
usethis::use_pkgdown()
pkgdown::build_site()
usethis::use_pkgdown_github_pages()

###
### Set ghostscript location
### If this causes issues with CRAN, just remove all .pdf vignettes
###
Sys.setenv(R_GSCMD = "C:\\Program Files\\gs\\gs10.01.1\\bin\\gswin64c.exe")

###
### Install package
###
devtools::install()

###
### Add CRAN comments
###
usethis::use_cran_comments()

###
### Run R-CMD-CHECK to be checked with CRANs win-builder service. Check email response for this.
###
devtools::check()
devtools::check(remote = TRUE, manual = TRUE)
devtools::check_win_devel()

###
### Update version
###
usethis::use_release_issue()
usethis::use_version('major')

