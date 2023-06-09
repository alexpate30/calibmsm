#' European Group for Blood and Marrow Transplantation data
#'
#' ebmt data in format ready to assess calibration
#'
#' @format ## 'ebmtcal'
#' A data frame with 2,279 rows and 17 columns:
#' \describe{
#'   \item{id}{Patient indentifier}
#'   \item{rec, rec.s}{Time until and event indicator for recovery variable}
#'   \item{ae, ae.s}{Time until and event indicator for adverse event variable}
#'   \item{recae, recae.s}{Time until and event indicator for recovery + adverse event variable}
#'   \item{rel, rel.s}{Time until and event indicator for relapse variable}
#'   \item{srv, srv.s}{Time until and event indicator for death variable}
#'   \item{year}{Year of transplant}
#'   \item{agecl}{Age at transplant}
#'   \item{proph}{Prophylaxis}
#'   \item{match}{Donor-recipient match}
#'   \item{dtcens}{Time of censoring}
#'   \item{dtcens.s}{Event indicator, 1:censoring occured, 0: absorbing state entered before censoring occured}
#' }
#' @source <https://search.r-project.org/CRAN/refmans/mstate/html/EBMT-data.htmlE>
"ebmtcal"

#' European Group for Blood and Marrow Transplantation data
#'
#' ebmt data in 'msdata' format, derived from the ebmt dataset using software from the mstate package
#'
#' @format ## 'msebmtcal'
#' A data frame with 15,512 rows and 8 columns:
#' \describe{
#'   \item{id}{Patient indentifier}
#'   \item{from}{transition from state}
#'   \item{to}{transition to state}
#'   \item{trans}{transition number}
#'   \item{Tstart}{time entered state 'from'}
#'   \item{Tstop}{time leaving state 'from'}
#'   \item{time}{time in state 'from'}
#'   \item{status}{event indicator, 1:transitioned to state 'to'}
#' }
"msebmtcal"

#' Predicted transition probabilities out of every state at time s = 0
#'
#'
#' @format ## 'tps0'
#' A data frame with 13,674 rows and 14 columns:
#' \describe{
#'   \item{id}{Patient indentifier}
#'   \item{pstate1, pstate2, pstate3, pstate4, pstate5, pstate6}{Predicted transition probabilities of transitions into states 1 to 6}
#'   \item{se1, se2, se3, se4, se5, se6}{Standard error of the predicted transition probabilities of transitions into states 1 to 6}
#'   \item{j}{State from which the predicted transition probabilities are estimated from}
#' }
"tps0"

#' Predicted transition probabilities out of every state at time s = 100
#'
#'
#' @format ## 'tps100'
#' A data frame with 13,674 rows and 14 columns:
#' \describe{
#'   \item{id}{Patient indentifier}
#'   \item{pstate1, pstate2, pstate3, pstate4, pstate5, pstate6}{Predicted transition probabilities of transitions into states 1 to 6}
#'   \item{se1, se2, se3, se4, se5, se6}{Standard error of the predicted transition probabilities of transitions into states 1 to 6}
#'   \item{j}{State from which the predicted transition probabilities are estimated from}
#' }
"tps100"

#' Predicted risks for a competing risks model out of state j = 0
#'
#' Used in vignette: Comparison-with-graphical-calibration-curves-in-competing-risks-setting. The predicted transition probabilities
#' are made at out of state j = 1 at time s = 0, and treat all states as absorbing states.
#'
#' @format ## 'tp.cmprsk.j0'
#' A data frame with 2,279 rows and 13 columns:
#' \describe{
#'   \item{id}{Patient indentifier}
#'   \item{pstate1, pstate2, pstate3, pstate4, pstate5, pstate6}{Predicted transition probabilities of transitions into states 1 to 6}
#'   \item{se1, se2, se3, se4, se5, se6}{Standard error of the predicted transition probabilities of transitions into states 1 to 6}
#' }
"tp.cmprsk.j0"

#' European Group for Blood and Marrow Transplantation data in competing risks format, for transitions out of the initial state only
#'
#' Used in vignette: Comparison-with-graphical-calibration-curves-in-competing-risks-setting.
#' This is the ebmt data in 'msdata' format, derived from the ebmt dataset using software from the mstate package.
#' However, only transitions out of the initial state are considered, making this data for competing risks model out of the initial state.
#'
#' @format ## 'msebmtcal.cmprsk'
#' A data frame with 9,116 rows and 8 columns:
#' \describe{
#'   \item{id}{Patient indentifier}
#'   \item{from}{transition from state}
#'   \item{to}{transition to state}
#'   \item{trans}{transition number}
#'   \item{Tstart}{time entered state 'from'}
#'   \item{Tstop}{time leaving state 'from'}
#'   \item{time}{time in state 'from'}
#'   \item{status}{event indicator, 1:transitioned to state 'to'}
#' }
"msebmtcal.cmprsk"
