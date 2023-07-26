#' Identify patids for individuals in state j at time t
#'
#' @description
#' Extract patids for individuals in state j at time t from a dataset in 'msdata'
#' format. Used internally in calib_blr, calib_mlr and calib_pv.
#'
#' @param data.mstate Validation data in `msdata` format
#' @param tmat Transition probability matrix
#' @param j State j
#' @param t Follow up time
extract_ids_states <- function(data.mstate, tmat, j, t){

  ### Define maximum state number
  max.state <- max(data.mstate$to)

  ### Identify which states are absorbing states
  absorbing.states <- which(apply(tmat, 1, function(x) {sum(!is.na(x))}) == 0)

  ### For non-absorbing states, to be in state j at time t, you must have an observations from state j, where Tstart <= t < Tstop
  if (!(j %in% absorbing.states)){
    ## Extract ids
    ids.state.j <- base::subset(data.mstate, from == j & Tstart <= t & t < Tstop) |>
      dplyr::select(id) |>
      dplyr::distinct(id)
    ## Put into numeric vector
    ids.state.j <- as.numeric(ids.state.j$id)
  } else if (j %in% absorbing.states){
    ### For absorbing state, just have to have moved into it
    ids.state.j <- base::subset(data.mstate, to == j & t >= Tstop & status == 1) |>
      dplyr::select(id) |>
      dplyr::distinct(id)
    ## Put into numeric vector
    ids.state.j <- as.numeric(ids.state.j$id)
  }

  return(ids.state.j)
}
