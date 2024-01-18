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
#'
#' @noRd
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


#' Apply landmarking
#'
#' @description
#' Reduce cohort to individual who are uncensored and in state `j` at time `s`. Choosing
#' `exclude.cens.t = FALSE` (default) will also remove individuals who are censored by time `t`.
#' Choosing `exclude.cens.t = TRUE` will allow individuals censored at time `t` to remain
#' in the cohort (required for pseudo-value approach). This step applied through the variable
#' 'state.poly' which will be NA if an individual is censored at time `t`.
#'
#' @noRd
apply_landmark <- function(data.raw, data.mstate, j, s, t, exclude.cens.t = FALSE, data.return = "data.raw"){

  ### Extract transition matrix from msdata object
  tmat <- attributes(data.mstate)$trans

  ## Identify individuals in state j at time s
  ids.state.js <- extract_ids_states(data.mstate = data.mstate, tmat = tmat, j = j, t = s)

  ## Create landmarked dataset
  if (data.return == "data.raw"){
    if (exclude.cens.t == FALSE){
      data.lmk.js <-  data.raw |> base::subset(id %in% ids.state.js)
    } else if (exclude.cens.t == TRUE){
      data.lmk.js <-  data.raw |> base::subset(id %in% ids.state.js) |> base::subset(dtcens > t | dtcens <= t & dtcens.s == 0)
    }
  } else if (data.return == "data.mstate"){
    data.lmk.js <-  data.mstate |> base::subset(id %in% ids.state.js)
  }

  return(data.lmk.js)

}


#' Identify valid transitions
#'
#' @description
#' Identify states which can be entered when in state j at time s
#'
#' Note that this returns which states have people in at time t, amongst those that
#' were in state j at time s.
#'
#' @noRd
identify_valid_transitions <- function(data.raw, data.mstate, j, s, t){

  ### Landmark the dataset, retaining only individuals uncensored at time t
  data.raw.lmk.js <- apply_landmark(data.raw = data.raw, data.mstate, j = j, s = s, t = t, exclude.cens.t = TRUE)

  ### Assign the maximum state an individual may enter
  max.state <- max(data.mstate$to)

  ### Extract transition matrix from msdata object
  tmat <- attributes(data.mstate)$trans

  ### Identify states these individuals are in at time t
  ids.state.list <- vector("list", max.state)
  for (k in 1:max.state){
    ids.state.list[[k]] <- extract_ids_states(data.mstate, tmat, k, t)
  }

  ### Create a variable to say which state an individual was in at the time of interest
  ## Create list containing the relevant data
  v1 <- data.raw.lmk.js$id
  m1 <- outer(v1, ids.state.list, FUN = Vectorize('%in%'))
  state.poly <- lapply(split(m1, row(m1)), function(x) (1:max.state)[x])

  ## Change integer(0) values to NA's
  idx <- !sapply(state.poly, length)
  state.poly[idx] <- NA

  ## Add to data.raw
  valid.transitions <- as.numeric(names(table(unlist(state.poly))))

  return(valid.transitions)

}
