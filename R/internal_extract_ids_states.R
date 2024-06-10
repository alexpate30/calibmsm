#' Identify patids for individuals in state j at time t
#'
#' @description
#' Extract patids for individuals in state j at time t from a dataset in 'msdata'
#' format. Used internally in calib_blr, calib_mlr and calib_pv.
#'
#' @param data.ms Validation data in `msdata` format
#' @param tmat Transition probability matrix
#' @param j State j
#' @param t Follow up time
#'
#' @noRd
extract_ids_states <- function(data.ms, tmat, j, t){

  ### Define maximum state number
  max.state <- max(data.ms$to)

  ### Identify which states are absorbing states
  absorbing.states <- which(apply(tmat, 1, function(x) {sum(!is.na(x))}) == 0)

  ### For non-absorbing states, to be in state j at time t, you must have an observations from state j, where Tstart <= t < Tstop
  if (!(j %in% absorbing.states)){
    ## Extract ids
    ids.state.j <- base::subset(data.ms, from == j & Tstart <= t & t < Tstop) |>
      dplyr::select(id) |>
      dplyr::distinct(id)
    ## Put into numeric vector
    ids.state.j <- as.numeric(ids.state.j$id)
  } else if (j %in% absorbing.states){
    ### For absorbing state, just have to have moved into it
    ids.state.j <- base::subset(data.ms, to == j & t >= Tstop & status == 1) |>
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
#' `exclude.cens.t = FALSE` (default) will allow individuals censored at time `t` to remain
#' in the cohort (required for pseudo-value approach).
#' Choosing `exclude.cens.t = TRUE` will also remove individuals who are censored by time `t`.
#'
#' Note that the `exclude.cens.t` argument will not work for `data.return = "data.ms"`.
#' @noRd
apply_landmark <- function(data.raw, data.ms, j, s, t, exclude.cens.t = FALSE, data.return = "data.raw"){

  ### Extract transition matrix from msdata object
  tmat <- attributes(data.ms)$trans

  ## Identify individuals in state j at time s
  ids.state.js <- extract_ids_states(data.ms = data.ms, tmat = tmat, j = j, t = s)

  ## Create landmarked dataset
  if (data.return == "data.raw"){
    if (exclude.cens.t == FALSE){
      data.lmk.js <-  data.raw |> base::subset(id %in% ids.state.js)
    } else if (exclude.cens.t == TRUE){
      data.lmk.js <-  data.raw |> base::subset(id %in% ids.state.js) |> base::subset(dtcens > t | dtcens <= t & dtcens.s == 0)
    }
  } else if (data.return == "data.ms"){
    data.lmk.js <-  data.ms |> base::subset(id %in% ids.state.js)
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
identify_valid_transitions <- function(data.raw, data.ms, j, s, t){

  ### Landmark the dataset, retaining only individuals uncensored at time t
  data.raw.lmk.js <- apply_landmark(data.raw = data.raw, data.ms, j = j, s = s, t = t, exclude.cens.t = TRUE)

  ### Assign the maximum state an individual may enter
  max.state <- max(data.ms$to)

  ### Extract transition matrix from msdata object
  tmat <- attributes(data.ms)$trans

  ### Identify states these individuals are in at time t
  ids.state.list <- vector("list", max.state)
  for (k in 1:max.state){
    ids.state.list[[k]] <- extract_ids_states(data.ms, tmat, k, t)
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


#' Apply bootstrapping to a dataset of class `msdata`
#'
#' @description
#' Apply bootstrapping to datasets of class `msdata` (i.e. `data.ms`). This is non-trivial
#' because there is more than one row per individual, and a new `id` variable `id2` must be
#' assigned.
#'
#' @noRd
apply_bootstrap_msdata <- function(data.ms, indices){

  ### Break up data.ms by id
  data.ms.list <- split(data.ms, data.ms$id)

  ### Extract the relevant list elements based on id
  data.ms.boot <- lapply(1:length(indices),
                             function(x) {
                               data.frame(data.ms.list[[indices[x]]], "id2" = x)
                             }
  )
  names(data.ms.boot) <- names(data.ms.list)

  ### Combine into a single dataset and give appropriate class
  data.ms.boot <- do.call("rbind", data.ms.boot)
  rownames(data.ms.boot) <- NULL
  class(data.ms.boot) <- c("msdata", "data.frame")

  ### Return
  return(data.ms.boot)

}
