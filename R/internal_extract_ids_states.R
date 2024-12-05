#' Identify patids for individuals in state j at time t
#'
#' @description
#' Extract patids for individuals in state j at time t from a dataset in 'msdata'
#' format. Used internally in calib_blr, calib_mlr and calib_pv.
#'
#' @param data_ms Validation data in `msdata` format
#' @param tmat Transition probability matrix
#' @param j State j
#' @param t Follow up time
#'
#' @noRd
extract_ids_states <- function(data_ms, tmat, j, t){

  ### Define maximum state number
  max_state <- max(data_ms$to)

  ### Identify which states are absorbing states
  absorbing_states <- which(apply(tmat, 1, function(x) {sum(!is.na(x))}) == 0)

  ### For non-absorbing states, to be in state j at time t, you must have an observations from state j, where Tstart <= t < Tstop
  if (!(j %in% absorbing_states)){
    ## Extract ids
    ids_state_j <- base::subset(data_ms, from == j & Tstart <= t & t < Tstop) |>
      dplyr::select(id) |>
      dplyr::distinct(id)

  } else if (j %in% absorbing_states){
    ### For absorbing state, just have to have moved into it
    ids_state_j <- base::subset(data_ms, to == j & t >= Tstop & status == 1) |>
      dplyr::select(id) |>
      dplyr::distinct(id)

  }

  ## Put into numeric vector
  if (is.factor(data_ms$id) == FALSE){
    ids_state_j <- as.numeric(ids_state_j$id)
  } else if (is.factor(data_ms$id) == TRUE){
    ids_state_j <- as.numeric(as.character(ids_state_j$id))
  }

  return(ids_state_j)

}


#' Apply landmarking
#'
#' @description
#' Reduce cohort to individual who are uncensored and in state `j` at time `s`. Choosing
#' `exclude_cens_t = FALSE` (default) will allow individuals censored at time `t` to remain
#' in the cohort (required for pseudo-value approach).
#' Choosing `exclude_cens_t = TRUE` will also remove individuals who are censored by time `t`.
#'
#' Note that the `exclude_cens_t` argument will not work for `data_return = "data_ms"`.
#' @noRd
apply_landmark <- function(data_raw, data_ms, j, s, t, exclude_cens_t = FALSE, data_return = "data_raw"){

  ### Extract transition matrix from msdata object
  tmat <- attributes(data_ms)$trans

  ## Identify individuals in state j at time s
  ids_state_js <- extract_ids_states(data_ms = data_ms, tmat = tmat, j = j, t = s)

  ## Create landmarked dataset
  if (data_return == "data_raw"){
    if (exclude_cens_t == FALSE){
      data_lmk_js <-  data_raw |> base::subset(id %in% ids_state_js)
    } else if (exclude_cens_t == TRUE){
      data_lmk_js <-  data_raw |> base::subset(id %in% ids_state_js) |> base::subset(dtcens > t | dtcens <= t & dtcens_s == 0)
    }
  } else if (data_return == "data_ms"){
    data_lmk_js <-  data_ms |> base::subset(id %in% ids_state_js)
  }

  return(data_lmk_js)

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
identify_valid_transitions <- function(data_raw, data_ms, j, s, t){

  ### Landmark the dataset, retaining only individuals uncensored at time t
  data_raw_lmk_js <- apply_landmark(data_raw = data_raw, data_ms, j = j, s = s, t = t, exclude_cens_t = TRUE)

  ### Assign the maximum state an individual may enter
  max_state <- max(data_ms$to)

  ### Extract transition matrix from msdata object
  tmat <- attributes(data_ms)$trans

  ### Identify states these individuals are in at time t
  ids_state_list <- vector("list", max_state)
  for (k in 1:max_state){
    ids_state_list[[k]] <- extract_ids_states(data_ms, tmat, k, t)
  }

  ### Create a variable to say which state an individual was in at the time of interest
  ## Create list containing the relevant data
  v1 <- data_raw_lmk_js$id
  m1 <- outer(v1, ids_state_list, FUN = Vectorize('%in%'))
  state_poly <- lapply(split(m1, row(m1)), function(x) (1:max_state)[x])

  ## Change integer(0) values to NA's
  idx <- !sapply(state_poly, length)
  state_poly[idx] <- NA

  ## Add to data_raw
  valid_transitions <- as.numeric(names(table(unlist(state_poly))))

  return(valid_transitions)

}


#' Apply bootstrapping to a dataset of class `msdata`
#'
#' @description
#' Apply bootstrapping to datasets of class `msdata` (i.e. `data_ms`). This is non-trivial
#' because there is more than one row per individual, and a new `id` variable `id2` must be
#' assigned.
#'
#' @noRd
apply_bootstrap_msdata <- function(data_ms, indices){

  ### Break up data_ms by id
  data_ms_list <- split(data_ms, data_ms$id)

  ### Extract the relevant list elements based on id
  data_ms_boot <- lapply(1:length(indices),
                             function(x) {
                               data.frame(data_ms_list[[indices[x]]], "id2" = x)
                             }
  )
  names(data_ms_boot) <- names(data_ms_list)

  ### Combine into a single dataset and give appropriate class
  data_ms_boot <- do.call("rbind", data_ms_boot)
  rownames(data_ms_boot) <- NULL
  class(data_ms_boot) <- c("msdata", "data.frame")

  ### Return
  return(data_ms_boot)

}
