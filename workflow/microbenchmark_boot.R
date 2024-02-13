### Create a new id for these individuals (calc_pv_aj relies on each individual having a unique identifier),
### meaning the duplicate values in the bootstrapped datasets will cause problems
load_all()
data.raw <- ebmtcal
data.mstate <- msebmtcal

indices <- 1:nrow(data.raw)
indices <- sample(1:nrow(data.raw), nrow(data.raw))
### Create bootstrapped dataset
data.raw.boot <- data.raw[indices, ]

data.raw.boot$id2 <- 1:nrow(data.raw.boot)

### Create bootstrapped data.mstate (we replicate the choice of patients that was chosen in data.raw)
myfunc1 <- function(data.mstate, data.raw.boot){
  out <- do.call("rbind",
                 lapply(1:nrow(data.raw.boot),
                        function(x) {
                          base::subset(data.mstate, id == data.raw.boot$id[x]) |>
                            dplyr::mutate(id2 = data.raw.boot$id2[x])
                        }
                 )
  )

  rownames(out) <- NULL
  return(out)
}


myfunc2 <- function(data.mstate, indices){

  ### Break up data.mstate by id
  data.mstate.list <- split(data.mstate, data.mstate$id)

  ### Extract the relevant list elements based on id
  data.mstate.boot <- lapply(1:length(indices),
                                  function(x) {
                                    data.frame(data.mstate.list[[indices[x]]], "id2" = x)
                                  }
  )
  names(data.mstate.boot) <- names(data.mstate.list)

  ### Combine into a single dataset and give appropriate class
  data.mstate.boot <- do.call("rbind", data.mstate.boot)
  rownames(data.mstate.boot) <- NULL
  class(data.mstate.boot) <- c("msdata", "data.frame")

  ### Return
  return(data.mstate.boot)

}

test1 <- myfunc1(data.mstate, data.raw.boot)
test2 <- myfunc2(data.mstate, indices)
testthat::expect_equal(test1, test2)

str(data.mstate.list)
str(data.mstate.list.boot)
names(data.mstate.list.boot)
names(data.mstate.list)
testthat::expect_equal(data.mstate.list, data.mstate.list.boot)

str(data.mstate.list[[1]])


microbenchmark::microbenchmark(myfunc1(data.mstate, data.raw.boot), myfunc2(data.mstate, data.raw.boot, indices), times = 10)

