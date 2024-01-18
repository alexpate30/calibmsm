test_that("check weights, stabilised == FALSE", {

  ### Create weights for j = 1 and s = 0
  weights.manual <-
    calc_weights(data.mstate = msebmtcal,
                 data.raw = ebmtcal,
                 covs = c("year", "agecl", "proph", "match"),
                 t = 1826,
                 s = 0,
                 landmark.type = "state",
                 j = 1,
                 max.weight = 10,
                 stabilised = FALSE)

  ### This should be of legnth 2279
  expect_equal(length(weights.manual$ipcw), 2279)

  ### There should be 501 NA's
  expect_equal(sum(is.na(weights.manual$ipcw)), 501)

  ### Create weights for j = 3 and s = 100
  weights.manual <-
    calc_weights(data.mstate = msebmtcal,
                 data.raw = ebmtcal,
                 covs = c("year", "agecl", "proph", "match"),
                 t = 1826,
                 s = 100,
                 landmark.type = "state",
                 j = 3,
                 max.weight = 10,
                 stabilised = FALSE)

  ### This should be of legnth 2279
  expect_equal(length(weights.manual$ipcw), 2279)

  ### There should be 501 NA's
  expect_equal(sum(is.na(weights.manual$ipcw)), 802)

})

test_that("check weights, stabilised == TRUE", {

  ### Create weights for j = 1 and s = 0
  weights.manual <-
    calc_weights(data.mstate = msebmtcal,
                 data.raw = ebmtcal,
                 covs = c("year", "agecl", "proph", "match"),
                 t = 1826,
                 s = 0,
                 landmark.type = "state",
                 j = 1,
                 max.weight = 10,
                 stabilised = TRUE)

  ### This should be of legnth 2279
  expect_equal(length(weights.manual$ipcw), 2279)

  ### There should be 501 NA's
  expect_equal(sum(is.na(weights.manual$ipcw)), 501)


  ### Create weights for j = 3 and s = 100
  weights.manual <-
    calc_weights(data.mstate = msebmtcal,
                 data.raw = ebmtcal,
                 covs = c("year", "agecl", "proph", "match"),
                 t = 1826,
                 s = 100,
                 landmark.type = "state",
                 j = 3,
                 max.weight = 10,
                 stabilised = FALSE)

  ### This should be of legnth 2279
  expect_equal(length(weights.manual$ipcw), 2279)

  ### There should be 501 NA's
  expect_equal(sum(is.na(weights.manual$ipcw)), 802)

})

test_that("check weights, stabilised == FALSE, add max.follow", {

  ### Create weights for j = 1 and s = 0
  weights.manual <-
    calc_weights(data.mstate = msebmtcal,
                 data.raw = ebmtcal,
                 covs = c("year", "agecl", "proph", "match"),
                 t = 1826,
                 s = 0,
                 landmark.type = "state",
                 j = 1,
                 max.weight = 10,
                 stabilised = FALSE,
                 max.follow = 1827)

  ### This should be of legnth 2279
  expect_equal(length(weights.manual$ipcw), 2279)

  ### There should be 501 NA's
  expect_equal(sum(is.na(weights.manual$ipcw)), 501)


  ### Create weights for j = 3 and s = 100
  weights.manual <-
    calc_weights(data.mstate = msebmtcal,
                 data.raw = ebmtcal,
                 covs = c("year", "agecl", "proph", "match"),
                 t = 1826,
                 s = 100,
                 landmark.type = "state",
                 j = 3,
                 max.weight = 10,
                 stabilised = FALSE,
                 max.follow = 1827)

  ### This should be of legnth 2279
  expect_equal(length(weights.manual$ipcw), 2279)

  ### There should be 501 NA's
  expect_equal(sum(is.na(weights.manual$ipcw)), 802)

  ### Expect error if max.follow < t
  expect_error(
      calc_weights(data.mstate = msebmtcal,
                   data.raw = ebmtcal,
                   covs = c("year", "agecl", "proph", "match"),
                   t = 1826,
                   s = 0,
                   landmark.type = "state",
                   j = 1,
                   max.weight = 10,
                   stabilised = FALSE,
                   max.follow = 1825)
      )

})

test_that("check weights, stabilised == FALSE, covs = null", {

  ### Create weights for j = 1 and s = 0
  weights.manual <-
    calc_weights(data.mstate = msebmtcal,
                 data.raw = ebmtcal,
                 t = 1826,
                 s = 0,
                 landmark.type = "state",
                 j = 1,
                 max.weight = 10,
                 stabilised = FALSE)

  ### This should be of legnth 2279
  expect_equal(length(weights.manual$ipcw), 2279)

  ### There should be 501 NA's
  expect_equal(sum(is.na(weights.manual$ipcw)), 501)


  ### Create weights for j = 3 and s = 100
  weights.manual <-
    calc_weights(data.mstate = msebmtcal,
                 data.raw = ebmtcal,
                 t = 1826,
                 s = 100,
                 landmark.type = "state",
                 j = 3,
                 max.weight = 10,
                 stabilised = FALSE)

  ### This should be of legnth 2279
  expect_equal(length(weights.manual$ipcw), 2279)

  ### There should be 501 NA's
  expect_equal(sum(is.na(weights.manual$ipcw)), 802)

})

test_that("check weights at j = 1 and s= 100, landmark.type = state/all", {

  ### Create weights for j = 1 and s = 100
  weights.manual.j1 <-
    calc_weights(data.mstate = msebmtcal,
                 data.raw = ebmtcal,
                 covs = c("year", "agecl", "proph", "match"),
                 t = 1826,
                 s = 100,
                 landmark.type = "state",
                 j = 1,
                 max.weight = 10,
                 stabilised = FALSE)

  ### This should be of legnth 2279
  expect_equal(length(weights.manual.j1$ipcw), 2279)

  ### There should be 501 NA's
  expect_equal(sum(is.na(weights.manual.j1$ipcw)), 802)

  ### Create weights for j = 3 and s = 0
  weights.manual.j3 <-
    calc_weights(data.mstate = msebmtcal,
                 data.raw = ebmtcal,
                 covs = c("year", "agecl", "proph", "match"),
                 t = 1826,
                 s = 100,
                 landmark.type = "state",
                 j = 3,
                 max.weight = 10,
                 stabilised = FALSE)

  ### This should be of legnth 2279
  expect_equal(length(weights.manual.j3$ipcw), 2279)

  ### There should be 501 NA's
  expect_equal(sum(is.na(weights.manual.j3$ipcw)), 802)

  ### The weights should not be the same for j = 1 and j = 3, because landmark.type == "state"
  expect_false(sum(weights.manual.j1$ipcw, na.rm = TRUE) - sum(weights.manual.j3$ipcw, na.rm = TRUE) == 0)

  ###
  ### Now repeat with landmark.type = "all", and weights should be the same
  ###

  ### Create weights for j = 1 and s = 100
  suppressWarnings(weights.manual.j1 <-
                       calc_weights(data.mstate = msebmtcal,
                                    data.raw = ebmtcal,
                                    covs = c("year", "agecl", "proph", "match"),
                                    t = 1826,
                                    s = 100,
                                    landmark.type = "all",
                                    j = 1,
                                    max.weight = 10,
                                    stabilised = FALSE))

  ### This should be of legnth 2279
  expect_equal(length(weights.manual.j1$ipcw), 2279)

  ### There should be 501 NA's
  expect_equal(sum(is.na(weights.manual.j1$ipcw)), 802)

  ### Create weights for j = 3 and s = 0
  suppressWarnings(weights.manual.j3 <-
                     calc_weights(data.mstate = msebmtcal,
                                  data.raw = ebmtcal,
                                  covs = c("year", "agecl", "proph", "match"),
                                  t = 1826,
                                  s = 100,
                                  landmark.type = "all",
                                  j = 3,
                                  max.weight = 10,
                                  stabilised = FALSE))

  ### This should be of legnth 2279
  expect_equal(length(weights.manual.j3$ipcw), 2279)

  ### There should be 501 NA's
  expect_equal(sum(is.na(weights.manual.j3$ipcw)), 802)

  ### The weights should not be the same for j = 1 and j = 3, because landmark.type == "state"
  expect_true(sum(weights.manual.j1$ipcw, na.rm = TRUE) - sum(weights.manual.j3$ipcw, na.rm = TRUE) == 0)

})

