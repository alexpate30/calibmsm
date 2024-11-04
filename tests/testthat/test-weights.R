test_that("check weights, stabilised == FALSE", {

  ### Create weights for j = 1 and s = 0
  weights_manual <-
    calc_weights(data_ms = msebmtcal,
                 data_raw = ebmtcal,
                 covs = c("year", "agecl", "proph", "match"),
                 t = 1826,
                 s = 0,
                 landmark_type = "state",
                 j = 1,
                 max_weight = 10,
                 stabilised = FALSE)

  ### This should be of legnth 2279
  expect_equal(length(weights_manual$ipcw), 2279)

  ### There should be 501 NA's
  expect_equal(sum(is.na(weights_manual$ipcw)), 501)

  ### Create weights for j = 3 and s = 100
  weights_manual <-
    calc_weights(data_ms = msebmtcal,
                 data_raw = ebmtcal,
                 covs = c("year", "agecl", "proph", "match"),
                 t = 1826,
                 s = 100,
                 landmark_type = "state",
                 j = 3,
                 max_weight = 10,
                 stabilised = FALSE)

  ### This should be of legnth 2279
  expect_equal(length(weights_manual$ipcw), 2279)

  ### There should be 501 NA's
  expect_equal(sum(is.na(weights_manual$ipcw)), 802)

})

test_that("check weights, stabilised == TRUE", {

  ### Create weights for j = 1 and s = 0
  weights_manual <-
    calc_weights(data_ms = msebmtcal,
                 data_raw = ebmtcal,
                 covs = c("year", "agecl", "proph", "match"),
                 t = 1826,
                 s = 0,
                 landmark_type = "state",
                 j = 1,
                 max_weight = 10,
                 stabilised = TRUE)

  ### This should be of legnth 2279
  expect_equal(length(weights_manual$ipcw), 2279)

  ### There should be 501 NA's
  expect_equal(sum(is.na(weights_manual$ipcw)), 501)


  ### Create weights for j = 3 and s = 100
  weights_manual <-
    calc_weights(data_ms = msebmtcal,
                 data_raw = ebmtcal,
                 covs = c("year", "agecl", "proph", "match"),
                 t = 1826,
                 s = 100,
                 landmark_type = "state",
                 j = 3,
                 max_weight = 10,
                 stabilised = FALSE)

  ### This should be of legnth 2279
  expect_equal(length(weights_manual$ipcw), 2279)

  ### There should be 501 NA's
  expect_equal(sum(is.na(weights_manual$ipcw)), 802)

})

test_that("check weights, stabilised == FALSE, add max_follow", {

  ### Create weights for j = 1 and s = 0
  weights_manual <-
    calc_weights(data_ms = msebmtcal,
                 data_raw = ebmtcal,
                 covs = c("year", "agecl", "proph", "match"),
                 t = 1826,
                 s = 0,
                 landmark_type = "state",
                 j = 1,
                 max_weight = 10,
                 stabilised = FALSE,
                 max_follow = 1827)

  ### This should be of legnth 2279
  expect_equal(length(weights_manual$ipcw), 2279)

  ### There should be 501 NA's
  expect_equal(sum(is.na(weights_manual$ipcw)), 501)


  ### Create weights for j = 3 and s = 100
  weights_manual <-
    calc_weights(data_ms = msebmtcal,
                 data_raw = ebmtcal,
                 covs = c("year", "agecl", "proph", "match"),
                 t = 1826,
                 s = 100,
                 landmark_type = "state",
                 j = 3,
                 max_weight = 10,
                 stabilised = FALSE,
                 max_follow = 1827)

  ### This should be of legnth 2279
  expect_equal(length(weights_manual$ipcw), 2279)

  ### There should be 501 NA's
  expect_equal(sum(is.na(weights_manual$ipcw)), 802)

  ### Expect error if max_follow < t
  expect_error(
      calc_weights(data_ms = msebmtcal,
                   data_raw = ebmtcal,
                   covs = c("year", "agecl", "proph", "match"),
                   t = 1826,
                   s = 0,
                   landmark_type = "state",
                   j = 1,
                   max_weight = 10,
                   stabilised = FALSE,
                   max_follow = 1825)
      )

})

test_that("check weights, stabilised == FALSE, covs = null", {

  ### Create weights for j = 1 and s = 0
  weights_manual <-
    calc_weights(data_ms = msebmtcal,
                 data_raw = ebmtcal,
                 t = 1826,
                 s = 0,
                 landmark_type = "state",
                 j = 1,
                 max_weight = 10,
                 stabilised = FALSE)

  ### This should be of legnth 2279
  expect_equal(length(weights_manual$ipcw), 2279)

  ### There should be 501 NA's
  expect_equal(sum(is.na(weights_manual$ipcw)), 501)


  ### Create weights for j = 3 and s = 100
  weights_manual <-
    calc_weights(data_ms = msebmtcal,
                 data_raw = ebmtcal,
                 t = 1826,
                 s = 100,
                 landmark_type = "state",
                 j = 3,
                 max_weight = 10,
                 stabilised = FALSE)

  ### This should be of legnth 2279
  expect_equal(length(weights_manual$ipcw), 2279)

  ### There should be 501 NA's
  expect_equal(sum(is.na(weights_manual$ipcw)), 802)

})

test_that("check weights at j = 1 and s= 100, landmark_type = state/all", {

  ### Create weights for j = 1 and s = 100
  weights_manual_j1 <-
    calc_weights(data_ms = msebmtcal,
                 data_raw = ebmtcal,
                 covs = c("year", "agecl", "proph", "match"),
                 t = 1826,
                 s = 100,
                 landmark_type = "state",
                 j = 1,
                 max_weight = 10,
                 stabilised = FALSE)

  ### This should be of legnth 2279
  expect_equal(length(weights_manual_j1$ipcw), 2279)

  ### There should be 501 NA's
  expect_equal(sum(is.na(weights_manual_j1$ipcw)), 802)

  ### Create weights for j = 3 and s = 0
  weights_manual_j3 <-
    calc_weights(data_ms = msebmtcal,
                 data_raw = ebmtcal,
                 covs = c("year", "agecl", "proph", "match"),
                 t = 1826,
                 s = 100,
                 landmark_type = "state",
                 j = 3,
                 max_weight = 10,
                 stabilised = FALSE)

  ### This should be of legnth 2279
  expect_equal(length(weights_manual_j3$ipcw), 2279)

  ### There should be 501 NA's
  expect_equal(sum(is.na(weights_manual_j3$ipcw)), 802)

  ### The weights should not be the same for j = 1 and j = 3, because landmark_type == "state"
  expect_false(sum(weights_manual_j1$ipcw, na.rm = TRUE) - sum(weights_manual_j3$ipcw, na.rm = TRUE) == 0)

  ###
  ### Now repeat with landmark_type = "all", and weights should be the same
  ###

  ### Create weights for j = 1 and s = 100
  suppressWarnings(weights_manual_j1 <-
                       calc_weights(data_ms = msebmtcal,
                                    data_raw = ebmtcal,
                                    covs = c("year", "agecl", "proph", "match"),
                                    t = 1826,
                                    s = 100,
                                    landmark_type = "all",
                                    j = 1,
                                    max_weight = 10,
                                    stabilised = FALSE))

  ### This should be of legnth 2279
  expect_equal(length(weights_manual_j1$ipcw), 2279)

  ### There should be 501 NA's
  expect_equal(sum(is.na(weights_manual_j1$ipcw)), 802)

  ### Create weights for j = 3 and s = 0
  suppressWarnings(weights_manual_j3 <-
                     calc_weights(data_ms = msebmtcal,
                                  data_raw = ebmtcal,
                                  covs = c("year", "agecl", "proph", "match"),
                                  t = 1826,
                                  s = 100,
                                  landmark_type = "all",
                                  j = 3,
                                  max_weight = 10,
                                  stabilised = FALSE))

  ### This should be of legnth 2279
  expect_equal(length(weights_manual_j3$ipcw), 2279)

  ### There should be 501 NA's
  expect_equal(sum(is.na(weights_manual_j3$ipcw)), 802)

  ### The weights should not be the same for j = 1 and j = 3, because landmark_type == "state"
  expect_true(sum(weights_manual_j1$ipcw, na.rm = TRUE) - sum(weights_manual_j3$ipcw, na.rm = TRUE) == 0)

})

