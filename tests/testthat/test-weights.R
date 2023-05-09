test_that("check weights output has correct length and number of NA's, stabilised == FALSE", {

  ### Create weights for j = 1 and s = 0
  weights.manual <-
    calc_weights(data.mstate = msebmtcal,
                 data.raw = ebmtcal,
                 covs = c("year", "agecl", "proph", "match"),
                 t.eval = 1826,
                 s = 0,
                 landmark.type = "state",
                 j = 1,
                 max.weight = 10,
                 stabilised = FALSE)

  ### This should be of legnth 2279
  expect_equal(length(weights.manual$ipcw), 2279)

  ### There should be 501 NA's
  expect_equal(sum(is.na(weights.manual$ipcw)), 501)
  expect_equal(sum(is.na(weights.manual$pcw)), 501)

  ### Create weights for j = 3 and s = 100
  weights.manual <-
    calc_weights(data.mstate = msebmtcal,
                 data.raw = ebmtcal,
                 covs = c("year", "agecl", "proph", "match"),
                 t.eval = 1826,
                 s = 100,
                 landmark.type = "state",
                 j = 3,
                 max.weight = 10,
                 stabilised = FALSE)

  ### This should be of legnth 2279
  expect_equal(length(weights.manual$ipcw), 2279)

  ### There should be 501 NA's
  expect_equal(sum(is.na(weights.manual$ipcw)), 802)
  expect_equal(sum(is.na(weights.manual$pcw)), 802)

})

test_that("check weights output has correct length and number of NA's, stabilised == TRUE", {

  ### Create weights for j = 1 and s = 0
  weights.manual <-
    calc_weights(data.mstate = msebmtcal,
                 data.raw = ebmtcal,
                 covs = c("year", "agecl", "proph", "match"),
                 t.eval = 1826,
                 s = 0,
                 landmark.type = "state",
                 j = 1,
                 max.weight = 10,
                 stabilised = TRUE)

  ### This should be of legnth 2279
  expect_equal(length(weights.manual$ipcw), 2279)

  ### There should be 501 NA's
  expect_equal(sum(is.na(weights.manual$ipcw)), 501)
  expect_equal(sum(is.na(weights.manual$pcw)), 501)

  ### Create weights for j = 3 and s = 100
  weights.manual <-
    calc_weights(data.mstate = msebmtcal,
                 data.raw = ebmtcal,
                 covs = c("year", "agecl", "proph", "match"),
                 t.eval = 1826,
                 s = 100,
                 landmark.type = "state",
                 j = 3,
                 max.weight = 10,
                 stabilised = FALSE)

  ### This should be of legnth 2279
  expect_equal(length(weights.manual$ipcw), 2279)

  ### There should be 501 NA's
  expect_equal(sum(is.na(weights.manual$ipcw)), 802)
  expect_equal(sum(is.na(weights.manual$pcw)), 802)

})

