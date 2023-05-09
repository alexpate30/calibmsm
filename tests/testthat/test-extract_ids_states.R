test_that("check extract_ids_states assigns correct states for each individual", {

  ### Extract ids in each state at time t.eval
  states <- lapply(1:6, extract_ids_states, data.mstate = msebmtcal, tmat = attributes(msebmtcal)$trans, t.eval = 1826)

  expect_equal(sum(is.na(states[[1]])), 0)
  expect_equal(sum(is.na(states[[2]])), 0)
  expect_equal(sum(is.na(states[[3]])), 0)
  expect_equal(sum(is.na(states[[4]])), 0)
  expect_equal(sum(is.na(states[[5]])), 0)
  expect_equal(sum(is.na(states[[6]])), 0)

  expect_equal(length(states[[1]]), 213)
  expect_equal(length(states[[2]]), 246)
  expect_equal(length(states[[3]]), 190)
  expect_equal(length(states[[4]]), 268)
  expect_equal(length(states[[5]]), 357)
  expect_equal(length(states[[6]]), 504)

})
