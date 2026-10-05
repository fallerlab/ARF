test_that("ARF_check_organism returns TRUE for valid organism codes", {
  expect_true(ARF_check_organism("hs"))
  expect_true(ARF_check_organism("mm"))
  expect_true(ARF_check_organism("sc"))
})

test_that("ARF_check_organism returns FALSE for unknown organism codes", {
  expect_false(suppressMessages(ARF_check_organism("at")))
  expect_false(suppressMessages(ARF_check_organism("ec")))
  expect_false(suppressMessages(ARF_check_organism("human")))
  expect_false(suppressMessages(ARF_check_organism("")))
})

test_that("ARF_check_organism emits a message for invalid codes", {
  expect_message(ARF_check_organism("at"))
})
