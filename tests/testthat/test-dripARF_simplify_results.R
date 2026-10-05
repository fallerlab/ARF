make_mock_results <- function() {
  data.frame(
    comp           = c("A_vs_B", "A_vs_B", "A_vs_B"),
    Description    = c("eL22", "uL1", "eS25"),
    C1.avg.read.c  = c(100, 200, 150),
    C2.avg.read.c  = c(150, 180, 120),
    RPSEA.NES      = c(1.5, 0.8, 1.2),
    RPSEA.NES_randZ = c(1.2, 0.5, 0.9),
    RPSEA.padj     = c(0.01, 0.20, 0.06),
    ORA.padj       = c(0.03, 0.40, 0.04),
    ORA.overlap    = c(5L, 2L, 3L),
    stringsAsFactors = FALSE
  )
}

test_that("dripARF_simplify_results returns a data.frame with required columns", {
  result <- dripARF_simplify_results(make_mock_results())

  expect_s3_class(result, "data.frame")
  expected_cols <- c("comp", "Description", "top", "ES1", "ES2",
                     "ORA", "RPSEA.padj", "ES1.pass", "ES2.pass", "ORA.sigN",
                     "C1.avg.read.c", "C2.avg.read.c")
  expect_named(result, expected_cols)
})

test_that("dripARF_simplify_results 'top' flag uses all four filters", {
  result <- dripARF_simplify_results(make_mock_results())

  # eL22: padj=0.01 (<0.05), randZ=1.2 (>=1), ORA.padj=0.03 (<0.05), overlap=5 (>=0) -> TRUE
  expect_true(result$top[result$Description == "eL22"])
  # uL1: padj=0.20 fails RPSEA threshold -> FALSE
  expect_false(result$top[result$Description == "uL1"])
  # eS25: RPSEA.padj=0.06 fails -> FALSE
  expect_false(result$top[result$Description == "eS25"])
})

test_that("dripARF_simplify_results respects custom thresholds", {
  # Relax RPSEA_adjP_thr to 0.1; eS25 (padj=0.06) should still fail randZ < 1
  result <- dripARF_simplify_results(make_mock_results(),
                                     RPSEA_adjP_thr = 0.1, randZscore_thr = 1)
  expect_false(result$top[result$Description == "eS25"])

  # Set randZscore_thr very low so eS25 passes
  result2 <- dripARF_simplify_results(make_mock_results(),
                                      RPSEA_adjP_thr = 0.1, randZscore_thr = 0)
  expect_true(result2$top[result2$Description == "eS25"])
})

test_that("dripARF_simplify_results preserves row count", {
  mock <- make_mock_results()
  result <- dripARF_simplify_results(mock)
  expect_equal(nrow(result), nrow(mock))
})
