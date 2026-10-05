make_mini_proximity_df <- function(n = 200, seed = 42) {
  set.seed(seed)
  data.frame(
    rRNA  = rep("rRNA_28S", n),
    resno = seq_len(n),
    eL22  = c(runif(40, 1, 15), runif(n - 40, 25, 80)),
    uL1   = c(runif(n - 40, 25, 80), runif(40, 1, 15)),
    stringsAsFactors = FALSE
  )
}

test_that("dripARF_get_RP_proximity_sets returns ont/gene data.frame", {
  result <- dripARF_get_RP_proximity_sets(make_mini_proximity_df())

  expect_s3_class(result, "data.frame")
  expect_named(result, c("ont", "gene"))
})

test_that("dripARF_get_RP_proximity_sets includes both RPs and their Rand sets", {
  result <- dripARF_get_RP_proximity_sets(make_mini_proximity_df())

  expect_true("eL22" %in% result$ont)
  expect_true("uL1" %in% result$ont)
  expect_true(any(grepl("^Rand[0-9]+_eL22$", result$ont)))
  expect_true(any(grepl("^Rand[0-9]+_uL1$", result$ont)))
})

test_that("dripARF_get_RP_proximity_sets gene labels are rRNA_resno format", {
  result <- dripARF_get_RP_proximity_sets(make_mini_proximity_df())

  expect_true(all(grepl("^rRNA_28S_[0-9]+$", result$gene)))
})

test_that("dripARF_get_RP_proximity_sets accepts explicit thresholds", {
  prox <- make_mini_proximity_df()
  result <- dripARF_get_RP_proximity_sets(prox, thresholds = c(15, 20))

  expect_s3_class(result, "data.frame")
  expect_named(result, c("ont", "gene"))
  # With dist threshold=15 and set-size cap=20, eL22 set should be at most 20 genes
  eL22_genes <- result$gene[result$ont == "eL22"]
  expect_lte(length(eL22_genes), 20)
})

test_that("dripARF_get_RP_proximity_sets returns NULL when no positions pass threshold", {
  prox <- make_mini_proximity_df()
  # Threshold of 0 Angstrom: nothing should pass
  result <- dripARF_get_RP_proximity_sets(prox, thresholds = c(0, 5))

  expect_null(result)
})
