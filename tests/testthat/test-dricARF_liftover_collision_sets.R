ecoli_fa <- system.file("extdata", "EscherichiaColi.fa", package = "ARF")

# dricARF_liftover_collision_sets uses RP_proximity_df from the calling scope
# (a design constraint of the function). Supply a minimal stub so the function
# can run standalone in tests.
make_stub_rp_df <- function(rRNA_name, n = 120) {
  data.frame(
    rRNA  = rep(rRNA_name, n),
    resno = seq_len(n),
    uL18  = runif(n, 1, 100),
    stringsAsFactors = FALSE
  )
}

test_that("dricARF_liftover_collision_sets returns ont/gene data.frame", {
  skip_if_not(nzchar(ecoli_fa) && file.exists(ecoli_fa), "EscherichiaColi.fa not found")
  skip_on_cran()

  RP_proximity_df <<- make_stub_rp_df("ec_5S_rRNA")
  on.exit(rm("RP_proximity_df", envir = .GlobalEnv), add = TRUE)

  result <- suppressMessages(
    dricARF_liftover_collision_sets(
      target_species     = "ec",
      target_rRNAs_fasta = ecoli_fa,
      rRNA_pairs         = list(c("5S", "ec_5S_rRNA"))
    )
  )

  expect_s3_class(result, "data.frame")
  expect_named(result, c("ont", "gene"))
})

test_that("dricARF_liftover_collision_sets gene labels carry target_species prefix", {
  skip_if_not(nzchar(ecoli_fa) && file.exists(ecoli_fa), "EscherichiaColi.fa not found")
  skip_on_cran()

  RP_proximity_df <<- make_stub_rp_df("ec_5S_rRNA")
  on.exit(rm("RP_proximity_df", envir = .GlobalEnv), add = TRUE)

  result <- suppressMessages(
    dricARF_liftover_collision_sets(
      target_species     = "ec",
      target_rRNAs_fasta = ecoli_fa,
      rRNA_pairs         = list(c("5S", "ec_5S_rRNA"))
    )
  )

  expect_true(any(grepl("^ec_5S_rRNA_", result$gene)))
})

test_that("dricARF_liftover_collision_sets includes Rand sets alongside real sets", {
  skip_if_not(nzchar(ecoli_fa) && file.exists(ecoli_fa), "EscherichiaColi.fa not found")
  skip_on_cran()

  RP_proximity_df <<- make_stub_rp_df("ec_5S_rRNA")
  on.exit(rm("RP_proximity_df", envir = .GlobalEnv), add = TRUE)

  result <- suppressMessages(
    dricARF_liftover_collision_sets(
      target_species     = "ec",
      target_rRNAs_fasta = ecoli_fa,
      rRNA_pairs         = list(c("5S", "ec_5S_rRNA"))
    )
  )

  ont_vals <- unique(result$ont)
  expect_true(any(grepl("^Rand[0-9]+_", ont_vals)))
})
