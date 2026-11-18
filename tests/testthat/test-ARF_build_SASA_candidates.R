make_sasa_tsv <- function(path, rRNA_name = "23S_ribosomal_RNA",
                          n_per_chain = 10, chains = c("A", "B")) {
  # Chain A: constant low; chain B: gradually increasing with one large outlier
  # at the last position.  Per-residue diffs become 1,2,...,n-1,2n, so the
  # outlier's z-score exceeds 1.64 and passes the selection threshold.
  rows <- do.call(rbind, lapply(chains, function(ch) {
    data.frame(
      structure_id       = "TEST",
      model              = "1",
      entity_id          = "1",
      entity_type        = "polymer",
      entity_description = rRNA_name,
      auth_chain         = ch,
      auth_comp_id       = "A",
      auth_seq_id        = as.character(seq_len(n_per_chain)),
      sasa_A2            = if (ch == chains[1])
                             rep(1.0, n_per_chain)
                           else
                             c(seq(2, n_per_chain, length.out = n_per_chain - 1L),
                               as.numeric(n_per_chain) * 2 + 1),
      stringsAsFactors   = FALSE
    )
  }))
  utils::write.table(rows, path, sep = "\t", row.names = FALSE, quote = FALSE)
  invisible(path)
}

test_that("ARF_build_SASA_candidates returns a data.frame with expected columns", {
  tsv <- tempfile(fileext = ".tsv")
  make_sasa_tsv(tsv)

  result <- ARF_build_SASA_candidates(tsv, rRNAs = "23S_ribosomal_RNA",
                                      target_species = "ec")

  expect_s3_class(result, "data.frame")
  expect_true(all(c("rRNA", "auth_seq_id", "sasa_A2_diff") %in% colnames(result)))
  expect_gt(nrow(result), 0L)
})

test_that("ARF_build_SASA_candidates prefixes rRNA column with target_species", {
  tsv <- tempfile(fileext = ".tsv")
  make_sasa_tsv(tsv)

  result <- ARF_build_SASA_candidates(tsv, rRNAs = "23S_ribosomal_RNA",
                                      target_species = "ec")

  expect_true(all(startsWith(result$rRNA, "ec_")))
})

test_that("ARF_build_SASA_candidates keeps only residues with sasa_A2_diff > 0", {
  tsv <- tempfile(fileext = ".tsv")
  # equal SASA in both chains → all diffs = 0 → nothing passes
  rows <- data.frame(
    structure_id = "TEST", model = "1", entity_id = "1",
    entity_type = "polymer", entity_description = "23S_ribosomal_RNA",
    auth_chain = rep(c("A", "B"), each = 5),
    auth_comp_id = "A",
    auth_seq_id = as.character(rep(1:5, 2)),
    sasa_A2 = rep(5.0, 10),
    stringsAsFactors = FALSE
  )
  utils::write.table(rows, tsv, sep = "\t", row.names = FALSE, quote = FALSE)

  result <- ARF_build_SASA_candidates(tsv, rRNAs = "23S_ribosomal_RNA",
                                      target_species = "ec")

  expect_equal(nrow(result), 0L)
})

test_that("ARF_build_SASA_candidates ignores entity_descriptions not in rRNAs", {
  tsv <- tempfile(fileext = ".tsv")
  rows <- data.frame(
    structure_id = "TEST", model = "1", entity_id = "1",
    entity_type = "polymer", entity_description = "ribosomal_protein_L1",
    auth_chain = rep(c("A", "B"), each = 5),
    auth_comp_id = "A",
    auth_seq_id = as.character(rep(1:5, 2)),
    sasa_A2 = c(seq(1, 5), seq(5, 1)),
    stringsAsFactors = FALSE
  )
  utils::write.table(rows, tsv, sep = "\t", row.names = FALSE, quote = FALSE)

  result <- ARF_build_SASA_candidates(tsv, rRNAs = "23S_ribosomal_RNA",
                                      target_species = "ec")

  expect_equal(nrow(result), 0L)
})

test_that("ARF_build_SASA_candidates errors on missing file", {
  expect_error(
    ARF_build_SASA_candidates("no_such_file.tsv", rRNAs = "23S_rRNA", target_species = "ec"),
    "SASA output file not found"
  )
})
