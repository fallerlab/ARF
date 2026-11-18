rcsb_fasta <- test_path("../../test_data/rcsb_pdb_8B2L.fasta")

test_that("ARF_parse_structure_fasta returns a named list with conversion_table and rRNA_seqs", {
  skip_if_not(file.exists(rcsb_fasta), "8B2L FASTA not found")

  result <- ARF_parse_structure_fasta(rcsb_fasta, target_species = "AT")

  expect_type(result, "list")
  expect_named(result, c("conversion_table", "rRNA_seqs"))
})

test_that("ARF_parse_structure_fasta conversion_table has required columns", {
  skip_if_not(file.exists(rcsb_fasta), "8B2L FASTA not found")

  ct <- ARF_parse_structure_fasta(rcsb_fasta, target_species = "AT")$conversion_table

  expect_s3_class(ct, "data.frame")
  expect_true(all(c("names", "chainID", "RP_name", "organism", "RP_new") %in% colnames(ct)))
  expect_gt(nrow(ct), 0)
})

test_that("ARF_parse_structure_fasta extracts only rRNA sequences", {
  skip_if_not(file.exists(rcsb_fasta), "8B2L FASTA not found")

  rRNA_seqs <- ARF_parse_structure_fasta(rcsb_fasta, target_species = "AT")$rRNA_seqs

  expect_s4_class(rRNA_seqs, "BStringSet")
  expect_gt(length(rRNA_seqs), 0)
  # All names should contain "rRNA" in the original FASTA description
  raw_names <- names(Biostrings::readBStringSet(rcsb_fasta))
  rRNA_raw <- raw_names[grepl("[0-9]S", raw_names) & grepl("RNA", raw_names)]
  expect_equal(length(rRNA_seqs), length(rRNA_raw))
})

test_that("ARF_parse_structure_fasta prefixes rRNA names with target_species", {
  skip_if_not(file.exists(rcsb_fasta), "8B2L FASTA not found")

  rRNA_seqs <- ARF_parse_structure_fasta(rcsb_fasta, target_species = "AT")$rRNA_seqs

  expect_true(all(startsWith(names(rRNA_seqs), "AT_")))
})

test_that("ARF_parse_structure_fasta writes rRNA FASTA when rRNAs_output_file is given", {
  skip_if_not(file.exists(rcsb_fasta), "8B2L FASTA not found")

  out <- tempfile(fileext = ".fa")
  result <- ARF_parse_structure_fasta(rcsb_fasta, target_species = "AT",
                                      rRNAs_output_file = out)

  expect_true(file.exists(out))
  written <- Biostrings::readBStringSet(out)
  expect_equal(length(written), length(result$rRNA_seqs))
})

test_that("ARF_parse_structure_fasta errors on missing file", {
  expect_error(ARF_parse_structure_fasta("does_not_exist.fasta", "AT"),
               "does not exist")
})
