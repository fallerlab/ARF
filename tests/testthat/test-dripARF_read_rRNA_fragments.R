test_that("dripARF_read_rRNA_fragments reads mouse bedgraphs into a count matrix", {
  samples_file <- test_path("../../test_data/Rpl10aRps25Rpl22_samples_mouse.tsv")
  rRNAs_file   <- test_path("../../rRNAs/mouse_rRNAs.fa")

  skip_if_not(file.exists(samples_file), "mouse samples TSV not found")
  skip_if_not(file.exists(rRNAs_file),   "mouse rRNAs FASTA not found")

  samples <- read_ARF_samples_file(samples_file)
  # Use first two samples only to keep the test fast
  mini_samples <- samples[1:2, ]
  mini_samples$bedGraphFile <- test_path("../../", mini_samples$bedGraphFile)

  counts <- dripARF_read_rRNA_fragments(mini_samples, rRNAs_file)

  expect_s3_class(counts, "data.frame")
  expect_equal(ncol(counts), 2)
  expect_equal(colnames(counts), mini_samples$sampleName)
  # Row names should be rRNA_accession_position (e.g. NR_003278.3_1)
  expect_true(all(grepl("^NR_[0-9]+\\.[0-9]+_[0-9]+$", rownames(counts))))
  # All counts must be non-negative
  expect_true(all(counts >= 0, na.rm = TRUE))
})

test_that("dripARF_read_rRNA_fragments returns NULL for unsupported file type", {
  samples_file <- test_path("../../test_data/Rpl10aRps25Rpl22_samples_mouse.tsv")
  rRNAs_file   <- test_path("../../rRNAs/mouse_rRNAs.fa")

  skip_if_not(file.exists(samples_file), "mouse samples TSV not found")
  skip_if_not(file.exists(rRNAs_file),   "mouse rRNAs FASTA not found")

  samples <- read_ARF_samples_file(samples_file)
  bad_samples        <- samples[1, , drop = FALSE]
  bad_samples$bedGraphFile <- "notafile.bam.txt"

  result <- suppressMessages(dripARF_read_rRNA_fragments(bad_samples, rRNAs_file))
  expect_null(result)
})
