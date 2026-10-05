test_that("read_ARF_samples_file reads mouse samples TSV correctly", {
  samples_file <- test_path("../../test_data/Rpl10aRps25Rpl22_samples_mouse.tsv")
  skip_if_not(file.exists(samples_file), "mouse samples TSV not found")

  df <- read_ARF_samples_file(samples_file)

  expect_s3_class(df, "data.frame")
  expect_equal(nrow(df), 12)
  expect_true(all(c("sampleName", "bedGraphFile", "group", "taggedRP", "ribosomeIsolation") %in% colnames(df)))
  expect_equal(rownames(df), df$sampleName)
})

test_that("read_ARF_samples_file rownames equal first column", {
  samples_file <- test_path("../../test_data/Rpl10aRps25Rpl22_samples_mouse.tsv")
  skip_if_not(file.exists(samples_file), "mouse samples TSV not found")

  df <- read_ARF_samples_file(samples_file)
  expect_equal(rownames(df), df[, 1])
})

test_that("read_ARF_samples_file reads E. coli samples TSV with extra columns", {
  samples_file <- test_path("../../test_data/Escherichia_PRJNA906596/riboseq/samples.tsv")
  skip_if_not(file.exists(samples_file), "E. coli samples TSV not found")

  df <- read_ARF_samples_file(samples_file)

  expect_s3_class(df, "data.frame")
  expect_equal(nrow(df), 9)
  expect_true("taggedRP" %in% colnames(df))
  expect_true("ribosomeIsolation" %in% colnames(df))
})
