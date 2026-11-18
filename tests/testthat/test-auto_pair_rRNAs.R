write_rna_fasta <- function(seqs, nms, path) {
  xs <- Biostrings::RNAStringSet(seqs)
  names(xs) <- nms
  Biostrings::writeXStringSet(xs, path)
}

test_that("auto_pair_rRNAs returns a list of two-element character vectors", {
  src <- tempfile(fileext = ".fa")
  tgt <- tempfile(fileext = ".fa")
  seq_a <- paste(rep("ACGU", 8), collapse = "")   # 32 nt
  seq_b <- paste(rep("GCUA", 16), collapse = "")  # 64 nt
  write_rna_fasta(c(seq_a, seq_b), c("src_5S", "src_16S"), src)
  write_rna_fasta(c(seq_a, seq_b), c("tgt_5S", "tgt_16S"), tgt)

  pairs <- suppressMessages(auto_pair_rRNAs(src, tgt))

  expect_type(pairs, "list")
  expect_length(pairs, 2)
  expect_true(all(lengths(pairs) == 2))
  expect_true(all(vapply(pairs, is.character, logical(1))))
})

test_that("auto_pair_rRNAs matches each source to the best-scoring target", {
  src <- tempfile(fileext = ".fa")
  tgt <- tempfile(fileext = ".fa")
  seq_a <- paste(rep("ACGU", 8), collapse = "")   # 32 nt
  seq_b <- paste(rep("GCUA", 16), collapse = "")  # 64 nt — clearly distinct by length
  write_rna_fasta(c(seq_a, seq_b), c("src_a", "src_b"), src)
  write_rna_fasta(c(seq_a, seq_b), c("tgt_a", "tgt_b"), tgt)

  pairs <- suppressMessages(auto_pair_rRNAs(src, tgt))

  expect_equal(pairs[[1]], c("src_a", "tgt_a"))
  expect_equal(pairs[[2]], c("src_b", "tgt_b"))
})

test_that("auto_pair_rRNAs prints pairing messages", {
  src <- tempfile(fileext = ".fa")
  tgt <- tempfile(fileext = ".fa")
  write_rna_fasta("ACGUACGUACGU", "src_5S", src)
  write_rna_fasta("ACGUACGUACGU", "tgt_5S", tgt)

  expect_message(auto_pair_rRNAs(src, tgt), "Auto-paired")
})

test_that("auto_pair_rRNAs strips FASTA header spaces (uses first word only)", {
  src <- tempfile(fileext = ".fa")
  tgt <- tempfile(fileext = ".fa")
  writeLines(c(">src_5S some description text", "ACGUACGUACGU"), src)
  writeLines(c(">tgt_5S rRNA molecule", "ACGUACGUACGU"), tgt)

  pairs <- suppressMessages(auto_pair_rRNAs(src, tgt))

  expect_equal(pairs[[1]][1], "src_5S")
  expect_equal(pairs[[1]][2], "tgt_5S")
})

test_that("auto_pair_rRNAs handles DNA-alphabet input (T→U conversion)", {
  src <- tempfile(fileext = ".fa")
  tgt <- tempfile(fileext = ".fa")
  # write with T (DNA) — function should convert to U before alignment
  writeLines(c(">src_5S", "ACGTACGTACGT"), src)
  writeLines(c(">tgt_5S", "ACGUACGUACGU"), tgt)

  expect_no_error(suppressMessages(auto_pair_rRNAs(src, tgt)))
})
