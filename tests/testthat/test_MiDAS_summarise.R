context("Summarising MiDAS results")

test_that("summarize amino acid position", {
  aa_sum <- summariseAAPosition(MiDAS_tut_HLA, "DRA_2")
  aa_sum_test <- data.frame(
    `HLA-DRA (2)` = c("*", "K"),
    `HLA-DRA alleles` = c("*01:02", "*01:01"),
    count = c(711L, 1289L),
    frequency = formattable::percent(c(0.3555, 0.6445)),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  expect_equal(aa_sum, aa_sum_test)
  
  aa_sum <- summariseAAPosition(MiDAS_tut_HLA, "DRA_-25")
  aa_sum_test <- data.frame(
    `HLA-DRA (-25)` = c("*", "M"),
    `HLA-DRA alleles` = c("*01:02", "*01:01"),
    count = c(711L, 1289L),
    frequency = formattable::percent(c(0.3555, 0.6445)),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  expect_equal(aa_sum, aa_sum_test)

  # checkHlaCalls Format test is ommited here
  expect_error(
    summariseAAPosition(MiDAS_tut_HLA, NA),
    "pos is not a string \\(a length one character vector\\)."
  )

  expect_error(
    summariseAAPosition(MiDAS_tut_HLA, "1"),
    "amino acid position should be formatted like: A_9."
  )

  expect_error(
    summariseAAPosition(MiDAS_tut_HLA, "A_1", aln = "foo"),
    "aln is not a matrix or NULL."
  )

  expect_error(
    summariseAAPosition(MiDAS_tut_HLA, "A_1", na.rm = "foo"),
    "na.rm is not a flag \\(a length one logical vector\\)."
  )

  hla_calls <- MiDAS_tut_HLA
  hla_calls[, c("DMA_1", "DMA_2")] <- NA
  expect_error(
    summariseAAPosition(hla_calls, "DMA_100"),
    "hla_calls for given gene contains only NA."
  )

  expect_error(
    summariseAAPosition(MiDAS_tut_HLA, "DRB1_10000"),
    "amino acid position 10000 was not found in amino acid sequence."
  )

  hla_calls <- MiDAS_tut_HLA
  hla_calls[1, "DRB1_1"] <- "DRB1*00:00"
  expect_error(
    summariseAAPosition(hla_calls, "DRB1_100"),
    "allele DRB1\\*00:00 could not be found in the alignment file."
  )
})
