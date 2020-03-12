context("Summarising MiDAS results")

test_that("summarize amino acid position", {
  file <- system.file("extdata", "HLAHD_output_example.txt", package = "MiDAS")
  hla_calls <- readHlaCalls(file)
  aa_sum <- summariseAAPosition(hla_calls, "DMA_135")

  aa_sum_test <- data.frame(
    `HLA-DMA (135)` = c("K"),
    `HLA-DMA alleles` = "*01:01, *01:02, *01:03",
    count = 40,
    frequency = formattable::percent(1),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  expect_equal(aa_sum, aa_sum_test)

  # checkHlaCalls Format test is ommited here
  expect_error(
    summariseAAPosition(hla_calls, NA),
    "pos is not a string \\(a length one character vector\\)."
  )

  expect_error(
    summariseAAPosition(hla_calls, "1"),
    "amino acid position should be formatted like: A_9."
  )

  expect_error(
    summariseAAPosition(hla_calls, "A_1", aln = "foo"),
    "aln is not a matrix or NULL."
  )

  expect_error(
    summariseAAPosition(hla_calls, "A_1", na.rm = "foo"),
    "na.rm is not a flag \\(a length one logical vector\\)."
  )

  hla_calls[, c("DMA_1", "DMA_2")] <- NA
  expect_error(
    summariseAAPosition(hla_calls, "DMA_100"),
    "hla_calls for given gene contains only NA."
  )

  expect_error(
    summariseAAPosition(hla_calls, "DRB1_10000"),
    "amino acid position 10000 is higher than amino acid sequence length."
  )

  hla_calls[1, "DRB1_1"] <- "DRB1*00:00"
  expect_error(
    summariseAAPosition(hla_calls, "DRB1_100"),
    "allele DRB1\\*00:00 could not be found in the nucleotide alignment file."
  )
})
