context("Transforming MiDAS objects")

test_that("Amino acids variability is infered correctly", {
  hla_calls <- system.file("extdata/HLAHD_output_example.txt", package = "MiDAS")
  hla_calls <- readHlaCalls(hla_calls)
  aa_variation <-
    hlaToAAVariation(hla_calls,
                     indels = TRUE,
                     unkchar = TRUE,
                     as_df = FALSE)
  load(system.file("extdata", "test_aa_variation.Rdata", package = "MiDAS"))
  expect_equal(aa_variation, test_aa_variation)

  expect_error(hlaToAAVariation(hla_calls, indels = "foo"),
               "indels is not a flag \\(a length one logical vector\\)."
  )

  expect_error(hlaToAAVariation(hla_calls, unkchar = "foo"),
               "unkchar is not a flag \\(a length one logical vector\\)."
  )

  expect_error(hlaToAAVariation(hla_calls, alnpath = "foo/bar/foo/bar"),
               "Path 'foo/bar/foo/bar' does not exist"
  )

  expect_error(hlaToAAVariation(hla_calls, alnpath = system.file(package = "MiDAS")),
               sprintf("no alignment files was found in path %s", system.file("inst", package = "MiDAS"))
  )

  expect_error(hlaToAAVariation(hla_calls, as_df = 1),
               "as_df is not a flag \\(a length one logical vector\\)."
  )
})

test_that("HLA calls table is converted to additional variables", {
  hla_calls <- system.file("extdata/HLAHD_output_example.txt",
                           package = "MiDAS"
  )
  hla_calls <- readHlaCalls(hla_calls)
  hla_supertypes <- hlaToVariable(hla_calls, dictionary = "4digit_supertype")
  load(system.file("extdata", "test_hla_supertypes.RData", package = "MiDAS"))
  expect_equal(hla_supertypes, test_hla_supertypes)

  expect_error(
    hlaToVariable(c("A*01:01", "A*02:01"), dictionary = "4digit_supertype"),
    "hla_calls is not a data frame"
  )

  expect_error(
    hlaToVariable(
      hla_calls = hla_calls,
      dictionary = "4digit_supertype",
      nacols.rm = "yes"
    ),
    "nacols.rm is not a flag \\(a length one logical vector\\)."
  )
})

test_that("HLA calls table is converted to counts table", {
  hla_calls <- system.file("extdata/HLAHD_output_example.txt",
                           package = "MiDAS"
  )
  hla_calls <- readHlaCalls(hla_calls)
  load(system.file("extdata", "test_hla_counts.RData", package = "MiDAS"))

  hla_counts <- hlaCallsToCounts(hla_calls, inheritance_model = "dominant")
  expect_equal(hla_counts, test_hla_counts[["dominant"]])

  hla_counts <- hlaCallsToCounts(hla_calls, inheritance_model = "recessive")
  expect_equal(hla_counts, test_hla_counts[["recessive"]])

  hla_counts <- hlaCallsToCounts(hla_calls, inheritance_model = "additive")
  expect_equal(hla_counts, test_hla_counts[["additive"]])

  expect_error(
    hlaCallsToCounts(c("A*01:01", "A*02:01"), inheritance_model = "additive"),
    "hla_calls is not a data frame"
  )

  expect_error(
    hlaCallsToCounts(hla_calls, inheritance_model = 123),
    "inheritance_model is not a string \\(a length one character vector\\)."
  )

  expect_error(
    hlaCallsToCounts(hla_calls, inheritance_model = "foo"),
    "inheritance_model should be one of 'dominant', 'recessive', 'additive'"
  )
})

test_that("hla frequencies are calculated properly", {
  minimal_hla_calls <- data.frame(
    ID = c("P1", "P2"),
    A_1 = c("A*01:01", "A*02:01"),
    A_2 = c("A*02:01", "A*01:01"),
    stringsAsFactors = FALSE
  )
  hla_freq <- getHlaFrequencies(minimal_hla_calls)
  test_hla_freq <- data.frame(
    allele = c("A*01:01", "A*02:01"),
    Freq = c(0.5, 0.5),
    stringsAsFactors = FALSE
  )
  expect_equal(hla_freq, test_hla_freq)

  # checkHlaCallsFormat tests are ommited here
})

test_that("amino acids variation data frame is converted to counts table", {
  minimal_hla_calls <- data.frame(
    ID = c("P1", "P2"),
    A_1 = c("A*01:01", "A*02:01"),
    A_2 = c("A*02:01", "A*01:01"),
    stringsAsFactors = FALSE
  )
  aa_var <- hlaToAAVariation(minimal_hla_calls)[, 1:5]

  aa_counts <- aaVariationToCounts(aa_var, inheritance_model = "additive")
  test_aa_counts <- data.frame(
    ID = c("P1", "P2"),
    A_44_K = c(1, 1),
    A_44_R = c(1, 1),
    A_62_Q = c(1, 1),
    A_62_G = c(1, 1),
    stringsAsFactors = FALSE
  )
  expect_equal(aa_counts, test_aa_counts)

  aa_counts <- aaVariationToCounts(aa_var, inheritance_model = "dominant")
  test_aa_counts <- data.frame(
    ID = c("P1", "P2"),
    A_44_K = c(1, 1),
    A_44_R = c(1, 1),
    A_62_Q = c(1, 1),
    A_62_G = c(1, 1),
    stringsAsFactors = FALSE
  )
  expect_equal(aa_counts, test_aa_counts)

  aa_counts <- aaVariationToCounts(aa_var, inheritance_model = "recessive")
  test_aa_counts <- data.frame(
    ID = c("P1", "P2"),
    A_44_K = c(0, 0),
    A_44_R = c(0, 0),
    A_62_Q = c(0, 0),
    A_62_G = c(0, 0),
    stringsAsFactors = FALSE
  )
  expect_equal(aa_counts, test_aa_counts)

  expect_error(
    aaVariationToCounts(c("x", "c"), inheritance_model = "additive"),
    "aa_variation is not a data frame"
  )

  expect_error(
    aaVariationToCounts(aa_var[, -1], inheritance_model = "additive"),
    "first column of aa_variation must be named ID"
  )

  expect_error(
    aaVariationToCounts(aa_var, inheritance_model = 123),
    "inheritance_model is not a string \\(a length one character vector\\)."
  )

  expect_error(
    aaVariationToCounts(aa_var, inheritance_model = "foo"),
    "inheritance_model should be one of 'dominant', 'recessive', 'additive'"
  )
})

test_that("amino acids frequencies are calculated properly", {
  minimal_hla_calls <- data.frame(
    ID = c("P1", "P2"),
    A_1 = c("A*01:01", "A*02:01"),
    A_2 = c("A*02:01", "A*01:01"),
    stringsAsFactors = FALSE
  )
  aa_var <- hlaToAAVariation(minimal_hla_calls)[, 1:5]
  aa_freq <- getAAFrequencies(aa_var)
  test_aa_freq <- data.frame(
    aa_pos = c("A_44_K", "A_44_R", "A_62_G", "A_62_Q"),
    Freq = c(0.5, 0.5, 0.5, 0.5),
    stringsAsFactors = FALSE
  )
  expect_equal(aa_freq, test_aa_freq)

  expect_error(
    aaVariationToCounts(c("x", "c"), inheritance_model = "additive"),
    "aa_variation is not a data frame"
  )

  expect_error(
    aaVariationToCounts(aa_var[, -1], inheritance_model = "additive"),
    "first column of aa_variation must be named ID"
  )
})

test_that("hla counts table can be rreverted to hla calls", {
  hla_calls <- data.frame(ID = c("PAT1", "PAT2", "PAT3"),
                          A_1 = c("A*02:01", "A*02:01", "A*01:01"),
                          A_2 = c("A*02:01", "A*02:06", "A*24:02"),
                          B_1 = c("B*13:02", "B*15:01", "B*13:02"),
                          B_2 = c("B*40:01", "B*57:01", "B*57:01"),
                          stringsAsFactors = FALSE
  )
  hla_counts <- hlaCallsToCounts(hla_calls, inheritance_model = "additive")
  expect_equal(countsToHlaCalls(hla_counts), hla_calls)

  err_hla_counts <- hla_counts
  colnames(err_hla_counts) <- NULL
  expect_error(
    countsToHlaCalls(err_hla_counts),
    "count table has no column names"
  )

  err_hla_counts <- hla_counts
  colnames(err_hla_counts)[2] <- NA
  expect_error(
    countsToHlaCalls(err_hla_counts),
    "column names contains NA values"
  )

  err_hla_counts <- hla_counts
  colnames(err_hla_counts) <- make.names(colnames(err_hla_counts))
  expect_error(
    countsToHlaCalls(err_hla_counts),
    "counts table column names contains improperly formated HLA alleles numbers"
  )

  err_hla_counts <- hla_counts
  err_hla_counts[2, 2] <- 1.5
  expect_error(
    countsToHlaCalls(err_hla_counts),
    "counts can only take values 0, 1 or 2"
  )

  err_hla_counts <- hla_counts
  err_hla_counts[2, 2:3] <- 2
  expect_error(
    countsToHlaCalls(err_hla_counts),
    "some samples have more than two alleles per gene"
  )
})
