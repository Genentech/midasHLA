context("Transforming MiDAS objects")

test_that("Amino acids variability is infered correctly", {
  hla_calls <- system.file("extdata/HLAHD_output_example.txt", package = "MiDAS")
  hla_calls <- readHlaCalls(hla_calls)
  aa_variation <- hlaToAAVariation(hla_calls, indels = TRUE, unkchar = TRUE)
  load(system.file("extdata", "test_aa_variation.Rdata", package = "MiDAS"))
  expect_equal(aa_variation, test_aa_variation)

  expect_error(hlaToAAVariation(aa_variation), "hla_calls is not a data frame")

  expect_error(hlaToAAVariation(data.frame()), "input data frame have to have at least 1 rows and 2 columns")

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
  hla_counts <- hlaCallsToCounts(hla_calls)
  load(system.file("extdata", "test_hla_counts.RData", package = "MiDAS"))
  expect_equal(hla_counts, test_hla_counts)

  expect_error(
    hlaCallsToCounts(c("A*01:01", "A*02:01")),
    "hla_calls is not a data frame"
  )
})
