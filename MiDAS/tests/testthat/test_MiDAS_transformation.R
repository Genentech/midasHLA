context("Transforming MiDAS objects")

test_that("Amino acids variability is infered correctly", {
  hla_calls <- system.file("extdata/HLAHD_output_example.txt", package = "MiDAS")
  hla_calls <- readHlaCalls(hla_calls)
  aa_variation <- hlaToAAVariation(hla_calls, indels = TRUE, unkchar = TRUE)
  load(system.file("extdata", "test_aa_variation.Rdata", package = "MiDAS"))
  expect_equal(aa_variation, test_aa_variation)

  expect_error(hlaToAAVariation(aa_variation), "hla_calls is not a data frame")

  expect_error(hlaToAAVariation(data.frame()), "input data frame have to have at least 1 rows and 2 columns")

  expect_error(hlaToAAVariation(hla_calls[, -1],
                                "first column of input data frame should specify samples id"
               )
  )

  expect_error(hlaToAAVariation(hla_calls[, rep(1, 5)],
                                "values in input data frame doesn't follow HLA numbers specification"
               )
  )

  expect_error(hlaToAAVariation(hla_calls, indels = "foo"),
               "indels is not a flag \\(a length one logical vector\\)."
  )

  expect_error(hlaToAAVariation(hla_calls, unkchar = "foo"),
               "unkchar is not a flag \\(a length one logical vector\\)."
  )

  expect_error(hlaToAAVariation(hla_calls, alnpath = "foo/bar/foo/bar"),
               "Path 'foo/bar/foo/bar' does not exist"
  )
})
