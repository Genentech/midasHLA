test_that("checkHlaCallsFormat", {
  file <- system.file("extdata", "MiDAS_tut_HLA.txt", package = "MiDAS")
  hla_calls <- readHlaCalls(file)
  expect_equal(checkHlaCallsFormat(hla_calls), TRUE)

  expect_error(checkHlaCallsFormat("A"), "hla_calls is not a data frame")
  expect_error(checkHlaCallsFormat(data.frame()),
               "data.frame\\(\\) have to have at least 1 rows and 2 columns"
  )
  hla_calls[, 1] <- as.factor(hla_calls[, 1])
  expect_error(checkHlaCallsFormat(hla_calls),
               "hla_calls can't contain factors"
  )
  fake_calls <- data.frame(ID = c("Sample1", "Sample2", "Sample3"),
                           A_1 = c("A*01", "A*02", "A*03"),
                           A_2 = c("A*01", "B*02", "C*03"),
                           stringsAsFactors = FALSE
  )
  expect_error(checkHlaCallsFormat(fake_calls[, c(2, 1, 3)]),
               "first column of fake_calls\\[, c\\(2, 1, 3\\)\\] should specify samples id"
  )

  expect_error(checkHlaCallsFormat(fake_calls[, c(1, 1, 3)]),
               "values: Sample1, Sample2, Sample3 in fake_calls\\[, c\\(1, 1, 3\\)\\] doesn't follow HLA numbers specification"
  )
})

test_that("checkKirCallsFormat", {
  expect_equal(checkKirCallsFormat(MiDAS_tut_KIR), TRUE)

  fake_kir_counts <- MiDAS_tut_KIR
  fake_kir_counts[, 1] <- as.factor(fake_kir_counts[, 1, drop = TRUE])
  expect_error(
    checkKirCallsFormat(fake_kir_counts),
    "kir_calls can't contain factors"
  )


  expect_error(
    checkKirCallsFormat(MiDAS_tut_KIR[, 1, drop = FALSE]),
    "Number of columns in kir_calls must equal 17."
  )

  fake_kir_counts <- MiDAS_tut_KIR
  colnames(fake_kir_counts) <- c("FOO", colnames(fake_kir_counts)[-1])
  expect_error(
    checkKirCallsFormat(fake_kir_counts),
    "Columns: 'FOO' in kir_calls should be named 'ID'"
  )
})

test_that("isExperimentCountsOrZeros", {
  hla_calls_file <- system.file("extdata", "HLAHD_output_example.txt", package = "MiDAS")
  hla_calls <- readHlaCalls(hla_calls_file)[1:5, 1:5]

  pheno_file <- system.file("extdata", "pheno_example.txt", package = "MiDAS")
  pheno <- read.table(pheno_file, header = TRUE, stringsAsFactors = FALSE)[1:5, ]

  midas <- prepareMiDAS(
    hla_calls = hla_calls,
    colData = pheno,
    experiment = c("hla_alleles", "hla_aa")
  )

  expect_equal(isExperimentCountsOrZeros(midas[["hla_alleles"]]), TRUE)

  expect_equal(isExperimentCountsOrZeros(midas[["hla_aa"]]), TRUE)

  expect_equal(isExperimentCountsOrZeros(matrix(runif(15), nrow = 3)), FALSE)

  expect_equal(isExperimentCountsOrZeros(LETTERS), FALSE)
})

test_that("checkStatisticalModel", {
  midas <- prepareMiDAS(
    kir_calls = MiDAS_tut_KIR,
    colData = MiDAS_tut_pheno,
    experiment = "kir_genes"
  )

  object <- lm(disease ~ term, data = midas)
  expect_equal(checkStatisticalModel(object), TRUE)

  expect_error(checkStatisticalModel(list(1)),
               "list\\(1\\) was not recognized as a fit from a model function \\(such as lm, glm and many others\\)."
  )

  expect_error(checkStatisticalModel(speed ~ cars),
               "speed ~ cars was not recognized as a fit from a model function \\(such as lm, glm and many others\\). speed ~ cars does not have 'call' attribute."
  )

  fake_model <- list(call = list(formula = "foo"))
  class(fake_model) <- "fake"
  expect_error(checkStatisticalModel(fake_model),
               "fake_model was not recognized as a fit from a model function \\(such as lm, glm and many others\\). fake_model does not have 'formula' attribute."
  )

  fake_model <- list(call = list(formula = 1 ~ 1))
  class(fake_model) <- "fake"
  expect_error(checkStatisticalModel(fake_model),
               "fake_model was not recognized as a fit from a model function \\(such as lm, glm and many others\\). fake_model does not have 'data' attribute."
  )
})

test_that("hasTidyMethod", {
  expect_equal(hasTidyMethod("lm"), TRUE)
  expect_equal(hasTidyMethod("foo"), FALSE)

  expect_error(
    assertthat::assert_that(hasTidyMethod("bar")),
    "Could not find 'tidy' function for statistical model 'bar'. Please ensure that 'tidy' for selected model is available. See the 'broom' package for more information on 'tidy' function."
  )
})

test_that("isCountsOrZeros", {
  expect_equal(isCountsOrZeros(c(1, 0, 2, NA)), TRUE)

  expect_error(
    assertthat::assert_that(isCountsOrZeros(c(1, 0, 2, NA, 1.5))),
    "values in c\\(1, 0, 2, NA, 1.5\\) are not counts \\(a positive integers\\) or zeros."
  )
})

test_that("isCharacterOrNULL", {
  expect_equal(isCharacterOrNULL(LETTERS), TRUE)
  expect_equal(isCharacterOrNULL(NULL), TRUE)

  expect_error(
    assertthat::assert_that(isCharacterOrNULL(1)),
    "1 is not a character vector or NULL."
  )
})

test_that("isNumberOrNULL", {
  expect_equal(isNumberOrNULL(1), TRUE)
  expect_equal(isNumberOrNULL(NULL), TRUE)

  expect_error(
    assertthat::assert_that(isNumberOrNULL("a")),
    "\"a\" is not a number \\(a length one numeric vector\\) or NULL."
  )
})

test_that("isStringOrNULL", {
  expect_equal(isStringOrNULL("foo"), TRUE)
  expect_equal(isStringOrNULL(NULL), TRUE)

  expect_error(
    assertthat::assert_that(isStringOrNULL(1)),
    "1 is not a string \\(a length one character vector\\) or NULL."
  )
})

test_that("stringMatches", {
  expect_equal(stringMatches("foo", c("foo", "bar")), TRUE)

  expect_error(
    assertthat::assert_that(stringMatches("foo", c("bar", "Foo"))),
    "\"foo\" should be one of \"bar\", \"Foo\"."
  )
})

test_that("isFlagOrNULL", {
  expect_equal(isFlagOrNULL(TRUE), TRUE)
  expect_equal(isFlagOrNULL(NULL), TRUE)
  expect_equal(isFlagOrNULL(NA), FALSE)

  expect_error(
    assertthat::assert_that(isFlagOrNULL(1)),
    "1 is not a flag \\(a length one logical vector\\) or NULL."
  )
})

test_that("characterMatches", {
  expect_equal(characterMatches("foo", c("foo", "bar")), TRUE)

  expect_error(
    assertthat::assert_that(characterMatches("foo", "bar")),
    '"foo" should match values "bar".'
  )
})

test_that("isClassOrNULL", {
  expect_equal(isClassOrNULL("foo", "character"), TRUE)
  expect_equal(isClassOrNULL(NULL, "character"), TRUE)

  expect_error(
    assertthat::assert_that(isClassOrNULL("foo", "bar")),
    "\"foo\" must be an instance of \"bar\" or NULL."
  )
})

test_that("colnamesMatches", {
  df <- data.frame(a = 1:5, b = 1:5)
  expect_equal(colnamesMatches(df, c("a", "b")), TRUE)

  expect_error(colnamesMatches(1:2, c("foo", "bar")), "x is not a data frame")

  expect_error(colnamesMatches(data.frame(one = 1:2), c("foo", "bar")),
               "Number of columns in data.frame\\(one = 1:2\\) must equal 2."
  )

  expect_error(
    assertthat::assert_that(colnamesMatches(df, c("foo", "bar"))),
    "Columns: 'a', 'b' in df should be named 'foo', 'bar'"
  )
})

test_that("isCountOrNULL", {
  expect_equal(isCountOrNULL(1), TRUE)
  expect_equal(isCountOrNULL(NULL), TRUE)

  expect_error(
    assertthat::assert_that(isCountOrNULL(1.5)),
    "1.5 is not a count \\(a single positive integer\\) or NULL."
  )
})

test_that("isTRUEorFALSE", {
  expect_equal(isTRUEorFALSE(TRUE), TRUE)
  expect_equal(isTRUEorFALSE(FALSE), TRUE)
  expect_equal(isTRUEorFALSE(NA), FALSE)

  expect_error(
    assertthat::assert_that(isTRUEorFALSE(1.5)),
    "1.5 is not a flag \\(a length one logical vector\\)."
  )
})

test_that("objectHasPlaceholder", {
  object <- lm(speed ~ dist, data = cars)
  expect_equal(objectHasPlaceholder(object, "dist"), TRUE)
  expect_equal(objectHasPlaceholder(object, "foo"), FALSE)

  expect_error(
    assertthat::assert_that(objectHasPlaceholder(object, "foo")),
    "placeholder 'foo' could not be found in object's formula"
  )
})

test_that("checkColDataFormat", {
  pheno <- data.frame(
    ID = 1:5,
    letter = LETTERS[1:5]
  )

  expect_equal(checkColDataFormat(pheno), TRUE)

  expect_error(
    checkColDataFormat(LETTERS),
    "LETTERS have to be a data frame"
  )

  expect_error(
    checkColDataFormat(data.frame()),
    "data.frame\\(\\) have to have at least 1 row and 2 columns"
  )

  expect_error(
    checkColDataFormat(pheno[, 2, drop = FALSE]),
    "pheno\\[, 2, drop = FALSE\\] have to have at least 1 row and 2 columns"
  )
})

test_that("isClass", {
  expect_equal(isClass("foo", "character"), TRUE)

  expect_error(
    assertthat::assert_that(isClassOrNULL("foo", "bar")),
    "\"foo\" must be an instance of \"bar\"."
  )
})

test_that("Frequency cutoffs validation", {
  # lower_frequency_cutof must be a number
  lower_frequency_cutoff <- "foo"
  upper_frequency_cutoff <- 0.5
  expect_error(
    validateFrequencyCutoffs(lower_frequency_cutoff, upper_frequency_cutoff),
    "lower_frequency_cutoff is not a number \\(a length one numeric vector\\)."
  )

  # lower_frequency_cutof must be positive
  lower_frequency_cutoff <- -1
  upper_frequency_cutoff <- 0.5
  expect_error(
    validateFrequencyCutoffs(lower_frequency_cutoff, upper_frequency_cutoff),
    "lower_frequency_cutoff must be a number greater than 0."
  )

  # upper_frequency_cutoff must be a number
  lower_frequency_cutoff <- 0.5
  upper_frequency_cutoff <- "foo"
  expect_error(
    validateFrequencyCutoffs(lower_frequency_cutoff, upper_frequency_cutoff),
    "upper_frequency_cutoff is not a number \\(a length one numeric vector\\)."
  )

  # upper_frequency_cutoff must be positive
  lower_frequency_cutoff <- 0
  upper_frequency_cutoff <- -1
  expect_error(
    validateFrequencyCutoffs(lower_frequency_cutoff, upper_frequency_cutoff),
    "upper_frequency_cutoff must be a number greater than 0."
  )

  # lower_frequency_cutoff is lower than upper_frequency_cutoff
  lower_frequency_cutoff <- 5
  upper_frequency_cutoff <- 1
  expect_error(
    validateFrequencyCutoffs(lower_frequency_cutoff, upper_frequency_cutoff),
    "lower_frequency_cutoff cannot be higher than upper_frequency_cutoff."
  )

  # Both lower_frequency_cutoff and upper_frequency_cutoff have to be either frequencies or counts
  lower_frequency_cutoff <- 0.5
  upper_frequency_cutoff <- 2
  expect_error(
    validateFrequencyCutoffs(lower_frequency_cutoff, upper_frequency_cutoff),
    "Both lower_frequency_cutoff and upper_frequency_cutoff have to be either frequencies or counts."
  )
})
