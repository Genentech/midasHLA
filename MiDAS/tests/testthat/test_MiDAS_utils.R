context("HLA allele numbers")

test_that("HLA allele numbers have proper format", {
  expect_equal(checkAlleleFormat(c("A*01", "A*01:24", "B*01:25:22",
                                   "C*050:24:55:54")
               ),
               c(TRUE, TRUE, TRUE, TRUE)
  )
  expect_equal(checkAlleleFormat(c("01", "01:24", "01:25:22",
                                   "050:24:55:54")
               ),
               c(FALSE, FALSE, FALSE, FALSE)
  )
  expect_equal(checkAlleleFormat(c("*01", "A*:22", "C*05:24:55:54:89")
               ),
               c(FALSE, FALSE, FALSE)
  )
})

test_that("HLA allele resolution is number of sets of digits * 2", {
  expect_equal(getAlleleResolution(c("A*01", "A*01:24", "B*01:25:22",
                                     "C*05:24:55:54")
               ),
               c(2, 4, 6, 8)
  )
  expect_error(getAlleleResolution("word"),
               "allele have to be a valid HLA allele number"
  )
})

test_that("Reduced HLA allele have desired resoulution", {
  expect_equal(reduceAlleleResolution(c("A*01", "A*01:24", "B*01:25:22",
                                        "C*05:24:55:54", "C*05:24:55:54N"), 2
               ),
               c("A*01", "A*01", "B*01", "C*05", "C*05:24:55:54N")
  )
  expect_error(getAlleleResolution("word"),
               "allele have to be a valid HLA allele number"
  )
  expect_error(reduceAlleleResolution("C*05:24:55:54", resolution = "four"),
               "resolution is not a count \\(a single positive integer\\)"
  )
})

test_that("HLA allels are converted to additional variables", {
  path <- system.file("extdata", "Match_allele_HLA_supertype.txt", package = "MiDAS")
  addvar <- convertAlleleToVariable(c("A*01:01", "A*02:01", "B*01", NA), dictionary = path)
  expect_equal(addvar, c("A01", "A02", NA, NA))
  dictionary <- read.table(path, header = TRUE, stringsAsFactors = FALSE)
  addvar <- convertAlleleToVariable(c("A*01:01", "A*02:01", "B*01", NA), dictionary = dictionary)
  expect_equal(addvar, c("A01", "A02", NA, NA))

  expect_error(convertAlleleToVariable(c("a", "b", "c"), dictionary = path,
                                       "allele have to be a valid HLA allele number"
  )
  )

  expect_error(convertAlleleToVariable(c("A*01", "A*02", "A*03"), dictionary = c("foo", "bar"),
                                       "dictionary have to be either path or data.frame"
  )
  )

  expect_error(
    convertAlleleToVariable(
      allele = c("A*01", "A*02", "A*03"),
      dictionary = file.path("foo", "bar")
    ),
    sprintf("Path '%s' does not exist", file.path("foo", "bar"))
  )

  expect_error(convertAlleleToVariable(c("A*01", "A*02", "A*03"), dictionary = dictionary[, 1],
                                       "match table have to consist out of two columns"
  )
  )

  expect_error(convertAlleleToVariable(c("A*01", "A*02", "A*03"), dictionary = dictionary[, c(2, 2)],
                                       "first column of match table must contain valid HLA allele numbers"
  )
  )

  expect_error(convertAlleleToVariable(c("A*01", "A*02", "A*03"), dictionary = dictionary[c(1, 1), ],
                                       "match table contains duplicated allele numbers"
  )
  )
})

test_that("HLA calls data frame have proper format", {
  file <- system.file("extdata", "HLAHD_output_example.txt", package = "MiDAS")
  hla_calls <- readHlaCalls(file)
  expect_equal(checkHlaCallsFormat(hla_calls), TRUE)

  expect_error(checkHlaCallsFormat("A"), "hla_calls is not a data frame")
  expect_error(checkHlaCallsFormat(data.frame()),
               "hla_calls have to have at least 1 rows and 2 columns"
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
               "first column of hla_calls should specify samples id"
  )

  expect_error(checkHlaCallsFormat(fake_calls[, c(1, 1, 3)]),
               "values: Sample1, Sample2, Sample3 in hla_calls doesn't follow HLA numbers specification"
  )
})

test_that("HLA allele are backquoted properly", {
  expect_equal(backquote(c("A:01:01", "A:02:01")), c("`A:01:01`", "`A:02:01`"))
  expect_error(backquote(1), "x is not a character vector")
})

context("HLA allele alignments")

test_that("Variable amino acids positions are detected properly", {
  hlaa_calls <- c("A*01:01", "A*01:02")
  hlaa_res <- 4
  hlaa_aln <- readHlaAlignments(system.file("extdata",
                                            "A_prot.txt",
                                            package = "MiDAS")
  )
  four_dig_numbers <- reduceAlleleResolution(rownames(hlaa_aln), resolution = 4)
  hlaa_aln <- hlaa_aln[!duplicated(four_dig_numbers), ]
  rownames(hlaa_aln) <- four_dig_numbers[!duplicated(four_dig_numbers)]
  hlaa_aln <- hlaa_aln[hlaa_calls, ]

  expect_equal(getVariableAAPos(hlaa_aln), c(9, 17))

  expect_error(getVariableAAPos(hlaa_calls), "alignment is not a matrix")
})

context("HLA statistical models handling")

test_that("HLA statistical models are updated properly", {
  library("survival")
  hla_calls_file <- system.file("extdata",
                                "HLAHD_output_example.txt",
                                package = "MiDAS"
  )
  hla_calls <- readHlaCalls(hla_calls_file)
  pheno_file <- system.file("extdata", "pheno_example.txt", package = "MiDAS")
  pheno <- read.table(pheno_file, header = TRUE, stringsAsFactors = FALSE)
  covar_file <- system.file("extdata", "covar_example.txt", package = "MiDAS")
  covar <- read.table(covar_file, header = TRUE, stringsAsFactors = FALSE)
  midas_data <-
    prepareMiDAS(hla_calls,
                     pheno,
                     covar,
                     analysis_type = "hla_allele",
                     inheritance_model = "additive")
  coxmod <- coxph(Surv(OS, OS_DIED) ~ 1, data = midas_data)
  coxmod$call$data <- midas_data
  coxmod_test <- coxph(Surv(OS, OS_DIED) ~ `A*01:01`, data = midas_data)
  coxmod_test$call$data <- midas_data
  expect_equal(updateModel(coxmod, "A*01:01"),
               coxmod_test
  )

  expect_error(updateModel(coxmod, 1),
               "x is not a character vector or formula"
  )

  expect_error(updateModel(coxmod, x = "A*01:01", backquote = 1),
               "backquote is not a flag \\(a length one logical vector\\)."
  )

  expect_error(updateModel(coxmod, x = "A*01:01", collapse = 1),
               "collapse is not a string \\(a length one character vector\\)."
  )
})


test_that("statistical models are statistical model", {
  object <- lm(speed ~ dist, data = cars)
  expect_equal(checkStatisticalModel(object), TRUE)

  expect_error(checkStatisticalModel(list(1)),
               "object is required to have the internal OBJECT bit set"
  )

  expect_error(updateModel(speed ~ cars),
               "object have to have an attribute 'call'"
  )

  fake_model <- list(call = list(formula = "foo"))
  class(fake_model) <- "fake"
  expect_error(updateModel(fake_model),
               "object have to be a model with defined formula"
  )

  fake_model <- list(call = list(formula = 1 ~ 1))
  class(fake_model) <- "fake"
  expect_error(updateModel(fake_model),
               "object need to have data attribute defined"
  )

  fake_model <- list(call = list(formula = 1 ~ 1, data = "bigData"))
  class(fake_model) <- "fake"
  expect_error(updateModel(fake_model),
               "object need to have data attribute defined"
  )
})

test_that("is counts or zeros", {
  expect_equal(isCountsOrZeros(c(1, 0, 2, NA)), TRUE)

  expect_error(
    assertthat::assert_that(isCountsOrZeros(c(1, 0, 2, NA, 1.5))),
    "values in c\\(1, 0, 2, NA, 1.5\\) are not counts \\(a positive integers\\) or zeros."
  )
})

test_that("is character or null", {
  expect_equal(isCharacterOrNULL(LETTERS), TRUE)
  expect_equal(isCharacterOrNULL(NULL), TRUE)

  expect_error(
    assertthat::assert_that(isCharacterOrNULL(1)),
    "1 is not a character vector or NULL."
  )
})

test_that("is number or null", {
  expect_equal(isNumberOrNULL(1), TRUE)
  expect_equal(isNumberOrNULL(NULL), TRUE)

  expect_error(
    assertthat::assert_that(isNumberOrNULL("a")),
    "\"a\" is not a number \\(a length one numeric vector\\) or NULL."
  )
})

test_that("is string or null", {
  expect_equal(isStringOrNULL("foo"), TRUE)
  expect_equal(isStringOrNULL(NULL), TRUE)

  expect_error(
    assertthat::assert_that(isStringOrNULL(1)),
    "1 is not a string \\(a length one character vector\\) or NULL."
  )
})

test_that("string matches", {
  expect_equal(stringMatches("foo", c("foo", "bar")), TRUE)

  expect_error(
    assertthat::assert_that(stringMatches("foo", c("bar", "Foo"))),
    "\"foo\" should be one of \"bar\", \"Foo\"."
  )
})

test_that("is flag or null", {
  expect_equal(isFlagOrNULL(TRUE), TRUE)
  expect_equal(isFlagOrNULL(NULL), TRUE)
  expect_equal(isFlagOrNULL(NA), FALSE)

  expect_error(
    assertthat::assert_that(isFlagOrNULL(1)),
    "1 is not a flag \\(a length one logical vector\\) or NULL."
  )
})

test_that("character maches choices", {
  expect_equal(characterMatches("foo", c("foo", "bar")), TRUE)

  expect_error(
    assertthat::assert_that(characterMatches("foo", "bar")),
    '"foo" should match values "bar".'
  )
})

test_that("is class or null", {
  expect_equal(isClassOrNULL("foo", "character"), TRUE)
  expect_equal(isClassOrNULL(NULL, "character"), TRUE)

  expect_error(
    assertthat::assert_that(isClassOrNULL("foo", "bar")),
    "\"foo\" must be an instance of \"bar\" or NULL."
  )
})

test_that("KIR haplotypes are converted to gene counts", {
  x <- c("1+3|16+3", "1+1")
  kir_hap <- kirHaplotypeToCounts(x)

  hap_dict <- system.file("extdata", "Match_kir_haplotype_gene.txt", package = "MiDAS")
  hap_dict <- read.table(hap_dict)
  hap1 <- colSums(hap_dict[c("1", "3"), ])
  hap1 <- ifelse(hap1 > 1, 1, hap1)
  hap2 <- colSums(hap_dict[c("1", "1"), ])
  hap2 <- ifelse(hap2 > 1, 1, hap2)
  test_kir_hap <- rbind(hap1, hap2)
  test_kir_hap <-
    as.data.frame(test_kir_hap,
                  optional = TRUE,
                  stringsAsFactors = FALSE)
  test_kir_hap <- cbind(haplotypes = x, test_kir_hap, stringsAsFactors = FALSE)
  rownames(test_kir_hap) <- NULL

  expect_equal(kir_hap, test_kir_hap)

  expect_error(kirHaplotypeToCounts(1), "x is not a character vector")
  expect_error(kirHaplotypeToCounts(x, hap_dict = "foo"),
               "Path 'foo' does not exist")
  expect_error(kirHaplotypeToCounts(x, binary = "yes"),
               "binary is not a flag \\(a length one logical vector\\).")
})

test_that("column names matches", {
  df <- data.frame(a = 1:5, b = 1:5)
  expect_equal(colnamesMatches(df, c("a", "b")), TRUE)

  expect_error(colnamesMatches(1:2, c("foo", "bar")), "x is not a data frame")

  expect_error(colnamesMatches(data.frame(one = 1:2), c("foo", "bar")),
               "Number of columns in data.frame\\(one = 1:2\\) must equal 2."
  )

  expect_error(
    assertthat::assert_that(colnamesMatches(df, c("foo", "bar"))),
    "Columns a, b in df should be named foo, bar"
  )
})

test_that("KIR counts have proper format", {
  file <- system.file("extdata", "KIP_output_example.txt", package = "MiDAS")
  kir_counts <- readKirCalls(file)
  expect_equal(checkKirCountsFormat(kir_counts), TRUE)

  expect_equal(checkKirCountsFormat(NULL, accept.null = TRUE), TRUE)

  expect_error(
    checkKirCountsFormat(kir_counts[, 1, drop = FALSE]),
    "kir_counts\\[, 1, drop = FALSE\\] have to have at least 1 rows and 2 columns"
  )

  fake_kir_counts <- kir_counts
  fake_kir_counts[, 1] <- as.factor(fake_kir_counts[, 1, drop = TRUE])
  expect_error(
    checkKirCountsFormat(fake_kir_counts),
    "fake_kir_counts can't contain factors"
  )

  fake_kir_counts <- kir_counts
  colnames(fake_kir_counts) <- c("FOO", colnames(fake_kir_counts)[-1])
  expect_error(
    checkKirCountsFormat(fake_kir_counts),
    "Columns FOO in kir_counts should be named ID"
  )
})

test_that("is count or null", {
  expect_equal(isCountOrNULL(1), TRUE)
  expect_equal(isCountOrNULL(NULL), TRUE)

  expect_error(
    assertthat::assert_that(isCountOrNULL(1.5)),
    "1.5 is not a count \\(a single positive integer\\) or NULL."
  )
})

test_that("is true or false", {
  expect_equal(isTRUEorFALSE(TRUE), TRUE)
  expect_equal(isTRUEorFALSE(FALSE), TRUE)
  expect_equal(isTRUEorFALSE(NA), FALSE)

  expect_error(
    assertthat::assert_that(isTRUEorFALSE(1.5)),
    "1.5 is not a flag \\(a length one logical vector\\)."
  )
})

test_that("tidy method exists", {
  expect_equal(hasTidyMethod("lm"), TRUE)
  expect_equal(hasTidyMethod("foo"), FALSE)

  expect_error(
    assertthat::assert_that(hasTidyMethod("bar")),
    "tidy function for object of class \"bar\" could not be found."
  )
})

test_that("likelihood ratio test", {
  df <- data.frame(OS = c(20, 30, 40), AGE = c(50, 60, 70))
  mod0 <- lm(OS ~ 1, data = df)
  mod1 <- lm(OS ~ AGE, data = df)
  lrt_res <- LRTest(mod0, mod1)
  expect_equal(
    lrt_res,
    data.frame(
      term = "AGE",
      dof = 1,
      logLik = 109.8401115921340078785,
      statistic = 219.680223184268015757,
      p.value = 1.062026e-49,
      stringsAsFactors = FALSE
    )
  )

  expect_error(LRTest(mod1, mod0), "variables AGE were not found in mod1")
})
