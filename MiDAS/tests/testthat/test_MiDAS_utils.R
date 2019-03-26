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
  expect_error(checkAlleleFormat(1), "allele is not a character vector")
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
  path <- system.file("extdata", "Match_4digit_supertype.txt", package = "MiDAS")
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
               "values in hla_calls doesn't follow HLA numbers specification"
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
  pheno <- read.table(pheno_file, header = TRUE)
  covar_file <- system.file("extdata", "covar_example.txt", package = "MiDAS")
  covar <- read.table(covar_file, header = TRUE)
  hla_data <- prepareHlaData(hla_calls, pheno, covar)
  coxmod <- coxph(Surv(OS, OS_DIED) ~ 1, data = hla_data$data)

  expect_equal(updateModel(coxmod, "A*01:01"),
               coxph(Surv(OS, OS_DIED) ~ `A*01:01`, data = hla_data$data)
  )

  expect_error(updateModel(list(1), "A*01:01"),
               "object have to have the internal OBJECT bit set"
  )

  expect_error(updateModel(data.frame(), "A*01:01"),
               "object have to have an attribue 'call'"
  )

  fake_model <- list(call = list(formula = "foo"))
  class(fake_model) <- "fake"
  expect_error(updateModel(fake_model, "A*01:01"),
               "object have to be a model with defined formula"
  )

  expect_error(updateModel(coxmod, 1),
               "x have to be a string \\(a length one character vector\\) or formula"
  )

  # expect_error(updateModel(coxmod, "foo"),
  #             "variable foo could not be found in object data"
  # ) TODO
})
