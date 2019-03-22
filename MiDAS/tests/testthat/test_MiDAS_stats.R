context("HLA allele statistical methods")

test_that("HLA allele associations are analyzed properly", {
  hla_calls_file <- system.file(
    "extdata", "HLAHD_output_example.txt", package = "MiDAS"
  )
  hla_calls <- readHlaCalls(hla_calls_file)
  pheno_file <- system.file("extdata", "pheno_example.txt", package = "MiDAS")
  pheno <- read.table(pheno_file, header = TRUE)
  covar_file <- system.file("extdata", "covar_example.txt", package = "MiDAS")
  covar <- read.table(covar_file, header = TRUE)
  res <- analyzeHlaAssociations(model = "coxph",
                                hla_calls,
                                pheno,
                                covar,
                                zygo = FALSE,
                                reduce_counts = FALSE,
                                correction = "BH"
  )
  load(system.file("extdata", "test_hla_analyze.Rdata", package = "MiDAS"))
  expect_equal(res, test_hla_analyze)

  expect_error(analyzeHlaAssociations(model = 1),
               "model is not a string \\(a length one character vector\\)."
  )

  expect_error(analyzeHlaAssociations(model = "foo"),
               "model foo is not implemented"
  )

  expect_error(analyzeHlaAssociations(model = "lm", hla_calls = hla_calls[, 1]), # other errors for hla_calls format are not checked here
               "hla_calls is not a data frame"
  )

  expect_error(
    analyzeHlaAssociations(model = "lm", hla_calls = hla_calls, pheno = 1),
    "pheno have to be a data frame"
  )

  expect_error(
    analyzeHlaAssociations(model = "lm",
                           hla_calls = hla_calls,
                           pheno = pheno[, 1, drop = FALSE]
    ),
    "pheno have to have at least 1 rows and 2 columns"
  )

  expect_error(
    analyzeHlaAssociations(model = "lm",
                           hla_calls = hla_calls,
                           pheno = pheno[, c(2, 2)]
    ),
    "first column in pheno must be named as first column in hla_calls"
  )

  pheno2 <- pheno
  pheno2$ID <- 1:nrow(pheno)
  expect_error(
    analyzeHlaAssociations(model = "lm",
                           hla_calls = hla_calls,
                           pheno = pheno2
    ),
    "IDs in pheno doesn't match IDs in hla_calls"
  )

  expect_error(
    analyzeHlaAssociations(model = "lm",
                           hla_calls = hla_calls,
                           pheno = pheno,
                           covar = 1
    ),
    "covar have to be a data frame"
  )

  expect_error(
    analyzeHlaAssociations(model = "lm",
                           hla_calls = hla_calls,
                           pheno = pheno,
                           covar = covar[, 1, drop = FALSE]
    ),
    "covar have to have at least 1 rows and 2 columns"
  )

  expect_error(
    analyzeHlaAssociations(model = "lm",
                           hla_calls = hla_calls,
                           pheno = pheno,
                           covar = covar[, c(2, 2)]
    ),
    "first column in covar must be named as first column in hla_calls"
  )

  covar2 <- covar
  covar2$ID <- 1:nrow(covar)
  expect_error(
    analyzeHlaAssociations(model = "lm",
                           hla_calls = hla_calls,
                           pheno = pheno,
                           covar = covar2
    ),
    "IDs in covar doesn't match IDs in hla_calls"
  )

  expect_error(
    analyzeHlaAssociations(model = "lm",
                           hla_calls = hla_calls,
                           pheno = pheno,
                           covar = covar,
                           zygo = "yes"
    ),
    "zygo is not a flag \\(a length one logical vector\\)."
  )

  expect_error(
    analyzeHlaAssociations(model = "lm",
                           hla_calls = hla_calls,
                           pheno = pheno,
                           covar = covar,
                           reduce_counts = "yes"
    ),
    "reduce_counts is not a flag \\(a length one logical vector\\)."
  )

  expect_error(
    analyzeHlaAssociations(model = "lm",
                           hla_calls = hla_calls,
                           pheno = pheno,
                           covar = covar,
                           correction = 12
    ),
    "correction is not a string \\(a length one character vector\\)."
  )

  covar <- dplyr::rename(covar, AGE_zygosity = AGE, `A*01:01` = SEX)
  expect_error(
    analyzeHlaAssociations(model = "lm",
                           hla_calls = hla_calls,
                           pheno = pheno,
                           covar = covar
    ),
    "some colnames in hla_calls and pheno and covar and zygosity are duplicated"
  )
})

test_that("HLA statistical models are defined properly", {
  hla_calls_file <- system.file(
    "extdata", "HLAHD_output_example.txt", package = "MiDAS"
  )
  hla_calls <- readHlaCalls(hla_calls_file)
  pheno_file <- system.file("extdata", "pheno_example.txt", package = "MiDAS")
  pheno <- read.table(pheno_file, header = TRUE)
  covar_file <- system.file("extdata", "covar_example.txt", package = "MiDAS")
  covar <- read.table(covar_file, header = TRUE)
  hla_data <- prepareHlaData(hla_calls, pheno, covar)

  fun <- hlaAssocModel(model = "glm",
                       response = hla_data$response[2],
                       variable = hla_data$covariate,
                       data = hla_data$data,
                       family = binomial(link = "logit")
  )
  expect_equal(
    fun,
    glm(
      OS_DIED ~ AGE + SEX,
      data = hla_data$data ,
      family = binomial(link = "logit")
    )
  )

  expect_error(hlaAssocModel(model = 1),
               "model have to be a string \\(a length one character vector\\) or a function"
  )

  expect_error(hlaAssocModel(model = "l"),
               "could not find function l"
  )

  expect_error(hlaAssocModel(model = "lm", response = 1),
               "response have to be a string or formula"
  )

  expect_error(hlaAssocModel(model = "lm",
                             response = hla_data$response[2],
                             variable = 1),
               "variable have to be a character or formula"
  )

  expect_error(hlaAssocModel(model = "lm",
                              response = hla_data$response[2],
                              variable = hla_data$covariate,
                              data = 1),
               "data is not a data.frame"
  )
})

test_that("Stepwise conditional alleles subset selection", {
  hla_calls_file <- system.file("extdata", "HLAHD_output_example.txt", package = "MiDAS")
  hla_calls <- readHlaCalls(hla_calls_file)
  pheno_file <- system.file("extdata", "pheno_example.txt", package = "MiDAS")
  pheno <- read.table(pheno_file, header = TRUE)
  covar_file <- system.file("extdata", "covar_example.txt", package = "MiDAS")
  covar <- read.table(covar_file, header = TRUE)
  hla_data <<- prepareHlaData(hla_calls, pheno, covar) # TODO there is a scope error
  object <- forwardConditionalSelection("coxph", hla_calls, pheno, covar, th = 0.02)
  test_object <- coxph(
    Surv(OS, OS_DIED) ~ AGE + SEX + `B*14:02` + `DRB1*11:01` + `DRA*01:02`,
    data = hla_data$data
  )

  expect_equal(object, test_object)

  expect_error(forwardConditionalSelection(model = 2,
                                           hla_calls = hla_calls,
                                           pheno = pheno,
                                           covar = covar,
                                           th = 0.05
               ),
               "model is not a string \\(a length one character vector\\)."
  )

  # assert tests with checkHlaCallsFormat & checkAdditionalData are omited here

  expect_error(forwardConditionalSelection(model = "coxph",
                                           hla_calls = hla_calls,
                                           pheno = pheno,
                                           covar = covar,
                                           th = "foo"
               ),
               "th is not a number \\(a length one numeric vector\\)."
  )

  expect_error(forwardConditionalSelection(model = "coxph",
                                           hla_calls = hla_calls,
                                           pheno = pheno,
                                           covar = covar,
                                           th = 0.05,
                                           keep = "yes"
              ),
              "keep is not a flag \\(a length one logical vector\\)."
  )

  expect_error(forwardConditionalSelection(model = "coxph",
                                           hla_calls = hla_calls,
                                           pheno = pheno,
                                           covar = covar,
                                           th = 0.05,
                                           rss_th = "foo"
              ),
              "rss_th is not a number \\(a length one numeric vector\\)."
  )
})
