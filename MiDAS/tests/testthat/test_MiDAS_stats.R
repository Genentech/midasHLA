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
    "pheno is not a data frame"
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
    "IDs in hla_calls doesn't match IDs in pheno"
  )

  expect_error(
    analyzeHlaAssociations(model = "lm",
                           hla_calls = hla_calls,
                           pheno = pheno,
                           covar = 1
    ),
    "covar have to be a data frame or NULL"
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
    "IDs in hla_calls doesn't match IDs in covar"
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
  expect_equal(hlaAssocModels(), c("coxph", "lm", "glm.logit", "glm.nb"))

  hla_calls_file <- system.file(
    "extdata", "HLAHD_output_example.txt", package = "MiDAS"
  )
  hla_calls <- readHlaCalls(hla_calls_file)
  hla_counts <- hlaCallsToCounts(hla_calls)
  pheno_file <- system.file("extdata", "pheno_example.txt", package = "MiDAS")
  pheno <- read.table(pheno_file, header = TRUE)
  covar_file <- system.file("extdata", "covar_example.txt", package = "MiDAS")
  covar <- read.table(covar_file, header = TRUE)
  data <- dplyr::left_join(hla_counts, pheno, by = "ID")
  data <- dplyr::left_join(data, covar, by = "ID")
  response <- colnames(pheno[, -1])
  covariate <- colnames(covar[, -1])

  fun <- hlaAssocModels(model = "coxph",
                 response = response,
                 covariate = covariate,
                 data = data
  )
  res <- fun("`A*01:01`")
  expect_equal(
    as.character(res$call),
    c("coxph", "as.formula(form)", "data")
  )

  response <- colnames(pheno[, 3, drop = FALSE])
  fun <- hlaAssocModels(model = "lm",
                        response = response,
                        covariate = covariate,
                        data = data
  )
  res <- fun("`A*01:01`")
  expect_equal(
    as.character(res$call),
    c("lm", "as.formula(form)", "data")
  )

  fun <- hlaAssocModels(model = "glm.logit",
                        response = response,
                        covariate = covariate,
                        data = data
  )
  res <- fun("`A*01:01`")
  expect_equal(
    as.character(res$call),
    c("glm", "as.formula(form)", "binomial(link = \"logit\")", "data")
  )

  fun <- hlaAssocModels(model = "glm.nb",
                        response = response,
                        covariate = covariate,
                        data = data
  )
  res <- fun("`A*01:01`")
  expect_equal(
    as.character(res$call)[-4],
    c("glm.nb", "as.formula(form)", "data", "log")
  )

  # check if null covariate is accepted
  fun <- hlaAssocModels(model = "glm.nb",
                        response = response,
                        covariate = NULL,
                        data = data
  )
  res <- fun("`A*01:01`")
  expect_equal(
    colnames(res$model),
    c("OS_DIED", "A*01:01")
  )

  expect_error(hlaAssocModels(model = 1),
               "model is not a string \\(a length one character vector\\)."
  )

  expect_error(hlaAssocModels(model = "lm", response = 1),
               "response is not a character vector"
  )

  expect_error(hlaAssocModels(model = "lm", response = character()),
               "response can not be empty"
  )

  expect_error(hlaAssocModels(model = "lm", response = response, covariate = 1),
               "covariate have to be a character or NULL"
  )

  expect_error(hlaAssocModels(model = "lm",
                              response = response,
                              covariate = covariate,
                              data = 1),
               "data is not a data.frame"
  )
})
