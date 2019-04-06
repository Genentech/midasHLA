context("HLA allele statistical methods")

test_that("HLA allele associations are analyzed properly", {
  library("survival")

  hla_calls_file <- system.file(
    "extdata", "HLAHD_output_example.txt", package = "MiDAS"
  )
  hla_calls <- readHlaCalls(hla_calls_file)
  pheno_file <- system.file("extdata", "pheno_example.txt", package = "MiDAS")
  pheno <- read.table(pheno_file, header = TRUE)
  covar_file <- system.file("extdata", "covar_example.txt", package = "MiDAS")
  covar <- read.table(covar_file, header = TRUE)
  midas_data <- prepareHlaData(hla_calls, pheno, covar, inheritance_model = "additive")

  load(system.file("extdata", "test_hla_analyze.Rdata", package = "MiDAS"))

  res <- analyzeHlaAssociations(model = "coxph",
                                response = c("OS", "OS_DIED"),
                                covariate = c("AGE", "SEX"),
                                midas_data,
                                correction = "BH"
  )
  expect_equal(as.data.frame(res), as.data.frame(test_hla_analyze[["base"]])) # Tibble doesn't respect tollerance https://github.com/tidyverse/tibble/issues/287 or something related mby

  res <- analyzeHlaAssociations(model = "coxph",
                                response = c("OS", "OS_DIED"),
                                covariate = NULL,
                                midas_data,
                                correction = "BH"
  )
  expect_equal(
    as.data.frame(res), as.data.frame(test_hla_analyze[["null_covar"]])
  )

  expect_error(analyzeHlaAssociations(model = 1),
               "model have to be a string \\(a length one character vector\\) or a function"
  )

  expect_error(analyzeHlaAssociations(model = "foo"),
               "could not find function foo"
  )

  expect_error(
    analyzeHlaAssociations(model = "lm",
                           response = 1,
                           covariate = c("AGE", "SEX"),
                           midas_data = midas_data,
                           correction = 12
    ),
    "response is not a character vector"
  )

  expect_error(
    analyzeHlaAssociations(model = "coxph",
                           response = "OS",
                           covariate = c("AGE", "SEX"),
                           midas_data = midas_data,
                           correction = 12
    ),
    "cox survival analysis requires response to be a character vector of length 2"
  )

  expect_error(
    analyzeHlaAssociations(model = "lm",
                           response = c("OS", "OS_DIED"),
                           covariate = c("AGE", "SEX"),
                           midas_data = midas_data,
                           correction = 12
    ),
    "response is not a string \\(a length one character vector\\)."
  )

  expect_error(
    analyzeHlaAssociations(model = "lm",
                           response = "foo",
                           covariate = c("AGE", "SEX"),
                           midas_data = midas_data,
                           correction = 12
    ),
    "response variables can not be found in midas_data"
  )

  expect_error(
    analyzeHlaAssociations(model = "lm",
                           response = "OS",
                           covariate = 1,
                           midas_data = midas_data,
                           correction = 12
    ),
    "covariate have to be a character or NULL"
  )

  expect_error(
    analyzeHlaAssociations(model = "lm",
                           response = "OS",
                           covariate = "foo",
                           midas_data = midas_data,
                           correction = 12
    ),
    "covariate variables can not be found in midas_data"
  )

  expect_error(
    analyzeHlaAssociations(model = "lm",
                           response = c("OS", "OS_DIED"),
                           covariate = c("AGE", "SEX"),
                           midas_data = midas_data,
                           correction = 12
    ),
    "response is not a string \\(a length one character vector\\)."
  )

  expect_error(
    analyzeHlaAssociations(model = "lm",
                           midas_data = midas_data,
                           response = "OS",
                           covariate = c("AGE", "SEX"),
                           correction = 1
    ),
    "correction is not a string \\(a length one character vector\\)."
  )

  expect_error(
    analyzeHlaAssociations(model = "lm",
                           midas_data = midas_data,
                           response = "OS",
                           covariate = c("AGE", "SEX"),
                           exponentiate = 1
    ),
    "exponentiate is not a flag \\(a length one logical vector\\)."
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
  midas_data <- prepareHlaData(hla_calls, pheno, covar, inheritance_model = "additive")

  fun <- hlaAssocModel(model = "glm",
                       response = attr(midas_data, "response")[2],
                       variable = attr(midas_data, "covariate"),
                       data = midas_data,
                       family = binomial(link = "logit")
  )
  expect_equal(
    fun,
    glm(
      OS_DIED ~ AGE + SEX,
      data = midas_data,
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
                             response = attr(midas_data, "response")[2],
                             variable = 1),
               "variable have to be a character or formula"
  )

  expect_error(hlaAssocModel(model = "lm",
                              response = attr(midas_data, "response")[2],
                              variable = attr(midas_data, "covariate"),
                              data = 1),
               "data is not a data.frame"
  )

  expect_error(hlaAssocModel(model = "list",
                             response = attr(midas_data, "response")[2],
                             variable = attr(midas_data, "covariate"),
                             data = midas_data,
                             family = binomial(link = "logit")),
               "object returned by list doesn't have OBJECT bit set"
  )
})

test_that("Stepwise conditional alleles subset selection", {
  library("survival")

  hla_calls_file <- system.file("extdata", "HLAHD_output_example.txt", package = "MiDAS")
  hla_calls <- readHlaCalls(hla_calls_file)
  pheno_file <- system.file("extdata", "pheno_example.txt", package = "MiDAS")
  pheno <- read.table(pheno_file, header = TRUE)
  covar_file <- system.file("extdata", "covar_example.txt", package = "MiDAS")
  covar <- read.table(covar_file, header = TRUE)
  midas_data <- prepareHlaData(hla_calls, pheno, covar, inheritance_model = "additive")

  object <- forwardConditionalSelection("coxph",
                                        midas_data,
                                        response = c("OS", "OS_DIED"),
                                        covariate = c("AGE", "SEX"),
                                        th = 0.02
  )

  test_object <- coxph(
    Surv(OS, OS_DIED) ~ AGE + SEX + `B*14:02` + `DRB1*11:01` + `DRA*01:02`,
    data = midas_data
  )
  expect_equal(object, test_object)

  # test if NULL covariate passes without errors
  object <- forwardConditionalSelection("coxph",
                                        midas_data,
                                        response = c("OS", "OS_DIED"),
                                        covariate = NULL,
                                        th = 0.04
  )

  test_object <- coxph(
    Surv(OS, OS_DIED) ~ `B*14:02` + `DPB1*14:01` + `B*15:01` + `DRA*01:01` +
      `DQA1*06:01` + `F*01:01`,
    data = midas_data
  )
  expect_equal(object, test_object)


  expect_error(forwardConditionalSelection(model = 2),
               "model have to be a string \\(a length one character vector\\) or a function"
  )

  expect_error(
    forwardConditionalSelection(model = "midas_data",
                                midas_data = midas_data
    ),
               "could not find function midas_data"
  )

  expect_error(
    forwardConditionalSelection(model = "coxph",
                                midas_data = midas_data,
                                response = 1
    ),
    "response is not a character vector"
  )

  expect_error(
    forwardConditionalSelection(model = "coxph",
                                midas_data = midas_data,
                                response = "OS"
    ),
    "cox survival analysis requires response to be a character vector of length 2"
  )

  expect_error(
    forwardConditionalSelection(model = "lm",
                                midas_data = midas_data,
                                response = c("OS", "OS_DIED")
    ),
    "response is not a string \\(a length one character vector\\)."
  )

  expect_error(
    forwardConditionalSelection(model = "lm",
                                midas_data = midas_data,
                                response = "foo"
    ),
    "response variables can not be found in midas_data"
  )

  expect_error(
    forwardConditionalSelection(model = "lm",
                                midas_data = midas_data,
                                response = "OS",
                                covariate = 1

    ),
    "covariate have to be a character or NULL"
  )

  expect_error(
    forwardConditionalSelection(model = "lm",
                                midas_data = midas_data,
                                response = "OS",
                                covariate = "foo"

    ),
    "covariate variables can not be found in midas_data"
  )

  expect_error(forwardConditionalSelection(model = "lm",
                                           midas_data = midas_data,
                                           response = "OS",
                                           covariate = NULL,
                                           th = "foo"
               ),
               "th is not a number \\(a length one numeric vector\\)."
  )

  expect_error(forwardConditionalSelection(model = "lm",
                                           midas_data = midas_data,
                                           response = "OS",
                                           covariate = NULL,
                                           th = 0.05,
                                           keep = "yes"
              ),
              "keep is not a flag \\(a length one logical vector\\)."
  )

  expect_error(forwardConditionalSelection(model = "lm",
                                           midas_data = midas_data,
                                           response = "OS",
                                           covariate = NULL,
                                           th = 0.05,
                                           rss_th = "foo"
              ),
              "rss_th is not a number \\(a length one numeric vector\\)."
  )
})

test_that("HLA data is properly formatted", {
  small_hla_calls <- data.frame(ID = 1:2,
                                A_1 = c("A*01:01", "A*01:02"),
                                A_2 = c("A*01:02", "A*01:01"),
                                stringsAsFactors = FALSE
  )
  small_pheno <- data.frame(ID = 1:2, OS = c(123, 321), OS_DIED = c(0, 0))
  small_covar <- data.frame(ID = 1:2, AGE = c(23, 24))
  midas_data <- prepareHlaData(small_hla_calls, small_pheno, small_covar, inheritance_model = "additive")
  expect_equal(midas_data,
               structure(
                 data.frame(ID = 1:2,
                            "A*01:01" = c(1, 1),
                            "A*01:02" = c(1, 1),
                            OS = c(123, 321),
                            OS_DIED = c(0, 0),
                            AGE = c(23, 24),
                            stringsAsFactors = FALSE,
                            check.names = FALSE
                 ),
                 response = c("OS", "OS_DIED"),
                 covariate = "AGE",
                 alleles = c("A*01:01", "A*01:02")
               )
  )

  expect_error(
    prepareHlaData(hla_calls = "foo",
                   pheno = small_pheno,
                   covar = small_covar,
                   inheritance_model = "additive"
    ),
    "hla_calls is not a data frame"
  )

  expect_error(
    prepareHlaData(hla_calls = data.frame(),
                   pheno = small_pheno,
                   covar = small_covar,
                   inheritance_model = "additive"
    ),
    "hla_calls have to have at least 1 rows and 2 columns"
  )

  small_hla_calls_fac <- small_hla_calls
  small_hla_calls_fac$A_1 <- as.factor(small_hla_calls_fac$A_1)
  expect_error(
    prepareHlaData(hla_calls = small_hla_calls_fac,
                   pheno = small_pheno,
                   covar = small_covar,
                   inheritance_model = "additive"
    ),
    "hla_calls can't contain factors"
  )

  expect_error(
    prepareHlaData(hla_calls = small_hla_calls[, 2:3],
                   pheno = small_pheno,
                   covar = small_covar,
                   inheritance_model = "additive"
    ),
    "first column of hla_calls should specify samples id"
  )

  expect_error(
    prepareHlaData(hla_calls = small_hla_calls[, c(1, 1, 1)],
                   pheno = small_pheno,
                   covar = small_covar,
                   inheritance_model = "additive"
    ),
    "values in hla_calls doesn't follow HLA numbers specification"
  )

  expect_error(
    prepareHlaData(hla_calls = small_hla_calls,
                   pheno = "foo",
                   covar = small_covar,
                   inheritance_model = "additive"
    ),
    "pheno have to be a data frame"
  )

  bad_pheno <- small_pheno
  colnames(bad_pheno) <- LETTERS[1:ncol(bad_pheno)]
  expect_error(
    prepareHlaData(hla_calls = small_hla_calls,
                   pheno = bad_pheno,
                   covar = small_covar,
                   inheritance_model = "additive"
    ),
    "first column in pheno must be named as first column in hla_calls"
  )

  bad_pheno <- small_pheno
  bad_pheno$ID <- paste0("bad", bad_pheno$ID)
  expect_error(
    prepareHlaData(hla_calls = small_hla_calls,
                   pheno = bad_pheno,
                   covar = small_covar,
                   inheritance_model = "additive"
    ),
    "IDs in pheno doesn't match IDs in hla_calls"
  )

  expect_error(
    prepareHlaData(hla_calls = small_hla_calls,
                   pheno = small_pheno,
                   covar = "foo",
                   inheritance_model = "additive"
    ),
    "covar have to be a data frame"
  )

  expect_error(
    prepareHlaData(hla_calls = small_hla_calls,
                   pheno = small_pheno,
                   covar = small_covar[, 1, drop = FALSE],
                   inheritance_model = "additive"
    ),
    "covar have to have at least 1 rows and 2 columns"
  )

  expect_error(
    prepareHlaData(hla_calls = small_hla_calls,
                   pheno = small_pheno,
                   covar = small_covar[, c(2, 2)],
                   inheritance_model = "additive"
    ),
    "first column in covar must be named as first column in hla_calls"
  )

  bad_covar <- small_covar
  bad_covar$ID <- paste0("bad", bad_covar$ID)
  expect_error(
    prepareHlaData(hla_calls = small_hla_calls,
                   pheno = small_pheno,
                   covar = bad_covar,
                   inheritance_model = "additive"
    ),
    "IDs in covar doesn't match IDs in hla_calls"
  )


  expect_error(
    prepareHlaData(small_hla_calls,
                   pheno = small_pheno,
                   covar = small_covar,
                   inheritance_model = 1
    ),
    "inheritance_model is not a string \\(a length one character vector\\)."
  )

  expect_error(
    prepareHlaData(small_hla_calls,
                   pheno = small_pheno,
                   covar = small_covar,
                   inheritance_model = "foo"
    ),
    "inheritance_model should be one of 'dominant', 'recessive', 'additive'"
  )
})
