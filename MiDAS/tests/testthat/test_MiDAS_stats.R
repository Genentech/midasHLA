context("HLA allele statistical methods")

test_that("HLA allele associations are analyzed properly", {
  hla_calls_file <-
    system.file("extdata", "HLAHD_output_example.txt", package = "MiDAS")
  hla_calls <- readHlaCalls(hla_calls_file)
  pheno_file <- system.file("extdata", "pheno_example.txt", package = "MiDAS")
  pheno <- read.table(pheno_file, header = TRUE, stringsAsFactors = FALSE)
  covar_file <- system.file("extdata", "covar_example.txt", package = "MiDAS")
  covar <- read.table(covar_file, header = TRUE, stringsAsFactors = FALSE)
  coldata <- dplyr::left_join(pheno, covar, by = "ID")
  midas <-
    prepareMiDAS(hla_calls,
                 colData = coldata,
                 analysis_type = "hla_allele",
                 inheritance_model = "additive")

  midas_data <- midasToWide(midas, analysis_type = "hla_allele")
  object <- lm(OS_DIED ~ AGE + SEX + term, data = midas_data)

  res <- analyzeAssociations(object,
                             variables = c("A*01:01", "A*02:01"),
                             correction = "BH"
  )

  test_res <- list(
    lm(OS_DIED ~ AGE + SEX + `A*01:01`, data = midas_data),
    lm(OS_DIED ~ AGE + SEX + `A*02:01`, data = midas_data)
  )
  test_res <- do.call("rbind", lapply(test_res, tidy))
  test_res$term <- gsub("`", "", test_res$term)
  test_res <- test_res[test_res$term %in% c("A*01:01", "A*02:01"), ]
  test_res$p.adjusted <- p.adjust(test_res$p.value, "BH")

  expect_equal(as.data.frame(res), as.data.frame(test_res)) # Tibble doesn't respect tollerance https://github.com/tidyverse/tibble/issues/287 or something related mby

  # Tests for checkStatisticalModel errors are ommitted here

  expect_error(analyzeAssociations(object, variables = 1),
               "variables is not a character vector"
  )

  expect_error(analyzeAssociations(object, variables = "thief"),
               "thief can not be found in object data"
  )

  expect_error(
    analyzeAssociations(object, variables = "A*01:01", placeholder = 1),
    "placeholder is not a string \\(a length one character vector\\)."
  )

  expect_error(
    analyzeAssociations(object, variables = "A*01:01", placeholder = "foo"),
    "placeholder 'foo' could not be found in object's formula"
  )

  expect_error(
    analyzeAssociations(object, variables = "A*01:01", n_correction = 1.5),
    "n_correction is not a count \\(a single positive integer\\) or NULL."
  )

  expect_error(
    analyzeAssociations(object, variables = "A*01:01", correction = 1),
    "correction is not a string \\(a length one character vector\\)."
  )

  expect_error(
    analyzeAssociations(object, variables = "A*01:01", exponentiate = 1),
    "exponentiate is not a flag \\(a length one logical vector\\)."
  )

  expect_error(
    analyzeAssociations(
      object,
      variables = c("A*01:01", "A*02:01"),
      n_correction = 1
    ),
    "n_correction must be at least 2."
  )
})

test_that("Stepwise conditional alleles subset selection", {
  library("survival")

  hla_calls_file <-
    system.file("extdata", "HLAHD_output_example.txt", package = "MiDAS")
  hla_calls <- readHlaCalls(hla_calls_file)
  pheno_file <- system.file("extdata", "pheno_example.txt", package = "MiDAS")
  pheno <- read.table(pheno_file, header = TRUE, stringsAsFactors = FALSE)
  covar_file <- system.file("extdata", "covar_example.txt", package = "MiDAS")
  covar <- read.table(covar_file, header = TRUE, stringsAsFactors = FALSE)
  coldata <- dplyr::left_join(pheno, covar, by = "ID")
  midas <-
    prepareMiDAS(
      hla_calls,
      colData = coldata,
      analysis_type = "hla_allele",
      inheritance_model = "additive"
    )
  midas_data <- midasToWide(midas, analysis_type = "hla_allele")

  object <- coxph(Surv(OS, OS_DIED) ~ AGE + SEX + term, data = midas_data)

  # keep = FALSE
  res <- analyzeConditionalAssociations(object,
                                        variables = c("B*14:02", "DRB1*11:01"),
                                        th = 0.05,
                                        keep = FALSE)
  res <- rapply(res, classes = "numeric", how = "replace", round, digits = 3)

  test_res <- tibble(term = c("B*14:02", "DRB1*11:01"),
                     estimate = c(3.72, 2.612),
                     std.error = c(1.59, 1.069),
                     statistic = c(2.339, 2.442),
                     p.value = c(0.019, 0.015),
                     conf.low = c(0.603, 0.516),
                     conf.high = c(6.838, 4.707),
                     p.adjusted = c(0.039, 0.015),
                     covariates = c("", "B*14:02")
  )

  expect_equal(res, test_res)

  # keep = TRUE
  res <- analyzeConditionalAssociations(object,
                                        variables = c("B*14:02", "DRB1*11:01"),
                                        th = 0.05,
                                        keep = TRUE)
  res <- rapply(res, classes = "numeric", how = "replace", round, digits = 3)
  test_res <- list(
    tibble(term = c("B*14:02", "DRB1*11:01"),
           estimate = c(3.72, 1.956),
           std.error = c(1.59, 0.960),
           statistic = c(2.339, 2.038),
           p.value = c(0.019, 0.042),
           conf.low = c(0.603, 0.075),
           conf.high = c(6.838, 3.838),
           p.adjusted = c(0.039, 0.083),
           covariates = c("", "")
    ),
    tibble(term = "DRB1*11:01",
           estimate = 2.612,
           std.error = 1.069,
           statistic = 2.442,
           p.value = 0.015,
           conf.low = 0.516,
           conf.high = 4.707,
           p.adjusted = 0.015,
           covariates = "B*14:02"
    )
  )
  expect_equal(res, test_res)

  # Tests for checkStatisticalModel errors are ommitted here

  expect_error(
    analyzeConditionalAssociations(object, variables = 1),
    "variables is not a character vector"
  )

  expect_error(analyzeConditionalAssociations(object, variables = "thief"),
               "thief can not be found in object data"
  )

  expect_error(
    analyzeConditionalAssociations(object, variables = "A*01:01", placeholder = 1),
    "placeholder is not a string \\(a length one character vector\\)."
  )

  expect_error(
    analyzeConditionalAssociations(object, variables = "A*01:01", placeholder = "foo"),
    "placeholder 'foo' could not be found in object's formula"
  )

  expect_error(
    analyzeConditionalAssociations(object, variables =  "A*01:01", correction = 1),
    "correction is not a string \\(a length one character vector\\)."
  )

  expect_error(
    analyzeConditionalAssociations(object, variables =  "A*01:01", n_correction = "foo"),
    "n_correction is not a count \\(a single positive integer\\) or NULL."
  )

  expect_error(
    analyzeConditionalAssociations(object, variables =  "A*01:01", th = "bar"),
    "th is not a number \\(a length one numeric vector\\)."
  )

  expect_error(
    analyzeConditionalAssociations(
      object,
      variables =  "A*01:01",
      th = 0.05,
      rss_th = "foo"
    ),
    "rss_th is not a number \\(a length one numeric vector\\)."
  )

  expect_error(
    analyzeConditionalAssociations(
      object,
      variables =  "A*01:01",
      th = 0.05,
      exponentiate = "yes"
    ),
    "exponentiate is not a flag \\(a length one logical vector\\)."
  )

  expect_error(
    analyzeConditionalAssociations(
      object,
      variables =  c("B*14:02", "DRB1*11:01"),
      th = 1,
      n_correction = 1
    ),
    "n_correction must be at least 2."
  )
})

test_that("MiDAS associations are analyzed properly", {
  hla_calls_file <-
    system.file("extdata", "HLAHD_output_example.txt", package = "MiDAS")
  hla_calls <- readHlaCalls(hla_calls_file)
  kir_file <-
    system.file("extdata", "KIP_output_example.txt", package = "MiDAS")
  kir_calls <- readKirCalls(kir_file, counts = TRUE)
  pheno_file <- system.file("extdata", "pheno_example.txt", package = "MiDAS")
  pheno <- read.table(pheno_file, header = TRUE, stringsAsFactors = FALSE)
  covar_file <- system.file("extdata", "covar_example.txt", package = "MiDAS")
  covar <- read.table(covar_file, header = TRUE, stringsAsFactors = FALSE)
  coldata <- dplyr::left_join(pheno, covar, by = "ID")

  midas <-
    prepareMiDAS(
      hla_calls = hla_calls,
      kir_call = kir_calls,
      colData = coldata,
      analysis_type = c(
        "hla_allele",
        "aa_level",
        "allele_g_group",
        "allele_supertype",
        "allele_group",
        "kir_genes",
        "hla_kir_interactions"
      ),
      inheritance_model = "additive"
    )



  #
  mode <- "linear"
  analysis_type_choice <-
    c(
      "hla_allele",
      "aa_level",
      "allele_g_group",
      "allele_supertype",
      "allele_group",
      "kir_genes",
      "hla_kir_interactions"
    )
  for (analysis_type in analysis_type_choice) {
    object <- lm(OS_DIED ~ AGE + SEX + term, data = midas)
    res <- runMiDAS(object,
                    mode = mode,
                    analysis_type = analysis_type,
                    exponentiate = FALSE
    )

    midas_data <- midasToWide(midas, analysis_type = analysis_type)
    object$call$data <- midas_data
    test_variables <- rownames(midas[[analysis_type]])
    test_res <-
      analyzeAssociations(object, variables = test_variables, exponentiate = FALSE)
    variables_freq <-
      MiDAS:::runMiDASGetVarsFreq(
        midas = midas,
        analysis_type = analysis_type,
        test_covar = all.vars(formula(object))[1]
      )
    test_res <- dplyr::left_join(test_res, variables_freq, by = "term")

    term_name <- switch (analysis_type,
                         "hla_allele" = "allele",
                         "aa_level" = "aa",
                         "expression_level" = "allele",
                         "allele_g_group" = "g.group",
                         "allele_supertype" = "supertype",
                         "allele_group" = "allele.group",
                         "kir_genes" = "kir.gene",
                         "hla_kir_interactions" = "hla.kir.interaction",
                         "term"
    )
    test_res <- dplyr::rename(test_res, !!term_name := term)


    expect_equal(as.data.frame(res), as.data.frame(test_res)) # Tibble doesn't respect tollerance https://github.com/tidyverse/tibble/issues/287 or something related mby
  }

  #
  mode <- "conditional"
  th <- 0.1
  keep <- FALSE
  analysis_type_choice <-
    c(
      "hla_allele",
      "aa_level",
      "allele_g_group",
      "allele_supertype",
      "allele_group",
      "kir_genes",
      "hla_kir_interactions"
    )
  for (analysis_type in analysis_type_choice) {
    object <- lm(OS_DIED ~ AGE + SEX + term, data = midas)
    res <- runMiDAS(object,
                    mode = mode,
                    analysis_type = analysis_type,
                    exponentiate = FALSE,
                    th = th,
                    keep = keep
    )

    midas_data <- midasToWide(midas, analysis_type = analysis_type)
    object$call$data <- midas_data
    test_variables <- rownames(midas[[analysis_type]])
    test_res <-
      analyzeConditionalAssociations(
        object,
        variables = test_variables,
        exponentiate = FALSE,
        th = th,
        keep = keep
      )
    variables_freq <-
      MiDAS:::runMiDASGetVarsFreq(
        midas = midas,
        analysis_type = analysis_type,
        test_covar = all.vars(formula(object))[1]
      )
    test_res <- dplyr::left_join(test_res, variables_freq, by = "term")

    term_name <- switch (analysis_type,
                         "hla_allele" = "allele",
                         "aa_level" = "aa",
                         "expression_level" = "allele",
                         "allele_g_group" = "g.group",
                         "allele_supertype" = "supertype",
                         "allele_group" = "allele.group",
                         "kir_genes" = "kir.gene",
                         "hla_kir_interactions" = "hla.kir.interaction",
                         "term"
    )
    test_res <- dplyr::rename(test_res, !!term_name := term)


    expect_equal(as.data.frame(res), as.data.frame(test_res))
  }

  #
  mode <- "conditional"
  th <- 0.1
  keep <- TRUE
  analysis_type_choice <-
    c(
      "hla_allele",
      "aa_level",
      "allele_g_group",
      "allele_supertype",
      "allele_group",
      "kir_genes",
      "hla_kir_interactions"
    )
  for (analysis_type in analysis_type_choice) {
    object <- lm(OS_DIED ~ AGE + SEX + term, data = midas)
    res <- runMiDAS(object,
                    mode = mode,
                    analysis_type = analysis_type,
                    exponentiate = FALSE,
                    th = th,
                    keep = keep
    )

    midas_data <- midasToWide(midas, analysis_type = analysis_type)
    object$call$data <- midas_data
    test_variables <- rownames(midas[[analysis_type]])
    test_res <-
      analyzeConditionalAssociations(
        object,
        variables = test_variables,
        exponentiate = FALSE,
        th = th,
        keep = keep
      )
    variables_freq <-
      MiDAS:::runMiDASGetVarsFreq(
        midas = midas,
        analysis_type = analysis_type,
        test_covar = all.vars(formula(object))[1]
      )
    test_res <- lapply(test_res, function(x) dplyr::left_join(x, variables_freq, by = "term"))

    term_name <- switch (analysis_type,
                         "hla_allele" = "allele",
                         "aa_level" = "aa",
                         "expression_level" = "allele",
                         "allele_g_group" = "g.group",
                         "allele_supertype" = "supertype",
                         "allele_group" = "allele.group",
                         "kir_genes" = "kir.gene",
                         "hla_kir_interactions" = "hla.kir.interaction",
                         "term"
    )
    test_res <- lapply(test_res, function(x) dplyr::rename(x, !!term_name := term))

    expect_equal(lapply(res, as.data.frame), lapply(test_res, as.data.frame))
  }

  #
  object <- lm(OS_DIED ~ AGE + SEX + term, data = midas)

  expect_error(runMiDAS(list()),
               "object is required to have the internal OBJECT bit set"
  )

  fake_object <- list(call = list(formula = 1 ~ 1, data = 1:5))
  class(fake_object) <- "foo"
  expect_error(runMiDAS(fake_object),
               "tidy function for object of class class\\(object\\)\\[1L\\] could not be found." #TODO
  )

  fake_object <- object
  fake_midas <- midas
  class(fake_midas) <- structure("MultiAssayExperiment", package = "MultiAssayExperiment")
  fake_object$call$data <- fake_midas
  expect_error(runMiDAS(fake_object),
               "object_details\\$data must be an instance of \"MiDAS\"." #TODO
  )

  # validObject test is ommited here

  fake_object <- object
  fake_object$call$formula <- OS ~ AGE + SEX
  expect_error(
    runMiDAS(fake_object),
    "placeholder 'term' could not be found in object's formula"
  )

  expect_error(runMiDAS(object, mode = 1),
               "mode is not a string \\(a length one character vector\\)."
  )

  expect_error(runMiDAS(object, mode = "foo"),
               "mode should be one of \"linear\", \"conditional\"."
  )

  expect_error(runMiDAS(object, mode = "linear", analysis_type = 1),
               "analysis_type is not a string \\(a length one character vector\\)."
  )

  expect_error(runMiDAS(object, mode = "linear", analysis_type = "foo"),
               "analysis_type should be one of \"hla_allele\", \"aa_level\", \"allele_g_group\", \"allele_supertype\", \"allele_group\", \"kir_genes\", \"hla_kir_interactions\"."
  )

  expect_error(runMiDAS(object, mode = "linear", analysis_type = "hla_allele", correction = 1),
               "correction is not a string \\(a length one character vector\\)."
  )

  expect_error(runMiDAS(object, mode = "linear", analysis_type = "hla_allele", n_correction = "foo"),
               "n_correction is not a count \\(a single positive integer\\) or NULL."
  )

  expect_error(runMiDAS(object, mode = "linear", analysis_type = "hla_allele", exponentiate = "foo"),
               "exponentiate is not a flag \\(a length one logical vector\\)."
  )
})

test_that("amino acid omnibus test works fine", {
  hla_calls_file <- system.file("extdata", "HLAHD_output_example.txt", package = "MiDAS")
  hla_calls <- readHlaCalls(hla_calls_file)
  pheno_file <- system.file("extdata", "pheno_example.txt", package = "MiDAS")
  pheno <- read.table(pheno_file, header = TRUE)
  covar_file <- system.file("extdata", "covar_example.txt", package = "MiDAS")
  covar <- read.table(covar_file, header = TRUE)
  midas_data <- prepareMiDAS(hla_calls, pheno, covar, analysis_type = "aa_level")
  object <- lm(OS ~ AGE + SEX + term, data = midas_data)
  omnibus_res <- aaPosOmnibusTest(object, aa_pos = c("B_11", "E_107", "A_246"))

  obj_B11 <- lm(OS ~ AGE + SEX + B_11_A + B_11_S + term, data = midas_data)
  obj_E107 <- lm(OS ~ AGE + SEX + E_107_R + E_107_G + term, data = midas_data)
  obj_A246 <- lm(OS ~ AGE + SEX + A_246_A + A_246_S + term, data = midas_data)
  LRT <- lapply(list(obj_B11, obj_E107, obj_A246), LRTest, mod0 = object)
  omnibus_res_test <- data.frame(
    aa_pos = c("B_11", "E_107", "A_246"),
    residues = c("A, S", "R, G", "A, S"),
    d.f. = c(1, 1, 1),
    statistic = sapply(LRT, `[[`, "statistic"),
    p.value = sapply(LRT, `[[`, "p.value"),
    p.adjusted = p.adjust(sapply(LRT, `[[`, "p.value"), method = "bonferroni"),
    stringsAsFactors = FALSE
  )
  expect_equal(omnibus_res, omnibus_res_test)

  # Tests for checkStatisticalModel errors are ommitted here

  expect_error(aaPosOmnibusTest(object, aa_pos = 1:5),
               "aa_pos is not a character vector"
  )

  expect_error(
    aaPosOmnibusTest(object, c("B_11", "E_107", "A_246"), correction = 1),
    "correction is not a string \\(a length one character vector\\)."
  )

  expect_error(
    aaPosOmnibusTest(object, c("B_11", "E_107", "A_246"), n_correction = 1.5),
    "n_correction is not a count \\(a single positive integer\\) or NULL."
  )

  expect_error(
    aaPosOmnibusTest(object, "FOO_2"),
    "amino acid position FOO_2 could not be found."
  )

  expect_error(
    aaPosOmnibusTest(object, c("B_11", "E_107", "A_246"), n_correction = 1),
    "n_correction must be at least 3."
  )
})
