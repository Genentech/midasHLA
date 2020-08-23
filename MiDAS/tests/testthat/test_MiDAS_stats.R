context("stats")

test_that("analyzeAssociations", {
  midas <-
    prepareMiDAS(hla_calls = MiDAS_tut_HLA,
                 colData = MiDAS_tut_pheno,
                 experiment = "hla_alleles")

  midas_data <- midasToWide(midas, experiment = "hla_alleles")
  object <- lm(disease ~ term, data = midas_data)

  res <- analyzeAssociations(object,
                             variables = c("A*01:01", "A*02:01"),
                             correction = "BH"
  )

  test_res <- list(
    lm(disease ~ `A*01:01`, data = midas),
    lm(disease ~ `A*02:01`, data = midas)
  )
  test_res <- do.call("rbind", lapply(test_res, tidy, conf.int = TRUE))
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

test_that("analyzeConditionalAssociations", {
  midas <-
    prepareMiDAS(hla_calls = MiDAS_tut_HLA,
                 colData = MiDAS_tut_pheno,
                 experiment = "hla_alleles")
  midas_data <- midasToWide(midas, experiment = "hla_alleles")

  object <- lm(disease ~ term, data = midas_data)

  # keep = FALSE
  res <- analyzeConditionalAssociations(object,
                                        variables = c("B*14:02", "DRB1*11:01"),
                                        th = 0.05,
                                        keep = FALSE)
  res <- rapply(res, classes = "numeric", how = "replace", round, digits = 3)

  test_res <- tibble(term = c("B*14:02", "DRB1*11:01"),
                     estimate = c(0.196, -0.092),
                     std.error = c(0.072, 0.044),
                     statistic = c(2.739, -2.099),
                     p.value = c(0.006, 0.036),
                     conf.low = c(0.056, -0.177),
                     conf.high = c(0.337, -0.006),
                     p.adjusted = c(0.013, 0.036),
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
           estimate = c(0.196, -0.096),
           std.error = c(0.072, 0.044),
           statistic = c(2.739, -2.194),
           p.value =  c(0.006, 0.028),
           conf.low = c(0.056, -0.182),
           conf.high = c(0.337, -0.01),
           p.adjusted = c(0.013, 0.057),
           covariates = c("", "")
    ),
    tibble(term = "DRB1*11:01",
           estimate = -0.092,
           std.error = 0.044,
           statistic = -2.099,
           p.value = 0.036,
           conf.low = -0.177,
           conf.high = -0.006,
           p.adjusted = 0.036,
           covariates = "B*14:02"
    )
  )
  expect_equal(res, test_res)

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

test_that("omnibusTest", {
  midas <-
    prepareMiDAS(
      hla_calls = MiDAS_tut_HLA,
      colData = MiDAS_tut_pheno,
      experiment = "hla_aa"
    )
  midas_data <- midasToWide(midas, experiment = "hla_aa")
  object <- lm(disease ~ term, data = midas_data)
  omnibus_groups <- list(
    A_77 = c("A_77_D", "A_77_N", "A_77_S"),
    A_79 = c("A_79_G", "A_79_R")
  )
  omnibus_res <- omnibusTest(object, omnibus_groups)

  obj_A77 <- lm(disease ~ A_77_D + A_77_N + A_77_S, data = midas_data)
  obj_A79 <- lm(disease ~ A_79_G + A_79_R , data = midas_data)
  LRT <-
    lapply(list(obj_A77, obj_A79),
           LRTest,
           mod0 = lm(disease ~ 1, data = midas_data))
  omnibus_res_test <- data.frame(
    group = c("A_77", "A_79"),
    term = c("A_77_D, A_77_N, A_77_S", "A_79_G, A_79_R"),
    dof = sapply(LRT, `[[`, "dof"),
    logLik = sapply(LRT, `[[`, "logLik"),
    statistic = sapply(LRT, `[[`, "statistic"),
    p.value = sapply(LRT, `[[`, "p.value"),
    p.adjusted = p.adjust(sapply(LRT, `[[`, "p.value"), method = "bonferroni"),
    stringsAsFactors = FALSE
  )
  expect_equal(omnibus_res, omnibus_res_test)
})

test_that("runMiDAS", {
  midas <-
    prepareMiDAS(
      hla_calls = MiDAS_tut_HLA,
      kir_call = MiDAS_tut_KIR,
      colData = MiDAS_tut_pheno,
      experiment = c(
        "hla_alleles",
        "hla_aa",
        "hla_g_groups",
        "hla_supertypes",
        "hla_NK_ligands",
        "kir_genes",
        "hla_kir_interactions",
        "hla_divergence"
      )
    )

  # linear
  conditional = FALSE
  omnibus = FALSE
  experiment_choice <-
    c(
      "hla_alleles",
      "hla_aa",
      "hla_g_groups",
      "hla_supertypes",
      "hla_NK_ligands",
      "kir_genes",
      "hla_kir_interactions",
      "hla_divergence"
    )
  for (experiment in experiment_choice) {
    object <- lm(disease ~ term, data = midas)
    res <- runMiDAS(object,
                    conditional = conditional,
                    omnibus = omnibus,
                    experiment = experiment,
                    exponentiate = FALSE
    )

    midas_data <- midasToWide(midas, experiment = experiment)
    object$call$data <- midas_data
    test_variables <- rownames(midas[[experiment]])
    test_res <-
      analyzeAssociations(object, variables = test_variables, exponentiate = FALSE)

    ex <- midas[[experiment]]
    if (is(ex, "SummarizedExperiment")) {ex <- assay(ex)}
    if (typeof(ex) == "integer") {
      variables_freq <-
        MiDAS:::runMiDASGetVarsFreq(
          midas = midas,
          experiment = experiment,
          test_covar = all.vars(formula(object))[1]
        )
      test_res <- dplyr::left_join(test_res, variables_freq, by = "term")
    }

    term_name <- switch (experiment,
                         "hla_alleles" = "allele",
                         "hla_aa" = "aa",
                         "expression_level" = "allele",
                         "hla_g_groups" = "g.group",
                         "hla_supertypes" = "supertype",
                         "hla_NK_ligands" = "allele.group",
                         "kir_genes" = "kir.gene",
                         "hla_kir_interactions" = "hla.kir.interaction",
                         "term"
    )
    test_res <- dplyr::rename(test_res, !!term_name := term)
    test_res <- dplyr::arrange(test_res, p.value)
    test_res <- test_res[, colnames(res)] # order columns


    expect_equal(as.data.frame(res), as.data.frame(test_res)) # Tibble doesn't respect tollerance https://github.com/tidyverse/tibble/issues/287 or something related mby
  }

  # conditional
  conditional <- TRUE
  omnibus <- FALSE
  th <- 0.2
  keep <- FALSE
  experiment_choice <-
    c(
      "hla_alleles",
      "hla_aa",
      "hla_g_groups",
      "hla_supertypes",
      "hla_NK_ligands",
      "kir_genes",
      "hla_kir_interactions",
      "hla_divergence"
    )
  for (experiment in experiment_choice) {
    object <- lm(disease ~ term, data = midas)
    res <- runMiDAS(object,
                    conditional = conditional,
                    omnibus = omnibus,
                    experiment = experiment,
                    exponentiate = FALSE,
                    th = th,
                    keep = keep
    )

    midas_data <- midasToWide(midas, experiment = experiment)
    object$call$data <- midas_data
    test_variables <- rownames(midas[[experiment]])
    test_res <-
      analyzeConditionalAssociations(
        object,
        variables = test_variables,
        exponentiate = FALSE,
        th = th,
        keep = keep
      )

    ex <- midas[[experiment]]
    if (is(ex, "SummarizedExperiment")) {ex <- assay(ex)}
    if (typeof(ex) == "integer") {
      variables_freq <-
        MiDAS:::runMiDASGetVarsFreq(
          midas = midas,
          experiment = experiment,
          test_covar = all.vars(formula(object))[1]
        )
      test_res <- dplyr::left_join(test_res, variables_freq, by = "term")
    }

    term_name <- switch (experiment,
                         "hla_alleles" = "allele",
                         "hla_aa" = "aa",
                         "expression_level" = "allele",
                         "hla_g_groups" = "g.group",
                         "hla_supertypes" = "supertype",
                         "hla_NK_ligands" = "allele.group",
                         "kir_genes" = "kir.gene",
                         "hla_kir_interactions" = "hla.kir.interaction",
                         "term"
    )
    test_res <- dplyr::rename(test_res, !!term_name := term)
    test_res <- dplyr::arrange(test_res, p.value)
    test_res <- test_res[, colnames(res)] # order columns


    expect_equal(as.data.frame(res), as.data.frame(test_res))
  }

  # frequency filtration
  conditional <- FALSE
  omnibus <- FALSE
  experiment_choice <-
    c(
      "hla_alleles",
      "hla_supertypes"
    )
  lower_frequency_cutoff <- 0.02
  upper_frequency_cutoff <- 0.06
  for (experiment in experiment_choice) {
    object <- lm(disease ~ term, data = midas)
    res <- runMiDAS(object,
                    conditional = conditional,
                    omnibus = omnibus,
                    experiment = experiment,
                    lower_frequency_cutoff = lower_frequency_cutoff,
                    upper_frequency_cutoff = upper_frequency_cutoff,
                    exponentiate = FALSE
    )

    midas_filtered <-
      filterByFrequency(
        object = midas,
        experiment = experiment,
        lower_frequency_cutoff = lower_frequency_cutoff,
        upper_frequency_cutoff = upper_frequency_cutoff
      )
    test_variables <- rownames(midas_filtered[[experiment]])
    midas_data <- midasToWide(midas_filtered, experiment = experiment)
    object$call$data <- midas_data
    test_res <- analyzeAssociations(
      object = object,
      variables = test_variables,
      exponentiate = FALSE
    )

    ex <- midas[[experiment]]
    if (is(ex, "SummarizedExperiment")) { ex <- assay(ex) }
    if (typeof(ex) == "integer") {
      variables_freq <-
        runMiDASGetVarsFreq(
          midas = midas_filtered,
          experiment = experiment,
          test_covar = all.vars(formula(object))[1]
        )
      test_res <-
        dplyr::left_join(test_res, variables_freq, by = "term")
    }

    term_name <- switch (experiment,
                         "hla_alleles" = "allele",
                         "hla_aa" = "aa",
                         "expression_level" = "allele",
                         "hla_g_groups" = "g.group",
                         "hla_supertypes" = "supertype",
                         "hla_NK_ligands" = "allele.group",
                         "kir_genes" = "kir.gene",
                         "hla_kir_interactions" = "hla.kir.interaction",
                         "term"
    )
    test_res <- dplyr::rename(test_res, !!term_name := term)
    test_res <- dplyr::arrange(test_res, p.value)
    test_res <- test_res[, colnames(res)] # order columns

    expect_equal(lapply(res, as.data.frame), lapply(test_res, as.data.frame))
  }

  #
  object <- lm(disease ~ term, data = midas)

  expect_error(runMiDAS(list()),
               "object was not recognized as a fit from a model function \\(such as lm, glm and many others\\)."
  )

  fake_object <- list(call = list(formula = 1 ~ 1, data = 1:5))
  class(fake_object) <- "foo"
  expect_error(runMiDAS(fake_object),
               "Could not find 'tidy' function for statistical model 'foo'. Please ensure that 'tidy' for selected model is available. See the 'broom' package for more information on 'tidy' function."
  )

  fake_object <- object
  fake_midas <- midas
  class(fake_midas) <- structure("MultiAssayExperiment", package = "MultiAssayExperiment")
  fake_object$call$data <- fake_midas
  expect_error(runMiDAS(fake_object),
               "data associated with statistical model must be an instance of MiDAS class."
  )

  fake_object <- object
  fake_object$call$formula <- disease ~ 1
  expect_error(
    runMiDAS(fake_object),
    "placeholder 'term' could not be found in object's formula"
  )

  expect_error(runMiDAS(object, experiment = 1),
               "experiment is not a string \\(a length one character vector\\)."
  )

  expect_error(runMiDAS(object, experiment = "foo"),
               "experiment should be one of \"hla_alleles\", \"hla_aa\", \"hla_g_groups\", \"hla_supertypes\", \"hla_NK_ligands\", \"kir_genes\", \"hla_kir_interactions\", \"hla_divergence\"."
  )

  expect_error(runMiDAS(object, experiment = "hla_alleles", inheritance_model = 1),
               "inheritance_model is not a string \\(a length one character vector\\) or NULL."
  )

  expect_error(runMiDAS(object, experiment = "hla_alleles", inheritance_model = "foo"),
               "inheritance_model should match values \"dominant\", \"recessive\", \"additive\"."
  )

  expect_error(runMiDAS(object, experiment = "hla_alleles", conditional = 1),
               "conditional is not a flag \\(a length one logical vector\\)."
  )

  expect_error(
    runMiDAS(
      object,
      experiment = "hla_alleles",
      conditional = TRUE,
      omnibus = 1
    ),
    "omnibus is not a flag \\(a length one logical vector\\)."
  )

  expect_error(
    runMiDAS(
      object,
      experiment = "hla_alleles",
      omnibus_groups_filter = 1
    ),
    "omnibus_groups_filter is not a character."
  )

  expect_error(runMiDAS(object, experiment = "hla_alleles", correction = 1),
               "correction is not a string \\(a length one character vector\\)."
  )

  expect_error(
    runMiDAS(
      object,
      experiment = "hla_alleles",
      conditional = TRUE,
      omnibus = FALSE,
      n_correction = "foo"
    ),
    "n_correction is not a count \\(a single positive integer\\) or NULL."
  )

  expect_error(
    runMiDAS(
      object,
      experiment = "hla_alleles",
      conditional = TRUE,
      omnibus = FALSE,
      exponentiate = "foo"
    ),
    "exponentiate is not a flag \\(a length one logical vector\\)."
  )

  expect_error(
    runMiDAS(
      object,
      experiment = "hla_divergence",
      conditional = FALSE,
      lower_frequency_cutoff = 0.1,
      upper_frequency_cutoff = 0.8
    ),
    "Frequency filtration does not support experiment 'hla_divergence'"
  )

  expect_error(
    runMiDAS(
      object,
      experiment = "hla_alleles",
      lower_frequency_cutoff = 0.53,
      upper_frequency_cutoff = 0.56
    ),
    "No variables available for analysis, please revisit your filtration criteria."
  )
})
