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
      correction = "BH",
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
                                        variables = c("DQB1*06:02", "B*57:01"),
                                        th = 0.05,
                                        keep = FALSE)
  res <- rapply(res, classes = "numeric", how = "replace", round, digits = 3)

  test_res <- tibble(term = c("DQB1*06:02", "B*57:01"),
                     estimate = c(0.15, 0.235),
                     std.error = c(0.035, 0.056),
                     statistic = c(4.319, 4.228),
                     p.value = c(0, 0),
                     conf.low = c(0.082, 0.126),
                     conf.high = c(0.219, 0.344),
                     p.adjusted = c(0, 0),
                     covariates = c("", "DQB1*06:02")
  )

  expect_equal(res, test_res)

  # keep = TRUE
  res <- analyzeConditionalAssociations(object,
                                        variables = c("DQB1*06:02", "B*57:01"),
                                        th = 0.05,
                                        keep = TRUE)
  res <- rapply(res, classes = "numeric", how = "replace", round, digits = 3)
  test_res <- list(
    tibble(term = c("DQB1*06:02", "B*57:01"),
           estimate = c(0.15, 0.236),
           std.error = c(0.035, 0.056),
           statistic = c(4.319, 4.21),
           p.value =  c(0, 0),
           conf.low = c(0.082, 0.126),
           conf.high = c(0.219, 0.346),
           p.adjusted = c(0, 0),
           covariates = c("", "")
    ),
    tibble(term = "B*57:01",
           estimate = 0.235,
           std.error = 0.056,
           statistic = 4.228,
           p.value = 0,
           conf.low = 0.126,
           conf.high = 0.344,
           p.adjusted = 0,
           covariates = "DQB1*06:02"
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
      correction = "BH",
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
    df = sapply(LRT, `[[`, "df"),
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
                    inheritance_model = "additive",
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
                    inheritance_model = "additive",
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
                    inheritance_model = "additive",
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

    expect_equal(as.data.frame(res), as.data.frame(test_res))
  }

  # linear omnibus
  conditional <- FALSE
  omnibus <- TRUE
  experiment <- "hla_aa"
  object <- lm(disease ~ term, data = midas)
  res <- runMiDAS(
    object,
    inheritance_model = "additive",
    conditional = conditional,
    omnibus = omnibus,
    omnibus_groups_filter = c("A_6", "A_12"),
    experiment = experiment,
    exponentiate = FALSE
  )
  test_res <- dplyr::tibble(
    aa_pos = c("A_6", "A_12"),
    residues = c("R, G", "V, M"),
    df = c(1, 1),
    statistic = c(1.00150233709155, 0.43786032955245),
    p.value = c(0.316947259129448, 0.508156995160654),
    p.adjusted = c(0.633894518258896, 1)
  )

  expect_equal(as.data.frame(res), as.data.frame(test_res))

  expect_error(
    runMiDAS(
      object,
      inheritance_model = "additive",
      conditional = FALSE,
      omnibus = TRUE,
      experiment = "hla_alleles",
      exponentiate = FALSE
    ),
    "Omnibus test does not support experiment hla_alleles"
  )

  # conditional omnibus
  conditional <- TRUE
  omnibus <- TRUE
  experiment <- "hla_aa"
  object <- lm(disease ~ term, data = midas)
  res <- runMiDAS(
    object,
    inheritance_model = "additive",
    conditional = conditional,
    omnibus = omnibus,
    omnibus_groups_filter = c("DRA_217", "DQA1_34"),
    experiment = experiment,
    exponentiate = FALSE
  )
  test_res <- dplyr::tibble(
    aa_pos = c("DRA_217", "DQA1_34"),
    residues = c("L, V", "Q, E"),
    df = c(1, 1),
    statistic = c(11.7515320865891, 8.92028215923278),
    p.value = c(0.000607931755127601, 0.00282020937328977),
    p.adjusted = c(0.00121586351025461, 0.00282020937329118),
    covariates = c("", "DRA_217")
  )
  expect_equal(as.data.frame(res), as.data.frame(test_res))

  res <- runMiDAS(
    object,
    inheritance_model = "additive",
    conditional = conditional,
    omnibus = omnibus,
    omnibus_groups_filter = c("B_178", "DRA_217"),
    experiment = experiment,
    exponentiate = FALSE,
    keep = TRUE
  )
  test_res <- list(
    dplyr::tibble(
      aa_pos = c("DRA_217", "B_178"),
      residues = c("L, V", "K, T"),
      df = c(1, 1),
      statistic = c(11.7515320865891, 5.63638399577462),
      p.value = c(0.000607931755127601, 0.0175914523736528),
      p.adjusted = c(0.0012158635102552, 0.0351829047473056),
      covariates = c("", "")
    )
  )
  expect_equal(lapply(res, as.data.frame), lapply(test_res, as.data.frame))

  expect_error(
    runMiDAS(
      object,
      inheritance_model = "additive",
      conditional = FALSE,
      omnibus = TRUE,
      experiment = "hla_alleles",
      exponentiate = FALSE
    ),
    "Omnibus test does not support experiment hla_alleles"
  )

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
               "inheritance_model should match values \"dominant\", \"recessive\", \"additive\", \"overdominant\"."
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
      inheritance_model = "additive",
      experiment = "hla_alleles",
      lower_frequency_cutoff = 0.53,
      upper_frequency_cutoff = 0.56
    ),
    "No variables available for analysis, please revisit your filtration criteria."
  )
})

test_that("HWETest", {
  midas <- prepareMiDAS(hla_calls = MiDAS_tut_HLA,
                        colData = MiDAS_tut_pheno,
                        experiment = "hla_alleles")
  res <- HWETest(midas, experiment = "hla_alleles", HWE_cutoff = 0.988)
  rownames(res) <- NULL
  test_res <- data.frame(
    var = c("A*02:01", "DRB1*13:02"),
    p.value = c(0.993279078042026, 0.988280722372832),
    stringsAsFactors = FALSE
  )
  expect_equal(as.data.frame(res), as.data.frame(test_res))

  res <-
    HWETest(
      midas,
      experiment = "hla_alleles",
      HWE_cutoff = 0.988,
      as.MiDAS = TRUE
    )
  test_res <- filterByVariables(midas, "hla_alleles", c("A*02:01", "DRB1*13:02"))
  expect_equal(res, test_res)

  expect_error(HWETest("foo"), "object must be an instance of \"MiDAS\".")

  expect_error(HWETest(midas, experiment = 1), "experiment is not a string \\(a length one character vector\\).")

  expect_error(HWETest(midas, experiment = "foo"), "experiment should be one of \"hla_alleles\", \"hla_aa\", \"hla_g_groups\", \"hla_supertypes\", \"hla_NK_ligands\".")

  expect_error(HWETest(midas, experiment = "hla_alleles", HWE_cutoff = "foo"), "HWE_cutoff is not a number \\(a length one numeric vector\\) or NULL.")

  expect_error(HWETest(midas, experiment = "hla_alleles", as.MiDAS = "foo"), "as.MiDAS is not a flag \\(a length one logical vector\\).")
})
