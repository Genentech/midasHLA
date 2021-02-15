context("stats_longtests")

test_that("runMiDAS", {
  midas <- MiDAS_tut_object
  midas_data <- midasToWide(
    MiDAS_tut_object,
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

    object$call$data <- midas_data
    test_variables <- rownames(midas[[experiment]])
    test_res <-
      analyzeAssociations(object, variables = test_variables, exponentiate = FALSE)

    ex <- midas[[experiment]]
    if (is(ex, "SummarizedExperiment")) {ex <- assay(ex)}
    if (typeof(ex) == "integer") {
      variables_freq <-
        midasHLA:::runMiDASGetVarsFreq(
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
        midasHLA:::runMiDASGetVarsFreq(
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
    residues = c("R", "V"),
    df = c(1, 1),
    statistic = c(1.00150233709155, 0.43786032955245),
    p.value = c(0.316947259129448, 0.508156995160654),
    p.adjusted = c(0.633894518258896, 1)
  )

  expect_equal(as.data.frame(res), as.data.frame(test_res))

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
    residues = c("L", "Q"),
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
      residues = c("L", "K"),
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
               "experiment should be one of \"hla_alleles\", \"hla_aa\", \"hla_g_groups\", \"hla_supertypes\", \"hla_NK_ligands\", \"kir_genes\", \"kir_haplotypes\", \"hla_kir_interactions\", \"hla_divergence\", \"hla_het\"."
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
