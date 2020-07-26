context("MiDAS class")

test_that("MiDAS object is valid", {
  pheno <- MiDAS_tut_pheno
  rownames(pheno) <- pheno$ID
  pheno$term <- runif(nrow(pheno))

  midas <-
    MultiAssayExperiment(
      experiments = ExperimentList(
        list(
          hla_calls = dfToExperimentMat(MiDAS_tut_HLA),
          kir_calls = dfToExperimentMat(MiDAS_tut_KIR)
        )
      ),
      colData = pheno,
      metadata = list(placeholder = "term")
    )
  class(midas) <- structure("MiDAS", package = "MiDAS")

  bad_placeholder <- midas
  S4Vectors::metadata(bad_placeholder)$placeholder <- 1
  expect_error(
    validObject(bad_placeholder),
    "placeholder is not a string \\(a length one character vector\\)."
  )

  bad_placeholder <- midas
  S4Vectors::metadata(bad_placeholder)$placeholder <- "A_1"
  expect_error(
    validObject(bad_placeholder),
    "Placeholder 'A_1' is used in one of the experiments"
  )

  bad_placeholder <- midas
  S4Vectors::metadata(bad_placeholder)$placeholder <- "foo"
  expect_error(
    validObject(bad_placeholder),
    "Placeholder 'foo' can not be found in object's colData"
  )
})

test_that("MiDAS object's experiment is extracted correctly", {
  midas <- prepareMiDAS(
    hla_calls = MiDAS_tut_HLA,
    colData = MiDAS_tut_pheno,
    experiment = "hla_alleles"
  )

  expect_equal(getExperiments(midas), "hla_alleles")
})

test_that("MiDAS object's hla_calls is extracted correctly", {
  midas <- prepareMiDAS(
    hla_calls = MiDAS_tut_HLA,
    colData = MiDAS_tut_pheno,
    experiment = "hla_alleles"
  )

  expect_equal(getHlaCalls(midas), MiDAS_tut_HLA)
})

test_that("MiDAS object's kir_calls is extracted correctly", {
  midas <- prepareMiDAS(
    kir_calls = MiDAS_tut_KIR,
    colData = MiDAS_tut_pheno,
    experiment = "kir_genes"
  )

  expect_equal(getKirCalls(midas), MiDAS_tut_KIR)
})

test_that("MiDAS object's placeholder is extracted correctly", {
  midas <- prepareMiDAS(
    kir_calls = MiDAS_tut_KIR,
    colData = MiDAS_tut_pheno,
    experiment = "kir_genes"
  )

  expect_equal(getPlaceholder(midas), "term")
})

test_that("MiDAS object's omnibus groups are extracted correctly", {
  midas <- prepareMiDAS(
    hla_calls = MiDAS_tut_HLA[, 1:3],
    colData = MiDAS_tut_pheno,
    experiment = "hla_aa"
  )

  expect_equal(getOmnibusGroups(midas, "hla_aa")[1:3],
               list(
                 `A_-22` = c("A_-22_*", "A_-22_V", "A_-22_I"),
                 `A_-21` = c("A_-21_*", "A_-21_M", "A_-21_V"),
                 `A_-18` = c("A_-18_*", "A_-18_R", "A_-18_P")
               ))

  expect_error(
    getOmnibusGroups(midas, 1),
    "experiment is not a string \\(a length one character vector\\)."
  )
})

test_that("MiDAS object's frequencies are extracted correctly", {
  hla_calls <- reduceHlaCalls(MiDAS_tut_HLA, resolution = 4) # allels have to be in same res. for 'hla_divergence'
  midas <- prepareMiDAS(
    hla_calls = hla_calls,
    colData = MiDAS_tut_pheno,
    experiment = c("hla_alleles", "hla_divergence")
  )

  alleles_subset <-
    c("A*01:01", "A*02:01", "A*02:06", "A*26:01", "B*07:02", "B*08:01", "B*13:02", "B*15:01", "B*27:05", "B*40:01", "B*57:01")
  midas <- midas[alleles_subset, ]
  freq <- getFrequencies(midas, "hla_alleles")
  test_freq <- data.frame(
    term = alleles_subset,
    Counts = c(236, 486, 22, 90, 179, 151, 66, 86, 59, 58, 44),
    Freq = formattable::percent(
      c(0.118, 0.243, 0.011, 0.045, 0.0895, 0.0755, 0.033, 0.043, 0.0295, 0.029, 0.022)
    ),
    stringsAsFactors = FALSE
  )
  rownames(test_freq) <- alleles_subset

  expect_equal(freq, test_freq)

  expect_error(
    getFrequencies(midas, 1),
    "experiment is not a string \\(a length one character vector\\)."
  )

  expect_error(getFrequencies(midas, "foo"),
               "experiment should be one of \"hla_alleles\", \"hla_divergence\"."
  )

  expect_error(getFrequencies(midas, "hla_divergence"),
               "Frequencies can not be calculated for experiment 'hla_divergence'"
  )
})

test_that("MiDAS's as.data.frame method works properly", {
  midas <- prepareMiDAS(
    kir_calls = MiDAS_tut_KIR,
    colData = MiDAS_tut_pheno,
    experiment = "kir_genes"
  )

  midas_df <- as.data.frame(midas)
  test_midas_df <-
    midasToWide(midas,
                experiment = getExperiments(midas))

  expect_equal(midas_df, test_midas_df)
})

test_that("MiDAS object is prepared properly", {
  hla_calls <- reduceHlaCalls(MiDAS_tut_HLA, resolution = 4)
  kir_calls <- MiDAS_tut_KIR
  phenotype <- MiDAS_tut_pheno

  args_c <-
    expand.grid(
      lower_frequency_cutoff = c(0, 0.54),
      upper_frequency_cutoff = c(1, 0.63),
      stringsAsFactors = FALSE
    )

  for (i in 1:nrow(args_c)) {
    args <- as.list(args_c[i, , drop = FALSE])

    midas <- prepareMiDAS(
      hla_calls = hla_calls,
      kir_calls = kir_calls,
      colData = phenotype,
      experiment = c(
        "hla_alleles",
        # "hla_aa",
        "hla_g_groups",
        "hla_supertypes",
        "hla_NK_ligands",
        "kir_genes",
        "hla_kir_interactions",
        "hla_divergence"
      ),
      lower_frequency_cutoff = args$lower_frequency_cutoff,
      upper_frequency_cutoff = args$upper_frequency_cutoff
    )
    validObject(midas)
  }

  expect_error(
    prepareMiDAS(),
    "hla_calls or kir_calls argument has to be specified"
  )


  expect_error(
    prepareMiDAS(hla_calls = cars[1:2, ]),
    "values: 2, 10 in hla_calls doesn't follow HLA numbers specification"
  )

  # TODO
  # expect_error(
  #   prepareMiDAS(kir_calls = cars[1:2, ]),
  #   "values: 2, 10 in kir_calls doesn't follow HLA numbers specification"
  # )

  expect_error(
    prepareMiDAS(hla_calls = hla_calls, colData = cars[1:2, ]),
    "first column in colData must be named 'ID'"
  )

  expect_error(
    prepareMiDAS(
      hla_calls = hla_calls,
      colData = phenotype,
      experiment = 1
    ),
    "experiment is not a character vector"
  )

  expect_error(
    prepareMiDAS(
      hla_calls = hla_calls,
      colData = phenotype,
      experiment = "foo"
    ),
    "experiment should match values \"hla_alleles\", \"hla_aa\", \"hla_g_groups\", \"hla_supertypes\", \"hla_NK_ligands\", \"kir_genes\", \"hla_kir_interactions\"."
  )

  expect_error(
    prepareMiDAS(
      hla_calls = hla_calls,
      colData = phenotype,
      experiment = "hla_alleles",
      placeholder = 1
    ),
    "placeholder is not a string \\(a length one character vector\\)."
  )

  expect_error(
    prepareMiDAS(
      hla_calls = hla_calls,
      colData = phenotype,
      experiment = "hla_alleles",
      placeholder = "outcome"
    ),
    "Placeholder 'outcome' can not be used, it is alredy used as column name in one of the inputs."
  )
})

test_that("MiDAS data for hla_alleles analysis is prepared properly", {
  args <- list(hla_calls = MiDAS_tut_HLA)
  experiment <- do.call(prepareMiDAS_hla_alleles, args)
  experiment_test <- do.call(hlaCallsToCounts, args)
  experiment_test <- dfToExperimentMat(experiment_test)
  expect_equal(experiment, experiment_test)
})

test_that("MiDAS data for hla_aa analysis is prepared properly", { # TODO
  args_c <- expand.grid(
    indels = c(TRUE, FALSE),
    unkchar = c(TRUE, FALSE),
    stringsAsFactors = FALSE
  )

  for (i in 1:nrow(args_c)) {
    args <- as.list(args_c[i, , drop = FALSE])
    args$hla_calls <- MiDAS_tut_HLA
    experiment <- do.call(prepareMiDAS_hla_aa, args)

    counts <-
      hlaToAAVariation(
        hla_calls = args$hla_calls,
        indels = args$indels,
        unkchar = args$unkchar
      )
    counts <-
      aaVariationToCounts(counts)
    counts <- dfToExperimentMat(counts)
    pos <- gsub("_.{1}$", "", rownames(counts))
    pos <- unique(pos)
    omnibus_groups <- lapply(pos, function(p) {
      grep(pattern = paste0("^", p, "_"),
           x = rownames(counts),
           value = TRUE)
    })
    names(omnibus_groups) <- pos
    experiment_test <-
      SummarizedExperiment(assays = counts,
                           metadata = list(omnibus_groups = omnibus_groups))

    expect_equal(experiment, experiment_test)
  }
})

test_that("MiDAS data for hla_g_groups analysis is prepared properly", {
  args <- list(hla_calls = MiDAS_tut_HLA)

  experiment <- do.call(prepareMiDAS_hla_g_groups, args)

  experiment_test <-
    hlaToVariable(
      hla_calls = args$hla_calls,
      dictionary = "allele_HLA_Ggroup",
      na.value = 0
    )
  experiment_test <-
    hlaCallsToCounts(experiment_test,
                     check_hla_format = FALSE)
  experiment_test <- dfToExperimentMat(experiment_test)

  expect_equal(experiment, experiment_test)

})

test_that("MiDAS data for hla_supertypes analysis is prepared properly", {
  args <- list(hla_calls = MiDAS_tut_HLA)

  experiment <- do.call(prepareMiDAS_hla_supertypes, args)

  experiment_test <-
    hlaToVariable(
      hla_calls = args$hla_calls,
      dictionary = "allele_HLA_supertype",
      na.value = 0
    )
  experiment_test <-
    hlaCallsToCounts(experiment_test,
                     check_hla_format = FALSE)
  experiment_test <- dfToExperimentMat(experiment_test)
  experiment_test <-
    experiment_test[rownames(experiment_test) != "Unclassified",]

  expect_equal(experiment, experiment_test)

})

test_that("MiDAS data for hla_NK_ligands analysis is prepared properly", {
  args <- list(hla_calls = MiDAS_tut_HLA)

  experiment <- do.call(prepareMiDAS_hla_NK_ligands, args)

  experiment_test <-
    Reduce(
      f = function(...)
        left_join(..., by = "ID"),
      x = lapply(
        c(
          "allele_HLA-B_Bw",
          "allele_HLA_Bw4+A23+A24+A32",
          "allele_HLA-C_C1-2"
        ),
        hlaToVariable,
        hla_calls = args$hla_calls,
        na.value = 0
      )
    )
  experiment_test <-
    hlaCallsToCounts(experiment_test,
                     check_hla_format = FALSE)
  experiment_test <- dfToExperimentMat(experiment_test)

  expect_equal(experiment, experiment_test)

})

test_that("MiDAS data for kir_genes analysis is prepared properly", {
  args <- list(kir_calls = MiDAS_tut_KIR)
  experiment <- do.call(prepareMiDAS_kir_genes, args)

  experiment_test <- dfToExperimentMat(args$kir_calls)

  expect_equal(experiment, experiment_test)
})

test_that("MiDAS data for hla_kir_interactions analysis is prepared properly", {
  args <- list(hla_calls = MiDAS_tut_HLA, kir_calls = MiDAS_tut_KIR)
  experiment <- do.call(prepareMiDAS_kir_genes, args)

  experiment_test <- getHlaKirInteractions(
    hla_calls = args$hla_calls,
    kir_counts = args$kir_calls
  )
  experiment_test <- dfToExperimentMat(args$kir_calls)

  expect_equal(experiment, experiment_test)
})

test_that("MiDAS data for hla_divergence analysis is prepared properly", {
  hla_calls <- reduceHlaCalls(MiDAS_tut_HLA, 4)
  experiment <- do.call(prepareMiDAS_hla_divergence, list(hla_calls))

  experiment_test <- hlaCallsGranthamDistance(hla_calls = hla_calls,
                                              genes = c("A", "B", "C"))
  experiment_test$ABC_avg <- rowMeans(experiment_test[-1])
  experiment_test <- dfToExperimentMat(experiment_test)

  expect_equal(experiment, experiment_test)


  hla_calls <- hla_calls[, c(1, 8:9)]
  expect_error(
    prepareMiDAS_hla_divergence(hla_calls),
    "Grantham distance can be calculated only for class I HLA alleles \\(A, B, C\\)."
  )
})

test_that("MiDAS is filtered correctly", {
  hla_calls <- reduceHlaCalls(MiDAS_tut_HLA, 4)
  midas <- prepareMiDAS(
    hla_calls = hla_calls,
    kir_calls = MiDAS_tut_KIR,
    colData = MiDAS_tut_pheno,
    experiment = c("hla_alleles", "hla_supertypes", "kir_genes", "hla_divergence")
  )

  # filtration works as expected
  experiment <- "hla_alleles"
  lower_frequency_cutoff <- 0.53
  upper_frequency_cutoff <- 0.56
  filtered_midas <-
    filterByFrequency(
      object = midas,
      experiment = experiment,
      lower_frequency_cutoff = lower_frequency_cutoff,
      upper_frequency_cutoff = upper_frequency_cutoff
    )
  test_filtered_midas <- midas
  test_filtered_midas[[experiment]] <-
    filterExperimentByFrequency(
      experiment = midas[[experiment]],
      lower_frequency_cutoff = lower_frequency_cutoff,
      upper_frequency_cutoff = upper_frequency_cutoff
    )
  expect_equal(filtered_midas, test_filtered_midas)

  # unsported experiments raise error
  experiment <- "hla_divergence"
  lower_frequency_cutoff <- 0.53
  upper_frequency_cutoff <- 0.56
  expect_error(
    filterByFrequency(
      object = midas,
      experiment = experiment,
      lower_frequency_cutoff = lower_frequency_cutoff,
      upper_frequency_cutoff = upper_frequency_cutoff
    ),
    "Frequency filtration does not support experiment 'hla_divergence'"
  )
})

test_that("filterByOmnibusGroups", {
  midas <- prepareMiDAS(
    hla_calls = MiDAS_tut_HLA,
    colData = MiDAS_tut_pheno,
    experiment = c("hla_alleles", "hla_aa")
  )

  # filtration works as expected
  experiment <- "hla_aa"
  mask <- c("A_83", "A_90")
  filtered_midas <-
    filterByOmnibusGroups(
      object = midas,
      experiment = experiment,
      groups = mask
    )
  test_filtered_midas <- midas
  vars <- c("A_83_G", "A_83_R",  "A_90_D", "A_90_A")
  test_filtered_midas[[experiment]] <- test_filtered_midas[[experiment]][vars, ]
  metadata(test_filtered_midas[[experiment]])$omnibus_groups <-
    metadata(test_filtered_midas[[experiment]])$omnibus_groups[mask]
  expect_equal(filtered_midas, test_filtered_midas)

  # unsported experiments raise error
  experiment <- "hla_alleles"
  expect_error(
    filterByOmnibusGroups(
      object = midas,
      experiment = experiment,
      groups = mask
    ),
    "Omnibus groups filtration does not support experiment 'hla_alleles'"
  )
})
