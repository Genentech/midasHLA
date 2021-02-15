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
  class(midas) <- structure("MiDAS", package = "midasHLA")

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
    "Placeholder 'A_1' is used in one of object's experiments"
  )

  bad_placeholder <- midas
  S4Vectors::metadata(bad_placeholder)$placeholder <- "foo"
  expect_error(
    validObject(bad_placeholder),
    "Placeholder 'foo' can not be found in object's colData"
  )

  duplicated_features <- midas
  rownames(duplicated_features[["hla_calls"]])[1] <- "ID"
  expect_error(
    validObject(duplicated_features),
    "Object contain duplicated features: ID"
  )
})

test_that("getExperiments", {
  expect_equal(
    getExperiments(MiDAS_tut_object),
    c(
      "hla_alleles",
      "hla_aa",
      "hla_g_groups",
      "hla_supertypes",
      "hla_NK_ligands",
      "kir_genes",
      "kir_haplotypes",
      "hla_kir_interactions",
      "hla_divergence",
      "hla_het"
    )
  )
})

test_that("getHlaCalls", {
  expect_equal(getHlaCalls(MiDAS_tut_object), MiDAS_tut_HLA)
})

test_that("getKirCalls", {
  expect_equal(getKirCalls(MiDAS_tut_object), MiDAS_tut_KIR)
})

test_that("getPlaceholder", {
  expect_equal(getPlaceholder(MiDAS_tut_object), "term")
})

test_that("getOmnibusGroups", {
  midas <-
    filterByOmnibusGroups(MiDAS_tut_object, "hla_aa", c("A_-22", "A_-21", "A_-18"))
  
  expect_equal(getOmnibusGroups(midas, "hla_aa"),
               list(
                 `A_-22` = c("A_-22_V", "A_-22_*", "A_-22_I"),
                 `A_-21` = c("A_-21_M", "A_-21_V", "A_-21_*"),
                 `A_-18` = c("A_-18_R", "A_-18_*", "A_-18_P")
               ))

  expect_error(
    getOmnibusGroups(midas, 1),
    "experiment is not a string \\(a length one character vector\\)."
  )
})

test_that("getFrequencies", {
  alleles_subset <-
    c("A*01:01", "A*02:01", "A*02:06", "A*26:01", "B*07:02", "B*08:01", "B*13:02", "B*15:01", "B*27:05", "B*40:01", "B*57:01")
  midas_sub <- MiDAS_tut_object[alleles_subset, ]
  freq <- getFrequencies(midas_sub, "hla_alleles")
  test_freq <- data.frame(
    term = alleles_subset,
    Counts = c(236, 486, 22, 90, 152, 128, 74, 86, 59, 58, 81),
    Freq = c(0.118, 0.243, 0.011, 0.045, 0.076, 0.064, 0.037, 0.043, 0.0295, 0.029, 0.0405),
    stringsAsFactors = FALSE
  )

  expect_equal(freq, test_freq)

  expect_error(
    getFrequencies(MiDAS_tut_object, 1),
    "experiment is not a string \\(a length one character vector\\)."
  )

  expect_error(getFrequencies(MiDAS_tut_object, "foo"),
               "experiment should be one of \"hla_alleles\", \"hla_aa\", \"hla_g_groups\", \"hla_supertypes\", \"hla_NK_ligands\", \"kir_genes\", \"kir_haplotypes\", \"hla_kir_interactions\", \"hla_divergence\", \"hla_het\"."
  )

  expect_error(getFrequencies(MiDAS_tut_object, "hla_divergence"),
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

test_that("filterByFrequency", {
  # filtration works as expected
  experiment <- "hla_alleles"
  lower_frequency_cutoff <- 0.53
  upper_frequency_cutoff <- 0.56
  filtered_midas <-
    filterByFrequency(
      object = MiDAS_tut_object,
      experiment = experiment,
      lower_frequency_cutoff = lower_frequency_cutoff,
      upper_frequency_cutoff = upper_frequency_cutoff
    )
  test_filtered_midas <- MiDAS_tut_object
  test_filtered_midas[[experiment]] <-
    filterExperimentByFrequency(
      experiment = test_filtered_midas[[experiment]],
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
      object = MiDAS_tut_object,
      experiment = experiment,
      lower_frequency_cutoff = lower_frequency_cutoff,
      upper_frequency_cutoff = upper_frequency_cutoff
    ),
    "Frequency filtration does not support experiment 'hla_divergence'"
  )
})

test_that("filterByOmnibusGroups", {
  # filtration works as expected
  experiment <- "hla_aa"
  mask <- c("A_83", "A_90")
  filtered_midas <-
    filterByOmnibusGroups(
      object = MiDAS_tut_object,
      experiment = experiment,
      groups = mask
    )
  test_filtered_midas <- MiDAS_tut_object
  vars <- c("A_83_R", "A_83_G", "A_90_D", "A_90_A")
  test_filtered_midas[[experiment]] <- test_filtered_midas[[experiment]][vars, ]

  expect_equal(filtered_midas, test_filtered_midas)

  # unsported experiments raise error
  experiment <- "hla_alleles"
  expect_error(
    filterByOmnibusGroups(
      object = MiDAS_tut_object,
      experiment = experiment,
      groups = mask
    ),
    "Omnibus groups filtration does not support experiment 'hla_alleles'"
  )
})

test_that("getAllelesForAA", {
  aa <- getAllelesForAA(MiDAS_tut_object, "DRA_2")
  aa_test <- data.frame(
    `HLA-DRA (2)` = c("*", "A"),
    `HLA-DRA alleles` = c("*01:02", "*01:01"),
    count = c(711L, 1289L),
    frequency = formattable::percent(c(0.3555, 0.6445)),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  expect_error(getAllelesForAA(MiDAS_tut_object, 2),
               "aa_pos is not a string \\(a length one character vector\\).")

  expect_error(getAllelesForAA(MiDAS_tut_object, "A"),
               "amino acid position should be formatted like: A_9.")

  MultiAssayExperiment::metadata(MiDAS_tut_object)[["hla_calls"]] <- NULL
  expect_error(getAllelesForAA(MiDAS_tut_object, "A_9"),
               "Could not find HLA calls associated with MiDAS object. Make sure to use prepareMiDAS for MiDAS object creation.")
})

test_that("prepareMiDAS_hla_alleles", {
  args <- list(hla_calls = MiDAS_tut_HLA)
  experiment <- do.call(prepareMiDAS_hla_alleles, args)
  experiment_test <- do.call(hlaCallsToCounts, args)
  experiment_test <- dfToExperimentMat(experiment_test)
  experiment_test <- SummarizedExperiment::SummarizedExperiment(
    assays = list(experiment_test),
    metadata = list(
      inheritance_model_applicable = TRUE,
      pop_mul = 2,
      omnibus_groups = NULL
    )
  )
  expect_equal(experiment, experiment_test)
})

test_that("prepareMiDAS_hla_g_groups", {
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
  experiment_test <- SummarizedExperiment::SummarizedExperiment(
    assays = list(experiment_test),
    metadata = list(
      inheritance_model_applicable = TRUE,
      pop_mul = 2,
      omnibus_groups = NULL
    )
  )

  expect_equal(experiment, experiment_test)

})

test_that("prepareMiDAS_hla_supertypes", {
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
  experiment_test <- SummarizedExperiment::SummarizedExperiment(
    assays = list(experiment_test),
    metadata = list(
      inheritance_model_applicable = TRUE,
      pop_mul = 2,
      omnibus_groups = NULL
    )
  )

  expect_equal(experiment, experiment_test)

})

test_that("prepareMiDAS_hla_NK_ligands", {
  args <- list(hla_calls = MiDAS_tut_HLA)

  experiment <- do.call(prepareMiDAS_hla_NK_ligands, args)

  allele_HLA_Bw <- read.table(
    file = system.file("extdata", "Match_allele_HLA_Bw.txt", package = "midasHLA"),
    header = TRUE,
    sep = "\t",
    stringsAsFactors = FALSE
  )
  allele_HLA_Bw <- allele_HLA_Bw[allele_HLA_Bw$group != "Bw6", ]
  experiment_test <-
    Reduce(
      f = function(...)
        left_join(..., by = "ID"),
      x = lapply(
        list(
          allele_HLA_Bw,
          "allele_HLA-Bw_only_B",
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
  experiment_test <- SummarizedExperiment::SummarizedExperiment(
    assays = list(experiment_test),
    metadata = list(
      inheritance_model_applicable = TRUE,
      pop_mul = 2,
      omnibus_groups = NULL
    )
  )

  expect_equal(experiment, experiment_test)
})

test_that("prepareMiDAS_kir_genes", {
  args <- list(kir_calls = MiDAS_tut_KIR)
  experiment <- do.call(prepareMiDAS_kir_genes, args)

  experiment_test <- dfToExperimentMat(args$kir_calls)
  experiment_test <- SummarizedExperiment::SummarizedExperiment(
    assays = list(experiment_test),
    metadata = list(
      inheritance_model_applicable = FALSE,
      pop_mul = 1,
      omnibus_groups = NULL
    )
  )

  expect_equal(experiment, experiment_test)
})

test_that("prepareMiDAS_kir_haplotypes", {
  args <- list(kir_calls = MiDAS_tut_KIR)
  experiment <- do.call(prepareMiDAS_kir_haplotypes, args)

  experiment_test <- countsToVariables(args$kir_calls, "kir_haplotypes")
  experiment_test <- dfToExperimentMat(experiment_test)
  experiment_test <- SummarizedExperiment::SummarizedExperiment(
    assays = list(experiment_test),
    metadata = list(
      inheritance_model_applicable = FALSE,
      pop_mul = 1,
      omnibus_groups = NULL
    )
  )

  expect_equal(experiment, experiment_test)
})

test_that("prepareMiDAS_hla_kir_interactions", {
  args <- list(hla_calls = MiDAS_tut_HLA, kir_calls = MiDAS_tut_KIR)
  experiment <- do.call(prepareMiDAS_kir_genes, args)

  experiment_test <- getHlaKirInteractions(
    hla_calls = args$hla_calls,
    kir_calls = args$kir_calls
  )
  experiment_test <- dfToExperimentMat(args$kir_calls)
  experiment_test <- SummarizedExperiment::SummarizedExperiment(
    assays = list(experiment_test),
    metadata = list(
      inheritance_model_applicable = FALSE,
      pop_mul = 1,
      omnibus_groups = NULL
    )
  )

  expect_equal(experiment, experiment_test)
})

test_that("prepareMiDAS_hla_divergence", {
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

test_that("prepareMiDAS_hla_het", {
  hla_calls <- data.frame(
    ID = c("SAM1", "SAM2", "SAM3", "SAM4"),
    A_1 = c("A*01:01:01", "A*01:01:02", "A*01:01:03", "A*01:01:04"),
    A_2 = c("A*01:01:01", "A*01:01:03", "A*01:01:02", "A*01:01:04"),
    stringsAsFactors = FALSE
  )

  experiment <- do.call(prepareMiDAS_hla_het, list(hla_calls))
  experiment_test <- SummarizedExperiment::SummarizedExperiment(
    assay = matrix(c(0L, 1L, 1L, 0L),
                   ncol = 4,
                   dimnames = list("A_het", c("SAM1", "SAM2", "SAM3", "SAM4"))
    ),
    metadata = list(
      inheritance_model_applicable = FALSE,
      pop_mul = 1,
      omnibus_groups = NULL
    )
  ) 
  
  expect_equal(experiment, experiment_test)

  experiment <- do.call(prepareMiDAS_hla_het, list(hla_calls, hla_het_resolution = 4))
  experiment_test <- SummarizedExperiment::SummarizedExperiment(
    assay = matrix(c(0L, 0L, 0L, 0L), 
                   ncol = 4, 
                   dimnames = list("A_het", c("SAM1", "SAM2", "SAM3", "SAM4"))
    ),
    metadata = list(
      inheritance_model_applicable = FALSE,
      pop_mul = 1,
      omnibus_groups = NULL
    )
  )

  expect_equal(experiment, experiment_test)

  expect_error(
    prepareMiDAS_hla_het(hla_calls, hla_het_resolution = "a"),
    "hla_het_resolution is not a count \\(a single positive integer\\)"
  )

  hla_calls <- data.frame(
    ID = c("SAM1", "SAM2", "SAM3", "SAM4"),
    F_1 = c("A*01:01:01", "A*01:01:02", "A*01:01:03", "A*01:01:04"),
    F_2 = c("A*01:01:01", "A*01:01:03", "A*01:01:02", "A*01:01:04"),
    stringsAsFactors = FALSE
  )
  expect_error(
    prepareMiDAS_hla_het(hla_calls),
    "Heterozygosity status can be calculated only for classical genes \\(A, B, C, DQA1, DQB1, DRA, DRB1, DPA1, DPB1\\)."
  )
})

test_that("prepareMiDAS_hla_custom", {
  hla_dictionary <- read.table(
    file = system.file("extdata", "Match_allele_HLA_Ggroup.txt", package = "midasHLA"),
    header = TRUE,
    sep = "\t",
    quote = "",
    stringsAsFactors = FALSE
  )
  midas <- prepareMiDAS(
    hla_calls = MiDAS_tut_HLA,
    colData = MiDAS_tut_pheno,
    experiment = c("hla_g_groups", "hla_custom"),
    hla_dictionary = hla_dictionary
  )
  expect_equal(midas[["hla_g_groups"]], midas[["hla_custom"]])
})

test_that("prepareMiDAS_kir_custom", {
  kir_dictionary <- read.table(
    file = system.file("extdata", "Match_counts_kir_haplotypes.txt", package = "midasHLA"),
    header = TRUE,
    sep = "\t",
    quote = "",
    stringsAsFactors = FALSE
  )
  midas <- prepareMiDAS(
    kir_calls = MiDAS_tut_KIR,
    colData = MiDAS_tut_pheno,
    experiment = c("kir_haplotypes", "kir_custom"),
    kir_dictionary = kir_dictionary
  )
  expect_equal(midas[["kir_haplotypes"]], midas[["kir_custom"]])
})
