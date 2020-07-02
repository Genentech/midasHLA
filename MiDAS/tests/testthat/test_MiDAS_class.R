context("MiDAS class")

test_that("MiDAS object is valid", {
  rleft_join <- function(init, ..., by = "ID") {
    df <- Reduce(function(...)
      dplyr::left_join(..., by = by),
      x = list(...),
      init = init)
    rownames(df) <- df[[by]]
    df
  }

  hla_calls_file <- system.file("extdata", "HLAHD_output_example.txt", package = "MiDAS")
  hla_calls <- readHlaCalls(hla_calls_file)

  kir_calls_file <- system.file("extdata", "KPI_output_example.txt", package = "MiDAS")
  kir_calls <- readKPICalls(kir_calls_file)

  pheno_file <- system.file("extdata", "pheno_example.txt", package = "MiDAS")
  pheno <- read.table(pheno_file, header = TRUE, stringsAsFactors = FALSE)
  covar_file <- system.file("extdata", "covar_example.txt", package = "MiDAS")
  covar <- read.table(covar_file, header = TRUE, stringsAsFactors = FALSE)
  phenotype <- rleft_join(pheno, covar)

  midas <-
    MultiAssayExperiment(
      experiments = ExperimentList(
        list(
          hla_calls = dfToExperimentMat(hla_calls),
          kir_calls = dfToExperimentMat(kir_calls)
        )
      ),
      colData = phenotype,
      metadata = list(inheritance_model = "dominant")
    )
  class(midas) <- structure("MiDAS", package = "MiDAS")

  empty_midas <- midas
  experiments(empty_midas) <- MultiAssayExperiment::ExperimentList()
  expect_error(
    validObject(empty_midas),
    "MiDAS object must contain hla_calls or kir_calls experiment"
  )

  # unit tests of checkHlaCalls and checkKirCalls functions are ommited here

  no_inheritance_midas <- midas
  metadata(no_inheritance_midas)$inheritance_model <- NULL
  expect_error(
    validObject(no_inheritance_midas),
    "inheritance_model is not a string \\(a length one character vector\\)."
  )

  bad_inheritance_midas <- midas
  metadata(bad_inheritance_midas)$inheritance_model <- "foo"
  expect_error(
    validObject(bad_inheritance_midas),
    "inheritance_model should be one of \"additive\", \"dominant\", \"recessive\"."
  )
})

test_that("MiDAS object's inheritance model is extracted correctly", {
  hla_calls_file <- system.file("extdata", "HLAHD_output_example.txt", package = "MiDAS")
  hla_calls <- readHlaCalls(hla_calls_file)

  pheno_file <- system.file("extdata", "pheno_example.txt", package = "MiDAS")
  pheno <- read.table(pheno_file, header = TRUE, stringsAsFactors = FALSE)

  midas <- prepareMiDAS(
    hla_calls = hla_calls,
    colData = pheno,
    inheritance_model = "additive",
    experiment = "hla_alleles"
  )

  expect_equal(getInheritanceModel(midas), "additive")
})

test_that("MiDAS object's experiment is extracted correctly", {
  hla_calls_file <- system.file("extdata", "HLAHD_output_example.txt", package = "MiDAS")
  hla_calls <- readHlaCalls(hla_calls_file)

  pheno_file <- system.file("extdata", "pheno_example.txt", package = "MiDAS")
  pheno <- read.table(pheno_file, header = TRUE, stringsAsFactors = FALSE)

  midas <- prepareMiDAS(
    hla_calls = hla_calls,
    colData = pheno,
    inheritance_model = "additive",
    experiment = "hla_alleles"
  )

  expect_equal(getExperiments(midas), "hla_alleles")
})

test_that("MiDAS object's hla_calls is extracted correctly", {
  hla_calls_file <- system.file("extdata", "HLAHD_output_example.txt", package = "MiDAS")
  hla_calls <- readHlaCalls(hla_calls_file)

  pheno_file <- system.file("extdata", "pheno_example.txt", package = "MiDAS")
  pheno <- read.table(pheno_file, header = TRUE, stringsAsFactors = FALSE)

  midas <- prepareMiDAS(
    hla_calls = hla_calls,
    colData = pheno,
    inheritance_model = "additive",
    experiment = "hla_alleles"
  )

  expect_equal(getHlaCalls(midas), hla_calls)
})

test_that("MiDAS object's kir_calls is extracted correctly", {
  kir_calls_file <- system.file("extdata", "KPI_output_example.txt", package = "MiDAS")
  kir_calls <- readKPICalls(kir_calls_file)
  kir_calls <- kir_calls
  pheno_file <- system.file("extdata", "pheno_example.txt", package = "MiDAS")
  pheno <- read.table(pheno_file, header = TRUE, stringsAsFactors = FALSE)

  midas <- prepareMiDAS(
    kir_calls = kir_calls,
    colData = pheno,
    inheritance_model = "additive",
    experiment = "kir_genes"
  )

  expect_equal(getKirCalls(midas), kir_calls)
})

test_that("MiDAS object's placeholder is extracted correctly", {
  kir_calls_file <- system.file("extdata", "KPI_output_example.txt", package = "MiDAS")
  kir_calls <- readKPICalls(kir_calls_file)
  kir_calls <- kir_calls[1:20, ]

  pheno_file <- system.file("extdata", "pheno_example.txt", package = "MiDAS")
  pheno <- read.table(pheno_file, header = TRUE, stringsAsFactors = FALSE)

  midas <- prepareMiDAS(
    kir_calls = kir_calls,
    colData = pheno,
    inheritance_model = "additive",
    experiment = "kir_genes"
  )

  expect_equal(getPlaceholder(midas), "term")
})

test_that("MiDAS object's omnibus groups are extracted correctly", {
  hla_calls_file <- system.file("extdata", "HLAHD_output_example.txt", package = "MiDAS")
  hla_calls <- readHlaCalls(hla_calls_file)

  pheno_file <- system.file("extdata", "pheno_example.txt", package = "MiDAS")
  pheno <- read.table(pheno_file, header = TRUE, stringsAsFactors = FALSE)

  midas <- prepareMiDAS(
    hla_calls = hla_calls,
    colData = pheno,
    inheritance_model = "additive",
    experiment = "hla_aa"
  )

  expect_equal(getOmnibusGroups(midas, "hla_aa")[1:3],
               list(
                 `A_-15` = c("A_-15_V", "A_-15_L"),
                 `A_-11` = c("A_-11_S", "A_-11_L"),
                 A_9 = c("A_9_F", "A_9_Y", "A_9_S")
               ))

  expect_error(
    getOmnibusGroups(midas, 1),
    "experiment is not a string \\(a length one character vector\\)."
  )
})

test_that("MiDAS object's frequencies are extracted correctly", {
  hla_calls_file <- system.file("extdata", "HLAHD_output_example.txt", package = "MiDAS")
  hla_calls <- readHlaCalls(hla_calls_file)[1:5, 1:5]

  pheno_file <- system.file("extdata", "pheno_example.txt", package = "MiDAS")
  pheno <- read.table(pheno_file, header = TRUE, stringsAsFactors = FALSE)[1:5, ]

  midas <- prepareMiDAS(
    hla_calls = hla_calls,
    colData = pheno,
    inheritance_model = "additive",
    experiment = c("hla_alleles", "hla_divergence")
  )

  freq <- getFrequencies(midas, "hla_alleles")
  test_freq <- data.frame(
    term = c(
      "A*01:01",
      "A*02:01",
      "A*02:06",
      "A*26:01",
      "B*07:02",
      "B*08:01",
      "B*13:02",
      "B*15:01",
      "B*27:05",
      "B*40:01",
      "B*57:01"
    ),
    Counts = c(2, 5, 1, 2, 2, 2, 1, 1, 1, 2, 1),
    Freq = formattable::percent(c(
      0.2, 0.5, 0.1, 0.2, 0.2, 0.2, 0.1, 0.1, 0.1, 0.2, 0.1
    )),
    stringsAsFactors = FALSE
  )
  rownames(test_freq) <-
    c(
      "A*01:01",
      "A*02:01",
      "A*02:06",
      "A*26:01",
      "B*07:02",
      "B*08:01",
      "B*13:02",
      "B*15:01",
      "B*27:05",
      "B*40:01",
      "B*57:01"
    )

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
  kir_calls_file <- system.file("extdata", "KPI_output_example.txt", package = "MiDAS")
  kir_calls <- readKPICalls(kir_calls_file)
  kir_calls <- kir_calls[1:20, ]

  pheno_file <- system.file("extdata", "pheno_example.txt", package = "MiDAS")
  pheno <- read.table(pheno_file, header = TRUE, stringsAsFactors = FALSE)

  midas <- prepareMiDAS(
    kir_calls = kir_calls,
    colData = pheno,
    inheritance_model = "additive",
    experiment = "kir_genes"
  )

  midas_df <- as.data.frame(midas)
  test_midas_df <-
    midasToWide(midas,
                experiment = getExperiments(midas))

  expect_equal(midas_df, test_midas_df)
})

test_that("MiDAS object is prepared properly", {
  rleft_join <- function(init, ..., by = "ID") {
    df <- Reduce(function(...)
      dplyr::left_join(..., by = by),
      x = list(...),
      init = init)
    rownames(df) <- df[[by]]
    df
  }

  hla_calls_file <- system.file("extdata", "HLAHD_output_example.txt", package = "MiDAS")
  hla_calls <- readHlaCalls(hla_calls_file)

  kir_calls_file <- system.file("extdata", "KPI_output_example.txt", package = "MiDAS")
  kir_calls <- readKPICalls(kir_calls_file)

  pheno_file <- system.file("extdata", "pheno_example.txt", package = "MiDAS")
  pheno <- read.table(pheno_file, header = TRUE, stringsAsFactors = FALSE)
  covar_file <- system.file("extdata", "covar_example.txt", package = "MiDAS")
  covar <- read.table(covar_file, header = TRUE, stringsAsFactors = FALSE)
  phenotype <- rleft_join(pheno, covar)

  args_c <-
    expand.grid(
      inheritance_model = c("additive", "dominant", "recessive"),
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
      inheritance_model = args$inheritance_model,
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
    prepareMiDAS(hla_calls = hla_calls, colData = phenotype, inheritance_model = 1),
    "inheritance_model is not a string \\(a length one character vector\\)."
  )

  expect_error(
    prepareMiDAS(hla_calls = hla_calls, colData = phenotype, inheritance_model = "foo"),
    "inheritance_model should be one of \"additive\", \"dominant\", \"recessive\"."
  )

  expect_error(
    prepareMiDAS(
      hla_calls = hla_calls,
      colData = phenotype,
      inheritance_model = "additive",
      experiment = 1
    ),
    "experiment is not a character vector"
  )

  expect_error(
    prepareMiDAS(
      hla_calls = hla_calls,
      colData = phenotype,
      inheritance_model = "additive",
      experiment = "foo"
    ),
    "experiment should match values \"hla_alleles\", \"hla_aa\", \"hla_g_groups\", \"hla_supertypes\", \"hla_NK_ligands\", \"kir_genes\", \"hla_kir_interactions\"."
  )

  expect_error(
    prepareMiDAS(
      hla_calls = hla_calls,
      colData = phenotype,
      inheritance_model = "additive",
      experiment = "hla_alleles",
      placeholder = 1
    ),
    "placeholder is not a string \\(a length one character vector\\)."
  )

  expect_error(
    prepareMiDAS(
      hla_calls = hla_calls,
      colData = phenotype,
      inheritance_model = "additive",
      experiment = "hla_alleles",
      placeholder = "AGE"
    ),
    "Placeholder 'AGE' can not be used, it is alredy used as column name in one of the inputs."
  )
})

test_that("MiDAS data for hla_alleles analysis is prepared properly", {
  hla_calls_file <- system.file("extdata", "HLAHD_output_example.txt", package = "MiDAS")
  hla_calls <- readHlaCalls(hla_calls_file)

  args_c <- expand.grid(
    inheritance_model = c("dominant", "recessive", "additive"),
    stringsAsFactors = FALSE
  )

  invisible(apply(
    X = args_c,
    MARGIN = 1,
    FUN = function(args) {
      args <- as.list(args)
      args$hla_calls <- hla_calls
      experiment <- do.call(prepareMiDAS_hla_alleles, args)

      experiment_test <- do.call(hlaCallsToCounts, args)
      experiment_test <- dfToExperimentMat(experiment_test)

      expect_equal(experiment, experiment_test)
    }
  ))
})

test_that("MiDAS data for hla_aa analysis is prepared properly", { # TODO
  hla_calls_file <- system.file("extdata", "HLAHD_output_example.txt", package = "MiDAS")
  hla_calls <- readHlaCalls(hla_calls_file)

  args_c <- expand.grid(
    inheritance_model = c("dominant", "recessive", "additive"),
    indels = c(TRUE, FALSE),
    unkchar = c(TRUE, FALSE),
    stringsAsFactors = FALSE
  )

  for (i in 1:nrow(args_c)) {
    args <- as.list(args_c[i, , drop = FALSE])
    args$hla_calls <- hla_calls
    experiment <- do.call(prepareMiDAS_hla_aa, args)

    counts <-
      hlaToAAVariation(
        hla_calls = args$hla_calls,
        indels = args$indels,
        unkchar = args$unkchar
      )
    counts <-
      aaVariationToCounts(counts, inheritance_model = args$inheritance_model)
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
  hla_calls_file <- system.file("extdata", "HLAHD_output_example.txt", package = "MiDAS")
  hla_calls <- readHlaCalls(hla_calls_file)

  args_c <- expand.grid(
    inheritance_model = c("dominant", "recessive", "additive"),
    stringsAsFactors = FALSE
  )

  for (i in 1:nrow(args_c)) {
    args <- as.list(args_c[i, , drop = FALSE])
    args$hla_calls <- hla_calls
    experiment <- do.call(prepareMiDAS_hla_g_groups, args)

    experiment_test <-
      hlaToVariable(hla_calls = args$hla_calls,
                    dictionary = "allele_HLA_Ggroup",
                    na.value = 0
      )
    experiment_test <-
      hlaCallsToCounts(
        experiment_test,
        inheritance_model = args$inheritance_model,
        check_hla_format = FALSE
      )
    experiment_test <- dfToExperimentMat(experiment_test)

    expect_equal(experiment, experiment_test)
  }
})

test_that("MiDAS data for hla_supertypes analysis is prepared properly", {
  hla_calls_file <- system.file("extdata", "HLAHD_output_example.txt", package = "MiDAS")
  hla_calls <- readHlaCalls(hla_calls_file)

  args_c <- expand.grid(
    inheritance_model = c("dominant", "recessive", "additive"),
    stringsAsFactors = FALSE
  )

  for (i in 1:nrow(args_c)) {
    args <- as.list(args_c[i, , drop = FALSE])
    args$hla_calls <- hla_calls
    experiment <- do.call(prepareMiDAS_hla_supertypes, args)

    experiment_test <-
      hlaToVariable(hla_calls = args$hla_calls,
                    dictionary = "allele_HLA_supertype",
                    na.value = 0
      )
    experiment_test <-
      hlaCallsToCounts(
        experiment_test,
        inheritance_model = args$inheritance_model,
        check_hla_format = FALSE
      )
    experiment_test <- dfToExperimentMat(experiment_test)
    experiment_test <-
      experiment_test[rownames(experiment_test) != "Unclassified", ]

    expect_equal(experiment, experiment_test)
  }
})

test_that("MiDAS data for hla_NK_ligands analysis is prepared properly", {
  hla_calls_file <- system.file("extdata", "HLAHD_output_example.txt", package = "MiDAS")
  hla_calls <- readHlaCalls(hla_calls_file)

  args_c <- expand.grid(
    inheritance_model = c("dominant", "recessive", "additive"),
    stringsAsFactors = FALSE
  )

  for (i in 1:nrow(args_c)) {
    args <- as.list(args_c[i, , drop = FALSE])
    args$hla_calls <- hla_calls
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
      hlaCallsToCounts(
        experiment_test,
        inheritance_model = args$inheritance_model,
        check_hla_format = FALSE
      )
    experiment_test <- dfToExperimentMat(experiment_test)

    expect_equal(experiment, experiment_test)
  }
})

test_that("MiDAS data for kir_genes analysis is prepared properly", {
  kir_calls_file <- system.file("extdata", "KPI_output_example.txt", package = "MiDAS")
  kir_calls <- readKPICalls(kir_calls_file)

  args_c <- expand.grid(
    "",
    stringsAsFactors = FALSE
  )

  for (i in 1:nrow(args_c)) {
    args <- as.list(args_c[i, , drop = FALSE])
    args$kir_calls <- kir_calls
    experiment <- do.call(prepareMiDAS_kir_genes, args)

    experiment_test <- dfToExperimentMat(kir_calls)

    expect_equal(experiment, experiment_test)
  }
})

test_that("MiDAS data for hla_kir_interactions analysis is prepared properly", {
  hla_calls_file <- system.file("extdata", "HLAHD_output_example.txt", package = "MiDAS")
  hla_calls <- readHlaCalls(hla_calls_file)
  kir_calls_file <- system.file("extdata", "KPI_output_example.txt", package = "MiDAS")
  kir_calls <- readKPICalls(kir_calls_file)

  args_c <- expand.grid(
    "",
    stringsAsFactors = FALSE
  )

  for (i in 1:nrow(args_c)) {
    args <- as.list(args_c[i, , drop = FALSE])
    args$hla_calls <- hla_calls
    args$kir_calls <- kir_calls
    experiment <- do.call(prepareMiDAS_kir_genes, args)

    experiment_test <- getHlaKirInteractions(
      hla_calls = args$hla_calls,
      kir_counts = args$kir_calls
    )
    experiment_test <- dfToExperimentMat(kir_calls)

    expect_equal(experiment, experiment_test)
  }
})

test_that("MiDAS data for hla_divergence analysis is prepared properly", {
  hla_calls_file <- system.file("extdata", "HLAHD_output_example.txt", package = "MiDAS")
  hla_calls <- readHlaCalls(hla_calls_file)

  args_c <- expand.grid(
    "",
    stringsAsFactors = FALSE
  )

  for (i in 1:nrow(args_c)) {
    args <- as.list(args_c[i, , drop = FALSE])
    args$hla_calls <- hla_calls
    experiment <- do.call(prepareMiDAS_hla_divergence, args)

    experiment_test <- hlaCallsGranthamDistance(
      hla_calls = hla_calls,
      genes = c("A", "B", "C")
    )
    experiment_test$ABC_avg <- rowMeans(experiment_test[-1])
    experiment_test <- dfToExperimentMat(experiment_test)

    expect_equal(experiment, experiment_test)
  }

  hla_calls <- hla_calls[, c(1, 8:9)]
  expect_error(prepareMiDAS_hla_divergence(hla_calls),
               "Grantham distance can be calculated only for class I HLA alleles \\(A, B, C\\)."
  )
})

test_that("MiDAS is filtered correctly", {
  hla_calls_file <- system.file("extdata", "HLAHD_output_example.txt", package = "MiDAS")
  hla_calls <- readHlaCalls(hla_calls_file)
  kir_calls_file <- system.file("extdata", "KPI_output_example.txt", package = "MiDAS")
  kir_calls <- readKPICalls(kir_calls_file)
  pheno_file <- system.file("extdata", "pheno_example.txt", package = "MiDAS")
  pheno <- read.table(pheno_file, header = TRUE, stringsAsFactors = FALSE)
  midas <- prepareMiDAS(
    hla_calls = hla_calls,
    kir_calls = kir_calls,
    colData = pheno,
    inheritance_model = "additive",
    experiment = c("hla_alleles", "hla_supertypes", "kir_genes", "hla_divergence")
  )

  # filtration works as expected
  experiment <- "hla_alleles"
  inheritance_model <- "additive"
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
      inheritance_model = inheritance_model,
      lower_frequency_cutoff = lower_frequency_cutoff,
      upper_frequency_cutoff = upper_frequency_cutoff
    )
  expect_equal(filtered_midas, test_filtered_midas)

  # unsported experiments raise error
  experiment <- "hla_divergence"
  inheritance_model <- "additive"
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
  hla_calls_file <- system.file("extdata", "HLAHD_output_example.txt", package = "MiDAS")
  hla_calls <- readHlaCalls(hla_calls_file)
  pheno_file <- system.file("extdata", "pheno_example.txt", package = "MiDAS")
  pheno <- read.table(pheno_file, header = TRUE, stringsAsFactors = FALSE)
  midas <- prepareMiDAS(
    hla_calls = hla_calls,
    colData = pheno,
    inheritance_model = "additive",
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
  vars <- c("A_83_G", "A_83_R", "A_90_A", "A_90_D")
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
