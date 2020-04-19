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

  kir_calls_file <- system.file("extdata", "KIP_output_example.txt", package = "MiDAS")
  kir_calls <- readKirCalls(kir_calls_file, counts = TRUE)

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
    phenotype = pheno,
    inheritance_model = "additive",
    analysis_type = character()
  )

  expect_equal(getInheritanceModel(midas), "additive")
})

test_that("MiDAS object's analysis_type is extracted correctly", {
  hla_calls_file <- system.file("extdata", "HLAHD_output_example.txt", package = "MiDAS")
  hla_calls <- readHlaCalls(hla_calls_file)

  pheno_file <- system.file("extdata", "pheno_example.txt", package = "MiDAS")
  pheno <- read.table(pheno_file, header = TRUE, stringsAsFactors = FALSE)

  midas <- prepareMiDAS(
    hla_calls = hla_calls,
    phenotype = pheno,
    inheritance_model = "additive",
    analysis_type = "hla_allele"
  )

  expect_equal(getAnalysisType(midas), "hla_allele")
})

test_that("MiDAS object's hla_calls is extracted correctly", {
  hla_calls_file <- system.file("extdata", "HLAHD_output_example.txt", package = "MiDAS")
  hla_calls <- readHlaCalls(hla_calls_file)

  pheno_file <- system.file("extdata", "pheno_example.txt", package = "MiDAS")
  pheno <- read.table(pheno_file, header = TRUE, stringsAsFactors = FALSE)

  midas <- prepareMiDAS(
    hla_calls = hla_calls,
    phenotype = pheno,
    inheritance_model = "additive",
    analysis_type = character()
  )

  expect_equal(getHlaCalls(midas), hla_calls)
})

test_that("MiDAS object's kir_calls is extracted correctly", {
  kir_calls_file <- system.file("extdata", "KIP_output_example.txt", package = "MiDAS")
  kir_calls <- readKirCalls(kir_calls_file, counts = TRUE)
  kir_calls <- kir_calls[1:20, ]

  pheno_file <- system.file("extdata", "pheno_example.txt", package = "MiDAS")
  pheno <- read.table(pheno_file, header = TRUE, stringsAsFactors = FALSE)

  midas <- prepareMiDAS(
    kir_calls = kir_calls,
    phenotype = pheno,
    inheritance_model = "additive",
    analysis_type = character()
  )

  expect_equal(getKirCalls(midas), kir_calls)
})

test_that("MiDAS's as.data.frame method works properly", {
  kir_calls_file <- system.file("extdata", "KIP_output_example.txt", package = "MiDAS")
  kir_calls <- readKirCalls(kir_calls_file, counts = TRUE)
  kir_calls <- kir_calls[1:20, ]

  pheno_file <- system.file("extdata", "pheno_example.txt", package = "MiDAS")
  pheno <- read.table(pheno_file, header = TRUE, stringsAsFactors = FALSE)

  midas <- prepareMiDAS(
    kir_calls = kir_calls,
    colData = pheno,
    inheritance_model = "additive",
    analysis_type = character()
  )

  midas_df <- as.data.frame(midas)
  test_midas_df <- colData(midas)
  test_midas_df <- as.data.frame(test_midas_df)

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

  kir_calls_file <- system.file("extdata", "KIP_output_example.txt", package = "MiDAS")
  kir_calls <- readKirCalls(kir_calls_file, counts = TRUE)

  pheno_file <- system.file("extdata", "pheno_example.txt", package = "MiDAS")
  pheno <- read.table(pheno_file, header = TRUE, stringsAsFactors = FALSE)
  covar_file <- system.file("extdata", "covar_example.txt", package = "MiDAS")
  covar <- read.table(covar_file, header = TRUE, stringsAsFactors = FALSE)
  phenotype <- rleft_join(pheno, covar)

  args_c <-
    expand.grid(
      inheritance_model = c("additive", "dominant", "recessive"),
      stringsAsFactors = FALSE
    )

  for (i in 1:nrow(args_c)) {
    args <- as.list(args_c[i, , drop = FALSE])

    midas <- prepareMiDAS(
      hla_calls = hla_calls,
      kir_calls = kir_calls,
      phenotype = phenotype,
      inheritance_model = args$inheritance_model,
      analysis_type = c(
        "hla_allele",
        "aa_level",
        "allele_g_group",
        "allele_supertype",
        "allele_group",
        "kir_genes",
        "hla_kir_interactions"
      )
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
    prepareMiDAS(hla_calls = hla_calls, phenotype = cars[1:2, ]),
    "first column in phenotype must be named 'ID'"
  )

  expect_error(
    prepareMiDAS(hla_calls = hla_calls, phenotype = phenotype, inheritance_model = 1),
    "inheritance_model is not a string \\(a length one character vector\\)."
  )

  expect_error(
    prepareMiDAS(hla_calls = hla_calls, phenotype = phenotype, inheritance_model = "foo"),
    "inheritance_model should be one of \"additive\", \"dominant\", \"recessive\"."
  )

  expect_error(
    prepareMiDAS(hla_calls = hla_calls, phenotype = phenotype, inheritance_model = "additive", analysis_type = 1),
    "analysis_type is not a character vector"
  )

  expect_error(
    prepareMiDAS(hla_calls = hla_calls, phenotype = phenotype, inheritance_model = "additive", analysis_type = "foo"),
    "analysis_type should match values \"hla_allele\", \"aa_level\", \"allele_g_group\", \"allele_supertype\", \"allele_group\", \"kir_genes\", \"hla_kir_interactions\"."
  )
})

test_that("MiDAS data for hla_allele analysis is prepared properly", {
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
      experiment <- do.call(prepareMiDAS_hla_allele, args)

      experiment_test <- do.call(hlaCallsToCounts, args)
      experiment_test <- dfToExperimentMat(experiment_test)

      expect_equal(experiment, experiment_test)
    }
  ))
})

test_that("MiDAS data for aa_level analysis is prepared properly", {
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
    experiment <- do.call(prepareMiDAS_aa_level, args)

    experiment_test <-
      hlaToAAVariation(
        hla_calls = args$hla_calls,
        indels = args$indels,
        unkchar = args$unkchar
      )
    experiment_test <-
      aaVariationToCounts(experiment_test, inheritance_model = args$inheritance_model)
    experiment_test <- dfToExperimentMat(experiment_test)

    expect_equal(experiment, experiment_test)
  }
})

test_that("MiDAS data for allele_g_group analysis is prepared properly", {
  hla_calls_file <- system.file("extdata", "HLAHD_output_example.txt", package = "MiDAS")
  hla_calls <- readHlaCalls(hla_calls_file)

  args_c <- expand.grid(
    inheritance_model = c("dominant", "recessive", "additive"),
    stringsAsFactors = FALSE
  )

  for (i in 1:nrow(args_c)) {
    args <- as.list(args_c[i, , drop = FALSE])
    args$hla_calls <- hla_calls
    experiment <- do.call(prepareMiDAS_allele_g_group, args)

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

test_that("MiDAS data for allele_supertype analysis is prepared properly", {
  hla_calls_file <- system.file("extdata", "HLAHD_output_example.txt", package = "MiDAS")
  hla_calls <- readHlaCalls(hla_calls_file)

  args_c <- expand.grid(
    inheritance_model = c("dominant", "recessive", "additive"),
    stringsAsFactors = FALSE
  )

  for (i in 1:nrow(args_c)) {
    args <- as.list(args_c[i, , drop = FALSE])
    args$hla_calls <- hla_calls
    experiment <- do.call(prepareMiDAS_allele_supertype, args)

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

test_that("MiDAS data for allele_group analysis is prepared properly", {
  hla_calls_file <- system.file("extdata", "HLAHD_output_example.txt", package = "MiDAS")
  hla_calls <- readHlaCalls(hla_calls_file)

  args_c <- expand.grid(
    inheritance_model = c("dominant", "recessive", "additive"),
    stringsAsFactors = FALSE
  )

  for (i in 1:nrow(args_c)) {
    args <- as.list(args_c[i, , drop = FALSE])
    args$hla_calls <- hla_calls
    experiment <- do.call(prepareMiDAS_allele_group, args)

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
  kir_calls_file <- system.file("extdata", "KIP_output_example.txt", package = "MiDAS")
  kir_calls <- readKirCalls(kir_calls_file, counts = TRUE)

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
  kir_calls_file <- system.file("extdata", "KIP_output_example.txt", package = "MiDAS")
  kir_calls <- readKirCalls(kir_calls_file, counts = TRUE)

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
