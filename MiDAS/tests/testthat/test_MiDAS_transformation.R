context("Transforming MiDAS objects")

test_that("Amino acids variability is infered correctly", {
  hla_calls <- system.file("extdata/HLAHD_output_example.txt", package = "MiDAS")
  hla_calls <- readHlaCalls(hla_calls)
  # hla_calls <- MiDAS_tut_HLA
  aa_variation <-
    hlaToAAVariation(hla_calls,
                     indels = TRUE,
                     unkchar = TRUE,
                     as_df = FALSE)
  load(system.file("extdata", "test_aa_variation.Rdata", package = "MiDAS"))
  expect_equal(aa_variation, test_aa_variation)

  expect_error(hlaToAAVariation(hla_calls, indels = "foo"),
               "indels is not a flag \\(a length one logical vector\\)."
  )

  expect_error(hlaToAAVariation(hla_calls, unkchar = "foo"),
               "unkchar is not a flag \\(a length one logical vector\\)."
  )

  expect_error(hlaToAAVariation(hla_calls, as_df = 1),
               "as_df is not a flag \\(a length one logical vector\\)."
  )
})

test_that("HLA calls table is converted to additional variables", {
  hla_calls <- system.file("extdata/HLAHD_output_example.txt",
                           package = "MiDAS"
  )
  hla_calls <- readHlaCalls(hla_calls)
  hla_supertypes <- hlaToVariable(hla_calls, dictionary = "allele_HLA_supertype")
  test_hla_supertypes <-
    lapply(
      hla_calls[,-1],
      convertAlleleToVariable,
      dictionary = system.file("extdata", "Match_allele_HLA_supertype.txt", package = "MiDAS")
    )
  na_mask <- vapply(test_hla_supertypes, function(x) all(is.na(x)), FUN.VALUE = logical(1))
  test_hla_supertypes <- test_hla_supertypes[! na_mask]
  test_hla_supertypes <- do.call(cbind, test_hla_supertypes)
  colnames(test_hla_supertypes) <-
    paste0("supertype_", colnames(test_hla_supertypes))
  test_hla_supertypes <-
    cbind(hla_calls[, 1, drop = FALSE], test_hla_supertypes, stringsAsFactors = FALSE)
  expect_equal(hla_supertypes, test_hla_supertypes)

  expect_error(
    hlaToVariable(c("A*01:01", "A*02:01"), dictionary = "4digit_supertype"),
    "hla_calls is not a data frame"
  )

  expect_error(
    hlaToVariable(
      hla_calls = hla_calls,
      dictionary = "4digit_supertype",
      reduce = "yes"
    ),
    "reduce is not a flag \\(a length one logical vector\\)."
  )

  expect_error(
    hlaToVariable(
      hla_calls = hla_calls,
      dictionary = "4digit_supertype",
      na.value = 1:5
    ),
    "na.value length must equal 1."
  )

  expect_error(
    hlaToVariable(
      hla_calls = hla_calls,
      dictionary = "4digit_supertype",
      nacols.rm = "yes"
    ),
    "nacols.rm is not a flag \\(a length one logical vector\\)."
  )
})

test_that("HLA calls table is converted to counts table", {
  hla_calls <- system.file("extdata/HLAHD_output_example.txt",
                           package = "MiDAS"
  )
  hla_calls <- readHlaCalls(hla_calls)
  load(system.file("extdata", "test_hla_counts.RData", package = "MiDAS"))

  hla_counts <- hlaCallsToCounts(hla_calls, inheritance_model = "dominant")
  expect_equal(hla_counts, test_hla_counts[["dominant"]])

  hla_counts <- hlaCallsToCounts(hla_calls, inheritance_model = "recessive")
  expect_equal(hla_counts, test_hla_counts[["recessive"]])

  hla_counts <- hlaCallsToCounts(hla_calls, inheritance_model = "additive")
  expect_equal(hla_counts, test_hla_counts[["additive"]])

  # check if non HLA data frame can be properly processed
  nonhla <-
    data.frame(
      ID = 1:2,
      A_1 = c("a1", "a2"),
      A_2 = c("a2", "a2"),
      stringsAsFactors = FALSE
    )
  nonhla_counts <-
    hlaCallsToCounts(nonhla,
                     inheritance_model = "recessive",
                     check_hla_format = FALSE)
  expect_equal(nonhla_counts,
               data.frame(
                 ID = 1:2,
                 a1 = c(0, 0),
                 a2 = c(0, 1),
                 stringsAsFactors = FALSE
               ))

  expect_error(
    hlaCallsToCounts(c("A*01:01", "A*02:01"), inheritance_model = "additive"),
    "hla_calls is not a data frame"
  )

  expect_error(
    hlaCallsToCounts(hla_calls, inheritance_model = 123),
    "inheritance_model is not a string \\(a length one character vector\\)."
  )

  expect_error(
    hlaCallsToCounts(hla_calls, inheritance_model = "foo"),
    "inheritance_model should be one of 'dominant', 'recessive', 'additive'"
  )

  expect_error(
    hlaCallsToCounts(
      hla_calls,
      inheritance_model = "additive",
      check_hla_format = 1
    ),
    "check_hla_format is not a flag \\(a length one logical vector\\)."
  )
})

test_that("hla frequencies are calculated properly", {
  minimal_hla_calls <- data.frame(
    ID = c("P1", "P2"),
    A_1 = c("A*01:01", "A*02:01"),
    A_2 = c("A*02:01", "A*01:01"),
    stringsAsFactors = FALSE
  )
  hla_freq <- getHlaFrequencies(minimal_hla_calls)
  test_hla_freq <- data.frame(
    allele = c("A*01:01", "A*02:01"),
    Freq = c(0.5, 0.5),
    stringsAsFactors = FALSE
  )
  expect_equal(hla_freq, test_hla_freq)

  # checkHlaCallsFormat tests are ommited here
})

test_that("amino acids variation data frame is converted to counts table", {
  minimal_hla_calls <- data.frame(
    ID = c("P1", "P2"),
    A_1 = c("A*01:01", "A*02:01"),
    A_2 = c("A*02:01", "A*01:01"),
    stringsAsFactors = FALSE
  )
  aa_var <- hlaToAAVariation(minimal_hla_calls)[, 1:5]

  aa_counts <- aaVariationToCounts(aa_var, inheritance_model = "additive")
  test_aa_counts <- data.frame(
    ID = c("P1", "P2"),
    `A_-15_L` = c(1, 1),
    `A_-15_V` = c(1, 1),
    A_44_K = c(1, 1),
    A_44_R = c(1, 1),
    stringsAsFactors = FALSE, check.names = FALSE
  )
  expect_equal(aa_counts, test_aa_counts)

  aa_counts <- aaVariationToCounts(aa_var, inheritance_model = "dominant")
  test_aa_counts <- data.frame(
    ID = c("P1", "P2"),
    `A_-15_L` = c(1, 1),
    `A_-15_V` = c(1, 1),
    A_44_K = c(1, 1),
    A_44_R = c(1, 1),
    stringsAsFactors = FALSE, check.names = FALSE
  )
  expect_equal(aa_counts, test_aa_counts)

  aa_counts <- aaVariationToCounts(aa_var, inheritance_model = "recessive")
  test_aa_counts <- data.frame(
    ID = c("P1", "P2"),
    `A_-15_L` = c(0, 0),
    `A_-15_V` = c(0, 0),
    A_44_K = c(0, 0),
    A_44_R = c(0, 0),
    stringsAsFactors = FALSE, check.names = FALSE
  )
  expect_equal(aa_counts, test_aa_counts)

  expect_error(
    aaVariationToCounts(c("x", "c"), inheritance_model = "additive"),
    "aa_variation is not a data frame"
  )

  expect_error(
    aaVariationToCounts(aa_var[, -1], inheritance_model = "additive"),
    "first column of aa_variation must be named ID"
  )

  expect_error(
    aaVariationToCounts(aa_var, inheritance_model = 123),
    "inheritance_model is not a string \\(a length one character vector\\)."
  )

  expect_error(
    aaVariationToCounts(aa_var, inheritance_model = "foo"),
    "inheritance_model should be one of 'dominant', 'recessive', 'additive'"
  )
})

test_that("amino acids frequencies are calculated properly", {
  minimal_hla_calls <- data.frame(
    ID = c("P1", "P2"),
    A_1 = c("A*01:01", "A*02:01"),
    A_2 = c("A*02:01", "A*01:01"),
    stringsAsFactors = FALSE
  )
  aa_var <- hlaToAAVariation(minimal_hla_calls)[, 1:5]
  aa_freq <- getAAFrequencies(aa_var)
  test_aa_freq <- data.frame(
    aa_pos = c("A_-15_L", "A_-15_V", "A_44_K", "A_44_R"),
    Freq = c(0.5, 0.5, 0.5, 0.5),
    stringsAsFactors = FALSE
  )
  expect_equal(aa_freq, test_aa_freq)

  expect_error(
    aaVariationToCounts(c("x", "c"), inheritance_model = "additive"),
    "aa_variation is not a data frame"
  )

  expect_error(
    aaVariationToCounts(aa_var[, -1], inheritance_model = "additive"),
    "first column of aa_variation must be named ID"
  )
})

test_that("hla counts table can be reverted to hla calls", {
  hla_calls <- data.frame(ID = c("PAT1", "PAT2", "PAT3"),
                          A_1 = c("A*02:01", "A*02:01", "A*01:01"),
                          A_2 = c("A*02:01", "A*02:06", "A*24:02"),
                          B_1 = c("B*13:02", "B*15:01", "B*13:02"),
                          B_2 = c("B*40:01", "B*57:01", "B*57:01"),
                          stringsAsFactors = FALSE
  )
  hla_counts <- hlaCallsToCounts(hla_calls, inheritance_model = "additive")
  expect_equal(countsToHlaCalls(hla_counts), hla_calls)

  err_hla_counts <- hla_counts
  colnames(err_hla_counts) <- NULL
  expect_error(
    countsToHlaCalls(err_hla_counts),
    "count table has no column names"
  )

  err_hla_counts <- hla_counts
  colnames(err_hla_counts)[2] <- NA
  expect_error(
    countsToHlaCalls(err_hla_counts),
    "column names contains NA values"
  )

  err_hla_counts <- hla_counts
  colnames(err_hla_counts) <- make.names(colnames(err_hla_counts))
  expect_error(
    countsToHlaCalls(err_hla_counts),
    "counts table column names contains improperly formated HLA alleles numbers"
  )

  err_hla_counts <- hla_counts
  err_hla_counts[2, 2] <- 1.5
  expect_error(
    countsToHlaCalls(err_hla_counts),
    "counts can only take values 0, 1 or 2"
  )

  err_hla_counts <- hla_counts
  err_hla_counts[2, 2:3] <- 2
  expect_error(
    countsToHlaCalls(err_hla_counts),
    "some samples have more than two alleles per gene"
  )
})

test_that("results are formatted properly", {
  res <- tibble(term = c("A*01:01", "A*02:01"),
                estimate = c(5, -6),
                p.value = c(0.05, 0.1))
  restab <- formatResults(res,
                          filter_by = c("p.value <= 0.05", "estimate > 0"),
                          arrange_by = c("p.value * estimate"),
                          select_cols = c("allele" = "term", "p.value"),
                          format = "html",
                          header = "informative header")
  correct_tab <-
    dplyr::filter(res, .data[["p.value"]] <= 0.05, .data[["estimate"]] > 0)
  correct_tab <-
    dplyr::arrange(correct_tab, .data[["p.value"]] * .data[["estimate"]])
  correct_tab <-
    dplyr::select(correct_tab, "allele" = .data[["term"]], .data[["p.value"]])
  correct_tab <-
    knitr::kable(correct_tab, format = "html", digits = 50)
  correct_tab <-
    kableExtra::add_header_above(correct_tab, header = c("informative header" = 2))
  correct_tab <-
    kableExtra::kable_styling(correct_tab,
                              bootstrap_options = c("striped", "hover", "condensed"))
  correct_tab <-
    kableExtra::scroll_box(correct_tab, width = "100%", height = "200px")

  expect_equal(restab, correct_tab)

  expect_error(formatResults(res, filter_by = 1),
               "filter_by is not a character vector")

  expect_error(formatResults(res, arrange_by = 1),
               "arrange_by is not a character vector")

  expect_error(formatResults(res, select_cols = 1),
               "select_cols is not a character vector")

  expect_error(formatResults(res, format = 1),
               "format is not a string \\(a length one character vector\\).")

  expect_error(formatResults(res, format = "markdown"),
               "format should be one of \"html\", \"latex\".")

  expect_error(formatResults(res, format ="html", header = 1),
               "header is not a string \\(a length one character vector\\) or NULL.")
})

# TODO
# test_that("results are formatted properly with preselected args", {
#   hla_calls_file <- system.file("extdata", "HLAHD_output_example.txt", package = "MiDAS")
#   hla_calls <- readHlaCalls(hla_calls_file)
#   pheno_file <- system.file("extdata", "pheno_example.txt", package = "MiDAS")
#   pheno <- read.table(pheno_file, header = TRUE, stringsAsFactors = FALSE)
#   covar_file <- system.file("extdata", "covar_example.txt", package = "MiDAS")
#   covar <- read.table(covar_file, header = TRUE, stringsAsFactors = FALSE)
#   phenotype <- dplyr::left_join(pheno, covar, by = "ID")
#   midas <- prepareMiDAS(
#     hla_calls = hla_calls,
#     colData = phenotype,
#     inheritance_model = "dominant",
#     analysis_type = c(
#       "hla_alleles"
#     )
#   )
#
#   object <- stats::glm(OS_DIED ~ 1 + term, data = midas, family = stats::binomial)
#   res <- runMiDAS(object, mode = "linear", analysis_type = "hla_alleles")
#   res_kable <- kableResults(res)
#   res_kable_test <- formatResults(res,
#                                   filter_by = "p.adjusted <= 1",
#                                   arrange_by = "p.value",
#                                   select_cols = c(
#                                     "allele",
#                                     "estimate",
#                                     "std.error",
#                                     "p.value",
#                                     "p.adjusted"
#   ),
#                                   format = "html",
#                                   header = "MiDAS analysis results"
#   )
#
#   expect_equal(res_kable, res_kable_test)
#
#   expect_error(kableResults(res, pvalue_cutoff = "a"),
#                "pvalue_cutoff is not a number \\(a length one numeric vector\\) or NULL."
#   )
#
#   expect_error(kableResults(res, format = 1),
#                "format is not a string \\(a length one character vector\\)."
#   )
#
#   expect_error(kableResults(res, format = "foo"),
#                "format should be one of \"html\", \"latex\"."
#   )
# })

test_that("counts to variables conversion", {
  file <- system.file("extdata", "KPI_output_example.txt", package = "MiDAS")
  kir_counts <- readKPICalls(file)[1:2, ]
  kir_haplotypes <- countsToVariables(kir_counts, "kir_haplotypes")
  kir_haplotypes_test <- data.frame(
    ID = c("PAT1", "PAT2"),
    cenAA = c(1, 0),
    cenBB = c(0, 0),
    cenAB = c(0, 1),
    telAA = c(0, 1),
    telBB = c(0, 0),
    telAB = c(1, 0),
    stringsAsFactors = FALSE
  )
  expect_equal(kir_haplotypes, kir_haplotypes_test)

  # works with a dictionary data frame
  dictionary <- data.frame(
    Name = c("cenAA", "cenAB", "telAA"),
    Expression = c(
      "KIR2DL3 & ! KIR2DL2 & ! KIR2DS2",
      "! (KIR2DL2 & ! KIR2DL3) & ! (KIR2DL3 & ! KIR2DL2 & ! KIR2DS2)",
      "! KIR2DS1 & ! KIR3DS1"
   ),
   stringsAsFactors = FALSE
  )
  kir_haplotypes <- countsToVariables(kir_counts, dictionary)
  expect_equal(kir_haplotypes,
               kir_haplotypes_test[, c("ID", "cenAA", "cenAB", "telAA")]
  )

  expect_error(
    countsToVariables(iris),
    "first column in counts must be named 'ID'"
  )

  expect_error(
    countsToVariables(kir_counts, na.value = 1:2),
    "na.value length must equal 1."
  )

  expect_error(
    countsToVariables(kir_counts, nacols.rm = 1),
    "nacols.rm is not a flag \\(a length one logical vector\\)."
  )

  expect_error(
    countsToVariables(kir_counts, dictionary = "foo"),
    "Path 'foo' does not exist"
  )

  expect_error(
    countsToVariables(kir_counts, dictionary = c("foo", "bar")),
    "dictionary is not a data frame"
  )
})

test_that("HLA - KIR interactions are infered correctly", {
  hla_file <- system.file("extdata", "HLAHD_output_example.txt", package = "MiDAS")
  hla_calls <- readHlaCalls(hla_file)
  kir_file <- system.file("extdata", "KPI_output_example.txt", package = "MiDAS")
  kir_counts <- readKPICalls(kir_file)
  hla_kir <- getHlaKirInteractions(hla_calls, kir_counts)
  load(system.file("extdata", "test_hla_kir_interactions.Rdata", package = "MiDAS"))
  expect_equal(hla_kir, test_hla_kir)

  # checkHlaCallsFormat are omitted here
  # checkKirCallsFormat are omitted here

  expect_error(
    getHlaKirInteractions(hla_calls, kir_counts, interactions_dict = 1),
    "interactions_dict is not a string \\(a length one character vector\\)."
  )

 fake_kir_counts <- kir_counts
 fake_kir_counts[, 1] <- paste0("foo", 1:nrow(kir_counts))
 expect_error(
   getHlaKirInteractions(hla_calls, fake_kir_counts),
   "IDs in hla_calls doesn't match IDs in kir_counts"
 )

 fake_kir_counts <- kir_counts
 fake_kir_counts[1:5, 1] <- paste0("foo", 1:5)
 expect_warning(
  getHlaKirInteractions(hla_calls, fake_kir_counts),
  "14 IDs in hla_calls matched IDs in kir_counts"
 )
})

test_that("Experiments are filtered correctly", {
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

  # filtering works as expected for fractions
  experiment <- midas[["hla_alleles"]]
  inheritance_model <- "additive"
  lower_frequency_cutoff <- 0.54
  upper_frequency_cutoff <- 0.56
  experiment_filtered <- filterExperimentByFrequency(
    experiment = experiment,
    inheritance_model = inheritance_model,
    lower_frequency_cutoff = lower_frequency_cutoff,
    upper_frequency_cutoff = upper_frequency_cutoff
  )
  expected_vars <- c("DPB1*04:01", "DRB4*01:03", "H*01:01", "K*01:02")
  expect_equal(experiment_filtered, experiment[expected_vars, ])

  # filtering works as expected for counts
  experiment <- midas[["hla_alleles"]]
  inheritance_model <- "additive"
  lower_frequency_cutoff <- 6
  upper_frequency_cutoff <- 8
  experiment_filtered <- filterExperimentByFrequency(
    experiment = experiment,
    inheritance_model = inheritance_model,
    lower_frequency_cutoff = lower_frequency_cutoff,
    upper_frequency_cutoff = upper_frequency_cutoff
  )
  expected_vars <-
    c("C*07:02",
      "DQB1*02:01")
  expect_equal(experiment_filtered, experiment[expected_vars, ])

  # filtering works as expected for boundry conditions NULL, NULL
  experiment <- midas[["hla_alleles"]]
  inheritance_model <- "additive"
  lower_frequency_cutoff <- NULL
  upper_frequency_cutoff <- NULL
  experiment_filtered <- filterExperimentByFrequency(
    experiment = experiment,
    inheritance_model = inheritance_model,
    lower_frequency_cutoff = lower_frequency_cutoff,
    upper_frequency_cutoff = upper_frequency_cutoff
  )
  expected_vars <- rownames(experiment)
  expect_equal(experiment_filtered, experiment[expected_vars, ])

  # filtering works as expected for boundry conditions 0, 0
  experiment <- midas[["hla_alleles"]]
  inheritance_model <- "additive"
  lower_frequency_cutoff <- 0
  upper_frequency_cutoff <- 0
  experiment_filtered <- filterExperimentByFrequency(
    experiment = experiment,
    inheritance_model = inheritance_model,
    lower_frequency_cutoff = lower_frequency_cutoff,
    upper_frequency_cutoff = upper_frequency_cutoff
  )
  expected_vars <- character(0L)
  expect_equal(experiment_filtered, experiment[expected_vars, ])

  # filtering works as expected for boundry conditions 1, 1
  experiment <- midas[["hla_alleles"]]
  inheritance_model <- "additive"
  lower_frequency_cutoff <- 1
  upper_frequency_cutoff <- 1
  experiment_filtered <- filterExperimentByFrequency(
    experiment = experiment,
    inheritance_model = inheritance_model,
    lower_frequency_cutoff = lower_frequency_cutoff,
    upper_frequency_cutoff = upper_frequency_cutoff
  )
  expected_vars <- character(0L)
  expect_equal(experiment_filtered, experiment[expected_vars, ])

  # experiment must be a matrix
  experiment <- LETTERS
  inheritance_model <- "additive"
  lower_frequency_cutoff <- 0.1
  upper_frequency_cutoff <- 0.5
  expect_error(
    filterExperimentByFrequency(
      experiment = experiment,
      inheritance_model = inheritance_model,
      lower_frequency_cutoff = lower_frequency_cutoff,
      upper_frequency_cutoff = upper_frequency_cutoff
    ),
    "experiment is not a matrix"
  )

  # experiment must be of type integer
  experiment <- matrix(LETTERS)
  inheritance_model <- "additive"
  lower_frequency_cutoff <- 0.1
  upper_frequency_cutoff <- 0.5
  expect_error(
    filterExperimentByFrequency(
      experiment = experiment,
      inheritance_model = inheritance_model,
      lower_frequency_cutoff = lower_frequency_cutoff,
      upper_frequency_cutoff = upper_frequency_cutoff
    ),
    "values in experiment are not counts \\(a positive integers\\) or zeros."
  )

  # inheritance_model must be a string
  experiment <- midas[["hla_alleles"]]
  inheritance_model <- 1
  lower_frequency_cutoff <- "foo"
  upper_frequency_cutoff <- 0.5
  expect_error(
    filterExperimentByFrequency(
      experiment = experiment,
      inheritance_model = inheritance_model,
      lower_frequency_cutoff = lower_frequency_cutoff,
      upper_frequency_cutoff = upper_frequency_cutoff
    ),
    "inheritance_model is not a string \\(a length one character vector\\)."
  )

  # inheritance_model must match allowed values
  experiment <- midas[["hla_alleles"]]
  inheritance_model <- "foo"
  lower_frequency_cutoff <- 0
  upper_frequency_cutoff <- 1
  inheritance_model_vals <- formals(filterExperimentByFrequency)
  inheritance_model_vals <- inheritance_model_vals[["inheritance_model"]]
  inheritance_model_vals <- eval(inheritance_model_vals)
  inheritance_model_vals <- paste(inheritance_model_vals, collapse = "\", \"")
  inheritance_model_vals <- paste0("\"", inheritance_model_vals, "\"")
  expect_error(
    filterExperimentByFrequency(
      experiment = experiment,
      inheritance_model = inheritance_model,
      lower_frequency_cutoff = lower_frequency_cutoff,
      upper_frequency_cutoff = upper_frequency_cutoff
    ),
    sprintf("inheritance_model should be one of %s.", inheritance_model_vals)
  )

  # lower_frequency_cutof must be a number
  experiment <- midas[["hla_alleles"]]
  inheritance_model <- "additive"
  lower_frequency_cutoff <- "foo"
  upper_frequency_cutoff <- 0.5
  expect_error(
    filterExperimentByFrequency(
      experiment = experiment,
      inheritance_model = inheritance_model,
      lower_frequency_cutoff = lower_frequency_cutoff,
      upper_frequency_cutoff = upper_frequency_cutoff
    ),
    "lower_frequency_cutoff is not a number \\(a length one numeric vector\\)."
  )

  # lower_frequency_cutof must be positive
  experiment <- midas[["hla_alleles"]]
  inheritance_model <- "additive"
  lower_frequency_cutoff <- -1
  upper_frequency_cutoff <- 0.5
  expect_error(
    filterExperimentByFrequency(
      experiment = experiment,
      inheritance_model = inheritance_model,
      lower_frequency_cutoff = lower_frequency_cutoff,
      upper_frequency_cutoff = upper_frequency_cutoff
    ),
    "lower_frequency_cutoff must be a number greater than 0."
  )

  # upper_frequency_cutoff must be a number
  experiment <- midas[["hla_alleles"]]
  inheritance_model <- "additive"
  lower_frequency_cutoff <- 0.5
  upper_frequency_cutoff <- "foo"
  expect_error(
    filterExperimentByFrequency(
      experiment = experiment,
      inheritance_model = inheritance_model,
      lower_frequency_cutoff = lower_frequency_cutoff,
      upper_frequency_cutoff = upper_frequency_cutoff
    ),
    "upper_frequency_cutoff is not a number \\(a length one numeric vector\\)."
  )

  # upper_frequency_cutoff must be positive
  experiment <- midas[["hla_alleles"]]
  inheritance_model <- "additive"
  lower_frequency_cutoff <- 0
  upper_frequency_cutoff <- -1
  expect_error(
    filterExperimentByFrequency(
      experiment = experiment,
      inheritance_model = inheritance_model,
      lower_frequency_cutoff = lower_frequency_cutoff,
      upper_frequency_cutoff = upper_frequency_cutoff
    ),
    "upper_frequency_cutoff must be a number greater than 0."
  )

  # lower_frequency_cutoff is lower than upper_frequency_cutoff
  experiment <- midas[["hla_alleles"]]
  inheritance_model <- "additive"
  lower_frequency_cutoff <- 5
  upper_frequency_cutoff <- 1
  expect_error(
    filterExperimentByFrequency(
      experiment = experiment,
      inheritance_model = inheritance_model,
      lower_frequency_cutoff = lower_frequency_cutoff,
      upper_frequency_cutoff = upper_frequency_cutoff
    ),
    "lower_frequency_cutoff cannot be higher than upper_frequency_cutoff."
  )

  # Both lower_frequency_cutoff and upper_frequency_cutoff have to be either frequencies or counts
  experiment <- midas[["hla_alleles"]]
  inheritance_model <- "additive"
  lower_frequency_cutoff <- 0.5
  upper_frequency_cutoff <- 2
  expect_error(
    filterExperimentByFrequency(
      experiment = experiment,
      inheritance_model = inheritance_model,
      lower_frequency_cutoff = lower_frequency_cutoff,
      upper_frequency_cutoff = upper_frequency_cutoff
    ),
    "Both lower_frequency_cutoff and upper_frequency_cutoff have to be either frequencies or counts."
  )
})

test_that("getExperimentFrequencies", {
  # matirx
  experiment_matrix <- matrix(
    data = c(0, 2, 0, 0, 0,
             0, 2, 0, 0, 0,
             2, 0, 0, 0, 0,
             0, 1, 1, 0, 0,
             0, 0, 0, 0, 0
    ),
    nrow = 5,
    ncol = 5,
    dimnames = list(
      c("A*01:01", "A*02:01", "A*02:06", "A*03:01", "A*23:01"),
      c("PAT1", "PAT2", "PAT3", "PAT4", "PAT5")
    )
  )

  # additive
  experiment_matrix_freq <-
    getExperimentFrequencies(experiment_matrix, "additive")
  experiment_matrix_freq_test <- data.frame(
    term = c("A*01:01", "A*02:01", "A*02:06", "A*03:01", "A*23:01"),
    Counts = c(2, 5, 1, 0, 0),
    Freq = formattable::percent(c(0.2, 0.5, 0.1, 0, 0), 2L),
    row.names = c("A*01:01", "A*02:01", "A*02:06", "A*03:01", "A*23:01"),
    stringsAsFactors = FALSE
  )
  expect_equal(experiment_matrix_freq, experiment_matrix_freq_test)

  # dominant
  experiment_matrix[experiment_matrix == 2] <- 1
  experiment_matrix_freq <-
    getExperimentFrequencies(experiment_matrix, "dominant")
  experiment_matrix_freq_test <- data.frame(
    term = c("A*01:01", "A*02:01", "A*02:06", "A*03:01", "A*23:01"),
    Counts = c(1, 3, 1, 0, 0),
    Freq = formattable::percent(c(0.2, 0.6, 0.2, 0, 0), 2L),
    row.names = c("A*01:01", "A*02:01", "A*02:06", "A*03:01", "A*23:01"),
    stringsAsFactors = FALSE
  )
  expect_equal(experiment_matrix_freq, experiment_matrix_freq_test)

  # SummarizedExperiment
  experiment_se <- SummarizedExperiment::SummarizedExperiment(experiment_matrix)
  experiment_se_freq <- getExperimentFrequencies(experiment_se, "dominant")
  expect_equal(experiment_matrix_freq, experiment_matrix_freq_test)

  expect_error(getExperimentFrequencies(as.matrix(LETTERS)),
               "values in experiment are not counts \\(a positive integers\\) or zeros.")

  expect_error(getExperimentFrequencies(experiment_matrix, 1),
               "inheritance_model is not a string \\(a length one character vector\\).")

  expect_error(getExperimentFrequencies(experiment_matrix, "foo"),
               "inheritance_model should be one of \"dominant\", \"recessive\", \"additive\".")
})
