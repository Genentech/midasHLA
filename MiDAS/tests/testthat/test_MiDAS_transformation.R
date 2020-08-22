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

test_that("hlaCallsToCounts", {
  hla_calls <- data.frame(
    ID = c("PAT1", "PAT2", "PAT3", "PAT4"),
    A_1 = c("A*01", "A*02", NA, "A*01"),
    A_2 = c("A*02", "A*01", "A*02", "A*01"),
    stringsAsFactors = FALSE
  )
  hla_counts <- hlaCallsToCounts(hla_calls)
  test_hla_counts <- data.frame(
    ID = c("PAT1", "PAT2", "PAT3", "PAT4"),
    "A*01" = c(1, 1, 0, 2),
    "A*02" = c(1, 1, 1, 0),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  expect_equal(hla_counts, test_hla_counts)

  # check if non HLA data frame can be properly processed
  nonhla <-
    data.frame(
      ID = 1:2,
      A_1 = c("a1", "a2"),
      A_2 = c("a2", "a2"),
      stringsAsFactors = FALSE
    )
  nonhla_counts <-  hlaCallsToCounts(nonhla, check_hla_format = FALSE)
  expect_equal(nonhla_counts,
               data.frame(
                 ID = 1:2,
                 a1 = c(1, 0),
                 a2 = c(1, 2),
                 stringsAsFactors = FALSE
               ))

  expect_error(
    hlaCallsToCounts(c("A*01:01", "A*02:01")),
    "hla_calls is not a data frame"
  )

  expect_error(
    hlaCallsToCounts(
      hla_calls,
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
    Counts = c(2, 2),
    Freq = formattable::percent(c(0.5, 0.5)),
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

  aa_counts <- aaVariationToCounts(aa_var)
  test_aa_counts <- data.frame(
    ID = c("P1", "P2"),
    `A_-15_L` = c(1, 1),
    `A_-15_V` = c(1, 1),
    A_44_K = c(1, 1),
    A_44_R = c(1, 1),
    stringsAsFactors = FALSE, check.names = FALSE
  )
  expect_equal(aa_counts, test_aa_counts)

  expect_error(
    aaVariationToCounts(c("x", "c")),
    "aa_variation is not a data frame"
  )

  expect_error(
    aaVariationToCounts(aa_var[, -1]),
    "first column of aa_variation must be named ID"
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
    aaVariationToCounts(c("x", "c")),
    "aa_variation is not a data frame"
  )

  expect_error(
    aaVariationToCounts(aa_var[, -1]),
    "first column of aa_variation must be named ID"
  )
})

test_that("formatResults", {
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
    knitr::kable(correct_tab, format = "html", format.args = list(digits = 4, scientific = -5))
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

test_that("kableResults", {
  midas <- prepareMiDAS(
    hla_calls = MiDAS_tut_HLA,
    colData = MiDAS_tut_pheno,
    experiment = "hla_alleles"
  )
  object <- lm(disease ~ term, data = midas)
  res <- runMiDAS(object, experiment = "hla_alleles")
  res_kable <- kableResults(res)
  res_kable_test <- formatResults(res,
                                  filter_by = "p.value <= 1",
                                  arrange_by = "p.value",
                                  select_cols = c("allele",
                                                  "p.value",
                                                  "p.adjusted",
                                                  "estimate",
                                                  "std.error",
                                                  "conf.low",
                                                  "conf.high",
                                                  "statistic",
                                                  "Ntotal",
                                                  "Ntotal [%]" = "Ntotal.percent",
                                                  "N(disease=0)",
                                                  "N(disease=0) [%]" = "N(disease=0).percent",
                                                  "N(disease=1)",
                                                  "N(disease=1) [%]" = "N(disease=1).percent"
                                  ),
                                  format = "html",
                                  header = "MiDAS analysis results"
  )

  expect_equal(res_kable, res_kable_test)

  expect_error(kableResults(LETTERS), "results is not a data frame")

  expect_error(kableResults(res, 1:5), "colnames is not a character vector or NULL.")

  expect_error(kableResults(res, pvalue_cutoff = "a"),
               "pvalue_cutoff is not a number \\(a length one numeric vector\\) or NULL."
  )

  expect_error(kableResults(res, format = 1),
               "format is not a string \\(a length one character vector\\)."
  )

  expect_error(kableResults(res, format = "foo"),
               "format should be one of \"html\", \"latex\"."
  )

  expect_error(
    kableResults(res, c("foo", "bar")),
    "colnames should match values \"allele\", \"p.value\", \"p.adjusted\", \"estimate\", \"std.error\", \"conf.low\", \"conf.high\", \"statistic\", \"Ntotal\", \"Ntotal.percent\", \"N\\(disease=0\\)\", \"N\\(disease=0\\).percent\", \"N\\(disease=1\\)\", \"N\\(disease=1\\).percent\"."
  )
})

test_that("countsToVariables", {
  kir_haplotypes <- countsToVariables(MiDAS_tut_KIR[1:2, ], "kir_haplotypes")
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
  kir_counts <- readKIRCalls(kir_file)
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
  getHlaKirInteractions(MiDAS_tut_HLA, fake_kir_counts),
  "995 IDs in hla_calls matched IDs in kir_calls"
 )
})

test_that("Experiments are filtered correctly", {
  hla_calls <- reduceHlaCalls(MiDAS_tut_HLA, 4)
  midas <- prepareMiDAS(
    hla_calls = hla_calls,
    kir_calls = MiDAS_tut_KIR,
    colData = MiDAS_tut_pheno,
    inheritance_model = "additive",
    experiment = c("hla_alleles", "hla_supertypes", "kir_genes", "hla_divergence")
  )

  # filtering works as expected for fractions
  experiment <- midas[["hla_alleles"]]
  lower_frequency_cutoff <- 0.50
  upper_frequency_cutoff <- 0.90
  experiment_filtered <- filterExperimentByFrequency(
    experiment = experiment,
    lower_frequency_cutoff = lower_frequency_cutoff,
    upper_frequency_cutoff = upper_frequency_cutoff
  )
  expected_vars <- c("DPA1*01:03", "DRA*01:01")
  expect_equal(experiment_filtered, experiment[expected_vars, ])

  # filtering works as expected for counts
  experiment <- midas[["hla_alleles"]]
  lower_frequency_cutoff <- 8
  upper_frequency_cutoff <- 10
  experiment_filtered <- filterExperimentByFrequency(
    experiment = experiment,
    lower_frequency_cutoff = lower_frequency_cutoff,
    upper_frequency_cutoff = upper_frequency_cutoff
  )
  expected_vars <-
    c("DRB1*01:03",
      "B*41:01",
      "C*08:03",
      "DRB1*08:04",
      "B*27:02",
      "B*41:02",
      "DPA1*02:06")
  expect_equal(experiment_filtered, experiment[expected_vars, ])

  # filtering works as expected for boundry conditions NULL, NULL
  experiment <- midas[["hla_alleles"]]
  lower_frequency_cutoff <- NULL
  upper_frequency_cutoff <- NULL
  experiment_filtered <- filterExperimentByFrequency(
    experiment = experiment,
    lower_frequency_cutoff = lower_frequency_cutoff,
    upper_frequency_cutoff = upper_frequency_cutoff
  )
  expected_vars <- rownames(experiment)
  expect_equal(experiment_filtered, experiment[expected_vars, ])

  # filtering works as expected for boundry conditions 0, 0
  experiment <- midas[["hla_alleles"]]
  lower_frequency_cutoff <- 0
  upper_frequency_cutoff <- 0
  experiment_filtered <- filterExperimentByFrequency(
    experiment = experiment,
    lower_frequency_cutoff = lower_frequency_cutoff,
    upper_frequency_cutoff = upper_frequency_cutoff
  )
  expected_vars <- character(0L)
  expect_equal(experiment_filtered, experiment[expected_vars, ])

  # filtering works as expected for boundry conditions 1, 1
  experiment <- midas[["hla_alleles"]]
  lower_frequency_cutoff <- 1
  upper_frequency_cutoff <- 1
  experiment_filtered <- filterExperimentByFrequency(
    experiment = experiment,
    lower_frequency_cutoff = lower_frequency_cutoff,
    upper_frequency_cutoff = upper_frequency_cutoff
  )
  expected_vars <- character(0L)
  expect_equal(experiment_filtered, experiment[expected_vars, ])

  # experiment must be a matrix
  experiment <- LETTERS
  lower_frequency_cutoff <- 0.1
  upper_frequency_cutoff <- 0.5
  expect_error(
    filterExperimentByFrequency(
      experiment = experiment,
      lower_frequency_cutoff = lower_frequency_cutoff,
      upper_frequency_cutoff = upper_frequency_cutoff
    ),
    "experiment is not a matrix"
  )

  # experiment must be of type integer
  experiment <- matrix(LETTERS)
  lower_frequency_cutoff <- 0.1
  upper_frequency_cutoff <- 0.5
  expect_error(
    filterExperimentByFrequency(
      experiment = experiment,
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
      carrier_frequency = 1,
      lower_frequency_cutoff = lower_frequency_cutoff,
      upper_frequency_cutoff = upper_frequency_cutoff
    ),
    "carrier_frequency is not a flag \\(a length one logical vector\\)."
  )

  # lower_frequency_cutof must be a number
  experiment <- midas[["hla_alleles"]]
  lower_frequency_cutoff <- "foo"
  upper_frequency_cutoff <- 0.5
  expect_error(
    filterExperimentByFrequency(
      experiment = experiment,
      lower_frequency_cutoff = lower_frequency_cutoff,
      upper_frequency_cutoff = upper_frequency_cutoff
    ),
    "lower_frequency_cutoff is not a number \\(a length one numeric vector\\)."
  )

  # lower_frequency_cutof must be positive
  experiment <- midas[["hla_alleles"]]
  lower_frequency_cutoff <- -1
  upper_frequency_cutoff <- 0.5
  expect_error(
    filterExperimentByFrequency(
      experiment = experiment,
      lower_frequency_cutoff = lower_frequency_cutoff,
      upper_frequency_cutoff = upper_frequency_cutoff
    ),
    "lower_frequency_cutoff must be a number greater than 0."
  )

  # upper_frequency_cutoff must be a number
  experiment <- midas[["hla_alleles"]]
  lower_frequency_cutoff <- 0.5
  upper_frequency_cutoff <- "foo"
  expect_error(
    filterExperimentByFrequency(
      experiment = experiment,
      lower_frequency_cutoff = lower_frequency_cutoff,
      upper_frequency_cutoff = upper_frequency_cutoff
    ),
    "upper_frequency_cutoff is not a number \\(a length one numeric vector\\)."
  )

  # upper_frequency_cutoff must be positive
  experiment <- midas[["hla_alleles"]]
  lower_frequency_cutoff <- 0
  upper_frequency_cutoff <- -1
  expect_error(
    filterExperimentByFrequency(
      experiment = experiment,
      lower_frequency_cutoff = lower_frequency_cutoff,
      upper_frequency_cutoff = upper_frequency_cutoff
    ),
    "upper_frequency_cutoff must be a number greater than 0."
  )

  # lower_frequency_cutoff is lower than upper_frequency_cutoff
  experiment <- midas[["hla_alleles"]]
  lower_frequency_cutoff <- 5
  upper_frequency_cutoff <- 1
  expect_error(
    filterExperimentByFrequency(
      experiment = experiment,
      lower_frequency_cutoff = lower_frequency_cutoff,
      upper_frequency_cutoff = upper_frequency_cutoff
    ),
    "lower_frequency_cutoff cannot be higher than upper_frequency_cutoff."
  )

  # Both lower_frequency_cutoff and upper_frequency_cutoff have to be either frequencies or counts
  experiment <- midas[["hla_alleles"]]
  lower_frequency_cutoff <- 0.5
  upper_frequency_cutoff <- 2
  expect_error(
    filterExperimentByFrequency(
      experiment = experiment,
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

  # allele frequecny
  experiment_matrix_freq <-
    getExperimentFrequencies(experiment_matrix)
  experiment_matrix_freq_test <- data.frame(
    term = c("A*01:01", "A*02:01", "A*02:06", "A*03:01", "A*23:01"),
    Counts = c(2, 5, 1, 0, 0),
    Freq = formattable::percent(c(0.2, 0.5, 0.1, 0, 0), 2L),
    row.names = c("A*01:01", "A*02:01", "A*02:06", "A*03:01", "A*23:01"),
    stringsAsFactors = FALSE
  )
  expect_equal(experiment_matrix_freq, experiment_matrix_freq_test)

  # carrier frequency
  experiment_matrix_freq <-
    getExperimentFrequencies(experiment_matrix, carrier_frequency =  TRUE)
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
  experiment_se_freq <- getExperimentFrequencies(experiment_se)
  expect_equal(experiment_matrix_freq, experiment_matrix_freq_test)

  expect_error(getExperimentFrequencies(as.matrix(LETTERS)),
               "values in experiment are not counts \\(a positive integers\\) or zeros.")

  expect_error(getExperimentFrequencies(experiment_matrix, 1),
               "carrier_frequency is not a flag \\(a length one logical vector\\).")
})

test_that("applyInheritanceModel", {
  experiment <- matrix(c(2, 0, 1, 1, 2, 1, 0, 2, 2), nrow = 3)

  inheritance_model <- "additive"
  additive <- applyInheritanceModel(experiment, inheritance_model)
  expect_equal(experiment, additive)

  inheritance_model <- "dominant"
  dominant <- applyInheritanceModel(experiment, inheritance_model)
  test_dominant <- matrix(c(1, 0, 1, 1, 1, 1, 0, 1, 1), nrow = 3)
  expect_equal(dominant, test_dominant)

  inheritance_model <- "recessive"
  recessive <- applyInheritanceModel(experiment, inheritance_model)
  test_recessive <- matrix(c(1, 0, 0, 0, 1, 0, 0, 1, 1), nrow = 3)
  expect_equal(recessive, test_recessive)

  se <- SummarizedExperiment(
    assays = list(matrix(c(2, 0, 1, 1, 2, 1, 0, 2, 2), nrow = 3)),
    colData = data.frame(foo = 1:3, row.names = 1:3)
  )
  se_dominant <- applyInheritanceModel(se, "dominant")
  test_se <- SummarizedExperiment(
    assays = list(matrix(c(1, 0, 1, 1, 1, 1, 0, 1, 1), nrow = 3)),
    colData = data.frame(foo = 1:3, row.names = 1:3)
  )
  expect_equal(se, se_dominant)
})
