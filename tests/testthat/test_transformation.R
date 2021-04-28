context("Transforming MiDAS objects")

test_that("hlaToAAVariation", {
  aa_variation <-
    hlaToAAVariation(MiDAS_tut_HLA,
                     indels = TRUE,
                     unkchar = TRUE,
                     as_df = FALSE)
  load(system.file("extdata", "test_aa_variation.Rdata", package = "midasHLA"))
  expect_equal(aa_variation, test_aa_variation)

  expect_error(hlaToAAVariation(MiDAS_tut_HLA, indels = "foo"),
               "indels is not a flag \\(a length one logical vector\\)."
  )

  expect_error(hlaToAAVariation(MiDAS_tut_HLA, unkchar = "foo"),
               "unkchar is not a flag \\(a length one logical vector\\)."
  )

  expect_error(hlaToAAVariation(MiDAS_tut_HLA, as_df = 1),
               "as_df is not a flag \\(a length one logical vector\\)."
  )
})

test_that("hlaToVariable", {
  hla_supertypes <- hlaToVariable(MiDAS_tut_HLA, dictionary = "allele_HLA_supertype")
  test_hla_supertypes <-
    lapply(
      MiDAS_tut_HLA[,-1],
      convertAlleleToVariable,
      dictionary = system.file("extdata", "Match_allele_HLA_supertype.txt", package = "midasHLA")
    )
  na_mask <- vapply(test_hla_supertypes, function(x) all(is.na(x)), FUN.VALUE = logical(1))
  test_hla_supertypes <- test_hla_supertypes[! na_mask]
  test_hla_supertypes <- do.call(cbind, test_hla_supertypes)
  colnames(test_hla_supertypes) <-
    paste0("supertype_", colnames(test_hla_supertypes))
  test_hla_supertypes <-
    cbind(MiDAS_tut_HLA[, 1, drop = FALSE], test_hla_supertypes, stringsAsFactors = FALSE)
  expect_equal(hla_supertypes, test_hla_supertypes)

  expect_error(
    hlaToVariable(c("A*01:01", "A*02:01"), dictionary = "4digit_supertype"),
    "hla_calls is not a data frame"
  )

  expect_error(
    hlaToVariable(
      hla_calls = MiDAS_tut_HLA,
      dictionary = "4digit_supertype",
      reduce = "yes"
    ),
    "reduce is not a flag \\(a length one logical vector\\)."
  )

  expect_error(
    hlaToVariable(
      hla_calls = MiDAS_tut_HLA,
      dictionary = "4digit_supertype",
      na.value = 1:5
    ),
    "na.value length must equal 1."
  )

  expect_error(
    hlaToVariable(
      hla_calls = MiDAS_tut_HLA,
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

test_that("getHlaFrequencies", {
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
    Freq = c(0.5, 0.5), 
    row.names = NULL,
    stringsAsFactors = FALSE
  )
  expect_equal(hla_freq, test_hla_freq)
})

test_that("getKIRFrequencies", {
  kir_freq <- getKIRFrequencies(MiDAS_tut_KIR)
  test_kir_freq <- data.frame(
    gene = c("KIR3DL3", "KIR2DS2", "KIR2DL2", "KIR2DL3", "KIR2DP1", "KIR2DL1", "KIR3DP1", "KIR2DL4", 
             "KIR3DL1", "KIR3DS1", "KIR2DL5", "KIR2DS3", "KIR2DS5", "KIR2DS4", "KIR2DS1", "KIR3DL2"),
    Counts = c(935, 455, 449, 849, 925, 918, 935, 935, 853, 365, 485, 294, 288, 853, 371, 935),
    Freq = c(1, 0.486631016042781, 0.480213903743316, 0.908021390374332, 0.989304812834225, 0.981818181818182,
             1, 1, 0.912299465240642, 0.390374331550802, 0.518716577540107, 0.314438502673797, 
             0.308021390374332, 0.912299465240642, 0.396791443850267, 1),
    row.names = NULL,
    stringsAsFactors = FALSE
  )
  expect_equal(kir_freq, test_kir_freq)
})

test_that("hlaToAAVariation", {
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
    A_44_K = c(1, 1),
    A_44_R = c(1, 1),
    A_62_Q = c(1, 1),
    A_62_G = c(1, 1),
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

test_that("getAAFrequencies", {
  minimal_hla_calls <- data.frame(
    ID = c("P1", "P2"),
    A_1 = c("A*01:01", "A*02:01"),
    A_2 = c("A*02:01", "A*01:01"),
    stringsAsFactors = FALSE
  )
  aa_var <- hlaToAAVariation(minimal_hla_calls)[, 1:5]
  aa_freq <- getAAFrequencies(aa_var)
  test_aa_freq <- data.frame(
    aa_pos = c("A_44_K", "A_44_R", "A_62_Q", "A_62_G"),
    Counts = c(2, 2, 2, 2),
    Freq = c(0.5, 0.5, 0.5, 0.5),
    row.names = NULL,
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
    knitr::kable(correct_tab, format = "html", format.args = list(digits = 4, scientific = -3))
  correct_tab <-
    kableExtra::add_header_above(correct_tab, header = c("informative header" = 2))
  correct_tab <-
    kableExtra::kable_styling(correct_tab,
                              bootstrap_options = c("striped", "hover", "condensed"))
  correct_tab <-
    kableExtra::scroll_box(correct_tab, width = "100%", height = "400px")

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

  expect_error(formatResults(res, format ="html", scroll_box_height = 1),
               "scroll_box_height is not a string \\(a length one character vector\\).")
})

test_that("kableResults", {
  midas <- prepareMiDAS(hla_calls = MiDAS_tut_HLA,
                        colData = MiDAS_tut_pheno,
                        experiment = "hla_alleles")
  object <- lm(disease ~ term, data = midas)
  res <- runMiDAS(object, inheritance_model = "additive", experiment = "hla_alleles")
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

  expect_error(kableResults(res, scroll_box_height = 1),
               "scroll_box_height is not a string \\(a length one character vector\\)."
  )
})

test_that("countsToVariables", {
  kir_haplotypes <- countsToVariables(MiDAS_tut_KIR[1:2, ], "kir_haplotypes")
  kir_haplotypes_test <- data.frame(
    ID = c("C001", "C002"),
    cenAA = c(1, 1),
    cenBB = c(0, 0),
    cenAB = c(0, 0),
    telAA = c(1, 1),
    telBB = c(0, 0),
    telAB = c(0, 0),
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
  kir_haplotypes <- countsToVariables(MiDAS_tut_KIR[1:2, ], dictionary)
  expect_equal(kir_haplotypes,
               kir_haplotypes_test[, c("ID", "cenAA", "cenAB", "telAA")]
  )

  expect_error(
    countsToVariables(iris),
    "first column in counts must be named 'ID'"
  )

  expect_error(
    countsToVariables(MiDAS_tut_KIR, na.value = 1:2),
    "na.value length must equal 1."
  )

  expect_error(
    countsToVariables(MiDAS_tut_KIR, nacols.rm = 1),
    "nacols.rm is not a flag \\(a length one logical vector\\)."
  )

  expect_error(
    countsToVariables(MiDAS_tut_KIR, dictionary = "foo"),
    "Path 'foo' does not exist"
  )

  expect_error(
    countsToVariables(MiDAS_tut_KIR, dictionary = c("foo", "bar")),
    "dictionary is not a data frame"
  )
})

test_that("getHlaKirInteractions", {
  hla_kir <- getHlaKirInteractions(MiDAS_tut_HLA, MiDAS_tut_KIR)
  load(system.file("extdata", "test_hla_kir_interactions.Rdata", package = "midasHLA"))
  expect_equal(hla_kir, test_hla_kir_interactions)

  expect_error(
    getHlaKirInteractions(MiDAS_tut_HLA, MiDAS_tut_KIR, interactions_dict = 1),
    "interactions_dict is not a string \\(a length one character vector\\)."
  )

 fake_kir_counts <- MiDAS_tut_KIR
 fake_kir_counts[, 1] <- paste0("foo", 1:nrow(fake_kir_counts))
 expect_error(
   getHlaKirInteractions(MiDAS_tut_HLA, fake_kir_counts),
   "IDs in hla_calls doesn't match IDs in kir_calls"
 )

 fake_kir_counts <- MiDAS_tut_KIR
 fake_kir_counts[1:5, 1] <- paste0("foo", 1:5)
 expect_warning(
  getHlaKirInteractions(MiDAS_tut_HLA, fake_kir_counts),
  "995 IDs in hla_calls matched IDs in kir_calls"
 )
})

test_that("filterExperimentByFrequency", {
  # filtering works as expected for fractions
  experiment <- MiDAS_tut_object[["hla_alleles"]]
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
  experiment <- MiDAS_tut_object[["hla_alleles"]]
  lower_frequency_cutoff <- 8
  upper_frequency_cutoff <- 10
  experiment_filtered <- filterExperimentByFrequency(
    experiment = experiment,
    lower_frequency_cutoff = lower_frequency_cutoff,
    upper_frequency_cutoff = upper_frequency_cutoff
  )
  
  # order of rows somehow depends on the R version...
  expected_vars <-
    c("B*27:02",
      "B*41:01",
      "B*41:02",
      "C*08:03",
      "DPA1*02:06",
      "DRB1*01:03",
      "DRB1*08:04")
  expected_vars <- expected_vars[order(expected_vars)]
  experiment_filtered <-
    experiment_filtered[order(rownames(experiment_filtered)),]

  expect_equal(experiment_filtered, experiment[expected_vars, ])

  # filtering works as expected for boundry conditions NULL, NULL
  experiment <- MiDAS_tut_object[["hla_alleles"]]
  lower_frequency_cutoff <- NULL
  upper_frequency_cutoff <- NULL
  experiment_filtered <- filterExperimentByFrequency(
    experiment = experiment,
    lower_frequency_cutoff = lower_frequency_cutoff,
    upper_frequency_cutoff = upper_frequency_cutoff
  )
  expect_equal(experiment_filtered, experiment)

  # filtering works as expected for boundry conditions 0, 0
  experiment <- MiDAS_tut_object[["hla_alleles"]]
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
  experiment <- MiDAS_tut_object[["hla_alleles"]]
  lower_frequency_cutoff <- 1
  upper_frequency_cutoff <- 1
  experiment_filtered <- filterExperimentByFrequency(
    experiment = experiment,
    lower_frequency_cutoff = lower_frequency_cutoff,
    upper_frequency_cutoff = upper_frequency_cutoff
  )
  expected_vars <- character(0L)
  expect_equal(experiment_filtered, experiment[expected_vars, ])
  
  # omnibus groups are filtred correctly
  experiment <- MiDAS_tut_object[["hla_aa"]]
  lower_frequency_cutoff <- 0.86
  upper_frequency_cutoff <- 0.87
  experiment_filtered <- filterExperimentByFrequency(
    experiment = experiment,
    lower_frequency_cutoff = lower_frequency_cutoff,
    upper_frequency_cutoff = upper_frequency_cutoff
  )
  test_experiment <- experiment[c("C_35_R", "C_309_V"), ]

  expect_equal(experiment_filtered, test_experiment)

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
    "Frequency filtration does not support provided experiment."
  )

  # carrier_frequency must be a string
  experiment <- MiDAS_tut_object[["hla_alleles"]]
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
  experiment <- MiDAS_tut_object[["hla_alleles"]]
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
  experiment <- MiDAS_tut_object[["hla_alleles"]]
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
  experiment <- MiDAS_tut_object[["hla_alleles"]]
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
  experiment <- MiDAS_tut_object[["hla_alleles"]]
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
  experiment <- MiDAS_tut_object[["hla_alleles"]]
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
  experiment <- MiDAS_tut_object[["hla_alleles"]]
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
  experiment_se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(experiment_matrix),
    metadata = list(
      inheritance_model_applicable = TRUE,
      pop_mul = 2,
      omnibus_groups = NULL
    )
  )

  # allele frequecny
  experiment_freq_test <- data.frame(
    term = c("A*01:01", "A*02:01", "A*02:06", "A*03:01", "A*23:01"),
    Counts = c(2, 5, 1, 0, 0),
    Freq = c(0.2, 0.5, 0.1, 0, 0),
    row.names = NULL,
    stringsAsFactors = FALSE
  )

  # matrix
  expect_equal(getExperimentFrequencies(experiment_matrix, 2), experiment_freq_test)

  # se
  experiment_freq <-
    getExperimentFrequencies(experiment_se)
  expect_equal(experiment_freq, experiment_freq_test)

  # carrier frequency
  experiment_freq <-
    getExperimentFrequencies(experiment_se, carrier_frequency =  TRUE)
  experiment_freq_test <- data.frame(
    term = c("A*01:01", "A*02:01", "A*02:06", "A*03:01", "A*23:01"),
    Counts = c(1, 3, 1, 0, 0),
    Freq = c(0.2, 0.6, 0.2, 0, 0),
    row.names = NULL,
    stringsAsFactors = FALSE
  )
  expect_equal(experiment_freq, experiment_freq_test)

  expect_error(
    getExperimentFrequencies(matrix(LETTERS)),
    "values in experiment are not counts \\(a positive integers\\) or zeros."
  )

  expect_error(getExperimentFrequencies(experiment_matrix),
               "pop_mul is not a number \\(a length one numeric vector\\)."
  )

  expect_error(getExperimentFrequencies(experiment_matrix, 1, carrier_frequency = 1),
               "carrier_frequency is not a flag \\(a length one logical vector\\).")

  expect_error(getExperimentFrequencies(experiment_matrix, 1, ref = "foo"),
               "ref is not a data frame")
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

  inheritance_model <- "overdominant"
  recessive <- applyInheritanceModel(experiment, inheritance_model)
  test_recessive <- matrix(c(0, 0, 1, 1, 0, 1, 0, 0, 0), nrow = 3)
  expect_equal(recessive, test_recessive)

  # works on summarized experiment
  se <- SummarizedExperiment(
    assays = list(matrix(c(2, 0, 1, 1, 2, 1, 0, 2, 2), nrow = 3)),
    colData = data.frame(foo = 1:3, row.names = 1:3)
  )
  se_dominant <- applyInheritanceModel(se, "dominant")
  assay(se) <- ifelse(assay(se) == 0, 0, 1)
  expect_equal(se_dominant, se)
})

test_that("getFrequencyMask", {
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
  experiment_matrix_freq <-
    getExperimentFrequencies(experiment_matrix, 2)
  mask <-
    getFrequencyMask(
      df = experiment_matrix_freq,
      lower_frequency_cutoff = 0.1,
      upper_frequency_cutoff = 0.5
    )
  expect_equal(mask, "A*01:01")

  mask <-
    getFrequencyMask(
      df = experiment_matrix_freq,
      lower_frequency_cutoff = 0.1,
      upper_frequency_cutoff = 0.51
    )
  expect_equal(mask, c("A*01:01", "A*02:01"))
})

test_that("filterExperimentByVariables", {
  experiment <- matrix(
    data = c(0, 2, 0, 0, 0,
             0, 2, 0, 0, 0,
             2, 0, 0, 0, 0,
             0, 1, 1, 0, 0,
             0, 0, 0, 0, 0),
    nrow = 5,
    ncol = 5,
    dimnames = list(
      c("A*01:01", "A*02:01", "A*02:06", "A*03:01", "A*23:01"),
      c("PAT1", "PAT2", "PAT3", "PAT4", "PAT5")
    )
  )
  experiment_filtered <- filterExperimentByVariables(experiment, c("A*01:01", "A*02:01"))
  expect_equal(experiment_filtered, experiment[1:2, ])

  experiment <- SummarizedExperiment::SummarizedExperiment(experiment)
  metadata(experiment)$omnibus_groups <-
    list(A = c("A*01:01", "A*02:01", "A*02:06", "A*03:01", "A*23:01"))
  experiment_filtered <- filterExperimentByVariables(experiment, c("A*01:01", "A*02:01"))
  test_experiment <- experiment[c("A*01:01", "A*02:01"), ]

  expect_equal(experiment_filtered, test_experiment)
})

test_that("getExperimentPopulationMultiplicator", {
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(matrix(1:5)),
    metadata = list(pop_mul = 2)
  )
  expect_equal(
    getExperimentPopulationMultiplicator(se),
    2
  )

  mat <- matrix(1:5)
  expect_equal(
    getExperimentPopulationMultiplicator(mat),
    NULL
  )
})

test_that("isExperimentInheritanceModelApplicable", {
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(matrix(1:5)),
    metadata = list(inheritance_model_applicable = TRUE)
  )
  expect_equal(
    isExperimentInheritanceModelApplicable(se),
    TRUE
  )

  mat <- matrix(1:5)
  expect_equal(
    isExperimentInheritanceModelApplicable(mat),
    FALSE
  )
})
