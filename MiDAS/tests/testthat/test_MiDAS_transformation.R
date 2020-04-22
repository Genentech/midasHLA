context("Transforming MiDAS objects")

test_that("Amino acids variability is infered correctly", {
  hla_calls <- system.file("extdata/HLAHD_output_example.txt", package = "MiDAS")
  hla_calls <- readHlaCalls(hla_calls)
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
    A_10_L = c(1, 1),
    A_10_V = c(1, 1),
    A_68_K = c(1, 1),
    A_68_R = c(1, 1),
    stringsAsFactors = FALSE
  )
  expect_equal(aa_counts, test_aa_counts)

  aa_counts <- aaVariationToCounts(aa_var, inheritance_model = "dominant")
  test_aa_counts <- data.frame(
    ID = c("P1", "P2"),
    A_10_L = c(1, 1),
    A_10_V = c(1, 1),
    A_68_K = c(1, 1),
    A_68_R = c(1, 1),
    stringsAsFactors = FALSE
  )
  expect_equal(aa_counts, test_aa_counts)

  aa_counts <- aaVariationToCounts(aa_var, inheritance_model = "recessive")
  test_aa_counts <- data.frame(
    ID = c("P1", "P2"),
    A_10_L = c(0, 0),
    A_10_V = c(0, 0),
    A_68_K = c(0, 0),
    A_68_R = c(0, 0),
    stringsAsFactors = FALSE
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
    aa_pos = c("A_10_L", "A_10_V", "A_68_K", "A_68_R"),
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

test_that("counts are conveerted into frequencies", {
  file <- system.file("extdata", "HLAHD_output_example.txt", package = "MiDAS")
  hla_calls <- readHlaCalls(file)
  hla_counts <- hlaCallsToCounts(hla_calls, inheritance_model = "additive")
  hla_freq <- getCountsFrequencies(hla_counts, inheritance_model = "additive")
  hla_freq <- hla_freq[, c("term", "Freq")]
  colnames(hla_freq) <- c("allele", "Freq")
  rownames(hla_freq) <- NULL
  test_hla_freq <- getHlaFrequencies(hla_calls)
  if (all(class(test_hla_freq$Freq) == c("formattable", "numeric"))) {
    stop("other freq functions are already updated delete this fix")
  }
  test_hla_freq$Freq <- formattable::percent(test_hla_freq$Freq)
  expect_equal(hla_freq, test_hla_freq)

  expect_error(getCountsFrequencies("foo"), "counts_table is not a data frame")

  expect_error(getCountsFrequencies(hla_counts[-1]),
               "first column of counts_table must be named ID")

  expect_error(getCountsFrequencies(hla_calls, inheritance_model = 1),
               "inheritance_model is not a string \\(a length one character vector\\).")

  expect_error(getCountsFrequencies(hla_calls, inheritance_model = "foo"),
               "inheritance_model should be one of \"dominant\", \"recessive\", \"additive\".")

  expect_error(getCountsFrequencies(hla_calls, inheritance_model = "additive"),
               "values in counts_table are not counts \\(a positive integers\\) or zeros.")
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
#       "hla_allele"
#     )
#   )
#
#   object <- stats::glm(OS_DIED ~ 1 + term, data = midas, family = stats::binomial)
#   res <- runMiDAS(object, mode = "linear", analysis_type = "hla_allele")
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
  file <- system.file("extdata", "KIP_output_example.txt", package = "MiDAS")
  kir_counts <- readKirCalls(file)[1:2, ]
  kir_haplotypes <- countsToVariables(kir_counts, "kir_haplotypes")
  kir_haplotypes_test <- data.frame(
    ID = c("PAT1", "PAT2"),
    cenAA = c(0, 1),
    cenBB = c(0, 0),
    cenAB = c(1, 0),
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
  kir_haplotypes <- countsToVariables(kir_counts, dictionary)
  expect_equal(kir_haplotypes,
               kir_haplotypes_test[, c("ID", "cenAA", "cenAB", "telAA")]
  )

  expect_error(
    countsToVariables(iris),
    "counts can't contain factors"
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
  kir_file <- system.file("extdata", "KIP_output_example.txt", package = "MiDAS")
  kir_counts <- readKirCalls(kir_file, counts = TRUE)
  hla_kir <- getHlaKirInteractions(hla_calls, kir_counts)
  load(system.file("extdata", "test_hla_kir_interactions.Rdata", package = "MiDAS"))
  expect_equal(hla_kir, test_hla_kir)

  # checkHlaCallsFormat are omitted here
  # checkKirCountsFormat are omitted here

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
  "15 IDs in hla_calls matched IDs in kir_counts"
 )
})
