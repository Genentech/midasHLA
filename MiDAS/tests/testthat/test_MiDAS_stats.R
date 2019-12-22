context("HLA allele statistical methods")

test_that("HLA allele associations are analyzed properly", {
  hla_calls_file <-
    system.file("extdata", "HLAHD_output_example.txt", package = "MiDAS")
  hla_calls <- readHlaCalls(hla_calls_file)
  pheno_file <- system.file("extdata", "pheno_example.txt", package = "MiDAS")
  pheno <- read.table(pheno_file, header = TRUE, stringsAsFactors = FALSE)
  covar_file <- system.file("extdata", "covar_example.txt", package = "MiDAS")
  covar <- read.table(covar_file, header = TRUE, stringsAsFactors = FALSE)
  midas_data <-
    prepareMiDAS(hla_calls,
                     pheno,
                     covar,
                     analysis_type = "hla_allele",
                     inheritance_model = "additive")

  object <- lm(OS_DIED ~ AGE + SEX + term, data = midas_data)
  #object$call$data <- midas_data
  res <- analyzeAssociations(object,
                             variables = c("A*01:01", "A*02:01"),
                             correction = "BH"
  )

  test_res <- list(
    lm(OS_DIED ~ AGE + SEX + `A*01:01`, data = midas_data),
    lm(OS_DIED ~ AGE + SEX + `A*02:01`, data = midas_data)
  )
  test_res <- do.call("rbind", lapply(test_res, tidy))
  test_res$term <- gsub("`", "", test_res$term)
  test_res <- test_res[test_res$term %in% c("A*01:01", "A*02:01"), ]
  test_res$p.adjusted <- p.adjust(test_res$p.value, "BH")

  expect_equal(as.data.frame(res), as.data.frame(test_res)) # Tibble doesn't respect tollerance https://github.com/tidyverse/tibble/issues/287 or something related mby

  # Tests for checkStatisticalModel errors are ommitted here

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
      n_correction = 1
    ),
    "n_correction must be at least 2."
  )
})

test_that("Stepwise conditional alleles subset selection", {
  library("survival")

  hla_calls_file <-
    system.file("extdata", "HLAHD_output_example.txt", package = "MiDAS")
  hla_calls <- readHlaCalls(hla_calls_file)
  pheno_file <- system.file("extdata", "pheno_example.txt", package = "MiDAS")
  pheno <- read.table(pheno_file, header = TRUE, stringsAsFactors = FALSE)
  covar_file <- system.file("extdata", "covar_example.txt", package = "MiDAS")
  covar <- read.table(covar_file, header = TRUE, stringsAsFactors = FALSE)
  midas_data <-
    prepareMiDAS(hla_calls,
                     pheno,
                     covar,
                     analysis_type = "hla_allele",
                     inheritance_model = "additive")

  object <- coxph(Surv(OS, OS_DIED) ~ AGE + SEX + term, data = midas_data)

  # keep = FALSE
  res <- analyzeConditionalAssociations(object,
                                        variables = c("B*14:02", "DRB1*11:01"),
                                        th = 0.05,
                                        keep = FALSE)
  res <- rapply(res, classes = "numeric", how = "replace", round, digits = 3)

  test_res <- tibble(term = c("B*14:02", "DRB1*11:01"),
                     estimate = c(3.72, 2.612),
                     std.error = c(1.59, 1.069),
                     statistic = c(2.339, 2.442),
                     p.value = c(0.019, 0.015),
                     conf.low = c(0.603, 0.516),
                     conf.high = c(6.838, 4.707),
                     p.adjusted = c(0.039, 0.015),
                     covariates = c("", "B*14:02")
  )

  expect_equal(res, test_res)

  # keep = TRUE
  res <- analyzeConditionalAssociations(object,
                                        variables = c("B*14:02", "DRB1*11:01"),
                                        th = 0.05,
                                        keep = TRUE)
  res <- rapply(res, classes = "numeric", how = "replace", round, digits = 3)
  test_res <- list(
    tibble(term = c("B*14:02", "DRB1*11:01"),
           estimate = c(3.72, 1.956),
           std.error = c(1.59, 0.960),
           statistic = c(2.339, 2.038),
           p.value = c(0.019, 0.042),
           conf.low = c(0.603, 0.075),
           conf.high = c(6.838, 3.838),
           p.adjusted = c(0.039, 0.083),
           covariates = c("", "")
    ),
    tibble(term = "DRB1*11:01",
           estimate = 2.612,
           std.error = 1.069,
           statistic = 2.442,
           p.value = 0.015,
           conf.low = 0.516,
           conf.high = 4.707,
           p.adjusted = 0.015,
           covariates = "B*14:02"
    )
  )
  expect_equal(res, test_res)

  # Tests for checkStatisticalModel errors are ommitted here

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
      n_correction = 1
    ),
    "n_correction must be at least 2."
  )
})

test_that("MiDAS associations are analyzed properly", {
  hla_calls_file <-
    system.file("extdata", "HLAHD_output_example.txt", package = "MiDAS")
  hla_calls <- readHlaCalls(hla_calls_file)
  pheno_file <- system.file("extdata", "pheno_example.txt", package = "MiDAS")
  pheno <- read.table(pheno_file, header = TRUE, stringsAsFactors = FALSE)
  covar_file <- system.file("extdata", "covar_example.txt", package = "MiDAS")
  covar <- read.table(covar_file, header = TRUE, stringsAsFactors = FALSE)
  kir_file <-
    system.file("extdata", "KIP_output_example.txt", package = "MiDAS")
  kir_counts <- readKirCalls(kir_file, counts = TRUE)
  midas_data <-
    prepareMiDAS(
      hla_calls,
      pheno,
      covar,
      kir_counts = kir_counts,
      analysis_type = c(
        "hla_allele",
        "aa_level",
        "expression_level",
        "allele_g_group",
        "allele_supertype",
        "allele_group",
        "kir_genes",
        "hla_kir_interactions",
        "custom"
      ),
      inheritance_model = "additive"
    )

  object <- lm(OS_DIED ~ AGE + SEX + term, data = midas_data)

  # conditional FALSE, analysis_type = "hla_allele", extra variables
  res <- runMiDAS(object,
                  analysis_type = "hla_allele",
                  variables = c("expression_A", "expression_C"),
                  exponentiate = FALSE
  )

  test_variables <- colnames(midas_data[, label(midas_data) == "hla_allele"])
  test_res <- analyzeAssociations(object, variables = c("expression_A", "expression_C", test_variables), exponentiate = FALSE)
  test_variables <- test_res$term[-1:-2] # constant variables are discarded
  test_res <- dplyr::rename(test_res, allele = term)
  variables_freq <- getCountsFrequencies(midas_data[, c("ID", test_variables)])
  test_res$Ntotal <- c(NA, NA, variables_freq$Counts)
  test_res$Ntotal.frequency <- formattable::percent(c(NA, NA, variables_freq$Freq))
  pos <- midas_data$OS_DIED == 1
  pos_freq <- getCountsFrequencies(midas_data[pos, c("ID", test_variables)])
  test_res$Npositive <- c(NA, NA, pos_freq$Counts)
  test_res$Npositive.frequency <- formattable::percent(c(NA, NA, pos_freq$Freq))
  neg_freq <- getCountsFrequencies(midas_data[! pos, c("ID", test_variables)])
  test_res$Nnegative <- c(NA, NA, neg_freq$Counts)
  test_res$Nnegative.frequency <- formattable::percent(c(NA, NA, neg_freq$Freq))

  expect_equal(as.data.frame(res), as.data.frame(test_res)) # Tibble doesn't respect tollerance https://github.com/tidyverse/tibble/issues/287 or something related mby

  # conditional FALSE, analysis_type = "hla_allele", pattern = "^A"
  res <- runMiDAS(object,
                  analysis_type = "hla_allele",
                  pattern = "^A",
                  exponentiate = FALSE
  )

  test_variables <- colnames(midas_data[, label(midas_data) == "hla_allele"])
  test_variables <- grep("^A", test_variables, value = TRUE)
  test_res <- analyzeAssociations(object, variables = test_variables, exponentiate = FALSE)
  test_res <- dplyr::rename(test_res, allele = term)
  variables_freq <- getCountsFrequencies(midas_data[, c("ID", test_variables)])
  test_res$Ntotal <- variables_freq$Counts
  test_res$Ntotal.frequency <- formattable::percent(variables_freq$Freq)
  pos <- midas_data$OS_DIED == 1
  pos_freq <- getCountsFrequencies(midas_data[pos, c("ID", test_variables)])
  test_res$Npositive <- pos_freq$Counts
  test_res$Npositive.frequency <- formattable::percent(pos_freq$Freq)
  neg_freq <- getCountsFrequencies(midas_data[! pos, c("ID", test_variables)])
  test_res$Nnegative <- neg_freq$Counts
  test_res$Nnegative.frequency <- formattable::percent(neg_freq$Freq)

  expect_equal(as.data.frame(res), as.data.frame(test_res))

  ## analysis_type = "hla_allele" variables = NULL
  res <- runMiDAS(object,
                  analysis_type = "hla_allele",
                  variables = NULL,
                  exponentiate = FALSE
  )

  test_variables <- colnames(midas_data[, label(midas_data) == "hla_allele"])
  test_res <- analyzeAssociations(object, variables = test_variables, exponentiate = FALSE)
  test_variables <- test_res$term # constant variables are discarded
  test_res <- dplyr::rename(test_res, allele = term)
  variables_freq <- getCountsFrequencies(midas_data[, c("ID", test_variables)])
  test_res$Ntotal <- variables_freq$Counts
  test_res$Ntotal.frequency <- variables_freq$Freq
  pos <- midas_data$OS_DIED == 1
  pos_freq <- getCountsFrequencies(midas_data[pos, c("ID", test_variables)])
  test_res$Npositive <- pos_freq$Counts
  test_res$Npositive.frequency <- pos_freq$Freq
  neg_freq <- getCountsFrequencies(midas_data[! pos, c("ID", test_variables)])
  test_res$Nnegative <- neg_freq$Counts
  test_res$Nnegative.frequency <- neg_freq$Freq

  expect_equal(as.data.frame(res), as.data.frame(test_res))


  ## analysis_type = "aa_level" variables = NULL
  res <- runMiDAS(object,
                  analysis_type = "aa_level",
                  variables = NULL,
                  exponentiate = FALSE
  )

  test_variables <- colnames(midas_data[, label(midas_data) == "aa_level"])
  test_res <- analyzeAssociations(object, variables = test_variables, exponentiate = FALSE)
  test_variables <- test_res$term # constant variables are discarded
  test_res <- dplyr::rename(test_res, aa = term)
  variables_freq <- getCountsFrequencies(midas_data[, c("ID", test_variables)])
  test_res$Ntotal <- variables_freq$Counts
  test_res$Ntotal.frequency <- variables_freq$Freq
  pos <- midas_data$OS_DIED == 1
  pos_freq <- getCountsFrequencies(midas_data[pos, c("ID", test_variables)])
  test_res$Npositive <- pos_freq$Counts
  test_res$Npositive.frequency <- pos_freq$Freq
  neg_freq <- getCountsFrequencies(midas_data[! pos, c("ID", test_variables)])
  test_res$Nnegative <- neg_freq$Counts
  test_res$Nnegative.frequency <- neg_freq$Freq

  expect_equal(as.data.frame(res), as.data.frame(test_res))

  ## analysis_type = "expression_level" variables = NULL
  res <- runMiDAS(object,
                  analysis_type = "expression_level",
                  variables = NULL,
                  exponentiate = FALSE
  )

  test_variables <-
    colnames(midas_data[, label(midas_data) == "expression_level", drop = FALSE])
  test_res <- analyzeAssociations(object, variables = test_variables, exponentiate = FALSE)
  test_variables <- test_res$term # constant variables are discarded
  test_res <- dplyr::rename(test_res, allele = term)

  expect_equal(as.data.frame(res), as.data.frame(test_res))

  ## analysis_type = "allele_g_group" variables = NULL
  res <- runMiDAS(object,
                  analysis_type = "allele_g_group",
                  variables = NULL,
                  exponentiate = FALSE
  )

  test_variables <- colnames(midas_data[, label(midas_data) == "allele_g_group"])
  test_res <- analyzeAssociations(object, variables = test_variables)
  test_variables <- test_res$term # constant variables are discarded
  test_res <- dplyr::rename(test_res, g.group = term)
  variables_freq <- getCountsFrequencies(midas_data[, c("ID", test_variables)])
  test_res$Ntotal <- variables_freq$Counts
  test_res$Ntotal.frequency <- variables_freq$Freq
  pos <- midas_data$OS_DIED == 1
  pos_freq <- getCountsFrequencies(midas_data[pos, c("ID", test_variables)])
  test_res$Npositive <- pos_freq$Counts
  test_res$Npositive.frequency <- pos_freq$Freq
  neg_freq <- getCountsFrequencies(midas_data[! pos, c("ID", test_variables)])
  test_res$Nnegative <- neg_freq$Counts
  test_res$Nnegative.frequency <- neg_freq$Freq

  expect_equal(as.data.frame(res), as.data.frame(test_res))

  ## analysis_type = "allele_supertype" variables = NULL
  res <- runMiDAS(object,
                  analysis_type = "allele_supertype",
                  variables = NULL,
                  exponentiate = FALSE
  )

  test_variables <- colnames(midas_data[, label(midas_data) == "allele_supertype"])
  test_res <- analyzeAssociations(object, variables = test_variables)
  test_variables <- test_res$term # constant variables are discarded
  test_res <- dplyr::rename(test_res, supertype = term)
  variables_freq <- getCountsFrequencies(midas_data[, c("ID", test_variables)])
  test_res$Ntotal <- variables_freq$Counts
  test_res$Ntotal.frequency <- variables_freq$Freq
  pos <- midas_data$OS_DIED == 1
  pos_freq <- getCountsFrequencies(midas_data[pos, c("ID", test_variables)])
  test_res$Npositive <- pos_freq$Counts
  test_res$Npositive.frequency <- pos_freq$Freq
  neg_freq <- getCountsFrequencies(midas_data[! pos, c("ID", test_variables)])
  test_res$Nnegative <- neg_freq$Counts
  test_res$Nnegative.frequency <- neg_freq$Freq

  expect_equal(as.data.frame(res), as.data.frame(test_res))

  ## analysis_type = "allele_group" variables = NULL
  res <- runMiDAS(object,
                  analysis_type = "allele_group",
                  variables = NULL,
                  exponentiate = FALSE
  )

  test_variables <- colnames(midas_data[, label(midas_data) == "allele_group"])
  test_res <- analyzeAssociations(object, variables = test_variables)
  test_variables <- test_res$term # constant variables are discarded
  test_res <- dplyr::rename(test_res, allele.group = term)
  variables_freq <- getCountsFrequencies(midas_data[, c("ID", test_variables)])
  test_res$Ntotal <- variables_freq$Counts
  test_res$Ntotal.frequency <- variables_freq$Freq
  pos <- midas_data$OS_DIED == 1
  pos_freq <- getCountsFrequencies(midas_data[pos, c("ID", test_variables)])
  test_res$Npositive <- pos_freq$Counts
  test_res$Npositive.frequency <- pos_freq$Freq
  neg_freq <- getCountsFrequencies(midas_data[! pos, c("ID", test_variables)])
  test_res$Nnegative <- neg_freq$Counts
  test_res$Nnegative.frequency <- neg_freq$Freq

  expect_equal(as.data.frame(res), as.data.frame(test_res))

  ## analysis_type = "kir_genes" variables = NULL
  res <- runMiDAS(object,
                  analysis_type = "kir_genes",
                  variables = NULL,
                  exponentiate = FALSE
  )
  test_variables <- colnames(midas_data[, label(midas_data) == "kir_genes"])
  test_res <- analyzeAssociations(object, variables = test_variables)
  test_variables <- test_res$term # constant variables are discarded
  test_res <- dplyr::rename(test_res, kir.gene = term)
  variables_freq <- getCountsFrequencies(midas_data[, c("ID", test_variables)])
  test_res$Ntotal <- variables_freq$Counts
  test_res$Ntotal.frequency <- variables_freq$Freq
  pos <- midas_data$OS_DIED == 1
  pos_freq <- getCountsFrequencies(midas_data[pos, c("ID", test_variables)])
  test_res$Npositive <- pos_freq$Counts
  test_res$Npositive.frequency <- pos_freq$Freq
  neg_freq <- getCountsFrequencies(midas_data[! pos, c("ID", test_variables)])
  test_res$Nnegative <- neg_freq$Counts
  test_res$Nnegative.frequency <- neg_freq$Freq

  expect_equal(as.data.frame(res), as.data.frame(test_res))

  ## analysis_type = "hla_kir_interactions" variables = NULL
  res <- runMiDAS(object,
                  analysis_type = "hla_kir_interactions",
                  variables = NULL,
                  exponentiate = FALSE
  )
  test_variables <-
    colnames(midas_data[, label(midas_data) == "hla_kir_interactions"])
  test_res <- analyzeAssociations(object, variables = test_variables)
  test_variables <- test_res$term # constant variables are discarded
  test_res <- dplyr::rename(test_res, hla.kir.interaction = term)
  variables_freq <- getCountsFrequencies(midas_data[, c("ID", test_variables)])
  test_res$Ntotal <- variables_freq$Counts
  test_res$Ntotal.frequency <- variables_freq$Freq
  pos <- midas_data$OS_DIED == 1
  pos_freq <- getCountsFrequencies(midas_data[pos, c("ID", test_variables)])
  test_res$Npositive <- pos_freq$Counts
  test_res$Npositive.frequency <- pos_freq$Freq
  neg_freq <- getCountsFrequencies(midas_data[! pos, c("ID", test_variables)])
  test_res$Nnegative <- neg_freq$Counts
  test_res$Nnegative.frequency <- neg_freq$Freq

  expect_equal(as.data.frame(res), as.data.frame(test_res))

  # conditional TRUE keep TRUE
  res <- runMiDAS(object,
                  analysis_type = "hla_allele",
                  conditional = TRUE,
                  keep = TRUE,
                  exponentiate = FALSE
  )

  test_variables <- colnames(midas_data[, label(midas_data) == "hla_allele"])
  test_res <- analyzeConditionalAssociations(object,
                                             variables = test_variables,
                                             th = 0.05,
                                             keep = TRUE
  )
  test_res <- lapply(test_res, dplyr::rename, allele = term)
  alleles <- lapply(test_res, `[`, "allele")
  alleles <- unique(unlist(alleles))
  variables_freq <- getCountsFrequencies(midas_data[, c("ID", alleles)])
  colnames(variables_freq) <- c("allele", "Ntotal", "Ntotal.frequency")
  test_res <-
    lapply(test_res, dplyr::left_join, y = variables_freq, by = "allele")
  pos <- midas_data$OS_DIED == 1
  pos_freq <- getCountsFrequencies(midas_data[pos, c("ID", alleles)])
  colnames(pos_freq) <- c("allele", "Npositive", "Npositive.frequency")
  test_res <-
    lapply(test_res, dplyr::left_join, y = pos_freq, by = "allele")
  neg_freq <- getCountsFrequencies(midas_data[! pos, c("ID", alleles)])
  colnames(neg_freq) <- c("allele", "Nnegative", "Nnegative.frequency")
  test_res <-
    lapply(test_res, dplyr::left_join, y = neg_freq, by = "allele")

  res <- lapply(res, as.data.frame)
  test_res <- lapply(test_res, as.data.frame)
  expect_equal(res, test_res) # Tibble doesn't respect tollerance https://github.com/tidyverse/tibble/issues/287 or something related mby

  # conditional TRUE keep FALSE
  res <- runMiDAS(object,
                  analysis_type = "hla_allele",
                  conditional = TRUE,
                  keep = FALSE,
                  exponentiate = FALSE
  )

  test_variables <- colnames(midas_data[, label(midas_data) == "hla_allele"])
  test_res <- analyzeConditionalAssociations(object,
                                             variables = test_variables,
                                             th = 0.05,
                                             keep = FALSE,
                                             exponentiate = FALSE
  )
  test_variables <- test_res$term # constant variables are discarded
  test_res <- dplyr::rename(test_res, allele = term)
  variables_freq <- getCountsFrequencies(midas_data[, c("ID", test_variables)])
  test_res$Ntotal <- variables_freq$Counts
  test_res$Ntotal.frequency <- variables_freq$Freq
  pos <- midas_data$OS_DIED == 1
  pos_freq <- getCountsFrequencies(midas_data[pos, c("ID", test_variables)])
  test_res$Npositive <- pos_freq$Counts
  test_res$Npositive.frequency <- pos_freq$Freq
  neg_freq <- getCountsFrequencies(midas_data[! pos, c("ID", test_variables)])
  test_res$Nnegative <- neg_freq$Counts
  test_res$Nnegative.frequency <- neg_freq$Freq

  expect_equal(as.data.frame(res), as.data.frame(test_res))

  # Test lower and upper frequency thresholds
  # %
  res <- runMiDAS(object, analysis_type = "hla_allele", lower_frequency_cutoff = 0.85)
  freqs <- getHlaFrequencies(hla_calls)
  expect_equal(res$allele, freqs$allele[freqs$Freq > 0.85 & freqs$Freq != 1])

  res <- runMiDAS(object, analysis_type = "hla_allele", upper_frequency_cutoff = 0.03)
  expect_equal(res$allele, freqs$allele[freqs$Freq < 0.03 & freqs$Freq != 1])

  # counts
  counts <- prepareMiDAS(hla_calls, analysis_type = "hla_allele")
  counts <- colSums(counts[-1], na.rm = TRUE)
  res <- runMiDAS(object, analysis_type = "hla_allele", lower_frequency_cutoff = 34)
  expect_equal(res$allele, names(counts)[counts > 34 & freqs$Freq != 1])

  res <- runMiDAS(object, analysis_type = "hla_allele", upper_frequency_cutoff = 2)
  expect_equal(res$allele, names(counts)[counts < 2 & freqs$Freq != 1])

  # Tests for checkStatisticalModel errors are ommitted here

  expect_error(runMiDAS(object, analysis_type = 1),
               "analysis_type is not a string \\(a length one character vector\\)."
  )

  expect_error(runMiDAS(object, analysis_type = "a"),
               "analysis_type should be one of \"hla_allele\", \"aa_level\", \"expression_level\", \"allele_g_group\", \"allele_supertype\", \"allele_group\", \"kir_genes\", \"hla_kir_interactions\"."
  )

  expect_error(runMiDAS(object, analysis_type = "hla_allele", pattern = 1),
               "pattern is not a string \\(a length one character vector\\) or NULL."
  )

  expect_error(runMiDAS(object, analysis_type = "hla_allele", variables = 1),
               "variables is not a character vector or NULL."
  )

  expect_error(runMiDAS(object, analysis_type = "hla_allele", placeholder = 1),
               "placeholder is not a string \\(a length one character vector\\)."
  )

  expect_error(runMiDAS(object, analysis_type = "hla_allele", conditional = 1),
               "conditional is not a flag \\(a length one logical vector\\)."
  )

  expect_error(runMiDAS(object, analysis_type = "hla_allele", keep = 1),
               "keep is not a flag \\(a length one logical vector\\)."
  )

  expect_error(runMiDAS(object, analysis_type = "hla_allele", variables = "thief"),
               "thief can not be found in object data"
  )

  expect_error(runMiDAS(object, analysis_type = "hla_allele", lower_frequency_cutoff = "foo"),
               "lower_frequency_cutoff is not a number \\(a length one numeric vector\\) or NULL."
  )

  expect_error(runMiDAS(object, analysis_type = "hla_allele", upper_frequency_cutoff = "foo"),
               "upper_frequency_cutoff is not a number \\(a length one numeric vector\\) or NULL."
  )

  expect_error(runMiDAS(object, analysis_type = "hla_allele", pvalue_cutoff = "foo"),
               "pvalue_cutoff is not a number \\(a length one numeric vector\\) or NULL."
  )

  expect_error(runMiDAS(object, analysis_type = "hla_allele", correction = NA),
               "correction is not a string \\(a length one character vector\\)."
  )

  expect_error(runMiDAS(object, analysis_type = "hla_allele", n_correction = "foo"),
               "n_correction is not a count \\(a single positive integer\\) or NULL."
  )

  expect_error(runMiDAS(object, analysis_type = "hla_allele", exponentiate = "NA"),
               "exponentiate is not a flag \\(a length one logical vector\\)."
  )

  expect_error(runMiDAS(object, analysis_type = "hla_allele", th = "NA"),
               "th is not a number \\(a length one numeric vector\\)."
  )

  expect_error(runMiDAS(object, analysis_type = "hla_allele", rss_th = "NA"),
               "rss_th is not a number \\(a length one numeric vector\\)."
  )
})

test_that("MiDAS data is prepared properly", {
  rleft_join <- function(init, ..., by = "ID") {
    Reduce(function(...)
      dplyr::left_join(..., by = by),
      x = list(...),
      init = init)
  }

  hla_calls_file <- system.file("extdata", "HLAHD_output_example.txt", package = "MiDAS")
  hla_calls <- readHlaCalls(hla_calls_file)
  pheno_file <- system.file("extdata", "pheno_example.txt", package = "MiDAS")
  pheno <- read.table(pheno_file, header = TRUE, stringsAsFactors = FALSE)
  covar_file <- system.file("extdata", "covar_example.txt", package = "MiDAS")
  covar <- read.table(covar_file, header = TRUE, stringsAsFactors = FALSE)

  # hla_allele
  midas_hla_allele <-
    prepareMiDAS(hla_calls,
                     pheno,
                     covar,
                     analysis_type = "hla_allele",
                     inheritance_model = "additive")
  midas_hla_allele_test <-
    hlaCallsToCounts(hla_calls, inheritance_model = "additive")
  Hmisc::label(midas_hla_allele_test[-1], self = FALSE) <-
    rep("hla_allele", ncol(midas_hla_allele_test) - 1)
  midas_hla_allele_test <- rleft_join(midas_hla_allele_test, pheno, covar)
  midas_hla_allele_test$term <- 1
  expect_equal(midas_hla_allele, midas_hla_allele_test)

  # aa_level
  midas_aa_level <- prepareMiDAS(hla_calls,
                                     pheno,
                                     covar,
                                     analysis_type = "aa_level",
                                     inheritance_model = "additive")
  midas_aa_level_test <- hlaToAAVariation(hla_calls)
  midas_aa_level_test <-
    aaVariationToCounts(midas_aa_level_test, inheritance_model = "additive")
  Hmisc::label(midas_aa_level_test[-1], self = FALSE) <-
    rep("aa_level", ncol(midas_aa_level_test) - 1)
  midas_aa_level_test <- rleft_join(midas_aa_level_test, pheno, covar)
  midas_aa_level_test$term <- 1
  expect_equal(midas_aa_level, midas_aa_level_test)

  # expression_levels
  midas_expression_levels <- prepareMiDAS(hla_calls,
                                              pheno,
                                              covar,
                                              analysis_type = "expression_level",
                                              inheritance_model = "additive")
  expression_dicts <- grep("expression", listMiDASDictionaries(), value = TRUE)
  midas_expression_levels_test <- Reduce(
    f = function(...) dplyr::left_join(..., by = "ID"),
    x = lapply(expression_dicts, function(x) {
      expr <- hlaToVariable(hla_calls = hla_calls, dictionary = x, na.value = NA)
      expr$sum <- rowSums(expr[, -1, drop = FALSE])
      gene <- gsub("_1", "", colnames(expr)[2])
      expr <- expr[, c("ID", "sum")]
      colnames(expr) <- c("ID", gene)
      expr
    })
  )
  Hmisc::label(midas_expression_levels_test[-1], self = FALSE) <-
    rep("expression_level", ncol(midas_expression_levels_test) - 1)
  midas_expression_levels_test <-
    rleft_join(midas_expression_levels_test, pheno, covar)
  midas_expression_levels_test$term <- 1
  expect_equal(midas_expression_levels, midas_expression_levels_test)

  # allele_g_group
  midas_allele_g_group <- prepareMiDAS(hla_calls,
                                          pheno,
                                          covar,
                                          analysis_type = "allele_g_group",
                                          inheritance_model = "additive")
  midas_allele_g_group_test <-
    hlaToVariable(hla_calls, dictionary = "allele_HLA_Ggroup")
  midas_allele_g_group_test <-
    hlaCallsToCounts(
      midas_allele_g_group_test,
      inheritance_model = "additive",
      check_hla_format = FALSE
    )
  Hmisc::label(midas_allele_g_group_test[-1], self = FALSE) <-
    rep("allele_g_group", ncol(midas_allele_g_group_test) - 1)
  midas_allele_g_group_test <- rleft_join(midas_allele_g_group_test, pheno, covar)
  midas_allele_g_group_test$term <- 1
  expect_equal(midas_allele_g_group, midas_allele_g_group_test)

  # allele_supertype
  midas_allele_supertype <- prepareMiDAS(hla_calls,
                                          pheno,
                                          covar,
                                          analysis_type = "allele_supertype",
                                          inheritance_model = "additive")
  test_midas_allele_supertype <-
    hlaToVariable(hla_calls, dictionary = "allele_HLA_supertype")
  test_midas_allele_supertype <-
    hlaCallsToCounts(
      test_midas_allele_supertype,
      inheritance_model = "additive",
      check_hla_format = FALSE
    )
  Hmisc::label(test_midas_allele_supertype[-1], self = FALSE) <-
    rep("allele_supertype", ncol(test_midas_allele_supertype) - 1)
  test_midas_allele_supertype <-
    rleft_join(test_midas_allele_supertype, pheno, covar)
  test_midas_allele_supertype <-
    subset(test_midas_allele_supertype, select = - Unclassified)
  test_midas_allele_supertype$term <- 1
  expect_equal(midas_allele_supertype, test_midas_allele_supertype)

  # allele_groups
  midas_allele_groups <- prepareMiDAS(hla_calls,
                                              pheno,
                                              covar,
                                              analysis_type = "allele_group",
                                              inheritance_model = "additive")
  allele_groups_lib <- c("allele_HLA-B_Bw", "allele_HLA_Bw4+A23+A24+A32", "allele_HLA-C_C1-2")
  test_midas_allele_group <- Reduce(
    f = function(...) left_join(..., by = "ID"),
    x = lapply(allele_groups_lib, hlaToVariable, hla_calls = hla_calls)
  )
  test_midas_allele_group <-
    hlaCallsToCounts(
      test_midas_allele_group,
      inheritance_model = "additive",
      check_hla_format = FALSE
    )
  Hmisc::label(test_midas_allele_group[-1], self = FALSE) <-
    rep("allele_group", ncol(test_midas_allele_group) - 1)
  test_midas_allele_group <-
    rleft_join(test_midas_allele_group, pheno, covar)
  test_midas_allele_group$term <- 1
  expect_equal(midas_allele_groups, test_midas_allele_group)

  # kir_genes
  kir_path <- system.file("extdata", "KIP_output_example.txt", package = "MiDAS")
  kir_counts <- readKirCalls(kir_path, counts = TRUE)
  midas_kir_genes <- prepareMiDAS(hla_calls,
                                      pheno,
                                      covar,
                                      kir_counts = kir_counts,
                                      analysis_type = "kir_genes",
                                      inheritance_model = "additive")
  test_midas_kir_genes <- kir_counts
  Hmisc::label(test_midas_kir_genes[-1], self = FALSE) <-
    rep("kir_genes", ncol(kir_counts) - 1)
  test_midas_kir_genes <-
    rleft_join(hla_calls[, 1, drop = FALSE], test_midas_kir_genes, pheno, covar)
  test_midas_kir_genes$term <- 1
  expect_equal(midas_kir_genes, test_midas_kir_genes)

  # hla_kir_interactions
  midas_hla_kir_interactions <- prepareMiDAS(
    hla_calls,
    pheno,
    covar,
    kir_counts = kir_counts,
    analysis_type = "hla_kir_interactions",
    inheritance_model = "additive"
  )
  test_midas_hla_kir_interactions <-
    getHlaKirInteractions(hla_calls = hla_calls, kir_counts = kir_counts)
  Hmisc::label(test_midas_hla_kir_interactions[-1], self = FALSE) <-
    rep("hla_kir_interactions", ncol(test_midas_hla_kir_interactions) - 1)
  test_midas_hla_kir_interactions <-
    rleft_join(hla_calls[, 1, drop = FALSE],
               test_midas_hla_kir_interactions,
               pheno,
               covar
    )
  test_midas_hla_kir_interactions$term <- 1
  expect_equal(midas_hla_kir_interactions, test_midas_hla_kir_interactions)

  # custom
  midas_custom <- prepareMiDAS(hla_calls,
                                   pheno,
                                   covar,
                                   analysis_type = "custom",
                                   inheritance_model = "additive")
  midas_custom_test <- rleft_join(hla_calls, pheno, covar)
  gene_idx <-
    ! colnames(midas_custom_test) %in% c("ID", "OS", "OS_DIED", "AGE", "SEX")
  Hmisc::label(midas_custom_test[, gene_idx], self = FALSE) <-
    rep("custom", sum(gene_idx))
  midas_custom_test$term <- 1
  expect_equal(midas_custom, midas_custom_test)

  # check more analysis types at once
  midas_multiple <- prepareMiDAS(hla_calls,
                                     pheno,
                                     covar,
                                     kir_counts = kir_counts,
                                     analysis_type = c("hla_allele", "aa_level", "expression_level", "allele_g_group", "allele_supertype", "allele_group", "kir_genes", "custom"),
                                     inheritance_model = "additive")
  midas_multiple_test <-
    rleft_join(
      midas_hla_allele,
      midas_aa_level,
      midas_expression_levels,
      midas_allele_g_group,
      midas_allele_supertype,
      midas_allele_groups,
      midas_kir_genes,
      midas_custom,
      by = c("ID", "OS", "OS_DIED", "AGE", "SEX")
    )
  # order of columns is different a bit expensive but just sort them
  midas_multiple <- midas_multiple[, order(colnames(midas_multiple))]
  midas_multiple_test <-
    midas_multiple_test[, order(colnames(midas_multiple_test))]
  midas_multiple_test[, grepl("term", colnames(midas_multiple_test))] <- NULL
  midas_multiple_test$term <- 1
  expect_equal(midas_multiple, midas_multiple_test)

  # test for checkHlaCallsFormat are ommitted here

  expect_error(
    prepareMiDAS(hla_calls, analysis_type = 1),
    "analysis_type is not a character vector"
  )

  expect_error(
    prepareMiDAS(hla_calls, analysis_type = "foo"),
    "analysis_type should match values \"hla_allele\", \"aa_level\", \"expression_level\", \"allele_g_group\", \"allele_supertype\", \"allele_group\", \"kir_genes\", \"hla_kir_interactions\", \"custom\"."
  )

  expect_error(
    prepareMiDAS(hla_calls, analysis_type = "hla_allele", inheritance_model = 1),
    "inheritance_model is not a string \\(a length one character vector\\)."
  )

  expect_error(
    prepareMiDAS(hla_calls, analysis_type = "hla_allele", inheritance_model = "bar"),
    "inheritance_model should be one of \"dominant\", \"recessive\", \"additive\"."
  )

  expect_error(
    prepareMiDAS(hla_calls, analysis_type = "hla_allele", indels = "no"),
    "indels is not a flag \\(a length one logical vector\\)."
  )

  expect_error(
    prepareMiDAS(hla_calls, analysis_type = "hla_allele", unkchar = "nope"),
    "unkchar is not a flag \\(a length one logical vector\\)."
  )

  # checkAdditionalData on ... argument are ommitted here

  expect_error(
    prepareMiDAS(hla_calls[, c("ID", "DMA_1", "DMA_2")], analysis_type = "expression_level"),
    "no expression levels were found for input hla_calls"
  )

  expect_error(
    prepareMiDAS(hla_calls[, c("ID", "DOB_1", "DOB_2")], analysis_type = "allele_group"),
    "no allele could be assigned to allele groups for input hla_calls"
  )

  expect_error(
    prepareMiDAS(hla_calls, analysis_type = "kir_genes"),
    "\"kir_genes\" analysis type requires kir_counts argument to be specified"
  )

  expect_error(
    prepareMiDAS(hla_calls, analysis_type = "hla_kir_interactions"),
    "\"hla_kir_interactions\" analysis type requires kir_counts argument to be specified"
  )
})

test_that("amino acid omnibus test works fine", {
  hla_calls_file <- system.file("extdata", "HLAHD_output_example.txt", package = "MiDAS")
  hla_calls <- readHlaCalls(hla_calls_file)
  pheno_file <- system.file("extdata", "pheno_example.txt", package = "MiDAS")
  pheno <- read.table(pheno_file, header = TRUE)
  covar_file <- system.file("extdata", "covar_example.txt", package = "MiDAS")
  covar <- read.table(covar_file, header = TRUE)
  midas_data <- prepareMiDAS(hla_calls, pheno, covar, analysis_type = "aa_level")
  object <- lm(OS ~ AGE + SEX + term, data = midas_data)
  omnibus_res <- aaPosOmnibusTest(object, aa_pos = c("B_11", "E_107", "A_246"))

  obj_B11 <- lm(OS ~ AGE + SEX + B_11_A + B_11_S + term, data = midas_data)
  obj_E107 <- lm(OS ~ AGE + SEX + E_107_R + E_107_G + term, data = midas_data)
  obj_A246 <- lm(OS ~ AGE + SEX + A_246_A + A_246_S + term, data = midas_data)
  LRT <- lapply(list(obj_B11, obj_E107, obj_A246), LRTest, mod0 = object)
  omnibus_res_test <- data.frame(
    aa_pos = c("B_11", "E_107", "A_246"),
    residues = c("A, S", "R, G", "A, S"),
    d.f. = c(1, 1, 1),
    statistic = sapply(LRT, `[[`, "statistic"),
    p.value = sapply(LRT, `[[`, "p.value"),
    p.adjusted = p.adjust(sapply(LRT, `[[`, "p.value"), method = "bonferroni"),
    stringsAsFactors = FALSE
  )
  expect_equal(omnibus_res, omnibus_res_test)

  # Tests for checkStatisticalModel errors are ommitted here

  expect_error(aaPosOmnibusTest(object, aa_pos = 1:5),
               "aa_pos is not a character vector"
  )

  expect_error(
    aaPosOmnibusTest(object, c("B_11", "E_107", "A_246"), correction = 1),
    "correction is not a string \\(a length one character vector\\)."
  )

  expect_error(
    aaPosOmnibusTest(object, c("B_11", "E_107", "A_246"), n_correction = 1.5),
    "n_correction is not a count \\(a single positive integer\\) or NULL."
  )

  expect_error(
    aaPosOmnibusTest(object, "FOO_2"),
    "amino acid position FOO_2 could not be found."
  )

  expect_error(
    aaPosOmnibusTest(object, c("B_11", "E_107", "A_246"), n_correction = 1),
    "n_correction must be at least 3."
  )
})
