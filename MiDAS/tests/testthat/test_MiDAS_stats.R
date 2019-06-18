context("HLA allele statistical methods")

test_that("HLA allele associations are analyzed properly", {
  hla_calls_file <-
    system.file("extdata", "HLAHD_output_example.txt", package = "MiDAS")
  hla_calls <- readHlaCalls(hla_calls_file)
  pheno_file <- system.file("extdata", "pheno_example.txt", package = "MiDAS")
  pheno <- read.table(pheno_file, header = TRUE, stringsAsFactors = FALSE)
  covar_file <- system.file("extdata", "covar_example.txt", package = "MiDAS")
  covar <- read.table(covar_file, header = TRUE, stringsAsFactors = FALSE)
  midas_data <<-
    prepareHlaData(hla_calls, pheno, covar, inheritance_model = "additive")

  object <- lm(OS_DIED ~ AGE + SEX, data = midas_data)
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
    analyzeAssociations(object, variables = "A*01:01", correction = 1),
    "correction is not a string \\(a length one character vector\\)."
  )

  expect_error(
    analyzeAssociations(object, variables = "A*01:01", exponentiate = 1),
    "exponentiate is not a flag \\(a length one logical vector\\)."
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
  midas_data <<-
    prepareHlaData(hla_calls, pheno, covar, inheritance_model = "additive")

  object <- coxph(Surv(OS, OS_DIED) ~ AGE + SEX, data = midas_data)
  res <- analyzeConditionalAssociations(object,
                                        variables = c("B*14:02", "DRB1*11:01"),
                                        th = 0.05)
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

  # Tests for checkStatisticalModel errors are ommitted here

  expect_error(
    analyzeConditionalAssociations(object, variables = 1),
    "variables is not a character vector"
  )

  expect_error(analyzeConditionalAssociations(object, variables = "thief"),
               "thief can not be found in object data"
  )

  expect_error(
    analyzeConditionalAssociations(object, variables =  "A*01:01", correction = 1),
    "correction is not a string \\(a length one character vector\\)."
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
})

test_that("HLA data is properly formatted", {
  small_hla_calls <- data.frame(ID = 1:2,
                                A_1 = c("A*01:01", "A*01:02"),
                                A_2 = c("A*01:02", "A*01:01"),
                                stringsAsFactors = FALSE
  )
  small_pheno <- data.frame(ID = 1:2, OS = c(123, 321), OS_DIED = c(0, 0))
  small_covar <- data.frame(ID = 1:2, AGE = c(23, 24))
  midas_data <- prepareHlaData(small_hla_calls, small_pheno, small_covar, inheritance_model = "additive")
  expect_equal(midas_data,
               data.frame(
                 ID = c(1, 2),
                 "A*01:01" = c(1, 1),
                 "A*01:02" = c(1, 1),
                 OS = c(123, 321),
                 OS_DIED = c(0, 0),
                 AGE = c(23, 24),
                 check.names = FALSE
               )
  )

  expect_error(
    prepareHlaData(hla_calls = "foo",
                   pheno = small_pheno,
                   covar = small_covar,
                   inheritance_model = "additive"
    ),
    "hla_calls is not a data frame"
  )

  expect_error(
    prepareHlaData(hla_calls = data.frame(),
                   pheno = small_pheno,
                   covar = small_covar,
                   inheritance_model = "additive"
    ),
    "hla_calls have to have at least 1 rows and 2 columns"
  )

  small_hla_calls_fac <- small_hla_calls
  small_hla_calls_fac$A_1 <- as.factor(small_hla_calls_fac$A_1)
  expect_error(
    prepareHlaData(hla_calls = small_hla_calls_fac,
                   pheno = small_pheno,
                   covar = small_covar,
                   inheritance_model = "additive"
    ),
    "hla_calls can't contain factors"
  )

  expect_error(
    prepareHlaData(hla_calls = small_hla_calls[, 2:3],
                   pheno = small_pheno,
                   covar = small_covar,
                   inheritance_model = "additive"
    ),
    "first column of hla_calls should specify samples id"
  )

  expect_error(
    prepareHlaData(hla_calls = small_hla_calls[, c(1, 1, 1)],
                   pheno = small_pheno,
                   covar = small_covar,
                   inheritance_model = "additive"
    ),
    "values in hla_calls doesn't follow HLA numbers specification"
  )

  expect_error(
    prepareHlaData(hla_calls = small_hla_calls,
                   pheno = "foo",
                   covar = small_covar,
                   inheritance_model = "additive"
    ),
    "additional_data_frame have to be a data frame"
  )

  bad_pheno <- small_pheno
  colnames(bad_pheno) <- LETTERS[1:ncol(bad_pheno)]
  expect_error(
    prepareHlaData(hla_calls = small_hla_calls,
                   pheno = bad_pheno,
                   covar = small_covar,
                   inheritance_model = "additive"
    ),
    "first column in additional_data_frame must be named as first column in hla_calls"
  )

  bad_pheno <- small_pheno
  bad_pheno$ID <- paste0("bad", bad_pheno$ID)
  expect_error(
    prepareHlaData(hla_calls = small_hla_calls,
                   pheno = bad_pheno,
                   covar = small_covar,
                   inheritance_model = "additive"
    ),
    "IDs in additional_data_frame doesn't match IDs in hla_calls"
  )

  expect_error(
    prepareHlaData(small_hla_calls,
                   pheno = small_pheno,
                   covar = small_covar,
                   inheritance_model = 1
    ),
    "inheritance_model is not a string \\(a length one character vector\\)."
  )

  expect_error(
    prepareHlaData(small_hla_calls,
                   pheno = small_pheno,
                   covar = small_covar,
                   inheritance_model = "foo"
    ),
    "inheritance_model should be one of 'dominant', 'recessive', 'additive'"
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
  midas_data <<-
    prepareHlaData(hla_calls, pheno, covar, inheritance_model = "additive")

  object <- lm(OS_DIED ~ AGE + SEX, data = midas_data)
  res <- analyzeMiDASData(object,
                          analysis_type = "hla_alleles",
                          variables = c("A*01:01", "A*02:01"),
                          kable_output = FALSE
  )

  test_res <- analyzeAssociations(object, variables = c("A*01:01", "A*02:01"))
  test_res <- dplyr::rename(test_res, allele = term)
  variables_freq <- getCountsFrequencies(midas_data[, c("ID", "A*01:01", "A*02:01")])
  test_res$Ntotal <- variables_freq$Counts
  test_res$Ntotal.frequency <- variables_freq$Freq
  pos <- midas_data$OS_DIED == 1
  pos_freq <- getCountsFrequencies(midas_data[pos, c("ID", "A*01:01", "A*02:01")])
  test_res$Npositive <- pos_freq$Counts
  test_res$Npositive.frequency <- pos_freq$Freq
  neg_freq <- getCountsFrequencies(midas_data[! pos, c("ID", "A*01:01", "A*02:01")])
  test_res$Nnegative <- neg_freq$Counts
  test_res$Nnegative.frequency <- neg_freq$Freq

  expect_equal(as.data.frame(res), as.data.frame(test_res)) # Tibble doesn't respect tollerance https://github.com/tidyverse/tibble/issues/287 or something related mby

  res <- analyzeMiDASData(object,
                          analysis_type = "hla_alleles",
                          conditional = TRUE,
                          kable_output = FALSE
  )

  test_res <-
    analyzeConditionalAssociations(object, variables = colnames(midas_data)[-1], th = 0.05)
  test_res <- dplyr::rename(test_res, allele = term)
  alleles <- test_res$allele
  variables_freq <- getCountsFrequencies(midas_data[, c("ID", alleles)])
  test_res$Ntotal <- variables_freq$Counts
  test_res$Ntotal.frequency <- variables_freq$Freq
  pos <- midas_data$OS_DIED == 1
  pos_freq <- getCountsFrequencies(midas_data[pos, c("ID", alleles)])
  test_res$Npositive <- pos_freq$Counts
  test_res$Npositive.frequency <- pos_freq$Freq
  neg_freq <- getCountsFrequencies(midas_data[! pos, c("ID", alleles)])
  test_res$Nnegative <- neg_freq$Counts
  test_res$Nnegative.frequency <- neg_freq$Freq

  expect_equal(as.data.frame(res), as.data.frame(test_res)) # Tibble doesn't respect tollerance https://github.com/tidyverse/tibble/issues/287 or something related mby

  # Tests for checkStatisticalModel errors are ommitted here

  expect_error(analyzeMiDASData(object, analysis_type = 1),
               "analysis_type is not a string \\(a length one character vector\\)."
  )

  expect_error(analyzeMiDASData(object, analysis_type = "a"),
               "analysis_type should be one of \"hla_alleles\", \"aa_level\", \"expression_levels\", \"allele_groups\", \"custom\"."
  )

  expect_error(analyzeMiDASData(object, analysis_type = "hla_alleles", conditional = 1),
               "conditional is not a flag \\(a length one logical vector\\)."
  )

  expect_error(analyzeMiDASData(object, analysis_type = "hla_alleles", variables = 1),
               "variables is not a character vector or NULL."
  )

  expect_error(analyzeMiDASData(object, analysis_type = "hla_alleles", variables = "thief"),
               "thief can not be found in object data"
  )

  expect_error(analyzeMiDASData(object, analysis_type = "hla_alleles", frequency_cutoff = "foo"),
               "frequency_cutoff is not number \\(a length one numeric vector\\) or NULL."
  )

  expect_error(analyzeMiDASData(object, analysis_type = "hla_alleles", pvalue_cutoff = "foo"),
               "pvalue_cutoff is not number \\(a length one numeric vector\\) or NULL."
  )

  expect_error(analyzeMiDASData(object, analysis_type = "hla_alleles", correction = NA),
               "correction is not a string \\(a length one character vector\\)."
  )

  expect_error(analyzeMiDASData(object, analysis_type = "hla_alleles", logistic = "NA"),
               "logistic is not a flag \\(a length one logical vector\\) or NULL."
  )

  expect_error(analyzeMiDASData(object, analysis_type = "hla_alleles", binary_phenotype = "NA"),
               "binary_phenotype is not a flag \\(a length one logical vector\\) or NULL."
  )

  expect_error(analyzeMiDASData(object, analysis_type = "hla_alleles", th = "NA"),
               "th is not a number \\(a length one numeric vector\\)."
  )

  expect_error(analyzeMiDASData(object, analysis_type = "hla_alleles", rss_th = "NA"),
               "rss_th is not a number \\(a length one numeric vector\\)."
  )

  expect_error(analyzeMiDASData(object, analysis_type = "hla_alleles", kable_output = "NA"),
               "kable_output is not a flag \\(a length one logical vector\\)."
  )

  expect_error(analyzeMiDASData(object, analysis_type = "hla_alleles", format = 1),
               "format is not a string \\(a length one character vector\\)."
  )

  expect_error(analyzeMiDASData(object, analysis_type = "hla_alleles", format = "pdf"),
               "format should be one of \"html\", \"latex\"."
  )
})

test_that("MiDAS data is prepared properly", {
  rleft_join <- function(init, ...) {
    Reduce(function(...)
      dplyr::left_join(..., by = "ID"),
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
    prepareMiDASData(hla_calls,
                     pheno,
                     covar,
                     analysis_type = "hla_allele",
                     inheritance_model = "additive")
  midas_hla_allele_test <-
    hlaCallsToCounts(hla_calls, inheritance_model = "additive")
  midas_hla_allele_test <- rleft_join(midas_hla_allele_test, pheno, covar)
  expect_equal(midas_hla_allele, midas_hla_allele_test)

  # aa_level
  midas_aa_level <- prepareMiDASData(hla_calls,
                                     pheno,
                                     covar,
                                     analysis_type = "aa_level",
                                     inheritance_model = "additive")
  midas_aa_level_test <- hlaToAAVariation(hla_calls)
  midas_aa_level_test <-
    aaVariationToCounts(midas_aa_level_test, inheritance_model = "additive")
  midas_aa_level_test <- rleft_join(midas_aa_level_test, pheno, covar)
  expect_equal(midas_aa_level, midas_aa_level_test)

  # expression_levels
  midas_expression_levels <- prepareMiDASData(hla_calls,
                                              pheno,
                                              covar,
                                              analysis_type = "expression_levels",
                                              inheritance_model = "additive")
  expression_dicts <- grep("expression", listMiDASDictionaries(), value = TRUE)
  midas_expression_levels_test <- Reduce(
    f = function(...) dplyr::left_join(..., by = "ID"),
    x = lapply(expression_dicts, function(x) {
      expr <- hlaToVariable(hla_calls = hla_calls, dictionary = x)
      expr$sum <- rowSums(expr[, -1, drop = FALSE])
      gene <- gsub("_1", "", colnames(expr)[2])
      expr <- expr[, c("ID", "sum")]
      colnames(expr) <- c("ID", paste0("expression_", gene))
      expr
    })
  )
  midas_expression_levels_test <-
    rleft_join(midas_expression_levels_test, pheno, covar)
  expect_equal(midas_expression_levels, midas_expression_levels_test)

  # allele_groups
  midas_allele_groups <- prepareMiDASData(hla_calls,
                                          pheno,
                                          covar,
                                          analysis_type = "allele_groups",
                                          inheritance_model = "additive")
  groups_dicts <-
    grep("expression",
         listMiDASDictionaries(),
         value = TRUE,
         invert = TRUE)
  midas_allele_groups_test <- Reduce(
    f = function(...) dplyr::left_join(..., by = "ID"),
    x = lapply(groups_dicts, hlaToVariable, hla_calls = hla_calls)
  )
  midas_allele_groups_test <- rleft_join(midas_allele_groups_test, pheno, covar)
  expect_equal(midas_allele_groups, midas_allele_groups_test)

  # custom
  midas_custom <- prepareMiDASData(hla_calls,
                                   pheno,
                                   covar,
                                   analysis_type = "custom",
                                   inheritance_model = "additive")
  midas_custom_test <- rleft_join(hla_calls, pheno, covar)
  expect_equal(midas_custom, midas_custom_test)

  # test for checkHlaCallsFormat are ommitted here

  expect_error(
    prepareMiDASData(hla_calls, analysis_type = 1),
    "analysis_type is not a string \\(a length one character vector\\)."
  )

  expect_error(
    prepareMiDASData(hla_calls, analysis_type = "foo"),
    "analysis_type should be one of \"hla_allele\", \"aa_level\", \"expression_levels\", \"allele_groups\", \"custom\"."
  )

  expect_error(
    prepareMiDASData(hla_calls, analysis_type = "hla_allele", inheritance_model = 1),
    "inheritance_model is not a string \\(a length one character vector\\)."
  )

  expect_error(
    prepareMiDASData(hla_calls, analysis_type = "hla_allele", inheritance_model = "bar"),
    "inheritance_model should be one of \"dominant\", \"recessive\", \"additive\"."
  )

  expect_error(
    prepareMiDASData(hla_calls, analysis_type = "hla_allele", indels = "no"),
    "indels is not a flag \\(a length one logical vector\\)."
  )

  expect_error(
    prepareMiDASData(hla_calls, analysis_type = "hla_allele", unkchar = "nope"),
    "unkchar is not a flag \\(a length one logical vector\\)."
  )

  # checkAdditionalData on ... argument are ommitted here

  expect_error(
    prepareMiDASData(hla_calls[, c("ID", "DMA_1", "DMA_2")], analysis_type = "expression_levels"),
    "no expression levels were found for input hla_calls"
  )

  # this is ill due to ggroups matches problem there is one more above, when groups are fixed this will start to fail so uncomment and remove linies as needed
  expect_error(
    prepareMiDASData(hla_calls[, c("ID", "DMA_1", "DMA_2")], analysis_type = "allele_groups"),
   # "no expression levels were found for input hla_calls"
    "match table contains duplicated allele numbers"
  )
})
