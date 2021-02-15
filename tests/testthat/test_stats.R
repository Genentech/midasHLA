context("stats")

test_that("analyzeAssociations", {
  midas_data <- midasToWide(MiDAS_tut_object, experiment = "hla_alleles")
  object <- lm(disease ~ term, data = midas_data)

  res <- analyzeAssociations(object,
                             variables = c("A*01:01", "A*02:01"),
                             correction = "BH"
  )

  test_res <- list(
    lm(disease ~ `A*01:01`, data = midas_data),
    lm(disease ~ `A*02:01`, data = midas_data)
  )
  test_res <- do.call("rbind", lapply(test_res, tidy, conf.int = TRUE))
  test_res$term <- gsub("`", "", test_res$term)
  test_res <- test_res[test_res$term %in% c("A*01:01", "A*02:01"), ]
  test_res$p.adjusted <- p.adjust(test_res$p.value, "BH")

  expect_equal(as.data.frame(res), as.data.frame(test_res)) # Tibble doesn't respect tollerance https://github.com/tidyverse/tibble/issues/287 or something related mby

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
    analyzeAssociations(object, variables = "A*01:01", placeholder = "foo"),
    "placeholder 'foo' could not be found in object's formula"
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
      correction = "BH",
      n_correction = 1
    ),
    "n_correction must be at least 2."
  )
})

test_that("analyzeConditionalAssociations", {
  midas_data <- midasToWide(MiDAS_tut_object, experiment = "hla_alleles")

  object <- lm(disease ~ term, data = midas_data)

  # keep = FALSE
  res <- analyzeConditionalAssociations(object,
                                        variables = c("DQB1*06:02", "B*57:01"),
                                        th = 0.05,
                                        keep = FALSE)
  res <- rapply(res, classes = "numeric", how = "replace", round, digits = 3)

  test_res <- tibble(term = c("DQB1*06:02", "B*57:01"),
                     estimate = c(0.15, 0.235),
                     std.error = c(0.035, 0.056),
                     statistic = c(4.319, 4.228),
                     p.value = c(0, 0),
                     conf.low = c(0.082, 0.126),
                     conf.high = c(0.219, 0.344),
                     p.adjusted = c(0, 0),
                     covariates = c("", "DQB1*06:02")
  )

  expect_equal(res, test_res)

  # keep = TRUE
  res <- analyzeConditionalAssociations(object,
                                        variables = c("DQB1*06:02", "B*57:01"),
                                        th = 0.05,
                                        keep = TRUE)
  res <- rapply(res, classes = "numeric", how = "replace", round, digits = 3)
  test_res <- list(
    tibble(term = c("DQB1*06:02", "B*57:01"),
           estimate = c(0.15, 0.236),
           std.error = c(0.035, 0.056),
           statistic = c(4.319, 4.21),
           p.value =  c(0, 0),
           conf.low = c(0.082, 0.126),
           conf.high = c(0.219, 0.346),
           p.adjusted = c(0, 0),
           covariates = c("", "")
    ),
    tibble(term = "B*57:01",
           estimate = 0.235,
           std.error = 0.056,
           statistic = 4.228,
           p.value = 0,
           conf.low = 0.126,
           conf.high = 0.344,
           p.adjusted = 0,
           covariates = "DQB1*06:02"
    )
  )
  expect_equal(res, test_res)

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
    analyzeConditionalAssociations(object, variables = "A*01:01", placeholder = "foo"),
    "placeholder 'foo' could not be found in object's formula"
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
      correction = "BH",
      n_correction = 1
    ),
    "n_correction must be at least 2."
  )
})

test_that("omnibusTest", {
  midas_data <- midasToWide(MiDAS_tut_object, experiment = "hla_aa")
  object <- lm(disease ~ term, data = midas_data)
  omnibus_groups <- list(
    A_77 = c("A_77_D", "A_77_N", "A_77_S"),
    A_79 = c("A_79_G", "A_79_R")
  )
  omnibus_res <- omnibusTest(object, omnibus_groups)

  obj_A77 <- lm(disease ~ A_77_D + A_77_N + A_77_S, data = midas_data)
  obj_A79 <- lm(disease ~ A_79_G + A_79_R , data = midas_data)
  LRT <-
    lapply(list(obj_A77, obj_A79),
           LRTest,
           mod0 = lm(disease ~ 1, data = midas_data))
  omnibus_res_test <- data.frame(
    group = c("A_77", "A_79"),
    term = c("A_77_D, A_77_N", "A_79_G"),
    df = sapply(LRT, `[[`, "df"),
    logLik = sapply(LRT, `[[`, "logLik"),
    statistic = sapply(LRT, `[[`, "statistic"),
    p.value = sapply(LRT, `[[`, "p.value"),
    p.adjusted = p.adjust(sapply(LRT, `[[`, "p.value"), method = "bonferroni"),
    stringsAsFactors = FALSE
  )
  expect_equal(omnibus_res, omnibus_res_test)
})

test_that("HWETest", {
  midas <- MiDAS_tut_object
  res <- HWETest(midas, experiment = "hla_alleles", HWE_cutoff = 0.988)
  rownames(res) <- NULL
  test_res <- data.frame(
    var = c("A*02:01", "DRB1*13:02"),
    p.value = c(0.993279078042026, 0.988280722372832),
    stringsAsFactors = FALSE
  )
  expect_equal(as.data.frame(res), as.data.frame(test_res))

  res <-
    HWETest(
      midas,
      experiment = "hla_alleles",
      HWE_cutoff = 0.988,
      as.MiDAS = TRUE
    )
  test_res <- filterByVariables(midas, "hla_alleles", c("A*02:01", "DRB1*13:02"))
  expect_equal(res, test_res)

  expect_error(HWETest("foo"), "object must be an instance of \"MiDAS\".")

  expect_error(HWETest(midas, experiment = 1), "experiment is not a string \\(a length one character vector\\).")

  expect_error(HWETest(midas, experiment = "foo"), "experiment should be one of \"hla_alleles\", \"hla_aa\", \"hla_g_groups\", \"hla_supertypes\", \"hla_NK_ligands\".")

  expect_error(HWETest(midas, experiment = "hla_alleles", HWE_cutoff = "foo"), "HWE_cutoff is not a number \\(a length one numeric vector\\) or NULL.")

  expect_error(HWETest(midas, experiment = "hla_alleles", as.MiDAS = "foo"), "as.MiDAS is not a flag \\(a length one logical vector\\).")
})
