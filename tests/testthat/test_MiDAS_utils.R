context("UTILS")

test_that("checkAlleleFormat", {
  expect_equal(checkAlleleFormat(c("A*01", "A*01:24", "B*01:25:22",
                                   "C*050:24:55:54")
               ),
               c(TRUE, TRUE, TRUE, TRUE)
  )
  expect_equal(checkAlleleFormat(c("01", "01:24", "01:25:22",
                                   "050:24:55:54")
               ),
               c(FALSE, FALSE, FALSE, FALSE)
  )
  expect_equal(checkAlleleFormat(c("*01", "A*:22", "C*05:24:55:54:89")
               ),
               c(FALSE, FALSE, FALSE)
  )
  expect_equal(checkAlleleFormat(c("*01", "A*:22", NA)
               ),
               c(FALSE, FALSE, NA)
  )
})

test_that("getAlleleResolution", {
  expect_equal(getAlleleResolution(c("A*01", "A*01:24", "B*01:25:22",
                                     "C*05:24:55:54", NA)
               ),
               c(2, 4, 6, 8, NA)
  )
  expect_error(getAlleleResolution("word"),
               "allele have to be a valid HLA allele number"
  )
})

test_that("reduceAlleleResolution", {
  expect_equal(reduceAlleleResolution(c("A*01", "A*01:24", "B*01:25:22",
                                        "C*05:24:55:54", "C*05:24:55:54N"), 2
               ),
               c("A*01", "A*01", "B*01", "C*05", "C*05:24:55:54N")
  )
  expect_error(getAlleleResolution("word"),
               "allele have to be a valid HLA allele number"
  )
  expect_error(reduceAlleleResolution("C*05:24:55:54", resolution = "four"),
               "resolution is not a count \\(a single positive integer\\)"
  )
})

test_that("getVariableAAPos", {
  allele <- c("TAP1*01:01:01:01", "TAP1*02:01:02")
  aln <- readHlaAlignments(system.file("extdata",
                                       "TAP1_prot.txt",
                                       package = "MiDAS"))
  aln <- aln[allele, ]
  expect_equal(getVariableAAPos(aln), c(`333` = 333, `637` = 637))

  expect_error(getVariableAAPos(c()), "alignment is not a matrix")
})

test_that("convertAlleleToVariable", {
  path <- system.file("extdata", "Match_allele_HLA_supertype.txt", package = "MiDAS")
  addvar <- convertAlleleToVariable(c("A*01:01", "A*02:01", "B*01", NA), dictionary = path)
  expect_equal(addvar, c("A01", "A02", NA, NA))
  dictionary <- read.table(path, header = TRUE, stringsAsFactors = FALSE)
  addvar <- convertAlleleToVariable(c("A*01:01", "A*02:01", "B*01", NA), dictionary = dictionary)
  expect_equal(addvar, c("A01", "A02", NA, NA))

  expect_error(
    convertAlleleToVariable(c("a", "b", "c"), dictionary = path),
    "allele have to be a valid HLA allele number"
  )

  expect_error(
    convertAlleleToVariable(c("A*01", "A*02", "A*03"), dictionary = c("foo", "bar")),
    "dictionary have to be either a path or a data.frame"
  )

  expect_error(
    convertAlleleToVariable(
      allele = c("A*01", "A*02", "A*03"),
      dictionary = file.path("foo", "bar")
    ),
    sprintf("Path '%s' does not exist", file.path("foo", "bar"))
  )

  expect_error(
    convertAlleleToVariable(c("A*01", "A*02", "A*03"), dictionary = dictionary[, 1, drop=FALSE]),
    "dictionary must have two columns"
  )

  expect_error(
    convertAlleleToVariable(c("A*01", "A*02", "A*03"), dictionary = dictionary[, c(2, 2)]),
    "first column in dictionary must contain valid HLA allele numbers"
  )

  expect_error(
    convertAlleleToVariable(c("A*01", "A*02", "A*03"), dictionary = dictionary[c(1, 1),]),
    "dictionary contains duplicated allele numbers"
  )
})

test_that("backquote", {
  expect_equal(backquote(c("A:01:01", "A:02:01")), c("`A:01:01`", "`A:02:01`"))
  expect_error(backquote(1), "x is not a character vector")
})

test_that("updateModel", {
  midas <-
    prepareMiDAS(
      hla_calls = MiDAS_tut_HLA,
      colData = MiDAS_tut_pheno,
      experiment = "hla_alleles"
    )
  midas_data <- midasToWide(midas, experiment = "hla_alleles")
  mod <- lm(disease ~ 1, data = midas_data)
  mod$call$data <- midas_data
  mod_test <- lm(disease ~ `A*01:01`, data = midas_data)
  mod_test$call$data <- midas_data
  expect_equal(updateModel(mod, "A*01:01"),
               mod_test
  )

  expect_error(updateModel(mod, 1),
               "x is not a character vector"
  )

  expect_error(updateModel(mod, x = "A*01:01", placeholder = 1),
               "placeholder is not a string \\(a length one character vector\\)."
  )

  expect_error(updateModel(mod, x = "A*01:01", backquote = 1),
               "backquote is not a flag \\(a length one logical vector\\)."
  )

  expect_error(updateModel(mod, x = "A*01:01", collapse = 1),
               "collapse is not a string \\(a length one character vector\\)."
  )

  expect_error(
    updateModel(
      mod,
      x = "A*01:01",
      placeholder = "foo"
    ),
    "placeholder 'foo' could not be found in object's formula"
  )
})

test_that("listMiDASDictionaries", {
  # NOTE ordering is somehow OS dependent
  expect_equal(
    sort(listMiDASDictionaries(pattern = ".*")),
    sort(c(
      "allele_HLA_Bw",
      "allele_HLA_Ggroup",
      "allele_HLA_supertype",
      "allele_HLA-Bw_only_B",
      "allele_HLA-C_C1-2",
      "counts_hla_kir_interactions",
      "counts_kir_haplotypes",
      "kir_haplotype_gene",
      "kir_nomenclature_gene"
    ))
  )
})

test_that("LRTest", {
  df <- data.frame(OS = c(20, 30, 40), AGE = c(50, 60, 70))
  mod0 <- lm(OS ~ 1, data = df)
  mod1 <- lm(OS ~ AGE, data = df)
  lrt_res <- LRTest(mod0, mod1)
  expect_equal(
    lrt_res,
    data.frame(
      term = "AGE",
      df = 1,
      logLik = 109.840111592134,
      statistic = 219.680223184268,
      p.value = 1.06202635429558e-49,
      stringsAsFactors = FALSE
    )
  )

  expect_error(LRTest(mod1, mod0), "variables AGE were not found in mod1")
})

test_that("getObjectDetails", {
  midas <- prepareMiDAS(
    kir_calls = MiDAS_tut_KIR,
    colData = MiDAS_tut_pheno,
    experiment = c("kir_genes"))
  obj <- lm(disease ~ term, data = midas)
  obj_det <- getObjectDetails(obj)
  expect_equal(
    obj_det,
    list(
      call = quote(lm(formula = disease ~ term, data = midas)),
      formula_vars = c("disease", "term"),
      data = midas
    )
  )
})

test_that("runMiDASGetVarsFreq", {
  midas <- prepareMiDAS(
    kir_calls = MiDAS_tut_KIR,
    colData = MiDAS_tut_pheno,
    experiment = c("kir_genes"))
  freq <- runMiDASGetVarsFreq(midas, "kir_genes", "disease")[1, ]
  freq_test <- data.frame(
    term = "KIR3DL3",
    Ntotal = 935,
    Ntotal.percent = formattable::percent(1),
    `N(disease=0)` = 467,
    `N(disease=0).percent` = formattable::percent(1),
    `N(disease=1)` = 468,
    `N(disease=1).percent` = formattable::percent(1),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  expect_equal(freq, freq_test)
})

test_that("distGrantham", {
  aa1 <- c("A", "S", "W")
  aa2 <- c("A", "S", "V")
  d <- distGrantham(aa1, aa2)
  d_test <- sum(dict_dist_grantham[paste0(aa1, aa2)]) / length(aa1)
  expect_equal(d, d_test)

  expect_error(distGrantham(1, aa2), "aa1 is not a character vector")

  expect_error(distGrantham(aa1, 1), "aa2 is not a character vector")

  expect_error(distGrantham(aa1, aa2[-3]),
               "aa1 and aa2 must have equal lengths.")

  expect_error(distGrantham(aa1, c("F", "O", "O")),
               "SO, WO are not valid amino acids pairs")
})

test_that("hlaCallsGranthamDistance", {
  gdist <- hlaCallsGranthamDistance(MiDAS_tut_HLA[1:5, ], genes = c("A", "B", "C"))
  gdist_test <- structure(list(
    ID = c("C001", "C002", "C003", "C004", "C005"),
    A = c(3.8121546961326, 0.87292817679558, 3.55801104972376, 10.6077348066298, 9.17127071823204),
    B = c(8.8121546961326, 6.64640883977901, 11.4088397790055, 9.41436464088398, 7.68508287292818),
    C = c(4.74585635359116, 2.90055248618785, 6.23204419889503, 7.05524861878453, 4.46961325966851)
    ),
    class = "data.frame",
    row.names = c(NA,-5L)
  )
  expect_equal(gdist, gdist_test)

  gdist <-
    hlaCallsGranthamDistance(MiDAS_tut_HLA[1:5,],
                             genes = c("A", "B", "C"),
                             aa_selection = "B_pocket")
  gdist_test <- structure(list(
    ID = c("C001", "C002", "C003", "C004", "C005"),
    A = c(0, 2, 0, 26.5454545454545, 7.72727272727273),
    B = c(30.7272727272727, 16.1818181818182, 48, 52.9090909090909, 27.3636363636364),
    C = c(36.6363636363636, 23.5454545454545, 37.1818181818182, 32.0909090909091, 32.0909090909091)
    ),
    class = "data.frame",
    row.names = c(NA,-5L)
  )
  expect_equal(gdist, gdist_test)

  # checkHlaCallsFormat test is ommitted here

  expect_error(hlaCallsGranthamDistance(MiDAS_tut_HLA, genes = 1),
               "genes is not a character vector")

  expect_error(hlaCallsGranthamDistance(MiDAS_tut_HLA, genes = c("A", NA)),
               "genes contains 1 missing values")

  hla_calls_bad <- MiDAS_tut_HLA
  hla_calls_bad[2, 2] <- "A*01"
  expect_error(
    hlaCallsGranthamDistance(hla_calls_bad, genes = "A"),
    "Allele resolutions for gene A are not equal"
  )

  expect_error(hlaCallsGranthamDistance(MiDAS_tut_HLA, genes = "A", aa_selection = "Z"),
               "aa_selection should be one of \"binding_groove\", \"B_pocket\", \"F_pocket\".")
})

test_that("hlaAlignmentGrantham", {
  aln <- hlaAlignmentGrantham("TAP1", 2, 2:182)
  aln_test <- readHlaAlignments(gene = "TAP1", resolution = 2)
  aln_test <- aln_test[, 2:182]
  mask <- apply(aln_test, 1, function(x) any(x == "" | x == "X" | x == "."))
  aln_test <- aln_test[! mask, ]
  expect_equal(aln, aln_test)
})

test_that("getHlaCallsGenes", {
  genes <- getHlaCallsGenes(MiDAS_tut_HLA)
  expect_equal(genes, c("A", "B", "C", "DPA1", "DPB1", "DQA1", "DQB1", "DRA", "DRB1"))
})

test_that("dfToExperimentMat", {
  mat <- dfToExperimentMat(MiDAS_tut_HLA)
  ids <- MiDAS_tut_HLA[["ID"]]
  test_mat <- MiDAS_tut_HLA[, -1]
  test_mat <- t(test_mat)
  colnames(test_mat) <- ids
  expect_equal(mat, test_mat)
})

test_that("experimentMatToDf", {
  ids <- MiDAS_tut_HLA[["ID"]]
  mat <- MiDAS_tut_HLA[, -1]
  mat <- t(mat)
  colnames(mat) <- ids
  expect_equal(experimentMatToDf(mat), MiDAS_tut_HLA)
})

test_that("midasToWide", {
  midas <- prepareMiDAS(
    hla_calls = MiDAS_tut_HLA[MiDAS_tut_HLA$ID %in% c("P001", "P002"), 1:2],
    colData = MiDAS_tut_pheno[MiDAS_tut_pheno$ID %in% c("P001", "P002"), ],
    experiment = "hla_alleles"
  )

  wide <- midasToWide(midas, "hla_alleles")
  test_wide <- data.frame(
    ID = c("P001", "P002"),
    `A*02:01` = 0:1,
    `A*11:88` = 1:0,
    disease = c(1L, 1L),
    lab_value = c(-0.85, -1.6),
    outcome = c(1L, 1L),
    term = colData(midas)$term,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  test_wide <- test_wide[, colnames(wide)] # somehow order of columns is R version dependent...
  expect_equal(wide, test_wide)

  expect_error(midasToWide(midas, 1), "experiment is not a character vector")

  expect_error(
    midasToWide(midas, "foo"),
    "experiment should match values \"hla_alleles\"."
  )
})

test_that("checkKirGenesFormat", {
  genes <- c("KIR3DL3", "KIR2DS4")
  expect_equal(checkKirGenesFormat(genes), c(TRUE, TRUE))
  expect_equal(checkKirGenesFormat(LETTERS), rep(FALSE, length(LETTERS)))
})

test_that("iterativeLRT", {
  MiDASdat <-
    prepareMiDAS(
      hla_calls = MiDAS_tut_HLA[, 1:3],
      colData = MiDAS_tut_pheno,
      experiment = "hla_aa"
    )
  MiDASdat <- filterByOmnibusGroups(MiDASdat, "hla_aa", c("A_29", "A_44", "A_65"))
  omnibus_groups <- getOmnibusGroups(MiDASdat, "hla_aa")
  placeholder <- getPlaceholder(MiDASdat)
  MiDASdat <- as.data.frame(MiDASdat)
  object <- lm(disease ~ outcome + term, data = MiDASdat)

  res <- iterativeLRT(object, placeholder, omnibus_groups)
  test_res <- lapply(
    X = omnibus_groups, 
    FUN = function (gr) {
      mod0 <- updateModel(
        object = object,
        x = "1",
        placeholder = placeholder,
        backquote = FALSE
      )
      test <- LRTest(
        mod0,
        updateModel(
          object = object,
          x = gr,
          placeholder = placeholder,
          collapse = " + ",
          backquote = TRUE
        )
      )
    })
  test_res <- dplyr::bind_rows(test_res, .id = "group")
  expect_equal(res, test_res)

  MiDASdat$A_29_D <- NA
  MiDASdat$A_29_A <- NA
  res <- iterativeLRT(object, placeholder, omnibus_groups)
  test_res <- lapply(
    X = omnibus_groups[c("A_44", "A_65")], 
    FUN = function (gr) {
      mod0 <- updateModel(
        object = object,
        x = "1",
        placeholder = placeholder,
        backquote = FALSE
      )
      test <- LRTest(
        mod0,
        updateModel(
          object = object,
          x = gr,
          placeholder = placeholder,
          collapse = " + ",
          backquote = TRUE
        )
      )
    })
  test_res[["A_29"]] <- data.frame(
    term = "A_29_D, A_29_A",
    df = NA,
    logLik = NA,
    statistic = NA,
    p.value = NA
  )
  test_res <- test_res[c("A_29", "A_44", "A_65")]
  test_res <- dplyr::bind_rows(test_res, .id = "group")
  expect_equal(res, test_res)
})

test_that("iterativeModel", {
  MiDASdat <-
    prepareMiDAS(
      hla_calls = MiDAS_tut_HLA[, 1:3],
      colData = MiDAS_tut_pheno,
      experiment = "hla_alleles"
    )
  placeholder <- getPlaceholder(MiDASdat)
  variables <- c("A*01:01", "A*01:02", "A*01:234")
  MiDASdat <- as.data.frame(MiDASdat)
  object <- lm(disease ~ outcome + term, data = MiDASdat)

  res <- iterativeModel(object, placeholder, variables)
  res_test <- lapply(
    X = variables,
    FUN = function(x) {
        obj <- updateModel(
          object = object,
          x = x,
          placeholder = placeholder,
          backquote = TRUE,
          collapse = " + "
        )
        r <- tidy(x = obj, conf.int = TRUE, exponentiate = exponentiate)
        r$term <- gsub("`", "", r$term)
        r <- r[r$term %in% variables, ]
    }
  )
  res_test <- dplyr::bind_rows(res_test)

  expect_equal(as.data.frame(res), as.data.frame(res_test))
})

test_that("getReferenceFrequencies", {
  freq <-
    getReferenceFrequencies(ref = allele_frequencies, pop = "USA NMDP Chinese")[1:3,]
  attr(freq, "reshapeWide") <- NULL
  freq_test <- data.frame(
    var = c("A*01:01", "A*01:03", "A*02:01"),
    `USA NMDP Chinese` = c(0.0145, 7e-05, 0.0946),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  expect_equal(freq, freq_test)

  freq <-
    getReferenceFrequencies(ref = allele_frequencies,
                            pop = "USA NMDP Chinese",
                            carrier_frequency = TRUE)[1:3,]
  attr(freq, "reshapeWide") <- NULL
  freq_test <- data.frame(
    var = c("A*01:01", "A*01:03", "A*02:01"),
    `USA NMDP Chinese` = c(0.02878975, 0.0001399951, 0.18025084),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  expect_equal(freq, freq_test)
})

test_that("adjustPValues", {
  p_val <- c(0.1, 0.001, 0.01)
  p_adj <- adjustPValues(p = p_val, method = "bonferroni", n = 3)
  expect_equal(p_adj, p_val * 3)

  p_adj <- adjustPValues(p = p_val, method = "bonferroni", n = 1)
  expect_equal(p_adj, p_val)

  expect_error(adjustPValues(p = p_val, method = "bonferroni", n = 0),
               "n must be >= 1")
})

test_that("filterListByElements", {
  A <- list(
    A = c("A", "B"),
    B = c("B", "C")
  )
  
  a <- filterListByElements(list = A, elements = c("B"))
  test_a <- list(
    A = c("B"),
    B = c("B")
  )
  expect_equal(a, test_a)
  
  a <- filterListByElements(list = A, elements = c("C"))
  test_a <- list(B = c("C"))
  expect_equal(a, test_a)
  
  a <- filterListByElements(list = A, elements = c())
  expect_equal(a, structure(list(), .Names = character(0)))
})
