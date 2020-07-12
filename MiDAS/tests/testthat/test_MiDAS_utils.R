context("HLA allele numbers")

test_that("HLA allele numbers have proper format", {
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
})

test_that("HLA allele resolution is number of sets of digits * 2", {
  expect_equal(getAlleleResolution(c("A*01", "A*01:24", "B*01:25:22",
                                     "C*05:24:55:54")
               ),
               c(2, 4, 6, 8)
  )
  expect_error(getAlleleResolution("word"),
               "allele have to be a valid HLA allele number"
  )
})

test_that("Reduced HLA allele have desired resoulution", {
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

test_that("HLA allels are converted to additional variables", {
  path <- system.file("extdata", "Match_allele_HLA_supertype.txt", package = "MiDAS")
  addvar <- convertAlleleToVariable(c("A*01:01", "A*02:01", "B*01", NA), dictionary = path)
  expect_equal(addvar, c("A01", "A02", NA, NA))
  dictionary <- read.table(path, header = TRUE, stringsAsFactors = FALSE)
  addvar <- convertAlleleToVariable(c("A*01:01", "A*02:01", "B*01", NA), dictionary = dictionary)
  expect_equal(addvar, c("A01", "A02", NA, NA))

  expect_error(convertAlleleToVariable(c("a", "b", "c"), dictionary = path,
                                       "allele have to be a valid HLA allele number"
  )
  )

  expect_error(convertAlleleToVariable(c("A*01", "A*02", "A*03"), dictionary = c("foo", "bar"),
                                       "dictionary have to be either path or data.frame"
  )
  )

  expect_error(
    convertAlleleToVariable(
      allele = c("A*01", "A*02", "A*03"),
      dictionary = file.path("foo", "bar")
    ),
    sprintf("Path '%s' does not exist", file.path("foo", "bar"))
  )

  expect_error(convertAlleleToVariable(c("A*01", "A*02", "A*03"), dictionary = dictionary[, 1],
                                       "match table have to consist out of two columns"
  )
  )

  expect_error(convertAlleleToVariable(c("A*01", "A*02", "A*03"), dictionary = dictionary[, c(2, 2)],
                                       "first column of match table must contain valid HLA allele numbers"
  )
  )

  expect_error(convertAlleleToVariable(c("A*01", "A*02", "A*03"), dictionary = dictionary[c(1, 1), ],
                                       "match table contains duplicated allele numbers"
  )
  )
})

test_that("HLA calls data frame have proper format", {
  file <- system.file("extdata", "HLAHD_output_example.txt", package = "MiDAS")
  hla_calls <- readHlaCalls(file)
  expect_equal(checkHlaCallsFormat(hla_calls), TRUE)

  expect_error(checkHlaCallsFormat("A"), "hla_calls is not a data frame")
  expect_error(checkHlaCallsFormat(data.frame()),
               "hla_calls have to have at least 1 rows and 2 columns"
  )
  hla_calls[, 1] <- as.factor(hla_calls[, 1])
  expect_error(checkHlaCallsFormat(hla_calls),
               "hla_calls can't contain factors"
  )
  fake_calls <- data.frame(ID = c("Sample1", "Sample2", "Sample3"),
                           A_1 = c("A*01", "A*02", "A*03"),
                           A_2 = c("A*01", "B*02", "C*03"),
                           stringsAsFactors = FALSE
  )
  expect_error(checkHlaCallsFormat(fake_calls[, c(2, 1, 3)]),
               "first column of hla_calls should specify samples id"
  )

  expect_error(checkHlaCallsFormat(fake_calls[, c(1, 1, 3)]),
               "values: Sample1, Sample2, Sample3 in hla_calls doesn't follow HLA numbers specification"
  )
})

test_that("HLA allele are backquoted properly", {
  expect_equal(backquote(c("A:01:01", "A:02:01")), c("`A:01:01`", "`A:02:01`"))
  expect_error(backquote(1), "x is not a character vector")
})

context("HLA allele alignments")

test_that("Variable amino acids positions are detected properly", {
  hlaa_calls <- c("TAP1*01:01", "TAP1*02:01")
  hlaa_res <- 4
  hlaa_aln <- readHlaAlignments(system.file("extdata",
                                            "TAP1_prot.txt",
                                            package = "MiDAS")
  )
  four_dig_numbers <- reduceAlleleResolution(rownames(hlaa_aln), resolution = 4)
  hlaa_aln <- hlaa_aln[!duplicated(four_dig_numbers), ]
  rownames(hlaa_aln) <- four_dig_numbers[!duplicated(four_dig_numbers)]
  hlaa_aln <- hlaa_aln[hlaa_calls, ]

  expect_equal(getVariableAAPos(hlaa_aln), c(`333` = 333, `637` = 637))

  expect_error(getVariableAAPos(hlaa_calls), "alignment is not a matrix")
})

context("HLA statistical models handling")

test_that("HLA statistical models are updated properly", {
  library("survival")
  hla_calls_file <- system.file("extdata",
                                "HLAHD_output_example.txt",
                                package = "MiDAS"
  )
  hla_calls <- readHlaCalls(hla_calls_file)
  pheno_file <- system.file("extdata", "pheno_example.txt", package = "MiDAS")
  pheno <- read.table(pheno_file, header = TRUE, stringsAsFactors = FALSE)
  covar_file <- system.file("extdata", "covar_example.txt", package = "MiDAS")
  covar <- read.table(covar_file, header = TRUE, stringsAsFactors = FALSE)
  coldata <- dplyr::left_join(pheno, covar, by = "ID")
  midas <-
    prepareMiDAS(
      hla_calls,
      colData = coldata,
      experiment = "hla_alleles",
      inheritance_model = "additive"
    )
  midas_data <- midasToWide(midas, experiment = "hla_alleles")
  coxmod <- coxph(Surv(OS, OS_DIED) ~ 1, data = midas_data)
  coxmod$call$data <- midas_data
  coxmod_test <- coxph(Surv(OS, OS_DIED) ~ `A*01:01`, data = midas_data)
  coxmod_test$call$data <- midas_data
  expect_equal(updateModel(coxmod, "A*01:01"),
               coxmod_test
  )

  expect_error(updateModel(coxmod, 1),
               "x is not a character vector"
  )

  expect_error(updateModel(coxmod, x = "A*01:01", placeholder = 1),
               "placeholder is not a string \\(a length one character vector\\)."
  )

  expect_error(updateModel(coxmod, x = "A*01:01", backquote = 1),
               "backquote is not a flag \\(a length one logical vector\\)."
  )

  expect_error(updateModel(coxmod, x = "A*01:01", collapse = 1),
               "collapse is not a string \\(a length one character vector\\)."
  )

  expect_error(
    updateModel(
      coxmod,
      x = "A*01:01",
      placeholder = "foo"
    ),
    "placeholder 'foo' could not be found in object's formula"
  )
})


test_that("statistical models are statistical model", {
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

  object <- lm(OS ~ OS_DIED, data = midas)
  expect_equal(checkStatisticalModel(object), TRUE)

  expect_error(checkStatisticalModel(list(1)),
               "object is required to have the internal OBJECT bit set"
  )

  expect_error(checkStatisticalModel(speed ~ cars),
               "object have to have an attribute 'call'"
  )

  fake_model <- list(call = list(formula = "foo"))
  class(fake_model) <- "fake"
  expect_error(checkStatisticalModel(fake_model),
               "object have to be a model with defined formula"
  )

  fake_model <- list(call = list(formula = 1 ~ 1))
  class(fake_model) <- "fake"
  expect_error(checkStatisticalModel(fake_model),
               "object need to have data attribute defined"
  )
})

test_that("is counts or zeros", {
  expect_equal(isCountsOrZeros(c(1, 0, 2, NA)), TRUE)

  expect_error(
    assertthat::assert_that(isCountsOrZeros(c(1, 0, 2, NA, 1.5))),
    "values in c\\(1, 0, 2, NA, 1.5\\) are not counts \\(a positive integers\\) or zeros."
  )
})

test_that("is character or null", {
  expect_equal(isCharacterOrNULL(LETTERS), TRUE)
  expect_equal(isCharacterOrNULL(NULL), TRUE)

  expect_error(
    assertthat::assert_that(isCharacterOrNULL(1)),
    "1 is not a character vector or NULL."
  )
})

test_that("is number or null", {
  expect_equal(isNumberOrNULL(1), TRUE)
  expect_equal(isNumberOrNULL(NULL), TRUE)

  expect_error(
    assertthat::assert_that(isNumberOrNULL("a")),
    "\"a\" is not a number \\(a length one numeric vector\\) or NULL."
  )
})

test_that("is string or null", {
  expect_equal(isStringOrNULL("foo"), TRUE)
  expect_equal(isStringOrNULL(NULL), TRUE)

  expect_error(
    assertthat::assert_that(isStringOrNULL(1)),
    "1 is not a string \\(a length one character vector\\) or NULL."
  )
})

test_that("string matches", {
  expect_equal(stringMatches("foo", c("foo", "bar")), TRUE)

  expect_error(
    assertthat::assert_that(stringMatches("foo", c("bar", "Foo"))),
    "\"foo\" should be one of \"bar\", \"Foo\"."
  )
})

test_that("is flag or null", {
  expect_equal(isFlagOrNULL(TRUE), TRUE)
  expect_equal(isFlagOrNULL(NULL), TRUE)
  expect_equal(isFlagOrNULL(NA), FALSE)

  expect_error(
    assertthat::assert_that(isFlagOrNULL(1)),
    "1 is not a flag \\(a length one logical vector\\) or NULL."
  )
})

test_that("character maches choices", {
  expect_equal(characterMatches("foo", c("foo", "bar")), TRUE)

  expect_error(
    assertthat::assert_that(characterMatches("foo", "bar")),
    '"foo" should match values "bar".'
  )
})

test_that("is class or null", {
  expect_equal(isClassOrNULL("foo", "character"), TRUE)
  expect_equal(isClassOrNULL(NULL, "character"), TRUE)

  expect_error(
    assertthat::assert_that(isClassOrNULL("foo", "bar")),
    "\"foo\" must be an instance of \"bar\" or NULL."
  )
})

test_that("KIR haplotypes are converted to gene counts", {
  x <- c("cA01~tA01+cB02~tA01", "cA01~tA01+cA01~tB01_2DS5")
  kir_hap <- kirHaplotypeToCounts(x)
  test_kir_hap <- data.frame(
    haplotypes = c("cA01~tA01+cB02~tA01", "cA01~tA01+cA01~tB01_2DS5"),
    KIR3DL3 = c(1, 1),
    KIR2DS2 = c(1, 0),
    KIR2DL2 = c(1, 0),
    KIR2DL3 = c(1, 1),
    KIR2DP1 = c(1, 1),
    KIR2DL1 = c(1, 1),
    KIR3DP1 = c(1, 1),
    KIR2DL4 = c(1, 1),
    KIR3DL1 = c(1, 1),
    KIR3DS1 = c(0, 1),
    KIR2DL5 = c(0, 1),
    KIR2DS3 = c(0, 0),
    KIR2DS5 = c(0, 1),
    KIR2DS4 = c(1, 1),
    KIR2DS1 = c(0, 1),
    KIR3DL2 = c(1, 1),
    stringsAsFactors = FALSE
  )
  expect_equal(kir_hap, test_kir_hap)
})

test_that("colnamesMatches", {
  df <- data.frame(a = 1:5, b = 1:5)
  expect_equal(colnamesMatches(df, c("a", "b")), TRUE)

  expect_error(colnamesMatches(1:2, c("foo", "bar")), "x is not a data frame")

  expect_error(colnamesMatches(data.frame(one = 1:2), c("foo", "bar")),
               "Number of columns in data.frame\\(one = 1:2\\) must equal 2."
  )

  expect_error(
    assertthat::assert_that(colnamesMatches(df, c("foo", "bar"))),
    "Columns: 'a', 'b' in df should be named 'foo', 'bar'"
  )
})

test_that("KIR counts have proper format", {
  file <- system.file("extdata", "KPI_output_example.txt", package = "MiDAS")
  kir_counts <- readKPICalls(file)
  expect_equal(checkKirCallsFormat(kir_counts), TRUE)

  expect_equal(checkKirCallsFormat(NULL, accept.null = TRUE), TRUE)

  fake_kir_counts <- kir_counts
  fake_kir_counts[, 1] <- as.factor(fake_kir_counts[, 1, drop = TRUE])
  expect_error(
    checkKirCallsFormat(fake_kir_counts),
    "fake_kir_counts can't contain factors"
  )


  expect_error(
    checkKirCallsFormat(kir_counts[, 1, drop = FALSE]),
    "Number of columns in kir_calls must equal 17."
  )

  fake_kir_counts <- kir_counts
  colnames(fake_kir_counts) <- c("FOO", colnames(fake_kir_counts)[-1])
  expect_error(
    checkKirCallsFormat(fake_kir_counts),
    "Columns: 'FOO' in kir_calls should be named 'ID'"
  )
})

test_that("is count or null", {
  expect_equal(isCountOrNULL(1), TRUE)
  expect_equal(isCountOrNULL(NULL), TRUE)

  expect_error(
    assertthat::assert_that(isCountOrNULL(1.5)),
    "1.5 is not a count \\(a single positive integer\\) or NULL."
  )
})

test_that("is true or false", {
  expect_equal(isTRUEorFALSE(TRUE), TRUE)
  expect_equal(isTRUEorFALSE(FALSE), TRUE)
  expect_equal(isTRUEorFALSE(NA), FALSE)

  expect_error(
    assertthat::assert_that(isTRUEorFALSE(1.5)),
    "1.5 is not a flag \\(a length one logical vector\\)."
  )
})

test_that("tidy method exists", {
  expect_equal(hasTidyMethod("lm"), TRUE)
  expect_equal(hasTidyMethod("foo"), FALSE)

  expect_error(
    assertthat::assert_that(hasTidyMethod("bar")),
    "tidy function for object of class \"bar\" could not be found."
  )
})

test_that("likelihood ratio test", {
  df <- data.frame(OS = c(20, 30, 40), AGE = c(50, 60, 70))
  mod0 <- lm(OS ~ 1, data = df)
  mod1 <- lm(OS ~ AGE, data = df)
  lrt_res <- LRTest(mod0, mod1)
  expect_equal(
    lrt_res,
    data.frame(
      term = "AGE",
      dof = 1,
      logLik = 109.8401115921340078785,
      statistic = 219.680223184268015757,
      p.value = 1.062026e-49,
      stringsAsFactors = FALSE
    )
  )

  expect_error(LRTest(mod1, mod0), "variables AGE were not found in mod1")
})

test_that("object has placeholder", {
  object <- lm(speed ~ dist, data = cars)
  expect_equal(objectHasPlaceholder(object, "dist"), TRUE)
  expect_equal(objectHasPlaceholder(object, "foo"), FALSE)

  expect_error(
    assertthat::assert_that(objectHasPlaceholder(object, "foo")),
    "placeholder 'foo' could not be found in object's formula"
  )
})

test_that("phenotype data is properly formatted", {
  pheno <- data.frame(
    ID = 1:5,
    letter = LETTERS[1:5]
  )

  expect_equal(checkColDataFormat(pheno), TRUE)

  expect_error(
    checkColDataFormat(LETTERS),
    "LETTERS have to be a data frame"
  )

  expect_error(
    checkColDataFormat(data.frame()),
    "data.frame\\(\\) have to have at least 1 row and 2 columns"
  )

  expect_error(
    checkColDataFormat(pheno[, 2, drop = FALSE]),
    "pheno\\[, 2, drop = FALSE\\] have to have at least 1 row and 2 columns"
  )
})

test_that("check if function exists", {
  expect_equal(functionExists("lm"), TRUE)

  expect_error(
    assertthat::assert_that(functionExists("foo")),
    "Function foo could not be found."
  )
})

test_that("is class", {
  expect_equal(isClass("foo", "character"), TRUE)

  expect_error(
    assertthat::assert_that(isClassOrNULL("foo", "bar")),
    "\"foo\" must be an instance of \"bar\"."
  )
})

test_that("Grantham distance is calculated properly", {
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

test_that("Between allele Grantham distance is calculated properly", {
  file <- system.file("extdata", "HLAHD_output_example.txt", package = "MiDAS")
  hla_calls <- readHlaCalls(file)[1:5, ]
  gdist <- hlaCallsGranthamDistance(hla_calls, genes = c("A", "B", "C"))
  gdist_test <- structure(list(
    ID = c("PAT1", "PAT2", "PAT3", "PAT4", "PAT5"),
    A = c(0, 0, 0, 0.121546961325967, 0),
    B = c(7.20441988950276, 13.9337016574586, 0, 10.2265193370166, 0),
    C = c(7.27624309392265, 6.58563535911602, 0, 0, 0)
  ),
  class = "data.frame",
  row.names = c(NA,-5L))
  expect_equal(gdist, gdist_test)

  # checkHlaCallsFormat test is ommitted here

  expect_error(hlaCallsGranthamDistance(hla_calls, genes = 1),
               "genes is not a character vector")

  expect_error(hlaCallsGranthamDistance(hla_calls, genes = c("A", NA)),
               "genes contains 1 missing values")

  hla_calls_bad <- hla_calls
  hla_calls_bad[2, 2] <- "A*01"
  expect_error(
    hlaCallsGranthamDistance(hla_calls_bad, genes = "A"),
    "Allele resolutions for gene A are not equal"
  )
})

test_that("Frequency cutoffs validation", {
  # lower_frequency_cutof must be a number
  lower_frequency_cutoff <- "foo"
  upper_frequency_cutoff <- 0.5
  expect_error(
    validateFrequencyCutoffs(lower_frequency_cutoff, upper_frequency_cutoff),
    "lower_frequency_cutoff is not a number \\(a length one numeric vector\\)."
  )

  # lower_frequency_cutof must be positive
  lower_frequency_cutoff <- -1
  upper_frequency_cutoff <- 0.5
  expect_error(
    validateFrequencyCutoffs(lower_frequency_cutoff, upper_frequency_cutoff),
    "lower_frequency_cutoff must be a number greater than 0."
  )

  # upper_frequency_cutoff must be a number
  lower_frequency_cutoff <- 0.5
  upper_frequency_cutoff <- "foo"
  expect_error(
    validateFrequencyCutoffs(lower_frequency_cutoff, upper_frequency_cutoff),
    "upper_frequency_cutoff is not a number \\(a length one numeric vector\\)."
  )

  # upper_frequency_cutoff must be positive
  lower_frequency_cutoff <- 0
  upper_frequency_cutoff <- -1
  expect_error(
    validateFrequencyCutoffs(lower_frequency_cutoff, upper_frequency_cutoff),
    "upper_frequency_cutoff must be a number greater than 0."
  )

  # lower_frequency_cutoff is lower than upper_frequency_cutoff
  lower_frequency_cutoff <- 5
  upper_frequency_cutoff <- 1
  expect_error(
    validateFrequencyCutoffs(lower_frequency_cutoff, upper_frequency_cutoff),
    "lower_frequency_cutoff cannot be higher than upper_frequency_cutoff."
  )

  # Both lower_frequency_cutoff and upper_frequency_cutoff have to be either frequencies or counts
  lower_frequency_cutoff <- 0.5
  upper_frequency_cutoff <- 2
  expect_error(
    validateFrequencyCutoffs(lower_frequency_cutoff, upper_frequency_cutoff),
    "Both lower_frequency_cutoff and upper_frequency_cutoff have to be either frequencies or counts."
  )
})

test_that("getHlaCallsGenes", {
  file <- system.file("extdata", "HLAHD_output_example.txt", package = "MiDAS")
  hla_calls <- readHlaCalls(file)[, 1:5]
  genes <- getHlaCallsGenes(hla_calls)
  expect_equal(genes, c("A", "B"))
})

test_that("dfToExperimentMat", {
  file <- system.file("extdata", "HLAHD_output_example.txt", package = "MiDAS")
  hla_calls <- readHlaCalls(file)
  mat <- dfToExperimentMat(hla_calls)
  ids <- hla_calls[["ID"]]
  test_mat <- hla_calls[, -1]
  test_mat <- t(test_mat)
  colnames(test_mat) <- ids
  expect_equal(mat, test_mat)
})

test_that("experimentMatToDf", {
  file <- system.file("extdata", "HLAHD_output_example.txt", package = "MiDAS")
  hla_calls <- readHlaCalls(file)
  ids <- hla_calls[["ID"]]
  mat <- hla_calls[, -1]
  mat <- t(mat)
  colnames(mat) <- ids
  expect_equal(experimentMatToDf(mat), hla_calls)
})

test_that("midasToWide", {
  hla_calls_file <- system.file("extdata", "HLAHD_output_example.txt", package = "MiDAS")
  hla_calls <- readHlaCalls(hla_calls_file)[1:5, 1:5]

  pheno_file <- system.file("extdata", "pheno_example.txt", package = "MiDAS")
  pheno <- read.table(pheno_file, header = TRUE, stringsAsFactors = FALSE)[1:5, ]

  midas <- prepareMiDAS(
    hla_calls = hla_calls,
    colData = pheno,
    inheritance_model = "additive",
    experiment = "hla_alleles"
  )

  wide <- midasToWide(midas, "hla_alleles")
  test_wide <- data.frame(
    primary = c("PAT1", "PAT2", "PAT3", "PAT4", "PAT5"),
    ID = c("PAT1", "PAT2", "PAT3", "PAT4", "PAT5"),
    OS = c(280L, 458L, 415L, 211L, 631L),
    OS_DIED = c(1L, 0L, 0L, 1L, 0L),
    term = wide$term, # there is a rounding error
    `A*01:01` = c(0L, 0L, 2L, 0L, 0L),
    `A*02:01` = c(2L, 2L, 0L, 1L, 0L),
    `A*02:06` = c(0L, 0L, 0L, 1L, 0L),
    `A*26:01` = c(0L, 0L, 0L, 0L, 2L),
    `B*07:02` = c(0L, 0L, 0L, 0L, 2L),
    `B*08:01` = c(0L, 0L, 2L, 0L, 0L),
    `B*13:02` = c(1L, 0L, 0L, 0L, 0L),
    `B*15:01` = c(1L, 0L, 0L, 0L, 0L),
    `B*27:05` = c(0L, 0L, 0L, 1L, 0L),
    `B*40:01` = c(0L, 1L, 0L, 1L, 0L),
    `B*57:01` = c(0L, 1L, 0L, 0L, 0L),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  expect_equal(wide, test_wide)

  expect_error(midasToWide(midas, 1), "experiment is not a character vector")

  expect_error(
    midasToWide(midas, "foo"),
    "experiment should match values \"hla_alleles\"."
  )
})

test_that("isExperimentCountsOrZeros", {
  hla_calls_file <- system.file("extdata", "HLAHD_output_example.txt", package = "MiDAS")
  hla_calls <- readHlaCalls(hla_calls_file)[1:5, 1:5]

  pheno_file <- system.file("extdata", "pheno_example.txt", package = "MiDAS")
  pheno <- read.table(pheno_file, header = TRUE, stringsAsFactors = FALSE)[1:5, ]

  midas <- prepareMiDAS(
    hla_calls = hla_calls,
    colData = pheno,
    inheritance_model = "additive",
    experiment = c("hla_alleles", "hla_aa")
  )

  expect_equal(isExperimentCountsOrZeros(midas[["hla_alleles"]]), TRUE)

  expect_equal(isExperimentCountsOrZeros(midas[["hla_aa"]]), TRUE)

  expect_equal(isExperimentCountsOrZeros(matrix(runif(15), nrow = 3)), FALSE)

  expect_equal(isExperimentCountsOrZeros(LETTERS), FALSE)
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
      inheritance_model = "dominant",
      experiment = "hla_aa"
    )
  MiDASdat <- filterByOmnibusGroups(MiDASdat, "hla_aa", c("A_29", "A_44", "A_65"))
  omnibus_groups <- getOmnibusGroups(MiDASdat, "hla_aa")
  placeholder <- getPlaceholder(MiDASdat)
  MiDASdat <- as.data.frame(MiDASdat)
  object <- lm(disease ~ outcome + term, data = MiDASdat)

  res <- iterativeLRT(object, placeholder, omnibus_groups)
  test_res <- data.frame(
    group = c("A_29", "A_44", "A_65"),
    term = c("A_29_D, A_29_A", "A_44_R, A_44_K", "A_65_R, A_65_G"),
    dof = c(1, 3, 3),
    logLik = c(16365.357381683, 16365.5053334948, 16366.7408713989),
    statistic = c(32730.714763366, 32731.0106669896, 32733.4817427978),
    p.value = c(0, 0, 0),
    stringsAsFactors = FALSE
  )
  expect_equal(res, test_res)

  MiDASdat$A_29_A <- NA
  res <- iterativeLRT(object, placeholder, omnibus_groups)
  test_res <- data.frame(
    group = c("A_29", "A_44", "A_65"),
    term = c("A_29_D, A_29_A", "A_44_R, A_44_K", "A_65_R, A_65_G"),
    dof = c(NA, 3, 3),
    logLik = c(NA, 16365.5053334948, 16366.7408713989),
    statistic = c(NA, 32731.0106669896, 32733.4817427978),
    p.value = c(NA, 0, 0),
    stringsAsFactors = FALSE
  )
  expect_equal(res, test_res)
})

test_that("iterativeModel", {
  MiDASdat <-
    prepareMiDAS(
      hla_calls = MiDAS_tut_HLA[, 1:3],
      colData = MiDAS_tut_pheno,
      inheritance_model = "dominant",
      experiment = "hla_alleles"
    )
  placeholder <- getPlaceholder(MiDASdat)
  variables <- rownames(MiDASdat)[["hla_alleles"]][1:3]
  MiDASdat <- as.data.frame(MiDASdat)
  object <- lm(disease ~ outcome + term, data = MiDASdat)

  res <- iterativeModel(object, placeholder, variables)
  res_test <- dplyr::tibble(
    term = c("A*01:01:01", "A*01:01:41", "A*01:01:47"),
    estimate = c(-2.50126952397595e-16, -4.22559004808231e-16, NA),
    std.error = c(4.79744241693974e-16, 4.45423469957002e-15, NA),
    statistic = c(-0.521375621965567, -0.0948668027863513, NA),
    p.value = c(0.60233737113293, 0.924458860404277, NA)
  )
  expect_equal(as.data.frame(res), as.data.frame(res_test))
})
