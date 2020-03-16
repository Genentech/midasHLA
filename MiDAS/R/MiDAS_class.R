#' MiDAS class definition
#'
#' @importFrom methods setClass
#' @importClassesFrom MultiAssayExperiment MultiAssayExperiment
#'
#' @export
.MiDAS <- setClass(
  Class = "MiDAS",
  slots = c(inheritance_model = "character"),
  contains = "MultiAssayExperiment"
)

#' MiDAS class constructor
#'
#' @importFrom MultiAssayExperiment MultiAssayExperiment
#' @importFrom magrittr %>% %<>%
#'
#' @export
MiDAS <- function(hla_calls = NULL,
                  kir_counts = NULL,
                  colData = NULL,
                  inheritance_model = c("additive", "dominant", "recesive"),
                  analysis_type = c("hla_allele", "aa_level", "expression_level", "allele_g_group", "allele_supertype", "allele_group", "kir_genes", "hla_kir_interactions"),
                  ...) {

  experiments <- list()

  if ("hla_allele" %in% analysis_type) {
    hla_allele <- hlaCallsToCounts(
      hla_calls = hla_calls,
      inheritance_model = inheritance_model
    )

    hla_allele <- t(hla_allele[, -1]) # TODO w/o ID column
    colnames(hla_allele) <- hla_calls$ID # TODO
    experiments[["hla_allele"]] <- hla_allele
  }

  if ("aa_level" %in% analysis_type) {
    aa_level <- hlaToAAVariation(
      hla_calls = hla_calls,
      indels = TRUE,
      unkchar = FALSE
    ) %>%
      aaVariationToCounts(inheritance_model = inheritance_model)

    aa_level <- t(aa_level[, -1]) # TODO w/o ID column
    colnames(aa_level) <- hla_calls$ID # TODO
    experiments[["aa_level"]] <- hla_allele
  }

  if ("expression_level" %in% analysis_type) {
    lib <- listMiDASDictionaries()
    lib <- grep("expression", lib, value = TRUE)
    expression_level <- Reduce(
      f = function(...) left_join(..., by = "ID"),
      x = lapply(lib, hlaToVariable, hla_calls = hla_calls, na.value = NA)
    )

    assert_that(
      ncol(expression_level) > 1,
      msg = "no expression levels were found for input hla_calls"
    )

    expression_level <- expression_level %>%
      gather("expression", "value", -c("ID")) %>%
      mutate(expression = gsub("_[12]", "", .data$expression)) %>%
      group_by(!!! syms(c("ID", "expression"))) %>%
      summarise_all(funs(sum)) %>%
      spread(.data$expression, .data$value, sep = NULL)

    expression_level <- t(expression_level[, -1]) # TODO w/o ID column
    colnames(expression_level) <- hla_calls$ID # TODO
    experiments[["expression_level"]] <- expression_level
  }

  if ("allele_g_group" %in% analysis_type) {
    lib <- "allele_HLA_Ggroup"
    allele_g_group <- hlaToVariable(hla_calls = hla_calls,
                                    dictionary = lib,
                                    na.value = 0
    )

    assert_that(
      ncol(allele_g_group) > 1,
      msg = "no allele could be assigned to G group for input hla_calls"
    )

    allele_g_group <- hlaCallsToCounts(
      hla_calls = allele_g_group,
      inheritance_model = inheritance_model,
      check_hla_format = FALSE
    )

    allele_g_group <- t(allele_g_group[, -1]) # TODO w/o ID column
    colnames(allele_g_group) <- hla_calls$ID # TODO
    experiments[["allele_g_group"]] <- allele_g_group
  }

  if ("allele_supertype" %in% analysis_type) {
    lib <- "allele_HLA_supertype"
    allele_supertype <- hlaToVariable(hla_calls = hla_calls,
                                      dictionary = lib,
                                      na.value = 0
    )

    assert_that(
      ncol(allele_supertype) > 1,
      msg = "no allele could be assigned to supertype for input hla_calls"
    )

    allele_supertype <- hlaCallsToCounts(
      hla_calls = allele_supertype,
      inheritance_model = inheritance_model,
      check_hla_format = FALSE
    ) %>%
      select(-"Unclassified")

    allele_supertype <- t(allele_supertype[, -1]) # TODO w/o ID column
    colnames(allele_supertype) <- hla_calls$ID # TODO
    experiments[["allele_supertype"]] <- allele_supertype # TODO w/o ID column
  }

  if ("allele_group" %in% analysis_type) {
    lib <- c(
      "allele_HLA-B_Bw",
      "allele_HLA_Bw4+A23+A24+A32",
      "allele_HLA-C_C1-2"
    )
    allele_group <- Reduce(
      f = function(...) left_join(..., by = "ID"),
      x = lapply(lib, hlaToVariable, hla_calls = hla_calls, na.value = 0)
    )

    assert_that(
      ncol(allele_group) > 1,
      msg = "no allele could be assigned to allele groups for input hla_calls"
    )

    allele_group <- hlaCallsToCounts(
      hla_calls = allele_group,
      inheritance_model = inheritance_model,
      check_hla_format = FALSE
    )

    allele_group <- t(allele_group[, -1]) # TODO w/o ID column
    colnames(allele_group) <- hla_calls$ID # TODO
    experiments[["allele_group"]] <- allele_group
  }

  if ("kir_genes" %in% analysis_type) {
    assert_that(
      ! missing(kir_counts),
      msg = "\"kir_genes\" analysis type requires kir_counts argument to be specified"
    )

    kir_counts <- t(kir_counts[, -1]) # TODO w/o ID column
    colnames(kir_counts) <- kir_counts$ID # TODO
    experiments[["kir_genes"]] <- kir_counts
  }

  if ("hla_kir_interactions" %in% analysis_type) {
    assert_that(
      ! missing(kir_counts),
      msg = "\"hla_kir_interactions\" analysis type requires kir_counts argument to be specified"
    )

    hla_kir_interactions <- getHlaKirInteractions(
      hla_calls = hla_calls,
      kir_counts = kir_counts
    )

    hla_kir_interactions <- t(hla_kir_interactions[, -1]) # TODO w/o ID column
    colnames(hla_kir_interactions) <- hla_kir_interactions$ID # TODO
    experiments[["hla_kir_interactions"]] <- hla_kir_interactions
  }

  midas <- MultiAssayExperiment(
    experiments = experiments,
    colData = colData
  )

  midas <- .MiDAS(midas, inheritance_model = inheritance_model)

  return(midas)
}

#' Transform MiDAS to wide format data.frame
#'
#' @importFrom tidyr spread
#'
midasToWide <- function(object) {
  longFormat(object, colDataCols = TRUE) %>%
    as.data.frame() %>%
    subset(select = -assay) %>%
    subset(select = -colname) %>%
    tidyr::spread(key = "rowname", value = "value") %>%
    mutate(term = 1)
}

library("survival")
hla_calls_file <- system.file("extdata", "HLAHD_output_example.txt", package = "MiDAS")
hla_calls <- readHlaCalls(hla_calls_file)
pheno_file <- system.file("extdata", "pheno_example.txt", package = "MiDAS")
pheno <- read.table(pheno_file, header = TRUE, stringsAsFactors = FALSE)
covar_file <- system.file("extdata", "covar_example.txt", package = "MiDAS")
covar <- read.table(covar_file, header = TRUE, stringsAsFactors = FALSE)
coldata <- left_join(pheno, covar, by ="ID")
rownames(coldata) <- coldata$ID
midas <- MiDAS(hla_calls = hla_calls,
                           colData = coldata,
                           analysis_type = "hla_allele",
                           inheritance_model = "additive"
)
midas_data <- midasToWide(midas)
object <- coxph(Surv(OS, OS_DIED) ~ AGE + SEX + term, data = midas_data)
runMiDAS(object, analysis_type = "hla_allele", variables = rownames(midas[[1]]))
