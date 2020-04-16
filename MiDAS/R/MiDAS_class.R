#' @include MiDAS_class.R

#' MiDAS class
#'
#' The \code{MiDAS} class is a \code{\link{MultiAssayExperiment}} object
#' containing all the data and metadata about a set of measurments and thier
#' transformations required for MiDAS analysis.
#'
#' A object of class \code{MiDAS} have at least one of
#' \link[=readHlaCalls]{HLA calls} or \link[=readKirCalls]{KIR calls} experiment
#' matrices defined (they can be accessed using \code{\link{getHlaCalls}} and
#' \code{\link{getKirCalls}} functions).
#'
#' Other experiments slots can be occupied by transformations of those matrices.
#' Those transformations can be further used by \code{\link{runMiDAS}} function
#' for performing statistical analyses. See \code{\link{prepareMiDAS}} and
#' \code{\link{runMiDAS}} for more informations.
#'
#' An inheritance model of \code{MiDAS} object is held in object's metadata. It
#' specifies inheritance model under which all experiments matrices has been
#' constructed. It can be accessed using \code{\link{getInheritanceModel}}
#' function.
#'
#' Experiments slots available for \code{\link{runMiDAS}} function are defined
#' in object's metadata under \code{analysis_type} variable. It can be accessed
#' using \code{\link{getAnalysisType}} function.
#'
#' @importFrom methods setClass
#' @importClassesFrom MultiAssayExperiment MultiAssayExperiment
#'
#' @export
MiDAS <- setClass(
  Class = "MiDAS",
  contains = "MultiAssayExperiment"
)

#' Validity method for class MiDAS
#'
#' Valid \code{MiDAS} object must contain hla_calls or kir_calls experiment in
#' proper \link[=readHlaCalls]{HLA calls} or \link[=readKirCalls]{KIR calls}
#' format. It also have to have \code{inheritance_model} metadata's variable
#' defined, accepted values are \code{"additive"}, \code{"dominant"},
#' \code{"recessive"}.
#'
#' \code{analysis_type} metadata's variable used to determine analyses available
#' via \code{\link{runMiDAS}} is optional and does not determine object
#' validity.
#'
#' @importFrom assertthat assert_that is.string see_if
#'
setValidity(Class = "MiDAS", method = function(object) {
  hla_calls <- getHlaCalls(object)
  kir_calls <- getKirCalls(object)
  inheritance_model <- getInheritanceModel(object)

  assert_that(
    see_if(! is.null(hla_calls) | ! is.null(kir_calls),
           msg = "MiDAS object must contain hla_calls or kir_calls experiment"
    ),
    if (! is.null(hla_calls)) { checkHlaCallsFormat(hla_calls) } else { TRUE },
    if (! is.null(kir_calls)) { checkKirCallsFormat(kir_calls) } else { TRUE },
    is.string(inheritance_model),
    stringMatches(
      x = inheritance_model,
      choice = c("additive", "dominant", "recessive")
    )
  )

  return(TRUE)
})

#' @name getInheritanceModel
#'
#' @title Extract inheritance model from MiDAS object.
#'
#' @param object \code{\link{MiDAS}} object
#'
#' @return String name of object's inheritance model.
#'
#' @export
setGeneric(
  name = "getInheritanceModel",
  def = function(object) standardGeneric("getInheritanceModel")
)

#' @rdname getInheritanceModel
#'
#' @importFrom S4Vectors metadata
#'
setMethod(
  f = "getInheritanceModel",
  signature = "MiDAS",
  definition = function (object) metadata(object)$inheritance_model
)

#' @name getAnalysisType
#'
#' @title Extract available analysis types from MiDAS object.
#'
#' @param object \code{\link{MiDAS}} object
#'
#' @return Character giving names of available analysis types.
#'
#' @export
setGeneric(
  name = "getAnalysisType",
  def = function(object) standardGeneric("getAnalysisType")
)

#' @rdname getAnalysisType
#'
#' @importFrom S4Vectors metadata
#'
setMethod(
  f = "getAnalysisType",
  signature = "MiDAS",
  definition = function (object) metadata(object)$analysis_type
)

#' @name getHlaCalls
#'
#' @title Extract HLA calls from MiDAS object.
#'
#' @param object \code{\link{MiDAS}} object
#'
#' @return Data frame with object's HLA calls.
#'
#' @export
setGeneric(
  name = "getHlaCalls",
  def = function(object) standardGeneric("getHlaCalls")
)

#' @rdname getHlaCalls
#'
setMethod(
  f = "getHlaCalls",
  signature = "MiDAS",
  definition = function (object) {
    hla_calls <- object[["hla_calls"]]
    if (! is.null(hla_calls)) {
      hla_calls <- experimentMatToDf(hla_calls)
    }

    return(hla_calls)
  }
)

#' @name getKirCalls
#'
#' @title Extract KIR calls from MiDAS object.
#'
#' @param object \code{\link{MiDAS}} object.
#'
#' @return Data frame with object's KIR calls.
#'
#' @export
setGeneric(
  name = "getKirCalls",
  def = function(object) standardGeneric("getKirCalls")
)

#' @rdname getKirCalls
#'
setMethod(
  f = "getKirCalls",
  signature = "MiDAS",
  definition = function (object) {
    kir_calls <- object[["kir_calls"]]
    if (! is.null(kir_calls)) {
      kir_calls <- experimentMatToDf(kir_calls)
    }

    return(kir_calls)
  }
)

#' Prepare MiDAS object for statistical analysis
#'
#' \code{prepareMiDAS} transform HLA alleles calls and KIR calls according
#' to selected analysis creating \code{\link{MiDAS}} object.
#'
#' @inheritParams checkHlaCallsFormat
#' @inheritParams checkKirCallsFormat
#' @inheritParams hlaCallsToCounts
#' @param phenotype Data frame holding additional variables like phenotypic
#'   observations or covariates. It have to contain \code{'ID'} column holding
#'   samples indentifires corresponding to indentifires in \code{hla_calls} and
#'   \code{kir_calls}. Importantly rows of \code{hla_calls} and
#'   \code{kir_calls} without corresponding phenotype are discarded.
#' @param analysis_type Character vector indicating analysis type for which data
#'   should be prepared. Valid choices are \code{"hla_allele"},
#'   \code{"aa_level"}, \code{"allele_g_group"}, \code{"allele_supertype"},
#'   \code{"allele_group"}, \code{"kir_genes"}, \code{"hla_kir_interactions"}.
#'   See details for further explanations.
#'
#' @param ... Attributes used in data transformation.
#'
#' @return Object of class \code{\link{MiDAS}}
#'
#' @examples
#' # read hla calls file
#' hla_calls_file <- system.file("extdata", "HLAHD_output_example.txt", package = "MiDAS")
#' hla_calls <- readHlaCalls(hla_calls_file)
#'
#' # read phenotypic data and covariates
#' pheno_file <- system.file("extdata", "pheno_example.txt", package = "MiDAS")
#' pheno <- read.table(pheno_file, header = TRUE, stringsAsFactors = FALSE)
#' covar_file <- system.file("extdata", "covar_example.txt", package = "MiDAS")
#' covar <- read.table(covar_file, header = TRUE, stringsAsFactors = FALSE)
#' phenotype <- left_join(pheno, covar, by ="ID")
#'
#' # create MiDAS object
#' midas <- MiDAS(hla_calls = hla_calls,
#'                phenotype = phenotype,
#'                inheritance_model = "additive",
#'                analysis_type = "hla_allele"
#' )
#'
#' @importFrom assertthat assert_that see_if
#' @importFrom MultiAssayExperiment ExperimentList MultiAssayExperiment
#' @importFrom S4Vectors DataFrame
#'
#' @export
prepareMiDAS <- function(hla_calls = NULL,
                         kir_calls = NULL,
                         phenotype,
                         inheritance_model = c("additive", "dominant", "recessive"),
                         analysis_type = c(
                           "hla_allele",
                           "aa_level",
                           "allele_g_group",
                           "allele_supertype",
                           "allele_group",
                           "kir_genes",
                           "hla_kir_interactions"
                         ),
                         ...
) {
  inheritance_model_choice <- eval(formals()[["inheritance_model"]])
  analysis_type_choice <- eval(formals()[["analysis_type"]])
  assert_that(
    see_if(
      ! is.null(hla_calls) | ! is.null(kir_calls),
      msg = "hla_calls or kir_calls argument has to be specified"
    ),
    if (! is.null(hla_calls)) {
      checkHlaCallsFormat(hla_calls)
    } else { TRUE },
    if (! is.null(kir_calls)) {
      checkKirCallsFormat(kir_calls)
    } else { TRUE },
    checkPhenotypeFormat(phenotype),
    is.string(inheritance_model),
    stringMatches(inheritance_model, inheritance_model_choice),
    is.character(analysis_type),
    characterMatches(analysis_type, analysis_type_choice)
  )

  dot.args <- if (...length()) {
    as.list(...)
  } else {
    list()
  }

  experiments <- ExperimentList()

  if (! is.null(hla_calls)) {
    experiments[["hla_calls"]] <- dfToExperimentMat(hla_calls)
  }

  if (! is.null(kir_calls)) {
    experiments[["kir_calls"]] <- dfToExperimentMat(kir_calls)
  }

  # prepare data for different analyses types
  for (at in analysis_type) {
    fun <- paste0("prepareMiDAS_", at)
    args <- list(
      hla_calls = hla_calls,
      kir_calls = kir_calls,
      inheritance_model = inheritance_model
    )
    args <- c(args, dot.args)
    experiment <- do.call(
      what = fun,
      args = args
    )
    experiments[[at]] <- experiment
  }

  phenotype <-
    DataFrame(phenotype, row.names = phenotype[["ID"]], check.names = TRUE)

  metadata <- list(
    inheritance_model = inheritance_model,
    analysis_type = analysis_type
  )

  mae <- MultiAssayExperiment(
    experiments = experiments,
    colData = phenotype,
    metadata = metadata
  )

  class(mae) <- structure("MiDAS", package = "MiDAS")

  return(mae)
}


#' @rdname prepareMiDAS
#'
#' @title Prepare MiDAS data on HLA allele level
#'
#' @details \code{'hla_allele'} - \code{hla_calls} are transformed into counts
#' under \code{inheritance_model} of choice (see \code{\link{hlaCallsToCounts}}
#' for  more details).
#'
#' @inheritParams prepareMiDAS
#'
#' @return Matrix
#'
prepareMiDAS_hla_allele <- function(hla_calls, inheritance_model, ...) {
  hla_allele <- hlaCallsToCounts(
    hla_calls = hla_calls,
    inheritance_model = inheritance_model
  ) %>%
    dfToExperimentMat()

  return(hla_allele)
}

#' Prepare MiDAS data on HLA amino acid level
#'
#' \code{hla_calls} are first converted to amino acid level, taking only
#' variable positions under consideration. Than variable amino acid positions
#' are transformed to counts under \code{inheritance_model} of choice (see
#' \code{\link{hlaToAAVariation}} and \code{\link{aaVariationToCounts}} for more
#' details).
#'
#' @inheritParams aaVariationToCounts
#' @param hla_calls Data frame
#' @param inheritance_model String
#' @param ... Not used
#'
#' @return Matrix
#'
#' @importFrom assertthat assert_that is.flag
#'
prepareMiDAS_aa_level <- function(hla_calls,
                                  inheritance_model,
                                  indels = TRUE,
                                  unkchar = FALSE,
                                  ...) {
  assert_that(
    is.flag(indels),
    is.flag(unkchar)
  )

  hla_aa <-
    hlaToAAVariation(
      hla_calls = hla_calls,
      indels = indels,
      unkchar = unkchar
    ) %>%
    aaVariationToCounts(inheritance_model = inheritance_model) %>%
    dfToExperimentMat()

  return(hla_aa)
}

#' Prepare MiDAS data on HLA allele's G groups level
#'
#' \code{hla_calls} are transformed to HLA alleles groups using G group
#' dictionary shipped with the package. Than they are transformed to counts and
#' \code{\link{hlaCallsToCounts}} for more details).
#'
#' @param hla_calls Data frame
#' @param inheritance_model String
#' @param ... Not used
#'
#' @return Matrix
#'
#' @importFrom assertthat assert_that
#'
prepareMiDAS_allele_g_group <- function(hla_calls,
                                        inheritance_model,
                                        ...) {
  lib <- "allele_HLA_Ggroup"
  allele_g_group <- hlaToVariable(hla_calls = hla_calls,
                                  dictionary = lib,
                                  na.value = 0
  )

  assert_that(
    ncol(allele_g_group) > 1,
    msg = "no allele could be assigned to G group for input hla_calls"
  )

  allele_g_group <-
    hlaCallsToCounts(
      hla_calls = allele_g_group,
      inheritance_model = inheritance_model,
      check_hla_format = FALSE
    ) %>%
    dfToExperimentMat()

  return(allele_g_group)
}

#' Prepare MiDAS data on HLA allele's supertypes level
#'
#' \code{hla_calls} are transformed to HLA alleles groups using supertypes
#' dictionary shipped with the package. Than they are transformed to counts
#' under \code{inheritance_model} of choice (see \code{\link{hlaToVariable}} and
#' \code{\link{hlaCallsToCounts}} for more details).
#'
#' @param hla_calls Data frame
#' @param inheritance_model String
#' @param ... Not used
#'
#' @return Matrix
#'
#' @importFrom assertthat assert_that
#'
prepareMiDAS_allele_supertype <- function(hla_calls, inheritance_model, ...) {
  lib <- "allele_HLA_supertype"
  allele_supertype <- hlaToVariable(hla_calls = hla_calls,
                                    dictionary = lib,
                                    na.value = 0
  )

  assert_that(
    ncol(allele_supertype) > 1,
    msg = "no allele could be assigned to supertype for input hla_calls"
  )

  allele_supertype <-
    hlaCallsToCounts(
      hla_calls = allele_supertype,
      inheritance_model = inheritance_model,
      check_hla_format = FALSE
    ) %>%
    subset(select = - Unclassified) %>%
    dfToExperimentMat()

  return(allele_supertype)
}

#' Prepare MiDAS data on HLA allele's groups level
#'
#' \code{hla_calls} are transformed to HLA alleles groups using Bw4/6, C1/2 and
#' Bw4+A23+A24+A32 dictionaries shipped with the package. Than they are
#' transformed to counts under \code{inheritance_model} of choice (see
#' \code{\link{hlaToVariable}} and \code{\link{hlaCallsToCounts}} for more
#' details).
#'
#' @param hla_calls Data frame
#' @param inheritance_model String
#' @param ... Not used
#'
#' @return Matrix
#'
#' @importFrom assertthat assert_that
#'
prepareMiDAS_allele_group <- function(hla_calls, inheritance_model, ...) {
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

  allele_group <-
    hlaCallsToCounts(
      hla_calls = allele_group,
      inheritance_model = inheritance_model,
      check_hla_format = FALSE
    ) %>%
    dfToExperimentMat()

  return(allele_group)
}

#' Prepare MiDAS data on KIR genes level
#'
#' \code{kir_counts} data frame is joined with other inputs.
#'
#' @param kir_calls Data frame
#' @param ... Not used
#'
#' @return Matrix
#'
prepareMiDAS_kir_genes <- function(kir_calls, ...) {
  kir_genes <-
    kir_calls %>%
    dfToExperimentMat()

  return(kir_genes)
}

#' Prepare MiDAS data on HLA - KIR interactions level
#'
#' \code{hla_calls} are processed with \code{kir_counts} into HLA - KIR
#' interactions variables (see \code{\link{getHlaKirInteractions}} for more
#' details).
#'
#' @param hla_calls Data frame
#' @param kir_calls Data frame
#' @param ... Not used
#'
#' @return Matrix
#'
prepareMiDAS_hla_kir_interactions <- function(hla_calls, kir_calls, ...) {
  hla_kir_interactions <-
    getHlaKirInteractions(
      hla_calls = hla_calls,
      kir_counts = kir_calls
    ) %>%
    dfToExperimentMat()

  return(hla_kir_interactions)
}

#' Filter midas object
#'
#' @inheritParams MiDAS
#' @param filter Character giving suffix of filter function.
#'
#' @return Object of class MiDAS.
#'
#' @importFrom assertthat assert_that
#'
filterMiDAS <- function(object, filter_by = c("hla_allele_frequency"), ...) {
  filter_by_choices <- eval(formals()[["filter_by"]])
  assert_that(
    validObject(object),
    is.character(filter_by),
    characterMatches(filter_by, filter_by_choices)
  )

  # get informations required for recreating MiDAS object
  midas_components <- list(
    inheritance_model = getInheritanceModel(object),
    hla_calls = getHlaCalls(object),
    kir_calls = getKirCalls(object),
    analysis_types = names(object)
  )

  # for (by in filter_by) {
  #
  # }
}

#' Transform MiDAS to wide format data.frame
#'
#' @param object Object of class MiDAS
#' @param analysis_type Character
#' @param placeholder String
#'
#' @return Data frame
#'
#' @importFrom assertthat assert_that
#' @importFrom dplyr filter mutate
#' @importFrom MultiAssayExperiment longFormat
#' @importFrom rlang !! :=
#' @importFrom tidyr spread
#'
midasToWide <- function(object, analysis_type, placeholder = "term") {
  assert_that(
    validObject(object),
    is.character(analysis_type),
    characterMatches(analysis_type, names(object))
  )

  wide_df <-
    object[, , analysis_type] %>%
    longFormat(colDataCols = TRUE) %>%
    as.data.frame() %>%
    subset(select = -assay) %>%
    subset(select = -colname) %>%
    spread(key = "rowname", value = "value") %>%
    mutate(!!placeholder := 1)

  return(wide_df)
}

#' Helper transform data frame to experiment matrix
#'
#' Function deletes 'ID' column from a \code{df}, then transpose it and sets
#' the columnames to values from deleted 'ID' column.
#'
#' @param df Data frame
#'
#' @return Matrix
#'
dfToExperimentMat <- function(df) {
  cols <- df[["ID"]]
  mat <- t(subset(df, select = -ID))
  colnames(mat) <- cols

  return(mat)
}

#' Helper transform experimentt matrix to data frame
#'
#' Function transpose \code{mat} and inserts columnames of input \code{mat} as
#' a 'ID' column.
#'
#' @param mat Matrix
#'
#' @return Data frame
#'
experimentMatToDf <- function(mat) {
  ID <- colnames(mat)
  df <-
    as.data.frame(
      t(mat),
      stringsAsFactors = FALSE,
      make.names = FALSE
    )
  rownames(df) <- NULL
  df <- cbind(ID, df, stringsAsFactors = FALSE)

  return(df)
}
