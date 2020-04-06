#' MiDAS class to hold all data and meta data about one MiDAS analysis.
#'
#' The \link{\code{MiDAS}} class is a \link{\code{MultiAssayExperiment}} object
#' containing all the data and metadata about a set of measurments and thier
#' transformations required for MiDAS analysis.
#'
#' @slot metadata A list that must at least contain a \code{inheritance_model}
#'   element.
#'
#' The MiDAS object must contain at least \code{hla_calls} (see
#' \link{\code{readHlaCalls}}) or \code{kir_calls} (see
#' \link{\code{readKirCalls}}) experiments and defined \code{inheritance_model}
#' element of \code{metadata}.
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
#' @importFrom assertthat assert_that is.string
#'
setValidity(Class = "MiDAS", method = function(object) {
  assert_that(
    if (! is.null(getHlaCalls(midas))) {
      checkHlaCallsFormat(getHlaCalls(midas))
    } else { TRUE },
    if (! is.null(getKirCalls(midas))) {
      checkKirCallsFormat(getKirCalls(midas))
    } else { TRUE },
    is.string(getInheritanceModel(midas)),
    stringMatches(
      x = getInheritanceModel(midas),
      choice = c("additive", "dominant", "recessive")
    )
  )

  return(TRUE)
})

#' Helper crating object of class MiDAS
#'
#' @examples
#' # read hla calls file
#' hla_calls_file <- system.file("extdata", "HLAHD_output_example.txt", package = "MiDAS")
#' hla_calls <- readHlaCalls(hla_calls_file)
#'
#' # read batch covariates
#' pheno_file <- system.file("extdata", "pheno_example.txt", package = "MiDAS")
#' pheno <- read.table(pheno_file, header = TRUE, stringsAsFactors = FALSE)
#' covar_file <- system.file("extdata", "covar_example.txt", package = "MiDAS")
#' covar <- read.table(covar_file, header = TRUE, stringsAsFactors = FALSE)
#' colData <- left_join(pheno, covar, by ="ID")
#' rownames(colData) <- colData$ID

#' # create MiDAS object
#' midas <- MiDAS(hla_calls = hla_calls,
#'                colData = colData,
#'                inheritance_model = "additive",
#'                analysis_type = "hla_allele"
#' )
#'
#'
#'
#' @importFrom assertthat assert_that
#' @importFrom methods as
#' @importFrom MultiAssayExperiment ExperimentList MultiAssayExperiment
#'
#' @export
MiDAS <- function(hla_calls = NULL,
                  kir_calls = NULL,
                  colData,
                  inheritance_model = c("additive", "dominant", "recesive"),
                  analysis_type = c("hla_allele", "aa_level", "expression_level",
                                    "allele_g_group", "allele_supertype",
                                    "allele_group", "kir_genes",
                                    "hla_kir_interactions", "custom"),
                  ...) {
  inheritance_model_choice <- eval(formals()[["inheritance_model"]])
  analysis_type_choice <- eval(formals()[["analysis_type"]])
  assert_that(
    if (! is.null(hla_calls)) {
      checkHlaCallsFormat(hla_calls)
    } else { TRUE },
    if (! is.null(kir_calls)) {
      checkKirCallsFormat(kir_calls)
    } else { TRUE },
    is.data.frame(colData),
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

  metadata <- list(
    hla_calls = hla_calls,
    kir_calls = kir_calls,
    inheritance_model = inheritance_model
  )

  mae <- MultiAssayExperiment(
    experiments = as(experiments, "ExperimentList"),
    colData = as(colData, "DataFrame"),
    metadata = metadata
  )

  class(mae) <- structure("MiDAS", package = "MiDAS")

  return(mae)
}

#' @name getInheritanceModel
#'
#' @title Extract inheritance model from MiDAS  object.
#'
#' @param object \link{\code{MiDAS}} object
#'
#' @return String giving name of object's inheritance model.
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

#' @name getHlaCalls
#'
#' @title Extract HLA calls from MiDAS  object.
#'
#' @param object \link{\code{MiDAS}} object
#'
#' @return String giving name of object's inheritance model.
#'
#' @export
setGeneric(
  name = "getHlaCalls",
  def = function(object) standardGeneric("getHlaCalls")
)

#' @rdname getHlaCalls
#'
#' @importFrom S4Vectors metadata
#'
setMethod(
  f = "getHlaCalls",
  signature = "MiDAS",
  definition = function (object) metadata(object)$hla_calls
)

#' @name getKirCalls
#'
#' @title Extract KIR calls from MiDAS  object.
#'
#' @param a
#'
#' @return String giving name of object's inheritance model.
#'
#' @export
setGeneric(
  name = "getKirCalls",
  def = function(object) standardGeneric("getKirCalls")
)

#' @rdname getKirCalls
#'
#' @importFrom S4Vectors metadata
#'
setMethod(
  f = "getKirCalls",
  signature = "MiDAS",
  definition = function (object) metadata(object)$kir_calls
)

#' Prepare MiDAS data on HLA allele level
#'
#' @param hla_calls Data frame
#' @param inheritance_model String
#' @param ... Not used
#'
#' @return Matrix
#'
prepareMiDAS_hla_allele <- function(hla_calls, inheritance_model, ...) {
  hla_allele <- hlaCallsToCounts(hla_calls = hla_calls,
                                 inheritance_model = inheritance_model)
  hla_allele <- t(hla_allele[, -1]) # TODO w/o ID column
  colnames(hla_allele) <- hla_calls$ID # TODO

  return(hla_allele)
}

#' Transform MiDAS to wide format data.frame
#'
#' @importFrom assertthat assert_that
#' @importFrom dplyr filter mutate
#' @importFrom MultiAssayExperiment longFormat
#' @importFrom tidyr spread
#'
midasToWide <- function(object, analysis_type) {
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
    mutate(term = 1)

  return(wide_df)
}
