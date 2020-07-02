#' @include MiDAS_class.R
NULL

#' MiDAS class
#'
#' The \code{MiDAS} class is a \code{\link{MultiAssayExperiment}} object
#' containing all the data and metadata about a set of measurments and thier
#' transformations required for MiDAS analysis.
#'
#' A object of class \code{MiDAS} have at least one of
#' \link[=readHlaCalls]{HLA calls} or \link[=readKPICalls]{KIR calls} experiment
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
#' in object's metadata under \code{experiment} variable. It can be accessed
#' using \code{\link{getExperiments}} function.
#'
#' @importFrom methods setClass new
#' @importClassesFrom MultiAssayExperiment MultiAssayExperiment
#'
#' @export
MiDAS <- setClass(
  Class = "MiDAS",
  contains = "MultiAssayExperiment"
)

# Validity method for class MiDAS
#
# Valid \code{MiDAS} object must contain hla_calls or kir_calls experiment in
# proper \link[=readHlaCalls]{HLA calls} or \link[=readKPICalls]{KIR calls}
# format. It also have to have \code{inheritance_model} metadata's variable
# defined, accepted values are \code{"additive"}, \code{"dominant"},
# \code{"recessive"}.
#
# \code{experiment} metadata's variable used to determine analyses available
# via \code{\link{runMiDAS}} is optional and does not determine object
# validity.
#
# @importFrom assertthat assert_that is.string see_if
#
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
    ),
    is.string(getPlaceholder(object)),
    see_if(
      ! getPlaceholder(object) %in% unlist(rownames(object)),
      msg = sprintf("Placeholder '%s' is used in one of the experiments", getPlaceholder(object))
    ),
    see_if(
      getPlaceholder(object) %in% colnames(colData(object)),
      msg = sprintf("Placeholder '%s' can not be found in object's colData", getPlaceholder(object))
    )
  )

  return(TRUE)
})

#' @rdname MiDAS-class
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

#' @rdname MiDAS-class
#'
#' @importFrom S4Vectors metadata
#'
setMethod(
  f = "getInheritanceModel",
  signature = "MiDAS",
  definition = function (object) metadata(object)$inheritance_model
)

#' @rdname MiDAS-class
#'
#' @title Extract available analysis types from MiDAS object.
#'
#' @param object \code{\link{MiDAS}} object
#'
#' @return Character giving names of available analysis types.
#'
#' @export
setGeneric(
  name = "getExperiments",
  def = function(object) standardGeneric("getExperiments")
)

#' @rdname MiDAS-class
#'
#' @importFrom S4Vectors metadata
#'
setMethod(
  f = "getExperiments",
  signature = "MiDAS",
  definition = function (object) metadata(object)$experiment
)

#' @rdname MiDAS-class
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

#' @rdname MiDAS-class
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

#' @rdname MiDAS-class
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

#' @rdname MiDAS-class
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

#' @rdname MiDAS-class
#'
#' @title Extract placeholder name from MiDAS object.
#'
#' @param object \code{\link{MiDAS}} object.
#'
#' @return String giving object's placeholder.
#'
#' @export
setGeneric(
  name = "getPlaceholder",
  def = function(object) standardGeneric("getPlaceholder")
)

#' @rdname MiDAS-class
#'
#' @importFrom S4Vectors metadata
#'
setMethod(
  f = "getPlaceholder",
  signature = "MiDAS",
  definition = function (object) {
    placeholder <- metadata(object)$placeholder

    return(placeholder)
  }
)

#' @rdname MiDAS-class
#'
#' @title Extract omnibus groups from MiDAS object.
#'
#' @param object \code{\link{MiDAS}} object.
#' @param experiment String specifying experiment.
#'
#' @return Character vector of omnibus groups for given experiment.
#'
#' @export
setGeneric(
  name = "getOmnibusGroups",
  def = function(object, experiment) standardGeneric("getOmnibusGroups")
)

#' @rdname MiDAS-class
#'
#' @importFrom assertthat assert_that is.string
#' @importFrom S4Vectors metadata
#'
setMethod(
  f = "getOmnibusGroups",
  signature = "MiDAS",
  definition = function (object, experiment) {
    assert_that(is.string(experiment))
    experiment <- object[[experiment]]
    omnibus_groups <- if (isClass(experiment, "SummarizedExperiment")) {
      metadata(experiment)$omnibus_groups
    } else {
      NULL
    }

    return(omnibus_groups)
  }
)

#' @rdname MiDAS-class
#'
#' @title Calculate variables frequencies for given experiment in MiDAS object.
#'
#' @param object \code{\link{MiDAS}} object.
#' @param experiment String specifying experiment.
#'
#' @return Data frame containing variables, its corresponding total counts
#'   and frequencies.
#'
#' Variables frequencies are counted in reference to sample size, depending on
#' the inheritance model under which the counts table has been generated one
#' might need to take under consideration both gene copies. Here sample size is
#' assumed to be depended on both gene copies for \code{"additive"} inheritance
#' model (`n / (2 * j)` where `n` is the number of term occurrences and `j`
#' is the sample size). For other models the sample size is taken as is
#' (`n / j`).
#'
#' @export
setGeneric(
  name = "getFrequencies",
  def = function(object, experiment) standardGeneric("getFrequencies")
)

#' @rdname MiDAS-class
#'
#' @importFrom assertthat assert_that is.string
#' @importFrom S4Vectors metadata
#'
setMethod(
  f = "getFrequencies",
  signature = "MiDAS",
  definition = function (object, experiment) {
    assert_that(
      is.string(experiment),
      stringMatches(experiment, getExperiments(object))
    )
    mat <- object[[experiment]]
    assert_that(
      is.matrix(mat) && isCountsOrZeros(mat),
      msg = sprintf("Frequencies can not be calculated for experiment '%s'",
                    experiment
            )
    )
    inheritance_model <- getInheritanceModel(object)
    freq <- getExperimentFrequencies(mat, inheritance_model)

    return(freq)
  }
)

#' @rdname MiDAS-class
#'
#' @title Filter MiDAS object by frqeuncy
#'
#' @param object \code{\link{MiDAS}} object
#' @param experiment Matrix
#' @param lower_frequency_cutoff Number
#' @param upper_frequency_cutoff Number
#'
#' @return Object of class MiDAS.
#'
#' @importFrom assertthat assert_that is.string see_if
#' @importFrom methods validObject
#'
#' @export
setGeneric(
  name = "filterByFrequency",
  def = function(object,
                 experiment,
                 lower_frequency_cutoff = NULL,
                 upper_frequency_cutoff = NULL) {
    standardGeneric("filterByFrequency")
  }
)

#' @rdname MiDAS-class
#'
#' @importFrom assertthat assert_that is.string
#' @importFrom S4Vectors metadata
#'
setMethod(
  f = "filterByFrequency",
  signature = "MiDAS",
  definition = function(object,
                        experiment,
                        lower_frequency_cutoff = NULL,
                        upper_frequency_cutoff = NULL) {
  assert_that(
    is.character(experiment),
    characterMatches(experiment, getExperiments(object)),
    validateFrequencyCutoffs(lower_frequency_cutoff, upper_frequency_cutoff)
  )

  inheritance_model <- getInheritanceModel(object)
  for (ex in experiment) {
    mat <- object[[ex]]
    assert_that(
      is.matrix(mat) && isCountsOrZeros(mat),
      msg = sprintf("Frequency filtration does not support experiment '%s'", ex)
    )
    object[[ex]] <- filterExperimentByFrequency(
      experiment = mat,
      inheritance_model = inheritance_model,
      lower_frequency_cutoff = lower_frequency_cutoff,
      upper_frequency_cutoff = upper_frequency_cutoff
    )
  }

  return(object)
})

#' @rdname MiDAS-class
#'
#' @title Filter MiDAS object by omnibus groups
#'
#' @param object \code{\link{MiDAS}} object
#' @param experiment Matrix
#' @param groups Character
#'
#' @return Object of class MiDAS.
#'
#' @importFrom assertthat assert_that is.string see_if
#' @importFrom methods validObject
#'
#' @export
setGeneric(
  name = "filterByOmnibusGroups",
  def = function(object,
                 experiment,
                 groups) {
    standardGeneric("filterByOmnibusGroups")
  }
)

#' @rdname MiDAS-class
#'
#' @importFrom assertthat assert_that is.string see_if
#' @importFrom S4Vectors metadata
#'
setMethod(
  f = "filterByOmnibusGroups",
  signature = "MiDAS",
  definition = function(object,
                        experiment,
                        groups) {
    assert_that(
      is.string(experiment),
      characterMatches(experiment, getExperiments(object)),
      is.character(groups)
    )
    omnibus_groups <- getOmnibusGroups(object, experiment)
    se <- object[[experiment]]
    assert_that(
      see_if(
        !is.null(omnibus_groups),
        msg = sprintf(
          "Omnibus groups filtration does not support experiment '%s'",
          experiment
        )
      ),
      see_if(
        isClass(se, "SummarizedExperiment"),
        msg = sprintf(
          "Omnibus groups filtration does not support experiment '%s'",
          experiment
        )
      ),
      characterMatches(groups, names(omnibus_groups))
    )

    mask <- unlist(omnibus_groups[groups])
    object[[experiment]] <- se[mask, ]
    metadata(object[[experiment]])$omnibus_groups <- omnibus_groups[groups]

    return(object)
  })

#' Coerce MiDAS to Data Frame
#'
#' @rdname as.data.frame
#' @method as.data.frame MiDAS
#'
#' @inheritParams as.data.frame
#'
#' @export
#'
as.data.frame.MiDAS <- function(x, ...) {
  midasToWide(object = x,
              experiment = getExperiments(x))
}

#' Prepare MiDAS object for statistical analysis
#'
#' \code{prepareMiDAS} transform HLA alleles calls and KIR calls according
#' to selected analysis creating \code{\link{MiDAS}} object.
#'
#' @inheritParams checkHlaCallsFormat
#' @inheritParams checkKirCallsFormat
#' @inheritParams hlaCallsToCounts
#' @param colData Data frame holding additional variables like phenotypic
#'   observations or covariates. It have to contain \code{'ID'} column holding
#'   samples identifiers corresponding to identifiers in \code{hla_calls} and
#'   \code{kir_calls}. Importantly rows of \code{hla_calls} and
#'   \code{kir_calls} without corresponding phenotype are discarded.
#' @param experiment Character vector indicating analysis type for which data
#'   should be prepared. Valid choices are \code{"hla_allele"},
#'   \code{"aa_level"}, \code{"allele_g_group"}, \code{"allele_supertype"},
#'   \code{"allele_group"}, \code{"kir_genes"}, \code{"hla_kir_interactions"},
#'   \code{"hla_divergence"}.
#'   See details for further explanations.
#' @param placeholder String
#' @param lower_frequency_cutoff Number
#' @param upper_frequency_cutoff Number
#' @param ... Attributes used in data transformation.
#'
#' @return Object of class \code{\link{MiDAS}}
#'
#' @examples
#' \dontrun{
#' # read hla calls file
#' hla_calls_file <- system.file("extdata", "HLAHD_output_example.txt", package = "MiDAS")
#' hla_calls <- readHlaCalls(hla_calls_file)
#'
#' # read kir calls file
#' kir_calls_file <- system.file("extdata", "KIP_output_example.txt", package = "MiDAS")
#' kir_calls <- readKPICalls(kir_calls_file, counts = TRUE)
#'
#' # read phenotypic data and covariates
#' pheno_file <- system.file("extdata", "pheno_example.txt", package = "MiDAS")
#' pheno <- read.table(pheno_file, header = TRUE, stringsAsFactors = FALSE)
#' covar_file <- system.file("extdata", "covar_example.txt", package = "MiDAS")
#' covar <- read.table(covar_file, header = TRUE, stringsAsFactors = FALSE)
#' phenotype <- left_join(pheno, covar, by ="ID")
#'
#' # create MiDAS object
#' midas <- prepareMiDAS(hla_calls = hla_calls,
#'                       kir_calls = kir_calls,
#'                       colData = phenotype,
#'                       inheritance_model = "additive",
#'                       experiment = "hla_allele"
#' )
#' }
#'
#' @importFrom assertthat assert_that see_if
#' @importFrom MultiAssayExperiment ExperimentList MultiAssayExperiment
#' @importFrom rlang is_empty
#' @importFrom stats runif
#' @importFrom S4Vectors DataFrame
#'
#' @export
prepareMiDAS <- function(hla_calls = NULL,
                         kir_calls = NULL,
                         colData,
                         inheritance_model = c("additive", "dominant", "recessive"),
                         experiment = c(
                           "hla_allele",
                           "aa_level",
                           "allele_g_group",
                           "allele_supertype",
                           "allele_group",
                           "kir_genes",
                           "hla_kir_interactions",
                           "hla_divergence"
                         ),
                         placeholder = "term",
                         lower_frequency_cutoff = NULL,
                         upper_frequency_cutoff = NULL,
                         ...
) {
  inheritance_model_choice <- eval(formals()[["inheritance_model"]])
  experiment_choice <- eval(formals()[["experiment"]])
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
    checkColDataFormat(colData),
    is.string(inheritance_model),
    stringMatches(inheritance_model, inheritance_model_choice),
    is.character(experiment),
    characterMatches(experiment, experiment_choice),
    is.string(placeholder),
    see_if(
      !placeholder %in% c(colnames(hla_calls), colnames(kir_calls), colnames(colData)),
      msg = sprintf(
        "Placeholder '%s' can not be used, it is alredy used as column name in one of the inputs.",
        placeholder
      )
    )
  )

  dot.args <- if (...length()) {
    list(...)
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
  for (e in experiment) {
    fun <- paste0("prepareMiDAS_", e)
    args <- list(
      hla_calls = hla_calls,
      kir_calls = kir_calls,
      inheritance_model = inheritance_model
    )
    args <- c(args, dot.args)
    E <- do.call(
      what = fun,
      args = args
    )

    # frequency filter
    if ((!is.null(lower_frequency_cutoff) ||
         !is.null(upper_frequency_cutoff))
        && is.integer(E)) {
      E <- filterExperimentByFrequency(
        experiment = E,
        inheritance_model = inheritance_model,
        lower_frequency_cutoff = lower_frequency_cutoff,
        upper_frequency_cutoff = upper_frequency_cutoff
      )
    }

    experiments[[e]] <- E
  }

  # insert placeholder
  colData[[placeholder]] <- runif(nrow(colData))

  colData <-
    DataFrame(colData, row.names = colData[["ID"]], check.names = TRUE)

  metadata <- list(
    inheritance_model = inheritance_model,
    experiment = experiment,
    placeholder = placeholder
  )

  mae <- MultiAssayExperiment(
    experiments = experiments,
    colData = colData,
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
  assert_that(
    checkHlaCallsFormat(hla_calls)
  )

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
#' @inheritParams hlaToAAVariation
#' @param hla_calls Data frame
#' @param inheritance_model String
#' @param ... Not used
#'
#' @return Matrix
#'
#' @importFrom assertthat assert_that is.flag
#' @importFrom SummarizedExperiment SummarizedExperiment
#'
prepareMiDAS_aa_level <- function(hla_calls,
                                  inheritance_model,
                                  indels = TRUE,
                                  unkchar = FALSE,
                                  ...) {
  assert_that(
    checkHlaCallsFormat(hla_calls),
    is.flag(indels),
    is.flag(unkchar)
  )

  counts <-
    hlaToAAVariation(
      hla_calls = hla_calls,
      indels = indels,
      unkchar = unkchar
    ) %>%
    aaVariationToCounts(inheritance_model = inheritance_model) %>%
    dfToExperimentMat()
  vars <- counts %>%
    rownames()
  pos <- vars %>%
    gsub(pattern = "_.{1}$", replacement = "") %>%
    unique()
  omnibus_groups <- lapply(
    pos,
    FUN = function(p) {
      grep(pattern = paste0("^", p, "_"),
           x = vars,
           value = TRUE)
    }
  )
  names(omnibus_groups) <- pos

  hla_aa <- SummarizedExperiment(
    assays = counts,
    metadata = list(omnibus_groups = omnibus_groups)
  )

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
  assert_that(
    checkHlaCallsFormat(hla_calls)
  )

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
  assert_that(
    checkHlaCallsFormat(hla_calls)
  )

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
  assert_that(
    checkHlaCallsFormat(hla_calls)
  )

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
#' @importFrom assertthat assert_that
#'
prepareMiDAS_kir_genes <- function(kir_calls, ...) {
  assert_that(
    checkKirCallsFormat(kir_calls)
  )

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
#' @importFrom magrittr %>%
#'
prepareMiDAS_hla_kir_interactions <- function(hla_calls, kir_calls, ...) {
  assert_that(
    checkHlaCallsFormat(hla_calls),
    checkKirCallsFormat(hla_calls)
  )

  hla_kir_interactions <-
    getHlaKirInteractions(
      hla_calls = hla_calls,
      kir_counts = kir_calls
    ) %>%
    dfToExperimentMat()

  return(hla_kir_interactions)
}

#' Prepare MiDAS data on HLA divergence level
#'
#' Distances between Class I alleles are calculated using Grantham distance,
#' as implemented in \code{hlaCallsGranthamDistance} function. Additionally
#' average distance in Class I genese is calculated.
#'
#' @param hla_calls Data frame
#' @param ... Not used
#'
#' @return Matrix
#'
prepareMiDAS_hla_divergence <- function(hla_calls, ...) {
  genes <- getHlaCallsGenes(hla_calls)
  assert_that(
    any(c("A", "B", "C") %in% genes),
    msg = "Grantham distance can be calculated only for class I HLA alleles (A, B, C)."
  )

  hla_divergence <-
    hlaCallsGranthamDistance(
      hla_calls = hla_calls,
      genes = genes[genes %in% c("A", "B", "C")]
    )
  hla_divergence$ABC_avg <- rowMeans(hla_divergence[-1])
  experiment_mat <- dfToExperimentMat(hla_divergence)

  return(experiment_mat)
}
