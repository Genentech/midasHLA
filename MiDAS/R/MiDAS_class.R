#' @include MiDAS_class.R
NULL

#' MiDAS class
#'
#' @description
#' The \code{MiDAS} class is a \code{\link{MultiAssayExperiment}} object
#' containing data and metadata required for MiDAS analysis.
#'
#' Valid \code{MiDAS} object must have unique features names across all
#' experiments and colData. It's metadata list need to have \code{placeholder}
#' element, which is a string specifying name of column in colData used when
#' defining statistical model for downstream analyses (see
#' \code{\link{runMiDAS}} for more details). Optionally the object's metadata
#' can also store \code{'hla_calls'} and \code{'kir_calls'} data frames (see
#' \code{\link{prepareMiDAS}} for more details).
#'
#' @importFrom methods setClass new
#' @importClassesFrom MultiAssayExperiment MultiAssayExperiment
#'
#' @export
MiDAS <- setClass(
  Class = "MiDAS",
  contains = "MultiAssayExperiment"
)

# MiDAS class initialize method
#
# @inheritParams methods::initialize
# @inheritParams MultiAssayExperiment
#
# @return MiDAS object
#
# @importFrom MultiAssayExperiment MultiAssayExperiment
#
setMethod("initialize", "MiDAS", function(.Object, experiments, colData, metadata) {
  midas <- MultiAssayExperiment(
    experiments = experiments,
    colData = colData,
    metadata = metadata
  )
  class(midas) <- structure("MiDAS", package = "MiDAS")

  return(midas)
})

# Validity method for class MiDAS
#
# A valid \code{MiDAS} object must have unique features names across all
# experiments and colData. It's metadata list need to have "placeholder"
# element, which is a string specifying name of column in colData that should
# be used in statistical model definitions for downstream analyses. Optionally
# the object's metadata can also store \code{'hla_calls'} and \code{'kir_calls'}
# data frames. See \code{\link{prepareMiDAS}} for more details.
#
# @importFrom assertthat assert_that is.string see_if
#
setValidity(Class = "MiDAS", method = function(object) {
  hla_calls <- getHlaCalls(object)
  kir_calls <- getKirCalls(object)
  placeholder <- getPlaceholder(object)

  assert_that(
    if (! is.null(hla_calls)) { checkHlaCallsFormat(hla_calls) } else { TRUE },
    if (! is.null(kir_calls)) { checkKirCallsFormat(kir_calls) } else { TRUE },
    is.string(placeholder),
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

#' Get available experiments in MiDAS object.
#'
#' @param object \code{\link{MiDAS}} object.
#'
#' @return Character vector giving names of experiments in \code{object}.
#'
#' @importFrom S4Vectors metadata
#' @export
setGeneric(
  name = "getExperiments",
  def = function(object) standardGeneric("getExperiments")
)

#' @rdname MiDAS-class
#'
#' @title Get available experiments in MiDAS object.
#'
#' @export
setMethod(
  f = "getExperiments",
  signature = "MiDAS",
  definition = function (object) metadata(object)$experiment
)

#' Get HLA calls from MiDAS object.
#'
#' @param object \code{\link{MiDAS}} object.
#'
#' @return HLA calls data frame.
#'
#' @importFrom S4Vectors metadata
#' @export
setGeneric(
  name = "getHlaCalls",
  def = function(object) standardGeneric("getHlaCalls")
)

#' @rdname MiDAS-class
#'
#' @title Get HLA calls from MiDAS object.
#'
#' @importFrom S4Vectors metadata
#' @export
setMethod(
  f = "getHlaCalls",
  signature = "MiDAS",
  definition = function (object) { metadata(object)[["hla_calls"]] }
)

#' Get KIR calls from MiDAS object.
#'
#' @param object \code{\link{MiDAS}} object.
#'
#' @return KIR calls data frame.
#'
#' @importFrom S4Vectors metadata
#' @export
setGeneric(
  name = "getKirCalls",
  def = function(object) standardGeneric("getKirCalls")
)

#' @rdname MiDAS-class
#'
#' @title Get KIR calls from MiDAS object.
#'
#' @importFrom S4Vectors metadata
#' @export
setMethod(
  f = "getKirCalls",
  signature = "MiDAS",
  definition = function (object) { metadata(object)[["kir_calls"]] }
)

#' Get placeholder name from MiDAS object.
#'
#' @param object \code{\link{MiDAS}} object.
#'
#' @return String giving name of placeholder.
#'
#' @importFrom S4Vectors metadata
#' @export
setGeneric(
  name = "getPlaceholder",
  def = function(object) standardGeneric("getPlaceholder")
)

#' @rdname MiDAS-class
#'
#' @title Get placeholder name from MiDAS object.
#'
#' @importFrom S4Vectors metadata
#' @export
setMethod(
  f = "getPlaceholder",
  signature = "MiDAS",
  definition = function (object) {
    placeholder <- metadata(object)[["placeholder"]]

    return(placeholder)
  }
)

#' Get omnibus groups from MiDAS object.
#'
#' @details For some experiments features can be naturally divided into groups
#' (here called omnibus groups). For example, in \code{'hla_aa'} experiment
#' features can be grouped by amino acid position (\code{"B_46_E"},
#' \code{"B_46_A"}) can be grouped into \code{B_46} group). Such groups can be
#' than used to perform omnibus test, see \code{\link{runMiDAS}} for more
#' details.
#'
#' @param object \code{\link{MiDAS}} object.
#' @param experiment String specifying experiment.
#'
#' @return List of omnibus groups for given experiment.
#'
#' @importFrom assertthat assert_that is.string
#' @importFrom S4Vectors metadata
#' @export
setGeneric(
  name = "getOmnibusGroups",
  def = function(object, experiment) standardGeneric("getOmnibusGroups")
)

#' Get omnibus groups from MiDAS object.
#'
#' @details For some experiments features can be naturally divided into groups
#' (here called omnibus groups). For example, in \code{'hla_aa'} experiment
#' features can be grouped by amino acid position (\code{"B_46_E"},
#' \code{"B_46_A"}) can be grouped into \code{B_46} group). Such groups can be
#' than used to perform omnibus test, see \code{\link{runMiDAS}} for more
#' details.
#'
#' @param object \code{\link{MiDAS}} object.
#' @param experiment String specifying experiment.
#'
#' @importFrom assertthat assert_that is.string
#' @importFrom S4Vectors metadata
#' @export
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

#' Calculate features frequencies for given experiment in MiDAS object.
#'
#' @inheritParams getExperimentFrequencies
#' @param object \code{\link{MiDAS}} object.
#' @param compare Logical flag indicating if \code{hla_calls} frequencies
#'   should be compared to reference frequencies given in \code{ref}.
#' @param ref_pop Character vector giving names of reference populations in
#'   \code{ref} to compare with. Optionally vector can be named, then those
#'   names will be used as population names.
#' @param ref Named list of reference frequencies data frames. See
#'   \code{\link{allele_frequencies}} for an example on how reference
#'   frequency data frame should be formatted.
#'
#' @return Data frame with features from selected experiment and their
#'   corresponding frequencies.
#'
#' @importFrom assertthat assert_that is.string
#' @importFrom dplyr select
#' @importFrom rlang !!
#' @importFrom stats reshape
#' @importFrom S4Vectors metadata
#' @export
setGeneric(
  name = "getFrequencies",
  def = function(object,
                 experiment,
                 carrier_frequency = FALSE,
                 compare = FALSE,
                 ref_pop = c("USA NMDP African American pop 2", "USA NMDP Chinese", "USA NMDP European Caucasian"),
                 ref = list(hla_alleles = allele_frequencies)
  ) {
    standardGeneric("getFrequencies")
  }
)

#' @rdname MiDAS-class
#'
#' @title Calculate features frequencies for given experiment in MiDAS object.
#'
#' @inheritParams getExperimentFrequencies
#' @param compare Logical flag indicating if \code{hla_calls} frequencies
#'   should be compared to reference frequencies given in \code{ref}.
#' @param ref_pop Character vector giving names of reference populations in
#'   \code{ref} to compare with. Optionally vector can be named, then those
#'   names will be used as population names.
#' @param ref Named list of reference frequencies data frames. See
#'   \code{\link{allele_frequencies}} for an example on how reference
#'   frequency data frame should be formatted.
#'
#' @importFrom assertthat assert_that is.string
#' @importFrom dplyr select
#' @importFrom rlang !!
#' @importFrom stats reshape
#' @importFrom S4Vectors metadata
#' @export
setMethod(
  f = "getFrequencies",
  signature = "MiDAS",
  definition = function (object,
                         experiment,
                         carrier_frequency = FALSE,
                         compare = FALSE,
                         ref_pop = c("USA NMDP African American pop 2", "USA NMDP Chinese", "USA NMDP European Caucasian"),
                         ref = list(hla_alleles = allele_frequencies)) {
    assert_that(
      is.string(experiment),
      stringMatches(experiment, getExperiments(object)),
      isTRUEorFALSE(compare),
      is.character(ref_pop),
      is.list(ref)
    )
    mat <- object[[experiment]]
    assert_that(
      is.matrix(mat) && isCountsOrZeros(mat),
      msg = sprintf("Frequencies can not be calculated for experiment '%s'",
                    experiment
            )
    )

    ref <- ref[[experiment]]
    if (compare && ! is.null(ref)) {
      assert_that(
        length(ref_pop) > 0,
        msg = "Please specify reference populations using 'ref_pop' argument."
      )
      ref <- getReferenceFrequencies(ref, ref_pop, carrier_frequency)
      freq <- getExperimentFrequencies(mat, carrier_frequency, ref)
    } else {
      if (compare) warn(sprintf("Could not find reference frequencies for experiment: '%s'", experiment))
      freq <- getExperimentFrequencies(mat, carrier_frequency)
    }

    return(freq)
  }
)

#' Filter MiDAS object by frequency
#'
#' @inheritParams getExperimentFrequencies
#' @param object \code{\link{MiDAS}} object.
#' @param experiment String specifying experiment.
#' @param lower_frequency_cutoff Number giving lower frequency threshold.
#'   Numbers greater than 1 are interpreted as number of feature occurrences,
#'   numbers between 0 and 1 as fractions.
#' @param upper_frequency_cutoff Number giving upper frequency threshold.
#'   Numbers greater than 1 are interpreted as number of feature occurrences,
#'   numbers between 0 and 1 as fractions.
#'
#' @return Filtered \code{\link{MiDAS}} object.
#'
#' @importFrom assertthat assert_that is.string
#' @importFrom S4Vectors metadata
#' @export
setGeneric(
  name = "filterByFrequency",
  def = function(object,
                 experiment,
                 lower_frequency_cutoff = NULL,
                 upper_frequency_cutoff = NULL,
                 carrier_frequency = FALSE) {
    standardGeneric("filterByFrequency")
  }
)

#' @rdname MiDAS-class
#'
#' @title Filter MiDAS object by frequency
#'
#' @param lower_frequency_cutoff Number giving lower frequency threshold.
#'   Numbers greater than 1 are interpreted as number of feature occurrences,
#'   numbers between 0 and 1 as fractions.
#' @param upper_frequency_cutoff Number giving upper frequency threshold.
#'   Numbers greater than 1 are interpreted as number of feature occurrences,
#'   numbers between 0 and 1 as fractions.
#'
#' @importFrom assertthat assert_that is.string
#' @importFrom S4Vectors metadata
#' @export
setMethod(
  f = "filterByFrequency",
  signature = "MiDAS",
  definition = function(object,
                        experiment,
                        lower_frequency_cutoff = NULL,
                        upper_frequency_cutoff = NULL,
                        carrier_frequency = FALSE) {
  assert_that(
    is.character(experiment),
    characterMatches(experiment, getExperiments(object)),
    validateFrequencyCutoffs(lower_frequency_cutoff, upper_frequency_cutoff)
  )

  for (ex in experiment) {
    mat <- object[[ex]]
    assert_that(
      is.matrix(mat) && isCountsOrZeros(mat),
      msg = sprintf("Frequency filtration does not support experiment '%s'", ex)
    )
    object[[ex]] <- filterExperimentByFrequency(
      experiment = mat,
      carrier_frequency = carrier_frequency,
      lower_frequency_cutoff = lower_frequency_cutoff,
      upper_frequency_cutoff = upper_frequency_cutoff
    )
  }

  return(object)
})


#' Filter MiDAS object by omnibus groups
#'
#' @param object \code{\link{MiDAS}} object.
#' @param experiment String specifying experiment.
#' @param groups Character vector specifying omnibus groups to select. See
#'   \code{\link{getOmnibusGroups}} for more details.
#'
#' @return Filtered \code{\link{MiDAS}} object.
#'
#' @importFrom assertthat assert_that is.string see_if
#' @importFrom S4Vectors metadata
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
#' @title Filter MiDAS object by omnibus groups
#'
#' @param object \code{\link{MiDAS}} object.
#' @param experiment String specifying experiment.
#' @param groups Character vector specifying omnibus groups to select. See
#'   \code{\link{getOmnibusGroups}} for more details.
#'
#' @importFrom assertthat assert_that is.string see_if
#' @importFrom S4Vectors metadata
#' @export
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

    mask <- unlist(omnibus_groups[groups]) # TODO what if some groups are missing due to smth
    object[[experiment]] <- se[mask, ]
    S4Vectors::metadata(object[[experiment]])$omnibus_groups <- omnibus_groups[groups]

    return(object)
  })


#' Filter MiDAS object by features
#'
#' @param object \code{\link{MiDAS}} object.
#' @param experiment String specifying experiment.
#' @param variables Character vector specifying features to select.
#'
#' @return Filtered \code{\link{MiDAS}} object.
#'
#' @export
setGeneric(
  name = "filterByVariables",
  def = function(object, experiment, variables) standardGeneric("filterByVariables")
)

#' @rdname MiDAS-class
#'
#' @title Filter MiDAS object by features
#'
#' @param variables Character vector specifying features to select.
#'
#' @importFrom S4Vectors metadata
#' @export
setMethod(
  f = "filterByVariables",
  signature = "MiDAS",
  definition = function (object, experiment, variables) {
    assert_that(
      validObject(object),
      is.string(experiment),
      stringMatches(experiment, getExperiments(object)),
      is.character(variables),
      characterMatches(variables, rownames(object))
    )
    object[[experiment]] <- filterExperimentByVariables(object[[experiment]], variables)

    return(object)
  }
)

#' Coerce MiDAS to Data Frame
#'
#' @method as.data.frame MiDAS
#'
#' @inheritParams base::as.data.frame
#' @export
as.data.frame.MiDAS <- function(x, ...) {
  midasToWide(object = x,
              experiment = getExperiments(x))
}

#' Construct a MiDAS object
#'
#' \code{prepareMiDAS} transform HLA alleles calls and KIR calls according
#' to selected \code{experiment}s creating a \code{\link{MiDAS}} object.
#'
#' \code{experiment} specifies analysis types for which \code{hla_calls} and
#' \code{kir_call} should be prepared.
#' \describe{
#'   \item{\code{'hla_alleles'}}{
#'     \code{hla_calls} are transformed to counts matrix describing number of
#'     allele occurrences for each sample. This experiment is used to test
#'     associations on HLA alleles level.
#'   }
#'   \item{\code{'hla_aa'}}{
#'     \code{hla_calls} are transformed to a matrix of variable amino acid
#'     positions. See \code{\link{hlaToAAVariation}} for more details. This
#'     experiment is used to test associations on amino acid level.
#'   }
#'   \item{\code{"hla_g_groups"}}{
#'     \code{hla_calls} are translated into HLA G groups and transformed to
#'     matrix describing number of G group occurrences for each sample. See
#'     \code{\link{hlaToVariable}} for more details. This experiment is used to
#'     test associations on HLA G groups level.
#'   }
#'   \item{\code{"hla_supertypes"}}{
#'     \code{hla_calls} are translated into HLA supertypes and transformed to
#'     matrix describing number of G group occurrences for each sample. See
#'     \code{\link{hlaToVariable}} for more details. This experiment is used to
#'     test associations on HLA supertypes level.
#'   }
#'   \item{\code{"hla_NK_ligands"}}{
#'     \code{hla_calls} are translated into NK ligands, which includes HLA
#'     Bw4/Bw6 and HLA C1/C2 groups and transformed to matrix describing number
#'     of their occurrences for each sample. See \code{\link{hlaToVariable}} for
#'     more details.This experiment is used to test associations on HLA NK
#'     ligands level.
#'   }
#'   \item{\code{"kir_genes"}}{
#'     \code{kir_calls} are transformed to counts matrix describing number of
#'     KIR gene occurrences for each sample. This experiment is used to test
#'     associations on KIR genes level.
#'   }
#'   \item{\code{"hla_kir_interactions"}}{
#'     \code{hla_calls} and \code{kir_calls} are translated to HLA - KIR
#'     interactions as defined in
#'     \href{https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6558367/}{Pende et al., 2019.}.
#'     See \code{\link{getHlaKirInteractions}} for more details. This experiment
#'     is used to test associations on HLA - KIR interactions level.
#'   }
#'   \item{\code{"hla_divergence"}}{
#'     Grantham distance for class I HLA alleles is calculated based on
#'     \code{hla_calls} using original formula by
#'     \href{http://www.sciencemag.org/content/185/4154/862.long}{Grantham R. 1974.}.
#'     See \code{\link{hlaCallsGranthamDistance}} for more details. This
#'     experiment is used to test associations on HLA divergence level measured
#'     by Grantham distance.
#'   }
#'   \item{\code{"hla_het"}}{
#'     \code{hla_calls} are transformed to heterozygosity status, where \code{1}
#'     designates a heterozygote and \code{0} homozygote. Heterozygosity status
#'     is calculated only for classical HLA genes (A, B, C, DQA1, DQB1, DRA,
#'     DRB1, DPA1, DPB1). This experiment is used to test associations on HLA
#'     divergence level measured by heterozygosity.
#'   }
#' }
#'
#' @inheritParams checkHlaCallsFormat
#' @inheritParams checkKirCallsFormat
#' @inheritParams hlaCallsToCounts
#' @inheritParams filterByFrequency
#' @inheritParams prepareMiDAS_hla_het
#' @param colData Data frame holding additional variables like phenotypic
#'   observations or covariates. It have to contain \code{'ID'} column holding
#'   samples identifiers corresponding to identifiers in \code{hla_calls} and
#'   \code{kir_calls}. Importantly rows of \code{hla_calls} and
#'   \code{kir_calls} without corresponding phenotype are discarded.
#' @param experiment Character vector indicating analysis type for which data
#'   should be prepared. Valid choices are \code{"hla_alleles"},
#'   \code{"hla_aa"}, \code{"hla_g_groups"}, \code{"hla_supertypes"},
#'   \code{"hla_NK_ligands"}, \code{"kir_genes"}, \code{"hla_kir_interactions"},
#'   \code{"hla_divergence"}, \code{"hla_het"}.
#'   See details for further explanations.
#' @param placeholder String giving name for dummy variable inserted to
#'   \code{colData}. This variable can be than used to define base statistical
#'   model used by \code{\link{runMiDAS}}.
#' @param indels Logical indicating whether indels should be considered when
#'   checking amino acid variability in \code{'hla_aa'} experiment.
#' @param unkchar Logical indicating whether unknown characters in the alignment
#'   should be considered when checking amino acid variability in
#'   \code{'hla_aa'} experiment.
#' @param hla_het_resolution Number specifying HLA alleles resolution used to
#'   calculate heterogeneity in \code{"hla_het"} experiment.
#'
#' @return Object of class \code{\link{MiDAS}}
#'
#' @examples
#' midas <- prepareMiDAS(hla_calls = MiDAS_tut_HLA,
#'                       kir_calls = MiDAS_tut_KIR,
#'                       colData = MiDAS_tut_pheno,
#'                       experiment = "hla_alleles")
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
                         experiment = c(
                           "hla_alleles",
                           "hla_aa",
                           "hla_g_groups",
                           "hla_supertypes",
                           "hla_NK_ligands",
                           "kir_genes",
                           "hla_kir_interactions",
                           "hla_divergence",
                           "hla_het"
                         ),
                         placeholder = "term",
                         lower_frequency_cutoff = NULL,
                         upper_frequency_cutoff = NULL,
                         ...
) {
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
      kir_calls = kir_calls
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
        lower_frequency_cutoff = lower_frequency_cutoff,
        upper_frequency_cutoff = upper_frequency_cutoff
      )
    }

    experiments[[e]] <- E
  }

  # insert placeholder
  colData[[placeholder]] <- runif(nrow(colData))

  colData <- DataFrame(colData[, ! colnames(colData) == "ID"],
                       row.names = colData[["ID"]],
                       check.names = TRUE)

  new("MiDAS", experiments = experiments, colData = colData, metadata = metadata)
}


#' Prepare MiDAS data on HLA allele level
#'
#' @inheritParams prepareMiDAS
#' @param ... Not used
#'
#' @return Matrix
#'
prepareMiDAS_hla_alleles <- function(hla_calls, ...) {
  assert_that(
    checkHlaCallsFormat(hla_calls)
  )

  hla_alleles <- hlaCallsToCounts(hla_calls = hla_calls) %>%
    dfToExperimentMat()

  return(hla_alleles)
}

#' Prepare MiDAS data on HLA amino acid level
#'
#' @inheritParams checkHlaCallsFormat
#' @inheritParams hlaToAAVariation
#' @param ... Not used
#'
#' @return Matrix
#'
#' @importFrom assertthat assert_that is.flag
#' @importFrom SummarizedExperiment SummarizedExperiment
#'
prepareMiDAS_hla_aa <- function(hla_calls,
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
    aaVariationToCounts() %>%
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
#' @param ... Not used
#'
#' @return Matrix
#'
#' @importFrom assertthat assert_that
#'
prepareMiDAS_hla_g_groups <- function(hla_calls,
                                      ...) {
  assert_that(
    checkHlaCallsFormat(hla_calls)
  )

  lib <- "allele_HLA_Ggroup"
  hla_g_groups <- hlaToVariable(hla_calls = hla_calls,
                                dictionary = lib,
                                na.value = 0
  )

  assert_that(
    ncol(hla_g_groups) > 1,
    msg = "no allele could be assigned to G group for input hla_calls"
  )

  hla_g_groups <-
    hlaCallsToCounts(
      hla_calls = hla_g_groups,
      check_hla_format = FALSE
    ) %>%
    dfToExperimentMat()

  return(hla_g_groups)
}

#' Prepare MiDAS data on HLA allele's supertypes level
#'
#' @inheritParams checkHlaCallsFormat
#' @param ... Not used
#'
#' @return Matrix
#'
#' @importFrom assertthat assert_that
#' @importFrom dplyr select
#' @importFrom rlang .data
#'
prepareMiDAS_hla_supertypes <- function(hla_calls, ...) {
  assert_that(
    checkHlaCallsFormat(hla_calls)
  )

  lib <- "allele_HLA_supertype"
  hla_supertypes <- hlaToVariable(hla_calls = hla_calls,
                                    dictionary = lib,
                                    na.value = 0
  )

  assert_that(
    ncol(hla_supertypes) > 1,
    msg = "no allele could be assigned to supertype for input hla_calls"
  )

  hla_supertypes <-
    hlaCallsToCounts(
      hla_calls = hla_supertypes,
      check_hla_format = FALSE
    ) %>%
    subset(select = - Unclassified) %>%
    dfToExperimentMat()

  return(hla_supertypes)
}

#' Prepare MiDAS data on HLA allele's groups level
#'
#' @inheritParams checkHlaCallsFormat
#' @param ... Not used
#'
#' @return Matrix
#'
#' @importFrom assertthat assert_that
#'
prepareMiDAS_hla_NK_ligands <- function(hla_calls, ...) {
  assert_that(
    checkHlaCallsFormat(hla_calls)
  )

  lib <- c(
    "allele_HLA_Bw4",
    "allele_HLA-B_only_Bw",
    "allele_HLA-C_C1-2"
  )
  hla_NK_ligands <- Reduce(
    f = function(...) left_join(..., by = "ID"),
    x = lapply(lib, hlaToVariable, hla_calls = hla_calls, na.value = 0)
  )

  assert_that(
    ncol(hla_NK_ligands) > 1,
    msg = "no allele could be assigned to allele groups for input hla_calls"
  )

  hla_NK_ligands <-
    hlaCallsToCounts(
      hla_calls = hla_NK_ligands,
      check_hla_format = FALSE
    ) %>%
    dfToExperimentMat()

  return(hla_NK_ligands)
}

#' Prepare MiDAS data on KIR genes level
#'
#' @inheritParams checkKirCallsFormat
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
#' @inheritParams checkHlaCallsFormat
#' @inheritParams checkKirCallsFormat
#' @param ... Not used
#'
#' @return Matrix
#'
#' @importFrom magrittr %>%
#'
prepareMiDAS_hla_kir_interactions <- function(hla_calls, kir_calls, ...) {
  assert_that(
    checkHlaCallsFormat(hla_calls),
    checkKirCallsFormat(kir_calls)
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
#' @inheritParams checkHlaCallsFormat
#' @param ... Not used
#'
#' @return Matrix
#'
prepareMiDAS_hla_divergence <- function(hla_calls, ...) {
  assert_that(
    checkHlaCallsFormat(hla_calls)
  )
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

#' Prepare MiDAS data on HLA heterozygosity level
#'
#' @inheritParams checkHlaCallsFormat
#' @param hla_het_resolution Number specifying HLA alleles resolution used to
#'   calculate heterogeneity.
#' @param ... Not used
#'
#' @return Matrix
#'
prepareMiDAS_hla_het <- function(hla_calls, hla_het_resolution = 8, ...) {
  assert_that(
    checkHlaCallsFormat(hla_calls),
    is.number(hla_het_resolution) # TODO see how this check was done elsewhere
  )

  genes <- getHlaCallsGenes(hla_calls)
  classical_genes <- c("A", "B", "C", "DQA1", "DQB1", "DRA", "DRB1", "DPA1", "DPB1")
  assert_that(
    any(classical_genes %in% genes),
    msg = "Heterozygosity status can be calculated only for classical genes (A, B, C, DQA1, DQB1, DRA, DRB1, DPA1, DPB1)."
  )

  # filter non-classical genes
  sel <- c("ID", paste0(rep(classical_genes, each = 2), "_", 1:2))
  hla_calls <- hla_calls[, colnames(hla_calls) %in% sel, drop = FALSE]
  genes <- getHlaCallsGenes(hla_calls)

  # resolution
  hla_calls <- reduceHlaCalls(hla_calls, hla_het_resolution)

  hla_het <- hla_calls[, "ID", drop = FALSE]
  for (g in genes) {
    i <- paste0(g, "_1")
    j <- paste0(g, "_2")
    het <- hla_calls[, i, drop = TRUE] != hla_calls[, j, drop = TRUE]
    nm <- paste0(g, "_het")
    hla_het[[nm]] <- as.integer(het)
  }
  hla_het <- dfToExperimentMat(hla_het)

  return(hla_het)
}
