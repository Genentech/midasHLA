#' @include MiDAS_class.R
NULL

#' MiDAS class
#'
#' @description
#' The \code{MiDAS} class is a \code{\link{MultiAssayExperiment}} object
#' containing data and metadata required for MiDAS analysis.
#'
#' Valid \code{MiDAS} object must have unique features names across all
#' experiments and colData. It's metadata list needs to have a \code{placeholder}
#' element, which is a string specifying name of column in colData used when
#' defining statistical model for downstream analyses (see
#' \code{\link{runMiDAS}} for more details). Optionally the object's metadata
#' can also store \code{'hla_calls'} and \code{'kir_calls'} data frames (see
#' \code{\link{prepareMiDAS}} for more details).
#'
#' @importFrom methods setClass new
#' @importClassesFrom MultiAssayExperiment MultiAssayExperiment
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
  variables <- c(unlist(rownames(object)), colnames(colData(object)))

  assert_that(
    if (! is.null(hla_calls)) { checkHlaCallsFormat(hla_calls) } else { TRUE },
    if (! is.null(kir_calls)) { checkKirCallsFormat(kir_calls) } else { TRUE },
    is.string(placeholder),
    see_if(
      ! getPlaceholder(object) %in% unlist(rownames(object)),
      msg = sprintf("Placeholder '%s' is used in one of object's experiments",
                    getPlaceholder(object))
    ),
    see_if(
      getPlaceholder(object) %in% colnames(colData(object)),
      msg = sprintf("Placeholder '%s' can not be found in object's colData",
                    getPlaceholder(object))
    ),
    see_if(
      getPlaceholder(object) %in% colnames(colData(object)),
      msg = sprintf("Placeholder '%s' can not be found in object's colData",
                    getPlaceholder(object))
    ),
    see_if(
      ! any(duplicated(variables)),
      msg = sprintf("Object contain duplicated features: %s",
                    paste(variables[duplicated(variables)], collapse = ", "))
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
  definition = function (object) names(object)
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
#' then used to perform omnibus test, see \code{\link{runMiDAS}} for more
#' details.
#'
#' @param object \code{\link{MiDAS}} object.
#' @param experiment String specifying experiment.
#'
#' @return List of omnibus groups for a given experiment.
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
#' then used to perform omnibus test, see \code{\link{runMiDAS}} for more
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
      og <- metadata(experiment)$omnibus_groups
      if (is.null(og)) {
        NULL
      } else {
        el <- rownames(experiment)
        filterListByElements(og, el)
      }
    } else {
      NULL
    }

    return(omnibus_groups)
  }
)

#' Calculate features frequencies for a given experiment in MiDAS object.
#'
#' @inheritParams getExperimentFrequencies
#' @param object \code{\link{MiDAS}} object.
#' @param compare Logical flag indicating if \code{hla_calls} frequencies
#'   should be compared to reference frequencies given in \code{ref}.
#' @param ref_pop Named list of character vectors giving names of reference
#'   populations in \code{ref} to compare with. Optionally vectors can be named,
#'   then those names will be used as population names. Each vector should
#'   correspond to a specific experiment.
#' @param ref Named list of reference frequencies data frames. Each element
#'   should give reference for a specific experiment. See
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
                 ref_pop = list(
                   hla_alleles = c(
                     "USA NMDP African American pop 2",
                     "USA NMDP Chinese",
                     "USA NMDP European Caucasian",
                     "USA NMDP Hispanic South or Central American",
                     "USA NMDP Japanese",
                     "USA NMDP North American Amerindian",
                     "USA NMDP South Asian Indian"
                   ),
                   kir_genes = c(
                     "USA California African American KIR",
                     "USA California Asian American KIR",
                     "USA California Caucasians KIR",
                     "USA California Hispanic KIR"
                   )
                 ),
                 ref = list(hla_alleles = allele_frequencies,
                            kir_genes = kir_frequencies)
  ) {
    standardGeneric("getFrequencies")
  }
)

#' @rdname MiDAS-class
#'
#' @title Calculate features frequencies for a given experiment in MiDAS object.
#'
#' @inheritParams getExperimentFrequencies
#' @param compare Logical flag indicating if \code{hla_calls} frequencies
#'   should be compared to reference frequencies given in \code{ref}.
#' @param ref_pop Named list of character vectors giving names of reference
#'   populations in \code{ref} to compare with. Optionally vectors can be named,
#'   then those names will be used as population names. Each vector should
#'   correspond to a specific experiment.
#' @param ref Named list of reference frequencies data frames. Each element
#'   should give reference for a specific experiment. See
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
                         ref_pop = list(
                           hla_alleles = c(
                             "USA NMDP African American pop 2",
                             "USA NMDP Chinese",
                             "USA NMDP European Caucasian",
                             "USA NMDP Hispanic South or Central American",
                             "USA NMDP Japanese",
                             "USA NMDP North American Amerindian",
                             "USA NMDP South Asian Indian"
                           ),
                           kir_genes = c(
                             "USA California African American KIR",
                             "USA California Asian American KIR",
                             "USA California Caucasians KIR",
                             "USA California Hispanic KIR"
                           )
                         ),
                         ref = list(hla_alleles = allele_frequencies,
                                    kir_genes = kir_frequencies)
  ) {
    assert_that(
      is.string(experiment),
      stringMatches(experiment, getExperiments(object)),
      isTRUEorFALSE(carrier_frequency),
      isTRUEorFALSE(compare),
      is.list(ref_pop),
      is.list(ref)
    )
    ex <- object[[experiment]]
    assert_that(
      ! is.null(getExperimentPopulationMultiplicator(ex)), # if pop_mul is not set frequency cannot be calculated
      msg = sprintf("Frequencies can not be calculated for experiment '%s'",
                    experiment
            )
    )

    ref_pop <- ref_pop[[experiment]]
    ref <- ref[[experiment]]
    if (compare && ! is.null(ref)) {
      assert_that(
        length(ref_pop) > 0,
        msg = "Please specify reference populations using 'ref_pop' argument."
      )
      ref <- getReferenceFrequencies(ref, ref_pop, carrier_frequency)
      freq <- getExperimentFrequencies(experiment = ex,
                                       carrier_frequency = carrier_frequency,
                                       ref = ref)
    } else {
      if (compare) warn(sprintf("Could not find reference frequencies for experiment: '%s'", experiment))
      freq <- getExperimentFrequencies(experiment = ex,
                                       carrier_frequency = carrier_frequency)
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
#'   Numbers greater than 1 are interpreted as the number of feature occurrences,
#'   numbers between 0 and 1 as fractions.
#' @param upper_frequency_cutoff Number giving upper frequency threshold.
#'   Numbers greater than 1 are interpreted as the number of feature occurrences,
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
#'   Numbers greater than 1 are interpreted as the number of feature occurrences,
#'   numbers between 0 and 1 as fractions.
#' @param upper_frequency_cutoff Number giving upper frequency threshold.
#'   Numbers greater than 1 are interpreted as the number of feature occurrences,
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
      ! is.null(getExperimentPopulationMultiplicator(mat)), # if pop_mul is not set frequency cannot be calculated
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
      characterMatches(variables, unlist(rownames(object)))
    )
    object[[experiment]] <- filterExperimentByVariables(object[[experiment]], variables)

    return(object)
  }
)

#' Get HLA alleles for amino acid position
#'
#' List HLA alleles and amino acid residues at a given position.
#'
#' @inheritParams summariseAAPosition
#' @param object \code{\link{MiDAS}} object.
#'
#' @return Data frame containing HLA alleles, their corresponding amino acid
#'   residues and frequencies at requested position.
#'
#' @export
setGeneric(
  name = "getAllelesForAA",
  def = function(object, aa_pos) standardGeneric("getAllelesForAA")
)

#' @rdname MiDAS-class
#'
#' @title Get HLA alleles for amino acid position
#'
#' @inheritParams summariseAAPosition
#'
#' @importFrom S4Vectors metadata
#' @export
setMethod(
  f = "getAllelesForAA",
  signature = "MiDAS",
  definition = function (object, aa_pos) {
    assert_that(
      validObject(object),
      is.string(aa_pos),
      see_if(grepl("^[A-Z]+[0-9]*_-*[0-9]+$", aa_pos),
             msg = "amino acid position should be formatted like: A_9."
      )
    )

    hla_calls <- getHlaCalls(object)
    assert_that(
      see_if(! is.null(hla_calls),
             msg = "Could not find HLA calls associated with MiDAS object. Make sure to use prepareMiDAS for MiDAS object creation."
      )
    )
    if (! is.null(object[["hla_alleles"]])) { # filter hla calls to match hla_alleles experiment
      alleles <- rownames(object[["hla_alleles"]])
      hla_calls[, -1] <- rapply(
        object = hla_calls[, -1],
        f = function(x) ifelse(x %in% alleles, x, as.character(NA)),
        how = "replace"
      )
    }

    alleles_for_aa <- summariseAAPosition(hla_calls, aa_pos)

    return(alleles_for_aa)
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
#' @param hla_divergence_aa_selection String specifying variable region in peptide binding
#'   groove which should be considered for Grantham distance calculation. Valid
#'   choices includes: \code{"binding_groove"}, \code{"B_pocket"},
#'   \code{"F_pocket"}. See details for more information.
#' @param hla_het_resolution Number specifying HLA alleles resolution used to
#'   calculate heterogeneity in \code{"hla_het"} experiment.
#' @param hla_dictionary Data frame giving HLA allele dictionary used in
#'   \code{'hla_custom'} experiment. See \code{\link{hlaToVariable}} for more
#'   details.
#' @param kir_dictionary Data frame giving KIR genes dictionary used in
#'   \code{'kir_custom'} experiment. See \code{\link{countsToVariables}} for more
#'   details.
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
                           "kir_haplotypes",
                           "hla_kir_interactions",
                           "hla_divergence",
                           "hla_het",
                           "hla_custom",
                           "kir_custom"
                         ),
                         placeholder = "term",
                         lower_frequency_cutoff = NULL,
                         upper_frequency_cutoff = NULL,
                         indels = TRUE,
                         unkchar = FALSE,
                         hla_divergence_aa_selection = "binding_groove",
                         hla_het_resolution = 8,
                         hla_dictionary = NULL,
                         kir_dictionary = NULL
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
        "Placeholder '%s' can not be used, it is already used as column name in one of the inputs.",
        placeholder
      )
    )
  )

  metadata <- list(placeholder = placeholder)
  if (! is.null(hla_calls)) { metadata[["hla_calls"]] <- hla_calls }
  if (! is.null(kir_calls)) { metadata[["kir_calls"]] <- kir_calls }

  # prepare data for different analyses types
  experiments <- ExperimentList()
  args <- list(
    hla_calls = hla_calls,
    kir_calls = kir_calls,
    indels = indels,
    unkchar = unkchar,
    hla_divergence_aa_selection = hla_divergence_aa_selection,
    hla_het_resolution = hla_het_resolution,
    hla_dictionary = hla_dictionary,
    kir_dictionary = kir_dictionary
  )
  for (e in experiment) {
    fun <- paste0("prepareMiDAS_", e)
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

# Note
# Experiment metadata can contain following elements:
# \code{inheritance_model_applicable} - logical flag indicating if inheritance
# model can be applied to this experiment, \code{pop_mul} - population
# multiplicator used for features frequency calcualtions (for experiments like
# \code{"hla_alleles"} population size have to take into account presence of
# two gene copies (2 * number of samples)), \code{omnibus_groups} - named list
# of chartacter vectors giving the features groupins used for omnibus test and
# omnibus groups filtering. Set them to \code{NULL} where not applicable.

#' Prepare MiDAS data on HLA allele level
#'
#' @inheritParams prepareMiDAS
#' @param ... Not used
#'
#' @return Matrix
#'
#' @importFrom assertthat assert_that is.flag
#' @importFrom SummarizedExperiment SummarizedExperiment
#'
prepareMiDAS_hla_alleles <- function(hla_calls, ...) {
  assert_that(
    checkHlaCallsFormat(hla_calls)
  )

  hla_alleles <- hlaCallsToCounts(hla_calls = hla_calls) %>%
    dfToExperimentMat()

  hla_alleles <- SummarizedExperiment(
    assays = hla_alleles,
    metadata = list(
      inheritance_model_applicable = TRUE,
      pop_mul = 2,
      omnibus_groups = NULL)
  )

  return(hla_alleles)
}

#' Prepare MiDAS data on HLA amino acid level
#'
#' @inheritParams checkHlaCallsFormat
#' @inheritParams hlaToAAVariation
#' @param ... Not used
#'
#' @return SummarizedExperiment
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
    metadata = list(
      inheritance_model_applicable = TRUE,
      pop_mul = 2,
      omnibus_groups = omnibus_groups
    )
  )

  return(hla_aa)
}

#' Prepare MiDAS data on HLA allele's G groups level
#'
#' @inheritParams checkHlaCallsFormat
#' @param ... Not used
#'
#' @return Matrix
#'
#' @importFrom assertthat assert_that
#' @importFrom SummarizedExperiment SummarizedExperiment
#'
prepareMiDAS_hla_g_groups <- function(hla_calls, ...) {
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

  hla_g_groups <- SummarizedExperiment(
    assays = hla_g_groups,
    metadata = list(
      inheritance_model_applicable = TRUE,
      pop_mul = 2,
      omnibus_groups = NULL)
  )

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
#' @importFrom SummarizedExperiment SummarizedExperiment
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
    )
  hla_supertypes <-
    hla_supertypes[, ! colnames(hla_supertypes) == "Unclassified"] %>%
    dfToExperimentMat()

  hla_supertypes <- SummarizedExperiment(
    assays = hla_supertypes,
    metadata = list(
      inheritance_model_applicable = TRUE,
      pop_mul = 2,
      omnibus_groups = NULL)
  )

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
#' @importFrom SummarizedExperiment SummarizedExperiment
#'
prepareMiDAS_hla_NK_ligands <- function(hla_calls, ...) {
  assert_that(
    checkHlaCallsFormat(hla_calls)
  )

  # Bw6 is defined both in Bw and Bw HLA-B only dictionaries, here this unambiguity is removed
  bw_with_A <-
    system.file("extdata", "Match_allele_HLA_Bw.txt", package = "MiDAS")
  bw_with_A <-
    read.table(
      file = bw_with_A,
      header = TRUE,
      sep = "\t",
      stringsAsFactors = FALSE
    )
  mask <- bw_with_A$group != "Bw6"
  bw_with_A <- bw_with_A[mask, ]
  
  lib <- list(bw_with_A, "allele_HLA-Bw_only_B", "allele_HLA-C_C1-2")
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

  hla_NK_ligands <- SummarizedExperiment(
    assays = hla_NK_ligands,
    metadata = list(
      inheritance_model_applicable = TRUE,
      pop_mul = 2,
      omnibus_groups = NULL)
  )

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
#' @importFrom SummarizedExperiment SummarizedExperiment
#'
prepareMiDAS_kir_genes <- function(kir_calls, ...) {
  assert_that(
    checkKirCallsFormat(kir_calls)
  )

  kir_genes <-
    kir_calls %>%
    dfToExperimentMat()

  kir_genes <- SummarizedExperiment(
    assays = kir_genes,
    metadata = list(
      inheritance_model_applicable = FALSE,
      pop_mul = 1,
      omnibus_groups = NULL)
  )

  return(kir_genes)
}

#' Prepare MiDAS data on KIR haplotypes level
#'
#' @inheritParams checkKirCallsFormat
#' @param ... Not used
#'
#' @return Matrix
#'
#' @importFrom assertthat assert_that
#'
prepareMiDAS_kir_haplotypes <- function(kir_calls, ...) {
  assert_that(
    checkKirCallsFormat(kir_calls)
  )

  kir_haplotypes <- countsToVariables(counts = kir_calls,
                                      dictionary = "kir_haplotypes")

  assert_that(
    ncol(kir_haplotypes) > 1,
    msg = "None of the KIR genes could be assigned to custom variables."
  )

  kir_haplotypes <- dfToExperimentMat(kir_haplotypes)

  kir_haplotypes <- SummarizedExperiment(
    assays = kir_haplotypes,
    metadata = list(
      inheritance_model_applicable = FALSE,
      pop_mul = 1,
      omnibus_groups = NULL)
  )

  return(kir_haplotypes)
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
#' @importFrom SummarizedExperiment SummarizedExperiment
#'
prepareMiDAS_hla_kir_interactions <- function(hla_calls, kir_calls, ...) {
  assert_that(
    checkHlaCallsFormat(hla_calls),
    checkKirCallsFormat(kir_calls)
  )

  hla_kir_interactions <-
    getHlaKirInteractions(
      hla_calls = hla_calls,
      kir_calls = kir_calls
    ) %>%
    dfToExperimentMat()

  hla_kir_interactions <- SummarizedExperiment(
    assays = hla_kir_interactions,
    metadata = list(
      inheritance_model_applicable = FALSE,
      pop_mul = 1,
      omnibus_groups = NULL)
  )

  return(hla_kir_interactions)
}

#' Prepare MiDAS data on HLA divergence level
#'
#' @inheritParams checkHlaCallsFormat
#' @param hla_divergence_aa_selection String specifying variable region in peptide binding
#'   groove which should be considered for Grantham distance calculation. Valid
#'   choices includes: \code{"binding_groove"}, \code{"B_pocket"},
#'   \code{"F_pocket"}. See details for more information.
#' @param ... Not used
#'
#' @return Matrix
#'
prepareMiDAS_hla_divergence <- function(hla_calls,
                                        hla_divergence_aa_selection = "binding_groove",
                                        ...) {
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
      genes = genes[genes %in% c("A", "B", "C")],
      aa_selection = hla_divergence_aa_selection
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
    is.count(hla_het_resolution)
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

  hla_het <- SummarizedExperiment(
    assays = dfToExperimentMat(hla_het),
    metadata = list(
      inheritance_model_applicable = FALSE,
      pop_mul = 1,
      omnibus_groups = NULL)
  )

  return(hla_het)
}

#' Prepare MiDAS data on custom HLA level
#'
#' @inheritParams checkHlaCallsFormat
#' @param hla_dictionary Data frame giving HLA allele dictionary. See
#'   \code{\link{hlaToVariable}} for more details.
#' @param ... Not used
#'
#' @return Matrix
#'
#' @importFrom assertthat assert_that
#'
prepareMiDAS_hla_custom <- function(hla_calls, hla_dictionary, ...) {
  assert_that(
    checkHlaCallsFormat(hla_calls),
    is.data.frame(hla_dictionary)
  )

  hla_custom <- hlaToVariable(hla_calls = hla_calls,
                              dictionary = hla_dictionary)

  assert_that(
    ncol(hla_custom) > 1,
    msg = "None of the HLA alleles could be assigned to custom variables."
  )

  # check variables types, character need to be summarized into counts
  if (is.character(unlist(hla_custom[, -1]))) {
    hla_custom <-
      hlaCallsToCounts(
        hla_calls = hla_custom,
        check_hla_format = FALSE
      ) %>%
      dfToExperimentMat()
    hla_custom <- SummarizedExperiment(
      assays = list(hla_custom),
      metadata = list(
        inheritance_model_applicable = TRUE,
        pop_mul = 2,
        omnibus_groups = NULL
      )
    )
  } else {
    hla_custom <- dfToExperimentMat(hla_custom)
    hla_custom <- SummarizedExperiment(
      assays = list(hla_custom),
      metadata = list(
        inheritance_model_applicable = FALSE,
        pop_mul = NULL,
        omnibus_groups = NULL
      )
    )
  }

  return(hla_custom)
}

#' Prepare MiDAS data on custom KIR level
#'
#' @inheritParams checkKirCallsFormat
#' @param kir_dictionary Data frame giving KIR genes dictionary. See
#'   \code{\link{countsToVariables}} for more details.
#' @param ... Not used
#'
#' @return Matrix
#'
#' @importFrom assertthat assert_that
#'
prepareMiDAS_kir_custom <- function(kir_calls, kir_dictionary, ...) {
  assert_that(
    checkKirCallsFormat(kir_calls),
    is.data.frame(kir_dictionary)
  )

  kir_custom <- countsToVariables(counts = kir_calls,
                                  dictionary = kir_dictionary)

  assert_that(
    ncol(kir_calls) > 1,
    msg = "None of the KIR genes could be assigned to custom variables."
  )

  kir_custom <- dfToExperimentMat(kir_custom)

  kir_custom <- SummarizedExperiment(
    assays = kir_custom,
    metadata = list(
      inheritance_model_applicable = FALSE,
      pop_mul = 1,
      omnibus_groups = NULL)
  )

  return(kir_custom)
}
