#' Generate amino acid variation matrix
#'
#' \code{hlaToAAVariation} convert HLA calls data frame to a matrix of variable
#'  amino acid positions.
#'
#' Variable amino acid positions are found by comparing elements of the
#' alignment column wise. Some of the values in alignment can be treated
#' specially using \code{indels} and \code{unkchar} arguments. Function
#' processes alignments for all HLA genes found in \code{hla_calls}.
#'
#' Variable amino acid position uses protein alignments from
#' \href{ftp://ftp.ebi.ac.uk/pub/databases/ipd/imgt/hla/alignments/}{EBI database}.
#'
#' @inheritParams checkHlaCallsFormat
#' @param indels Logical indicating whether indels should be considered when
#'   checking variability.
#' @param unkchar Logical indicating whether unknown characters in the alignment
#'   should be considered when checking variability.
#' @param as_df Logical indicating if data frame should be returned.
#'   Otherwise a matrix is returned.
#'
#' @return Matrix or data frame containing variable amino acid positions.
#'   Rownames corresponds to ID column in \code{hla_calls}, and colnames to
#'   alignment positions. If no variation is found one column matrix filled with
#'   \code{NA}'s is returned.
#'
#' @examples
#' hlaToAAVariation(MiDAS_tut_HLA)
#'
#' @importFrom assertthat assert_that see_if
#' @importFrom stringi stri_split_fixed
#' @export
hlaToAAVariation <- function(hla_calls,
                             indels = TRUE,
                             unkchar = FALSE,
                             as_df = TRUE){
  assert_that(
    checkHlaCallsFormat(hla_calls),
    isTRUEorFALSE(indels),
    isTRUEorFALSE(unkchar),
    isTRUEorFALSE(as_df)
  )
  ids <- hla_calls[, 1]
  hla_calls <- hla_calls[, -1]

  # get names of genes
  gene_names <- vapply(X = colnames(hla_calls),
                       FUN = function(x) stri_split_fixed(x, "_")[[1]][1], # colnames are of form A_1, A_2, B_1, ...
                       FUN.VALUE = character(length = 1)
  )
  gene_names_uniq <- unique(gene_names)

  # assert that alignment files are available
  available_genes <- list.files(
    path = system.file("extdata", package = "MiDAS"),
    pattern = "_prot.Rdata$"
  )
  available_genes <- vapply(
    X = stri_split_fixed(available_genes, "_prot.Rdata"),
    `[[`, 1,
    FUN.VALUE = character(length = 1)
  )
  av_genes_idx <- gene_names_uniq %in% available_genes
  assert_that(
    any(av_genes_idx),
    msg = sprintf("Alignments for genes %s are not available.",
                  paste(gene_names_uniq[! av_genes_idx], collapse = ", ")
    )
  )
  if (! all(av_genes_idx)) {
    warn(sprintf(
      "Alignments for genes %s are not available and will be omitted.",
      paste(gene_names_uniq[! av_genes_idx], collapse = ", ")
    ))
    gene_names_uniq <- gene_names_uniq[av_genes_idx]
  }

  # read alignment matrices in all resolutions
  hla_aln <- lapply(X = gene_names_uniq,
                    FUN = function(x) {
                      alns <- lapply(
                        X = c(2, 4, 6, 8),
                        FUN = function(res) {
                          readHlaAlignments(
                            gene = x,
                            resolution = res,
                            unkchar = "*"
                          )
                        }
                      )
                      aln <- do.call(rbind, alns)
                      aln <- aln[! duplicated(rownames(aln)), ]

                      return(aln)
                    }
  )

  # get aa variations for each gene
  aa_variation <- list()
  for (i in 1:length(gene_names_uniq)) {
    x_calls <- hla_calls[, gene_names == gene_names_uniq[i], drop = FALSE]

    # mark alleles w/o reference as NAs
    x_calls_unlist <- unlist(x_calls)
    ref_allele <- rownames(hla_aln[[i]])
    mask_alleles_wo_ref <- ! x_calls_unlist %in% ref_allele
    if (any(mask_alleles_wo_ref[! is.na(x_calls_unlist)], na.rm = TRUE)) {
      warn(sprintf(
        "Alignments for alleles %s are not available and will be omitted.",
        paste(x_calls_unlist[mask_alleles_wo_ref], collapse = ", ")
      ))
      x_calls_unlist[mask_alleles_wo_ref] <- NA
    }

    # check if there is possibility for variability
    x_calls_uniq <- na.omit(unique(x_calls_unlist))
    if (length(x_calls_uniq) <= 1) next()

    # get variable aa positions
    hla_aln[[i]] <- hla_aln[[i]][x_calls_uniq, ]
    var_pos <- getVariableAAPos(hla_aln[[i]],
                                varchar = sprintf("[A-Z%s%s]",
                                                  ifelse(indels, "\\.", ""),
                                                  ifelse(unkchar, "\\*", "")
                                )
    )
    var_aln <- lapply(colnames(x_calls), function(allele) {
      mask <- 1:nrow(hla_aln[[i]]) # NAs in character index gives oob error, so it is needed to refer to indexes
      names(mask) <- rownames(hla_aln[[i]])
      x <- hla_aln[[i]][mask[x_calls[, allele]], var_pos, drop = FALSE]
      colnames(x) <- paste0(allele, "_", "AA_", colnames(x))
      return(x)
    })
    var_aln <- do.call(cbind, var_aln)
    ord <- as.vector(vapply(1:length(var_pos),
                  function(j) {
                    c(j, j + length(var_pos))
                  },
                  FUN.VALUE = numeric(length = 2)
    ))
    var_aln <- var_aln[, ord, drop = FALSE]

    aa_variation[[length(aa_variation) + 1]] <- var_aln
  }

  if (length(aa_variation) > 1) {
    aa_variation <- do.call(cbind, aa_variation)
    rownames(aa_variation) <- ids
  } else if (length(aa_variation) == 1) {
    aa_variation <- aa_variation[[1]]
    rownames(aa_variation) <- ids
  } else {
    aa_variation <- matrix(nrow = length(ids))
    rownames(aa_variation) <- ids
  }

  if (as_df) {
    aa_variation <- as.data.frame(aa_variation,
                                  optional = TRUE,
                                  stringsAsFactors = FALSE
    )
    aa_variation <- cbind(ID = ids, aa_variation, stringsAsFactors = FALSE)
    rownames(aa_variation) <- NULL
  }

  return(aa_variation)
}

#' Convert HLA calls to variables
#'
#' \code{hlaToVariable} converts HLA calls data frame to additional variables.
#'
#' \code{dictionary} file should be a tsv format with header and two columns.
#' First column should hold allele numbers and second corresponding additional
#' variables. Optionally a data frame formatted in the same manner can be passed
#' instead.
#'
#' \code{dictionary} can be also used to access dictionaries shipped with the
#' package. They can be referred to by using one of the following strings:
#' \describe{
#'   \item{\code{"allele_HLA_Bw"}}{
#'     Translates HLA-B alleles together with A*23, A*24 and A*32 into Bw4 and
#'     Bw6 allele groups. In some cases HLA alleles containing Bw4 epitope, on
#'     nucleotide level actually carries a premature stop codon. Meaning that
#'     although on nucleotide level the allele would encode a Bw4 epitope it's
#'     not really there and it is assigned to Bw6 group. However in 4-digit
#'     resolution these alleles can not be distinguished from other Bw4 groups.
#'     Since alleles with premature stop codons are rare, Bw4 group is assigned.
#'   }
#'   \item{\code{"allele_HLA-B_only_Bw"}}{
#'     Translates HLA-B alleles (without A*23, A*24 and A*32) into Bw4 and Bw6
#'     allele groups.
#'   }
#'   \item{\code{"allele_HLA-C_C1-2"}}{
#'     Translates HLA-C alleles into C1 and C2 allele groups.
#'   }
#'   \item{\code{"allele_HLA_supertype"}}{
#'    Translates HLA-A and HLA-B alleles into supertypes, a classification that
#'    group HLA alleles based on peptide binding specificities.
#'   }
#'   \item{\code{"allele_HLA_Ggroup"}}{
#'     Translates HLA alleles into G groups, which defines amino acid identity
#'     only in the exons relevant for peptide binding. Note that alleles
#'     DRB1*01:01:01 and DRB1*01:16 match more than one G group, here this
#'     ambiguity was removed by deleting matching with DRB5*01:01:01G group.
#'   }
#' }
#'
#' \code{reduce} control if conversion should happen in a greedy way, such that
#' if some HLA number cannot be converted, it's resolution is reduced by 2 and
#' another attempt is taken. This process stops when alleles cannot be further
#' reduced or all have been successfully converted.
#'
#' @inheritParams checkHlaCallsFormat
#' @inheritParams convertAlleleToVariable
#' @param reduce Logical indicating if function should try to reduce allele
#'   resolution when no matching entry  in the dictionary is found. See details.
#' @param na.value Vector of length one speciyfing value for alleles with
#'   no matching entry in \code{dictionary}. Default is to use \code{0}.
#' @param nacols.rm Logical indicating if result columns that contain only
#'   \code{NA} should be removed.
#'
#' @return Data frame of HLA variables.
#'
#' @examples
#' hlaToVariable(MiDAS_tut_HLA, dictionary = "allele_HLA_supertype")
#'
#' @importFrom assertthat assert_that is.string see_if
#' @importFrom rlang warn
#' @export
hlaToVariable <- function(hla_calls,
                          dictionary,
                          reduce = TRUE,
                          na.value = 0,
                          nacols.rm = TRUE) {
  assert_that(
    checkHlaCallsFormat(hla_calls),
    isTRUEorFALSE(reduce),
    see_if(length(na.value) == 1, msg = "na.value length must equal 1."),
    isTRUEorFALSE(nacols.rm)
  )

  if (is.string(dictionary)) {
    lib <- listMiDASDictionaries()
    if (dictionary %in% lib) {
      if (dictionary %in% c("allele_HLA-B_Bw", "allele_HLA-Bw_only_B")) {
        warn("In ambiguous cases Bw4 will be assigned! See 'hlaToVariable' documentation for more details.")
      }
      dictionary <- system.file(
        "extdata",
        paste0("Match_", dictionary, ".txt"),
        package = "MiDAS"
      )
    }
  }

  variable <- as.data.frame(
    lapply(hla_calls[, -1], convertAlleleToVariable, dictionary = dictionary),
    stringsAsFactors = FALSE,
    row.names = NULL
  )

  if (reduce) {
    max_resolution <- getAlleleResolution(unlist(hla_calls[, -1]))
    max_resolution <- max(max_resolution, na.rm = TRUE)
    while (any(is.na(variable)) & max_resolution > 2) {
      max_resolution <- max_resolution - 2
      hla_calls <- reduceHlaCalls(hla_calls, resolution = max_resolution)
      variable[is.na(variable)] <- convertAlleleToVariable(
        allele = hla_calls[, -1][is.na(variable)],
        dictionary = dictionary
      )
    }
  }

  # add dictionary prefix to column names
  if (is.string(dictionary)) {
    dict_prefix <- gsub(".txt$", "", gsub("^.*_", "", dictionary))
  } else {
    dict_prefix <- colnames(dictionary)[2] # colnames are allele, name_of_variable
  }
  colnames(variable) <- paste0(dict_prefix, "_", colnames(variable))

  # get all na columns
  j <- vapply(variable, function(x) ! all(is.na(x)), logical(length = 1))

  if (nacols.rm) {
    variable <- variable[, j, drop = FALSE]
  }

  variable <- cbind(hla_calls[, 1, drop = FALSE], variable, stringsAsFactors = FALSE)
  colnames(variable) <- c("ID", colnames(variable[, -1]))

  if (ncol(variable) <= 1) {
    warn("HLA alleles could not be converted to any new variables.")
  }

  return(variable)
}

#' Reduce HLA calls resolution
#'
#' \code{reduceHlaCalls} reduces HLA calls data frame to specified resolution.
#'
#' Alleles with resolution greater than \code{resolution} or optional suffixes
#' are returned unchanged.
#'
#' @inheritParams checkHlaCallsFormat
#' @inheritParams reduceAlleleResolution
#'
#' @return HLA calls reduced to specified resolution.
#'
#' @examples
#' reduceHlaCalls(MiDAS_tut_HLA, resolution = 2)
#'
#' @export
reduceHlaCalls <- function(hla_calls, resolution = 4) {
  assert_that(checkHlaCallsFormat(hla_calls))
  hla_calls[, -1] <- as.data.frame(
    lapply(hla_calls[, -1], reduceAlleleResolution, resolution = resolution),
    stringsAsFactors = FALSE
  )

  return(hla_calls)
}

#' Transform HLA calls to counts table
#'
#' \code{hlaCallsToCounts} converts HLA calls data frame into a counts table.
#'
#' @inheritParams checkHlaCallsFormat
#' @param check_hla_format Logical indicating if \code{hla_calls} format should
#'   be checked. This is useful if one wants to use \code{hlaCallsToCounts} with
#'   input not adhering to HLA nomenclature standards. See examples.
#'
#' @return HLA allele counts data frame.
#'
#' @importFrom assertthat assert_that is.string
#' @importFrom qdapTools mtabulate
#'
hlaCallsToCounts <- function(hla_calls,
                             check_hla_format = TRUE) {
  assert_that(
    isTRUEorFALSE(check_hla_format),
    if (check_hla_format) {
      checkHlaCallsFormat(hla_calls)
    } else {
      TRUE
    }
  )

  hla_counts <- hla_calls[, -1, drop = FALSE]
  hla_counts <- mtabulate(as.data.frame(t(hla_counts)))
  rownames(hla_counts) <- NULL

  hla_counts <- cbind(ID = hla_calls[, 1, drop = FALSE],
                      hla_counts,
                      stringsAsFactors = FALSE
  )

  return(hla_counts)
}

#' Calculate HLA allele frequencies
#'
#' \code{getHlaFrequencies} calculates allele frequencies in HLA calls data
#' frame.
#'
#' Both gene copies are taken into consideration for frequencies calculation,
#' \code{frequency = n / (2 * j)} where \code{n} is the number of allele
#' occurrences and \code{j} is the number of samples in \code{hla_calls}.
#'
#' @inheritParams checkHlaCallsFormat
#' @inheritParams getExperimentFrequencies
#' @param compare Logical flag indicating if \code{hla_calls} frequencies
#'   should be compared to reference frequencies given in \code{ref}.
#' @param ref_pop Character vector giving names of reference populations in
#'   \code{ref} to compare with. Optionally vector can be named, then those
#'   names will be used as population names.
#' @param ref Data frame giving reference allele frequencies. See
#'   \code{\link{allele_frequencies}} for an example.
#'
#' @return Data frame containing HLA alleles and their corresponding frequencies.
#'
#' @examples
#' getHlaFrequencies(MiDAS_tut_HLA)
#'
#' @importFrom assertthat assert_that
#' @importFrom dplyr left_join
#' @importFrom formattable percent
#' @export
getHlaFrequencies <- function(hla_calls,
                              carrier_frequency = FALSE,
                              compare = FALSE,
                              ref_pop = c(
                                "USA NMDP African American pop 2",
                                "USA NMDP Chinese",
                                "USA NMDP European Caucasian",
                                "USA NMDP Hispanic South or Central American",
                                "USA NMDP Japanese",
                                "USA NMDP North American Amerindian",
                                "USA NMDP South Asian Indian"
                              ),
                              ref = allele_frequencies) {
  assert_that(
    checkHlaCallsFormat(hla_calls),
    isTRUEorFALSE(compare),
    is.data.frame(ref),
    colnamesMatches(ref, c("var", "population", "frequency")),
    is.character(ref_pop),
    characterMatches(ref_pop, unique(ref$population))
  )

  allele <- unlist(hla_calls[, -1])
  allele_counts <- table(allele, useNA = "no")
  allele_freq <- allele_counts / (2 * nrow(hla_calls))

  allele_freq <- data.frame(
    allele = names(allele_counts),
    Counts = as.vector(allele_counts),
    Freq = as.vector(allele_freq),
    stringsAsFactors = FALSE
  )

  if (compare) {
    ref <- getReferenceFrequencies(ref, ref_pop, carrier_frequency)
    allele_freq <- left_join(allele_freq, ref, by = c("allele" = "var"))
  }

  # format frequencies as percent
  allele_freq[, -c(1, 2)] <-
    rapply(
      object = allele_freq[, -c(1, 2), drop = FALSE],
      f = function(col) percent(col),
      how = "replace"
    )

  return(allele_freq)
}

#' Calculate KIR genes frequencies
#'
#' \code{getKIRFrequencies} calculates KIR genes frequencies in KIR calls data
#' frame.
#'
#' @inheritParams checkKirCallsFormat
#'
#' @return Data frame containing KIR genes and their corresponding frequencies.
#'
#' @examples
#' getKIRFrequencies(MiDAS_tut_KIR)
#'
#' @importFrom assertthat assert_that
#' @importFrom dplyr left_join
#' @importFrom formattable percent
#' @export
getKIRFrequencies <- function(kir_calls) {
  assert_that(
    checkKirCallsFormat(kir_calls)
  )

  kir_sums <- colSums(kir_calls[, -1, drop = FALSE], na.rm = TRUE)
  kir_freq <- kir_sums / nrow(kir_calls)

  kir_freq <- data.frame(
    gene = names(kir_sums),
    Counts = kir_sums,
    Freq = percent(kir_freq),
    stringsAsFactors = FALSE
  )

  return(kir_freq)
}

#' Transform amino acid variation data frame into counts table
#'
#' \code{aaVariationToCounts} convert amino acid variation data frame into
#' counts table.
#'
#' @param aa_variation Amino acid variation data frame as returned by
#'   \link{hlaToAAVariation}.
#'
#' @return Amino acid counts data frame.
#'
#' @importFrom assertthat assert_that is.string
#' @importFrom qdapTools mtabulate
#' @importFrom stats na.omit
#'
aaVariationToCounts <- function(aa_variation) {
  assert_that(
    is.data.frame(aa_variation),
    see_if(colnames(aa_variation)[1] == "ID",
           msg = "first column of aa_variation must be named ID"
    )
  )

  ids <- aa_variation[, 1]
  aa_counts <- aa_variation[, -1]
  aa_ids <- colnames(aa_variation[, -1])
  aa_ids <- gsub("_[12]_AA", "", aa_ids)
  aa_counts <- lapply(1:(ncol(aa_counts)),
                         function(i) {
                           x <- paste(aa_ids[i], aa_counts[, i], sep = "_")
                           x[is.na(aa_counts[, i])] <- NA
                           return(x)
                         }
  )
  ord <- na.omit(unique(unlist(aa_counts)))
  aa_counts <- do.call(rbind, aa_counts)
  aa_counts <- mtabulate(as.data.frame(aa_counts, stringsAsFactors = FALSE))
  rownames(aa_counts) <- NULL
  aa_counts <- aa_counts[, ord]
  aa_counts <- cbind(ID = aa_variation[, 1, drop = FALSE],
                     aa_counts,
                     stringsAsFactors = FALSE
  )

  return(aa_counts)
}

#' Calculate amino acid frequencies
#'
#' \code{getAAFrequencies} calculates amino acid frequencies in amino acid
#' data frame.
#'
#' Both gene copies are taken into consideration for frequencies calculation,
#' \code{frequency = n / (2 * j)} where \code{n} is the number of amino acid
#' occurrences and \code{j} is the number of samples in \code{aa_variation}.
#'
#' @inheritParams aaVariationToCounts
#'
#' @return Data frame containing amino acid positions and their corresponding
#'   frequencies.
#'
#' @examples
#' aa_variation <- hlaToAAVariation(MiDAS_tut_HLA)
#' getAAFrequencies(aa_variation)
#'
#' @importFrom assertthat assert_that
#' @export
getAAFrequencies <- function(aa_variation) {
  assert_that(
    is.data.frame(aa_variation),
    see_if(colnames(aa_variation)[1] == "ID",
           msg = "first column of aa_variation must be named ID"
    )
  )

  aa_pos <- aa_variation[, -1]
  aa_ids <- colnames(aa_variation[, -1])
  aa_ids <- gsub("_[12]_AA", "", aa_ids)
  aa_pos <- lapply(1:(ncol(aa_pos)),
                   function(i) {
                     paste(aa_ids[i], aa_pos[, i], sep = "_")
                   }
  )
  aa_pos <- unlist(aa_pos)

  aa_freq <- table(aa_pos, useNA = "no") / (2 * nrow(aa_variation))
  aa_freq <- as.data.frame(aa_freq, stringsAsFactors = FALSE)


  return(aa_freq)
}

#' Pretty format statistical analysis results helper
#'
#' \code{formatResults} format statistical analysis results table to html or
#' latex format.
#'
#' @param results Tibble as returned by \code{\link{runMiDAS}}.
#' @param filter_by Character vector specifying conditional expression used to
#'   filter \code{results}, this is equivalent to \code{...} argument passed to
#'   \code{\link[dplyr]{filter}}.
#' @param arrange_by Character vector specifying variable names to use for
#'   sorting. Equivalent to \code{...} argument passed to
#'   \code{\link[dplyr]{arrange}}.
#' @param select_cols Character vector specifying variable names that should be
#'   included in the output table. Can be also used to rename selected
#'   variables, see examples.
#' @param format String \code{"latex"} or \code{"html"}.
#' @param header String specifying header for result table. If \code{NULL}
#'   no header is added.
#'
#' @return Character vector of formatted table source code.
#'
#' @examples
#' \dontrun{
#' midas <- prepareMiDAS(hla_calls = MiDAS_tut_HLA,
#'                       colData = MiDAS_tut_pheno,
#'                       experiment = "hla_alleles")
#' object <- lm(disease ~ term, data = midas)
#' res <- runMiDAS(object, experiment = "hla_alleles")
#' formatResults(res,
#'               filter_by = c("p.value <= 0.05", "estimate > 0"),
#'               arrange_by = c("p.value * estimate"),
#'               select_cols = c("allele", "p-value" = "p.value"),
#'               format = "html",
#'               header = "HLA allelic associations")
#' }
#'
#' @importFrom assertthat assert_that
#' @importFrom dplyr arrange filter select
#' @importFrom knitr kable
#' @importFrom kableExtra add_header_above kable_styling scroll_box
#' @importFrom magrittr %>% %<>%
#' @importFrom stats setNames
#' @importFrom rlang parse_exprs .data
#'
formatResults <- function(results,
                          filter_by = "p.value <= 0.05",
                          arrange_by = "p.value",
                          select_cols = c("term", "estimate", "std.error", "p.value", "p.adjusted"),
                          format = c("html", "latex"),
                          header = NULL
                          ) {
  assert_that(
    is.character(filter_by),
    is.character(arrange_by),
    is.character(select_cols),
    is.string(format),
    stringMatches(format, choice = c("html", "latex")),
    isStringOrNULL(header)
  )

  filter_by <- parse_exprs(filter_by)
  arrange_by <- parse_exprs(arrange_by)

  results %<>%
    filter(!!! filter_by) %>%
    arrange(!!! arrange_by) %>%
    select(select_cols)

  if (format == "html" & isTRUE(getOption("knitr.in.progress"))) {
    results <-
      rapply(
        results,
        f = gsub,
        classes = "character",
        how = "replace",
        pattern = "(\\*)",
        replacement = "\\\\\\1"
      )
  }

  if (! (is.null(header) & format == "html")) {
    header <- setNames(ncol(results), header) # Still if format is 'latex' and header = NULL the result is not visualy appealing, and without it gives error. Issue created on github: https://github.com/haozhu233/kableExtra/issues/387
  }

  results %<>%
    kable(format = format, format.args = list(digits = 4, scientific = -5)) %>%
    add_header_above(header = header)

  if (format == "html") {
    results %<>%
      kable_styling(bootstrap_options = c("striped", "hover", "condensed")) %>%
      scroll_box(width = "100%", height = "200px")
  }

  return(results)
}

#' Create association analysis results table in HTML or LaTeX
#'
#' \code{kableResults} convert results table (\code{\link{runMiDAS}} output) to
#' HTML or LaTeX format.
#'
#' @inheritParams formatResults
#' @param colnames Character vector of form \code{c("new_name" = "old_name")},
#'   used to rename \code{results} colnames.
#' @param header String specifying results table header.
#' @param pvalue_cutoff Number specifying p-value cutoff for results to be
#'   included in output. If \code{NULL} no filtering is done.
#'
#' @return Association analysis results table in HTML or LaTeX.
#'
#' @examples
#' midas <- prepareMiDAS(hla_calls = MiDAS_tut_HLA,
#'                       colData = MiDAS_tut_pheno,
#'                       experiment = "hla_alleles")
#' object <- lm(disease ~ term, data = midas)
#' res <- runMiDAS(object, experiment = "hla_alleles")
#' kableResults(results = res,
#'              colnames = c("HLA allele" = "allele"))
#'
#' @importFrom assertthat assert_that is.number is.string see_if
#' @importFrom dplyr ends_with mutate_at vars
#' @importFrom magrittr %>% %<>%
#' @importFrom rlang has_name list2 parse_expr warn !! :=
#' @export
kableResults <- function(results,
                         colnames = NULL,
                         header = "MiDAS analysis results",
                         pvalue_cutoff = NULL,
                         format = getOption("knitr.table.format")) {
  assert_that(
    is.data.frame(results),
    isCharacterOrNULL(colnames),
    isNumberOrNULL(pvalue_cutoff),
    is.string(format),
    stringMatches(format, choice = c("html", "latex"))
  )
  if (! is.null(colnames)) {
    assert_that(
      characterMatches(colnames, choice = colnames(results))
    )
  }

  filter_by <- ifelse(
    test = is.null(pvalue_cutoff),
    yes = "p.value <= 1",
    no = sprintf("p.value < %f", pvalue_cutoff)
  )

  # create rename vector
  select_cols <- colnames(results)
  names(select_cols) <- select_cols
  i <- na.omit(match(x = select_cols, table = colnames))
  names(select_cols)[i] <- names(colnames)

  # replace .percent with %
  names(select_cols) <- gsub(".percent", " [%]", names(select_cols))

  results %<>%
    formatResults(
      filter_by = filter_by,
      arrange_by = "p.value",
      select_cols = select_cols,
      format = format,
      header = header
    )

  return(results)
}

#' Convert counts table to variables
#'
#' \code{countsToVariables} converts counts table to additional variables.
#'
#' \code{dictionary} file should be a tsv format with header and two columns.
#' First column should be named \code{"Name"} and hold variable name, second
#' should be named \code{"Expression"} and hold expression used to identify
#' variable (eg. \code{"KIR2DL3 & ! KIR2DL2"} will match all samples with
#' \code{KIR2DL3} and without \code{KIR2DL2}). Optionally a data frame formatted
#' in the same manner can be passed instead.
#'
#' Dictionaries shipped with the package:
#' \describe{
#'   \item{\code{kir_haplotypes}}{
#'     KIR genes to KIR haplotypes dictionary.
#'   }
#' }
#'
#' @inheritParams hlaToVariable
#' @param counts Data frame with counts, such as returned by
#'   \code{\link{hlaCallsToCounts}} function. First column should contain
#'   samples IDs, following columns should contain counts (natural numbers
#'   including zero).
#' @param dictionary Path to file containing variables dictionary or data
#'   frame. See details for further explanations.
#' @param na.value Vector of length one speciyfing value for variables with no
#'   matching entry in \code{dictionary}. Default is to use \code{0}.
#'
#' @return Data frame of indicators for new variables, with \code{1} and
#'   \code{0} signaling presence and  absence of a variable respectively.
#'
#' @examples
#' countsToVariables(MiDAS_tut_KIR, "kir_haplotypes")
#'
#' @importFrom assertthat assert_that is.string
#' @importFrom rlang parse_exprs
#' @export
countsToVariables <- function(counts,
                              dictionary,
                              na.value = NA,
                              nacols.rm = TRUE) {
  assert_that(
    checkColDataFormat(counts),
    see_if(length(na.value) == 1, msg = "na.value length must equal 1."),
    isTRUEorFALSE(nacols.rm)
  )

  if (is.string(dictionary)) {
    pattern <- paste0("counts_", dictionary)
    dict_path <- listMiDASDictionaries(pattern = pattern, file.names = TRUE)
    if(length(dict_path) == 0) {
      dict_path <- dictionary
    }
    assert_that(is.readable(dict_path))

    dictionary <- read.table(
      file = dict_path,
      header = TRUE,
      sep = "\t",
      quote = "",
      stringsAsFactors = FALSE
    )
  }

  assert_that(
    is.data.frame(dictionary),
    colnamesMatches(x = dictionary, cols = c("Name", "Expression"))
  )

  expressions <- dictionary[, "Expression", drop = TRUE]
  expressions <- parse_exprs(expressions)

  variables <- lapply(
    X = expressions,
    FUN = function(expr) {
      vars <- all.vars(expr)
      if (all(has_name(counts, vars))) {
        cl <- do.call(
          what = substitute,
          args = list(expr = expr, env = counts[, vars])
        )
        test <- eval(cl)
      } else {
        test <- rep(NA, nrow(counts))
      }

      test <- as.integer(test)
      return(test)
    }
  )

  res <- do.call(cbind, variables)
  colnames(res) <- dictionary[, "Name", drop = TRUE]

  # add ID column
  res <- cbind(counts[, 1, drop = FALSE], res)

  if (nacols.rm) {
    mask_na <- vapply(res, function(x) ! all(is.na(x)), logical(length = 1))
    res <- res[, mask_na, drop = FALSE]
  }

  return(res)
}

#' Get HLA - KIR interactions
#'
#' \code{getHlaKirInteractions} calculate presence-absence matrix of HLA - KIR
#' interactions.
#'
#' \code{hla_calls} are first reduced to all possible resolutions and converted
#' to additional variables, such as G groups, using dictionaries shipped with
#' the package.
#'
#' \code{interactions_dict} file should be a tsv format with header and two
#' columns. First column should be named \code{"Name"} and hold interactions
#' names, second should be named \code{"Expression"} and hold expression used to
#' identify interaction (eg. \code{"C2 & KIR2DL1"} will match all samples
#' with \code{C2} and \code{KIR2DL1}). The package is shipped with an interactions
#' file based on \href{https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6558367/}{Pende et al., 2019.}
#'
#' @inheritParams checkHlaCallsFormat
#' @inheritParams checkKirCallsFormat
#' @param interactions_dict Path to HLA - KIR interactions dictionary.
#'
#' @return Data frame with presence-absence indicators for HLA - KIR
#'   interactions.
#'
#' @importFrom assertthat assert_that is.string
#' @importFrom dplyr left_join
#' @importFrom magrittr %>%
#' @importFrom rlang warn
#' @importFrom stringi stri_detect_regex
#' @export
getHlaKirInteractions <- function(hla_calls,
                                  kir_calls,
                                  interactions_dict = system.file("extdata", "Match_counts_hla_kir_interactions.txt", package = "MiDAS")) {
  assert_that(
    checkHlaCallsFormat(hla_calls),
    checkKirCallsFormat(kir_calls),
    is.string(interactions_dict)
  )
  id_matches <- hla_calls[, 1, drop = TRUE] %in% kir_calls[, 1, drop = TRUE] %>%
    sum()
  assert_that(id_matches > 0,
              msg = "IDs in hla_calls doesn't match IDs in kir_calls"
  )
  if (nrow(hla_calls) != id_matches) {
    msg <- sprintf("%i IDs in hla_calls matched IDs in kir_calls", id_matches)
    warn(msg)
  }

  # transform hla_calls to all possible variables and resolutions
  ## in practice only subset of variables could be used but this should be more time proof
  midas_dicts <- listMiDASDictionaries(pattern = "allele_") %>%
    grep(pattern = "expression", value = TRUE, invert = TRUE)
  hla_variables <- Reduce(
    f = function(x, y) {
      left_join(x, hlaToVariable(hla_calls, dictionary = y), by = "ID")
    },
    x = midas_dicts,
    init = hla_calls
  )
  hla_max_resolution <- hla_calls[, -1] %>%
    unlist() %>%
    getAlleleResolution() %>%
    max(na.rm = TRUE)
  while (hla_max_resolution > 2) {
    hla_variables <- hla_calls %>%
      reduceHlaCalls(resolution = hla_max_resolution - 2) %>%
      left_join(x = hla_variables, by = "ID")
    hla_max_resolution <- hla_max_resolution - 2
  }

  hla_counts <- hlaCallsToCounts(hla_variables, check_hla_format = FALSE)
  hla_counts[, -1] <- ceiling(hla_counts[, -1] / 2) # reduce to presence / absence indicators

  counts <- left_join(hla_counts, kir_calls, by = "ID")
  interactions <- countsToVariables(counts, dictionary = interactions_dict)

  return(interactions)
}

#' Filter experiment by frequency
#'
#' Helper function for experiments filtering
#'
#' @inheritParams getExperimentFrequencies
#' @param lower_frequency_cutoff Positive number or \code{NULL}. Numbers greater
#'   than 1 are interpreted as number of feature occurrences, numbers between 0
#'   and 1 as fractions.
#' @param upper_frequency_cutoff Positive number or \code{NULL}. Numbers greater
#'   than 1 are interpreted as number of feature occurrences, numbers between 0
#'   and 1 as fractions.
#'
#' @return Filtered experiment matrix.
#'
#' @importFrom assertthat assert_that see_if is.number is.string
#' @importFrom magrittr %>%
#'
filterExperimentByFrequency <- function(experiment,
                                        carrier_frequency = FALSE,
                                        lower_frequency_cutoff = NULL,
                                        upper_frequency_cutoff = NULL) {
  inheritance_model_choice <- eval(formals()[["inheritance_model"]])
  assert_that(
    see_if(
     ! is.null(getExperimentPopulationMultiplicator(experiment)), # if pop_mul is not set frequency cannot be calculated
     msg = "Frequency filtration does not support provided experiment."
    ),
    isTRUEorFALSE(carrier_frequency),
    validateFrequencyCutoffs(lower_frequency_cutoff, upper_frequency_cutoff)
  )

  filtered_vars <- getExperimentFrequencies(
    experiment = experiment,
    carrier_frequency = carrier_frequency
  ) %>%
    getFrequencyMask(lower_frequency_cutoff = lower_frequency_cutoff,
                     upper_frequency_cutoff = upper_frequency_cutoff)
  mask <- rownames(experiment) %in% filtered_vars
  experiment <- experiment[mask, , drop = FALSE]

  return(experiment)
}

#' Calculate experiment's features frequencies
#'
#' \code{getExperimentFrequencies} calculate features frequencies.
#'
#' @param experiment Matrix or SummarizedExperiment object.
#' @param pop_mul Number by which number of samples should be multiplied to get
#'   the population size.
#' @param carrier_frequency Logical flag indicating if carrier frequency should
#'   be returned.
#' @param ref Wide format data frame with first column named "var" holding
#'   features matching \code{experiment} and specific populations frequencies in
#'   following columns. See \code{\link{getReferenceFrequencies}} for more
#'   details.
#'
#' @return Data frame containing variables and their corresponding frequencies.
#'
#' @importFrom assertthat assert_that is.string see_if
#' @importFrom formattable percent
#' @importFrom SummarizedExperiment assay
#'
getExperimentFrequencies <-
  function(experiment,
           pop_mul = NULL,
           carrier_frequency = FALSE,
           ref = NULL) {
    UseMethod("getExperimentFrequencies", experiment)
  }

#' @rdname getExperimentFrequencies
#' @method getExperimentFrequencies matrix
#'
getExperimentFrequencies.matrix <-
  function(experiment,
           pop_mul = NULL,
           carrier_frequency = FALSE,
           ref = NULL) {
    assert_that(
      isCountsOrZeros(experiment),
      is.number(pop_mul),
      isTRUEorFALSE(carrier_frequency)
    )
    if (! is.null(ref)) {
      assert_that(is.data.frame(ref))
    }

    if (carrier_frequency) {
      experiment <- applyInheritanceModel(experiment, "dominant")
      pop_mul <- 1 # carrier frequency does not take account of gene copies
    }

    counts_sums <- rowSums(experiment, na.rm = TRUE)
    allele_freq <- counts_sums / (pop_mul * ncol(experiment))

    counts_df <- data.frame(
      term = rownames(experiment),
      Counts = counts_sums,
      Freq = allele_freq,
      stringsAsFactors = FALSE
    )

    if (! is.null(ref)) {
      counts_df <-
        left_join(counts_df, ref, by = c("term" = "var"))
    }

    # format frequencies as percent
    counts_df[, -c(1, 2)] <-
      rapply(
        object = counts_df[, -c(1, 2), drop = FALSE],
        f = function(col) percent(col),
        how = "replace"
      )

    return(counts_df)
  }

#' @rdname getExperimentFrequencies
#' @method getExperimentFrequencies SummarizedExperiment
#'
getExperimentFrequencies.SummarizedExperiment <-
  function(experiment,
           pop_mul = NULL,
           carrier_frequency = FALSE,
           ref = NULL) {
    assert_that(
      isNumberOrNULL(pop_mul),
      isTRUEorFALSE(carrier_frequency)
    )
    if (! is.null(ref)) {
      assert_that(is.data.frame(ref))
    }

    counts <- assay(experiment)
    pop_mul <- getExperimentPopulationMultiplicator(experiment)
    getExperimentFrequencies(experiment = counts,
                             pop_mul = pop_mul,
                             carrier_frequency = carrier_frequency,
                             ref = ref
    )
  }

#' Apply inheritance model
#'
#' Helper function transforming experiment counts to selected
#' \code{inheritance_model}.
#'
#' Under \code{"dominant"} model homozygotes and heterozygotes are coded as
#' \code{1}. In \code{"recessive"} model homozygotes are coded as \code{1} and
#' other as \code{0}. In \code{"additive"} model homozygotes are coded as
#' \code{2} and heterozygotes as \code{1}.
#'
#' @param experiment Matrix or SummarizedExperiment object.
#' @param inheritance_model String specifying inheritance model to use.
#'  Available choices are \code{"dominant"}, \code{"recessive"},
#'  \code{"additive"}.
#'
#' @return \code{experiment} converted to specified inheritance model.
#'
applyInheritanceModel <-
  function(experiment,
           inheritance_model = c("dominant", "recessive", "additive")) {
    UseMethod("applyInheritanceModel", experiment)
  }

#' @rdname applyInheritanceModel
#' @method applyInheritanceModel matrix
#'
applyInheritanceModel.matrix <- function(experiment,
                                         inheritance_model =  c("dominant", "recessive", "additive")) {
  .classify <- function(x, val) {
    x <- x >= val
    mode(x) <- "integer"
    x
  }
  switch (
    inheritance_model,
    "additive" = experiment,
    "dominant" = .classify(experiment, 1), # ifelse(x >= 1, 1, 0)
    "recessive" = .classify(experiment, 2) # ifelse(x >= 2, 1, 0)
  )
}

#' @rdname applyInheritanceModel
#' @method applyInheritanceModel SummarizedExperiment
#'
applyInheritanceModel.SummarizedExperiment <- function(experiment,
                                                       inheritance_model =  c("dominant", "recessive", "additive")) {
  SummarizedExperiment::assay(experiment) <-
    applyInheritanceModel(SummarizedExperiment::assay(experiment), inheritance_model)

  return(experiment)
}

#' Helper function for filtering frequency data frame
#'
#' @inheritParams filterExperimentByFrequency
#' @param df Data frame as returned by \code{getExperimentFrequencies}.
#'
#' @return Character vector containing names of variables after filtration.
#'
#' @importFrom dplyr filter
#' @importFrom magrittr %>%
#'
getFrequencyMask <- function(df,
                             lower_frequency_cutoff = NULL,
                             upper_frequency_cutoff = NULL) {
  lower_frequency_cutoff <- ifelse(is.null(lower_frequency_cutoff), -Inf, lower_frequency_cutoff)
  upper_frequency_cutoff <- ifelse(is.null(upper_frequency_cutoff), Inf, upper_frequency_cutoff)
  freqs_are_float <- lower_frequency_cutoff <= 1 || upper_frequency_cutoff <= 1
  variables_freq <- df %>%
    filter(.data$Counts > lower_frequency_cutoff |
             freqs_are_float) %>%
    filter(.data$Freq > lower_frequency_cutoff |
             ! freqs_are_float) %>%
    filter(.data$Counts < upper_frequency_cutoff |
             freqs_are_float) %>%
    filter(.data$Freq < upper_frequency_cutoff |
             ! freqs_are_float)

   filtered_vars <- variables_freq$term

  return(filtered_vars)
}

#' Filter experiment by variable
#'
#' Helper function for experiments filtering
#'
#' @param experiment Matrix or SummarizedExperiment object.
#' @param variables Character vector specifying features to choose.
#'
#' @return Filtered \code{experiment} object.
#'
filterExperimentByVariables <-
  function(experiment, variables) {
    UseMethod("filterExperimentByVariables", experiment)
  }

#' @rdname filterExperimentByVariables
#' @method filterExperimentByVariables matrix
#'
filterExperimentByVariables.matrix <- function(experiment, variables) {
  return(experiment[variables, ])
}

#' @rdname filterExperimentByVariables
#' @method filterExperimentByVariables SummarizedExperiment
#'
filterExperimentByVariables.SummarizedExperiment <- function(experiment, variables) {
  experiment <- experiment[variables, ]
  S4Vectors::metadata(experiment)$omnibus_groups <-
    S4Vectors::metadata(experiment)$omnibus_groups[variables]

  return(experiment)
}

#' Get experiment's population multiplicator
#'
#' \code{getExperimentPopulationMultiplicator} extracts population multiplicator
#' from experiment's metadata.
#'
#' @param experiment Matrix or SummarizedExperiment object.
#'
#' @return Experiment's population multiplicator number.
#'
#' @importFrom S4Vectors metadata
#'
getExperimentPopulationMultiplicator <-
  function(experiment) {
    UseMethod("getExperimentPopulationMultiplicator", experiment)
  }

#' @rdname getExperimentPopulationMultiplicator
#' @method getExperimentPopulationMultiplicator matrix
#'
getExperimentPopulationMultiplicator.matrix <- function(experiment) return(NULL)

#' @rdname getExperimentPopulationMultiplicator
#' @method getExperimentPopulationMultiplicator SummarizedExperiment
#'
getExperimentPopulationMultiplicator.SummarizedExperiment <-
  function(experiment) {
      pop_mul <- metadata(experiment)[["pop_mul"]]
      return(pop_mul)
  }

#' Check if experiment is inheritance model applicable
#'
#' \code{isExperimentInheritanceModelApplicable} check experiment's metadata
#' for presence of \code{"inheritance_model_applicable"} flag, indicating if
#' inheritance model can be applied.
#'
#' @param experiment Matrix or SummarizedExperiment object.
#'
#' @return Logical flag.
#'
#' @importFrom S4Vectors metadata
#'
isExperimentInheritanceModelApplicable <-
  function(experiment) {
    UseMethod("isExperimentInheritanceModelApplicable", experiment)
  }

#' @rdname isExperimentInheritanceModelApplicable
#' @method isExperimentInheritanceModelApplicable matrix
#'
isExperimentInheritanceModelApplicable.matrix <- function(experiment) {
  return(FALSE)
}

#' @rdname isExperimentInheritanceModelApplicable
#' @method isExperimentInheritanceModelApplicable SummarizedExperiment
#'
isExperimentInheritanceModelApplicable.SummarizedExperiment <-
  function(experiment) {
    inheritance_model_applicable <-
      metadata(experiment)[["inheritance_model_applicable"]]
    return(inheritance_model_applicable)
  }
