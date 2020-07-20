#' Convert HLA allele numbers to amino acid variation matrix
#'
#' \code{hlaToAAVariation} convert HLA allele numbers data frame to a matrix
#' holding information on amino acid variation.
#'
#' Variable amino acid positions are found by comparing elements of the
#' alignment column wise. Some of the values in alignment can be treated
#' specially using \code{indels} and \code{unkchar} arguments. Function process
#' alignments for all HLA genes found in \code{hla_calls}.
#'
#' To infer variable amino acid position function uses protein alignment files
#' that are shipped with the package. Those files were downloaded from
#' \href{ftp://ftp.ebi.ac.uk/pub/databases/ipd/imgt/hla/alignments/}{EBI database}.
#'
#' @inheritParams checkHlaCallsFormat
#' @param indels Logical indicating whether indels should be considered when
#'   checking variability.
#' @param unkchar Logical indicating whether unknown characters in the alignment
#'   should be considered when checking variability.
#' @param as_df Logical indicating if data frame should be returned.
#'   Otherwise matrix is returned.
#'
#' @return Matrix or data frame containing variable amino acid positions. See
#'   \code{as_df} parameter.
#'
#'   Rownames corresponds to ID column of input data frame, and colnames to
#'   alignment positions for given genes. If no variation in amino acid
#'   alignments is found function return one column matrix filled with `NA`s.
#'
#' @examples
#' file <- system.file("extdata", "HLAHD_output_example.txt", package = "MiDAS")
#' hla_calls <- readHlaCalls(file)
#' hlaToAAVariation(hla_calls)
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

  # get names of genes and corresponding resolutions
  gene_names <- vapply(X = colnames(hla_calls),
                       FUN = function(x) stri_split_fixed(x, "_")[[1]][1],
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

  # hla_resolution <- vapply(X = gene_names_uniq,
  #                          FUN = function(x) {
  #                            x_numbers <- unlist(hla_calls[, gene_names == x])
  #                            x_res <- getAlleleResolution(na.omit(x_numbers))
  #                            return(min(x_res))
  #                          },
  #                          FUN.VALUE = numeric(length = 1),
  #                          USE.NAMES = TRUE
  # )

  # read alignment matrices and convert to desired resolution
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
      mask <- 1:nrow(hla_aln[[i]]) # This is tmp solution as NAs in character index gives oob error
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

#' Convert HLA calls data frame according to match table
#'
#' \code{hlaToVariable} convert HLA calls data frame to additional variables
#' based on match table (dictionary).
#'
#' \code{reduce} control if conversion should happen in a greedy way, such that
#' if some hla numbers cannot be converted, their resolution is reduced by 2 and
#' another attempt is taken. This iterative process stops when alleles cannot be
#' further reduced or all have been successfully converted.
#'
#' \code{dictionary} file should be a tsv format with header and two columns.
#' First column should hold allele numbers and second corresponding additional
#' variables. Optionally a data frame formatted in the same manner can be passed
#' instead.
#'
#' \code{dictionary} can be also used to access matching files shipped with the
#' package. They can be referred to by using one of the following strings (to
#' list available dictionaries use \code{\link{listMiDASDictionaries}}):
#' \describe{
#'   \item{\code{allele_HLA-A_expression}}{
#'     Reference data to impute expression levels for HLA-A alleles.
#'   }
#'   \item{\code{allele_HLA-B_Bw}}{
#'     B alleles can be grouped in allele groups Bw4 and Bw6. In some cases HLA
#'     alleles containing Bw4 epitope, on nucleotide level actually carries a
#'     premature stop codon. Meaning that although on nucleotide level the
#'     allele would encode a Bw4 epitope it's not really there and it is
#'     assigned to Bw6 group. However in 4-digit resolution these alleles can
#'     not be distinguished from other Bw4 groups. Since alleles with premature
#'     stop codons are rare in those ambiguous cases those are assigned to Bw4
#'     group.
#'   }
#'   \item{allele_HLA_Bw4+A23+A24+A32}{
#'     Extends \code{allele_HLA-B_Bw} dictionary by inclusion of A*23, A*24 and
#'     A*32 HLA alleles.
#'   }
#'   \item{\code{allele_HLA-C_C1-2}}{
#'     C alleles can be grouped in allele groups C1 and C2.
#'   }
#'   \item{\code{allele_HLA-C_expression}}{
#'     Reference data to impute expression levels for HLA-C alleles.
#'   }
#'   \item{\code{allele_HLA_supertype}}{
#'     A and B alleles can be assigned to so-called supertypes, a
#'     classification that group HLA alleles based on peptide binding
#'     specificities.
#'   }
#'   \item{\code{allele_HLA_Ggroup}}{
#'     HLA alleles can be re-coded in G groups, which defines amino acid
#'     identity only in the exons relevant for peptide binding. Note that
#'     alleles "DRB1*01:01:01" and "DRB1*01:16" were matched with more than one
#'     G group, this ambiguity was removed by deleting matching with
#'     "DRB5*01:01:01G" group. Moreover in the original match file there were
#'     alleles named "DPA*...", here they are renamed to "DPA1*..." to adhere
#'     with HLA nomenclature.
#'   }
#' }
#'
#' @inheritParams checkHlaCallsFormat
#' @inheritParams convertAlleleToVariable
#' @param reduce logical indicating if function should try to reduce alleles
#'   resolution when no matching is found. See details for more details.
#' @param na.value Vector of length one speciyfing value for alleles with
#'   no values in dictionary. Default behaviour is to mark such instances with
#'  \code{0}, however in some cases \code{NA} might be more appropriate.
#' @param nacols.rm logical indicating if result columns that contain only
#'   \code{NA} should be removed.
#'
#' @return Data frame of HLA numbers converted to additional variables according
#'   to match table.
#'
#' @examples
#' file <- system.file("extdata", "HLAHD_output_example.txt", package = "MiDAS")
#' hla_calls <- readHlaCalls(file)
#' hlaToVariable(hla_calls, dictionary = "allele_HLA_supertype")
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
    lib <- listMiDASDictionaries(pattern = "allele")
    if (dictionary %in% lib) {
      if (dictionary == "allele_HLA-B_Bw") {
        warn("In ambiguous cases Bw4 will be assigned! See documentation for more details.")
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
  dict_prefix <- gsub(".txt$", "", gsub("^.*_", "", dictionary))
  colnames(variable) <- paste0(dict_prefix, "_", colnames(variable))

  # # set non-original NAs to 0
  # i <- is.na(variable) & ! is.na(hla_calls[, -1, drop = FALSE])

  # get all na columns
  j <- vapply(variable, function(x) ! all(is.na(x)), logical(length = 1))

  # variable[i] <- na.value get all na columns
  if (nacols.rm) {
    variable <- variable[, j, drop = FALSE]
  }

  variable <- cbind(hla_calls[, 1, drop = FALSE], variable, stringsAsFactors = FALSE)
  colnames(variable) <- c("ID", colnames(variable[, -1]))

  if (ncol(variable) <= 1) {
    warn("No new variables colud be found.")
  }

  return(variable)
}

#' Reduce HLA calls data frame resolution
#'
#' \code{reduceHlaCalls} reduce HLA calls data frame to specified resolution.
#'
#' If \code{resolution} is greater than resolution of \code{hla_calls} elements,
#' those elements will be unchanged. Elements with optional suffixes are not
#' reduced.
#'
#' @inheritParams checkHlaCallsFormat
#' @inheritParams reduceAlleleResolution
#'
#' @return Data frame containing HLA allele calls reduced to required
#'   resolution.
#'
#' @examples
#' file <- system.file("extdata", "HLAHD_output_example.txt", package = "MiDAS")
#' hla_calls <- readHlaCalls(file)
#' reduceHlaCalls(hla_calls, resolution = 2)
#'
#' @export
reduceHlaCalls <- function(hla_calls,
                           resolution = 4) {
  assert_that(checkHlaCallsFormat(hla_calls))
  hla_calls[, -1] <- as.data.frame(
    lapply(hla_calls[, -1], reduceAlleleResolution, resolution = resolution),
    stringsAsFactors = FALSE
  )

  return(hla_calls)
}

#' Transform HLA calls to counts table
#'
#' \code{hlaCallsToCounts} convert HLA calls data frame into counts table.
#'
#' @inheritParams checkHlaCallsFormat
#' @param check_hla_format Logical indicating if \code{hla_calls} format should
#'   be checked. This is useful if one wants to use \code{hlaCallsToCounts} with
#'   input not adhering to HLA nomenclature standards. See examples.
#'
#' @return Data frame containing counts of HLA alleles according to specified
#'   inheritance model.
#'
#' @examples
#' file <- system.file("extdata", "HLAHD_output_example.txt", package = "MiDAS")
#' hla_calls <- readHlaCalls(file)
#' hlaCallsToCounts(hla_calls)
#'
#' # usage with non-HLA alleles numbers input
#' hla_vars <- hlaToVariable(hla_calls, dictionary = "allele_HLA_supertype")
#' hlaCallsToCounts(hla_calls, check_hla_format = FALSE)
#'
#' @importFrom assertthat assert_that is.string
#' @importFrom qdapTools mtabulate
#'
#' @export
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

#' Calculate alleles frequencies
#'
#' \code{getHlaFrequencies} calculates alleles frequencies in HLA calls data
#' frame.
#'
#' Allele frequencies are counted in reference to sample taking both gene copies
#' into consideration. `n / (2 * j)` where `n` is the number of allele
#' occurrences and `j` is the sample size.
#'
#' @inheritParams checkHlaCallsFormat
#'
#' @return Data frame containing alleles and thier corresponding frequencies.
#'
#' @examples
#' file <- system.file("extdata", "HLAHD_output_example.txt", package = "MiDAS")
#' hla_calls <- readHlaCalls(file)
#' getHlaFrequencies(hla_calls)
#'
#' @importFrom assertthat assert_that
#'
#' @export
getHlaFrequencies <- function(hla_calls) {
  assert_that(
    checkHlaCallsFormat(hla_calls)
  )

  allele <- unlist(hla_calls[, -1])
  allele_freq <- table(allele, useNA = "no") / (2 * nrow(hla_calls))
  allele_freq <- as.data.frame(allele_freq, stringsAsFactors = FALSE)

  return(allele_freq)
}

#' Transform amino acid variations data frame to counts table
#'
#' \code{aaVariationToCounts} converts amino acid variations data frame into
#' counts table.
#'
#' @inheritParams hlaCallsToCounts
#' @param aa_variation Data frame holding amino acid variation data as returned
#'   by \link{hlaToAAVariation}.
#'
#' @return Data frame containing counts of amino acid at specific positions
#'   according to inheritance specified model.
#'
#' @seealso \code{\link{hlaToAAVariation}}
#'
#' @examples
#' file <- system.file("extdata", "HLAHD_output_example.txt", package = "MiDAS")
#' hla_calls <- readHlaCalls(file)
#' aa_variation <- hlaToAAVariation(hla_calls)
#' aaVariationToCounts(aa_variation)
#'
#' @importFrom assertthat assert_that is.string
#' @importFrom qdapTools mtabulate
#' @importFrom stats na.omit
#'
#' @export
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

#' Calculate amino acid's frequencies
#'
#' \code{getAAFrequencies} calculates amino acid's frequencies in amino acid
#' variations data frame.
#'
#' Amino acid's frequencies are counted in reference to sample taking both gene
#' copies into consideration. `n / (2 * j)` where `n` is the number of amino
#' acid occurrences and `j` is the sample size.
#'
#' @inheritParams aaVariationToCounts
#'
#' @return Data frame containing the amino acid's positions and their
#'   corresponding frequencies.
#'
#' @seealso \code{\link{hlaToAAVariation}}
#'
#' @examples
#' file <- system.file("extdata", "HLAHD_output_example.txt", package = "MiDAS")
#' hla_calls <- readHlaCalls(file)
#' aa_variation <- hlaToAAVariation(hla_calls)
#' getAAFrequencies(aa_variation)
#'
#' @importFrom assertthat assert_that
#'
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

#' Convert HLA counts table to HLA calls
#'
#' \code{countsToHlaCalls} convert counts table to HLA calls data frame, this
#' is useful when working with data from UK Biobank.
#'
#' Note that proper HLA calls reconstruction from counts table is only possible
#' under additive inheritance model. This mode of operation is the only one
#' implemented so the function will always treat counts table as coming from
#' \code{hlaCallsToCounts(hla_calls, inheritance_model = 'additive')}.
#'
#' @param counts Data frame with HLA alleles counts, as returned by
#'   \code{\link{hlaCallsToCounts}} function. First column should contain
#'   samples IDs, following columns should be named with valid HLA alleles
#'   numbers.
#'
#' @return Data frame containing HLA allele calls.
#'
#' @examples
#' file <- system.file("extdata", "HLAHD_output_example.txt", package = "MiDAS")
#' hla_calls <- readHlaCalls(file)
#' hla_counts <- hlaCallsToCounts(hla_calls)
#' countsToHlaCalls(hla_counts)
#'
#' @importFrom assertthat assert_that see_if
#' @importFrom uniqtag cumcount make_unique_all
#'
#' @export
countsToHlaCalls <- function(counts) {
  assert_that(
    see_if(! is.null(colnames(counts)),
           msg = "count table has no column names"
    ),
    see_if(! any(is.na(colnames(counts))),
           msg = "column names contains NA values"
    ),
    see_if(all(checkAlleleFormat(colnames(counts)[-1])),
           msg = "counts table column names contains improperly formated HLA alleles numbers"
    ),
    see_if(
      all(counts[-1] == 0 | counts[-1] == 1 | counts[-1] == 2, na.rm = TRUE),
      msg = "counts can only take values 0, 1 or 2"
    )
  )

  ids <- counts[1]
  counts <- counts[-1]
  counts[is.na(counts)] <- 0

  alleles <- colnames(counts)
  genes <- gsub("\\*.*$", "", alleles)
  genes <- sort(unique(genes))
  genes <- make_unique_all(rep(genes, each = 2), sep = "_")

  haplotypes <- apply(counts, 1, function(row) {
    hap_ids <- row != 0
    hap <- alleles[hap_ids]
    hap <- rep(hap, times = row[hap_ids])
    row_genes <- gsub("\\*.*", "", hap)
    assert_that(
      see_if(! any(cumcount(row_genes) > 2),
             msg = "some samples have more than two alleles per gene"
      )
    )
    names(hap) <- make_unique_all(row_genes, sep = "_")
    hap[genes]
  })
  haplotypes <- t(haplotypes)
  colnames(haplotypes) <- genes

  new_df <- cbind(ids, haplotypes, stringsAsFactors = FALSE)

  return(new_df)
}

#' Helper function for pretty formating statistical analysis results
#'
#' \code{formatResults} format statistical analysis results table to html or
#' latex format.
#'
#' @param results Tibble as returned by \code{\link{analyzeAssociations}}.
#' @param filter_by Character vector specifying conditional expression used to
#'   filter \code{results}, this is equivalent to \code{...} argument passed to
#'   \code{\link[dplyr]{filter}} except it has to be a character vector.
#' @param arrange_by Character vector specifying variable names to use for
#'   sorting. Equivalent to \code{...} argument passed to
#'   \code{\link[dplyr]{arrange}}.
#' @param select_cols Character vector specifying variable names that should be
#'   included in the output table. Can be also used to rename selected
#'   variables, see examples.
#' @param format String with possible values \code{"latex"} and \code{"html"}.
#' @param header String specifying header for result table. If \code{NULL}
#'   no header is added.
#'
#' @return Character vector of formatted table source code.
#'
#' @seealso \code{\link{runMiDAS}}, \code{\link{analyzeAssociations}},
#'   \code{\link{analyzeConditionalAssociations}}.
#'
#' @examples
#' \dontrun{
#' hla_calls <- readHlaCalls(system.file("extdata", "HLAHD_output_example.txt", package = "MiDAS"))
#' pheno <- read.table(
#'   system.file("extdata", "pheno_example.txt", package = "MiDAS"),
#'   header = TRUE)
#' midas_data <- prepareMiDAS(hla_calls, pheno, analysis_type = "hla_alleles")
#' object <- lm(OS ~ 1 + term, data = midas_data)
#' res <- analyzeAssociations(object, variables = colnames(midas_data)[-1])
#' formatResults(res,
#'               filter_by = c("p.value <= 0.05", "estimate > 0"),
#'               arrange_by = c("p.value * estimate"),
#'               select_cols = c("allele" = "term", "p.value"),
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
#' @export
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
  select_cols_quo <- parse_exprs(select_cols)
  names(select_cols_quo) <- names(select_cols)

  results %<>%
    filter(!!! filter_by) %>%
    arrange(!!! arrange_by) %>%
    select(!!! select_cols_quo)

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
    kable(format = format, digits = 50) %>% # TODO empty results throw corrupted data.frame warning
    add_header_above(header = header)

  if (format == "html") {
    results %<>%
      kable_styling(bootstrap_options = c("striped", "hover", "condensed")) %>%
      scroll_box(width = "100%", height = "200px")
  }

  return(results)
}

#' Pretty format association analysis results
#'
#' \code{kableResults} formats results table to specified format. It uses
#' \code{\link{formatResults}} with pre specified arguments to return pretty
#' formatted table depending on the type of analysis and model type.
#'
#' @inheritParams formatResults
#' @param cols Character vector specifying columns to kable. Names can be used
#'   to rename columns.
#' @param header String specifying results table header.
#' @param pvalue_cutoff Number specifying p-value cutoff for results to be
#'   included in output. If \code{NULL} no filtering is done.
#'
#' @return A character vector with pretty formatted \code{results} table.
#'
#' @seealso \code{\link{formatResults}}, \code{\link{runMiDAS}}
#'
#' @importFrom assertthat assert_that is.number is.string see_if
#' @importFrom dplyr ends_with mutate_at vars
#' @importFrom magrittr %>% %<>%
#' @importFrom rlang has_name list2 parse_expr warn !! :=
#'
kableResults <- function(results,
                         cols = c("estimate", "std.error", "p.value", "p.adjusted", "Ntotal", "Ntotal (%)" = "Ntotal.frequency", Npositive = "N R=1", Npositive.frequency = "N R=1 (%%)", Nnegative = "N R=0", Nnegative.frequency = "N R=0 (%%)"),
                         header = "MiDAS analysis results",
                         pvalue_cutoff = NULL,
                         format = getOption("knitr.table.format")) {
  assert_that(
    is.data.frame(results),
    is.character(cols),
    isNumberOrNULL(pvalue_cutoff),
    is.string(format),
    stringMatches(format, choice = c("html", "latex"))
  )

  filter_by <- ifelse(
    test = is.null(pvalue_cutoff),
    yes = "p.adjusted <= 1",
    no = sprintf("p.value <= %f", pvalue_cutoff)
  )
  passed_filter <- eval(parse_expr(filter_by), envir = as.list(results))
  if (! any(passed_filter, na.rm = TRUE)) {
    warn(sprintf("None of the results meets filtering criteria: %s", filter_by))
  }

  term_name <- colnames(results)[1] # term name is always added..
  if (term_name %in% cols) {
    cols <- cols[cols != term_name]
  }
  select_cols <- c(
    unlist(list2(
      !! term_name := colnames(results)[1] # term name
    )),
    cols
  )

  present_cols <- has_name(results, select_cols)
  if (test_present_cols <- ! all(present_cols)) {
    warn(
      sprintf(
        "Columns %s could't be found in results, will be ommited.",
        ifelse(test_present_cols,
               paste0("\"",
                      paste(
                        select_cols[! present_cols],
                        collapse = "\", \""
                       ),
                      "\""
               ),
               ""
        )
      )
    )
  }
  assert_that(
    sum(present_cols) != 0,
    msg = sprintf("results does not contain any of the following columns: %s",
                  paste(select_cols, collapse = ", ")
    )
  )
  select_cols <- select_cols[present_cols]

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

#' Convert counts data frame according to match table
#'
#' \code{countsToVariables} convert counts data frame to variables based on
#' match table (dictionary).
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
#'   \item{\code{hla_kir_interactions}}{
#'     HLA - KIR interactions based on Pende et al., 2019.
#'   }
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
#' @param dictionary Path to the file containing variables matchings or data
#'   frame providing this information. See details for further explanations.
#' @param na.value Vector of length one speciyfing value for variables for which
#'   no matching is found in \code{counts}. Default behaviour is to mark such
#'   instances with \code{NA}.
#'
#' @return Data frame of indicators for new variables, with \code{1} signaling
#'   presence of variable and \code{0} absence.
#'
#' @examples
#' file <- system.file("extdata", "KPI_output_example.txt", package = "MiDAS")
#' kir_counts <- readKPICalls(file)
#' countsToVariables(kir_counts, "kir_haplotypes")
#'
#' @importFrom assertthat assert_that is.string
#' @importFrom rlang parse_exprs
#'
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
#' \code{getHlaKirInteractions} calculates binary presence-absence matrix of HLA
#' - KIR interactions.
#'
#' In order to be able to compare input data with \code{interactions_dict}
#' \code{hla_calls} are first converted to variables such as G groups, using
#' matching files shipped with the packages. Moreover \code{hla_calls} are also
#' reduced to all possible resolutions.
#'
#' \code{interactions_dict} file should be a tsv format with header and two
#' columns. First column should be named \code{"Name"} and hold interactions
#' names, second should be named \code{"Expression"} and hold expression used to
#' identify interaction (eg. \code{"C2 & KIR2DL1"} will match all samples
#' with \code{C2} and \code{KIR2DL1}). The package is shipped with interactions
#' file created based on Pende, et al. 2019.
#'
#' @inheritParams checkHlaCallsFormat
#' @param kir_counts Data frame containing KIR genes counts, as return by
#'   \code{\link{readKPICalls}}.
#' @param interactions_dict Path to the file containing HLA - KIR interactions
#'   matchings. See details for further details.
#'
#' @return Data frame with binary presence-absence indicators for HLA - KIR
#'   interactions.
#'
#' @examples
#' hla_file <- system.file("extdata", "HLAHD_output_example.txt", package = "MiDAS")
#' hla_calls <- readHlaCalls(hla_file)
#' kir_file <- system.file("extdata", "KPI_output_example.txt", package = "MiDAS")
#' kir_counts <- readKPICalls(kir_file)
#' getHlaKirInteractions(hla_calls, kir_counts)
#'
#' @importFrom assertthat assert_that is.string
#' @importFrom dplyr left_join
#' @importFrom magrittr %>%
#' @importFrom rlang warn
#' @importFrom stringi stri_detect_regex
#'
#' @export
getHlaKirInteractions <- function(hla_calls,
                                  kir_counts,
                                  interactions_dict = system.file("extdata", "Match_counts_hla_kir_interactions.txt", package = "MiDAS")) {
  assert_that(
    checkHlaCallsFormat(hla_calls),
    checkKirCallsFormat(kir_counts),
    is.string(interactions_dict)
  )
  id_matches <- hla_calls[, 1, drop = TRUE] %in% kir_counts[, 1, drop = TRUE] %>%
    sum()
  assert_that(id_matches > 0,
              msg = "IDs in hla_calls doesn't match IDs in kir_counts"
  )
  if (nrow(hla_calls) != id_matches) {
    msg <- sprintf("%i IDs in hla_calls matched IDs in kir_counts", id_matches)
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

  counts <- left_join(hla_counts, kir_counts, by = "ID")
  interactions <- countsToVariables(counts, dictionary = interactions_dict)

  return(interactions)
}

#' Filter experiment by frequency
#'
#' Helper function for experiments filtering
#'
#' @inheritParams getExperimentFrequencies
#' @param lower_frequency_cutoff Number of positive value or \code{NULL}
#' @param upper_frequency_cutoff Number of positive value or \code{NULL}
#'
#' @return Matrix
#'
#' @importFrom assertthat assert_that see_if is.number is.string
#' @importFrom dplyr filter
#' @importFrom magrittr %>%
#'
filterExperimentByFrequency <- function(experiment,
                                        carrier_frequency = FALSE,
                                        lower_frequency_cutoff = NULL,
                                        upper_frequency_cutoff = NULL) {
  inheritance_model_choice <- eval(formals()[["inheritance_model"]])
  assert_that(
    is.matrix(experiment),
    isCountsOrZeros(experiment),
    isTRUEorFALSE(carrier_frequency),
    validateFrequencyCutoffs(lower_frequency_cutoff, upper_frequency_cutoff)
  )

  lower_frequency_cutoff <- ifelse(is.null(lower_frequency_cutoff), -Inf, lower_frequency_cutoff)
  upper_frequency_cutoff <- ifelse(is.null(upper_frequency_cutoff), Inf, upper_frequency_cutoff)
  freqs_are_float <- lower_frequency_cutoff <= 1 || upper_frequency_cutoff <= 1
  variables_freq <- experiment %>%
    getExperimentFrequencies(carrier_frequency = carrier_frequency) %>%
    filter(.data$Counts > lower_frequency_cutoff |
             freqs_are_float) %>%
    filter(.data$Freq > lower_frequency_cutoff |
             ! freqs_are_float) %>%
    filter(.data$Counts < upper_frequency_cutoff |
             freqs_are_float) %>%
    filter(.data$Freq < upper_frequency_cutoff |
             ! freqs_are_float)
  variables <- variables_freq$term
  experiment <- experiment[variables, , drop = FALSE]

  return(experiment)
}

#' Calculate variables frequencies TODO
#'
#' \code{getExperimentFrequencies} calculate variables frequencies based on counts
#' table, such as produced by \code{\link{hlaCallsToCounts}}.
#'
#' Variables frequencies are counted in reference to sample size, depending on
#' the inheritance model under which the counts table has been generated one
#' might need to take under consideration both gene copies. Here sample size is
#' assumed to be depended on both gene copies for \code{"additive"} inheritance
#' model (`n / (2 * j)` where `n` is the number of term occurrences and `j`
#' is the sample size). For other models the sample size is taken as is
#' (`n / j`).
#'
#' @inheritParams hlaCallsToCounts
#' @param experiment Matrix
#' @param type String "allele_frequency" or "carrier_frequency"
#'
#' @return Data frame containing variables, its corresponding total counts
#'   and frequencies.
#'
#' @seealso \code{\link{hlaCallsToCounts}}
#'
#' @importFrom assertthat assert_that is.string
#' @importFrom formattable percent
#' @importFrom SummarizedExperiment assay
#'
#' @export
#'
getExperimentFrequencies <-
  function(experiment,
           carrier_frequency = FALSE,
           ref = NULL) {
    UseMethod("getExperimentFrequencies", experiment)
  }

#' @rdname getExperimentFrequencies
#' @method getExperimentFrequencies matrix
#'
getExperimentFrequencies.matrix <-
  function(experiment,
           carrier_frequency = FALSE,
           ref = NULL) {
    inheritance_model_choice <-  eval(formals()[["inheritance_model"]])
    assert_that(
      is.matrix(experiment),
      isCountsOrZeros(experiment),
      isTRUEorFALSE(carrier_frequency)
      # TODO ref
    )

    if (carrier_frequency) {
      experiment <- applyInheritanceModel(experiment, "dominant")
    }

    # Under additive inheritance model population size equals 2 * nrow(counts_table), in other cases it's 1 * nrow(counts_table)
    counts_sums <- rowSums(experiment, na.rm = TRUE)
    allele_freq <- counts_sums / (2 * ncol(experiment)) # the population size is 2x because genes comes in two copies

    counts_df <- data.frame(
      term = rownames(experiment),
      Counts = counts_sums,
      Freq = percent(allele_freq),
      stringsAsFactors = FALSE
    )

    if (! is.null(ref)) {
      if (carrier_frequency) {
         ref[, -1] <- lapply(ref[, -1, drop = FALSE], function (x) 2 * x * (1 - x) + x^2) # HWE 2qp + q^2
      }
      counts_df <-
        left_join(counts_df, ref, by = c("term" = "var"))
    }

    return(counts_df)
  }

#' @rdname getExperimentFrequencies
#' @method getExperimentFrequencies SummarizedExperiment
#'
getExperimentFrequencies.SummarizedExperiment <-
  function(experiment,
           carrier_frequency = TRUE) {
    counts <- assay(experiment)
    getExperimentFrequencies(experiment = counts,
                             carrier_frequency = carrier_frequency
    )
  }

#' Format experiment to inheritance model
#'
#' Helper function transforming experiments counts to selected
#' \code{inheritance_model}. Under \code{"additive"} inheritance model
#' experiment is returned unchanged. For \code{"dominant"} occurence of the
#' variable is reported as \code{1} and \code{0} otherwise. Under
#' \code{"recessive"} inheritence model occurence of the variable occurence of
#' two copies of the variable is reported as \code{1} and \code{0} otherwise.
#'
#' @param experiment Matrix
#' @param inheritance_model String specifying inheritance model to use.
#'   Available choices are \code{"dominant"}, \code{"recessive"},
#'   \code{"additive"}. In \code{"dominant"} model homozygotes and heterozygotes
#'   are coded as \code{1}. In \code{"recessive"} model homozygotes are coded as
#'   \code{1} and all other as \code{0}. In \code{"additive"} model homozygotes
#'   are coded as \code{2} and heterozygotes as \code{1}.
#'
#' @return Matrix
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
  switch (
    inheritance_model,
    "additive" = experiment,
    "dominant" = ceiling(experiment / 2),
    "recessive" = floor(experiment / 2)
  )
}

#' @rdname applyInheritanceModel
#' @method applyInheritanceModel SummarizedExperiment
#'
#' @importFrom SummarizedExperiment assay
#'
applyInheritanceModel.SummarizedExperiment <- function(experiment,
                                                       inheritance_model =  c("dominant", "recessive", "additive")) {
  assay(experiment) <- applyInheritanceModel(assay(experiment), inheritance_model)

  return(experiment)
}
