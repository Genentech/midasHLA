#' Converts HLA allele numbers to amino acid variation
#'
#' \code{hlaToAAVariation} converts HLA allele numbers data frame to a matrix
#' holding information on amino acid level variation.
#'
#' \code{alnpath} can be used to provide path to directory containing custom
#' alignment files. Each alignment file have to be named following EBI database
#' convention GENENAME_prot.txt. By default \code{alnpath} points to directory
#' containing alignments files shipped with the package.
#'
#' @inheritParams checkHlaCallsFormat
#' @param indels Logical indicating whether indels should be considered as
#'   variability.
#' @param unkchar Logical indicating whether unknown characters in the alignment
#'   should be treated as variability.
#' @param alnpath String providing optional path to directory containing HLA
#'   alignment files. See details for further explanations.
#' @param as_df Logical indicating if data frame should be returned.
#'   Otherwise function matrix is returned.
#'
#' @return Data frame or matrix containing variable amino acid positions. See
#'   \code{as_df} parameter.
#'
#'   Rownames corresponds to ID column of input data frame, and colnames to
#'   alignment positions for given genes. If no variation in amino acids
#'   alignments is found function return one column matrix filled with `NA`.
#'
#' @examples
#' file <- system.file("extdata", "HLAHD_output_example.txt", package = "MiDAS")
#' hla_calls <- readHlaCalls(file)
#' hlaToAAVariation(hla_calls)
#'
#' @importFrom assertthat assert_that see_if is.dir is.flag
#' @importFrom stringi stri_split_fixed
#' @export
hlaToAAVariation <- function(hla_calls,
                             indels = TRUE,
                             unkchar = FALSE,
                             alnpath = system.file("extdata", package = "MiDAS"),
                             as_df = TRUE){
  assert_that(
    checkHlaCallsFormat(hla_calls),
    is.flag(indels),
    is.flag(unkchar),
    is.dir(alnpath),
    is.flag(as_df)
  )
  ids <- hla_calls[, 1]
  hla_calls <- hla_calls[, -1]

  # get names of genes and corresponding resolutions
  gene_names <- vapply(X = colnames(hla_calls),
                       FUN = function(x) stri_split_fixed(x, "_")[[1]][1],
                       FUN.VALUE = character(length = 1)
  )
  gene_names_uniq <- unique(gene_names)
  # discard genes for which no alignment files are available
  available_genes <- list.files(
    path = alnpath,
    pattern = "_prot.txt$"
  )
  assert_that(
    length(available_genes) >= 1,
    msg = sprintf("no alignment files was found in path %s", alnpath)
  )
  alnfiles_readable <- vapply(
    X = file.path(alnpath, available_genes),
    FUN = is.readable,
    FUN.VALUE = logical(length = 1)
  )
  assert_that(
    all(alnfiles_readable),
    msg = sprintf("files: %s are not readable",
                  paste(available_genes[!alnfiles_readable], collapse = ", ")
    )
  )
  available_genes <- vapply(
    X = stri_split_fixed(available_genes, "_prot.txt"),
    `[[`, 1,
    FUN.VALUE = character(length = 1)
  )
  gene_names_uniq <- gene_names_uniq[gene_names_uniq %in% available_genes]
  hla_resolution <- vapply(X = gene_names_uniq,
                           FUN = function(x) {
                             x_numbers <- unlist(hla_calls[, gene_names == x])
                             x_res <- getAlleleResolution(na.omit(x_numbers))
                             return(min(x_res))
                           },
                           FUN.VALUE = numeric(length = 1),
                           USE.NAMES = TRUE
  )

  # read alignment matrices and convert to desired resolution
  hla_aln <- lapply(X = gene_names_uniq,
                    FUN = function(x) {
                      aln_file <- file.path(alnpath, paste0(x, "_prot.txt"))
                      aln <- readHlaAlignments(
                        file = aln_file,
                        resolution = hla_resolution[x],
                        unkchar = "*"
                      )
                      return(aln)
                    }
  )

  # get aa variations for each gene
  aa_variation <- list()
  for (i in 1:length(gene_names_uniq)) {
    x_calls <- hla_calls[, gene_names == gene_names_uniq[i]]

    # check if there is possibility for variability
    x_calls_uniq <- na.omit(unique(unlist(x_calls)))
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
      colnames(x) <- paste0(allele, "_", "AA_", var_pos)
      return(x)
    })
    var_aln <- do.call(cbind, var_aln)
    ord <- as.vector(vapply(1:length(var_pos),
                  function(j) {
                    c(j, j + length(var_pos))
                  },
                  FUN.VALUE = numeric(length = 2)
    ))
    var_aln <- var_aln[, ord]

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

#' Converts hla calls data frame according to match table
#'
#' \code{hlaToVariable} converts hla calls data frame to additional variables
#' based on match table (dictionary).
#'
#' \code{reduce} controls if conversion should happen in a greedy way, such that
#' if some hla numbers cannot be converted, their resolution is reduced by 2 and
#' another attempt is taken. This iterative process stops when alleles cannot be
#' further reduced or all have been successfully converted.
#'
#' \code{dictionary} file should be a tsv format with header and two columns.
#' First column should hold allele numbers and second corresponding additional
#' variables. Optionally a data frame formatted in the same manner can be passed
#' instead.
#'
#' \code{dictionary} can be also used to access matchings files shipped with the
#' package. They can be referred to by using one of the following strings (to
#' list available dictionaries use \link{listMiDASDictionaries}):
#'
#' \code{"2digit_A-allele_expression"} reference data to impute expression
#' levels for HLA-A alleles.
#'
#' \code{"4digit_B-allele_Bw"} B alleles can be grouped in allele groups Bw4 and
#' Bw6. In some cases HLA alleles containing Bw4 epitope, on nucleotide level
#' actually carries a premature stop codon. Meaning that although on nucleotide
#' level the allele would encode a Bw4 epitope it's not really there and it is
#' assigned to Bw6 group. However in 4-digit resolution these alleles can not be
#' distinguished from other Bw4 groups. Since alleles with premature stop codons
#' are rare in those ambiguous cases those are assigned to Bw4 group.
#'
#' \code{"4digit_C-allele_C1-2"} C alleles can be grouped in allele groups C1
#' and C2.
#'
#' \code{"2digit_C-allele_expression"} reference data to impute expression
#' levels for HLA-C alleles.
#'
#' \code{"4digit_supertype"} A and B alleles can be assigned to so-called
#' supertypes.
#'
#' \code{"4digit_allele_Ggroup"} HLA alleles can be re-coded in G groups,
#' which defines amino acid identity only in the exons relevant for peptide
#' binding. Note that alleles "DRB1*01:01:01" and "DRB1*01:16" were matched with
#' more than one G group, this ambiguity was removed by deleting matching with
#' "DRB5*01:01:01G" group. Moreover in the original match file there were alleles
#' named "DPA*...", here they are renamed to "DPA1*..." to adhere with HLA
#' nomenclature.
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
#' @return Data frame of hla numbers converted to additional variables according
#'   to match table.
#'
#' @examples
#' file <- system.file("extdata", "HLAHD_output_example.txt", package = "MiDAS")
#' hla_calls <- readHlaCalls(file)
#' hlaToVariable(hla_calls, dictionary = "4digit_supertype")
#'
#' @importFrom assertthat assert_that is.string is.flag see_if
#' @importFrom rlang warn
#' @export
hlaToVariable <- function(hla_calls,
                          dictionary,
                          reduce = TRUE,
                          na.value = 0,
                          nacols.rm = TRUE) {
  assert_that(
    checkHlaCallsFormat(hla_calls),
    is.flag(reduce),
    see_if(length(na.value) == 1, msg = "na.value length must equal 1."),
    is.flag(nacols.rm)
  )

  if (is.string(dictionary)) {
    lib <- listMiDASDictionaries()
    if (dictionary %in% lib) {
      if (dictionary == "4digit_B-allele_Bw") {
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

  # set non-original NAs to 0
  i <- is.na(variable) & ! is.na(hla_calls[, -1, drop = FALSE])

  # get all na columns
  j <- vapply(variable, function(x) ! all(is.na(x)), logical(length = 1))

  variable[i] <- na.value
  if (nacols.rm) {
    variable <- variable[, j]
  }

  variable <- cbind(hla_calls[, 1], variable, stringsAsFactors = FALSE)
  colnames(variable) <- c("ID", colnames(variable[, -1]))

  return(variable)
}

#' Reduce hla calls data frame resolution
#'
#' \code{reduceHlaCalls} reduces hla calls data frame to specified resolution.
#'
#' If \code{resolution} is greater than resolution of \code{hla_calls} elements,
#' those elements will be returned unchanged. Elements with optional suffixes
#' are not reduced.
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
#' \code{hlaCallsToCounts} converts HLA calls data frame into counts table.
#'
#' @inheritParams checkHlaCallsFormat
#' @param inheritance_model String specifying inheritance model to use.
#'   Available choices are \code{"dominant"}, \code{"recessive"},
#'   \code{"additive"}. In \code{"dominant"} model homozygotes and heterozygotes
#'   are coded as \code{1}. In \code{"recessive"} model homozygotes are coded as
#'   \code{1} and all other as \code{0}. In \code{"additive"} model homozygotes
#'   are coded as \code{2} and heterozygotes as \code{1}.
#' @param check_hla_format Logical indicating if \code{hla_calls} format should
#'   be checked. This is useful if one wants to use \code{hlaCallsToCounts} with
#'   input not adhering to HLA nomenclature standards.
#'
#' @return Data frame containing counts of HLA alleles counted according to
#'   specified model.
#'
#' @examples
#' file <- system.file("extdata", "HLAHD_output_example.txt", package = "MiDAS")
#' hla_calls <- readHlaCalls(file)
#' hlaCallsToCounts(hla_calls, inheritance_model = "additive")
#'
#' @importFrom assertthat assert_that is.flag is.string
#' @importFrom qdapTools mtabulate
#'
#' @export
hlaCallsToCounts <- function(hla_calls,
                             inheritance_model = c("dominant", "recessive", "additive"),
                             check_hla_format = TRUE) {
  assert_that(
    is.string(inheritance_model),
    see_if(
      pmatch(inheritance_model,
             table = c("dominant", "recessive", "additive"),
             nomatch = 0
      ) != 0,
      msg = "inheritance_model should be one of 'dominant', 'recessive', 'additive'"
    ),
    is.flag(check_hla_format)
  )

  if (check_hla_format) {
    assert_that(checkHlaCallsFormat(hla_calls))
  }

  inheritance_model <- match.arg(inheritance_model)

  hla_counts <- hla_calls[, -1, drop = FALSE]
  i <- do.call(
    what = "cbind",
    args = lapply(
      X = hla_counts,
      FUN = function(x) {
        vapply(
          X = x,
          FUN = function(x) {
            x <- suppressWarnings(as.numeric(x))
            test <- ! is.na(x)

            return(test)
          },
          FUN.VALUE = logical(length = 1)
        )
      }
    )
  )
  hla_counts[i] <- NA
  hla_counts <- mtabulate(as.data.frame(t(hla_counts)))
  rownames(hla_counts) <- NULL
  hla_counts <- hla_counts[, order(colnames(hla_counts)), drop = FALSE]

  hla_counts <- switch(inheritance_model,
                       "dominant" = as.data.frame(
                         lapply(hla_counts,
                                function(x) ifelse(x == 2, 1, x)
                         ),
                         stringsAsFactors = FALSE,
                         optional = TRUE
                       ),
                       "recessive" = as.data.frame(
                         lapply(hla_counts,
                                function(x) ifelse(x == 2, 1, 0)
                         ),
                         stringsAsFactors = FALSE,
                         optional = TRUE
                       ),
                       "additive" = hla_counts # Do nothing this is default res
  )

  # set 0 to NAs where appropiate
  genes <- colnames(hla_calls[, -1, drop = FALSE])
  origin_dict <- data.frame(
    allele = unlist(hla_calls[, -1, drop = FALSE]),
    gene = rep(genes, each = nrow(hla_calls)),
    stringsAsFactors = FALSE
  )
  origin_dict <- origin_dict[! is.na(origin_dict$allele), ]
  for (col in colnames(hla_counts)) {
    origin <- origin_dict$gene[origin_dict$allele == col]
    na_i <- hla_calls[, origin, drop = FALSE]
    na_i <- is.na(na_i)
    na_i <- rowSums(na_i) == ncol(na_i)
    hla_counts[na_i, col] <- NA
  }

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
#' @return Data frame containing the allele and its corresponding frequencies.
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

#' Transform amino acids variation data frame to counts table
#'
#' \code{aaVariationToCounts} converts variation data frame data frame into
#' counts table.
#'
#' @inheritParams hlaCallsToCounts
#' @param aa_variation Data frame holding amino acid variation data as returned
#'   by \link{hlaToAAVariation}.
#'
#' @return Data frame containing counts of amino acids at specific positions
#'   counted according to specified model.
#'
#' @examples
#' file <- system.file("extdata", "HLAHD_output_example.txt", package = "MiDAS")
#' hla_calls <- readHlaCalls(file)
#' aa_variation <- hlaToAAVariation(hla_calls)
#' aaVariationToCounts(aa_variation, inheritance_model = "additive")
#'
#' @importFrom assertthat assert_that is.string
#' @importFrom qdapTools mtabulate
#'
#' @export
aaVariationToCounts <- function(aa_variation,
                                inheritance_model = c("dominant", "recessive", "additive")) {
  assert_that(
    is.data.frame(aa_variation),
    see_if(colnames(aa_variation)[1] == "ID",
           msg = "first column of aa_variation must be named ID"
    ),
    is.string(inheritance_model),
    see_if(
      pmatch(inheritance_model,
             table = c("dominant", "recessive", "additive"),
             nomatch = 0
      ) != 0,
      msg = "inheritance_model should be one of 'dominant', 'recessive', 'additive'"
    )
  )

  inheritance_model <- match.arg(inheritance_model)

  ids <- aa_variation[, 1]
  aa_counts <- aa_variation[, -1]
  aa_ids <- colnames(aa_variation[, -1])
  aa_ids <- gsub("_[12]_AA", "", aa_ids)
  aa_counts <- lapply(1:(ncol(aa_counts)),
                         function(i) {
                           paste(aa_ids[i], aa_counts[, i], sep = "_")
                         }
  )
  ord <- unique(unlist(aa_counts))
  aa_counts <- do.call(rbind, aa_counts)
  aa_counts <- mtabulate(as.data.frame(aa_counts, stringsAsFactors = FALSE))
  rownames(aa_counts) <- NULL
  aa_counts <- aa_counts[, ord]

  aa_counts <- switch(
    inheritance_model,
    "dominant" = as.data.frame(
      lapply(aa_counts,
             function(x)
               ifelse(x == 2, 1, x)),
      stringsAsFactors = FALSE,
      optional = TRUE
    ),
    "recessive" = as.data.frame(
      lapply(aa_counts,
             function(x)
               ifelse(x == 2, 1, 0)),
      stringsAsFactors = FALSE,
      optional = TRUE
    ),
    "additive" = aa_counts # Do nothing this is default res
  )

  aa_counts <- cbind(ID = aa_variation[, 1, drop = FALSE],
                     aa_counts,
                     stringsAsFactors = FALSE
  )

  return(aa_counts)
}

#' Calculate amino acids frequencies
#'
#' \code{getAAFrequencies} calculates amino acids frequencies in amino acids
#' data frame.
#'
#' Amino acids frequencies are counted in reference to sample taking both gene
#' copies into consideration. `n / (2 * j)` where `n` is the number of amino
#' acid occurrences and `j` is the sample size.
#'
#' @inheritParams aaVariationToCounts
#'
#' @return Data frame containing the amino acids with positions and its
#'   corresponding frequencies.
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
#' \code{countsToHlaCalls} converts counts table to hla calls data frame, this
#' is useful when working with data from UK Biobank.
#'
#' Note that proper HLA calls reconstruction from counts table is only possible
#' under additive inheritance model. This mode of operation is the only one
#' implemented so the function will always treat counts table as coming from
#' \code{hlaCallsToCounts(hla_calls, inheritance_model = 'additive')}.
#'
#' @param counts Data frame with HLA alleles counts, as returned by
#'   \link{hlaCallsToCounts} function. First column should contain samples IDs,
#'   following columns should be named with valid HLA alleles numbers.
#'
#' @return Data frame containing HLA allele calls.
#'
#' @examples
#' file <- system.file("extdata", "HLAHD_output_example.txt", package = "MiDAS")
#' hla_calls <- readHlaCalls(file)
#' hla_counts <- hlaCallsToCounts(hla_calls, inheritance_model = "additive")
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
#' \code{formatResults} formats statistical analysis results table to html or
#' latex format.
#'
#' @param results Tibble as returned by \link{analyzeAssociations}.
#' @param filter_by Character specifying conditional expression used to filter
#'   \code{results}, this is equivalent to \code{...} argument passed to
#'   \link[dplyr]{filter} except it has to be a character vector.
#' @param arrange_by Character specifying variable names to use for sorting.
#'   Equivalent to \code{...} argument passed to \link[dplyr]{arrange}.
#' @param select_cols Character specifying variable names that should be
#'   included in the output table. Can be also used to rename selected
#'   variables, see examples.
#' @param format A character string. Possible values are \code{"latex"} and
#'   \code{"html"}.
#' @param header String specifying header for result table. If \code{NULL}
#'   header is omitted.
#'
#' @return A character vector of the table source code.
#'
#' @examples
#' hla_calls <- readHlaCalls(system.file("extdata", "HLAHD_output_example.txt", package = "MiDAS"))
#' hla_counts <- hlaCallsToCounts(hla_calls, inheritance_model = "additive")
#' midas_data <- read.table(
#'   system.file("extdata", "pheno_example.txt", package = "MiDAS"),
#'   header = TRUE)
#' midas_data <- dplyr::left_join(x = midas_data, y = hla_counts, by = "ID")
#' object <- lm(OS ~ 1, data = midas_data)
#' res <- analyzeAssociations(object, variables = colnames(midas_data)[-1])
#' formatResults(res,
#'               filter_by = c("p.value <= 0.05", "estimate > 0"),
#'               arrange_by = c("p.value * estimate"),
#'               select_cols = c("allele" = "term", "p.value"),
#'               format = "html",
#'               header = "HLA allelic associations")
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

#' Calculate variables frequencies
#'
#' \code{getCountsFrequencies} calculate variables frequencies based on counts
#' table, such as produced by \link{hlaCallsToCounts}.
#'
#' Variables frequencies are counted in reference to sample size, depending on
#' the inheritance model under which the counts table has been generated one
#' might need to take under consideration both gene copies. Here sample size is
#' assumed to be depended on both gene copies if any count is greater than
#' \code{1} (`n / (2 * j)` where `n` is the number of term occurrences and `j`
#' is the sample size). If this is not the case the sample size is taken as is
#' (`n / j`).
#'
#' @param counts_table Data frame containing variables counts, such as produced
#'   by \link{hlaCallsToCounts}.
#'
#' @return Data frame containing variables, its corresponding total counts
#'   and frequencies.
#'
#' @examples
#' file <- system.file("extdata", "HLAHD_output_example.txt", package = "MiDAS")
#' hla_calls <- readHlaCalls(file)
#' hla_counts <- hlaCallsToCounts(hla_calls, inheritance_model = "additive")
#' getCountsFrequencies(hla_counts)
#'
#' @importFrom assertthat assert_that
#' @importFrom formattable percent
#'
#' @export
getCountsFrequencies <- function(counts_table) {
  assert_that(
    is.data.frame(counts_table),
    see_if(colnames(counts_table)[1] == "ID",
           msg = "first column of counts_table must be named ID"
    )
  )

  counts_table <- counts_table[-1]
  assert_that(isCountsOrZeros(counts_table))
  counts_sums <- colSums(counts_table, na.rm = TRUE)

  # Under additive inheritance model population size equals 2 * nrow(counts_table), in other cases it's 1 * nrow(counts_table)
  pop_mul <- ifelse(max(counts_table, na.rm = TRUE) > 1, 2, 1)
  counts_freq <- counts_sums / (pop_mul * nrow(counts_table))

  counts_df <- data.frame(
    term = colnames(counts_table),
    Counts = counts_sums,
    Freq = percent(counts_freq),
    stringsAsFactors = FALSE
  )

  return(counts_df)
}

#' Pretty format association analysis results
#'
#' \link{formatAssociationsResults} formats results table to specified
#' format. It uses \link{formatResults} with prespecifed arguments to return
#' nice table depending on the type of analysis and model type. This function is
#' intended only to be used internally by \link{analyzeMiDASData}.
#'
#' @inheritParams formatResults
#' @param type String specifying type of analysis from which \code{results} were
#'   produced. Possible values includes \code{'hla_allele'}, \code{'aa_level'},
#'   \code{'expression_level'}, \code{'allele_g_group'},
#'   \code{'allele_supertype'}, \code{'allele_group'}, \code{'kir_genes'},
#'   \code{'hla_kir_interactions'}, \code{'custom'}.
#' @param response_variable String giving the name of response variable, it is
#'   used to produce binary phenotype column names.
#' @param logistic Logical indicating if statistical model is logistic. If set
#'  to \code{TRUE}, estimate will be renamed to odds ratio.
#' @param pvalue_cutoff Number specifying p-value cutoff for results to be
#'   included in output. If \code{NULL} cutoff of \code{0.05} on
#'   \code{p.adjusted} value is used instead.
#'
#' @return A character vector with pretty formatted \code{results} table.
#'
#' @importFrom assertthat assert_that is.flag is.number is.string see_if
#' @importFrom dplyr ends_with mutate_at vars
#' @importFrom magrittr %>% %<>%
#' @importFrom rlang has_name list2 parse_expr warn !! :=
#'
formatAssociationsResults <- function(results,
                                      type = "hla_allele",
                                      response_variable = "R",
                                      logistic = FALSE,
                                      pvalue_cutoff = NULL,
                                      format = getOption("knitr.table.format")) {
  assert_that(
    is.string(type),
    stringMatches(type,
                  choice = c("hla_allele", "aa_level", "expression_level", "allele_g_group", "allele_supertype", "allele_group", "kir_genes", "hla_kir_interactions", "custom")
    ),
    is.string(response_variable),
    is.flag(logistic),
    isNumberOrNULL(pvalue_cutoff),
    is.string(format),
    stringMatches(format, choice = c("html", "latex"))
  )

  filter_by <- ifelse(
    test = is.null(pvalue_cutoff),
    yes = "p.adjusted <= 0.05",
    no = sprintf("p.value <= %f", pvalue_cutoff)
  )
  passed_filter <- eval(parse_expr(filter_by), envir = as.list(results))
  if (! any(passed_filter, na.rm = TRUE)) {
    warn(sprintf("None of the results meets filtering criteria: %s", filter_by))
  }

  estimate_name <- ifelse(logistic, "odds ratio", "estimate")
  term_name <- switch (type,
                       "hla_allele" = "allele",
                       "aa_level" = "aa",
                       "expression_level" = "allele",
                       "allele_g_group" = "g group",
                       "allele_supertype" = "supertype",
                       "allele_group" = "allele group",
                       "kir_genes" = "kir gene",
                       "hla_kir_interactions" = "hla kir interaction",
                       "term"
  )
  select_cols <- unlist(list2(
    !! term_name := "term",
    !! estimate_name := "estimate",
    "std.error",
    "p.value",
    "p.adjusted",
    "Ntotal",
    "Ntotal (%)" = "Ntotal.frequency",
    !! sprintf("N %s=1", response_variable) := "Npositive",
    !! sprintf("N %s=1 (%%)", response_variable) := "Npositive.frequency",
    !! sprintf("N %s=0", response_variable) := "Nnegative",
    !! sprintf("N %s=0 (%%)", response_variable) := "Nnegative.frequency"
  ))
  present_cols <- has_name(results, select_cols)
  assert_that(
    sum(present_cols) != 0,
    msg = sprintf("results does not contain any of the following columns: %s",
                  paste(select_cols, collapse = ", ")
    )
  )
  select_cols <- select_cols[present_cols]


  header <- switch (type,
                    "hla_allele" = "HLA allelic associations",
                    "aa_level" = "HLA AA associations",
                    "expression_level" = "HLA expression level associations",
                    "allele_g_group" = "HLA alleles G groups associations",
                    "allele_supertype" = "HLA alleles supertypes associations",
                    "allele_group" = "HLA alleles groups associations",
                    "kir_genes" = "KIR genes associations",
                    "hla_kir_interactions" = "HLA - KIR interaction associations",
                    "Associations results"
  )

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

#' Get HLA - KIR interactions
#'
#' \code{getHlaKirInteractions} calculates binary presence-absence matrix of HLA
#' - KIR interactions.
#'
#' @inheritParams checkHlaCallsFormat
#' @param kir_counts Data frame containing KIR genes counts, as return by
#'   \link{readKirCalls}.
#' @param interactions_dict Path to the file containing HLA - KIR interactions
#'   matchings. See details for further details.
#'
#' @return Data frame with binary presence-absence indicators for HLA - KIR
#'   interactions.
#'
#' In order to be able to compare input data with \code{interactions_dict}
#' \code{hla_calls} are first converted to variables such as G groups, using
#' matching files shipped with the packages. Moreover \code{hla_calls} are also
#' reduced to all possible resolutions.
#'
#' \code{interactions_dict} should be a tsv formatted file with following
#' columns: HLA, KIR, Affinity, Type. The package is shipped with interactions
#' file created based on Pende et al., 2019.
#'
#' @examples
#' hla_file <- system.file("extdata", "HLAHD_output_example.txt", package = "MiDAS")
#' hla_calls <- readHlaCalls(hla_file)
#' kir_file <- system.file("extdata", "KIP_output_example.txt", package = "MiDAS")
#' kir_counts <- readKirCalls(kir_file, counts = TRUE)
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
                                  interactions_dict = system.file("extdata", "Interactions_hla_kir.tsv", package = "MiDAS")) {
  assert_that(
    checkHlaCallsFormat(hla_calls),
    checkKirCountsFormat(kir_counts),
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
  midas_dicts <- listMiDASDictionaries() %>%
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

  hla_counts <- hlaCallsToCounts(hla_variables,
                                 inheritance_model = "dominant",
                                 check_hla_format = FALSE
  )

  # find hla - kir interactions
  interaction_vars <- left_join(hla_counts, kir_counts, by = "ID")
  vars <- colnames(interaction_vars[, -1, drop = FALSE])
  interactions_dict <- read.table(
    file = interactions_dict,
    sep = "\t",
    header = TRUE,
    stringsAsFactors = FALSE
  )

  posible_interactions <- vapply(
    X = interactions_dict$HLA,
    FUN = function(pattern) {
      any(
        stri_detect_regex(str = vars, pattern = pattern, max_count = 1),
        na.rm = TRUE
      )
    },
    FUN.VALUE = logical(1)
  ) &
    vapply(
      X = interactions_dict$KIR,
      FUN = function(pattern) {
        any(
          stri_detect_regex(str = vars, pattern = pattern, max_count = 1),
          na.rm = TRUE
        )
      },
      FUN.VALUE = logical(1)
    )

  interactions_dict <- interactions_dict[posible_interactions, ]

  # assert that interactions_dict is not empty

  interactions <- apply(
    X = interactions_dict[, c("HLA", "KIR")],
    MARGIN = 1,
    FUN = function(row) {
      hla_col <- grep(pattern = paste0("^", row[[1]], "$"), x = colnames(interaction_vars))
      hla_counts <- interaction_vars[, hla_col, drop = FALSE]
      hla_na_i <- apply(hla_counts, 1, function(x) all(is.na(x)))
      kir_col <- grep(pattern = paste0("^", row[[2]], "$"), x = colnames(interaction_vars))
      kir_counts <- interaction_vars[, kir_col, drop = FALSE]
      kir_na_i <- apply(kir_counts, 1, function(x) all(is.na(x)))
      x <- rowSums(hla_counts, na.rm = TRUE) + rowSums(kir_counts, na.rm = TRUE)
      x[hla_na_i | kir_na_i] <- NA
      x <- ifelse(x == 1, 0, x)
      x <- ifelse(x >= 2, 1, x) # it is assumed that kir_calls are in dominant mode..., >= is used here to handle more complex interactions like KIR3DL1_Bw4 + KIR3DL1_A*23 + KIR3DL1_A*24 + KIR3DL1_A*32
    }
  )
  colnames(interactions) <- interactions_dict$Name

  interactions <- cbind(hla_calls[, 1, drop = FALSE], interactions)

  return(interactions)
}

#' Converts counts data frame according to match table
#'
#' \code{countsToVariables} converts counts data frame to variables based on
#' match table (dictionary).
#'
#' @inheritParams hlaToVariable
#' @param counts Data frame with counts, such as returned by
#'   \link{hlaCallsToCounts} function. First column should contain samples IDs,
#'   following columns should contain counts (natural numbers including zero).
#' @param dictionary Path to the file containing variable matchings or data frame providing
#'   this information. See details for further explanations.
#' @param na.value Vector of length one speciyfing value for variables for which no
#'   matching is found in \code{counts}. Default behaviour is to mark such
#'   instances with \code{NA}.
#'
#' @return Data frame of indicators for new variables, with \code{1} signaling
#'   presence of variable and \code{0} absence.
#'
#' \code{dictionary} file should be a tsv format with header and two columns.
#' First column should be named \code{"Name"} and hold variable name, second
#' should be named \code{"Expression"} and hold expression used to identify
#' variable (eg. ...). Optionally a data frame formatted in the same manner can
#' be passed instead.
#'
#' Dictionaries shipped with package:
#'
#' \code{hla_kir_interactions} interactions based on Pende et al., 2019.
#'
#' \code{kir_haplotypes} kir haplotypes.
#'
#' @examples
#' file <- system.file("extdata", "KIP_output_example.txt", package = "MiDAS")
#' kir_counts <- readKirCalls(file)
#' countsToVariables(kir_counts, "kir_haplotypes")
#'
#' @importFrom assertthat assert_that is.flag is.string
#' @importFrom rlang parse_exprs
countsToVariables <- function(counts,
                              dictionary,
                              na.value = NA,
                              nacols.rm = TRUE) {
  assert_that(
    checkKirCountsFormat(counts),
    see_if(length(na.value) == 1, msg = "na.value length must equal 1."),
    is.flag(nacols.rm)
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
