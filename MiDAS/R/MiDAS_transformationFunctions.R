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
#'
#' @return Matrix containing variable amino acid positions.
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
                             alnpath = system.file("extdata", package = "MiDAS")){
  assert_that(
    checkHlaCallsFormat(hla_calls),
    is.flag(indels),
    is.flag(unkchar),
    is.dir(alnpath)
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

  return(aa_variation)
}

#' Converts hla calls data frame according to match table
#'
#' \code{hlaToVariable} converts hla calls data frame to additional variables
#' based on match table (dictionary).
#'
#' \code{dictionary} file should be a tsv format with header and two columns.
#' First column should hold allele numbers and second corresponding additional
#' variables. Optionally a data frame formatted in the same manner can be passed
#' instead.
#'
#' \code{dictionary} optional parameter that can be also used to access
#' matchings files shipped with the package. They can be referred to by using
#' one of the following strings:
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
#' \code{"4digit_allele_Ggroup"} HLA alleles can be re-coded in G groups,
#' which defines amino acid identity only in the exons relevant for peptide
#' binding.
#'
#' \code{"4digit_supertype"} A and B alleles can be assigned to so-called
#' supertypes.
#'
#' @inheritParams checkHlaCallsFormat
#' @inheritParams convertAlleleToVariable
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
                          nacols.rm = TRUE) {
  assert_that(
    is.data.frame(hla_calls),
    checkHlaCallsFormat(hla_calls),
    is.flag(nacols.rm)
  )
  if (is.string(dictionary)) {
    lib <- list.files(
      path = system.file("extdata", package = "MiDAS"),
      pattern = "^Match_.*txt$"
    )
    lib <- gsub("^Match_", "", gsub(".txt$", "", lib))
    if (dictionary %in% lib) {
      if (dictionary == "4digit_B-allele_Bw") {
        warn("In ambigious cases Bw4 will be assigned! See documentation for more details.")
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
    row.names = 1:nrow(hla_calls)
  )
  if (nacols.rm) {
    variable <- Filter(function(x) ! all(is.na(x)), variable)
  }
  variable <- cbind(hla_calls[, 1], variable)
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
#'
#' @return Data frame containing counts of HLA alleles counted according to
#'   specified model.
#'
#' @examples
#' file <- system.file("extdata", "HLAHD_output_example.txt", package = "MiDAS")
#' hla_calls <- readHlaCalls(file)
#' hlaCallsToCounts(hla_calls, inheritance_model = "additive")
#'
#' @importFrom assertthat assert_that is.string
#' @importFrom qdapTools mtabulate
#'
#' @export
hlaCallsToCounts <- function(hla_calls,
                             inheritance_model = c("dominant", "recessive", "additive")) {
  assert_that(
    checkHlaCallsFormat(hla_calls),
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

  hla_counts <- hla_calls[, -1, drop = FALSE]
  hla_counts <- mtabulate(as.data.frame(t(hla_counts)))
  rownames(hla_counts) <- NULL
  hla_counts <- hla_counts[, order(colnames(hla_counts))]

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

  hla_counts <- cbind(ID = hla_calls[, 1], hla_counts)

  return(hla_counts)
}
