#' Reads data table with HLA allele calls
#'
#' \code{readHlaCalls} reads table with HLA allele calls from file.
#'
#' Input file have to be a tsv formatted table with header. First column should
#' contain samples IDs, further columns should hold corresponding HLA allele
#' numbers.
#'
#' \code{resolution} parameter can be used to reduce HLA allele numbers. If
#' reduction is not needed \code{resolution} can be set to 8. \code{resolution}
#' parameter can take following values: 2, 4, 6, 8. For more details
#' about HLA allele numbers resolution see
#' \url{http://hla.alleles.org/nomenclature/naming.html}.
#'
#' @inheritParams reduceAlleleResolution
#' @inheritParams utils::read.table
#' @param file Path to input file.
#'
#' @return Data frame containing HLA allele calls.
#'
#' \code{NA} values in input file are parsed unchanged.
#'
#' @examples
#' file <- system.file("extdata", "HLAHD_output_example.txt", package = "MiDAS")
#' hla_calls <- readHlaCalls(file)
#'
#' @importFrom assertthat assert_that is.readable see_if
#' @importFrom stats na.omit
#' @importFrom stringi stri_split_fixed
#' @importFrom utils read.table
#' @export
readHlaCalls <- function(file,
                         resolution = 4,
                         na.strings = c("Not typed", "-", "NA")) {
  assert_that(is.readable(file),
              is.count(resolution),
              is.character(na.strings)
  )
  hla_calls <- read.table(file,
                          header = TRUE,
                          sep = "\t",
                          stringsAsFactors = FALSE,
                          na.strings = na.strings
  )
  assert_that(checkHlaCallsFormat(hla_calls))

  # set colnames based on allele numbers
  gene_names <- vapply(X = 2:ncol(hla_calls),
                       FUN = function(i) {
                         names <- stri_split_fixed(hla_calls[, i], "*")
                         names <- vapply(X = names,
                                FUN = function(x) x[1],
                                FUN.VALUE = character(length = 1)
                         )
                         names <- unique(na.omit(names))
                         assert_that(
                           see_if(length(names) <= 1,
                                  msg = "Gene names in columns are not identical"
                           ),
                           see_if(length(names) != 0,
                                  msg = "One of the columns contains only NA"
                           )
                         )
                         return(names)
                       },
                       FUN.VALUE = character(length = 1)
  )
  ord <- order(gene_names)
  gene_names <- gene_names[ord]
  gene_names <- toupper(gene_names)
  gene_names_id <- unlist(lapply(table(gene_names), seq_len))
  gene_names <- paste(gene_names, gene_names_id, sep = "_")
  hla_calls <- hla_calls[, c(1, ord + 1)]
  colnames(hla_calls) <- c("ID", gene_names)

  hla_calls <- reduceHlaCalls(hla_calls, resolution = resolution)

  return(hla_calls)
}

#' Reads HLA allele alignments
#'
#' \code{readHlaAlignments} reads HLA allele alignments from file.
#'
#' HLA allele alignment file should follow format used EBI database, for details
#' see
#' \url{ftp://ftp.ebi.ac.uk/pub/databases/ipd/imgt/hla/alignments/README.md}.
#'
#' All protein alignment files from EBI database are shipped with the package.
#' They can be easily accessed using \code{gene} parameter. If \code{gene} is
#' set to \code{NULL} file parameter is used instead and alignment is read from
#' the provided file. In EBI database alignments for DRB1, DRB3, DRB4 and DRB5
#' genes are provided as a single file, here they are separated into separate
#' files for accordance with the package functionality.
#'
#' @inheritParams readHlaCalls
#' @inheritParams reduceAlleleResolution
#' @param gene Character vector of length one specifying the name of a gene for
#'   which alignment is required. See details for further explanations.
#' @param trim Logical indicating if alignment should be trimmed to start codon
#'   of the mature protein.
#' @param unkchar Character to be used to represent positions with unknown
#'   sequence.
#'
#' @return Matrix containing HLA allele alignments.
#'
#'   Rownames corresponds to allele numbers and columns to positions in the
#'   alignment. Sequences following the termination codon are marked as empty
#'   character (\code{""}). Unknown sequences are marked with a character of
#'   choice, by default \code{""}. Stop codons are represented by a hash (X).
#'   Insertion and deletions are marked with period (.).
#'
#' @examples
#' hla_alignments <- readHlaAlignments(gene = "A")
#'
#' @importFrom assertthat assert_that is.count is.readable is.string
#' @importFrom stringi stri_flatten stri_split_regex stri_sub
#' @importFrom stringi stri_subset_fixed stri_read_lines stri_detect_regex
#' @export
readHlaAlignments <- function(file,
                              gene = NULL,
                              trim = TRUE,
                              unkchar = "",
                              resolution = 8) {
  if (is.null(gene)) {
    assert_that(is.readable(file))
  } else {
    assert_that(is.string(gene))
    gene <- toupper(gene)
    file <- paste0(system.file("extdata", package = "MiDAS"),
                   "/",
                   gene,
                   "_prot.txt"
    )
    available_genes <- list.files(
      path = system.file("extdata", package = "MiDAS"),
      pattern = "_prot.txt$",
      full.names = TRUE
    )
    assert_that(
      see_if(file %in% available_genes,
             msg = sprintf("alignment for %s is not available", gene)
      ),
      is.readable(file)
    )
  }
  assert_that(
    isTRUEorFALSE(trim),
    is.string(unkchar),
    is.count(resolution)
  )

  aln_raw <- stri_read_lines(file)
  aln <- stri_split_regex(aln_raw, "\\s+")

  # extract lines containing alignments and omit empty alignment lines
  allele_lines <- vapply(X = aln,
                         FUN = function(x) any(checkAlleleFormat(x)),
                         FUN.VALUE = logical(length = 1)
  )
  assert_that(
    see_if(any(allele_lines),
           msg = "input file contains no correct HLA alignments"
    )
  )
  aln <- aln[allele_lines]
  aln <- aln[vapply(aln, length, integer(length = 1)) > 2]
  assert_that(
    see_if(all(stri_detect_regex(unlist(aln), "^[A-Z0-9:*.-]*$")),
                             msg = "alignments lines contain non standard characters"
    )
  )


  # convert aln into matrix and substitute for corresponding aa in ref
  allele_numbers <- vapply(X = aln,
                           FUN = function(x) x[2],
                           FUN.VALUE = character(length = 1)
  )
  aln <- vapply(X = aln,
                FUN = function(x) stri_flatten(x[-c(1, 2)]),
                FUN.VALUE = character(length = 1)
  )
  aln <- vapply(X = unique(allele_numbers),
                      FUN = function(x) {
                        stri_flatten(aln[allele_numbers == x])
                      },
                      FUN.VALUE = character(length = 1)
  )
  ref_seq <- stri_sub(aln[1],
                      seq(1, nchar(aln[1]), 1),
                      seq(1, nchar(aln[1]), 1)
  )
  aln <- do.call(rbind,
                 lapply(aln,
                        function(a) {
                          a <- stri_sub(a,
                                        seq(1, length(ref_seq), 1),
                                        seq(1, length(ref_seq), 1)
                          )
                          a[a == "-"] <- ref_seq[a == "-"]
                          a[a == "*"] <- unkchar
                          return(a)
                        }
                 )
  )

  # reduce alignment matrix to selected resolution
  allele_numbers <- reduceAlleleResolution(rownames(aln),
                                           resolution = resolution
  )
  unique_numbers <- ! duplicated(allele_numbers)
  aln <- aln[unique_numbers, ]
  rownames(aln) <- allele_numbers[unique_numbers]

  # discard aa '5 to start codon of mature protein
  if (trim) {
    raw_first_codon_idx <- nchar(stri_subset_fixed(aln_raw, "Prot")[1])
    raw_alignment_line <- stri_sub(aln_raw[allele_lines][1],
                                   1,
                                   raw_first_codon_idx
    )
    raw_alignment_seq <- stri_split_regex(raw_alignment_line, "\\s+")
    raw_alignment_seq <- unlist(raw_alignment_seq)[-c(1, 2)]
    first_codon_idx <- nchar(stri_flatten(raw_alignment_seq))
    assert_that(
      see_if(is.count(first_codon_idx),
             msg = "start codon is not marked properly in the input file"
      )
    )
    aln <- aln[, first_codon_idx:ncol(aln)]
  }

  return(aln)
}

#' Reads data table with KIR haplotypes calls
#'
#' \code{readKirCalls} reads table with KIR haplotypes calls from file.
#'
#' Input file have to be a tsv formatted table with two columns and header.
#' First column should contain samples IDs, second column should hold
#' corresponding KIR haplotypes.
#'
#' @inheritParams kirHaplotypeToCounts
#' @inheritParams utils::read.table
#' @param file Path to input file.
#' @param counts Logical flag indicating if KIR haplotypes should be converted
#'   to gene counts.
#'
#' @return Data frame containing KIR haplotypes calls or corresponding gene
#'   counts.
#'
#' @examples
#' file <- system.file("extdata", "KIP_output_example.txt", package = "MiDAS")
#' readKirCalls(file)
#'
#' @importFrom assertthat assert_that is.readable see_if
#' @importFrom dplyr left_join select
#' @importFrom stats na.omit setNames
#'
#' @export
readKirCalls <- function(file,
                         hap_dict = system.file("extdata", "Match_kir_haplotype_gene.txt", package = "MiDAS"),
                         counts = TRUE,
                         binary = TRUE,
                         na.strings = c("", "NA")) {
  assert_that(
    is.readable(file),
    is.readable(hap_dict),
    isTRUEorFALSE(counts),
    isTRUEorFALSE(binary),
    is.character(na.strings)
  )

  kir_calls <- read.table(file = file,
                          header = TRUE,
                          sep = "\t",
                          stringsAsFactors = FALSE,
                          na.strings = na.strings
  )
  assert_that(
    see_if(ncol(kir_calls) == 2,
           msg = sprintf(
             fmt = "KIR haplotypes calls table should have 2 columns, not %i",
             ncol(kir_calls)
           )
    ),
    see_if(
      all(is_hap <- grepl("^[0-9+|/]+$", na.omit(kir_calls[, 2, drop = TRUE]))),
      msg = sprintf(
        fmt = "rows %s of input file contains unexpected characters",
        paste(which(! is_hap), collapse = ", ")
      )
    )
  )

  if (counts) {
    haplotypes <- kir_calls[, 2, drop = TRUE]
    kir_counts <- kirHaplotypeToCounts(
      x = haplotypes,
      hap_dict = hap_dict,
      binary = binary
    )

    kir_calls <- left_join(
      x = kir_calls,
      y = kir_counts,
      by = setNames(colnames(kir_counts)[1], colnames(kir_calls)[2]),
      na_matches = "never"
    )

    # If there are multiple matches between x and y, all combinations of the matches are returned.
    kir_calls <- kir_calls[! duplicated(kir_calls[, 1, drop = TRUE]), ]
    rownames(kir_calls) <- NULL

    # discard haplotype designation from final table
    kir_calls <- select(kir_calls, - !!colnames(kir_calls)[2])
  }

  return(kir_calls)
}
