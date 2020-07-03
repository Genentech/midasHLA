#' Read HLA allele calls data
#'
#' \code{readHlaCalls} read HLA allele calls from file
#'
#' Input file have to be a tsv formatted table with header. First column should
#' contain sample IDs, further columns should hold corresponding HLA allele
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

#' Read HLA allele alignments
#'
#' \code{readHlaAlignments} read HLA allele alignments from file.
#'
#' HLA allele alignment file should follow EBI database format, for details
#' see
#' \url{ftp://ftp.ebi.ac.uk/pub/databases/ipd/imgt/hla/alignments/README.md}.
#'
#' All protein alignment files from EBI database are shipped with the package.
#' They can be easily accessed using \code{gene} parameter. If \code{gene} is
#' set to \code{NULL} file parameter is used instead and alignment is read from
#' the provided file. In EBI database alignments for DRB1, DRB3, DRB4 and DRB5
#' genes are provided as a single file, here they are separated.
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
                              trim = FALSE,
                              unkchar = "",
                              resolution = 8) {
  assert_that(
    isTRUEorFALSE(trim),
    is.string(unkchar),
    is.count(resolution)
  )

  if (is.null(gene)) {
    assert_that(is.readable(file))
    aln_raw <- stri_read_lines(file)
    aln <- stri_split_regex(aln_raw, "\\s+")

    # extract lines containing alignments and omit empty alignment lines
    nonempty_lines <- vapply(aln, length, integer(length = 1)) >= 2
    aln <- aln[nonempty_lines]
    allele_numbers <- vapply(aln, `[`, character(length = 1), 2)
    allele_lines <- checkAlleleFormat(allele_numbers)

    assert_that(
      see_if(any(allele_lines),
             msg = "could not find alleles numbers in the alignment file"
      )
    )
    allele_numbers <- allele_numbers[allele_lines]
    aln <- aln[allele_lines]
    assert_that(
      see_if(all(stri_detect_regex(unlist(aln), "^[A-Z0-9:*.-]*$")),
                               msg = "alignments lines contain non standard characters"
      )
    )

    tmp_aln_env <- new.env(size = 5000)
    for (i in 1:length(allele_numbers)) {
      assign(
        x = allele_numbers[i],
        value = append(
          x = get0(
            x = allele_numbers[i],
            envir = tmp_aln_env,
            ifnotfound = character(length = 0)
          ),
          values = aln[[i]][-c(1, 2)] # discard empty element and allele number
        ),
        envir = tmp_aln_env
      )
    }
    aln_list <- as.list(tmp_aln_env)[unique(allele_numbers)] # convert to list and sort

    ref_seq <- stri_flatten(aln_list[[1]])
    seq_along_ref <- seq(1, nchar(ref_seq), 1)
    ref_seq <- stri_sub(ref_seq,
                        seq_along_ref,
                        seq_along_ref
    )
    ref_len <- length(ref_seq)
    aln <- do.call(rbind,
                   lapply(aln_list,
                          function(a) {
                            a <- stri_flatten(a)
                            a <- stri_sub(a,
                                          seq_along_ref,
                                          seq_along_ref
                            )
                            i <- a == "-"
                            a[i] <- ref_seq[i]

                            return(a)
                          }
                   )
    )

    # find AA positions numbers
    aln_raw <- aln_raw[nonempty_lines]
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
    if (first_codon_idx > 1) {
      aln_colnames <- c(seq(1 - first_codon_idx, -1, 1),
                        seq(1, ncol(aln) + 1 - first_codon_idx, 1)
                      )
    } else {
      aln_colnames <- seq(1, ncol(aln) + 1 - first_codon_idx, 1)
    }
    colnames(aln) <- aln_colnames

    # discard aa '5 to start codon of mature protein
    if (trim) {
      aln <- aln[, first_codon_idx:ncol(aln)]
    }
  } else {
    assert_that(is.string(gene))
    gene <- toupper(gene)
    file <- paste0(system.file("extdata", package = "MiDAS"),
                   "/",
                   gene,
                   "_prot.Rdata"
    )
    available_genes <- list.files(
      path = system.file("extdata", package = "MiDAS"),
      pattern = "_prot.Rdata$",
      full.names = TRUE
    )
    assert_that(
      see_if(file %in% available_genes,
             msg = sprintf("alignment for %s is not available", gene)
      )
    )

    cached_aln_obj <- readRDS(file) # list(readHlaAlignments(file, trim = FALSE, unkchar = "*", resolution = 8), first_codon_idx)
    aln <- cached_aln_obj[[1]]

    # discard aa '5 to start codon of mature protein
    if (trim) {
      first_codon_idx <- cached_aln_obj[[2]]
      aln <- aln[, first_codon_idx:ncol(aln)]
    }
  }

  # substitute unkchar
  aln[aln == "*"] <- unkchar

  # reduce alignment matrix to selected resolution
  allele_numbers <- reduceAlleleResolution(rownames(aln),
                                           resolution = resolution
  )
  unique_numbers <- ! duplicated(allele_numbers)
  aln <- aln[unique_numbers, , drop = FALSE]
  rownames(aln) <- allele_numbers[unique_numbers]

  return(aln)
}

#' Reads data table with KIR haplotypes calls
#'
#' \code{readKPICalls} reads table with KIR haplotypes calls from file.
#'
#' Input file have to be a tsv formatted table with two columns and header.
#' First column should contain samples IDs, second column should hold
#' corresponding KIR haplotypes.
#'
#' @inheritParams utils::read.table
#'
#' @return Data frame containing KIR gene counts.
#'
#' @examples
#' file <- system.file("extdata", "KPI_output_example.txt", package = "MiDAS")
#' readKPICalls(file)
#'
#' @importFrom assertthat assert_that is.readable see_if
#' @importFrom dplyr left_join select
#' @importFrom stats na.omit setNames
#'
#' @export
readKPICalls <- function(file,
                         na.strings = c("", "NA", "uninterpretable")) {
  assert_that(
    is.readable(file),
    is.character(na.strings)
  )

  kpi_output <- read.table(
    file = file,
    header = TRUE,
    sep = "\t",
    stringsAsFactors = FALSE,
    na.strings = na.strings
  )
  kir_calls <- kpi_output[, -2, drop = FALSE]
  checkKirCallsFormat(kir_calls)

  return(kir_calls)
}
