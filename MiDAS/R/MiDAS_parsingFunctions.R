#' Reads data table with HLA allele calls
#'
#' Reads data table with HLA allele calls from file in tsv format.
#'
#' @param file Path to the file containing HLA allele calls.
#' @param resolution Integer specifying the resolution to which the HLA allele
#'        calls should be reduced to. Valid values should be one of
#'        `2, 4, 6, 8`. To disable this functionality see `reduce` parameter.
#' @param reduce Logical flag specifying whether HLA numbers should be reduced
#'        to provided resolution.
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
                         reduce = TRUE) {
  assert_that(is.readable(file),
              is.count(resolution),
              is.flag(reduce)
  )
  hla_calls <- read.table(file,
                          header = TRUE,
                          sep = "\t",
                          stringsAsFactors = FALSE
  )
  assert_that(
    see_if(! all(checkAlleleFormat(hla_calls[, 1]), na.rm = TRUE),
           msg = "First column of input file should specify samples id"
    ),
    see_if(all(checkAlleleFormat(unlist(hla_calls[, -1])), na.rm = TRUE),
           msg = "Values in input file doesn't follow HLA numbers specification"
    )
  )

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

  if (reduce) {
    hla_calls[, -1] <- as.data.frame(
      lapply(hla_calls[, -1], reduceAlleleResolution, resolution = resolution),
      stringsAsFactors = FALSE
    )
  }

  return(hla_calls)
}

#' Reads HLA allele alignments
#'
#' Reads HLA allele alignments from file in msf format.
#'
#' @param file Path to the file containing HLA allele alignments.
#' @param gene Character vector of length one specifying the name of a gene for
#' which alignment is required. All the protein alignment files from EBI
#' database are shipped with the package and this parameter can be used to
#' provide simpler access to those files. If it's set to \code{NULL} file
#' parameter is used instead.
#' @param trim Logical indicating if alignment should be trimmed to start codon
#' of the mature protein.
#' @param unkchar Character to be used to represent positions with unknown
#' sequence.
#' @param resolution Integer specifying the resoultion with which alignment
#'        matrix should be returned.
#'
#' @return Matrix containing HLA allele alignments. Rownames corresponds to
#' allele numbers and columns to positions in the alignment. Sequences
#' following the termination codon are marked as empty character. Unknown
#' sequences are marked with a character of choice, that defaults to empty
#' character (""). Stop codons are represented by a hash (X). Insertion and
#' deletions are marked with period (.).
#'
#' @examples
#' hla_alignments <- readHlaAlignments(gene = "A")
#'
#' @importFrom assertthat assert_that is.count is.flag is.readable is.string
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
    is.flag(trim),
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
