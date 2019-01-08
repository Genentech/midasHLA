#' Reads data table with HLA allele calls
#'
#' Reads data table with HLA allele calls from file in tsv format.
#'
#' @param file Path to the file containing HLA allele calls.
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
readHlaCalls <- function(file) {
  assert_that(is.readable(file))
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

  return(hla_calls)
}

#' Reads HLA allele alignments
#'
#' Reads HLA allele alignments from file in msf format.
#'
#' @param file Path to the file containing HLA allele alignments.
#' @param trim Logical indicating if alignment should be trimmed to start codon
#' of the mature protein.
#' @param unkchar Character to be used to represent positions with unknown
#' sequence.
#'
#' @return Matrix containing HLA allele alignments. Rownames corresponds to
#' allele numbers and columns to positions in the alignment. Sequences
#' following the termination codon are marked as empty character. Unknown
#' sequences are marked with a character of choice, that defaults to empty
#' character (""). Stop codons are represented by a hash (X). Insertion and
#' deletions are marked with period (.).
#'
#' @examples
#' file <- system.file("extdata", "A_prot.txt", package = "MiDAS")
#' hla_alignments <- readHlaAlignments(file)
#'
#' @importFrom assertthat assert_that is.count is.readable
#' @importFrom stringi stri_flatten stri_split_regex stri_sub
#' @importFrom stringi stri_subset_fixed stri_read_lines stri_detect_regex
#' @export
readHlaAlignments <- function(file,
                              trim = TRUE,
                              unkchar = "") {
  assert_that(is.readable(file))
  aln_raw <- stri_read_lines(file)
  aln <- stri_split_regex(aln_raw, "\\s+")

  # extract lines containing alignments
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
    assert_that(is.count(first_codon_idx))
    aln <- aln[, first_codon_idx:ncol(aln)]
  }

  return(aln)
}
