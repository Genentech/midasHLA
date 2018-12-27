#' Reads data table with HLA allele calls
#'
#' Reads data table with HLA allele calls from file in tsv format.
#'
#' @param file Path to the file containing HLA allele calls.
#'
#' @return Data frame containing HLA allele calls.
#'
#' @examples
#' file <- system.file("extdata/HLAHD_output_example.txt", package = "MiDAS")
#' hla_calls <- readHlaCalls(file)
readHlaCalls <- function(file,
                         sep = "\t",
                         header = TRUE,
                         verbose = TRUE) {
  if (! file.exists(file)) {
    stop("Error: File ", file, "doesn't exist.")
  }
  hla_calls <- read.table(file,
                          header = header,
                          sep = sep,
                          stringsAsFactors = FALSE
  )
  if (all(checkAlleleFormat(hla_calls[, 1]), na.rm = TRUE)) {
    stop("Error: First column of input file should specify samples id.")
  }
  if (! all(checkAlleleFormat(unlist(hla_calls[, -1])), na.rm = TRUE)) {
    stop("Error: Values in input file doesn't follow HLA numbers specification.")
  }

  # set colnames based on allele numbers
  gene_names <- vapply(X = 2:ncol(hla_calls),
                       FUN = function(i) {
                         names <- stringi::stri_split_fixed(hla_calls[, i], "*")
                         names <- vapply(X = names,
                                FUN = function(x) x[1],
                                FUN.VALUE = character(length = 1)
                         )
                         names <- unique(na.omit(names))
                         if (length(names) > 1) {
                           stop("Error: Gene names in columns are not identical.")
                         }
                         if (length(names) == 0) {
                           stop("Error: One of the columns contains only NA.")
                         }
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
#'
#' @return Matrix containing HLA allele alignments.
#'
#' @importFrom assertthat assert_that is.readable
#' @importFrom stringi stri_flatten stri_split_regex stri_sub
#' @importFrom stringi stri_subset_fixed stri_read_lines
readHlaAlignments <- function(file,
                              trim = TRUE) {
  assert_that(is.readable(file)) # TODO check if file follows alignment format
  aln_raw <- stri_read_lines(file)
  aln <- stri_split_regex(aln_raw, "\\s+")

  # extract lines containing alignments
  allele_lines <- vapply(X = aln,
                         FUN = function(x) any(checkAlleleFormat(x)),
                         FUN.VALUE = logical(length = 1)
  )
  aln <- aln[allele_lines]

  # locate index of start codon of the mature protein
  raw_first_codon_idx <- nchar(stri_subset_fixed(aln_raw, "Prot")[1])
  raw_alignment_line <- stri_sub(aln_raw[allele_lines][1],
                                 1,
                                 raw_first_codon_idx
  )
  raw_alignment_seq <- unlist(stri_split_regex(raw_alignment_line,
                                               "\\s+"
                              )
  )[-c(1, 2)]
  first_codon_idx <- nchar(stri_flatten(raw_alignment_seq))

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
  ref_seq <- stri_sub(aln[1], seq(1, nchar(aln[1]), 1), seq(1, nchar(aln[1]), 1))
  aln <- do.call(rbind,
                 lapply(aln,
                        function(a) {
                          a <- stri_sub(a,
                                        seq(1, length(ref_seq), 1),
                                        seq(1, length(ref_seq), 1)
                          )
                          a[a == "-"] <- ref_seq[a == "-"]
                          a[a == "*"] <- ""
                          return(a)
                        }
                 )
  )

  # discard aa '5 to start codon of mature protein
  if (trim) {
    aln <- aln[, first_codon_idx:ncol(aln)]
  }

  return(aln)
}
