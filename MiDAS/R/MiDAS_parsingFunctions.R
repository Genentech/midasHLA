#' Reads data table with HLA allele calls
#'
#' Reads data table with HLA allele calls from file in tsv format.
#'
#' @param file Path to the file containing HLA allele calls.
#'
#' @return Data frame containing HLA allele calls.
#'
#' @example
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
readHlaAlignments <- function(file,
                              trim = TRUE) {
  if (! file.exists(file)) {
    stop("Error: File ", file, "doesn't exist.")
  }
  aln_raw <- stringi::stri_read_lines(file)

  #TODO check if file follows alignment format
  # header <- stringi::stri_subset_regex(aln_raw, "^#")
  # !!!!!!!!!!!
  # instead check various semi products


  # extract lines containing alignments
  aln_split <- stringi::stri_split_regex(aln_raw, "\\s+")
  allele_lines <- vapply(X = aln_split,
                         FUN = function(x) any(checkAlleleFormat(x)),
                         FUN.VALUE = logical(length = 1)
  )
  aln_split <- aln_split[allele_lines]

  # locate index of start codon of the mature protein
  first_codon_idx <- nchar(stringi::stri_subset_fixed(aln_raw, "Prot")[1])
  first_codon_idx <- substr(aln_raw[allele_lines][1], 1, first_codon_idx)
  first_codon_idx <- unlist(stringi::stri_split_regex(first_codon_idx,
                                                      "\\s+"
  )
  )[-c(1, 2)]
  first_codon_idx <- nchar(paste(first_codon_idx, collapse = ""))

  #
  allele_numbers <- vapply(X = aln_split,
                           FUN = function(x) x[2],
                           FUN.VALUE = character(length = 1)
  )
  aln_split <- vapply(X = aln_split,
                      FUN = function(x) stringi::stri_flatten(x[-c(1, 2)]),
                      FUN.VALUE = character(length = 1)
  )
  aln_split <- vapply(X = unique(allele_numbers),
                      FUN = function(x) {
                        stringi::stri_flatten(aln_split[allele_numbers == x])
                      },
                      FUN.VALUE = character(length = 1)
  )

  # convert to matrix
  aln_split <- do.call(rbind,
                       lapply(aln_split,
                              function(a) {
                                stringi::stri_sub(a,
                                                  seq(1, nchar(aln_split[1]), 1),
                                                  seq(1, nchar(aln_split[1]), 1)
                                )
                              }
                       )
  )

  # substitute - for corresponding AA in reference sequnce
  hla_alignments <- t(apply(X = aln_split,
                            MARGIN = 1,
                            FUN = function(a) {
                              a[a == "-"] <- aln_split[1, ][a == "-"]
                              return(a)
                            }
  )
  ) # mby this could be better solved by extracting ref seq beforehand. Test it after wrapping on function

  # discard aa '5 to start codon of mature protein
  if (trim) {
    hla_alignments <- hla_alignments[, first_codon_idx:ncol(hla_alignments)]
  }

  return(hla_alignments)
}
