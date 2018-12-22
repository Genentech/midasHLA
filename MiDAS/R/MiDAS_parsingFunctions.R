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

file <- system.file("extdata/HLAHD_output_example.txt", package = "MiDAS")
hla_calls <- readHlaCalls(file)

#' Reads HLA allele alignments
#'
#' Reads HLA allele alignments from file in msf format.
#'
#' @param file Path to the file containing HLA allele alignments.
#'
#' @return Matrix containing HLA allele alignments.
#'
readHlaAlignments <- function(file) {
  return(hla_alignments)
}
