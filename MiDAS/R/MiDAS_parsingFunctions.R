#' Reads data table with HLA allele calls
#'
#' Reads data table with HLA allele calls from file in tsv format.
#'
#' @param file Path to the file containing HLA allele calls.
#'
#' @return Data frame containing HLA allele calls.
readHlaCalls <- function(file) {
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
readHlaAlignments <- function(file) {
  return(hla_alignments)
}
