#' Checks allele format
#'
#' Checks correctnes of the HLA allele number format and resolution.
#'
#' @param allele Character vector containing HLA allele numbers.
#' @param resolution Integer vector of length one specifiying desired
#' resolution.
#'
#' @return Logical vector of length specifying if \code{allele} follows HLA
#' alleles naming conventions and have desired resolution.
checkAlleleFormat <- function(allele, resolution=4) {
  return(isCorrect)
}

#' Infers HLA allele resolution
#'
#' Infers HLA allele resolution.
#'
#' @param allele Character vector containing HLA allele numbers.
#'
#' @return Integer vector specifiying alleles resolutions.
getAlleleResolution <- function(allele) {
  return(allele_resolution)
}

#' Reduce HLA allele number
#'
#' Reduces HLA allele number to lower level representation.
#'
#' @param allele Character vector containing HLA allele numbers.
#' @param resolution Integer vector of length one specifiying desired
#' resolution. Cannot be greater than \code{allele} resolution.
#'
#' @return Character vector containing reduced HLA allele numbers.
reduceAlleleResolution <- function(allele, resolution=4) {
  return(reduced_allele)
}
