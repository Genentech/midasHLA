#' Checks allele format
#'
#' @param allele Character vector containing HLA allele numbers.
#' @param resolution Integer vector of length one specifiying desired
#' resolution.
#'
#' @return Logical vector of length specifying if \code{allele} follows HLA
#' alleles naming conventions and have desired resolution.
#'
#' @example
#' allele <- c("A*01:01", "A*01:02")
#' checkAlleleFormat(allele)
checkAlleleFormat <- function(allele) {
  if (! typeof(allele) == "character") {
    stop("Error: Allele have to be of type character.")
  }
  pattern <- ".*[*][0-9][0-9](:[0-9][0-9])*[NLSCAQ]{0,1}$"
  is_correct <- stringi::stri_detect(pattern, allele)
  return(is_correct)
}

#' Infers HLA allele resolution
#'
#' @param allele Character vector containing HLA allele numbers.
#'
#' @return Integer vector specifiying alleles resolutions.
#'
#' @example
#' allele <- c("A*01:01", "A*01:02")
#' getAlleleResolution(allele)
getAlleleResolution <- function(allele) {
  return(allele_resolution)
}

#' Reduce HLA allele number
#'
#' @param allele Character vector containing HLA allele numbers.
#' @param resolution Integer vector of length one specifiying desired
#' resolution. Cannot be greater than \code{allele} resolution.
#'
#' @return Character vector containing reduced HLA allele numbers.
reduceAlleleResolution <- function(allele,
                                   resolution=4) {
  return(reduced_allele)
}
