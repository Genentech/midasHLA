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
  pattern <- "^[A-Z]+[*][0-9]+(:[0-9]+){0,3}[NLSCAQ]{0,1}$"
  is_correct <- stringi::stri_detect_regex(allele, pattern)
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
  if (! all(checkAlleleFormat(allele))) {
    stop("Error: Allele have to be a valid HLA allele number.")
  }
  allele_resolution <- 2 * (stringi::stri_count_fixed(allele, ":") + 1)
  return(allele_resolution)
}

#' Reduce HLA allele number
#'
#' @param allele Character vector containing HLA allele numbers.
#' @param resolution Numeric vector of length one specifiying desired
#' resolution. Cannot be greater than \code{allele} resolution.
#'
#' @return Character vector containing reduced HLA allele numbers.
#'
#' @example
#' reduceAlleleResolution(c("A*01", "A*01:24", "C*05:24:55:54"), 2)
reduceAlleleResolution <- function(allele,
                                   resolution=4) {
  if (! is.numeric(resolution) && length(resolution) == 1) {
    stop("Error: Resolution have to be of type numeric with length one.")
  }
  if (any(getAlleleResolution(allele) < resolution)) {
    stop("Error: Input resolution can't be lower than requested resolution.")
  }
  resolution <- floor(resolution) / 2
  reduced_allele <- vapply(X = stringi::stri_split_fixed(allele, ":"),
                           FUN = function(a) {
                             paste(a[1:resolution], collapse = ":")
                           },
                           FUN.VALUE = character(length = 1)
                          )
  return(reduced_allele)
}
