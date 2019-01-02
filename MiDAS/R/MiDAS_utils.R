#' Checks allele format
#'
#' @param allele Character vector containing HLA allele numbers.
#'
#' @return Logical vector of length specifying if \code{allele} follows HLA
#' alleles naming conventions and have desired resolution.
#'
#' @examples
#' allele <- c("A*01:01", "A*01:02")
#' checkAlleleFormat(allele)
#'
#' @importFrom assertthat assert_that
#' @importFrom stringi stri_detect_regex
#' @export
checkAlleleFormat <- function(allele) {
  assert_that(is.character(allele))
  pattern <- "^[A-Z]+[*][0-9]+(:[0-9]+){0,3}[NLSCAQ]{0,1}$"
  is_correct <- stri_detect_regex(allele, pattern)
  return(is_correct)
}

#' Infers HLA allele resolution
#'
#' @param allele Character vector containing HLA allele numbers.
#'
#' @return Integer vector specifiying alleles resolutions.
#'
#' @examples
#' allele <- c("A*01:01", "A*01:02")
#' getAlleleResolution(allele)
#'
#' @importFrom assertthat assert_that see_if
#' @importFrom stringi stri_count_fixed
#' @export
getAlleleResolution <- function(allele) {
  assert_that(
    see_if(all(checkAlleleFormat(allele)),
         msg = "allele have to be a valid HLA allele number"
    )
  )
  allele_resolution <- 2 * (stri_count_fixed(allele, ":") + 1)
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
#' @examples
#' reduceAlleleResolution(c("A*01", "A*01:24", "C*05:24:55:54"), 2)
#'
#' @importFrom assertthat assert_that is.count see_if
#' @importFrom stringi stri_split_fixed
#' @export
reduceAlleleResolution <- function(allele,
                                   resolution=4) {
  assert_that(
    is.count(resolution),
    see_if(all(getAlleleResolution(allele) >= resolution),
           msg = "input resolution can't be lower than requested resolution"
    )
  )
  resolution <- floor(resolution) / 2
  reduced_allele <- vapply(X = stri_split_fixed(allele, ":"),
                           FUN = function(a) {
                             paste(a[1:resolution], collapse = ":")
                           },
                           FUN.VALUE = character(length = 1)
                          )
  return(reduced_allele)
}
