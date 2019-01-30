#' Checks allele format
#'
#' @param allele Character vector containing HLA allele numbers.
#'
#' @return Logical vector of length 1 specifying if \code{allele} follows HLA
#'   alleles naming conventions and have desired resolution.
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
  pattern <- "^[A-Z0-9]+[*][0-9]+(:[0-9]+){0,3}[NLSCAQ]{0,1}$"
  is_correct <- stri_detect_regex(allele, pattern)
  return(is_correct)
}

#' Infers HLA allele resolution
#'
#' @inheritParams checkAlleleFormat
#'
#' @return Integer vector specifying alleles resolutions.
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
    see_if(all(checkAlleleFormat(allele), na.rm = TRUE),
         msg = "allele have to be a valid HLA allele number"
    )
  )
  allele_resolution <- 2 * (stri_count_fixed(allele, ":") + 1)
  return(allele_resolution)
}

#' Reduce HLA allele number
#'
#' @inheritParams checkAlleleFormat
#' @param resolution Numeric vector of length one specifying desired
#'   resolution. Note if `resolution` is greater than resolution of any
#'   \code{allele} elements, those elements will be returned unchanged.
#'
#' @return Character vector containing reduced HLA allele numbers, importantly
#'   alleles containing optional suffixes are omitted and returned unchanged.
#'
#' @examples
#' reduceAlleleResolution(c("A*01", "A*01:24", "C*05:24:55:54"), 2)
#'
#' @importFrom assertthat assert_that is.count see_if
#' @importFrom stringi stri_split_fixed stri_detect_regex
#' @export
reduceAlleleResolution <- function(allele,
                                   resolution=4) {
  assert_that(
    is.count(resolution)
  )
  na_idx <- is.na(allele)
  letter_alleles <- stri_detect_regex(allele, pattern = "[A-Z]{1}$")
  allele_res <- getAlleleResolution(allele)
  to_reduce <- allele_res >= resolution & ! letter_alleles & ! na_idx
  resolution <- floor(resolution) / 2
  allele[to_reduce] <- vapply(
    X = stri_split_fixed(allele[to_reduce], ":"),
    FUN = function(a) {
      paste(a[1:resolution], collapse = ":")
    },
    FUN.VALUE = character(length = 1)
  )
  allele[na_idx] <- NA
  return(allele)
}

#' Returns positions of variable amino acids in the alignment.
#'
#' @param alignment Matrix containing amino acids level alignment.
#' @param varchar Regex matching characters that should be considered when
#'   looking for variable amino acids positions. Eg. when varchar = "[A-Z]"
#'   occurence of deletion/insertion (".") will not be treated as variability.
#'   In order to detect this kind of variability \code{varchar = "[A-Z\\.]"}
#'   should be used.
#'
#' @return Integer vector specifying which alignment columns are variable.
#'
#' @examples
#' file <- system.file("extdata", "A_prot.txt", package = "MiDAS")
#' alignment <- readHlaAlignments(file)
#' getVariableAAPos(alignment)
#'
#' @importFrom assertthat assert_that is.count see_if
#' @export
getVariableAAPos <- function(alignment,
                             varchar = "[A-Z]") {
  assert_that(is.matrix(alignment))
  var_cols <- apply(alignment,
        2,
        function(col) {
          col <- col[grepl(sprintf("^%s$", varchar), col)]
          is.variable <- logical(length = 1)
          if (length(col) == 0) {
            is.variable <- FALSE
          } else {
            is.variable <- any(col != col[1])
          }
          return(is.variable)
        }
  )
  return(which(var_cols))
}

#' Converts allele numbers to additional variables based on match file.
#'
#' @inheritParams checkAlleleFormat
#' @param dictionary Path to the file containing HLA allele numbers matchings.
#'   The file should be in tsv format with header and two columns. First column
#'   should hold allele numbers and second corresponding additional variables.
#'   Function also accepts match table in form of data frame.
#'
#' @return Vector containing allele numbers converted to additional variables
#'   according to matching file.
#'
#' @examples
#' dictionary <- system.file("extdata", "Match_4digit_supertype.txt", package = "MiDAS")
#' convertAlleleToVariable(c("A*01:01", "A*02:01"), dictionary = dictionary)
#'
#' @importFrom assertthat assert_that is.readable see_if
#' @importFrom stats setNames
#' @importFrom utils type.convert
#' @export
convertAlleleToVariable <- function(allele,
                                    dictionary) {
  assert_that(
    see_if(all(checkAlleleFormat(allele), na.rm = TRUE),
           msg = "allele have to be a valid HLA allele number"
    ),
    see_if(is.string(dictionary) | is.data.frame(dictionary),
           msg = "dictionary have to be either path or data.frame"
    )
  )
  if (is.character(dictionary)) {
    assert_that(
      is.readable(dictionary)
    )
    dictionary <- read.table(
      file = dictionary,
      header = TRUE,
      stringsAsFactors = FALSE
    )
  }
  assert_that(
    see_if(ncol(dictionary) == 2,
           msg = "match table have to consist out of two columns"
    ),
    see_if(all(checkAlleleFormat(dictionary[, 1]), na.rm = TRUE),
           msg = "first column of match table must contain valid HLA allele numbers"
    ),
    see_if(! any(duplicated(dictionary[, 1]), na.rm = TRUE),
           msg = "match table contains duplicated allele numbers")
  )
  dictionary <- setNames(dictionary[, 2], dictionary[, 1])
  variable <- dictionary[allele]
  variable <- type.convert(variable, as.is = TRUE)
  return(variable)
}

#' Assert hla calls data frame format
#'
#' @param hla_calls Data frame containing HLA allele calls, in a format as
#'   return by \code{\link{readHlaCalls}} function.
#'
#' @return Logical indicating if hla_calls follows hla calls data frame format.
#'   Otherwise raise error.
#'
#' @examples
#' file <- system.file("extdata", "HLAHD_output_example.txt", package = "MiDAS")
#' hla_calls <- readHlaCalls(file)
#' checkHlaCallsFormat(hla_calls)
#'
#' @export
checkHlaCallsFormat <- function(hla_calls) {
  assert_that(
    see_if(! any(vapply(hla_calls, is.factor, logical(length = 1))),
           msg = "input can't contain factors"
    ),
    see_if(! all(checkAlleleFormat(hla_calls[, 1]), na.rm = TRUE),
           msg = "first column of input should specify samples id"
    ),
    see_if(all(checkAlleleFormat(unlist(hla_calls[, -1])), na.rm = TRUE),
           msg = "values in input doesn't follow HLA numbers specification"
    )
  )

  return(TRUE)
}
