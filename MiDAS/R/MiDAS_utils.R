#' Check allele format
#'
#' \code{checkAlleleFormat} test if the input character follows HLA nomenclature
#' specifications.
#'
#' Correct HLA number should consist of HLA gene name followed by "*" and sets
#' of digits separated with ":". Maximum number of sets of digits is 4 which
#' is termed 8-digit resolution. Optionally HLA numbers can be supplemented with
#' additional suffix indicating its expression status. See
#' \url{http://hla.alleles.org/nomenclature/naming.html} for more details.
#'
#' HLA alleles with identical sequences across exons encoding the peptide
#' binding domains might be designated with G group allele numbers. Those
#' numbers have additional G or GG suffix. See
#' \url{http://hla.alleles.org/alleles/g_groups.html} for more details.
#'
#' @param allele Character vector containing HLA allele numbers.
#'
#' @return Logical vector specifying if \code{allele} follows HLA alleles naming
#'   conventions.
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
  pattern <- "^[A-Z0-9]+[*][0-9]+(:[0-9]+){0,3}((?=G)(G|GG)|(N|L|S|C|A|Q)){0,1}$"
  is_correct <- stri_detect_regex(allele, pattern)
  return(is_correct)
}

#' Infers HLA allele resolution
#'
#' \code{getAlleleResolution} returns the resolution of input HLA allele
#' numbers.
#'
#' HLA allele resolution can take following values: 2, 4, 6, 8.See
#' \url{http://hla.alleles.org/nomenclature/naming.html} for more details.
#'
#' @inheritParams checkAlleleFormat
#'
#' @return Integer vector specifying alleles resolutions.
#'
#'   \code{NA} values are accepted and returned as \code{NA}.
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

#' Reduce HLA allele resolution
#'
#' \code{reduceAlleleResolution} reduces HLA allele numbers vector to specified
#' resolution.
#'
#' In cases when allele numbers contains additional suffix their resolution
#' can not be unambiguously reduced. These cases are returned unchanged.
#' Function behaves in the same manner if \code{resolution} is higher than
#' resolution of input HLA allele numbers.
#'
#' @inheritParams checkAlleleFormat
#' @param resolution Numeric vector of length one specifying output resolution.
#'
#' @return Character vector containing reduced HLA allele numbers.
#'
#'   \code{NA} values are accepted and returned as \code{NA}.
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
  letter_alleles <- stri_detect_regex(allele, pattern = "(N|L|S|C|A|Q){1}$")
  is_ggroup <- stri_detect_regex(allele, pattern = "(G|GG){1}$")
  allele_res <- getAlleleResolution(allele)
  to_reduce <- allele_res > resolution & ! letter_alleles & ! na_idx
  resolution <- floor(resolution) / 2
  if (any(is_ggroup & to_reduce)) {
    warning("Reducing G groups alleles, major allele gene name will be used.")
  }
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

#' Returns positions of variable amino acids in the alignment
#'
#' \code{getVariableAAPos} finds variable amino acids positions in the
#' alignment.
#'
#' The variable amino acid positions in the alignment are those at which
#' different amino acids can be found. As the alignments can also contain indels
#' and unknown characters, the user choice might be to consider those positions
#' also as variable. This can be achieved by passing appropriate regular
#' expression in \code{varchar}. Eg. when \code{varchar = "[A-Z]"} occurence of
#' deletion/insertion (".") will not be treated as variability. In order to
#' detect this kind of variability \code{varchar = "[A-Z\\\\.]"} should be used.
#'
#' @param alignment Matrix containing amino acids level alignment.
#' @param varchar Regex matching characters that should be considered when
#'   looking for variable amino acids positions. See details for further
#'   explanations.
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

#' Converts allele numbers to additional variables
#'
#' \code{convertAlleleToVariable} convert input HLA allele numbers to additional
#' variables based on the supplied match table (dictionary).
#'
#' \code{dictionary} file should be a tsv format with header and two columns.
#' First column should hold allele numbers and second corresponding additional
#' variables. Optionally a data frame formatted in the same manner can be passed
#' insted.
#'
#' @inheritParams checkAlleleFormat
#' @param dictionary Path to the file containing HLA allele numbers matchings or
#'   data frame providing this information. See details for further
#'   explanations.
#'
#' @return Vector containing HLA allele numbers converted to additional
#'   variables according to matching file.
#'
#'   Type of the returned vector depends on the type of the additional variable.
#'   For example for HLA alleles supertypes character vector is returned while
#'   for expression values numeric vector will be returned.
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
#' \code{checkHlaCallsFormat} asserts if hla calls data frame have proper
#' format.
#'
#' @param hla_calls Data frame containing HLA allele calls, as return by
#'   \code{\link{readHlaCalls}} function.
#'
#' @return Logical indicating if \code{hla_calls} follows hla calls data frame
#'   format. Otherwise raise error.
#'
#' @importFrom assertthat assert_that see_if
#' @examples
#' file <- system.file("extdata", "HLAHD_output_example.txt", package = "MiDAS")
#' hla_calls <- readHlaCalls(file)
#' checkHlaCallsFormat(hla_calls)
#'
#' @export
checkHlaCallsFormat <- function(hla_calls) {
  assert_that(
    is.data.frame(hla_calls),
    see_if(nrow(hla_calls) >= 1 & ncol(hla_calls) >= 2,
           msg = "hla_calls have to have at least 1 rows and 2 columns"
    ),
    see_if(! any(vapply(hla_calls, is.factor, logical(length = 1))),
           msg = "hla_calls can't contain factors"
    ),
    see_if(! all(checkAlleleFormat(hla_calls[, 1]), na.rm = TRUE),
           msg = "first column of hla_calls should specify samples id"
    ),
    see_if(all(checkAlleleFormat(unlist(hla_calls[, -1])), na.rm = TRUE),
           msg = "values in hla_calls doesn't follow HLA numbers specification"
    )
  )

  return(TRUE)
}

#' Backquote string
#'
#' \code{backquote} places backticks around string.
#'
#' \code{backquote} is usefull when using HLA allele numbers in fomulas, where
#' \code{'*'} and \code{':'} characters have special meanings.
#'
#' @param x Character vector.
#'
#' @return Character vector with its elements backticked.
#'
#' @examples
#' backquote("A*01:01")
#'
#' @importFrom assertthat assert_that
#' @export
backquote <- function(x) {
  assert_that(is.character(x))
  backquoted <- paste0("`", x, "`")
  return(backquoted)
}
