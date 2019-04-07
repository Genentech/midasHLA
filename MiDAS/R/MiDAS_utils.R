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
#' @importFrom rlang warn
#' @export
reduceAlleleResolution <- function(allele,
                                   resolution=4) {
  assert_that(
    see_if(all(checkAlleleFormat(allele), na.rm = TRUE),
           msg = "allele have to be a valid HLA allele number"
    ),
    is.count(resolution)
  )
  na_idx <- is.na(allele)
  letter_alleles <- stri_detect_regex(allele, pattern = "(N|L|S|C|A|Q){1}$")
  is_ggroup <- stri_detect_regex(allele, pattern = "(G|GG){1}$")
  allele_res <- getAlleleResolution(allele)
  to_reduce <- allele_res > resolution & ! letter_alleles & ! na_idx
  resolution <- floor(resolution) / 2
  if (any(is_ggroup & to_reduce)) {
    warn("Reducing G groups alleles, major allele gene name will be used.")
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
    see_if(! all(checkAlleleFormat(as.character(hla_calls[, 1])), na.rm = TRUE),
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
  x <- gsub("`", "", x)
  backquted <- paste0("`", x, "`")
  return(backquted)
}

#' Assert additional data
#'
#' \code{checkAdditionalData} asserts if phenotype or covariate data frame
#' has proper format.
#'
#' @inheritParams checkHlaCallsFormat
#' @param data_frame Data frame containing phenotype or covariate data
#'   corresponding to accompanying hla calls data frame.
#' @param accept.null Logical indicating if NULL data_frame should be accepted.
#'
#' @return Logical indicating if \code{data_frame} is properly formatted.
#'   Otherwise raise error.
#'
#' @importFrom assertthat assert_that see_if
#' @examples
#' hla_file <- system.file("extdata", "HLAHD_output_example.txt", package = "MiDAS")
#' pheno_file <- system.file("extdata", "pheno_example.txt", package = "MiDAS")
#' pheno <- read.table(pheno_file, header = TRUE)
#' hla_calls <- readHlaCalls(hla_file)
#' checkAdditionalData(pheno, hla_calls)
#'
#' @export
checkAdditionalData <- function(data_frame,
                                hla_calls,
                                accept.null = FALSE) {
  data_frame_name <- deparse(substitute(data_frame))
  hla_calls_name <- deparse(substitute(hla_calls))
  if (! (is.null(data_frame) & accept.null)) {
    assert_that(
      see_if(is.data.frame(data_frame),
             msg = sprintf("%s have to be a data frame",
                           data_frame_name
             )
      ),
      checkHlaCallsFormat(hla_calls),
      see_if(nrow(data_frame) >= 1 & ncol(data_frame) >= 2,
             msg = sprintf("%s have to have at least 1 rows and 2 columns",
                           data_frame_name
             )
      ),
      see_if(colnames(data_frame)[1] == colnames(hla_calls)[1],
             msg = sprintf(
               "first column in %s must be named as first column in %s",
               data_frame_name, hla_calls_name
             )
      ),
      see_if(any(hla_calls[, 1] %in% data_frame[, 1]),
             msg = sprintf(
               "IDs in %s doesn't match IDs in %s",
               data_frame_name, hla_calls_name
             )
      )
    )
  }

  return(TRUE)
}

#' Add new variables to statistical model
#'
#' \code{updateModel} will add new variables to model and re-fit it.
#'
#' @param object An existing fit from a model function such as lm, glm and many
#'   others.
#' @param x Character vector specifying variables to be added to model or a
#'   formula giving a template which specifies how to update.
#' @param backquote Logical indicating if added variables should be quoted.
#'   Longer than one element vectors are accepted as well, specifying which new
#'   variables should be backquoted. Only relevant if x is of type character.
#' @param collapse Character specifying how new characters should be added to
#'   old formula. Only relevant if x is of type character.
#'
#' @return Updated fit of input model.
#'
#' @importFrom assertthat assert_that is.flag is.string
#' @importFrom stats update
#' @importFrom purrr is_formula
#'
#' @examples
#' library("survival")
#' hla_calls_file <- system.file("extdata", "HLAHD_output_example.txt", package = "MiDAS")
#' hla_calls <- readHlaCalls(hla_calls_file)
#' pheno_file <- system.file("extdata", "pheno_example.txt", package = "MiDAS")
#' pheno <- read.table(pheno_file, header = TRUE)
#' covar_file <- system.file("extdata", "covar_example.txt", package = "MiDAS")
#' covar <- read.table(covar_file, header = TRUE)
#' midas_data <- prepareHlaData(hla_calls, pheno, covar, inheritance_model = "additive")
#' coxmod <- hlaAssocModel(model = "coxph",
#'                         response = "Surv(OS, OS_DIED)",
#'                         variable = "1",
#'                         data = midas_data
#' )
#' updateModel(coxmod, "A*01:01", backquote = TRUE, collapse = " + ")
#'
#' @export
updateModel <- function(object, x, backquote = TRUE, collapse = " + ") {
  assert_that(
    checkStatisticalModel(object),
    see_if(is.character(x) | is_formula(x),
           msg = "x is not a character vector or formula"
    ),
    is.flag(backquote),
    is.string(collapse)
  )

  if (is.character(x)) {
    x[backquote] <- backquote(x[backquote])
    x <- paste0(". ~ . + ", paste(x, collapse = collapse))
  }

  new_object <- update(object = object, x, evaluate = FALSE)
  new_object <- eval.parent(new_object)

  return(new_object)
}

#' Assert statistical model
#'
#' \code{checkStatisticalModel} asserts if object is an existing fit from a
#' model function such as lm, glm and many others.
#'
#' @inheritParams checkHlaCallsFormat
#'
#' @return Logical indicating if \code{data_frame} is an existing fit from a
#' model function such as lm, glm and many others. Otherwise raise error.
#'
#' @importFrom assertthat assert_that see_if
#' @examples
#' object <- lm(dist ~ speed, data = cars)
#' checkStatisticalModel(object)
#'
#' @export
checkStatisticalModel <- function(object) { # TODO simplyfy output of this function; or something like object is not a stat model: potential problem bla bla
  assert_that(
    see_if(is.object(object),
           msg = "object have to have the internal OBJECT bit set"
    ),
    {
      object_call <- getCall(object)
      if (! is.null(object_call)) {
        object_formula <- eval(substitute(formula, env = as.list(object_call)))
        see_if(is_formula(object_formula),
               msg = "object have to be a model with defined formula"
        )
      } else {
        structure(FALSE, msg = "object have to have an attribute 'call'")
      }
    }
  )
}
