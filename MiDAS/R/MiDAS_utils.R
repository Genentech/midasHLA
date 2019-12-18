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
  pattern <- "^[A-Z0-9]+[*][0-9]+(:[0-9]+){0,3}((?=G)(G|GG)|(N|L|S|C|A|Q)){0,1}$"
  is_correct <- stri_detect_regex(allele, pattern)
  return(is_correct)
}

#' Infers HLA allele resolution
#'
#' \code{getAlleleResolution} returns the resolution of input HLA allele
#' numbers.
#'
#' HLA allele resolution can take the following values: 2, 4, 6, 8.See
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
#' \code{getVariableAAPos} finds variable amino acid positions in the
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
#' @param alignment Matrix containing amino acid level alignment.
#' @param varchar Regex matching characters that should be considered when
#'   looking for variable amino acid positions. See details for further
#'   explanations.
#'
#' @return Integer vector specifying which alignment columns are variable.
#'
#' @examples
#' alignment <- readHlaAlignments(gene = "TAP1")
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
#' instead.
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
#' dictionary <- system.file("extdata", "Match_allele_HLA_supertype.txt", package = "MiDAS")
#' convertAlleleToVariable(c("A*01:01", "A*02:01"), dictionary = dictionary)
#'
#' @importFrom assertthat assert_that is.string is.readable see_if
#' @importFrom stats setNames
#' @importFrom utils read.table
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
      sep = "\t",
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
  variable <- dictionary[as.character(allele)] # for all NAs vectrors default type is logical; recycling mechanism leads to unwanted results
  names(variable) <- NULL

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
#' @family assert functions
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
           msg = "hla_calls have to have at least 1 rows and 2 columns. Make sure the input file is in tsv format."
    ),
    see_if(! any(vapply(hla_calls, is.factor, logical(length = 1))),
           msg = "hla_calls can't contain factors"
    ),
    see_if(! all(checkAlleleFormat(as.character(hla_calls[, 1])), na.rm = TRUE),
           msg = "first column of hla_calls should specify samples id"
    )
  )

  alleles <- unlist(hla_calls[, -1])
  test_values <- checkAlleleFormat(alleles)
  alleles <- alleles[! test_values & ! is.na(alleles)]
  assert_that(
      all(test_values, na.rm = TRUE),
      msg = sprintf(
        "values: %s in hla_calls doesn't follow HLA numbers specification",
        paste(alleles, collapse = ", ")
      )
  )

  return(TRUE)
}

#' Backquote string
#'
#' \code{backquote} places backticks around string.
#'
#' \code{backquote} is useful when using HLA allele numbers in formulas, where
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
  backquoted <- paste0("`", x, "`")
  return(backquoted)
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
#' @family assert functions
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
      see_if(any(hla_calls[, 1, drop = TRUE] %in% data_frame[, 1, drop = TRUE]),
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
#' \code{updateModel} adds new variables to model and re-fit it.
#'
#' @param object An existing fit from a model function such as lm, glm and many
#'   others.
#' @param x Character vector specifying variables to be added to model.
#' @param placeholder String specifying term to substitute with value from
#'   \code{x}. Can be used only if \code{x} is a string. Ignored if set to
#'   \code{NULL}.
#' @param backquote Logical indicating if added variables should be quoted.
#'   Elements of this vector are recycled over \code{x}. Only relevant if
#'   \code{x} is of type character.
#' @param collapse String specifying how new characters should be added to
#'   old formula. Only relevant if \code{placeholder} is specified.
#'
#' @return Updated fit of input model.
#'
#' @importFrom assertthat assert_that is.string
#' @importFrom magrittr %>%
#' @importFrom stats update
#' @importFrom purrr is_formula
#'
#' @examples
#' object <- lm(dist ~ 1, data = cars)
#' updateModel(object, "dist")
#'
#' @export
updateModel <- function(object,
                        x,
                        placeholder = NULL,
                        backquote = TRUE,
                        collapse = " + ") {
  assert_that(
    checkStatisticalModel(object),
    is.character(x),
    isStringOrNULL(placeholder),
    see_if(
      ! (! is.null(placeholder) && length(x) != 1),
      msg = "placeholder argument can be used only with one new variable x."
    ),
    isTRUEorFALSE(backquote),
    is.string(collapse)
  )

  object_env <- attr(object$terms, ".Environment")

  if (backquote) {
    x <- backquote(x)
  }

  if (is.null(placeholder)) {
    x <- paste0(". ~ . + ", paste(x, collapse = collapse))
  } else {
    object_call <- getCall(object)
    x <- object_call[["formula"]] %>%
      eval(envir = object_env) %>%
      deparse() %>%
      gsub(pattern = placeholder, replacement = x)
  }

  # print(x)
  new_object <- update(object = object, x, evaluate = FALSE)
  new_object <- eval(new_object, envir = object_env)

  return(new_object)
}

#' Assert statistical model
#'
#' \code{checkStatisticalModel} asserts if object is an existing fit from a
#' model function such as lm, glm and many others.
#'
#' @inheritParams updateModel
#'
#' @return Logical indicating if \code{object} is an existing fit from a
#' model function such as lm, glm and many others. Otherwise raise error.
#'
#' @family assert functions
#'
#' @importFrom assertthat assert_that see_if
#' @importFrom stats getCall
#' @examples
#' object <- lm(dist ~ speed, data = cars)
#' checkStatisticalModel(object)
#'
#' @export
checkStatisticalModel <- function(object) { # TODO simplyfy output of this function; or something like object is not a stat model: potential problem bla bla
  assert_that(
    is.object(object),
    msg = "object is required to have the internal OBJECT bit set"
  )

  object_call <- getCall(object)
  assert_that(
    ! is.null(object_call),
    msg = "object have to have an attribute 'call'"
  )

  object_env <- attr(object$terms, ".Environment")

  object_formula <- eval(object_call[["formula"]], envir = object_env)
  assert_that(
    is_formula(object_formula),
    msg = "object have to be a model with defined formula"
  )

  object_data <- eval(object_call[["data"]], envir = object_env)
  assert_that(
    ! is.null(object_data) & is.data.frame(object_data),
    msg = "object need to have data attribute defined"
  )
}

#' Check if vector contains only counts or zeros
#'
#' \code{isCountsOrZeros} checks if vector contains only positive integers or
#' zeros.
#'
#' @param x Numeric vector or object that can be \code{unlist} to numeric
#'   vector.
#' @param na.rm Logical indicating if \code{NA} values should be accepted.
#'
#' @return Logical indicating if provided vector contains only positive integers
#'   or zeros.
#'
#' @family assert functions
#'
#' @importFrom rlang is_integerish
#'
isCountsOrZeros <- function(x, na.rm = TRUE) {
    x <- unlist(x)
    test <- is_integerish(x) & x >= 0
    test <- all(test, na.rm = na.rm)

  return(test)
}

#' Error message for isCountsOrZeros
#'
#' @inheritParams assertthat::on_failure
#'
assertthat::on_failure(isCountsOrZeros) <- function(call, env) {
  paste0("values in ", deparse(call$x), " are not counts (a positive integers) or zeros.")
}

#' Check if object is character vector or NULL
#'
#' \code{isCharacterOrNULL} checks if object is character vector or NULL.
#'
#' @param x object to test.
#'
#' @return Logical indicating if object is character vector or NULL
#'
#' @family assert functions
#'
isCharacterOrNULL <- function(x) {
    test <- is.character(x) | is.null(x)

  return(test)
}

#' Error message for isCharacterOrNULL
#'
#' @inheritParams assertthat::on_failure
#'
assertthat::on_failure(isCharacterOrNULL) <- function(call, env) {
  paste0(deparse(call$x), " is not a character vector or NULL.")
}

#' Check if object is number or NULL
#'
#' \code{isNumberOrNULL} checks if object is number (a length one numeric
#' vector) or NULL.
#'
#' @param x object to test.
#'
#' @return Logical indicating if object is number or NULL
#'
#' @family assert functions
#'
#' @importFrom assertthat is.number
#'
isNumberOrNULL <- function(x) {
    test <- is.number(x) | is.null(x)

  return(test)
}

#' Error message for isNumberOrNULL
#'
#' @inheritParams assertthat::on_failure
#'
assertthat::on_failure(isNumberOrNULL) <- function(call, env) {
  paste0(deparse(call$x),
         " is not a number (a length one numeric vector) or NULL."
  )
}

#' Check if object is string or NULL
#'
#' \code{isStringOrNULL} checks if object is string (a length one character
#' vector) or NULL.
#'
#' @param x object to test.
#'
#' @return Logical indicating if object is string or NULL
#'
#' @family assert functions
#'
#' @importFrom assertthat is.string
#'
isStringOrNULL <- function(x) {
    test <- is.string(x) | is.null(x)

  return(test)
}

#' Error message for isStringOrNULL
#'
#' @inheritParams assertthat::on_failure
#'
assertthat::on_failure(isStringOrNULL) <- function(call, env) {
  paste0(deparse(call$x),
         " is not a string (a length one character vector) or NULL."
  )
}

#' Check if string matches one of possible values
#'
#' \code{stringMatches} checks if string is equal to one of the choices.
#'
#' @param x string to test.
#' @param choice Character vector with possible values for \code{x}.
#'
#' @return Logical indicating if \code{x} matches one of the strings in
#'   \code{choice}.
#'
#' @family assert functions
#'
stringMatches <- function(x, choice) {
    test <- x %in% choice

  return(test)
}

#' Error message for stringMatches
#'
#' @inheritParams assertthat::on_failure
#'
assertthat::on_failure(stringMatches) <- function(call, env) {
  paste0(deparse(call$x),
         ' should be one of "',
         paste(eval(call$choice), collapse = '", "'),
         '".'
  )
}

#' Check if object is flag or NULL
#'
#' \code{isFlagOrNULL} checks if object is flag (a length one logical vector) or
#' NULL.
#'
#' @param x object to test.
#'
#' @return Logical indicating if object is flag or NULL
#'
#' @family assert functions
#'
#' @importFrom assertthat is.flag
#'
isFlagOrNULL <- function(x) {
    test <- (is.flag(x) && ! is.na(x)) || is.null(x)

  return(test)
}

#' Error message for isFlagOrNULL
#'
#' @inheritParams assertthat::on_failure
#'
assertthat::on_failure(isFlagOrNULL) <- function(call, env) {
  paste0(deparse(call$x),
         " is not a flag (a length one logical vector) or NULL."
  )
}

#' List HLA alleles dictionaries
#'
#' \code{listMiDASDictionaries} lists dictionaries shipped with MiDAS package.
#'
#' @param file.names Logical value. If FALSE, only the names of dictionaries are
#' returned. If TRUE their paths are returned.
#' @param pattern String used to match dictionary names, it can be a regular
#'   expression. By default all names are matched.
#'
#' @return Character vector with names of HLA alleles dictionaries.
#'
#' @export
listMiDASDictionaries <- function(pattern = ".*",
                                  file.names = FALSE) {
  pattern <- paste0("^Match.*", pattern, ".*.txt$")
  lib <- list.files(
    path = system.file("extdata", package = "MiDAS"),
    pattern = pattern,
    full.names = file.names
  )

  if (! file.names) {
    lib <- gsub("^Match_", "", gsub(".txt$", "", lib))
  }

  return(lib)
}

#' Check if character matches one of possible values
#'
#' \code{characterMatches} checks if all elements of a character vector matches
#' values in choices.
#'
#' @param x Character vector to test.
#' @param choice Character vector with possible values for \code{x}.
#'
#' @return Logical indicating if \code{x}'s elements matches any of the values
#' in \code{choice}.
#'
#' @family assert functions
#'
#' @importFrom assertthat assert_that
characterMatches <- function(x, choice) {
  assert_that(is.character(x))
  test <- x %in% choice
  test <- all(test)

  return(test)
}

#' Error message for characterMatches
#'
#' @inheritParams assertthat::on_failure
#'
assertthat::on_failure(characterMatches) <- function(call, env) {
  paste0(deparse(call$x),
         ' should match values "',
         paste(eval(call$choice), collapse = '", "'),
         '".'
  )
}

#' Check if object is of class x or null
#'
#' \code{isClassOrNULL} checks if object is an instance of a specified class or
#' is null.
#'
#' @param x object to test.
#' @param class String specifying class to test.
#'
#' @return Logical indicating if \code{x} is an instance of \code{class}.
#'
#' @family assert functions
#'
#' @importFrom assertthat assert_that
#' @importFrom methods is
isClassOrNULL <- function(x, class) {
  test <- is(x, class) | is.null(x)

  return(test)
}

#' Error message for isClassOrNULL
#'
#' @inheritParams assertthat::on_failure
#'
assertthat::on_failure(isClassOrNULL) <- function(call, env) {
  paste0(deparse(call$x),
         " must be an instance of ",
         deparse(call$class),
         " or NULL."
  )
}

#' Convert KIR haplotypes to gene counts
#'
#' \code{kirHaplotypeToCounts} convert vector of KIR haplotypes to data frame
#' of KIR gene counts.
#'
#' @param x Character vector specifying KIR haplotypes.
#' @param hap_dict String specifying path to KIR haplotypes dictionary. By
#'   default file shipped together with package is being used. See details for
#'   more information.
#' @param binary Logical flag indicating if haplotypes should be converted only
#'   to gene presence / absence indicators (it is the only way that allows
#'   unambiguous conversion). This argument is currently ignored.
#'
#' \code{hap_dict} have to be a tab separated values formatted file with first
#' column holding KIR haplotypes and gene counts in others. File should have
#' header with first column unnamed and gene names in the others.
#'
#' @return Data frame with haplotypes and corresponding gene counts. \code{NA}'s
#'   in \code{x} are removed during conversion.
#'
#' @seealso \code{\link{readKirCalls}}, \code{\link{getHlaKirInteractions}},
#'   \code{\link{checkKirCountsFormat}}, \code{\link{prepareMiDAS}}.
#'
#' @examples
#' x <- c(NA, "1+3|16+3", "1+1", NA)
#' kirHaplotypeToCounts(x)
#'
#' @importFrom assertthat assert_that is.readable see_if
#' @importFrom stats na.omit
#' @importFrom stringi stri_split_fixed
#'
#' @export
kirHaplotypeToCounts <- function(x,
                                 hap_dict = system.file("extdata", "Match_kir_haplotype_gene.txt", package = "MiDAS"),
                                 binary = TRUE) {
  assert_that(
    is.character(x),
    is.readable(hap_dict),
    isTRUEorFALSE(binary)
  )
  binary <- TRUE # its not possible yet to not oparate on binary states - once available delete this line
  hap_dict <- read.table(hap_dict, stringsAsFactors = FALSE)

  # x <- na.omit(x)
  x_split <- stri_split_fixed(x, "|")
  x_split_unlist <- unlist(x_split)

  x_split_unlist_haps <- stri_split_fixed(x_split_unlist, pattern = "+")
  assert_that(
    all(haps_match <- vapply(
      X = x_split_unlist_haps,
      FUN = function(x) all(na.omit(x) %in% rownames(hap_dict)),
      FUN.VALUE = logical(1)
    )),
    msg = sprintf(
      fmt = "%s haplotype was not found in hap_dict",
      paste(x_split_unlist[! haps_match], collapse = ", ")
    )
  )

  counts <- do.call(cbind, x_split_unlist_haps)
  counts <- apply(counts, 2, function(i) colSums(hap_dict[i, ]))
  colnames(counts) <- x_split_unlist

  if (binary) {
    counts <- ifelse(counts > 1, 1, counts)
  }

  counts <- vapply(
    X = x_split,
    FUN = function(hap) {
      hap <- ifelse(is.na(hap), NA, hap) # class character NA cannot be used for columns indexing
      hap <- counts[, hap, drop = FALSE]
      if (ncol(hap) > 1) {
        assert_that(
          all(apply(hap, 2, function(col) all(col == hap[, 1])), na.rm = TRUE),
          msg = sprintf("haplotype %s can not be unambigously converted to counts", hap)
        )
      }
      hap <- hap[, 1, drop = TRUE]
      return(hap)
    },
    FUN.VALUE = numeric(length = ncol(hap_dict))
  )
  counts <- t(counts)
  counts <- as.data.frame(counts, optional = TRUE, stringsAsFactors = FALSE)
  counts <- cbind(haplotypes = x, counts, stringsAsFactors = FALSE)

  return(counts)
}

#' Check column names
#'
#' \code{colnamesMatches} check if  data frame's columns are named as specified
#'
#' @param x Data frame to test.
#' @param cols Ordered character vector to test against \code{x}'s colnames.
#'
#' @return Logical indicating if \code{x}'s colnames equals \code{choice}.
#'
#' @family assert functions
#'
#' @importFrom assertthat assert_that
colnamesMatches <- function(x, cols) {
  assert_that(
    is.data.frame(x),
    see_if(ncol(x) == length(cols),
           msg = sprintf(
             "Number of columns in %s must equal %i.",
             deparse(substitute(x)),
             length(cols)
           )
    )
  )

  columns_names <- colnames(x)
  columns_test <- columns_names == cols
  test <- all(columns_test)

  return(test)
}

#' Error message for colnamesMatches
#'
#' @inheritParams assertthat::on_failure
#'
assertthat::on_failure(colnamesMatches) <- function(call, env) {
  curr_colnames <- colnames(eval(call$x, envir = env))
  future_colnames <- eval(call$cols, envir = env)

  sprintf("Columns %s in %s should be named %s",
           paste(curr_colnames, collapse = ", "),
           deparse(call$x),
           paste(future_colnames, collapse = ", ")
  )
}

#' Assert KIR counts data frame format
#'
#' \code{checkKirCountsFormat} asserts if KIR counts data frame have proper
#' format.
#'
#' @param kir_counts Data frame containing KIR gene counts, as returned by
#'   \code{\link{readKirCalls}} function.
#' @param accept.null Logical indicating if NULL \code{kir_counts} should be
#'   accepted.
#'
#' @return Logical indicating if \code{kir_counts} follow KIR counts data frame
#'   format. Otherwise raise error.
#'
#' @family assert functions
#'
#' @seealso \code{\link{readKirCalls}}, \code{\link{getHlaKirInteractions}},
#'   \code{\link{kirHaplotypeToCounts}}, \code{\link{prepareMiDAS}}.
#'
#' @importFrom assertthat assert_that see_if
#' @examples
#' file <- system.file("extdata", "KIP_output_example.txt", package = "MiDAS")
#' kir_counts <- readKirCalls(file)
#' checkKirCountsFormat(kir_counts)
#'
#' @export
checkKirCountsFormat <- function(kir_counts,
                                 accept.null = FALSE) {
  if (! (is.null(kir_counts) & accept.null)) {
    kir_counts_name <- deparse(substitute(kir_counts))
    assert_that(
      is.data.frame(kir_counts),
      see_if(nrow(kir_counts) >= 1 & ncol(kir_counts) >= 2,
              msg = paste0(kir_counts_name,
                           " have to have at least 1 rows and 2 columns"
              )
      ),
      see_if(! any(vapply(kir_counts, is.factor, logical(length = 1))),
           msg = paste0(kir_counts_name, " can't contain factors")
      )
    )

    kir_counts <- kir_counts[, 1, drop = FALSE]
      assert_that(
        colnamesMatches(kir_counts, "ID")
      )
  }

  return(TRUE)
}

#' Check if object is count or NULL
#'
#' \code{isCountOrNULL} check if object is a count (a single positive integer)
#' or NULL.
#'
#' @param x object to test.
#'
#' @return Logical indicating if object is count or NULL
#'
#' @family assert functions
#'
#' @importFrom assertthat is.count
#'
isCountOrNULL <- function(x) {
  test <- is.count(x) | is.null(x)

  return(test)
}

#' Error message for isCountOrNULL
#'
#' @inheritParams assertthat::on_failure
#'
assertthat::on_failure(isCountOrNULL) <- function(call, env) {
  paste0(deparse(call$x),
         " is not a count (a single positive integer) or NULL."
  )
}

#' Check if object is TRUE or FALSE flag
#'
#' \code{isTRUEorFALSE} check if object is a flag (a length one logical vector)
#' except NA.
#'
#' @param x object to test.
#'
#' @return Logical indicating if object is TRUE or FALSE flag
#'
#' @family assert functions
#'
#' @importFrom assertthat is.flag
#'
isTRUEorFALSE <- function(x) {
  test <- is.flag(x) && ! is.na(x)

  return(test)
}

#' Error message for isTRUEorFALSE
#'
#' @inheritParams assertthat::on_failure
#'
assertthat::on_failure(isTRUEorFALSE) <- function(call, env) {
  paste0(deparse(call$x),
         " is not a flag (a length one logical vector)."
  )
}

#' Check if tidy method for class exist
#'
#' \code{hasTidyMethod} check if there is tidy method available for given class.
#'
#' @param class Object class.
#'
#' @return Logical indicating if there is tidy method for given class.
#'
#' @family assert functions
#'
#' @importFrom utils methods
#'
hasTidyMethod <- function(class) {
  cl_fun <- paste0("tidy.", class)
  test <- cl_fun %in% methods("tidy")

  return(test)
}

#' Error message for hasTidyMethod
#'
#' @inheritParams assertthat::on_failure
#'
assertthat::on_failure(hasTidyMethod) <- function(call, env) {
  paste0("tidy function for object of class ",
         deparse(call$class),
         " could not be found."
  )
}

#' Likelihood ratio test
#'
#' \code{LRTest} carry out asymptotic likelihood ratio test for two models.
#'
#' \code{mod0} have to be a reduced version of \code{mod1}. See examples.
#'
#' @param mod0 An existing fit from a model function such as lm, glm and many
#'   others.
#' @param mod1 Object of the same class as \code{mod0} with extra terms
#'   included.
#'
#' @return Data frame with the results of likelihood ratio test of the supplied
#'   models.
#'
#'   Column \code{term} holds new variables apearing in \code{mod1},
#'   \code{dof} difference in degrees of freedom between models, \code{logLik}
#'   difference in log likelihoods, \code{statistic} \code{Chisq} statistic and
#'   \code{p.value} corresponding p-value.
#'
#' @examples
#' df <- data.frame(OS = c(20, 30, 40), AGE = c(50, 60, 70))
#' mod0 <- lm(OS ~ 1, data = df)
#' mod1 <- lm(OS ~ AGE, data = df)
#' MiDAS:::LRTest(mod0, mod1)
#'
#' @importFrom assertthat assert_that
#' @importFrom stats logLik pchisq
#'
LRTest <- function(mod0, mod1) {
  formula0 <- formula(mod0)
  vars0 <- all.vars(formula0)
  formula1 <- formula(mod1)
  vars1 <- all.vars(formula1)
  assert_that(
    all(vars0 %in% vars1),
    msg = sprintf("variables %s were not found in mod1",
                  paste(vars0[! vars0 %in% vars1], collapse = ", ")
    )
  )

  ll0 <- logLik(mod0)
  ll1 <- logLik(mod1)
  dof <- attr(ll1, "df") - attr(ll0, "df")
  statistic <- 2 * (as.numeric(ll1) - as.numeric(ll0))
  p.value <- pchisq(statistic, df = dof, lower.tail=FALSE)

  new_vars <- paste(vars1[! vars1 %in% vars0], collapse = ", ")
  res <- data.frame(
    term = new_vars,
    dof = dof,
    logLik = ll1 - ll0,
    statistic = statistic,
    p.value = p.value,
    stringsAsFactors = FALSE
  )

  return(res)
}
