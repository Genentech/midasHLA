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
      msg = "placeholder argument can be used only with one new variable in x."
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
    object_form <- object_call[["formula"]] %>%
      eval(envir = object_env)
    assert_that(
      objectHasPlaceholder(object, placeholder = placeholder)
    )
    x <- gsub(pattern = placeholder, replacement = x, x = deparse(object_form))
  }

  # print(x)
  new_object <- update(object = object, x, evaluate = FALSE)
  new_object <- eval(new_object, envir = object_env)

  return(new_object)
}

#' Assert statistical model
#'
#' \code{checkStatisticalModel} asserts if object is an existing fit from a
#' model function such as lm, glm and many others. Containing MiDAS object as
#' its data atribute.
#'
#' @inheritParams updateModel
#'
#' @return Logical indicating if \code{object} is an existing fit from a
#' model function such as lm, glm and many others. ontaining MiDAS object as
#' its data atribute. Otherwise raise error.
#'
#' @family assert functions
#'
#' @importFrom assertthat assert_that see_if
#' @importFrom stats getCall
#'
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
  assert_that(!is.null(object_data),
              msg = "object need to have data attribute defined")
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

assertthat::on_failure(stringMatches) <- function(call, env) {
  paste0(deparse(call$x),
         ' should be one of "',
         paste(eval(call$choice, envir = env), collapse = '", "'),
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
  test <- if (length(x)) {
    all(x %in% choice)
  } else {
    FALSE
  }

  return(test)
}

assertthat::on_failure(characterMatches) <- function(call, env) {
  paste0(deparse(call$x),
         ' should match values "',
         paste(eval(call$choice, envir = env), collapse = '", "'),
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

#' Check if placeholder is present in object formula
#'
#' \code{isTRUEorFALSE} check if object is a flag (a length one logical vector)
#' except NA.
#'
#' @param object statistical model to test.
#' @param placeholder string specifying name of placeholder.
#'
#' @return Logical indicating if placeholder is present in object formula.
#'
#' @family assert functions
#'
objectHasPlaceholder <- function(object, placeholder) {
  object_env <- attr(object$terms, ".Environment")
  object_call <- getCall(object)
  object_form <- object_call[["formula"]] %>%
    eval(envir = object_env)
  test <- placeholder %in% all.vars(object_form)

  return(test)
}

assertthat::on_failure(objectHasPlaceholder) <- function(call, env) {
  paste0("placeholder '",
         eval(call$placeholder, envir = env),
         "' could not be found in object's formula"
  )
}

#' Get attributes of statistical model object
#'
#' \code{getObjectDetails} extracts some of the statistical model object
#' attributes that are needed for \code{runMiDAS} internal calculations.
#'
#' TODO write unit tests
#'
#' @inheritParams checkStatisticalModel
#'
#' @return List with following elemnts:
#' \describe{
#'   \item{call}{Object's call}
#'   \item{formula_vars}{Character containing names of variables in object
#'     formula}
#'   \item{data}{MiDAS object associated with model}
#' }
#'
#' @importFrom MultiAssayExperiment colData
#'
getObjectDetails <- function(object) {
  object_call <- getCall(object)
  object_env <- attr(object$terms, ".Environment")
  object_formula <- eval(object_call[["formula"]], envir = object_env)
  object_data <- eval(object_call[["data"]], envir = object_env)

  object_details <- list(
    call = object_call,
    formula_vars = all.vars(object_formula),
    data = object_data
  )

  return(object_details)
}

#' Assert colData data
#'
#' \code{checkColDataFormat} asserts if colData data frame has proper format.
#'
#' @param data_frame Data frame containing colData data used to construct
#'   \code{\link{MiDAS}} object.
#'
#' @return Logical indicating if \code{data_frame} is properly formatted.
#'   Otherwise raise error.
#'
#' @family assert functions
#'
#' @importFrom assertthat assert_that see_if
#' @examples
#' \dontrun{
#' pheno_file <- system.file("extdata", "pheno_example.txt", package = "MiDAS")
#' pheno <- read.table(pheno_file, header = TRUE)
#' checkColDataFormat(pheno)
#' }
#'
checkColDataFormat <- function(data_frame) {
  data_frame_name <- deparse(substitute(data_frame))
  assert_that(
    see_if(
      is.data.frame(data_frame),
      msg = sprintf("%s have to be a data frame",
                    data_frame_name)
    ),
    see_if(
      nrow(data_frame) >= 1 & ncol(data_frame) >= 2,
      msg = sprintf("%s have to have at least 1 row and 2 columns",
                    data_frame_name)
    ),
    see_if(
      colnames(data_frame)[1] == "ID",
      msg = sprintf(
        "first column in %s must be named 'ID'",
        data_frame_name
      )
    )
  )

  return(TRUE)
}

#' Assert kir calls data frame format
#'
#' \code{checkKirCallsFormat} asserts if kir calls data frame have proper
#' format.
#'
#' @param kir_calls Data frame containing KIR calls, as return by
#'   \code{\link{readKirCalls}} function.
#'
#' @return Logical indicating if \code{kir_calls} follows kir calls data frame
#'   format. Otherwise raise error.
#'
#' @family assert functions
#'
#' @importFrom assertthat assert_that see_if
#' @examples
#' file <- system.file("extdata", "KIP_output_example.txt", package = "MiDAS")
#' kir_calls <- readKirCalls(file)
#' checkKirCallsFormat(kir_calls)
#'
#' @export
checkKirCallsFormat <- function(kir_calls) { # TODO
  test <- is.data.frame(kir_calls)
  return(test)
}

#' Check if function can be found in environment
#'
#' \code{functionExists} check if function exists
#'
#' @param name String giving name of function to check.
#'
#' @return Logical indicating if function exists.
#'
#' @family assert functions
#'
functionExists <- function(name) {
  fun <- get0(name)
  test <- is.function(fun)

  return(test)
}

assertthat::on_failure(functionExists) <- function(call, env) {
  sprintf("Function %s could not be found.", call$name)
}


#' Add variables frequencies to runMiDAS results
#'
#' This is a helper function
#'
#' @param midas MiDAS object
#' @param experiment String
#' @param test_covar String giving name of test covariate
#'
#' @importFrom dplyr filter left_join rename select
#' @importFrom magrittr %>%
#' @importFrom rlang !! :=
#'
runMiDASGetVarsFreq <- function(midas, experiment, test_covar) {
  variables_freq <- midas[[experiment]] %>%
    getExperimentFrequencies(inheritance_model = getInheritanceModel(midas)) %>%
    rename(Ntotal = .data$Counts, Ntotal.frequency = .data$Freq)

  test_covar_vals <- factor(colData(midas)[[test_covar]])
  if (nlevels(test_covar_vals) == 2) {
    ids <- rownames(colData(midas))

    lvl1_ids <- ids[test_covar_vals == levels(test_covar_vals)[1]]
    lvl1_freq <- midas[, lvl1_ids] %>%
      `[[`(experiment) %>%
      getExperimentFrequencies(inheritance_model = getInheritanceModel(midas)) %>%
      rename(
        !!sprintf("N(%s=%s)", test_covar, levels(test_covar_vals)[1]) := .data$Counts,
        !!sprintf("N.frequency(%s=%s)", test_covar, levels(test_covar_vals)[1]) := .data$Freq
      )
    variables_freq <- left_join(variables_freq, lvl1_freq, by = "term")

    lvl2_ids <- ids[test_covar_vals == levels(test_covar_vals)[2]]
    lvl2_freq <- midas[, lvl2_ids] %>%
      `[[`(experiment) %>%
      getExperimentFrequencies(inheritance_model = getInheritanceModel(midas)) %>%
      rename(
        !!sprintf("N(%s=%s)", test_covar, levels(test_covar_vals)[2]) := .data$Counts,
        !!sprintf("N.frequency(%s=%s)", test_covar, levels(test_covar_vals)[2]) := .data$Freq
      )
    variables_freq <- left_join(variables_freq, lvl2_freq, by = "term")
  }

  return(variables_freq)
}

#' Check if object is of class x
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
#'
isClass <- function(x, class) {
  test <- is(x, class)

  return(test)
}

assertthat::on_failure(isClass) <- function(call, env) {
  paste0(deparse(call$x),
         " must be an instance of ",
         deparse(call$class),
         "."
  )
}

#' Grantham distance
#'
#' \code{distGrantham} calculates normalized Grantham distance between two
#' amino acid sequences.
#'
#' Distance between amino acid sequences is normalized by length of compared
#' sequences.
#'
#' Lengths of \code{aa1} and \code{aa2} must be equal.
#'
#' @param aa1 Character vector giving amino acid sequence using one letter
#'   codings. Each element must correspond to single amino acid.
#' @param aa2 Character vector giving amino acid sequence using one letter
#'   codings. Each element must correspond to single amino acid.
#'
#' @return Integer normalized Grantham distance between \code{aa1} and \
#'   code{aa2}.
#'
#' @examples
#' distGrantham(
#'   aa1 = c("A", "S", "W"),
#'   aa2 = c("A", "S", "V")
#' )
#'
#' @importFrom assertthat assert_that see_if
#'
#' @export
distGrantham <- function(aa1, aa2) {
  assert_that(
    is.character(aa1),
    is.character(aa2),
    see_if(
      length(aa1) == length(aa2),
      msg = "aa1 and aa2 must have equal lengths."
    )
  )

  idx <- paste(aa1, aa2, sep = "")
  assert_that(
    all(test <- idx %in% names(dict_dist_grantham)),
    msg = sprintf(
      fmt = "%s are not valid amino acids pairs",
      paste(idx[! test], collapse = ", ")
    )
  )

  d <- sum(dict_dist_grantham[idx]) / length(idx)

  return(d)
}

#' Calculate Grantham disttance between alleles
#'
#' \code{hlaCallsGranthamDistance} calculate distance between HLA alleles using
#' distance of choice.
#'
#' The Grantham distance is calculated for pairs of alleles for HLA gene.
#'
#' In case of missing alleles \code{NA} are returned.
#'
#' Grantham distance is calculated only for class I HLA alleles. First
#' exons forming the variable region in the peptide binding groove (i.e.,
#' exon 2 and 3) are selected (positions 1-182 in IMGT/HLA alignments, however
#' here we take 2-182 as many 1st positions are missing). Then all the alleles
#' containing gaps, stop codons or indels are discarded. Finally distance is
#' calculated for each pair.
#'
#' @inheritParams checkHlaCallsFormat
#' @param genes Character vector spcifying genes for which allelic distance
#'   should be calculated.
#'
#' @return Data frame of normalized Grantham distances between pairs of alleles
#'   for each specified HLA gene. First column (\code{ID}) is the same as in
#'   \code{hla_calls}, further columns are named as given by \code{genes} and
#'   contain numeric values.
#'
#' @importFrom assertthat assert_that is.string noNA see_if
#' @importFrom magrittr %>%
#' @importFrom rlang warn
#' @importFrom stats na.omit
#'
#' @examples
#' file <- system.file("extdata", "HLAHD_output_example.txt", package = "MiDAS")
#' hla_calls <- readHlaCalls(file)
#' hlaCallsGranthamDistance(hla_calls, genes = "A")
#'
#' @export
hlaCallsGranthamDistance <- function(hla_calls, genes = c("A", "B", "C")) {
  assert_that(
    checkHlaCallsFormat(hla_calls),
    is.character(genes),
    noNA(genes)
  )

  target_genes <- getHlaCallsGenes(hla_calls)
  assert_that(
    characterMatches(x = genes, choice = target_genes) # asset is buggy hlaAlleleDistance(hla_calls, genes = c("A", "B", "FA"))
  )
  target_genes <- genes

  d <- list(ID = hla_calls[, 1, drop = TRUE])
  for (gene in target_genes) {
    if (! gene %in% c("A", "B", "C")) {
      warn(sprintf("Grantham distance is calculated only for class I HLA alleles. Ommiting gene %s", gene))
      next()
    }

    sel <- paste0(gene, c("_1", "_2"))
    pairs <- hla_calls[, sel]
    resolution <- getAlleleResolution(unlist(pairs)) %>%
      na.omit()
    assert_that( # this should become obsolete
      see_if(
        all(resolution == resolution[1]),
        msg = sprintf("Allele resolutions for gene %s are not equal", gene)
      )
    )

    # process alignment
    alignment <-hlaAlignmentGrantham(gene, resolution[1])

    allele_numbers <- rownames(alignment)
    d[[gene]] <- vapply(
      X = 1:nrow(pairs),
      FUN = function(i) {
        allele1 <- pairs[i, 1]
        allele2 <- pairs[i, 2]
        if (allele1 %in% allele_numbers && allele2 %in% allele_numbers) {
          aa1 <- alignment[allele1, ]
          aa2 <- alignment[allele2, ]
          distGrantham(aa1, aa2)
        } else {
          warn(
            sprintf(
              fmt = "Alleles %s could not be found in the alignmnet coercing to NA",
              paste(allele1, allele1, sep = ", ")
            )
          )
          as.numeric(NA)
        }
      },
      FUN.VALUE = numeric(length = 1L)
    )
  }

  hla_dist <- data.frame(d, stringsAsFactors = FALSE)

  return(hla_dist)
}

#' Helper function alignment for Grantham distance calculations
#'
#' \code{hlaAlignmentGrantham} get HLA alignment processed so that grantham
#' distance between alleles can be calculated. Processing includes extracting
#' exons 1 and 2, masking indels, gaps and stop codons.
#'
#' @param gene Character
#' @param resolution Number
#'
hlaAlignmentGrantham <- function(gene, resolution) {
  alignment <- readHlaAlignments(
    gene = gene,
    resolution = resolution,
    trim = TRUE
  )
  alignment <- alignment[, 2:182] # select exons 2 and 3 w/o 1st position as it is biased towards missing data
  mask <- apply(alignment, 1, function(x) any(x == "" | x == "X" | x == ".")) # mask gaps, stop codons, indels
  alignment <- alignment[! mask, ]
}

#' Validate frequency cutoffs
#'
#' \code{validateFrequencyCutoffs} checks if \code{lower_frequency_cutoff} and
#' \code{upper_frequency_cutoff} are valid.
#'
#' \code{lower_frequency_cutoff} and \code{upper_frequency_cutoff} should be a
#' positive numbers, giving either frequency or counts.
#' \code{lower_frequency_cutoff} has to be lower than
#' \code{upper_frequency_cutoff}.
#'
#' @param lower_frequency_cutoff Number
#' @param upper_frequency_cutoff Number
#'
#' @return Logical indicating if \code{lower_frequency_cutoff} and
#'   \code{upper_frequency_cutoff} are valid.
#'
#' @family assert functions
#'
#' @importFrom assertthat assert_that see_if
#'
validateFrequencyCutoffs <- function(lower_frequency_cutoff, upper_frequency_cutoff) {
  assert_that(
    isNumberOrNULL(lower_frequency_cutoff),
    if (! is.null(lower_frequency_cutoff)) {
      see_if(lower_frequency_cutoff >= 0,
             msg = "lower_frequency_cutoff must be a number greater than 0."
      )
    } else TRUE,
    isNumberOrNULL(upper_frequency_cutoff),
    if (! is.null(upper_frequency_cutoff)) {
      see_if(upper_frequency_cutoff >= 0,
             msg = "upper_frequency_cutoff must be a number greater than 0."
      )
    } else TRUE,
    if (! is.null(lower_frequency_cutoff) && ! is.null(upper_frequency_cutoff)) {
      see_if(! lower_frequency_cutoff > upper_frequency_cutoff,
             msg = "lower_frequency_cutoff cannot be higher than upper_frequency_cutoff."
      )
    } else TRUE,
    if (! is.null(lower_frequency_cutoff) && ! is.null(upper_frequency_cutoff)) {
      see_if((lower_frequency_cutoff <= 1 && upper_frequency_cutoff <= 1) ||
               (lower_frequency_cutoff >= 1 && upper_frequency_cutoff >= 1),
             msg = "Both lower_frequency_cutoff and upper_frequency_cutoff have to be either frequencies or counts."
      )
    } else TRUE
  )
}

#' Get HLA calls genes
#'
#' \code{getHlaCallsGenes} get's genes found in HLA calls.
#'
#' @inheritParams checkHlaCallsFormat
#'
#' @return Character vector of genes in \code{hla_calls}.
#'
#' @importFrom assertthat assert_that
#'
getHlaCallsGenes <- function(hla_calls) {
  assert_that(
    checkHlaCallsFormat(hla_calls)
  )

  genes <- hla_calls %>%
    colnames() %>%
    `[`(-1) %>% # discard ID column
    gsub(pattern = "_[0-9]+", replacement = "") %>%
    unique()

  return(genes)
}

#' Helper transform data frame to experiment matrix
#'
#' Function deletes 'ID' column from a \code{df}, then transpose it and sets
#' the colum names to values from deleted 'ID' column.
#'
#' @param df Data frame
#'
#' @return Matrix
#'
#' @importFrom utils type.convert
#'
dfToExperimentMat <- function(df) {
  cols <- df[["ID"]]
  mat <- t(subset(df, select = -ID))
  colnames(mat) <- cols

  # convert to apropiate type
  mat <- type.convert(mat, as.is = TRUE)

  return(mat)
}

#' Helper transform experiment matrix to data frame
#'
#' Function transpose \code{mat} and inserts colum names of input \code{mat} as
#' a 'ID' column.
#'
#' @param mat Matrix
#'
#' @return Data frame
#'
experimentMatToDf <- function(mat) {
  ID <- colnames(mat)
  df <-
    as.data.frame(
      t(mat),
      stringsAsFactors = FALSE,
      make.names = FALSE
    )
  rownames(df) <- NULL
  df <- cbind(ID, df, stringsAsFactors = FALSE)

  return(df)
}

#' Transform MiDAS to wide format data.frame
#'
#' @param object Object of class MiDAS
#' @param experiment Character
#'
#' @return Data frame
#'
#' @importFrom assertthat assert_that
#' @importFrom methods validObject
#' @importFrom MultiAssayExperiment longFormat colData
#' @importFrom tidyr spread
#'
midasToWide <- function(object, experiment) {
  assert_that(
    validObject(object),
    is.character(experiment),
    characterMatches(experiment, getExperiments(object))
  )

  object <- object[, , experiment]
  wide_df <-
    object %>%
    longFormat(colDataCols = TRUE) %>%
    as.data.frame() %>%
    subset(select = -assay) %>%
    subset(select = -colname) %>%
    spread(key = "rowname", value = "value")

  return(wide_df)
}
