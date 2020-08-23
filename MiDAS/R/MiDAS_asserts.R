#' Assert hla calls data frame format
#'
#' \code{checkHlaCallsFormat} asserts if hla calls data frame have proper
#' format.
#'
#' @param hla_calls HLA calls data frame, as returned by
#'   \code{\link{readHlaCalls}} function.
#'
#' @return Logical indicating if \code{hla_calls} follows hla calls data frame
#'   format. Otherwise raise an error.
#'
#' @importFrom assertthat assert_that see_if
#'
checkHlaCallsFormat <- function(hla_calls) {
  an <- deparse(substitute(hla_calls))
  assert_that(
    is.data.frame(hla_calls),
    see_if(nrow(hla_calls) >= 1 & ncol(hla_calls) >= 2,
           msg = paste0(an, " have to have at least 1 rows and 2 columns.")
    ),
    see_if(! any(vapply(hla_calls, is.factor, logical(length = 1))),
           msg = paste0(an, " can't contain factors")
    ),
    see_if(! all(checkAlleleFormat(as.character(hla_calls[, 1])), na.rm = TRUE),
           msg = paste0("first column of ", an, " should specify samples id")
    )
  )

  alleles <- unlist(hla_calls[, -1])
  test_values <- checkAlleleFormat(alleles)
  alleles <- alleles[! test_values & ! is.na(alleles)]
  truncated <- ifelse(length(alleles) < 10, "", ", ...")
  msg <- sprintf("values: %s%s in %s doesn't follow HLA numbers specification", paste(alleles, collapse = ", "), truncated, an)
  assert_that(
    all(test_values, na.rm = TRUE),
    msg = msg
  )

  return(TRUE)
}

#' Assert KIR counts data frame format
#'
#' \code{checkKirCallsFormat} asserts if KIR counts data frame have proper
#' format.
#'
#' @param kir_calls KIR calls data frame, as returned by
#'   \code{\link{readKIRCalls}} function.
#'
#' @return Logical indicating if \code{kir_calls} follow KIR counts data frame
#'   format. Otherwise raise an error.
#'
#' @importFrom assertthat assert_that see_if
#'
checkKirCallsFormat <- function(kir_calls) {
  assert_that(
    is.data.frame(kir_calls),
    see_if(!any(vapply(
      kir_calls, is.factor, logical(length = 1)
    )),
    msg = paste0(
      "kir_calls can't contain factors"
    )),
    colnamesMatches(
      kir_calls,
      c(
        "ID",
        "KIR3DL3",
        "KIR2DS2",
        "KIR2DL2",
        "KIR2DL3",
        "KIR2DP1",
        "KIR2DL1",
        "KIR3DP1",
        "KIR2DL4",
        "KIR3DL1",
        "KIR3DS1",
        "KIR2DL5",
        "KIR2DS3",
        "KIR2DS5",
        "KIR2DS4",
        "KIR2DS1",
        "KIR3DL2"
      )
    )
  )

  return(TRUE)
}

#' Check if frequencies can be calculated for an experiment
#'
#' \code{isExperimentCountsOrZeros} checks if experiment contains only positive
#' integers or zeros.
#'
#' @param x Matrix or SummarizedExperiment object.
#' @param na.rm Logical indicating if \code{NA} values should be accepted.
#'
#' @return Logical indicating if \code{x} contains only positive integers or
#'   zeros.
#'
#' @family assert functions
#'
isExperimentCountsOrZeros <- function(x, na.rm = TRUE) {
  test <- if (is.matrix(x)) {
    isCountsOrZeros(x)
  } else if (isClass(x, "SummarizedExperiment")) {
    isCountsOrZeros(assay(x))
  } else {
    FALSE
  }

  return(test)
}

#' Assert statistical model
#'
#' \code{checkStatisticalModel} asserts if object is an existing fit from a
#' model functions such as lm, glm and many others. Containing MiDAS object as
#' its data atribute.
#'
#' @inheritParams updateModel
#'
#' @return Logical indicating if \code{object} is an existing fit from a
#' model functions such as lm, glm and many others. Containing MiDAS object as
#' its data attribute. Otherwise raise an error.
#'
#' @importFrom assertthat assert_that see_if
#' @importFrom stats getCall
#'
checkStatisticalModel <- function(object) {
  on <- deparse(substitute(object))
  msg <- " was not recognized as a fit from a model function (such as lm, glm and many others)."
  assert_that(
    is.object(object),
    msg = paste0(on, msg)
  )

  object_call <- getCall(object)
  assert_that(
    ! is.null(object_call),
    msg = paste0(on, msg, " ", on, " does not have 'call' attribute.")
  )

  # Preserve object environment; this is important because arguments are given
  # in specific environment, but it is useful to be able to evaluate it in
  # original context -- doh!
  object_env <- attr(object$terms, ".Environment")

  object_formula <- eval(object_call[["formula"]], envir = object_env)
  assert_that(
    is_formula(object_formula),
    msg = paste0(on, msg, " ", on, " does not have 'formula' attribute.")
  )

  object_data <- eval(object_call[["data"]], envir = object_env)
  assert_that(!is.null(object_data),
              msg = paste0(on, msg, " ", on, " does not have 'data' attribute.")
  )

  return(TRUE)
}

#' Check if tidy method for class exist
#'
#' \code{hasTidyMethod} check if there is a tidy method available for a given class.
#'
#' @param class String giving object class.
#'
#' @return Logical indicating if there is a tidy method for a given class.
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
  paste0("Could not find 'tidy' function for statistical model '",
         eval(expr = call$class, envir = env),
         "'. Please ensure that 'tidy' for selected model is available. See the 'broom' package for more information on 'tidy' function."
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
#' \code{isCharacterOrNULL} checks if the object is a character vector or NULL.
#'
#' @param x object to test.
#'
#' @return Logical indicating if object is character vector or NULL
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

#' Check column names
#'
#' \code{colnamesMatches} check if  data frame's columns are named as specified
#'
#' @param x Data frame to test.
#' @param cols Ordered character vector to test against \code{x}'s colnames.
#'
#' @return Logical indicating if \code{x}'s colnames equals \code{choice}.
#'
#' @importFrom assertthat assert_that
#'
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
  mask <- ! curr_colnames %in% future_colnames

  sprintf("Columns: '%s' in %s should be named '%s'",
          paste(curr_colnames[mask], collapse = "', '"),
          deparse(call$x),
          paste(future_colnames[mask], collapse = "', '")
  )
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

#' Assert colData data
#'
#' \code{checkColDataFormat} asserts if the colData data frame has proper format.
#'
#' @param data_frame Data frame containing colData data used to construct
#'   \code{\link{MiDAS}} object.
#'
#' @return Logical indicating if \code{data_frame} is properly formatted.
#'   Otherwise raise an error.
#'
#' @family assert functions
#'
#' @importFrom assertthat assert_that see_if
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
