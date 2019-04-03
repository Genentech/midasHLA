#' Association analysis of HLA allele calls
#'
#' \code{analyzeHlaAssociations} performs associations analysis on single HLA
#' allele level using statistical model of choice.
#'
#' @inheritParams checkHlaCallsFormat
#' @inheritParams hlaCallsToCounts
#' @param model String specifying statistical model to use.
#' @param hla_data hla_data object as returned by \code{\link{prepareHlaData}}
#'   function.
#' @param response Character specifying which variables should be treated as
#'   response variable.
#' @param covariate Character specifying which variables should be treated as
#'   covariates. Can be \code{NULL}, than no covariates are considered in the
#'   model.
#' @param correction String specifying multiple testing correction method. See
#'   details for further information.
#' @param exponentiate Logical indicating whether or not to exponentiate the
#'   coefficient estimates. This is typical for logistic and multinomial
#'   regressions, but a bad idea if there is no log or logit link. Defaults to
#'   FALSE.
#' @param ... Further arguments passed to \code{model}.
#'
#' \code{correction} specifies p-value adjustment method to use, common choice
#' is Benjamini & Hochberg (1995) (\code{"BH"}). Internally this is passed to
#' \link{p.adjust}. Check there to get more details.
#'
#' @return Tibble containing combined results for all alleles in
#'   \code{hla_calls}.
#'
#' @examples
#' library("survival")
#'
#' hla_calls_file <- system.file("extdata", "HLAHD_output_example.txt", package = "MiDAS")
#' hla_calls <- readHlaCalls(hla_calls_file)
#' pheno_file <- system.file("extdata", "pheno_example.txt", package = "MiDAS")
#' pheno <- read.table(pheno_file, header = TRUE)
#' covar_file <- system.file("extdata", "covar_example.txt", package = "MiDAS")
#' covar <- read.table(covar_file, header = TRUE)
#' hla_data <- prepareHlaData(hla_calls = hla_calls,
#'                            pheno = pheno,
#'                            covar = covar,
#'                            inheritance_model = "additive"
#' )
#'
#' # Cox proportional hazards regression model
#' analyzeHlaAssociations(model = "coxph",
#'                        hla_data = hla_data,
#'                        response = c("OS", "OS_DIED"),
#'                        covariate = c("AGE", "SEX"),
#'                        correction = "BH"
#' )
#'
#' @importFrom assertthat assert_that see_if is.flag is.string
#' @importFrom broom tidy
#' @importFrom dplyr left_join filter mutate rename
#' @importFrom stats p.adjust
#' @importFrom purrr map_dfr
#'
#' @export
analyzeHlaAssociations <- function(model = "coxph",
                                   hla_data,
                                   response,
                                   covariate,
                                   correction = "BH",
                                   exponentiate = FALSE,
                                   ...) {
  assert_that(
    see_if(is.string(model) | is.function(model),
           msg = "model have to be a string (a length one character vector) or a function"
    ),
    if (is.string(model)) {
      see_if(is.function(get0(model)),
             msg = sprintf("could not find function %s", model)
      )
    } else {
      TRUE
    },
    is.character(response),
    {
      response_len <- length(response)
      is_cox <- as.character(substitute(model)) %in% c("coxph", "cph")
      if (is_cox & response_len != 2) {
        structure(FALSE, msg = "cox survival analysis requires response to be a character vector of length 2")
      } else if (! is_cox & response_len != 1) {
        structure(FALSE, msg = "response is not a string (a length one character vector).")
      } else {
        TRUE
      }
    },
    see_if(all(response %in% colnames(hla_data)),
           msg = "response variables can not be found in hla_data"
    ),
    see_if(is.character(covariate) | is.null(covariate),
           msg = "covariate have to be a character or NULL"
    ),
    {
      if (is.character(covariate)) {
        see_if(all(covariate %in% colnames(hla_data)),
               msg = "covariate variables can not be found in hla_data"
        )
      } else {
        TRUE
      }
    },
    is.string(correction),
    is.flag(exponentiate)
  )

  alleles <- backquote(attr(hla_data, "alleles"))

  response <- backquote(response)
  if (substitute(model) %in% c("coxph", "cph")) {
    response <- paste0("Surv(", paste(response, collapse = ","), ")")
  }

  if (is.null(covariate)) {
    covariate <- . ~ .
  } else {
    covariate <- backquote(covariate)
  }

  model_function <- hlaAssocModel(model = model,
                                  response = response,
                                  variable = covariate,
                                  data = hla_data,
                                  ...
  )

  results <- map_dfr(
    .x = alleles,
    .f = ~ tidy(updateModel(object = model_function,
                            x = .,
                            backquote = FALSE,
                            collapse = " + "),
                exponentiate = exponentiate
    )
  )
  results <- mutate(results, term = gsub("`", "", term))
  results <- filter(results, checkAlleleFormat(term))
  results <- rename(results, allele = term)

  results <- mutate(results, p.adjusted = p.adjust(p.value, correction))

  return(results)
}

#' Association models for analysis of HLA alleles
#'
#' \code{hlaAssocModels} make creating statistical models for HLA association
#' analyses easier.
#'
#' @inheritParams analyzeHlaAssociations
#' @param response String specifying response variable in \code{data}.
#' @param variable Character specifying variables in \code{data}.
#' @param data Data frame containing the variables in the model.
#' @param ... Further arguments passed to \code{model} function.
#'
#' @return Fit from specified \code{model} function.
#'
#' @examples
#' library("survival")
#'
#' hla_calls_file <- system.file("extdata", "HLAHD_output_example.txt", package = "MiDAS")
#' hla_calls <- readHlaCalls(hla_calls_file)
#' pheno_file <- system.file("extdata", "pheno_example.txt", package = "MiDAS")
#' pheno <- read.table(pheno_file, header = TRUE)
#' covar_file <- system.file("extdata", "covar_example.txt", package = "MiDAS")
#' covar <- read.table(covar_file, header = TRUE)
#' hla_data <- prepareHlaData(hla_calls, pheno, covar, inheritance_model = "additive")
#' hlaAssocModel(model = "coxph",
#'               response = "Surv(OS, OS_DIED)",
#'               variable = c("AGE", "SEX"),
#'               data = hla_data
#' )
#'
#' @importFrom assertthat assert_that is.string see_if
#' @importFrom stats as.formula
#' @importFrom purrr is_formula
#'
#' @export
hlaAssocModel <- function(model,
                          response,
                          variable,
                          data,
                          ...) {
  assert_that(
    see_if(is.string(model) | is.function(model),
           msg = "model have to be a string (a length one character vector) or a function"
    ),
    if (is.string(model)) {
      see_if(is.function(get0(model)),
             msg = sprintf("could not find function %s", model)
      )
    } else {
      TRUE
    },
    see_if(is.string(response) | is_formula(response),
           msg = "response have to be a string or formula"
    ),
    see_if(is.character(variable) | is_formula(variable),
           msg = "variable have to be a character or formula"
    ),
    is.data.frame(data)
  )

  if (is.string(model)) {
    model <- as.name(model)
  }

  if (is.character(response)) {
    response <- paste0(response, " ~ 1")
    response <- as.formula(response)
  }

  if (is.character(variable)) {
    variable <- paste(variable, collapse = " + ")
    variable <- paste(". ~ .", variable, sep = " + ")
  }

  form <- update(response, variable)
  model_fun <- eval.parent(substitute(model(formula = form, data = data, ...)))

  # Assert that model_fun is in fact model fit object
  assert_that(
    see_if(is.object(model_fun),
           msg = sprintf("object returned by %s doesn't have OBJECT bit set",
                         deparse(substitute(model))
           )
    ),
    {
      object_call <- get0("call", envir = as.environment(model_fun))
      if (! is.null(object_call)) {
        object_formula <- eval(substitute(formula, env = as.list(object_call)))
        see_if(is_formula(object_formula),
               msg = sprintf("object returned by %s is not a model with defined formula",
                             deparse(substitute(model))
               )
        )
      } else {
        structure(FALSE,
                  msg = sprintf("object returned by %s doesn't have an attribute 'call'",
                                deparse(substitute(model))
                  )
        )
      }
    }
  )

  return(model_fun)
}

#' Stepwise conditional variables selection
#'
#' \code{stepwiseConditionalSelection} performs stepwise conditional testing
#' adding the previous top-associated variable as covariate, until there’s no
#' more significant variables based on a self-defined threshold.
#'
#' Selection criteria is the p-value from the test on coefficients values.
#'
#' @inheritParams checkHlaCallsFormat
#' @inheritParams analyzeHlaAssociations
#' @inheritParams hlaCallsToCounts
#' @param th number specifying p-value threshold for a term to be included into
#'   model.
#' @param keep logical flag indicating if the output should be a list of models
#'   resulting from each selection step. Default is to return only the final
#'   model.
#' @param rss_th number specifying residual sum of squares threshold at which
#'   function should stop adding additional terms.
#' @param ... Further arguments passed to \code{model} function.
#'
#' As the residual sum of squares approaches \code{0} the perfect fit is
#' obtained making further attempts at model selection nonsense, thus function
#' is stopped. This behavior can be controlled using \code{rss_th}.
#'
#' @return selected model or list of models. See \code{keep} parameter.
#'
#' @examples
#' library("survival")
#'
#' hla_calls_file <- system.file("extdata", "HLAHD_output_example.txt", package = "MiDAS")
#' hla_calls <- readHlaCalls(hla_calls_file)
#' pheno_file <- system.file("extdata", "pheno_example.txt", package = "MiDAS")
#' pheno <- read.table(pheno_file, header = TRUE)
#' covar_file <- system.file("extdata", "covar_example.txt", package = "MiDAS")
#' covar <- read.table(covar_file, header = TRUE)
#' hla_data <- prepareHlaData(hla_calls = hla_calls,
#'                            pheno = pheno,
#'                            covar = covar,
#'                            inheritance_model = "additive"
#' )
#' forwardConditionalSelection(model = "coxph",
#'                             hla_data = hla_data,
#'                             response = c("OS", "OS_DIED"),
#'                             covariate = c("AGE", "SEX"),
#'                             th = 0.05,
#'                             keep = FALSE,
#'                             rss_th = 1e-07
#' )
#'
#' @importFrom assertthat assert_that is.flag is.number is.string see_if
#' @importFrom MASS addterm
#' @importFrom purrr is_formula map_dfr
#' @importFrom rlang warn
#' @importFrom stats formula resid update
#' @export
forwardConditionalSelection <- function(model,
                                        hla_data,
                                        response,
                                        covariate,
                                        th,
                                        keep = FALSE,
                                        rss_th = 1e-07,
                                        ...) {

  assert_that(
    see_if(is.string(model) | is.function(model),
           msg = "model have to be a string (a length one character vector) or a function"
    ),
    if (is.string(model)) {
      see_if(is.function(get0(model)),
             msg = sprintf("could not find function %s", model)
      )
    } else {
      TRUE
    },
    is.character(response),
    {
      response_len <- length(response)
      is_cox <- as.character(substitute(model)) %in% c("coxph", "cph")
      if (is_cox & response_len != 2) {
        structure(FALSE, msg = "cox survival analysis requires response to be a character vector of length 2")
      } else if (! is_cox & response_len != 1) {
        structure(FALSE, msg = "response is not a string (a length one character vector).")
      } else {
        TRUE
      }
    },
    see_if(all(response %in% colnames(hla_data)),
           msg = "response variables can not be found in hla_data"
    ),
    see_if(is.character(covariate) | is.null(covariate),
           msg = "covariate have to be a character or NULL"
    ),
    {
      if (is.character(covariate)) {
        see_if(all(covariate %in% colnames(hla_data)),
               msg = "covariate variables can not be found in hla_data"
        )
      } else {
        TRUE
      }
    },
    is.number(th),
    is.flag(keep),
    is.number(rss_th)
  )

  alleles <- attr(hla_data, "alleles")

  response <- backquote(response)
  if (substitute(model) %in% c("coxph", "cph")) {
    response <- paste0("Surv(", paste(response, collapse = ","), ")")
  }

  if (is.null(covariate)) {
    covariate <- . ~ .
  } else {
    covariate <- backquote(covariate)
  }

  object <- hlaAssocModel(model = model,
                          response = response,
                          variable = covariate,
                          data = hla_data
  )

  vars <- alleles
  prev_formula <- formula(object)
  prev_vars <- all.vars(prev_formula)

  best <- list(object)
  i <- 2

  while (TRUE) {
    new_vars <- vars[! vars %in% prev_vars]

    results <- map_dfr(
      .x = new_vars,
      .f = ~ tidy(updateModel(object = object,
                              x = .,
                              backquote = TRUE,
                              collapse = " + "
      ))
    )
    results <- results[results$term %in% backquote(new_vars), ]
    results <- results[! is.infinite(results$p.value), ]

    i_min <- which.min(results$p.value)
    if (length(i_min) == 0) break
    if (results$p.value[i_min] > th) break

    object <- updateModel(object,
                          results$term[i_min],
                          backquote = FALSE,
                          collapse = " + "
    )

    if (sum(resid(object) ^ 2) <= rss_th) {
      warn("Perfect fit was reached attempting further model selection is nonsense.")
      break
    }
    prev_formula <- formula(object)
    prev_vars <- all.vars(prev_formula)
    best[[i]] <- object
    i <- i + 1
  }

  if (! keep) {
    best <- best[[length(best)]]
  }

  return(best)
}

#' Prepare data for statistical analysis
#'
#' \code{prepareHlaData} binds HLA alleles calls data frame with phenotypic
#' observations and covariates, creating an input data for further statistical
#' analysis.
#'
#' @inheritParams checkHlaCallsFormat
#' @inheritParams analyzeHlaAssociations
#' @inheritParams hlaCallsToCounts
#' @param pheno Data frame holding phenotypic response variables.
#' @param covar Data frame holding covariates or NULL.
#'
#' \code{pheno} and \code{covar} should be data frames with first column holding
#' samples IDs and named \code{ID}. Those should correspond to \code{ID} column
#' in \code{hla_calls}.
#'
#' @return Data frame with hla counts and pheno, covar. It also holds names of
#'   variables under attributes: 'alleles', 'response', 'covariate'.
#'
#' @examples
#' hla_calls_file <- system.file("extdata", "HLAHD_output_example.txt", package = "MiDAS")
#' hla_calls <- readHlaCalls(hla_calls_file)
#' pheno_file <- system.file("extdata", "pheno_example.txt", package = "MiDAS")
#' pheno <- read.table(pheno_file, header = TRUE)
#' covar_file <- system.file("extdata", "covar_example.txt", package = "MiDAS")
#' covar <- read.table(covar_file, header = TRUE)
#' prepareHlaData(hla_calls = hla_calls,
#'                pheno = pheno,
#'                covar = covar,
#'                inheritance_model = "additive"
#' )
#'
#' @importFrom assertthat assert_that is.string see_if
#'
#' @export
prepareHlaData <- function(hla_calls,
                           pheno,
                           covar = NULL,
                           inheritance_model = "additive") {

  assert_that(
    checkHlaCallsFormat(hla_calls),
    checkAdditionalData(pheno, hla_calls),
    checkAdditionalData(covar, hla_calls, accept.null = TRUE),
    is.string(inheritance_model),
    see_if(
      pmatch(inheritance_model,
             table = c("dominant", "recessive", "additive"),
             nomatch = 0
      ) != 0,
      msg = "inheritance_model should be one of 'dominant', 'recessive', 'additive'"
    )
  )

  hla_counts <- hlaCallsToCounts(hla_calls,
                                 inheritance_model = inheritance_model
  )

  assert_that(
    see_if(
      anyDuplicated(
        c(
          colnames(hla_counts[, -1]), colnames(pheno[, -1]), colnames(covar[, -1])
        )
      ) == 0,
      msg = "some colnames in hla_calls and pheno and covar duplicated"
    ))

  data <- left_join(hla_counts, pheno, by = "ID")
  if (! is.null(covar)) {
    data <- left_join(data, covar, by = "ID")
  }

  pheno_var <- colnames(pheno)[-1]
  covar_var <- colnames(covar)[-1]
  alleles_var <- colnames(hla_counts)[-1]
  alleles_freq <- getHlaFrequencies(hla_calls)

  hla_data <- structure(
    data,
    response = pheno_var,
    covariate = covar_var,
    alleles = alleles_var,
    alleles_freq = alleles_freq,
    hla_calls = hla_calls,
    inheritance_model = inheritance_model
  )
  class(hla_data) <- c("data.frame", "hla_data")

  return(hla_data)
}
