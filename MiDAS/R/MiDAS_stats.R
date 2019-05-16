#' Association analysis
#'
#' \code{analyzeAssociations} performs associations analysis on single variable
#' level using statistical model of choice.
#'
#' @inheritParams updateModel
#' @param variables Character specifying variables to use in association tests
#'   or \code{NULL}. If \code{NULL} all variables in object data are tested.
#'   See details for further information.
#' @param correction String specifying multiple testing correction method. See
#'   details for further information.
#' @param exponentiate Logical indicating whether or not to exponentiate the
#'   coefficient estimates. Internally this is passed to \link[broom]{tidy}.
#'   This is typical for logistic and multinomial regressions, but a bad idea if
#'   there is no log or logit link. Defaults to FALSE.
#'
#' \code{variables} takes \code{NULL} as a default value. When specifed as such
#' column names of data frame associated with the \code{object} are used as
#' variables for testing. This exludes first column which should corresponds
#' to samples IDs as well as covariates and response variables defined in
#' object formula.
#'
#' \code{correction} specifies p-value adjustment method to use, common choice
#' is Benjamini & Hochberg (1995) (\code{"BH"}). Internally this is passed to
#' \link[stats]{p.adjust}. Check there to get more details.
#'
#' @return Tibble containing combined results for all alleles in
#'   \code{hla_calls}.
#'
#' @examples
#' library("survival")
#' hla_calls_file <- system.file("extdata", "HLAHD_output_example.txt", package = "MiDAS")
#' hla_calls <- readHlaCalls(hla_calls_file)
#' pheno_file <- system.file("extdata", "pheno_example.txt", package = "MiDAS")
#' pheno <- read.table(pheno_file, header = TRUE)
#' covar_file <- system.file("extdata", "covar_example.txt", package = "MiDAS")
#' covar <- read.table(covar_file, header = TRUE)
#' midas_data <- prepareHlaData(hla_calls = hla_calls,
#'                              pheno = pheno,
#'                              covar = covar,
#'                              inheritance_model = "additive"
#' )
#'
#' # Cox proportional hazards regression model
#' ## define base model with response and covariates
#' object <- coxph(Surv(OS, OS_DIED) ~ AGE + SEX, data = midas_data)
#'
#' ## test for alleles associations
#' analyzeAssociations(object = object,
#'                     correction = "BH"
#' )
#'
#' @importFrom assertthat assert_that see_if is.flag is.string
#' @importFrom broom tidy
#' @importFrom dplyr bind_rows
#' @importFrom stats p.adjust
#'
#' @export
analyzeAssociations <- function(object,
                                variables = NULL,
                                correction = "BH",
                                exponentiate = FALSE) {
  assert_that(
    checkStatisticalModel(object)
  )
  object_call <- getCall(object)
  object_formula <- eval(object_call[["formula"]], envir = parent.frame())
  object_data <- eval(object_call[["data"]], envir = parent.frame())
  object_variables <- colnames(object_data)[-1]

  assert_that(
    see_if(
      is.character(variables) | is.null(variables),
      msg = "variables is not a character vector or NULL"
    ),
    see_if(
      all(test_vars <- variables %in% object_variables) | is.null(variables),
      msg = sprintf("%s can not be found in object data",
                    paste(variables[! test_vars], collapse = ", ")
      )
    ),
    is.string(correction),
    is.flag(exponentiate)
  )

  if (is.null(variables)) {
    mask <- ! object_variables %in% all.vars(object_formula)
    variables <- object_variables[mask]
  }

  results <- lapply(variables,
                    updateModel,
                    object = object,
                    backquote = TRUE,
                    collapse = " + "
  )

  results <- lapply(results, tidy, exponentiate = exponentiate)
  results <- bind_rows(results)
  results$term <- gsub("`", "", results$term)
  results <- results[results$term %in% variables, ]

  results$p.adjusted <- p.adjust(results$p.value, correction)

  covariates <- formula(object)[[3]]
  covariates <- deparse(covariates)
  results$covariates <- covariates

  return(results)
}

#' Stepwise conditional association analysis
#'
#' \code{analyzeConditionalAssociations} performs stepwise conditional testing
#' adding the previous top-associated variable as covariate, until there is no
#' more significant variables based on a self-defined threshold.
#'
#' Selection criteria is the p-value from the test on coefficients values.
#'
#' @inheritParams updateModel
#' @inheritParams analyzeAssociations
#' @param th number specifying p-value threshold for a variable to be considered
#'   significant.
#' @param rss_th number specifying residual sum of squares threshold at which
#'   function should stop adding additional variables. As the residual sum of
#'   squares approaches \code{0} the perfect fit is obtained making further
#'   attempts at variables selection nonsense, thus function is stopped. This
#'   behavior can be controlled using \code{rss_th}.
#'
#' @return tibble with stepwise conditional testsing results.
#'
#' @examples
#' library("survival")
#' hla_calls_file <- system.file("extdata", "HLAHD_output_example.txt", package = "MiDAS")
#' hla_calls <- readHlaCalls(hla_calls_file)
#' pheno_file <- system.file("extdata", "pheno_example.txt", package = "MiDAS")
#' pheno <- read.table(pheno_file, header = TRUE, stringsAsFactors = FALSE)
#' covar_file <- system.file("extdata", "covar_example.txt", package = "MiDAS")
#' covar <- read.table(covar_file, header = TRUE, stringsAsFactors = FALSE)
#' midas_data <- prepareHlaData(hla_calls = hla_calls,
#'                              pheno = pheno,
#'                              covar = covar,
#'                              inheritance_model = "additive"
#' )
#'
#' ## define base model with covariates only
#' object <- coxph(Surv(OS, OS_DIED) ~ AGE + SEX, data = midas_data)
#' analyzeConditionalAssociations(object,
#'                             variables = c("B*14:02", "DRB1*11:01"),
#'                             th = 0.05,
#'                             rss_th = 1e-07
#' )
#'
#' @importFrom assertthat assert_that is.number
#' @importFrom dplyr bind_rows tibble
#' @importFrom purrr map_dfr
#' @importFrom rlang warn
#' @importFrom stats formula resid
#'
#' @export
analyzeConditionalAssociations <- function(object,
                                           variables = NULL,
                                           th,
                                           rss_th = 1e-07,
                                           exponentiate = FALSE) {
  assert_that(
    checkStatisticalModel(object)
  )
  object_call <- getCall(object)
  object_formula <- eval(object_call[["formula"]], envir = parent.frame())
  object_data <- eval(object_call[["data"]], envir = parent.frame())
  object_variables <- colnames(object_data)[-1]

  assert_that(
    see_if(
      is.character(variables) | is.null(variables),
      msg = "variables is not a character vector or NULL"
    ),
    see_if(
      all(test_vars <- variables %in% object_variables) | is.null(variables),
      msg = sprintf("%s can not be found in object data",
                    paste(variables[! test_vars], collapse = ", ")
      )
    ),
    is.number(th),
    is.number(rss_th)
  )

  if (is.null(variables)) {
    mask <- ! object_variables %in% all.vars(object_formula)
    variables <- object_variables[mask]
  }

  prev_formula <- formula(object)
  prev_variables <- all.vars(prev_formula)

  best <- list()
  i <- 1

  while (TRUE) {
    new_variables <- variables[! variables %in% prev_variables]

    results <- map_dfr(
      .x = new_variables,
      .f = ~ tidy(updateModel(object = object,
                              x = .,
                              backquote = TRUE,
                              collapse = " + "
                  )
      )
    )
    results <- results[results[["term"]] %in% backquote(new_variables), ]
    eesults <- results[! is.infinite(results[["p.value"]]), ]

    i_min <- which.min(results[["p.value"]])
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
    prev_variables <- all.vars(prev_formula)
    best[[i]] <- object
    i <- i + 1
  }

  if (length(best) > 0) {
    results <- lapply(
      X = best,
      FUN = function(obj) {
        cov <- formula(obj)[[3]]
        cov <- all.vars(cov)
        cov <- cov[-length(cov)]
        obj_tidy <- tidy(obj, exponentiate = exponentiate)
        obj_tidy <- obj_tidy[length(cov) + 1, ]
        obj_tidy$covariates <- paste(cov, collapse = " + ")
        return(obj_tidy)
      }
    )
    results <- bind_rows(results)
    results$term <- gsub("`", "", results$term)
    results <- results[results$term %in% variables, ]
  } else {
    results <- tibble()
  }

  return(results)
}

#' Prepare data for statistical analysis
#'
#' \code{prepareHlaData} binds HLA alleles calls data frame with additional data
#' frames like phenotypic observations or covariates, creating an input data for
#' further statistical analysis.
#'
#' @inheritParams checkHlaCallsFormat
#' @inheritParams hlaCallsToCounts
#' @param ... Data frames holding additional variables like phenotypic
#'   observations or covariates.
#'
#' \code{...} should be data frames with first column holding samples IDs and
#' named \code{ID}. Those should correspond to \code{ID} column in
#' \code{hla_calls}.
#'
#' @return Data frame with hla counts and additional variables.
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
#' @importFrom dplyr left_join
#'
#' @export
prepareHlaData <- function(hla_calls,
                           ...,
                           inheritance_model = "additive") {

  additional_data <- list(...)
  if (length(additional_data) == 0) {
    additional_data <- NULL
  }

  assert_that(
    checkHlaCallsFormat(hla_calls),
    all(
      vapply(
        X = additional_data,
        FUN = function(x) {
          additional_data_frame <- x
          checkAdditionalData(additional_data_frame,
                              hla_calls = hla_calls,
                              accept.null = TRUE
          )
        },
        FUN.VALUE = logical(1)
      )
    ),
    see_if(
      anyDuplicated(
        c(
          unique(unlist(hla_calls[, -1])),
          unlist(lapply(
            X = additional_data,
            FUN = function(x) colnames(x)[-1]
          ))
        )
      ) == 0,
      msg = "column names in additional data are duplicated or overlap with alleles in hla_calls"
    ),
    is.string(inheritance_model),
    see_if(
      pmatch(inheritance_model,
             table = c("dominant", "recessive", "additive"),
             nomatch = 0
      ) != 0,
      msg = "inheritance_model should be one of 'dominant', 'recessive', 'additive'"
    )
  )

  data <- hlaCallsToCounts(hla_calls,
                           inheritance_model = inheritance_model
  )
  if (length(additional_data) != 0) {
    for (i in 1:length(additional_data)) {
      data <- left_join(data, additional_data[[i]], by = "ID")
    }
  }

  return(data)
}
