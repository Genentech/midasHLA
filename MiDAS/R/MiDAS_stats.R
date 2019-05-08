#' Association analysis
#'
#' \code{analyzeAssociations} performs associations analysis on single variable
#' level using statistical model of choice.
#'
#' @inheritParams updateModel
#' @param variables Character specifying variables to use in association tests.
#' @param correction String specifying multiple testing correction method. See
#'   details for further information.
#' @param exponentiate Logical indicating whether or not to exponentiate the
#'   coefficient estimates. Internally this is passed to \link[broom]{tidy}.
#'   This is typical for logistic and multinomial regressions, but a bad idea if
#'   there is no log or logit link. Defaults to FALSE.
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
#'                     variables = c("B*14:02", "DRB1*11:01"),
#'                     correction = "BH"
#' )
#'
#' @importFrom assertthat assert_that see_if is.flag is.string
#' @importFrom broom tidy
#' @importFrom dplyr bind_rows left_join filter mutate rename
#' @importFrom stats p.adjust
#'
#' @export
analyzeAssociations <- function(object,
                                variables, # TODO check if variables are specified in object,  # opcja all
                                correction = "BH",
                                exponentiate = FALSE) {
  assert_that(
    checkStatisticalModel(object),
    is.character(variables),
    is.string(correction),
    is.flag(exponentiate)
  )

  results <- lapply(variables,
                    updateModel,
                    object = object,
                    backquote = TRUE,
                    collapse = " + "
  )

  results <- lapply(results, tidy, exponentiate = exponentiate)
  results <- bind_rows(results)
  results <- mutate(results, term = gsub("`", "", term))
  results <- filter(results, term %in% variables)

  results <- mutate(results, p.adjusted = p.adjust(p.value, correction))

  return(results)
}

#' Stepwise conditional variables selection
#'
#' \code{stepwiseConditionalSelection} performs stepwise conditional testing
#' adding the previous top-associated variable as covariate, until there’s no
#' more significant variables based on a self-defined threshold.
#'
#' Selection criteria is the p-value from the test on coefficients values.
#'
#' @inheritParams updateModel
#' @param variables Character specifying variables to use in selection
#'   procedure.
#' @param th number specifying p-value threshold for a term to be included into
#'   model.
#' @param keep logical flag indicating if the output should be a list of models
#'   resulting from each selection step. Default is to return only the final
#'   model.
#' @param rss_th number specifying residual sum of squares threshold at which
#'   function should stop adding additional terms.
#'
#' As the residual sum of squares approaches \code{0} the perfect fit is
#' obtained making further attempts at model selection nonsense, thus function
#' is stopped. This behavior can be controlled using \code{rss_th}.
#'
#' @return selected model or list of models. See \code{keep} parameter.
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
#' ## define base model with covariates only
#' object <- coxph(Surv(OS, OS_DIED) ~ AGE + SEX, data = midas_data)
#' forwardConditionalSelection(object,
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
forwardConditionalSelection <- function(object,
                                        variables,
                                        th,
                                        rss_th = 1e-07) {
  assert_that(
    checkStatisticalModel(object),
    is.character(variables),
    is.number(th),
    is.number(rss_th)
  )

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
      ))
    )
    results <- results[results$term %in% backquote(new_variables), ]
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
        obj_tidy <- tidy(obj)
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
