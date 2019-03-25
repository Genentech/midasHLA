#' Association analysis of HLA allele calls
#'
#' \code{analyzeHlaAssociations} performs associations analysis on HLA alleles
#' calls using statistical model of choice.
#'
#' @inheritParams checkHlaCallsFormat
#' @param model String specifying statistical model to use or corresponding
#'   function.
#' @param pheno Data frame holding phenotypic response variables.
#' @param covar Data frame holding covariates.
#' @param zygo Flag indicating whether zygosity should be added to
#'   covariates. See details for further explanations.
#' @param reduce_counts Flag indicating whether allele counts should be reduced
#'   to presence / absence indicators. See details for further explanations.
#' @param correction String specifying multiple testing correction method. See
#'   details for further information.
#'
#' \code{pheno} and \code{covar} should be data frames with first column holding
#' samples IDs and named \code{ID}. Those should correspond to \code{ID} column
#' in \code{hla_calls}.
#'
#' \code{zygo} indicate if additional covariate, indicating sample zygosity
#' status, should be added to covariates. HLA allele counts for each sample
#' can take following values \code{0, 1, 2}. To avoid implying ordering on those
#' levels and effect size, this information can be split between two variables.
#' If \code{zygo} is set to \code{TRUE} zygosity variable is added during model
#' fitting, it specifies if sample is homozygous for an allele.
#'
#' If \code{reduce_counts} is set to \code{TRUE} HLA allele counts are reduced
#' to presence / absence indicators. This is done by setting counts for
#' homozygotes as \code{1}.
#'
#' \code{correction} specifies p-value adjustment method to use, common choice
#' is Benjamini & Hochberg (1995) (\code{"BH"}). Internally this is passed to
#' \link{p.adjust}. Check there to get more details.
#'
#' @return Tibble containing combined results for all alleles in
#'   \code{hla_calls}.
#'
#' @examples
#' hla_calls_file <- system.file("extdata", "HLAHD_output_example.txt", package = "MiDAS")
#' hla_calls <- readHlaCalls(hla_calls_file)
#' pheno_file <- system.file("extdata", "pheno_example.txt", package = "MiDAS")
#' pheno <- read.table(pheno_file, header = TRUE)
#' covar_file <- system.file("extdata", "covar_example.txt", package = "MiDAS")
#' covar <- read.table(covar_file, header = TRUE)
#'
#' # Cox proportional hazards regression model
#' analyzeHlaAssociations(model = "coxph",
#'                        hla_calls,
#'                        pheno,
#'                        covar,
#'                        zygo = FALSE,
#'                        reduce_counts = FALSE,
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
                                   hla_calls,
                                   pheno,
                                   covar,
                                   zygo = FALSE,
                                   reduce_counts = FALSE,
                                   correction = "BH",
                                   exponentiate = FALSE) {
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
    checkHlaCallsFormat(hla_calls),
    checkAdditionalData(pheno, hla_calls),
    checkAdditionalData(covar, hla_calls, accept.null = TRUE),
    is.flag(zygo),
    is.flag(reduce_counts),
    is.string(correction)
  )

  hla_data <- prepareHlaData(hla_calls,
                             pheno,
                             covar,
                             zygo = FALSE,
                             reduce_counts = FALSE
  )
  alleles <- backquote(hla_data$alleles)
  response <- backquote(hla_data$response)
  if (substitute(model) %in% c("coxph", "cph")) {
    response <- paste0("Surv(", paste(response, collapse = ","), ")")
  }
  covariate <- backquote(hla_data$covariate)
  model_function <- hlaAssocModel(model = model,
                                  response = response,
                                  variable = covariate,
                                  data = hla_data$data
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
#' @param data Data frame containing variables in the model.
#'
#' @return Fit from specified \code{model} function.
#'
#' @examples
#' hla_calls_file <- system.file("extdata", "HLAHD_output_example.txt", package = "MiDAS")
#' hla_calls <- readHlaCalls(hla_calls_file)
#' pheno_file <- system.file("extdata", "pheno_example.txt", package = "MiDAS")
#' pheno <- read.table(pheno_file, header = TRUE)
#' covar_file <- system.file("extdata", "covar_example.txt", package = "MiDAS")
#' covar <- read.table(covar_file, header = TRUE)
#' hla_data <- prepareHlaData(hla_calls, pheno, covar)
#' hlaAssocModel(model = "coxph",
#'               response = hla_data$response,
#'               variable = hla_data$covariate,
#'               data = hla_data$data
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
        see_if(is_formula(formula(object_call)),
               msg = sprintf("object returned by %s is not a model with defined formula",
                             deparse(substitute(model))
               )
        )
      } else {
        structure(FALSE,
                  msg = sprintf("object returned by %s doesn't have an attribue 'call'",
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
#' @return selected model of the same class as \code{object} or list of models.
#'   See \code{keep} parameter.
#'
#' @examples
#' hla_calls_file <- system.file("extdata", "HLAHD_output_example.txt", package = "MiDAS")
#' hla_calls <- readHlaCalls(hla_calls_file)
#' pheno_file <- system.file("extdata", "pheno_example.txt", package = "MiDAS")
#' pheno <- read.table(pheno_file, header = TRUE)
#' covar_file <- system.file("extdata", "covar_example.txt", package = "MiDAS")
#' covar <- read.table(covar_file, header = TRUE)
#' forwardConditionalSelection(model = "coxph",
#'                             hla_calls = hla_calls,
#'                             pheno = pheno,
#'                             covar = covar,
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
                                        hla_calls,
                                        pheno,
                                        covar,
                                        th,
                                        keep = FALSE,
                                        rss_th = 1e-07) {

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
    checkHlaCallsFormat(hla_calls),
    checkAdditionalData(pheno, hla_calls),
    checkAdditionalData(covar, hla_calls, accept.null = TRUE),
    is.number(th),
    is.flag(keep),
    is.number(rss_th)
  )

  hla_data <- prepareHlaData(hla_calls, # Perhaps it would be beneficial to take this object outside, it has relatively simple structure so if someone needs to hack it it should be easy. Plus checking the data you are putting into functions would be beneficial.
                             pheno,
                             covar,
                             zygo = FALSE,
                             reduce_counts = FALSE
  )

  alleles <- hla_data$alleles
  response <- backquote(hla_data$response)
  if (substitute(model) %in% c("coxph", "cph")) {
    response <- paste0("Surv(", paste(response, collapse = ","), ")")
  }
  covariate <- backquote(hla_data$covariate)
  object <- hlaAssocModel(model = model,
                          response = response,
                          variable = covariate,
                          data = hla_data$data
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
#' @inheritParams checkHlaCallsFormat
#' @inheritParams analyzeHlaAssociations
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
#'                zygo = FALSE,
#'                reduce_counts = FALSE
#' )
#'
#' @importFrom assertthat assert_that is.flag see_if
#'
#' @export
prepareHlaData <- function(hla_calls,
                           pheno,
                           covar = NULL,
                           zygo = FALSE,
                           reduce_counts = FALSE) {
  assert_that(
    checkHlaCallsFormat(hla_calls),
    checkAdditionalData(pheno, hla_calls),
    checkAdditionalData(covar, hla_calls, accept.null = TRUE),
    is.flag(zygo),
    is.flag(reduce_counts)
  )

  hla_counts <- hlaCallsToCounts(hla_calls)
  zygosity <- hla_counts

  assert_that(
    see_if(
      anyDuplicated(
        c(
          colnames(hla_counts[, -1]), colnames(pheno[, -1]),
          colnames(covar[, -1]), paste0(colnames(zygosity[, -1]), "_zygosity")
        )
      ) == 0,
      msg = "some colnames in hla_calls and pheno and covar and zygosity are duplicated"
    ))

  if (reduce_counts) {
    hla_counts[, -1] <- lapply(hla_counts[, -1],
                               function(x) ifelse(x == 2, 1, x)
    )
  }

  data <- left_join(hla_counts, pheno, by = "ID")
  data <- left_join(data, covar, by = "ID")

  pheno_var <- colnames(pheno)[-1]
  covar_var <- colnames(covar)[-1]
  alleles_var <- colnames(hla_counts)[-1]

  if (zygo) {
    zygosity[, -1] <- lapply(zygosity[, -1], function(x) ifelse(x == 2, 1, 0))
    colnames(zygosity) <- c("ID", paste0(colnames(zygosity[, -1]), "_zygosity"))
    data <- left_join(data, zygosity, by = "ID")
    zygo_var <- paste0(colnames(hla_counts)[-1], "_zygosity")
    covar_var <- append(covar_var, zygo_var)
  }

  hla_data <- list(
    data = data,
    response = pheno_var,
    covariate = covar_var,
    alleles = alleles_var
  )

  return(hla_data)
}
