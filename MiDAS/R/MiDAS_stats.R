#' Association analysis
#'
#' \code{analyzeAssociations} perform association analysis on single variable
#' level using statistical model of choice.
#'
#' \code{correction} specifies p-value adjustment method to use, common choice
#' is Benjamini & Hochberg (1995) (\code{"BH"}). Internally this is passed to
#' \link[stats]{p.adjust}.
#'
#' @inheritParams updateModel
#' @param variables Character vector specifying variables to use in association
#'   tests.
#' @param placeholder String specyfing term in \code{object}'s formula which
#'   should be substituted with an allele during analysis.
#' @param correction String specifying multiple testing correction method. See
#'   details for further information.
#' @param n_correction Integer specifying number of comparisons to consider
#'   during multiple testing correction calculations, must be at least equal to
#'   number of comparisons being made; only set this (to non-default) when you
#'   know what you are doing!
#' @param exponentiate Logical flag indicating whether or not to exponentiate
#'   the coefficient estimates. Internally this is passed to
#'   \code{\link[broom]{tidy}}. This is typical for logistic and multinomial
#'   regressions, but a bad idea if there is no log or logit link. Defaults to
#'   FALSE.
#'
#' @return Tibble containing combined results for all alleles in
#'   \code{hla_calls}.
#'
#' @family MiDAS statistical functions
#' @seealso \code{\link[stats]{p.adjust}}, \code{\link[broom]{tidy}}
#'
#' @examples
#' \dontrun{
#' library("survival")
#' hla_calls_file <- system.file("extdata", "HLAHD_output_example.txt", package = "MiDAS")
#' hla_calls <- readHlaCalls(hla_calls_file)
#' pheno_file <- system.file("extdata", "pheno_example.txt", package = "MiDAS")
#' pheno <- read.table(pheno_file, header = TRUE)
#' covar_file <- system.file("extdata", "covar_example.txt", package = "MiDAS")
#' covar <- read.table(covar_file, header = TRUE)
#' midas_data <- prepareMiDAS(hla_calls = hla_calls,
#'                                pheno = pheno,
#'                                covar = covar,
#'                                experiment = "hla_allele",
#'                                inheritance_model = "additive"
#' )
#'
#' # Cox proportional hazards regression model
#' ## define base model with response and covariates
#' object <- coxph(Surv(OS, OS_DIED) ~ AGE + SEX + term, data = midas_data)
#'
#' ## test for alleles associations
#' analyzeAssociations(object = object,
#'                     variables = c("B*14:02", "DRB1*11:01")
#' )
#' }
#'
#' @importFrom assertthat assert_that see_if is.string
#' @importFrom broom tidy
#' @importFrom dplyr bind_rows
#' @importFrom stats p.adjust
#'
#' @export
analyzeAssociations <- function(object,
                                variables,
                                placeholder = "term",
                                correction = "bonferroni",
                                n_correction = NULL,
                                exponentiate = FALSE) {
  assert_that(
    checkStatisticalModel(object),
    hasTidyMethod(class(object)[1L])
  )
  object_call <- getCall(object)
  object_env <- attr(object$terms, ".Environment")
  object_data <- eval(object_call[["data"]], envir = object_env)
  object_variables <- colnames(object_data)[-1]

  assert_that(
    is.character(variables),
    see_if(
      all(test_vars <- variables %in% object_variables) | is.null(variables),
      msg = sprintf("%s can not be found in object data",
                    paste(variables[! test_vars], collapse = ", ")
      )
    ),
    is.string(placeholder),
    objectHasPlaceholder(object, placeholder),
    is.string(correction),
    isCountOrNULL(n_correction),
    isTRUEorFALSE(exponentiate)
  )

  results <- lapply(
    X = variables,
    FUN = function(x) tryCatch(
      expr = updateModel(
        object = object,
        x = x,
        placeholder = placeholder,
        backquote = TRUE,
        collapse = " + "
      ),
      error = function(e) {
        msg <- sprintf(
          "Error occurred while processing variable %s:\n\t%s",
          x,
          conditionMessage(e)
        )
        warn(msg)

        return(object)
      }
    )
  )

  results <- lapply(results, tidy, exponentiate = exponentiate)
  results <- bind_rows(results)
  results$term <- gsub("`", "", results$term)
  results <- results[results$term %in% variables, ]

  nc <- ifelse(is.null(n_correction), length(results$p.value), n_correction)
  assert_that(
    nc >= length(results$p.value),
    msg = sprintf("n_correction must be at least %i.", length(results$p.value))
  )
  results$p.adjusted <- p.adjust(
    p = results$p.value,
    method = correction,
    n = nc
  )

  #  This covariates were added for consistiency with conditional analyze, now however that we are filtering covariates there it doesn't make much sense to keep those?
  #  covariates <- formula(object)[[3]]
  #  covariates <- deparse(covariates)
  #  results$covariates <- covariates

  if (nrow(results) == 0) {
    warn("None of the variables could be tested. Returning empty table.")
  }

  return(results)
}

#' Stepwise conditional association analysis
#'
#' \code{analyzeConditionalAssociations} perform stepwise conditional testing
#' adding the previous top-associated variable as covariate, until there is no
#' more significant variables based on a self-defined threshold.
#'
#' Selection criteria is the p-value from the test on coefficients values.
#'
#' @inheritParams updateModel
#' @inheritParams analyzeAssociations
#' @param th Number specifying p-value threshold for a variable to be considered
#'   significant.
#' @param keep Logical flag indicating if the output should be a list of results
#'   resulting from each selection step. Default is to return only the final
#'   result.
#' @param rss_th Number specifying residual sum of squares threshold at which
#'   function should stop adding additional variables. As the residual sum of
#'   squares approaches \code{0} the perfect fit is obtained making further
#'   attempts at variables selection nonsense, thus function is stopped. This
#'   behavior can be controlled using \code{rss_th}.
#'
#' @return Tibble with stepwise conditional testing results.
#'
#' @family MiDAS statistical functions
#' @seealso \code{\link[stats]{p.adjust}}, \code{\link[broom]{tidy}}
#'
#' @examples
#' \dontrun{
#' library("survival")
#' hla_calls_file <- system.file("extdata", "HLAHD_output_example.txt", package = "MiDAS")
#' hla_calls <- readHlaCalls(hla_calls_file)
#' pheno_file <- system.file("extdata", "pheno_example.txt", package = "MiDAS")
#' pheno <- read.table(pheno_file, header = TRUE, stringsAsFactors = FALSE)
#' covar_file <- system.file("extdata", "covar_example.txt", package = "MiDAS")
#' covar <- read.table(covar_file, header = TRUE, stringsAsFactors = FALSE)
#' midas_data <- prepareMiDAS(hla_calls = hla_calls,
#'                              pheno = pheno,
#'                              covar = covar,
#'                              experiment = "hla_allele",
#'                              inheritance_model = "additive"
#' )
#'
#' ## define base model with covariates only
#' object <- coxph(Surv(OS, OS_DIED) ~ AGE + SEX + term, data = midas_data)
#' analyzeConditionalAssociations(object,
#'                             variables = c("B*14:02", "DRB1*11:01"),
#'                             th = 0.05,
#'                             rss_th = 1e-07
#' )
#' }
#'
#' @importFrom assertthat assert_that is.number is.string
#' @importFrom dplyr bind_rows tibble
#' @importFrom purrr map_dfr
#' @importFrom rlang warn
#' @importFrom stats formula resid
#'
#' @export
analyzeConditionalAssociations <- function(object,
                                           variables,
                                           placeholder = "term",
                                           correction = "bonferroni",
                                           n_correction = NULL,
                                           th,
                                           keep = FALSE,
                                           rss_th = 1e-07,
                                           exponentiate = FALSE) {
  assert_that(
    checkStatisticalModel(object),
    hasTidyMethod(class(object)[1L])
  )
  object_call <- getCall(object)
  object_env <- attr(object$terms, ".Environment")
  object_formula <- eval(object_call[["formula"]], envir = object_env)
  object_data <- eval(object_call[["data"]], envir = object_env)
  object_variables <- colnames(object_data)[-1]

  assert_that(
    is.character(variables),
    see_if(
      all(test_vars <- variables %in% object_variables),
      msg = sprintf("%s can not be found in object data",
                    paste(variables[! test_vars], collapse = ", ")
      )
    ),
    is.string(placeholder),
    objectHasPlaceholder(object, placeholder),
    is.string(correction),
    isCountOrNULL(n_correction),
    is.number(th),
    isTRUEorFALSE(keep),
    is.number(rss_th),
    isTRUEorFALSE(exponentiate)
  )

  prev_formula <- object_formula
  first_variables <- all.vars(object_formula)
  prev_variables <- first_variables
  new_variables <- variables[! variables %in% prev_variables]

  best <- list()
  i <- 1

  while (length(new_variables) > 0) {
    results <- lapply(
      X = new_variables,
      FUN = function(x) tryCatch(
        expr = updateModel(
          object = object,
          x = x,
          placeholder = placeholder,
          backquote = TRUE,
          collapse = " + "
        ),
        error = function(e) {
          msg <- sprintf(
            "Error occurred while processing variable %s:\n\t%s",
            x,
            conditionMessage(e)
          )
          warn(msg)

          return(object)
        }
      )
    )
    results <- lapply(results, tidy, exponentiate = exponentiate)
    results <- bind_rows(results)

    mask_new_vars <- backquote(results[["term"]]) %in% backquote(new_variables)
    results <- results[mask_new_vars, ]

    nc <- ifelse(is.null(n_correction), length(results$p.value), n_correction)
    assert_that(
      nc >= length(results$p.value),
      msg = sprintf("n_correction must be at least %i.", length(results$p.value))
    )
    results$p.adjusted <- p.adjust(
      p = results$p.value,
      method = correction,
      n = nc
    )

    results <- results[! is.infinite(results[["p.value"]]), ]

    mask <- ! prev_variables %in% first_variables
    results$covariates <- paste(prev_variables[mask], collapse = " + ")

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
    new_variables <- variables[! variables %in% prev_variables]

    results$term <- gsub("`", "", results$term)
    best[[i]] <- results
    i <- i + 1
  }

  if (length(best) == 0) {
    warn("No significant variables found. Returning empty table.") # Tibble to be more precise?
  }

  if (keep) {
    results <- best
  } else {
    if (length(best) == 0) {
      results <- results[0, ]
    } else {
      results <- lapply(best, function(res) {
        i_min <- which.min(res[["p.value"]])
        res[i_min, ]
      })
      results <- bind_rows(results)
    }
  }

  return(results)
}

#' Amino acid position omnibus test
#'
#' \code{aaPosOmnibusTest} calculate overall p-value for amino acid position
#' using likelihood rato test.
#'
#' Likelihood ratio test is conducted by comparing model given in \code{object}
#' with extended model, that is created by including effect of amino acid
#' residues at each position in \code{aa_pos}. Amino acid effect is assumend to
#' be simply additive.
#'
#' @inheritParams analyzeAssociations
#' @inheritParams summariseAAPosition
#'
#' @return Data frame containing omnibus test results for specified amino acid
#'   positions.
#'
#' @importFrom assertthat assert_that see_if is.string
#' @importFrom broom tidy
#' @importFrom dplyr bind_cols mutate select
#' @importFrom stats p.adjust
#' @importFrom rlang .data
#'
#' @export
aaPosOmnibusTest <- function(object,
                             aa_pos,
                             correction = "bonferroni",
                             n_correction = NULL) {
  assert_that(
    checkStatisticalModel(object)
  )
  object_call <- getCall(object)
  object_env <- attr(object$terms, ".Environment")
  object_data <- eval(object_call[["data"]], envir = object_env)
  object_variables <- colnames(object_data)[-1]
  base_vars <- all.vars(object_call$formula)

  assert_that(
    is.character(aa_pos), # proper check on aa_pos format is missing!
    is.string(correction),
    isCountOrNULL(n_correction)
  )

  variables <- lapply(
    X = aa_pos,
    FUN = function(x) {
      pattern <- paste0(x, "_")
      resids <- grep(
        pattern = pattern,
        x = object_variables,
        fixed = TRUE,
        value = TRUE
      )
      resids <- resids[! resids %in% base_vars] # discard residues used as covariates
      assert_that(
        length(resids) != 0,
        msg = sprintf("amino acid position %s could not be found.", x)
      )

      return(resids)
    }
  )

  results <- lapply(
    X = variables,
    FUN = function(x) tryCatch(
      expr = LRTest(
        object,
        updateModel(
          object = object,
          x = x,
          backquote = TRUE,
          collapse = " + "
        )
      ),
      error = function(e) {
        msg <- sprintf(
          "Error occurred while processing variable %s:\n\t%s",
          gsub("_[A-Z]", "", x[1]), # output aa_pos in err message
          conditionMessage(e)
        )
        warn(msg)

        return(object)
      }
    )
  )

  results <- bind_rows(results) %>%
    mutate(
      residues = gsub("[A-Z]+_[0-9]+_", "", .data$term),
      term = aa_pos
    ) %>%
    select(aa_pos = .data$term,
           .data$residues,
           d.f. = .data$dof,
           .data$statistic,
           .data$p.value
    )

  nc <- ifelse(is.null(n_correction), length(results$p.value), n_correction)
  assert_that(
    nc >= length(results$p.value),
    msg = sprintf("n_correction must be at least %i.", length(results$p.value))
  )
  results$p.adjusted <- p.adjust(
    p = results$p.value,
    method = correction,
    n = nc
  )

  return(results)
}

#' Run MiDAS statistical analysis
#'
#' \code{runMiDAS} perform association analysis on MiDAS data using mode and
#' statistical model specified by user. Function is intended for use with
#' \code{\link{prepareMiDAS}}. See examples section.
#'
#' \code{experiment} is used to select an experiment from \code{MiDAS} object
#' associated with \code{object} using. In standard work flow data are first
#' processed using \code{\link{prepareMiDAS}} creating experiments with
#' transformed data for given \code{experiment}. Next user constructs the
#' statistical model using function of choice (eg. \code{lm}. \code{coxph}).
#' Than \code{runMiDAS} is used to evaluate specified model uder \code{mode}
#' of choice. See the details of different \code{mode}'s implmenetions for
#' more details.
#'
#' \code{correction} specifies p-value adjustment method to use, common choice
#' is Benjamini & Hochberg (1995) (\code{"BH"}). Internally this is passed to
#' \link[stats]{p.adjust}. Check there to get more details.
#'
#' @inheritParams analyzeAssociations
#' @inheritParams filterByFrequency
#' @param experiment String indicating the experiment associated with
#'   \code{object}'s \code{MiDAS} data to use. Valid values includes:
#'   \code{"hla_allele"}, \code{"aa_level"}, \code{"allele_g_group"},
#'   \code{"allele_supertype"}, \code{"allele_group"}, \code{"kir_genes"},
#'   \code{"hla_kir_interactions"}. See \code{link{prepareMiDAS}} for more
#'   informations.
#' @param conditional Logical flag,
#' @param omnibus Logical flag.
#' @param ... other arguments
#'
#' @return Tibble containing analysis results.
#'
#' @family MiDAS statistical functions
#' @seealso \code{\link[stats]{p.adjust}}, \code{\link[broom]{tidy}}
#'
#' @examples
#' \dontrun{
#' # read hla calls file
#' hla_calls_file <- system.file("extdata", "HLAHD_output_example.txt", package = "MiDAS")
#' hla_calls <- readHlaCalls(hla_calls_file)
#'
#' # read kir calls file
#' kir_calls_file <- system.file("extdata", "KIP_output_example.txt", package = "MiDAS")
#' kir_calls <- readKirCalls(kir_calls_file, counts = TRUE)
#'
#' # read phenotypic data and covariates
#' pheno_file <- system.file("extdata", "pheno_example.txt", package = "MiDAS")
#' pheno <- read.table(pheno_file, header = TRUE, stringsAsFactors = FALSE)
#' covar_file <- system.file("extdata", "covar_example.txt", package = "MiDAS")
#' covar <- read.table(covar_file, header = TRUE, stringsAsFactors = FALSE)
#' phenotype <- left_join(pheno, covar, by ="ID")
#'
#' # create MiDAS object
#' midas <- prepareMiDAS(hla_calls = hla_calls,
#'                       kir_calls = kir_calls,
#'                       colData = phenotype,
#'                       inheritance_model = "additive",
#'                       experiment = "hla_allele"
#' )
#'
#' # constructs statistical model
#' object <- lm(OS ~ AGE + SEX + term, data = midas)
#'
#' # run analysis
#' runMiDAS(object, mode = "linear", experiment = "hla_allele")
#' }
#'
#' @importFrom assertthat assert_that is.number is.string
#'
#' @export
runMiDAS <- function(object,
                     experiment,
                     conditional = FALSE,
                     omnibus = FALSE,
                     lower_frequency_cutoff = NULL,
                     upper_frequency_cutoff = NULL,
                     correction = "bonferroni",
                     n_correction = NULL,
                     exponentiate = FALSE,
                     ...
                    ) {
  assert_that(
    checkStatisticalModel(object),
    hasTidyMethod(class(object)[1L])
  )
  object_details <- getObjectDetails(object)

  assert_that(
    isClass(object_details$data, "MiDAS"),
    validObject(object_details$data),
    objectHasPlaceholder(object, getPlaceholder(object_details$data)),
    is.string(experiment),
    stringMatches(experiment, choice = getExperiments(object_details$data)),
    isTRUEorFALSE(conditional),
    isTRUEorFALSE(omnibus),
    validateFrequencyCutoffs(lower_frequency_cutoff, upper_frequency_cutoff),
    is.string(correction),
    isCountOrNULL(n_correction),
    isTRUEorFALSE(exponentiate)
  )

  if (! is.null(lower_frequency_cutoff) || ! is.null(upper_frequency_cutoff)) {
    object_details$data <-
      filterByFrequency(
        object = object_details$data,
        experiment = experiment,
        lower_frequency_cutoff = lower_frequency_cutoff,
        upper_frequency_cutoff = upper_frequency_cutoff
      )
    assert_that(
      length(object_details$data[, , experiment]) != 0,
      msg = "No variables available for analysis, please revisit your filtration criteria."
    )
  }

  args <- list(
    object_details = object_details,
    experiment = experiment,
    correction = correction,
    n_correction = n_correction,
    exponentiate = exponentiate,
    ...
  )

  results <- if (! conditional && ! omnibus) {
    do.call(runMiDAS_linear, args)
  } else if (conditional && ! omnibus) {
    do.call(runMiDAS_conditional, args)
  }

  return(results)
}

#' @rdname runMiDAS
#'
#' @title runMiDAS linear
#'
#' @details statistical analysis is performed iteratively
#'   on each variable in selected experiment. This is done by substituting
#'   \code{placeholder} in the \code{object}'s formula with each variable in the
#'   experiment.
#'
#' @param object_details TODO
#'
#' @importFrom assertthat assert_that
#' @importFrom dplyr rename
#' @importFrom rlang call_modify !! :=
#'
runMiDAS_linear <- function(object_details,
                            experiment,
                            correction = "bonferroni",
                            n_correction = NULL,
                            exponentiate = FALSE,
                            ...) {
  # only experiments of class matrix can be used here
  assert_that(
    isClass(object_details$data[[experiment]], "matrix"),
    msg = sprintf("Unconditional runMiDAS does not support experiment %s",
                  experiment)
  )

  test_var <- rownames(object_details$data[[experiment]])
  placeholder <- getPlaceholder(object_details$data)

  # insert data for analysis
  data <- midasToWide(object_details$data, experiment)
  call <- call_modify(object_details$call, data = data)
  object <- eval(call)

  # run analysis
  results <- analyzeAssociations(object,
                                 variables = test_var,
                                 placeholder = placeholder,
                                 correction = correction,
                                 n_correction = n_correction,
                                 exponentiate = exponentiate
  )

  assert_that(
    nrow(results) > 0,
    msg = "Could not process any variables. Please check warning messages for more informations (warnings())."
  )

  # format linear results
  ## add variables frequencies
  if (typeof(object_details$data[[experiment]]) == "integer") {
    results <-
      left_join(
        results,
        runMiDASGetVarsFreq(
          midas = object_details$data,
          experiment = experiment,
          test_covar = object_details$formula_vars[1]
        ),
        by = "term"
      )
  }

  ## rename term
  term_name <- switch (experiment,
                       "hla_allele" = "allele",
                       "aa_level" = "aa",
                       "expression_level" = "allele",
                       "allele_g_group" = "g.group",
                       "allele_supertype" = "supertype",
                       "allele_group" = "allele.group",
                       "kir_genes" = "kir.gene",
                       "hla_kir_interactions" = "hla.kir.interaction",
                       "term"
  )
  results <- rename(results, !! term_name := .data$term)

  return(results)
}

#' @rdname runMiDAS
#'
#' @title runMiDAS conditional
#'
#' @details statistical analysis is performed in a
#'   stepwise conditional testing manner, adding the previous top-associated
#'   variable as a covariate to \code{object}'s formula. The analysis stops
#'   when there is no more siginifcant variabls, based on self-defined
#'   threshold. The p-values of variables are used as the selection criteria.
#'   This proces is reapeated for each variable available in the selected
#'   experiment. This is done by substituting \code{placeholder} in the
#'   \code{object}'s formula with each variable in the experiment.
#'
#' @inheritParams analyzeConditionalAssociations
#'
#' @importFrom assertthat assert_that
#' @importFrom dplyr rename bind_rows
#' @importFrom rlang call_modify !! :=
#'
runMiDAS_conditional <- function(object_details,
                                 experiment,
                                 correction = "bonferroni",
                                 n_correction = NULL,
                                 exponentiate = FALSE,
                                 th = 0.05,
                                 keep = FALSE,
                                 rss_th = 1e-07,
                                 ...) {
  # only experiments of class matrix can be used here
  assert_that(
    isClass(object_details$data[[experiment]], "matrix"),
    msg = sprintf("Conditional runMiDAS does not supported experiment %s",
                  experiment)
  )

  # get test covariates names
  test_var <- rownames(object_details$data[[experiment]])
  placeholder <- getPlaceholder(object_details$data)

  # insert data for analysis
  data <- midasToWide(object_details$data, experiment)
  call <- call_modify(object_details$call, data = data)
  object <- eval(call)

  # run analysis
  results <- analyzeConditionalAssociations(
    object,
    variables = test_var,
    placeholder = placeholder,
    correction = correction,
    n_correction = n_correction,
    th = th,
    keep = keep,
    rss_th = rss_th,
    exponentiate = exponentiate
  )

  assert_that(
    ! (nrow(results) == 0 && ! keep) || ! (length(results) == 0 && keep),
    msg = "Could not process any variables. Please check warning messages for more informations (warnings())."
  )

  # format conditional results
  term_name <- switch (experiment,
                       "hla_allele" = "allele",
                       "aa_level" = "aa",
                       "expression_level" = "allele",
                       "allele_g_group" = "g.group",
                       "allele_supertype" = "supertype",
                       "allele_group" = "allele.group",
                       "kir_genes" = "kir.gene",
                       "hla_kir_interactions" = "hla.kir.interaction",
                       "term"
  )

  if (typeof(object_details$data[[experiment]]) == "integer") {
    if (keep) {
      results <- lapply(
        X = results,
        FUN = function(x) {
          x <-
            left_join(
              x,
              runMiDASGetVarsFreq(
                midas = object_details$data,
                experiment = experiment,
                test_covar = object_details$formula_vars[1]
              ),
              by = "term"
            )
          rename(.data = x, !! term_name := .data$term)
        }
      )
    } else {
      results <-
        left_join(
          results,
          runMiDASGetVarsFreq(
            midas = object_details$data,
            experiment = experiment,
            test_covar = object_details$formula_vars[1]
          ),
          by = "term"
        )
      results <- rename(results, !! term_name := .data$term)
    }
  }

  return(results)
}
