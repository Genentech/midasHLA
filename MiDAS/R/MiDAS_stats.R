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
#'                                analysis_type = "hla_allele",
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
#'                              analysis_type = "hla_allele",
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

    results <- results[results[["term"]] %in% backquote(new_variables), ]

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

  if (keep) {
    results <- best
  } else {
    if (length(best) == 0) {
      warn("No significant variables found. Returning empty table.") # Tibble to be more precise?
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

#' Analyze MiDAS data
#'
#' \code{runMiDAS} perform association analysis on MiDAS data using
#' statistical model specified by user. Function is intended for use with
#' \code{\link{prepareMiDAS}}. See examples section.
#'
#' \code{analysis_type} is used to select variables from data associated with
#' \code{object} using \code{\link[Hmisc]{label}}s. In standard work flow data
#' are first processed using \code{\link{prepareMiDAS}}, columns of its
#' output data frame are labeled with the type of analysis they can be used for
#' eg. \code{hla_allele}. By specifying \code{analysis_type} function will
#' select all variables with corresponding label. This choice can be further
#' refined by using \code{pattern} argument or extended with \code{variables}.
#' In cases where one would rather prefer to specify more custom set of
#' variables analysis type \code{"none"} can be used. In this mode function will
#' only use variables specifed by \code{variables} argument.
#'
#' \code{correction} specifies p-value adjustment method to use, common choice
#' is Benjamini & Hochberg (1995) (\code{"BH"}). Internally this is passed to
#' \link[stats]{p.adjust}. Check there to get more details.
#'
#' @inheritParams analyzeAssociations
#' @inheritParams analyzeConditionalAssociations
#' @inheritParams formatAssociationsResults
#' @param analysis_type String indicating the type of analysis to be performed,
#'   it's used to select appropriate variables for testing from the data
#'   associated with \code{object}. Valid values are \code{"hla_allele"},
#'   \code{"aa_level"}, \code{"expression_level"}, \code{"allele_g_group"},
#'   \code{"allele_supertype"}, \code{"allele_group"}, \code{"kir_genes"},
#'   \code{"hla_kir_interactions"}, \code{"none"}. See details for further
#'   explanations.
#' @param pattern String containing a regular expression that is used
#'   to further select variables selected by \code{analysis_type}.
#' @param conditional Logical flag indicating if the analysis should be
#'   performed using stepwise conditional test. See
#'   \code{\link{analyzeConditionalAssociations}} for more details.
#' @param variables Character vector specifying additional variables to use in
#'   association tests except those selected by \code{analysis_type}. By default
#'   \code{NULL}. Note that additional variables are not considered when
#'   applying lower or upper frequency cutoff.
#' @param lower_frequency_cutoff Number specifying lower threshold for inclusion
#'   of a variable. If it's a number between \code{0} and \code{1} variables
#'   with frequency below this number will not be considered during analysis. If
#'   it's greater or equal \code{1} variables with number of counts less that
#'   this will not be considered during analysis. Only applied to discrete
#'   variables.
#' @param upper_frequency_cutoff Number specifying upper threshold for inclusion
#'   of a variable. If it's a number between \code{0} and \code{1} variables
#'   with frequency above this number will not be considered during analysis. If
#'   it's greater or equal \code{1} variables with number of counts greater that
#'   this will not be considered during analysis.Only applied to discrete
#'   variables.
#' @param exponentiate Logical flag indicating if coefficient estimates
#'   should be exponentiated. This is typical for logistic and multinomial
#'   regressions, but a bad idea if there is no log or logit link. If
#'   \code{NULL} function will try to figure this out by testing if response is
#'   binary (\code{0} or \code{1}).
#'
#' @return Tibble containing results for tested variables.
#'
#' @family MiDAS statistical functions
#' @seealso \code{\link[stats]{p.adjust}}, \code{\link[broom]{tidy}}
#'
#' @examples
#' library("survival")
#' hla_calls_file <- system.file("extdata", "HLAHD_output_example.txt", package = "MiDAS")
#' hla_calls <- readHlaCalls(hla_calls_file)
#' pheno_file <- system.file("extdata", "pheno_example.txt", package = "MiDAS")
#' pheno <- read.table(pheno_file, header = TRUE, stringsAsFactors = FALSE)
#' covar_file <- system.file("extdata", "covar_example.txt", package = "MiDAS")
#' covar <- read.table(covar_file, header = TRUE, stringsAsFactors = FALSE)
#' midas_data <- prepareMiDAS(hla_calls = hla_calls,
#'                                pheno = pheno,
#'                                covar = covar,
#'                                analysis_type = "hla_allele",
#'                                inheritance_model = "additive"
#' )
#'
#' object <- coxph(Surv(OS, OS_DIED) ~ AGE + SEX + term, data = midas_data)
#' runMiDAS(object, analysis_type = "hla_allele")
#'
#' @importFrom assertthat assert_that is.number is.string
#' @importFrom dplyr bind_rows filter left_join select rename
#' @importFrom stats getCall
#' @importFrom rlang !! := .data
#' @importFrom magrittr %>% %<>%
#' @importFrom Hmisc label
#'
#' @export
runMiDAS <- function(object,
                     analysis_type = c("hla_allele", "aa_level", "expression_level", "allele_g_group", "allele_supertype", "allele_group", "kir_genes", "hla_kir_interactions", "hla_divergence", "none"),
                     pattern = NULL,
                     variables = NULL,
                     placeholder = "term",
                     conditional = FALSE,
                     keep = FALSE,
                     lower_frequency_cutoff = NULL,
                     upper_frequency_cutoff = NULL,
                     pvalue_cutoff = NULL,
                     correction = "bonferroni",
                     n_correction = NULL,
                     exponentiate = FALSE,
                     th = 0.05,
                     rss_th = 1e-07) {
  assert_that(
    checkStatisticalModel(object),
    hasTidyMethod(class(object)[1L])
  )
  object_call <- getCall(object)
  object_env <- attr(object$terms, ".Environment")
  object_formula <- eval(object_call[["formula"]], envir = object_env)
  object_data <- eval(object_call[["data"]], envir = object_env)
  # assert object data more than one column, etc
  object_variables <- colnames(object_data)[-1]
  variables_labels <- label(object_data[, -1])
  # what happens if data is not labled?

  assert_that(
    is.string(analysis_type),
    stringMatches(analysis_type,
                  choice = c("hla_allele", "aa_level", "expression_level", "allele_g_group", "allele_supertype", "allele_group", "kir_genes", "hla_kir_interactions", "hla_divergence", "none")
    ),
    isStringOrNULL(pattern),
    isCharacterOrNULL(variables),
    is.string(placeholder),
    objectHasPlaceholder(object, placeholder),
    isTRUEorFALSE(conditional),
    isTRUEorFALSE(keep),
    see_if(
      all(test_vars <- variables %in% object_variables) | is.null(variables),
      msg = sprintf("%s can not be found in object data",
                    paste(variables[! test_vars], collapse = ", ")
      )
    ),
    isNumberOrNULL(lower_frequency_cutoff),
    isNumberOrNULL(upper_frequency_cutoff),
    isNumberOrNULL(pvalue_cutoff),
    is.string(correction),
    isCountOrNULL(n_correction),
    isTRUEorFALSE(exponentiate),
    is.number(th),
    is.number(rss_th)
  )

  assert_that(
    ! (analysis_type == "none" && is.null(variables)),
    msg = "For analysis type \"none\" variables argument can not be NULL."
  )

  assert_that(
    any(variables_labels == analysis_type, na.rm = TRUE) || ! is.null(variables),
    msg = "Argument variables = NULL can be used only with labeled variables, make sure to use prepareMiDAS function for data preparation."
  )

  # select variables based on analysis_type labels
  labeled_vars <- object_variables[variables_labels == analysis_type]
  if (! is.null(pattern)) {
    labeled_vars <- grep(pattern = pattern, x = labeled_vars, value = TRUE)
  }

  # create set of variables for further testing
  test_var <- c(labeled_vars, variables) %>%
    unique()
  test_var <- test_var[! test_var %in% all.vars(object_formula)]
  assert_that(length(test_var) != 0,
              msg = "No new variables found in object data."
  )

  # # guess if model used is logistic type
  # if (is.null(logistic)) {
  #   model_fun <- deparse(object_call[[1]])
  #   model_family <- deparse(object_call[["family"]])
  #   logistic <- grepl("coxph", model_fun) |
  #     (grepl("glm", model_fun) & grepl("binomial", model_family))
  # }

  # Filter variables on frequency cutoff
  variables_labels <- variables_labels[test_var] # if test_var != NULL select only corresponding labels
  mask_counts <- variables_labels %in% c(
    "hla_allele",
    "aa_level",
    "allele_g_group",
    "allele_supertype",
    "allele_group",
    "kir_genes",
    "hla_kir_interactions"
  ) & ! test_var %in% variables
  cts_vars <- test_var[mask_counts]
  ncts_vars <- test_var[! mask_counts]

  if (length(cts_vars)) {
    lower_frequency_cutoff <- ifelse(is.null(lower_frequency_cutoff), 0, lower_frequency_cutoff)
    upper_frequency_cutoff <- ifelse(is.null(upper_frequency_cutoff), Inf, upper_frequency_cutoff)
    variables_freq <- object_data %>%
      select("ID",!! cts_vars) %>%
      getCountsFrequencies() %>%
      rename(Ntotal = .data$Counts, Ntotal.frequency = .data$Freq) %>%
      filter(.data$Ntotal > lower_frequency_cutoff |
               lower_frequency_cutoff < 1) %>%
      filter(.data$Ntotal.frequency > lower_frequency_cutoff |
               lower_frequency_cutoff >= 1) %>%
      filter(.data$Ntotal < upper_frequency_cutoff |
               upper_frequency_cutoff < 1) %>%
      filter(.data$Ntotal.frequency < upper_frequency_cutoff |
               upper_frequency_cutoff >= 1)

    assert_that(
      nrow(variables_freq) != 0,
      msg = "No observations passes filtering criteria. Revisit your choice of 'lower_frequency_cutoff' and 'upper_frequency_cutoff'."
    )

    test_var <- c(ncts_vars, variables_freq$term)
  }

  if (conditional) {
    results_iter <- analyzeConditionalAssociations(object,
                                                   variables = test_var,
                                                   placeholder = placeholder,
                                                   correction = correction,
                                                   n_correction = n_correction,
                                                   th = th,
                                                   keep = TRUE,
                                                   rss_th = rss_th,
                                                   exponentiate = exponentiate
    )
    results <- lapply(results_iter, function(res) {
      i_min <- which.min(res[["p.value"]])
      res[i_min, ]
    })
    results <- bind_rows(results)
  } else {
    results <- analyzeAssociations(object,
                                   variables = test_var,
                                   placeholder = placeholder,
                                   correction = correction,
                                   n_correction = n_correction,
                                   exponentiate = exponentiate
    )
  }

  assert_that(
    nrow(results) > 0,
    msg = "Could not process any variables. Please check warning messages for more informations."
  )

  # Add frequency information to results table if there were any count variables
  if (length(cts_vars)) {
    results <- left_join(x = results, y = variables_freq, by = "term")
    if (conditional) {
      results_iter <- lapply(
        X = results_iter,
        FUN = left_join,
        y = variables_freq,
        by = "term")
    }
  }

  pheno_var <- all.vars(object_formula)[1]
  bin_pheno <- object_data[, pheno_var] %in% c(0, 1)
  bin_pheno <- all(bin_pheno, na.rm = TRUE)
  if (length(cts_vars) != 0 & bin_pheno) {
    pos_freq <- object_data %>%
      filter(.data[[!! pheno_var]] == 1) %>%
      select("ID", !! cts_vars) %>%
      getCountsFrequencies() %>%
      rename(Npositive = .data$Counts,
             Npositive.frequency = .data$Freq
      )
    results <- left_join(x = results, y = pos_freq, by = "term")
    if (conditional) {
      results_iter <- lapply(results_iter, left_join, y = pos_freq, by = "term")
    }

    neg_freq <- object_data %>%
      filter(.data[[!! pheno_var]] != 1) %>%
      select("ID", !! cts_vars) %>%
      getCountsFrequencies() %>%
      rename(Nnegative = .data$Counts,
             Nnegative.frequency = .data$Freq
      )
    results <- left_join(x = results,   y = neg_freq, by = "term")
    if (conditional && keep) {
      results_iter <- lapply(results_iter, left_join, y = neg_freq, by = "term")
    }
  }

  # rename term and estimate to match preety_table
  # estimate_name <- ifelse(logistic, "odds.ratio", "estimate") in favor of new policy :D
  term_name <- switch (analysis_type,
                       "hla_allele" = "allele",
                       "aa_level" = "aa",
                       "expression_level" = "allele",
                       "allele_g_group" = "g.group",
                       "allele_supertype" = "supertype",
                       "allele_group" = "allele.group",
                       "kir_genes" = "kir.gene",
                       "hla_kir_interactions" = "hla.kir.interaction",
                       "hla_divergence" = "hla.divergence",
                       "term"
  )

  if (conditional && keep) {
    results <- lapply(
      X = results_iter,
      FUN = function(x) {
        rename(.data = x,
               !! term_name := .data$term
               # !! estimate_name := .data$estimate
        )
      }
    )
  } else {
    results %<>%
     rename(!! term_name := .data$term)# , !! estimate_name := .data$estimate)
  }

  return(results)
}

#' Prepare MiDAS data for statistical analysis
#'
#' \code{prepareMiDAS} transform HLA alleles calls and KIR calls according
#' to selected analysis type and join obtained transformation with additional
#' data like phenotypic observations or covariates.
#'
#' Data frames passed as arguments to the function are joined using
#' \code{hla_calls} as a reference. As a consequence all samples with HLA data
#' are kept, independent of whether there's actually phenotypes / covariates.
#' Moreover phenotypes / covariates without HLA data are discarded.
#'
#' \code{...} should be data frames with first column holding sample IDs and
#' named \code{ID}. Those should correspond to \code{ID} column in
#' \code{hla_calls} and \code{kir_counts}.
#'
#' Choices for \code{analysis_type}:
#' \describe{
#'   \item{\code{hla_allele}}{
#'     \code{hla_calls} are transformed into counts under
#'     \code{inheritance_model} of choice (see \code{\link{hlaCallsToCounts}}
#'     for more details).
#'   }
#'   \item{\code{aa_level}}{
#'     \code{hla_calls} are first converted to amino acid level, taking only
#'     variable positions under consideration. Than variable amino acid
#'     positions are transformed to counts under \code{inheritance_model} of
#'     choice (see \code{\link{hlaToAAVariation}} and
#'     \code{\link{aaVariationToCounts}} for more details).
#'   }
#'   \item{\code{expression_level}}{
#'     \code{hla_calls} are transformed to expression levels using expression
#'     dictionaries shipped with package (see \code{\link{hlaToVariable}} for
#'     more details). Expression levels from both alleles are summed into single
#'     variable for each HLA gene.
#'   }
#'   \item{\code{allele_g_group}}{
#'     \code{hla_calls} are transformed to HLA alleles groups using G group
#'     dictionary shipped with the package. Than they are transformed to counts
#'     under \code{inheritance_model} of choice (see \code{\link{hlaToVariable}}
#'     and \code{\link{hlaCallsToCounts}} for more details).
#'   }
#' . \item{\code{allele_supertype}}{
#'     \code{hla_calls} are transformed to HLA alleles groups using supertypes
#'     dictionary shipped with the package. Than they are transformed to counts
#'     under \code{inheritance_model} of choice (see \code{\link{hlaToVariable}}
#'     and \code{\link{hlaCallsToCounts}} for more details).
#'   }
#'   \item{\code{allele_group}}{
#'     \code{hla_calls} are transformed to HLA alleles groups using Bw4/6, C1/2
#'     and Bw4+A23+A24+A32 dictionaries shipped with the package. Than they are
#'     transformed to counts under \code{inheritance_model} of choice
#'     (see \code{\link{hlaToVariable}} and \code{\link{hlaCallsToCounts}} for
#'     more details).
#'   }
#'   \item{\code{kir_genes}}{
#'     \code{kir_counts} data frame is joined with other inputs.
#'   }
#'   \item{\code{hla_kir_interactions}}{
#'     \code{hla_calls} are processed with \code{kir_counts} into HLA - KIR
#'     interactions variables (see \code{\link{getHlaKirInteractions}} for more
#'     details).
#'   }
#'   \item{\code{hla_divergence}}{
#'     Distances between Class I alleles are calculated using Grantham distance,
#'     as implemented in \code{hlaCallsGranthamDistance} function. Additionally
#'     average distance in Class I genese is calculated.
#'   }
#'   \item{\code{custom}}{
#'     No data transformation is done. All inputs are joined together.
#'   }
#' }
#'
#' @inheritParams checkHlaCallsFormat
#' @inheritParams hlaCallsToCounts
#' @inheritParams hlaToAAVariation
#' @param kir_counts Data frame with KIR genes counts. Required for
#'   \code{"kir_genes"} analysis type.
#' @param ... Data frames holding additional variables like phenotypic
#'   observations or covariates.
#' @param analysis_type Character vector indicating analysis type for which data
#'   should be prepared. Valid choices are \code{"hla_allele"},
#'   \code{"aa_level"}, \code{"expression_level"}, \code{"allele_group"},
#'   \code{"hla_divergence"}, \code{"custom"}. Each prepared variable will
#'   be labeled with corresponding
#'   \code{analysis_type}. See details for further explanations.
#' @param placeholder String specifying name of dummy column added to result
#'   data frame.
#'
#' @return Data frame containing prepared data.
#'
#' @family MiDAS statistical functions
#'
#' @examples
#' hla_calls_file <- system.file("extdata", "HLAHD_output_example.txt", package = "MiDAS")
#' hla_calls <- readHlaCalls(hla_calls_file)
#' pheno_file <- system.file("extdata", "pheno_example.txt", package = "MiDAS")
#' pheno <- read.table(pheno_file, header = TRUE)
#' covar_file <- system.file("extdata", "covar_example.txt", package = "MiDAS")
#' covar <- read.table(covar_file, header = TRUE)
#' prepareMiDAS(hla_calls, pheno, covar, analysis_type = "expression_level")
#'
#' @importFrom assertthat assert_that is.string see_if
#' @importFrom dplyr funs group_by left_join mutate summarise_all syms
#' @importFrom magrittr %>% %<>%
#' @importFrom rlang .data !!!
#' @importFrom stats runif
#' @importFrom tidyr gather spread
#' @importFrom Hmisc label label<-
#'
#' @export
prepareMiDAS <- function(hla_calls,
                             ...,
                             kir_counts = NULL,
                             analysis_type = c("hla_allele", "aa_level", "expression_level", "allele_g_group", "allele_supertype", "allele_group", "kir_genes", "hla_kir_interactions", "hla_divergence", "custom"),
                             inheritance_model = "additive",
                             placeholder = "term",
                             indels = TRUE,
                             unkchar = FALSE
) {
  assert_that(
    checkHlaCallsFormat(hla_calls),
    checkKirCountsFormat(kir_counts, accept.null = TRUE),
    is.character(analysis_type),
    characterMatches(
      x = analysis_type,
      choice = c("hla_allele", "aa_level", "expression_level", "allele_g_group", "allele_supertype", "allele_group", "kir_genes", "hla_kir_interactions", "hla_divergence", "custom")
    ),
    is.string(inheritance_model),
    is.string(placeholder),
    stringMatches(
      x = inheritance_model,
      choice = c("dominant", "recessive", "additive")
    ),
    isTRUEorFALSE(indels),
    isTRUEorFALSE(unkchar)
  )

  additional_data <- list(...)
  for (additional_data_frame in additional_data) { # rewrite those tests to addhere to current standard
    assert_that(
      checkAdditionalData(additional_data_frame, hla_calls, accept.null = FALSE)
    )
  }

  midas_data <- hla_calls[, 1, drop = FALSE]

  # Process hla_calls based on analysis type
  if ("hla_allele" %in% analysis_type) {
    hla_allele <- hlaCallsToCounts(
      hla_calls = hla_calls,
      inheritance_model = inheritance_model
    )

    label(hla_allele[-1], self = FALSE) <- rep("hla_allele", ncol(hla_allele) - 1)
    midas_data <- left_join(midas_data, hla_allele, by = "ID")
  }

  if ("aa_level" %in% analysis_type) {
    aa_level <- hlaToAAVariation(
      hla_calls = hla_calls,
      indels = indels,
      unkchar = unkchar
    ) %>%
      aaVariationToCounts(inheritance_model = inheritance_model)

    label(aa_level[-1], self = FALSE) <- rep("aa_level", ncol(aa_level) - 1)
    midas_data <- left_join(midas_data, aa_level, by = "ID")
  }

  if ("expression_level" %in% analysis_type) {
    lib <- listMiDASDictionaries()
    lib <- grep("expression", lib, value = TRUE)
    expression_level <- Reduce(
      f = function(...) left_join(..., by = "ID"),
      x = lapply(lib, hlaToVariable, hla_calls = hla_calls, na.value = NA)
    )

    assert_that(
      ncol(expression_level) > 1,
      msg = "no expression levels were found for input hla_calls"
    )

    expression_level %<>%
      gather("expression", "value", -c("ID")) %>%
      mutate(expression = gsub("_[12]", "", .data$expression)) %>%
      group_by(!!! syms(c("ID", "expression"))) %>%
      summarise_all(funs(sum)) %>%
      spread(.data$expression, .data$value, sep = NULL)

    label(expression_level[-1], self = FALSE) <- rep(
      x = "expression_level",
      ncol(expression_level) - 1
    )
    midas_data <- left_join(midas_data, expression_level, by = "ID")
  }

  if ("allele_g_group" %in% analysis_type) {
    lib <- "allele_HLA_Ggroup"
    allele_g_group <- hlaToVariable(hla_calls = hla_calls,
                                    dictionary = lib,
                                    na.value = 0
    )

    assert_that(
      ncol(allele_g_group) > 1,
      msg = "no allele could be assigned to G group for input hla_calls"
    )

    allele_g_group <- hlaCallsToCounts(
      hla_calls = allele_g_group,
      inheritance_model = inheritance_model,
      check_hla_format = FALSE
    )

    label(allele_g_group[-1], self = FALSE) <- rep(
      x = "allele_g_group",
      ncol(allele_g_group) - 1
    )
    midas_data <- left_join(midas_data, allele_g_group, by = "ID")
  }

  if ("allele_supertype" %in% analysis_type) {
    lib <- "allele_HLA_supertype"
    allele_supertype <- hlaToVariable(hla_calls = hla_calls,
                                      dictionary = lib,
                                      na.value = 0
    )

    assert_that(
      ncol(allele_supertype) > 1,
      msg = "no allele could be assigned to supertype for input hla_calls"
    )

    allele_supertype <- hlaCallsToCounts(
      hla_calls = allele_supertype,
      inheritance_model = inheritance_model,
      check_hla_format = FALSE
    ) %>%
      select(-"Unclassified")

    label(allele_supertype[-1], self = FALSE) <- rep(
      x = "allele_supertype",
      ncol(allele_supertype) - 1
    )
    midas_data <- left_join(midas_data, allele_supertype, by = "ID")
  }

  if ("allele_group" %in% analysis_type) {
    lib <- c(
      "allele_HLA-B_Bw",
      "allele_HLA_Bw4+A23+A24+A32",
      "allele_HLA-C_C1-2"
    )
    allele_group <- Reduce(
      f = function(...) left_join(..., by = "ID"),
      x = lapply(lib, hlaToVariable, hla_calls = hla_calls, na.value = 0)
    )

    assert_that(
      ncol(allele_group) > 1,
      msg = "no allele could be assigned to allele groups for input hla_calls"
    )

    allele_group <- hlaCallsToCounts(
      hla_calls = allele_group,
      inheritance_model = inheritance_model,
      check_hla_format = FALSE
    )

    label(allele_group[-1], self = FALSE) <- rep(
      x = "allele_group",
      ncol(allele_group) - 1
    )
    midas_data <- left_join(midas_data, allele_group, by = "ID")
  }

  if ("kir_genes" %in% analysis_type) {
    assert_that(
      ! missing(kir_counts),
      msg = "\"kir_genes\" analysis type requires kir_counts argument to be specified"
    )

    label(kir_counts[-1], self = FALSE) <- rep(
      x = "kir_genes",
      ncol(kir_counts) - 1
    )
    midas_data <- left_join(midas_data, kir_counts, by = "ID")
  }

  if ("hla_kir_interactions" %in% analysis_type) {
    assert_that(
      ! missing(kir_counts),
      msg = "\"hla_kir_interactions\" analysis type requires kir_counts argument to be specified"
    )

    hla_kir_interactions <- getHlaKirInteractions(
      hla_calls = hla_calls,
      kir_counts = kir_counts
    )

    label(hla_kir_interactions[-1], self = FALSE) <- rep(
      x = "hla_kir_interactions",
      ncol(hla_kir_interactions) - 1
    )

    midas_data <- left_join(midas_data, hla_kir_interactions, by = "ID")
  }

  if ("hla_divergence" %in% analysis_type) {
    hla_divergence <- hlaCallsGranthamDistance(
      hla_calls = hla_calls,
      genes = c("A", "B", "C")
    )
    hla_divergence$ABC_avg <- rowMeans(hla_divergence[-1])
    colnames(hla_divergence)[-1] <-
      paste0(colnames(hla_divergence[-1]), "_divergence")
    label(hla_divergence[-1], self = FALSE) <- rep(
      x = "hla_divergence",
      ncol(hla_divergence) - 1
    )

    midas_data <- left_join(midas_data, hla_divergence, by = "ID")
  }

  if ("custom" %in% analysis_type) {
    label(hla_calls[-1], self = FALSE) <- rep(
      x = "custom",
      ncol(hla_calls) - 1
    )
    midas_data <- left_join(midas_data, hla_calls, by = "ID")
  }

  # join with additional_data
  if (length(additional_data)) {
    midas_data <- Reduce(
      f = function(...) left_join(..., by = "ID"),
      x = additional_data,
      init = midas_data
    )
  }

  # add dummy column
  if (! is.null(placeholder)) {
    midas_data[[placeholder]] <- runif(nrow(midas_data))
  }

  return(midas_data)
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
#' @examples
#' hla_calls_file <- system.file("extdata", "HLAHD_output_example.txt", package = "MiDAS")
#' hla_calls <- readHlaCalls(hla_calls_file)
#' pheno_file <- system.file("extdata", "pheno_example.txt", package = "MiDAS")
#' pheno <- read.table(pheno_file, header = TRUE)
#' covar_file <- system.file("extdata", "covar_example.txt", package = "MiDAS")
#' covar <- read.table(covar_file, header = TRUE)
#' midas_data <- prepareMiDAS(hla_calls, pheno, covar, analysis_type = "aa_level")
#' object <- lm(OS ~ AGE + SEX, data = midas_data)
#' aaPosOmnibusTest(object, aa_pos = c("B_35", "E_128", "A_270"))
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
    is.character(aa_pos),
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
