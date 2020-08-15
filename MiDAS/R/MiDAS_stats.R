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
#'                                experiment = "hla_alleles",
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

  results <- iterativeModel(object, placeholder, variables, exponentiate)

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
#' @param th_adj Logical
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
#'                              experiment = "hla_alleles",
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
                                           th_adj = TRUE,
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
    results <- iterativeModel(object, placeholder, new_variables, exponentiate)

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

    # select optimization criteria
    crit <- ifelse(th_adj, "p.adjusted", "p.value")

    i_min <- which.min(results[[crit]])
    if (length(i_min) == 0) break
    if (results[[crit]][i_min] > th) break

    object <- updateModel(object,
                          results$term[i_min],
                          backquote = TRUE,
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

#' Omnibus test
#'
#' \code{OmnibusTest} calculate overall p-value for linear combination of
#' vairables using likelihood rato test.
#'
#' Likelihood ratio test is conducted by comparing model given in \code{object}
#' with extended model, that is created by including effect of variables given
#' in \code{variables} as thie linear combination.
#'
#' @inheritParams analyzeAssociations
#' @param omnibus_groups List of character vectors
#'
#' @return Data frame containing omnibus test results for specified variables.
#'
#' @importFrom dplyr bind_cols
#' @importFrom stats p.adjust
#'
omnibusTest <- function(object,
                        omnibus_groups,
                        placeholder = "term",
                        correction = "bonferroni",
                        n_correction = NULL) {
  results <- iterativeLRT(object, placeholder, omnibus_groups)

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
#'   \code{"hla_alleles"}, \code{"hla_aa"}, \code{"hla_g_groups"},
#'   \code{"hla_supertypes"}, \code{"hla_NK_ligands"}, \code{"kir_genes"},
#'   \code{"hla_kir_interactions"}. See \code{link{prepareMiDAS}} for more
#'   informations.
#' @param conditional Logical flag,
#' @param omnibus Logical flag.
#' @param omnibus_groups_filter Character
#' @param ... other arguments
#'
#' @return Tibble containing analysis results.
#'
#' @family MiDAS statistical functions
#' @seealso \code{\link[stats]{p.adjust}}, \code{\link[broom]{tidy}}
#'
#' @examples
#' \dontrun{
#' # create MiDAS object
#' midas <- prepareMiDAS(hla_calls = MiDAS_tut_HLA,
#'                       kir_calls = MiDAS_tut_KIR,
#'                       colData = MiDAS_tut_pheno,
#'                       experiment = "hla_alleles"
#' )
#'
#' # constructs statistical model
#' object <- lm(disease ~ outcome + term, data = midas)
#'
#' # run analysis
#' runMiDAS(object, experiment = "hla_alleles", inheritance_model = "dominant")
#' }
#'
#' @importFrom assertthat assert_that is.number is.string
#' @importFrom rlang warn
#'
#' @export
runMiDAS <- function(object,
                     experiment,
                     inheritance_model = NULL,
                     conditional = FALSE,
                     omnibus = FALSE,
                     omnibus_groups_filter = NULL,
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
    isStringOrNULL(inheritance_model),
    isTRUEorFALSE(conditional),
    isTRUEorFALSE(omnibus),
    isCharacterOrNULL(omnibus_groups_filter),
    validateFrequencyCutoffs(lower_frequency_cutoff, upper_frequency_cutoff),
    is.string(correction),
    isCountOrNULL(n_correction),
    isTRUEorFALSE(exponentiate)
  )

  # convert experiment to specifed inheritance model
  if (! is.null(inheritance_model)) {
    if (isInheritanceModelApplicable(experiment)) {
      assert_that(
        characterMatches(inheritance_model, c("dominant", "recessive", "additive"))
      )
      object_details$data[[experiment]] <- applyInheritanceModel(
        experiment = object_details$data[[experiment]],
        inheritance_model = inheritance_model
      )
    } else {
      warn(sprintf("Inheritance model can not be applied to experiment %s. Continuing without it.", experiment))
    }
  }

  if (! is.null(omnibus_groups_filter)) {
    object_details$data <-
      filterByOmnibusGroups(
        object = object_details$data,
        experiment = experiment,
        groups = omnibus_groups_filter
      )
  }

  if (! is.null(lower_frequency_cutoff) || ! is.null(upper_frequency_cutoff)) {
    object_details$data <-
      filterByFrequency(
        object = object_details$data,
        experiment = experiment,
        lower_frequency_cutoff = lower_frequency_cutoff,
        upper_frequency_cutoff = upper_frequency_cutoff
      )
  }

  assert_that(
    length(object_details$data[, , experiment]) != 0,
    msg = "No variables available for analysis, please revisit your filtration criteria."
  )

  args <- list(
    call = object_details$call,
    midas = object_details$data,
    experiment = experiment,
    test_covar = object_details$formula_vars[1],
    correction = correction,
    n_correction = n_correction,
    exponentiate = exponentiate,
    ...
  )

  results <- if (! conditional && ! omnibus) {
    do.call(runMiDAS_linear, args)
  } else if (conditional && ! omnibus) {
    do.call(runMiDAS_conditional, args)
  } else if (! conditional && omnibus) {
    do.call(runMiDAS_linear_omnibus, args)
  } else if (conditional && omnibus) {
    do.call(runMiDAS_conditional_omnibus, args)
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
#' @importFrom assertthat assert_that
#' @importFrom dplyr arrange rename
#' @importFrom rlang call_modify !! :=
#'
runMiDAS_linear <- function(call,
                            midas,
                            experiment,
                            test_covar,
                            correction = "bonferroni",
                            n_correction = NULL,
                            exponentiate = FALSE,
                            ...) {
  test_var <- rownames(midas[[experiment]])
  placeholder <- getPlaceholder(midas)

  # insert data for analysis
  data <- as.data.frame(midas)
  call <- call_modify(substitute(call), data = data)
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
  if (isExperimentCountsOrZeros(midas[[experiment]])) {
    results <-
      left_join(
        results,
        runMiDASGetVarsFreq(
          midas = midas,
          experiment = experiment,
          test_covar = test_covar
        ),
        by = "term"
      )
  }

  ## rename term
  term_name <- switch (experiment,
                       "hla_alleles" = "allele",
                       "hla_aa" = "aa",
                       "expression_level" = "allele",
                       "hla_g_groups" = "g.group",
                       "hla_supertypes" = "supertype",
                       "hla_NK_ligands" = "allele.group",
                       "kir_genes" = "kir.gene",
                       "hla_kir_interactions" = "hla.kir.interaction",
                       "term"
  )
  results <- rename(results, !! term_name := .data$term) %>%
    arrange(.data$p.value)

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
#' @importFrom dplyr arrange rename bind_rows
#' @importFrom rlang call_modify !! :=
#'
runMiDAS_conditional <- function(call,
                                 midas,
                                 experiment,
                                 test_covar,
                                 correction = "bonferroni",
                                 n_correction = NULL,
                                 exponentiate = FALSE,
                                 th = 0.05,
                                 th_adj = TRUE,
                                 keep = FALSE,
                                 rss_th = 1e-07,
                                 ...) {
  # get test covariates names
  test_var <- rownames(midas[[experiment]])
  placeholder <- getPlaceholder(midas)

  # insert data for analysis
  data <- as.data.frame(midas)
  call <- call_modify(substitute(call), data = data)
  object <- eval(call)

  # run analysis
  results <- analyzeConditionalAssociations(
    object,
    variables = test_var,
    placeholder = placeholder,
    correction = correction,
    n_correction = n_correction,
    th = th,
    th_adj = th_adj,
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
                       "hla_alleles" = "allele",
                       "hla_aa" = "aa",
                       "expression_level" = "allele",
                       "hla_g_groups" = "g.group",
                       "hla_supertypes" = "supertype",
                       "hla_NK_ligands" = "allele.group",
                       "kir_genes" = "kir.gene",
                       "hla_kir_interactions" = "hla.kir.interaction",
                       "term"
  )

  if (isExperimentCountsOrZeros(midas[[experiment]])) {
    if (keep) {
      results <- lapply(
        X = results,
        FUN = function(x) {
          x <-
            left_join(
              x,
              runMiDASGetVarsFreq(
                midas = midas,
                experiment = experiment,
                test_covar = test_covar
              ),
              by = "term"
            )
          rename(.data = x, !! term_name := .data$term) %>%
            arrange(p.value)
        }
      )
    } else {
      results <-
        left_join(
          results,
          runMiDASGetVarsFreq(
            midas = midas,
            experiment = experiment,
            test_covar = test_covar
          ),
          by = "term"
        )
      results <- rename(results, !! term_name := .data$term) %>%
        arrange(.data$p.value)
    }
  }

  return(results)
}

#' @rdname runMiDAS
#'
#' @title runMiDAS linear omnibus
#'
#' @details statistical analysis is performed iteratively on groups of variables
#'   like residues at particular amino acid position, using likelyhood ratio
#'   test. This is done by substituting \code{placeholder} in the
#'   \code{object}'s formula with linear combination of variable in particular
#'   group.
#'
#' @param object_details TODO
#'
#' @importFrom assertthat assert_that
#' @importFrom dplyr as_tibble mutate select
#' @importFrom magrittr %>%
#' @importFrom rlang call_modify !! :=
#'
runMiDAS_linear_omnibus <- function(call,
                                    midas,
                                    experiment,
                                    test_covar,
                                    correction = "bonferroni",
                                    n_correction = NULL,
                                    exponentiate = FALSE,
                                    ...) {
  omnibus_groups <- getOmnibusGroups(midas, experiment)
  assert_that(
    ! is.null(omnibus_groups),
    msg = sprintf("Omnibus test does not support experiment %s", experiment)
  )
  test_var <- rownames(midas[[experiment]])
  placeholder <- getPlaceholder(midas)

  # insert data for analysis
  data <- as.data.frame(midas)
  call <- call_modify(substitute(call), data = data)
  object <- eval(call)

  # run analysis
  results <- omnibusTest(object = object,
                         omnibus_groups = omnibus_groups,
                         placeholder = placeholder,
                         correction = correction,
                         n_correction = n_correction
  )

  assert_that(
    nrow(results) > 0,
    msg = "Could not process any variables. Please check warning messages for more informations (warnings())."
  )

  # format linear omnibus results
  group_name <- switch (experiment,
                       "hla_aa" = "aa_pos",
                       "group"
  )
  term_prefix <- switch (experiment,
                         "hla_aa" = "[A-Z0-9]+_-*[0-9]+_",
                         ""
  )
  term_name <- switch (experiment,
                       "hla_aa" = "residue",
                       "term"
  )
  results <- results %>%
    as_tibble() %>%
    mutate(term = gsub(term_prefix, "", .data$term)) %>%
    select(
      !! group_name := .data$group,
      !! term_name := .data$term,
      .data$dof,
      .data$statistic,
      .data$p.value,
      .data$p.adjusted
    ) %>%
    arrange(.data$p.value)

  return(results)
}

#' @rdname runMiDAS
#'
#' @title runMiDAS conditonal omnibus
#'
#' @details statistical analysis is performed iteratively on groups of variables
#'   like residues at particular amino acid position, using likelyhood ratio
#'   test. This is done by substituting \code{placeholder} in the
#'   \code{object}'s formula with linear combination of variable in particular
#'   group. In each iteration best variable from previous one is added to
#'   formula.
#'
#' @param object_details TODO
#'
#' @importFrom assertthat assert_that is.number
#' @importFrom dplyr arrange as_tibble mutate select
#' @importFrom magrittr %>%
#' @importFrom stats formula
#' @importFrom rlang call_modify !! :=
#'
runMiDAS_conditional_omnibus <- function(call,
                                    midas,
                                    experiment,
                                    test_covar,
                                    correction = "bonferroni",
                                    n_correction = NULL,
                                    keep = FALSE,
                                    th,
                                    rss_th = 1e-07,
                                    exponentiate = FALSE,
                                    ...) {
  omnibus_groups <- getOmnibusGroups(midas, experiment)
  assert_that(
    is.flag(keep),
    is.number(th),
    is.number(rss_th),
    ! is.null(omnibus_groups),
    msg = sprintf("Omnibus test does not support experiment %s", experiment)
  )
  placeholder <- getPlaceholder(midas)

  # insert data for analysis
  data <- as.data.frame(midas)
  call <- call_modify(substitute(call), data = data)
  object <- eval(call)

  # run analysis
  prev_formula <- formula(object)
  first_variables <- all.vars(prev_formula)
  prev_omnibus_groups <- c()
  new_omnibus_groups <-
    omnibus_groups[! names(omnibus_groups) %in% first_variables]

  best <- list()
  i <- 1
  while (length(new_omnibus_groups) > 0) {
    results <- omnibusTest(
      object = object,
      omnibus_groups = new_omnibus_groups,
      placeholder = placeholder,
      correction = correction,
      n_correction = n_correction
    )

    results <- results[! is.infinite(results[["p.value"]]), ]

    mask <- ! prev_omnibus_groups %in% first_variables
    results$covariates <-
      paste(prev_omnibus_groups[mask], collapse = " + ")

    i_min <- which.min(results[["p.value"]])
    if (length(i_min) == 0) break
    if (results$p.value[i_min] > th) break

    gr_name <- results$group[i_min]
    prev_omnibus_groups <- c(prev_omnibus_groups, gr_name)
    object <- updateModel(object,
                          omnibus_groups[[gr_name]],
                          backquote = TRUE,
                          collapse = " + "
    )

    if (sum(resid(object) ^ 2) <= rss_th) {
      warn("Perfect fit was reached attempting further model selection is nonsense.")
      break
    }
    prev_formula <- formula(object)
    prev_variables <- all.vars(prev_formula)
    new_omnibus_groups <-
      new_omnibus_groups[! names(new_omnibus_groups) %in% prev_omnibus_groups]

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

  # format linear omnibus results
  group_name <- switch (experiment,
                        "hla_aa" = "aa_pos",
                        "group"
  )
  term_prefix <- switch (experiment,
                         "hla_aa" = "[A-Z0-9]+_-*[0-9]+_",
                         ""
  )
  term_name <- switch (experiment,
                       "hla_aa" = "residue",
                       "term"
  )

  .formatResults <- function(df) {
    df %>%
      as_tibble() %>%
      mutate(term = gsub(term_prefix, "", .data$term)) %>%
      select(
        !! group_name := .data$group,
        !! term_name := .data$term,
        .data$dof,
        .data$statistic,
        .data$p.value,
        .data$p.adjusted
      ) %>%
      arrange(.data$p.value)
  }
  if (is.list(results)) {
    results <- lapply(results, .formatResults)
  } else {
    results <- .formatResults(results)
  }

  return(results)
}

#' Test for Hardy Weinberg Equilibrium
#'
#' Test HLA calls experiment (matrix) for Hardy Weinberg Equilibrium. This
#' function \code{HWChisqStats} from \code{\link{HardyWeinberg}} package.
#'
#' @inheritParams stats::p.adjust
#' @param X Matrix with genotype counts, with samples as columns and alleles as
#'   rows.
#'
#' @return A vector of p-values
#'
#' @importFrom assertthat assert_that is.string
#' @importFrom stats p.adjust
#' @importFrom methods validObject
#' @importFrom MultiAssayExperiment colData
#' @importFrom HardyWeinberg HWChisqStats
#'
HWETest <-
  function(object,
           experiment = c("hla_alleles", "hla_aa", "hla_g_groups", "hla_supertypes", "hla_NK_ligands"),
           HWE_group = NULL,
           HWE_cutoff = NULL,
           as.MiDAS = FALSE) {

    experiment_choice <- eval(formals()[["experiment"]])
    assert_that(
      validObject(object),
      is.string(experiment),
      stringMatches(experiment, experiment_choice),
      isNumberOrNULL(HWE_cutoff),
      isTRUEorFALSE(as.MiDAS)
    )
    HWE_group <- substitute(HWE_group)
    if (is.null(HWE_group)) {
      X <- list(p.value = object[[experiment]])
    } else {
      colData <-
        do.call(subset,
                list(
                  x = colData(object),
                  subset = HWE_group
                )
        )
      subset_ids <- colData$ID
      x <- object[[experiment]]
      X <- list(
        x[, colnames(x) %in% subset_ids],
        x[, ! colnames(x) %in% subset_ids]
      )
      names(X) <- c(deparse(HWE_group), paste0("not ", deparse(HWE_group)))
    }

    HWE.result <- data.frame(
      var = rownames(object[[experiment]]),
      stringsAsFactors = FALSE
    )
    for (i in 1:length(X)) {
      HWE.pvalue <- apply(
        X = X[[i]],
        MARGIN = 1,
        FUN = function(x) {
          c(
            AA = sum(x == 0, na.rm = TRUE),
            AB = sum(x == 1, na.rm = TRUE),
            BB = sum(x == 2, na.rm = TRUE)
          )
        }
      ) %>%
        t() %>%
        HWChisqStats(x.linked = FALSE, pvalues = TRUE)
      nm <- names(X)[[i]]
      HWE.result[[nm]] <- HWE.pvalue[HWE.result$var] # NAs will be inserted on missing
    }

    if (! is.null(HWE_cutoff)) {
      mask <- lapply(1:(ncol(HWE.result) - 1), function(i) {
        m <- HWE.result[, i + 1] < HWE_cutoff[i]
        m[is.na(m)] <- FALSE
        m
      })
      mask <- Reduce(f = `&`, x = mask)
      HWE.result <- HWE.result[mask, ]
    }

    if(as.MiDAS) {
      HWE.result <- filterByVariables(object, experiment, HWE.result$var)
    }

    return(HWE.result)
  }
