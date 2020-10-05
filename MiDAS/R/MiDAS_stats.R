#' Association analysis
#'
#' \code{analyzeAssociations} perform association analysis on a single variable
#' level using a statistical model of choice.
#'
#' \code{correction} specifies p-value adjustment method to use, common choice
#' is Benjamini & Hochberg (1995) (\code{"BH"}). Internally this is passed to
#' \link[stats]{p.adjust}.
#'
#' @inheritParams updateModel
#' @param variables Character vector specifying variables to use in association
#'   tests.
#' @param placeholder String specifying term in \code{object}'s formula which
#'   should be substituted with variables during analysis.
#' @param correction String specifying multiple testing correction method. See
#'   details for further information.
#' @param n_correction Integer specifying number of comparisons to consider
#'   during multiple testing correction calculations. For Bonferroni correction
#'   it is possible to specify a number lower than the number of comparisons
#'   being made. This is useful in cases when knowledge about the biology or
#'   redundance of alleles reduces the need for correction. For other methods it
#'   must be at least equal to the number of comparisons being made; only set
#'   this (to non-default) when you know what you are doing!
#' @param exponentiate Logical flag indicating whether or not to exponentiate
#'   the coefficient estimates. Internally this is passed to
#'   \code{\link[broom]{tidy}}. This is typical for logistic and multinomial
#'   regressions, but a bad idea if there is no log or logit link. Defaults to
#'   FALSE.
#'
#' @return Tibble containing combined results for all \code{variables}.
#'
#' @examples
#' midas <- prepareMiDAS(hla_calls = MiDAS_tut_HLA,
#'                       colData = MiDAS_tut_pheno,
#'                       experiment = "hla_alleles")
#'
#' # analyzeAssociations expects model data to be a data.frame
#' midas_data <- as.data.frame(midas)
#'
#' # define base model
#' object <- lm(disease ~ term, data = midas_data)
#'
#' # test for alleles associations
#' analyzeAssociations(object = object,
#'                     variables = c("B*14:02", "DRB1*11:01"))
#'
#' @importFrom assertthat assert_that see_if is.string
#' @importFrom broom tidy
#' @importFrom dplyr bind_rows
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
    nc >= length(results$p.value) || correction == "bonferroni",
    msg = sprintf("n_correction must be at least %i.", length(results$p.value))
  )
  results$p.adjusted <- adjustPValues(
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
#' adding the previous top-associated variable as covariate, until there are no
#' more significant variables based on a self-defined threshold.
#'
#' @inheritParams updateModel
#' @inheritParams analyzeAssociations
#' @param th Number specifying threshold for a variable to be considered
#'   significant.
#' @param th_adj Logical flag indicating if adjusted p-value should be used as
#'   threshold criteria, otherwise unadjusted p-value is used.
#' @param keep Logical flag indicating if the output should be a list of results
#'   resulting from each selection step. Default is to return only the final
#'   result.
#' @param rss_th Number specifying residual sum of squares threshold at which
#'   function should stop adding additional variables. As the residual sum of
#'   squares approaches \code{0} the perfect fit is obtained making further
#'   attempts at variable selection nonsense. This behavior can be controlled
#'   using \code{rss_th}.
#'
#' @return Tibble with stepwise conditional testing results.
#'
#' @examples
#' midas <- prepareMiDAS(hla_calls = MiDAS_tut_HLA,
#'                       colData = MiDAS_tut_pheno,
#'                       experiment = "hla_alleles")
#'
#' # analyzeConditionalAssociations expects model data to be a data.frame
#' midas_data <- as.data.frame(midas)
#'
#' # define base model
#' object <- lm(disease ~ term, data = midas_data)
#' analyzeConditionalAssociations(object,
#'                             variables = c("B*14:02", "DRB1*11:01"),
#'                             th = 0.05)
#'
#' @importFrom assertthat assert_that is.number is.string
#' @importFrom dplyr bind_rows tibble
#' @importFrom rlang warn
#' @importFrom stats formula resid
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

  # select optimization criteria
  crit <- ifelse(th_adj, "p.adjusted", "p.value")

  prev_formula <- object_formula
  first_variables <- all.vars(object_formula)
  prev_variables <- first_variables
  new_variables <- variables[! variables %in% prev_variables]

  best <- list()
  i <- 1
  while (length(new_variables) > 0) {
    # reorder object formula such that 'term' is at the end, otherwise TL;DR
    # when the new covariates are added to the model formula start to look like
    # this: resp ~ term + A + B + C
    # as a result if any of the new variables have the same effect
    # as any old one (say X = c(0, 0, 1, 1), Z = c(0, 0, 1, 1)),
    # it won't be corrected for and we will get multiple vars with same estimate
    # instead formula should look like resp ~ A + B + C + term
    new_form <- paste0(". ~ . - ", placeholder, " + ", placeholder)
    object <- update(object, new_form, evaluate = FALSE)
    object <- eval(object, envir = object_env)

    results <- iterativeModel(object, placeholder, new_variables, exponentiate)

    nc <- ifelse(is.null(n_correction), length(results$p.value), n_correction)
    assert_that(
      nc >= length(results$p.value) || correction == "bonferroni",
      msg = sprintf("n_correction must be at least %i.", length(results$p.value))
    )
    results$p.adjusted <- adjustPValues(
      p = results$p.value,
      method = correction,
      n = nc
    )

    results <- results[! is.infinite(results[["p.value"]]), ]

    mask <- ! prev_variables %in% first_variables
    results$covariates <- paste(prev_variables[mask], collapse = " + ")

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
#' \code{OmnibusTest} calculates overall p-value for linear combination of
#' variables using likelihood ratio test.
#'
#' Likelihood ratio test is conducted by comparing a model given in an
#' \code{object} with an extended model, that is created by including the effect
#' of variables given in \code{variables} as their linear combination.
#'
#' @inheritParams analyzeAssociations
#' @param omnibus_groups List of character vectors giving sets of variables for
#'  which omnibus test should be applied.
#'
#' @return Data frame containing omnibus test results for specified variables.
#'
#' @examples
#' midas <- prepareMiDAS(hla_calls = MiDAS_tut_HLA,
#'                       colData = MiDAS_tut_pheno,
#'                       experiment = "hla_aa")
#'
#' # define base model
#' object <- lm(disease ~ term, data = midas)
#' omnibusTest(object,
#'             omnibus_groups = list(
#'               A_29 = c("A_29_D", "A_29_A"),
#'               A_43 = c("A_43_Q", "A_43_R")
#'             ))
#'
#' @importFrom dplyr bind_cols
#' @export
omnibusTest <- function(object,
                        omnibus_groups,
                        placeholder = "term",
                        correction = "bonferroni",
                        n_correction = NULL) {
  results <- iterativeLRT(object, placeholder, omnibus_groups)

  nc <- ifelse(is.null(n_correction), length(results$p.value), n_correction)
  assert_that(
    nc >= length(results$p.value) || correction == "bonferroni",
    msg = sprintf("n_correction must be at least %i.", length(results$p.value))
  )
  results$p.adjusted <- adjustPValues(
    p = results$p.value,
    method = correction,
    n = nc
  )

  return(results)
}

#' Run MiDAS statistical analysis
#'
#' \code{runMiDAS} perform association analysis on MiDAS data using statistical
#' model of choice. Function is intended for use with \code{\link{prepareMiDAS}}.
#' See examples section.
#'
#' By default statistical analysis is performed iteratively on each variable in
#' selected experiment. This is done by substituting \code{placeholder} in the
#' \code{object}'s formula with each variable in the experiment.
#'
#' Setting \code{conditional} argument to \code{TRUE} will cause the statistical
#' analysis to be performed in a stepwise conditional testing manner, adding the
#' previous top-associated variable as a covariate to \code{object}'s formula.
#' The analysis stops when there is no more significant variables, based on
#' self-defined threshold (\code{th} argument). Either adjusted or unadjusted
#' p-values can be used as the selection criteria, which is controlled using
#' \code{th_adj} argument.
#'
#' Setting \code{omnibus} argument to \code{TRUE} will cause the statistical
#' analysis to be performed iteratively on groups of variables (like residues at
#' particular amino acid position) using likelihood ratio test.
#'
#' \code{correction} specifies p-value adjustment method to use, common choice
#' is Benjamini & Hochberg (1995) (\code{"BH"}). Internally this is passed to
#' \link[stats]{p.adjust}.
#'
#' @inheritParams analyzeAssociations
#' @inheritParams analyzeConditionalAssociations
#' @inheritParams filterByFrequency
#' @inheritParams applyInheritanceModel
#' @param experiment String indicating the experiment associated with
#'   \code{object}'s \code{MiDAS} data to use. Valid values includes:
#'   \code{"hla_alleles"}, \code{"hla_aa"}, \code{"hla_g_groups"},
#'   \code{"hla_supertypes"}, \code{"hla_NK_ligands"}, \code{"kir_genes"},
#'   \code{"kir_haplotypes"}, \code{"hla_kir_interactions"},
#'   \code{"hla_divergence"}, \code{"hla_het"}, \code{"hla_custom"},
#'   \code{"kir_custom"}. See \code{\link{prepareMiDAS}} for more information.
#' @param conditional Logical flag indicating if conditional analysis should be
#'   performed.
#' @param omnibus Logical flag indicating if omnibus test should be used.
#' @param omnibus_groups_filter Character vector specifying omnibus groups to
#'   use.
#'
#' @return Tibble containing analysis results.
#'
#' @examples
#' # create MiDAS object
#' midas <- prepareMiDAS(hla_calls = MiDAS_tut_HLA,
#'                       colData = MiDAS_tut_pheno,
#'                       experiment = c("hla_alleles", "hla_aa")
#' )
#'
#' # construct statistical model
#' object <- lm(disease ~ term, data = midas)
#'
#' # run analysis
#' runMiDAS(object, experiment = "hla_alleles", inheritance_model = "dominant")
#'
#' # omnibus test
#' # omnibus_groups_filter argument can be used to restrict omnibus test only
#' # to selected variables groups, here we restrict the analysis to HLA-A
#' # positions 29 and 43.
#' runMiDAS(
#'   object,
#'   experiment = "hla_aa",
#'   inheritance_model = "dominant",
#'   omnibus = TRUE,
#'   omnibus_groups_filter = c("A_29", "A_43")
#' )
#'
#' @importFrom assertthat assert_that is.number is.number is.string
#' @importFrom dplyr arrange as_tibble bind_rows mutate rename select
#' @importFrom magrittr %>%
#' @importFrom stats formula
#' @importFrom rlang call_modify warn !! :=
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
                     th = 0.05,
                     th_adj = TRUE,
                     keep = FALSE,
                     rss_th = 1e-07
                    ) {
  assert_that(
    checkStatisticalModel(object),
    hasTidyMethod(class(object)[1L])
  )
  object_details <- getObjectDetails(object)

  assert_that(
    see_if(
      isClass(object_details$data, "MiDAS"),
      msg = "data associated with statistical model must be an instance of MiDAS class."
    ),
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
    if (isExperimentInheritanceModelApplicable(object_details$data[[experiment]])) {
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
    th = th,
    th_adj = th_adj,
    keep = keep,
    rss_th = rss_th
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

#
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
    msg = "Could not process any variables. Please check warning messages for more information (warnings())."
  )

  # oreder results columns
  col_ord <- c("term", "p.value", "p.adjusted", "estimate", "std.error",
               "conf.low", "conf.high", "statistic")
  if (all(col_ord  %in% colnames(results))) {
    col_ord <- c(col_ord, colnames(results)[! colnames(results) %in% col_ord]) # append columns not specified in col_ord
    results <- results[, col_ord]
  }

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

#
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
    msg = "Could not process any variables. Please check warning messages for more information (warnings())."
  )

 # oreder results columns
 col_ord <- c("term", "p.value", "p.adjusted", "estimate", "std.error",
              "conf.low", "conf.high", "statistic", "covariates")
 if (all(col_ord %in% colnames(results))) {
   col_ord <- c(col_ord, colnames(results)[! colnames(results) %in% col_ord]) # append columns not specified in col_ord
   results <- results[, col_ord]
 }

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
            arrange(.data$p.value)
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

#
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
    msg = "Could not process any variables. Please check warning messages for more information (warnings())."
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
                       "hla_aa" = "residues",
                       "term"
  )
  results <- results %>%
    as_tibble() %>%
    mutate(term = gsub(term_prefix, "", .data$term)) %>%
    select(
      !! group_name := .data$group,
      !! term_name := .data$term,
      .data$df,
      .data$statistic,
      .data$p.value,
      .data$p.adjusted
    ) %>%
    arrange(.data$p.value)

  return(results)
}

#
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
    # reorder object formula such that 'term' is at the end, otherwise TL;DR
    # when the new covariates are added to the model formula start to look like
    # this: resp ~ term + A + B + C
    # as a result if any of the new variables have the same effect
    # as any old one (say X = c(0, 0, 1, 1), Z = c(0, 0, 1, 1)),
    # it won't be corrected for and we will get multiple vars with same estimate
    # instead formula should look like resp ~ A + B + C + term
    new_form <- paste0(". ~ . - ", placeholder, " + ", placeholder)
    object <- update(object, new_form)

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
      results <- results[0,]
    } else {
      results <- lapply(best, function(res) {
        i_min <- which.min(res[["p.value"]])
        res[i_min,]
      })
      results <- bind_rows(results)
    }
  }

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
                       "hla_aa" = "residues",
                       "term"
  )

  .formatResults <- function(df) {
    df %>%
      as_tibble() %>%
      mutate(term = gsub(term_prefix, "", .data$term)) %>%
      select(
        !! group_name := .data$group,
        !! term_name := .data$term,
        .data$df,
        .data$statistic,
        .data$p.value,
        .data$p.adjusted,
        .data$covariates
      ) %>%
      arrange(.data$p.value)
  }
  if (keep) {
    results <- lapply(results, .formatResults)
  } else {
    results <- .formatResults(results)
  }

  return(results)
}

#' Test for Hardy Weinberg equilibrium
#'
#' Test experiment features for Hardy Weinberg equilibrium.
#'
#' Setting \code{as.MiDAS} to \code{TRUE} will filter MiDAS object based on
#' p-value cut-off given by \code{HWE_cutoff}.
#'
#' @param object \code{\link{MiDAS}} object.
#' @param experiment String specifying experiment to test. Valid values
#'   includes \code{"hla_alleles"}, \code{"hla_aa"}, \code{"hla_g_groups"},
#'   \code{"hla_supertypes"}, \code{"hla_NK_ligands"}.
#' @param HWE_group Expression defining samples grouping to test for Hardy
#'   Weinberg equilibrium. By default samples are not grouped.
#' @param HWE_cutoff Number specifying p-value threshold. When \code{HWE_group}
#'   is specified both groups are thresholded.
#' @param as.MiDAS Logical flag indicating if MiDAS object should be returned.
#'
#' @return Data frame with Hardy Weinberg Equilibrium test results or a filtered
#'   MiDAS object.
#'
#' @examples
#' # create MiDAS object
#' midas <- prepareMiDAS(hla_calls = MiDAS_tut_HLA,
#'                       colData = MiDAS_tut_pheno,
#'                       experiment = "hla_alleles"
#' )
#'
#' # get HWE p-values as data frame
#' HWETest(midas, experiment = "hla_alleles")
#'
#' # get HWE in groups defined by disease status
#' # grouping by `disease == 1` will divide samples into two groups:
#' # `disease == 1` and `not disease == 1`
#' HWETest(midas, experiment = "hla_alleles", HWE_group = disease == 1)
#'
#' # filter MiDAS object by HWE test p-value
#' HWETest(midas, experiment = "hla_alleles", HWE_cutoff = 0.05, as.MiDAS = TRUE)
#'
#' @importFrom HardyWeinberg HWChisqStats
#' @importFrom assertthat assert_that is.string
#' @importFrom methods validObject
#' @importFrom SummarizedExperiment assay
#' @importFrom MultiAssayExperiment colData
#' @export
HWETest <-
  function(object,
           experiment = c("hla_alleles", "hla_aa", "hla_g_groups", "hla_supertypes", "hla_NK_ligands"),
           HWE_group = NULL,
           HWE_cutoff = NULL,
           as.MiDAS = FALSE) {
    experiment_choice <- eval(formals()[["experiment"]])
    assert_that(
      isClass(object, "MiDAS"),
      validObject(object),
      is.string(experiment),
      stringMatches(experiment, experiment_choice),
      isNumberOrNULL(HWE_cutoff),
      isTRUEorFALSE(as.MiDAS)
    )

    if (is(object[[experiment]], "SummarizedExperiment")) {
      x <- assay(object[[experiment]])
    } else {
      x <- object[[experiment]]
    }
    HWE_group <- substitute(HWE_group)
    if (is.null(HWE_group)) {
      X <- list(p.value = x)
    } else {
      colData <-
        do.call(subset,
                list(
                  x = colData(object),
                  subset = HWE_group
                )
        )
      subset_ids <- rownames(colData)
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
      HWE.pvalue <- apply( # get haplotypes for each variable
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

    # filter based on p-value cut-off
    if (! is.null(HWE_cutoff)) {
      # get p-value mask over columns ie. p.value or optional grouping columns
      mask <- lapply(1:(ncol(HWE.result) - 1), function(i) {
        m <- HWE.result[, i + 1] > HWE_cutoff
        m[is.na(m)] <- TRUE
        m
      })
      mask <- Reduce(f = `&`, x = mask)
      HWE.result <- HWE.result[mask, ]
    }

    if(as.MiDAS) {
      vars <- rownames(object[[experiment]])
      vars <- vars[vars %in% HWE.result$var]
      HWE.result <- filterByVariables(object, experiment, vars)
    }

    return(HWE.result)
  }
