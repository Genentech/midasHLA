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
#'                     variables = c("B*14:02", "DRB1*11:01")
#' )
#'
#' @importFrom assertthat assert_that see_if is.flag is.string
#' @importFrom broom tidy
#' @importFrom dplyr bind_rows
#' @importFrom stats p.adjust
#'
#' @export
analyzeAssociations <- function(object,
                                variables,
                                correction = "BH",
                                exponentiate = FALSE) {
  assert_that(
    checkStatisticalModel(object)
  )
  object_call <- getCall(object)
  object_data <- eval(object_call[["data"]], envir = parent.frame())
  object_variables <- colnames(object_data)[-1]

  assert_that(
    is.character(variables),
    see_if(
      all(test_vars <- variables %in% object_variables) | is.null(variables),
      msg = sprintf("%s can not be found in object data",
                    paste(variables[! test_vars], collapse = ", ")
      )
    ),
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
  results$term <- gsub("`", "", results$term)
  results <- results[results$term %in% variables, ]

  results$p.adjusted <- p.adjust(results$p.value, correction)

#  This covariates were added for consistiency with conditional analyze, now however that we are filtering covariates there it doesn't make much sense to keep those?
#  covariates <- formula(object)[[3]]
#  covariates <- deparse(covariates)
#  results$covariates <- covariates

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
#' @return tibble with stepwise conditional testing results.
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
#' @importFrom assertthat assert_that is.flag is.number is.string
#' @importFrom dplyr bind_rows tibble
#' @importFrom purrr map_dfr
#' @importFrom rlang warn
#' @importFrom stats formula resid
#'
#' @export
analyzeConditionalAssociations <- function(object,
                                           variables,
                                           correction = "BH",
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
    is.character(variables),
    see_if(
      all(test_vars <- variables %in% object_variables),
      msg = sprintf("%s can not be found in object data",
                    paste(variables[! test_vars], collapse = ", ")
      )
    ),
    is.string(correction),
    is.number(th),
    is.number(rss_th),
    is.flag(exponentiate)
  )

  prev_formula <- object_formula
  first_variables <- all.vars(object_formula)
  prev_variables <- first_variables
  new_variables <- variables[! variables %in% prev_variables]

  best <- list()
  i <- 1

  while (length(new_variables) > 0) {
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
    results$p.adjusted <- p.adjust(results[["p.value"]], correction)
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

    results <- results[i_min, ]
    results$term <- gsub("`", "", results$term)
    best[[i]] <- results
    i <- i + 1
  }

  if (length(best) == 0) {
    warn("No significant variables found. Returning empty table.") # Tibble to be more precise?
    results <- results[0, ]
  } else {
    results <- bind_rows(best)
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

#' Analyze associations in MiDAS data
#'
#' \code{analyzeMiDASData} performs association analysis on MiDAS data using
#' statistical model specified by user.
#'
#' @inheritParams analyzeAssociations
#' @inheritParams analyzeConditionalAssociations
#' @inheritParams formatAssociationsResults
#' @param analysis_type String indicating the type of analysis being performed,
#'   at this point it is only used for results formatting. Valid values are
#'   \code{"hla_alleles"}, \code{"aa_level"}, \code{"expression_levels"},
#'   \code{"allele_groups"}, \code{"custom"}.
#' @param conditional Logical indicating if the analysis should be performed
#'   using stepwise conditional tests or not. See
#'   \link{analyzeConditionalAssociations} for more details.
#' @param variables Character specifying variables to use in association tests
#'   or \code{NULL}. If \code{NULL} all variables in object data are tested.
#'   See details for further information.
#' @param frequency_cutoff Number specifying threshold for inclusion of a
#'   variable. If it's a number between 0 and 1 variables with frequency below
#'   this number will not be considered during analysis. If it's greater or
#'   equal 1 variables with number of counts less that this will not be
#'   considered during analysis.
#' @param logistic Logical indicating if statistical model used is logistic (eg.
#'   \code{coxph}). If \code{NULL} function will try to figure this out. This is
#'   only used for results formatting.
#' @param binary_phenotype Logical indicating if coefficient estimates should be
#'   exponentiated. This is typical for logistic and multinomial regressions,
#'   but a bad idea if there is no log or logit link. If \code{NULL} function
#'   will try to figure this out by testing if response is binary (\code{0} or
#'   \code{1}).
#' @param kable_output Logical indicating if additionally results should be
#'   pretty printed in specified \code{format}.
#'
#' \code{variables} takes \code{NULL} as a default value, than column names of
#' data frame associated with the \code{object} are used as variables in
#' association tests. This excludes first column which should corresponds
#' to samples IDs as well as covariates and response variables defined in
#' object formula.
#'
#' \code{correction} specifies p-value adjustment method to use, common choice
#' is Benjamini & Hochberg (1995) (\code{"BH"}). Internally this is passed to
#' \link[stats]{p.adjust}. Check there to get more details.
#'
#' @return Tibble containing results for all tested variables.
#'
#' @examples
#' library("survival")
#' hla_calls_file <- system.file("extdata", "HLAHD_output_example.txt", package = "MiDAS")
#' hla_calls <- readHlaCalls(hla_calls_file)
#' pheno_file <- system.file("extdata", "pheno_example.txt", package = "MiDAS")
#' pheno <- read.table(pheno_file, header = TRUE, stringsAsFactors = FALSE)
#' covar_file <- system.file("extdata", "covar_example.txt", package = "MiDAS")
#' covar <- read.table(covar_file, header = TRUE, stringsAsFactors = FALSE)
#' midas_data <- prepareMiDASData(hla_calls = hla_calls,
#'                                pheno = pheno,
#'                                covar = covar,
#'                                analysis_type = "hla_allele",
#'                                inheritance_model = "additive"
#' )
#'
#' object <- coxph(Surv(OS, OS_DIED) ~ AGE + SEX, data = midas_data)
#' analyzeMiDASData(object, analysis_type = "hla_alleles")
#'
#' @importFrom assertthat assert_that is.flag is.number is.string
#' @importFrom dplyr filter left_join select rename
#' @importFrom stats getCall
#' @importFrom rlang !! := .data
#' @importFrom magrittr %>% %<>%
#'
#' @export
analyzeMiDASData <- function(object,
                             analysis_type = c("hla_alleles", "aa_level", "expression_levels",  "allele_groups", "custom"),
                             conditional = FALSE,
                             variables = NULL,
                             frequency_cutoff = NULL,
                             pvalue_cutoff = NULL,
                             correction = "BH",
                             logistic = NULL,
                             binary_phenotype = NULL,
                             th = 0.05,
                             rss_th = 1e-07,
                             kable_output = TRUE,
                             format = getOption("knitr.table.format")) {

  assert_that(
    checkStatisticalModel(object)
  )
  object_call <- getCall(object)
  object_formula <- eval(object_call[["formula"]], envir = parent.frame())
  object_data <- eval(object_call[["data"]], envir = parent.frame())
  object_variables <- colnames(object_data)[-1]

  assert_that(
    is.string(analysis_type),
    stringMatches(analysis_type,
                  choice = c("hla_alleles", "aa_level", "expression_levels",  "allele_groups", "custom")
    ),
    is.flag(conditional),
    isCharacterOrNULL(variables),
    see_if(
      all(test_vars <- variables %in% object_variables) | is.null(variables),
      msg = sprintf("%s can not be found in object data",
                    paste(variables[! test_vars], collapse = ", ")
      )
    ),
    isNumberOrNULL(frequency_cutoff),
    isNumberOrNULL(pvalue_cutoff),
    is.string(correction),
    isFlagOrNULL(logistic),
    isFlagOrNULL(binary_phenotype),
    is.number(th),
    is.number(rss_th),
    is.flag(kable_output),
    is.string(format),
    stringMatches(format, choice = c("html", "latex"))
  )

  if (is.null(variables)) {
    mask <- ! object_variables %in% all.vars(object_formula)
    variables <- object_variables[mask]
    assert_that(length(variables) != 0,
                msg = "No new variables found in object data."
    )
  }

  # guess if model used is logistic type
  if (is.null(logistic)) {
    model_fun <- deparse(object_call[[1]])
    model_family <- deparse(object_call[["family"]])
    logistic <- grepl("coxph", model_fun) |
      (grepl("glm", model_fun) & grepl("binomial", model_family))
  }

  # Filter variables on frequency cutoff
  frequency_cutoff <- ifelse(is.null(frequency_cutoff), 0, frequency_cutoff)
  variables_freq <- object_data %>%
    select("ID",!! variables) %>%
    getCountsFrequencies() %>%
    rename(Ntotal = .data$Counts, Ntotal.frequency = .data$Freq) %>%
    filter(.data$Ntotal > frequency_cutoff |
             frequency_cutoff < 1) %>%
    filter(.data$Ntotal.frequency > frequency_cutoff |
             frequency_cutoff >= 1)

  if (conditional) {
    results <- analyzeConditionalAssociations(object,
                                              variables = variables_freq$term,
                                              correction = correction,
                                              th = th,
                                              rss_th = rss_th,
                                              exponentiate = logistic
    )
  } else {
    results <- analyzeAssociations(object,
                                   variables = variables_freq$term,
                                   correction = correction,
                                   exponentiate = logistic
    )
  }

  # Add frequency information to results table
  results <- left_join(x = results, y = variables_freq, by = "term")

  pheno_var <- all.vars(object_formula)[1]
  if (is.null(binary_phenotype)) {
    binary_phenotype <- object_data[, pheno_var] %in% c(0, 1)
    binary_phenotype <- all(binary_phenotype, na.rm = TRUE)
  }

  if (binary_phenotype) {
    results <- object_data %>%
      filter(.data[[!! pheno_var]] == 1) %>%
      select("ID", !! variables) %>%
      getCountsFrequencies() %>%
      rename(Npositive = .data$Counts,
             Npositive.frequency = .data$Freq) %>%
      left_join(x = results, by = "term")

    results <- object_data %>%
      filter(.data[[!! pheno_var]] != 1) %>%
      select("ID", !! variables) %>%
      getCountsFrequencies() %>%
      rename(Nnegative = .data$Counts,
             Nnegative.frequency = .data$Freq) %>%
      left_join(x = results, by = "term")
  }

  if (kable_output) {
    response_variable <- deparse(object_formula[[2]])

    preety_table <- formatAssociationsResults(
      results,
      type = analysis_type,
      response_variable = response_variable,
      logistic = logistic,
      pvalue_cutoff = pvalue_cutoff,
      format = format
    )
    print(preety_table)
  }

  # rename term and estimate to match preety_table
  estimate_name <- ifelse(logistic, "odds.ratio", "estimate")
  term_name <- switch (analysis_type,
                       "hla_alleles" = "allele",
                       "aa_level" = "aa",
                       "expression_levels" = "allele",
                       "allele_groups" = "allele.group",
                       "term"
  )

  results %<>%
    rename(!! term_name := .data$term, !! estimate_name := .data$estimate)

  return(results)
}

#' Prepare HLA calls data for statistical analysis
#'
#' \code{prepareMiDASData} transform HLA alleles calls according to selected
#' analysis type and joins obtained transformation with additional data frames
#' like phenotypic observations or covariates, creating an input data for
#' further statistical analysis.
#'
#' @inheritParams checkHlaCallsFormat
#' @inheritParams hlaCallsToCounts
#' @inheritParams hlaToAAVariation
#' @param ... Data frames holding additional variables like phenotypic
#'   observations or covariates.
#' @param analysis_type String indicating analysis type for which data should be
#'   prepared. Valid choices are \code{"hla_allele"}, \code{"aa_level"},
#'   \code{"expression_levels"}, \code{"allele_groups"}, \code{"custom"}. See
#'   details for further explanations.
#'
#' \code{...} should be data frames with first column holding samples IDs and
#' named \code{ID}. Those should correspond to \code{ID} column in
#' \code{hla_calls}.
#'
#' \code{analysis_type} specifies type of analysis for which \code{hla_calls}
#' should be prepared:
#'
#' \code{"hla_allele"} - \code{hla_calls} are transformed into counts  under
#' \code{inheritance_model} of choice (this is done with
#' \link{hlaCallsToCounts}).
#'
#' \code{"aa_level"} - \code{hla_calls} are first converted to amino acid level,
#' taking only variable positions under consideration. Than variable amino acid
#' positions are transformed to counts under \code{inheritance_model} of choice
#' (this is done with \link{aaVariationToCounts}).
#'
#' \code{"expression_levels"} - \code{hla_calls} are transformed to expression
#' levels using expression dictionaries shipped with package (this is done using
#' \link{hlaToVariable}). The expression levels from both alleles are than
#' summed into single variable for each translated HLA gene.
#'
#' \code{"allele_g_group"} - \code{hla_calls} are transformed to HLA alleles
#' groups using G group dictionary shipped with package (this is done using
#' \link{hlaToVariable}). Than those are transformed to counts under
#' \code{inheritance_model} of choice (this is done with
#' \link{hlaCallsToCounts}).
#'
#' \code{"allele_supertypes"} - \code{hla_calls} are transformed to HLA alleles
#' groups using supertypes dictionary shipped with package (this is done using
#' \link{hlaToVariable}). Than those are transformed to counts under
#' \code{inheritance_model} of choice (this is done with
#' \link{hlaCallsToCounts}).
#'
#' \code{"allele_groups"} - \code{hla_calls} are transformed to HLA alleles
#' groups using Bw4/6 and C1/2 groups dictionaries shipped with package (this is
#' done using \link{hlaToVariable}). Than those are transformed to counts under
#' \code{inheritance_model} of choice (this is done with
#' \link{hlaCallsToCounts}).
#'
#' \code{"custom"} - will not transform \code{hla_calls} and only joins it with
#' additional data(\code{...}).
#'
#' @return Data frame containing prepared data.
#'
#' @examples
#' hla_calls_file <- system.file("extdata", "HLAHD_output_example.txt", package = "MiDAS")
#' hla_calls <- readHlaCalls(hla_calls_file)
#' pheno_file <- system.file("extdata", "pheno_example.txt", package = "MiDAS")
#' pheno <- read.table(pheno_file, header = TRUE)
#' covar_file <- system.file("extdata", "covar_example.txt", package = "MiDAS")
#' covar <- read.table(covar_file, header = TRUE)
#' prepareMiDASData(hla_calls, pheno, covar, analysis_type = "expression_levels")
#'
#' @importFrom assertthat assert_that is.flag is.string see_if
#' @importFrom dplyr funs group_by left_join mutate summarise_all syms
#' @importFrom magrittr %>% %<>%
#' @importFrom rlang .data !!!
#' @importFrom tidyr gather spread
#' @importFrom Hmisc label label<-
#'
#' @export
prepareMiDASData <- function(hla_calls,
                             ...,
                             analysis_type = c("hla_allele", "aa_level", "expression_levels", "allele_g_group", "allele_supertypes", "allele_groups", "custom"),
                             inheritance_model = "additive",
                             indels = TRUE,
                             unkchar = FALSE
) {
  assert_that(
    checkHlaCallsFormat(hla_calls),
    is.string(analysis_type),
    stringMatches(
      x = analysis_type,
      choice = c("hla_allele", "aa_level", "expression_levels", "allele_g_group", "allele_supertypes", "allele_groups", "custom")
    ),
    is.string(inheritance_model),
    stringMatches(
      x = inheritance_model,
      choice = c("dominant", "recessive", "additive")
    ),
    is.flag(indels),
    is.flag(unkchar)
  )

  additional_data <- list(...)
  for (additional_data_frame in additional_data) { # rewrite those tests to addhere to current standard
    assert_that(
      checkAdditionalData(additional_data_frame, hla_calls, accept.null = FALSE)
    )
  }

  # Process hla_calls based on analysis type
  if (analysis_type == "hla_allele") {
    midas_data <- hlaCallsToCounts(
      hla_calls = hla_calls,
      inheritance_model = inheritance_model
    )

    label(midas_data[-1], self = FALSE) <- rep("hla_allele", ncol(midas_data) - 1)

  } else if (analysis_type == "aa_level") {
    midas_data <- hlaToAAVariation(
      hla_calls = hla_calls,
      indels = indels,
      unkchar = unkchar
    ) %>%
      aaVariationToCounts(inheritance_model = inheritance_model)

    label(midas_data[-1], self = FALSE) <- rep("aa_level", ncol(midas_data) - 1)

  } else if (analysis_type == "expression_levels") {
    lib <- listMiDASDictionaries()
    lib <- grep("expression", lib, value = TRUE)
    midas_data <- Reduce(
      f = function(...) left_join(..., by = "ID"),
      x = lapply(lib, hlaToVariable, hla_calls = hla_calls)
    )

    assert_that(
      ncol(midas_data) > 1,
      msg = "no expression levels were found for input hla_calls"
    )

    midas_data %<>%
      gather("expression", "value", -c("ID")) %>%
      mutate(expression = gsub("_.*", "", .data$expression)) %>%
      group_by(!!! syms(c("ID", "expression"))) %>%
      summarise_all(funs(sum)) %>%
      spread(.data$expression, .data$value, sep = "_")

    label(midas_data[-1], self = FALSE) <- rep(
      x = "expression_levels",
      ncol(midas_data) - 1
    )

  } else if (analysis_type == "allele_g_group") {
    lib <- "4digit_allele_Ggroup"
    midas_data <- hlaToVariable(hla_calls = hla_calls, dictionary = lib)

    assert_that(
      ncol(midas_data) > 1,
      msg = "no allele could be assigned to G group for input hla_calls"
    )

    midas_data <- hlaCallsToCounts(
      hla_calls = midas_data,
      inheritance_model = inheritance_model
    )

    label(midas_data[-1], self = FALSE) <- rep(
      x = "allele_g_group",
      ncol(midas_data) - 1
    )

  } else if (analysis_type == "allele_supertypes") {
    lib <- "4digit_supertype"
    midas_data <- hlaToVariable(hla_calls = hla_calls, dictionary = lib)


    assert_that(
      ncol(midas_data) > 1,
      msg = "no allele could be assigned to supertype for input hla_calls"
    )

    midas_data <- hlaCallsToCounts(
      hla_calls = midas_data,
      inheritance_model = inheritance_model,
      check_hla_format = FALSE
    ) %>%
      select(-"Unclassified")

    label(midas_data[-1], self = FALSE) <- rep(
      x = "allele_supertypes",
      ncol(midas_data) - 1
    )

  } else if (analysis_type == "allele_groups") {
    lib <- c("4digit_B-allele_Bw", "4digit_C-allele_C1-2")
    midas_data <- Reduce(
      f = function(...) left_join(..., by = "ID"),
      x = lapply(lib, hlaToVariable, hla_calls = hla_calls)
    )

    assert_that(
      ncol(midas_data) > 1,
      msg = "no allele could be assigned to allele groups for input hla_calls"
    )

    midas_data <- hlaCallsToCounts(
      hla_calls = midas_data,
      inheritance_model = inheritance_model,
      check_hla_format = FALSE
    )

    label(midas_data[-1], self = FALSE) <- rep(
      x = "allele_groups",
      ncol(midas_data) - 1
    )

  } else {
    midas_data <- hla_calls
  }

  # join with additional_data
  if (length(additional_data) != 0) {
    midas_data <- Reduce(
      f = function(...) left_join(..., by = "ID"),
      x = additional_data,
      init = midas_data
    )
  }

  return(midas_data)
}
