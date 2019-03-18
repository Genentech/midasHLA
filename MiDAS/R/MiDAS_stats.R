#' Association analysis of HLA allele calls
#'
#' \code{analyzeHlaAssociations} performs associations analysis on HLA alleles
#' calls using statistical model of choice.
#'
#' @inheritParams checkHlaCallsFormat
#' @param model String specifying statistical model to use.
#' @param pheno Data frame holding phenotypic response variables.
#' @param covar Data frame holding covariates.
#' @param zygo Flag indicating whether zygosity should be added to
#'   covariates. See details for further explanations.
#' @param reduce_counts Flag indicating whether allele counts should be reduced
#'   to presence / absence indicators.
#' @param correction String specifying multiple testing correction method.
#'
#' Available choices for \code{model} include:
#' \code{"coxph"} - Cox survival analysis,
#' \code{"lm"} - Linear regression,
#' \code{"glm.logit"} - Logistic regression,
#' \code{"glm.nb"} - Negative binomial regression
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
#' # Logistic regression
#' pheno <- pheno[, c(1, 3)]
#' analyzeHlaAssociations(model = "glm.logit",
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
#' @export
analyzeHlaAssociations <- function(model = "coxph",
                                   hla_calls,
                                   pheno,
                                   covar,
                                   zygo = FALSE,
                                   reduce_counts = FALSE,
                                   correction = "BH") {
  assert_that(
    is.string(model),
    see_if(model %in% hlaAssocModels(),
           msg = sprintf("model %s is not implemented", model)),
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
  covariate <- backquote(hla_data$covariate)
  model_function <- hlaAssocModels(model = model,
                                   response = response,
                                   covariate = covariate,
                                   data = hla_data$data
  )

  results <- map_dfr(
    .x = alleles,
    .f = ~tidy(model_function(.), exponentiate = FALSE) # TODO exponentiate could be passed somehow
  )
  results <- mutate(results, term = gsub("`", "", term))
  results <- filter(results, checkAlleleFormat(term))
  results <- rename(results, allele = term)

  results <- mutate(results, p.adjusted = p.adjust(p.value, correction))

  return(results)
}

#' Association models for analysis of HLA alleles
#'
#' \code{hlaAssocModels} is a collection of pre-configured models for use with
#' HLA alleles count table.
#'
#' \code{hlaAssocModels} is not intended to use by basic user.
#'
#' @inheritParams analyzeHlaAssociations
#' @param response Character specifying response variables in \code{data}.
#' @param covariate Character specifying covariates in \code{data}.
#' @param data Data frame containing variables in the model.
#'
#' @return Function for fitting association model of choice, it takes allele
#'   number as an argument. Additional arguments can be passed as well.
#'
#'   Returned function takes one required argument: HLA allele number
#'   (\code{allele}). Additional parameters are passed to statistical model
#'   function. Due to fact that allele numbers contains characters that have
#'   special meanings in formulas, they should be backquoted. This can be
#'   easily done with \link{backquote}. See examples section for general usage
#'   case.
#'
#'   If model is equal to \code{NULL} names of available models are returned
#'   instead.
#'
#' @examples
#' hla_calls_file <- system.file("extdata", "HLAHD_output_example.txt", package = "MiDAS")
#' hla_calls <- readHlaCalls(hla_calls_file)
#' hla_counts <- hlaCallsToCounts(hla_calls)
#' pheno_file <- system.file("extdata", "pheno_example.txt", package = "MiDAS")
#' pheno <- read.table(pheno_file, header = TRUE)
#' covar_file <- system.file("extdata", "covar_example.txt", package = "MiDAS")
#' covar <- read.table(covar_file, header = TRUE)
#' data <- dplyr::left_join(hla_counts, pheno, by = "ID")
#' data <- dplyr::left_join(data, covar, by = "ID")
#' response <- paste(colnames(pheno[, -1]), collapse = ", ")
#' covariate <- paste(colnames(covar[, -1]), collapse = " + ")
#' fun <- hlaAssocModels(model = "coxph",
#'                        response = response,
#'                        covariate = covariate,
#'                        data = data
#' )
#' allele <- backquote("A*01:01")
#' fun(allele)
#'
#' @importFrom assertthat assert_that see_if
#' @importFrom MASS glm.nb
#' @importFrom stats as.formula binomial glm lm
#' @importFrom survival coxph Surv
#' @export
hlaAssocModels <- function(model = NULL,
                           response,
                           covariate,
                           data) {
  if (is.null(model)) {
    return(c("coxph", "lm", "glm.logit", "glm.nb"))
  }
  assert_that(
    is.string(model),
    is.character(response),
    see_if(length(response) != 0, msg = "response can not be empty"),
    see_if(is.character(covariate) | is.null(covariate),
           msg = "covariate have to be a character or NULL"
    ),
    is.data.frame(data)
  )

  response_var <- paste(response, collapse = ", ")
  if (length(covariate) != 0) {
    covariate_var <- paste(covariate, collapse = " + ")
  } else {
    covariate_var <- "0"
  }

  model_function <- switch(
    model,
    coxph = function(allele, ...) { # Cox survival analysis
      form <- sprintf("Surv(%s) ~ %s + %s", response_var, allele, covariate_var)
      result <- coxph(formula = as.formula(form), data = data, ...)
      return(result)
    },
    lm = function(allele, ...) { # Linear regression
      form <- sprintf("%s ~ %s + %s", response_var, allele, covariate_var)
      result <- lm(formula = as.formula(form), data = data, ...)
      return(result)
    },
    glm.logit = function(allele, ...) { # Logistic regression
      form <- sprintf("%s ~ %s + %s", response_var, allele, covariate_var)
      result <- glm(
        formula = as.formula(form),
        family = binomial(link = "logit"),
        data = data,
        ...
      )
      return(result)
    },
    glm.nb = function(allele, ...) { # Negative binomial regression
      form <- sprintf("%s ~ %s + %s", response_var, allele, covariate_var)
      result <- glm.nb(formula = as.formula(form), data = data, ...)
      return(result)
    }
  )
  return(model_function)
}

#' Stepwise forward alleles subset selection
#'
#' \code{forwardAllelesSelection} does stepwise conditional testing adding the
#' previous top-associated allele as covariate, until there’s no more
#' significant alleles using a self-defined threshold.
#'
#' @param object object fitted by some model-fitting function.
#' @param scope formula specifying a maximal model which should include the
#'   current one. All additional terms in the maximal model with all marginal
#'   terms in the original model are tried.
#' @param th number specifying p-value threshold for a term to be included into
#'   model.
#' @param test string indicating test statistic to use for p-value calculation.
#'   Can be either "F" or "Chisq".
#' @param rss_th number specifying residual sum of squares threshold at which
#'   function should stop adding additional terms.
#'
#' All the variables in the \code{scope} should be defined in the \code{object}.
#'
#' The F test is only appropriate for lm and aov models, and perhaps for some
#' over-dispersed glm models. The Chisq test can be an exact test (lm models
#' with known scale) or a likelihood-ratio test depending on the method.
#'
#' As the residual sum of squares approaches \code{0} the perfect fit is
#' obtained making further attempts at model selection nonsense, thus function
#' is stopped. This behavior can be controlled using \code{rss_th}.
#'
#' @return selected model of the same class as \code{object} is returned.
#'
#' @examples
#' hla_calls_file <- system.file("extdata", "HLAHD_output_example.txt", package = "MiDAS")
#' hla_calls <- readHlaCalls(hla_calls_file)
#' pheno_file <- system.file("extdata", "pheno_example.txt", package = "MiDAS")
#' pheno <- read.table(pheno_file, header = TRUE)
#' covar_file <- system.file("extdata", "covar_example.txt", package = "MiDAS")
#' covar <- read.table(covar_file, header = TRUE)
#' hla_data <- prepareHlaData(hla_calls, pheno, covar)
#' object <- survival::coxph(survival::Surv(OS, OS_DIED) ~ AGE + SEX, data = hla_data$data)
#' scope <- survival::Surv(OS, OS_DIED) ~ AGE + SEX + `B*57:01` + `C*07:02`
#' forwardAllelesSelection(object, scope, th = 0.05, test = "Chisq")
#'
#' @importFrom assertthat assert_that is.number is.string see_if
#' @importFrom MASS addterm
#' @importFrom stringi stri_startswith_fixed
#' @importFrom purrr is_formula
#' @export
forwardAllelesSelection <- function(object,
                                    scope,
                                    th,
                                    test = c("F", "Chisq"),
                                    rss_th = 1e-07,
                                    ...) {
  assert_that(
    see_if("formula" %in% attr(object, "names"),
           msg = "object have to be a model"
    ),
    see_if(is_formula(scope), msg = "scope is not a formula"),
    is.number(th),
    is.string(test),
    is.number(rss_th)
  )

  test <- match.arg(test)

  criterium_pattern <- "Pr("
  scope_vars <- all.vars(scope)
  cur_vars <- all.vars(object$call$formula)

  while (! all(scope_vars %in% cur_vars)) {
    temp <- addterm(object, scope, test = test, ...)
    criterium_column <- stri_startswith_fixed(colnames(temp), criterium_pattern)
    assert_that(
      see_if(sum(criterium_column) == 1,
             msg = "Test criterium column couldn't be found..."
      )
    )
    test_val <- temp[-1, criterium_column]
    i_min <- which.min(test_val)
    if (length(i_min) == 0) break()
    if (test_val[i_min] > th) break()
    new_var <- rownames(temp)[i_min + 1]
    new_formula <- formula(object)
    new_formula <- as.formula(
      paste(new_formula[2], "~", paste(new_formula[3], new_var, sep = " + "))
    )
    object <- update(object, new_formula)
    if (sum(resid(object) ^ 2) <= rss_th) {
      warning("Perfect fit was reached attempting further model selection is nonsense.")
      break()
    }
    cur_vars <- all.vars(object$call$formula)
  }

  return(object)
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
