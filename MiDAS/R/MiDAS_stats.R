#' Association analysis of HLA allele calls
#'
#' \code{analyzeHlaAssociations} performs associations analysis on HLA alleles
#' calls using statistical model of choice.
#'
#' @inheritParams checkHlaCallsFormat
#' @param model String specifying statistical model to use.
#' @param pheno Data frame holding phenotypic response variables.
#' @param covar Data frame holding covariates.
#' @param zygo Flag indicating wheater zygocity should be added to
#'   covariates. See details for further explanations.
#' @param reduce_counts Flag indicating wheater allele counts should be reduced
#'   to presence / absence indicators.
#' @param correction String specifying multiple testing correction method.
#' @param ... further arguments passed to statistical model function.
#'
#' Available choices for \code{model} can be checked in \link{hlaAssocModels}
#' documentation or by using \code{hlaAssocModels()} function.
#'
#' \code{pheno} and \code{covar} should be data frames with first column holding
#' samples IDs and named \code{ID}. Those should correspond to \code{ID} column
#' in \code{hla_calls}.
#'
#' \code{zygo} indicate if additional covariate, indicating sample zygocity
#' status, should be added to covariates. HLA allele counts for each sample
#' can take following values \code{0, 1, 2}. To avoid implying ordering on those
#' levels and effect size, this information can be split between two variables.
#' If \code{zygo} is set to \code{TRUE} zyocity variable is added during model
#' fitting, it specifies if sample is double homezygote for an allele.
#'
#' If \code{reduce_counts} is set to \code{TRUE} HLA allele counts are reduced
#' to presence / absence indicators. This is done by setting counts for double
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
#' hla_calls <- readHlaCalls("data/HLAexample.txt")
#' pheno <- read.table("data/pheno.txt",header=T)
#' covar <- read.table("data/covar.txt",header=T)
#'
#' # Cox proportional hazards regression model
#' analyzeHlaAssociations <- function(model = "coxph",
#'                                    hla_calls,
#'                                    pheno,
#'                                    covar,
#'                                    zygo = FALSE,
#'                                    reduce_counts = FALSE,
#'                                    correction = "BH"
#' )
#'
#' # Logistic regression
#' pheno <- pheno[, c(1, 3)]
#' analyzeHlaAssociations <- function(model = "glm.logit",
#'                                    hla_calls,
#'                                    pheno,
#'                                    covar,
#'                                    zygo = FALSE,
#'                                    reduce_counts = FALSE,
#'                                    correction = "BH"
#' )
#'
#' @importFrom assertthat assert_that see_if is.flag is.string
#' @importFrom broom tidy
#' @importFrom dplyr left_join filter mutate
#' @importFrom purrr map_dfr
#' @export
analyzeHlaAssociations <- function(model = "coxph",
                                   hla_calls,
                                   pheno,
                                   covar,
                                   zygo = TRUE,
                                   reduce_counts = TRUE,
                                   correction = "BH",
                                   ...) {
  assert_that(
    is.string(model),
    see_if(! is.null(hlaAssocModels(model)),
           msg = sprintf("%s model is not implemented.", model)),
    checkHlaCallsFormat(hla_calls),
    is.data.frame(pheno),
    see_if(nrow(pheno) >= 1 & ncol(pheno) >= 2,
           msg = "pheno have to have at least 1 rows and 2 columns"
    ),
    see_if(colnames(pheno)[1] == colnames(hla_calls)[1],
           msg = "first column in pheno must be named as first column in hla_calls"
    ),
    see_if(any(hla_calls[, 1] %in% pheno[, 1]),
           msg = "IDs in hla_calls doesn't match IDs in pheno"
    ),
    is.data.frame(covar),
    see_if(nrow(covar) >= 1 & ncol(covar) >= 2,
           msg = "covar have to have at least 1 rows and 2 columns"
    ),
    see_if(colnames(covar)[1] == colnames(hla_calls)[1],
           msg = "first column in covar must be named as first column in hla_calls"
    ),
    see_if(any(hla_calls[, 1] %in% covar[, 1]),
           msg = "IDs in hla_calls doesn't match IDs in covar"
    ),
    see_if(
      anyDuplicated(c(
        colnames(hla_counts[, -1]), colnames(pheno[, -1]),
        colnames(covar[, -1]), paste0(colnames(zygocity[, -1]), "_zygocity")
      )) == 0,
      msg = "some colnames in hla_calls and pheno and covar and zygocity are duplicated"
    ),
    is.flag(zygo),
    is.flag(reduce_counts),
    is.string("BH")
  )

  hla_counts <- hlaCallsToCounts(hla_calls)
  zygocity <- hla_counts

  if (reduce_counts) {
    hla_counts[, -1] <- lapply(hla_counts[, -1],
                               function(x) ifelse(x == 2, 1, x)
    )
  }

  data <- left_join(hla_counts, pheno, by="ID")
  data <- left_join(data, covar, by="ID")

  pheno_var <- paste(backquote(colnames(pheno)[-1]), collapse = ", ")
  covar_var <- paste(backquote(colnames(covar)[-1]), collapse = " + ")
  alleles_var <- backquote(colnames(hla_counts)[-1])

  if (zygo) {
    zygocity[, -1] <- lapply(zygocity[, -1], function(x) ifelse(x == 2, 1, 0))
    colnames(zygocity) <- c("ID", paste0(colnames(zygocity[, -1]), "_zygocity"))
    data <- left_join(data, zygocity, by = "ID")
    zygo_var <- paste0(colnames(hla_counts)[-1], "_zygocity")
    zygo_var <- backquote(zygo_var)
    alleles_var <- paste(alleles_var, zygo_var, sep = " + ")
  }

  model_function <- hlaAssocModels(model = model,
                                   response = pheno_var,
                                   covariate = covar_var,
                                   data = data
  )

  results <- map_dfr(
    .x = alleles_var,
    .f = ~tidy(model_function(., ...), exponentiate=TRUE) # this have to be handled somehow
  )

  results <- mutate(results, term = gsub("`", "", term))
  results <- filter(results, checkAlleleFormat(term))

  results <- mutate(results, p.adjusted = p.adjust(p.value, correction))

  return(results)
}

#' Association models for analysis of HLA alleles
#'
#' \code{hlaAssocModels} is a collection of preconfigured models for use with
#' HLA alleles count table.
#'
#' @inheritParams analyzeHlaAssociations
#' @param response Character specifying response variables in \code{data}.
#' @param covariate Chararcter specifying covariates in \code{data}.
#' @param data Data frame containing variables in the model.
#'
#' @return Function for fitting association model of choice, it takes allele
#'   number as an argument. Additional arguments can be passed as well.
#'
#' @examples
#'
#' @importFrom assertthat assert_that see_if
#' @export
hlaAssocModels <- function(model = NULL,
                           response,
                           covariate,
                           data) {
  if (is.null(model)) {
    print("Available models: coxph, lm, glm")
  }
  assert_that(
    is.string(model),
    is.character(response),
    is.character(covariate),
    is.data.frame(data)
  )
  model_function <- switch(
    model,
    coxph = function(allele, ...) {
      form <- sprintf("Surv(%s) ~ %s + %s", pheno_var, allele, covar_var)
      coxph(formula = as.formula(form), data = data, ...)
    },
    lm = function(allele, ...) {
      form <- sprintf("%s ~ %s + %s", pheno_var, allele, covar_var)
      lm(formula = as.formula(form), data = data, ...)
    },
    glm = function(allele, ...) {
      form <- sprintf("%s ~ %s + %s", pheno_var, allele, covar_var)
      glm(formula = as.formula(form), family = binomial(link = "logit"),  data = data, ...)
    }
  )
  return(model_function)
}
