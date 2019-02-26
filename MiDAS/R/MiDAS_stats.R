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
#' hla_calls_file <- system.file("extdata", "HLAHD_output_example.txt", package = "MiDAS")
#' hla_calls <- readHlaCalls("data/HLAexample.txt")
#' pheno_file <- system.file("extdata", "pheno_example.txt", package = "MiDAS")
#' pheno <- read.table(pheno_file, header = TRUE)
#' covar_file <- system.file("extdata", "covar_example.txt", package = "MiDAS")
#' covar <- read.table(covar_file, header = TRUE)
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
                                   zygo = FALSE,
                                   reduce_counts = FALSE,
                                   correction = "BH") {
  assert_that(
    is.string(model),
    see_if(model %in% hlaAssocModels(),
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
    is.flag(zygo),
    is.flag(reduce_counts),
    is.string(correction)
  )

  hla_counts <- hlaCallsToCounts(hla_calls)
  zygocity <- hla_counts

  assert_that(
    see_if(
    anyDuplicated(c(
      colnames(hla_counts[, -1]), colnames(pheno[, -1]),
      colnames(covar[, -1]), paste0(colnames(zygocity[, -1]), "_zygocity")
    )) == 0,
    msg = "some colnames in hla_calls and pheno and covar and zygocity are duplicated"
  ))

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
    .f = ~tidy(model_function(.), exponentiate=FALSE) # TODO exponentiate could be passed somehow
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
#' \code{hlaAssocModels} is not indended to use by basic user.
#'
#' @inheritParams analyzeHlaAssociations
#' @param response Character specifying response variables in \code{data}.
#' @param covariate Chararcter specifying covariates in \code{data}.
#' @param data Data frame containing variables in the model.
#'
#' @return Function for fitting association model of choice, it takes allele
#'   number as an argument. Additional arguments can be passed as well.
#'
#'   Returned function takes one required argument: HLA allele number
#'   (\code{allele}). Additional parameters are passed to statistical model
#'   function. Due to fact that allele numbers contains characters that have
#'   special meanings in formulas, they should be backquoted. This can be easly
#'   done with \link{backquote}. See examples section for general usage case.
#'
#'   If model is equal to \code{NULL} names of available models are returned
#'   instead.
#'
#' @examples
#' hla_calls_file <- system.file("extdata", "HLAHD_output_example.txt", package = "MiDAS")
#' hla_calls <- readHlaCalls("data/HLAexample.txt")
#' hla_counts <- hlaCallsToCounts(hla_calls)
#' pheno_file <- system.file("extdata", "pheno_example.txt", package = "MiDAS")
#' pheno <- read.table(pheno_file, header = TRUE)
#' covar_file <- system.file("extdata", "covar_example.txt", package = "MiDAS")
#' covar <- read.table(covar_file, header = TRUE)
#' data <- left_join(hla_counts, pheno, by="ID")
#' data <- left_join(data, covar, by="ID")
#' response <- paste(colnames(pheno[, -1]), collapse = ", ")
#' covariate <- paste(colnames(covar[, -1]), collapse = " + ")
#' func <- hlaAssocModels(model = "coxph",
#'                        response = response,
#'                        covariate = covariate,
#'                        data = data
#' )
#' allele <- backquote("A*01:01")
#' func(allele)
#'
#' @importFrom assertthat assert_that see_if
#' @importFrom stats binomial glm lm
#' @importFrom survival coxph Surv
#' @export
hlaAssocModels <- function(model = NULL,
                           response,
                           covariate,
                           data) {
  if (is.null(model)) {
    return(c("coxph", "lm", "glm.logit"))
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
      form <- sprintf("Surv(%s) ~ %s + %s", response, allele, covariate)
      result <- coxph(formula = as.formula(form), data = data, ...)
      return(result)
    },
    lm = function(allele, ...) {
      form <- sprintf("%s ~ %s + %s", response, allele, covariate)
      result <- lm(formula = as.formula(form), data = data, ...)
      return(result)
    },
    glm.logit = function(allele, ...) {
      form <- sprintf("%s ~ %s + %s", response, allele, covariate)
      result <- glm(
        formula = as.formula(form),
        family = binomial(link = "logit"),
        data = data,
        ...
      )
      return(result)
    }
  )
  return(model_function)
}




hla_calls <- readHlaCalls("inst/extdata/HLAHD_output_example.txt")
pheno_file <- system.file("extdata", "pheno_example.txt", package = "MiDAS")
pheno <- read.table(pheno_file, header = TRUE)
covar <- read.table("inst/extdata/covar_example.txt",header=T)
res <- analyzeHlaAssociations(hla_calls = hla_calls, pheno = pheno, covar = covar, zygo = F, correction = "none")

HLA_data <- hla_calls
tmp <- HLA_data[, -1]
HLA_toAnalyze <- cbind(tmp[0], mtabulate(as.data.frame(t(tmp))))
row.names(HLA_toAnalyze) <- NULL
HLA_toAnalyze <- HLA_toAnalyze[ , order(names(HLA_toAnalyze))]
HLA_toAnalyze <- cbind(ID=(HLA_data[,1]),HLA_toAnalyze)
tmp <- merge(HLA_toAnalyze, pheno, by="ID")
dat <- merge(tmp, covar, by="ID")
coxfun <- function(x) {
  coxph( Surv(OS, OS_DIED) ~ x + AGE + SEX, data=dat)
}
firstHLA <- which(names(dat)=="A*01:01")
lastHLA <- which(names(dat)=="L*01:01")
suppressWarnings(HLA_4digit_OS_res <- map_dfr(dat[,firstHLA:lastHLA], ~coxfun(.x) %>% tidy(exponentiate=TRUE)))
HLA_4digit_summary <- HLA_4digit_OS_res %>% filter(term=="x")
HLA_4digit_summary$allele <- names(dat[firstHLA:lastHLA])

res %>% filter(p.value<0.05)
HLA_4digit_summary[c(8,2,3,4,6,7,5)] %>% filter(p.value<0.05)

all(HLA_4digit_summary$p.value == res$p.value, na.rm = T)








