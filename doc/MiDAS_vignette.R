## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE)

library("dplyr")
library("kableExtra")
library("knitr")
library("MiDAS")

## ----load_pheno, echo=TRUE, warning=FALSE-------------------------------------
pheno_file <- system.file("extdata", "MiDAS_tut_pheno.txt", package = "MiDAS")
pheno <- read.table(pheno_file, header = TRUE, stringsAsFactors = FALSE)

## ----show_pheno, echo=FALSE, warning=FALSE------------------------------------
pheno %>%
  head(10) %>%
  kable() %>%
  kable_styling(bootstrap_options = c("striped", "hover", "condensed"))

## ----load_hla, echo=TRUE, warning=FALSE---------------------------------------
# HLA calls can be loaded using the readHlaCalls function with the desired resolution
hla_calls_file <- system.file("extdata", "MiDAS_tut_HLA.txt", package = "MiDAS")
hla_calls <- readHlaCalls(hla_calls_file, resolution = 4)

## ----show_hla, echo=FALSE, warning=FALSE--------------------------------------
hla_calls %>%
  head(10) %>%
  kable() %>%
  kable_styling(bootstrap_options = c("striped", "hover", "condensed")) %>%
  scroll_box(width = "100%", height = "200px")

## ----load_kir, echo=TRUE, warning=FALSE---------------------------------------
# KIR calls (currently presence/absence calls, no allele-level resolution) can be loaded using the readKirCalls function
kir_calls_file <- system.file("extdata", "MiDAS_tut_KIR.txt", package = "MiDAS")
kir_calls <- readKirCalls(kir_calls_file)

## ----show_kir, echo=FALSE, warning=FALSE--------------------------------------
kir_calls %>%
  head(10) %>%
  kable() %>%
  kable_styling(bootstrap_options = c("striped", "hover", "condensed")) %>%
  scroll_box(width = "100%", height = "200px")

## ----test_data, echo=TRUE, warning=FALSE--------------------------------------
HLA <- reduceHlaCalls(MiDAS_tut_HLA, resolution = 4)
KIR <- MiDAS_tut_KIR
pheno <- MiDAS_tut_pheno

## ----creating_midas_data_set, echo=TRUE, warning=FALSE------------------------
midas <- prepareMiDAS(
  hla_calls = HLA,
  kir_calls = KIR,
  colData = pheno,
  experiment = c("hla_alleles", "hla_aa")
  )

## ----get_freq, echo=TRUE, warning=FALSE---------------------------------------
freq <- getFrequencies(
  object = midas,
  carrier_frequency = FALSE,
  experiment = "hla_alleles",
  compare = TRUE
)

## ----show_freq, echo=FALSE, warning=FALSE-------------------------------------
kable(freq) %>%
  kable_styling(bootstrap_options = c("striped", "hover", "condensed")) %>%
  scroll_box(width = "100%", height = "200px")

## ----filter_freq, echo=TRUE, warning=FALSE------------------------------------
midas <-
  filterByFrequency(
    object = midas,
    experiment = "hla_alleles",
    lower_frequency_cutoff = 0.01
    )

## ----model definition, echo=TRUE, warning=FALSE-------------------------------
# Logistic regression
object <- glm(disease ~ term, data = midas, family = binomial(link = "logit"))

## ----association_analysis, echo=TRUE, warning=FALSE---------------------------
results <-
  runMiDAS(
    object = object,
    experiment = "hla_alleles",
    conditional = FALSE,
    omnibus = FALSE,
    lower_frequency_cutoff = 0.05
    )

kableResults(results)

