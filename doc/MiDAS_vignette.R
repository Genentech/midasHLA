## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

devtools::load_all()
library("dplyr")
library("knitr")
library("kableExtra")

## ----loading_data-------------------------------------------------------------
coldata_file <- system.file("extdata", "MiDAS_tut_pheno.txt", package = "MiDAS")
coldata <- read.table(coldata_file, header = TRUE, stringsAsFactors = FALSE)

kable(coldata) %>%
  kable_styling(bootstrap_options = c("striped", "hover", "condensed")) %>%
  scroll_box(width = "100%", height = "200px")

# hla calls can be loaded using readHlaCalls function
hla_calls_file <- system.file("extdata", "MiDAS_tut_HLA.txt", package = "MiDAS")
hla_calls <- readHlaCalls(hla_calls_file)

kable(hla_calls) %>%
  kable_styling(bootstrap_options = c("striped", "hover", "condensed")) %>%
  scroll_box(width = "100%", height = "200px")

## ----creating_midas_data_set--------------------------------------------------
midas <- prepareMiDAS(
  hla_calls = hla_calls,
  colData = coldata,
  experiment = "hla_alleles"
  )

## ----midas_obj----------------------------------------------------------------
midas

## ----model definition---------------------------------------------------------
# Logistic regression
object <- glm(disease ~ term, data = midas, family = binomial(link = "logit"))

## ----association_analysis-----------------------------------------------------
results <-
  runMiDAS(
    object = object,
    experiment = "hla_alleles",
    conditional = FALSE,
    omnibus = FALSE
    )

kableResults(results)

