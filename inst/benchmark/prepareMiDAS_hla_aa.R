#!/usr/bin/env R
# Migdal
# Run prepareMiDAS for hla_alleles

devtools::load_all()

argv = commandArgs(trailingOnly=TRUE)
input <- argv[1] # input hla_calls

hla_calls <- readHlaCalls(input)
col_data <- data.frame(
  ID = hla_calls$ID,
  val = rbinom(nrow(hla_calls), 1, 0.2))

midas <-
  prepareMiDAS(
    hla_calls = hla_calls,
    colData = col_data,
    experiment = "hla_aa"
  )

