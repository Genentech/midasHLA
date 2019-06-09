#!/usr/bin/env Rscript
# HLA calls to VCF conversion
# By Migdal 06/2019
usage <-
"

Usage:
Rscript vcf2hla.R input output

input input gziped vcf file
output output HLA calls file

Example:
Rscript vcf2hla.R HLAHD_output_example.vcf.gz HLAHD_output_example.txt

"

suppressMessages(library("dplyr"))
suppressMessages(library("magrittr"))
suppressMessages(library("MiDAS"))
suppressMessages(library("stats"))
suppressPackageStartupMessages(library("vcfR"))
suppressMessages(library("tidyr"))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  stop(usage)
}

vcf_in_file <- args[1]
hla_out_file <- args[2]

vcf <- read.vcfR(vcf_in_file, verbose = FALSE)

vcf_tidy <- vcfR::extract_gt_tidy(vcf, verbose = FALSE) %>%
  dplyr::select(ID = Indiv, gt_GT_alleles) %>%
  dplyr::filter(! is.na(gt_GT_alleles)) %>%
  dplyr::mutate(gt_GT_alleles = gsub("[<>]", "", .data$gt_GT_alleles)) %>%
  tidyr::separate(gt_GT_alleles, c("a1", "a2"), sep = "[|/]") %>%
  dplyr::filter(checkAlleleFormat(a1) | checkAlleleFormat(a2)) %>%
  tidyr::gather(a, allele, -ID) %>%
  dplyr::mutate(gene = paste0(gsub("\\*.*", "_", allele), gsub("a", "", a))) %>%
  dplyr::select(ID, gene, allele) %>%
  tidyr::spread(gene, allele) %>%
  dplyr::arrange(ID) %>%
  as.data.frame(stringsAsFactors = FALSE)

write.table(vcf_tidy, file = "test_out.txt", sep = "\t", quote = FALSE)
