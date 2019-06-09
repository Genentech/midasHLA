#!/usr/bin/env Rscript
# VCF to HLA calls conversion
# By Migdal 06/2019
usage <-
"

Usage:
Rscript vcf2hla.R input output

input input gziped vcf file; bgz files might not be fully supported
output output HLA calls file

Example:
Rscript vcf2hla.R HLAHD_output_example.vcf.bgz HLAHD_output_example.txt

"

suppressMessages(library("dplyr"))
suppressMessages(library("magrittr"))
suppressMessages(library("MiDAS"))
suppressMessages(library("stats"))
suppressMessages(library("tidyr"))
suppressMessages(library("vcfR")) # TODO rewrite using VariantAnnotation. vcfR doesn't support bgzip format which we would like to work on, this might have some impact on parsing vcf's however I've not observed it so far.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  stop(usage)
}

vcf_in_file <- args[1]
hla_out_file <- args[2]

vcf <- read.vcfR(vcf_in_file, verbose = FALSE)

vcf_tidy <- vcf %>%
  vcfR::extract_gt_tidy(verbose = FALSE) %>%
  dplyr::select(ID = Indiv, gt_GT_alleles) %>%
  dplyr::filter(! is.na(gt_GT_alleles)) %>%
  dplyr::mutate(gt_GT_alleles = gsub("[<>]", "", gt_GT_alleles)) %>%
  tidyr::separate(gt_GT_alleles, c("a1", "a2"), sep = "[|/]") %>%
  dplyr::filter(checkAlleleFormat(a1) | checkAlleleFormat(a2)) %>%
  tidyr::gather(a, allele, -ID) %>%
  dplyr::mutate(gene = paste0(gsub("\\*.*", "_", allele), gsub("a", "", a))) %>%
  dplyr::select(ID, gene, allele) %>%
  tidyr::spread(gene, allele) %>%
  as.data.frame(stringsAsFactors = FALSE)

write.table(vcf_tidy, file = hla_out_file, sep = "\t", quote = FALSE)
