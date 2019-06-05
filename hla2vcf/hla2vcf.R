#!/usr/bin/env Rscript
# HLA calls to VCF conversion
# By Migdal 06/2019
usage <-
"

Usage:
Rscript hla2vcf.R input resolution output

input input HLA calls file
resolution desierd HLA calls resolution
output output gziped vcf file

Example:
Rscript hla2vcf.R HLAHD_output_example.txt 8 output HLAHD_output_example.vcf.gz

"

suppressMessages(library("dplyr"))
suppressMessages(library("magrittr"))
suppressMessages(library("MiDAS"))
suppressMessages(library("stats"))
suppressPackageStartupMessages(library("vcfR"))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) {
  stop(usage)
}

hla_calls_path <- args[1]
resolution <- as.numeric(args[2])
vcf_out_file <- args[3]

hla_calls <- readHlaCalls(hla_calls_path, resolution = resolution)

genes <- hla_calls[-1] %>%
  colnames() %>%
  gsub("_.*", "", .) %>%
  unique()

alleles_by_gene <- lapply(genes, function(x) {
  x_pattern <- paste0(x, "_")
  hla_calls %>%
    select(starts_with(x_pattern)) %>%
    unlist() %>%
    unique() %>%
    na.omit() %>%
    as.character()
})
names(alleles_by_gene) <- genes

# filter empty genes
mask <- vapply(alleles_by_gene, function(x) length(x) > 0, FUN.VALUE = logical(1))
alleles_by_gene <- alleles_by_gene[mask]
genes <- names(alleles_by_gene)
if (length(genes) == 0) {
  stop("No alleles found in HLA calls file!")
}

# create metadata
meta <- c(
  "##fileformat=VCFv4.2",
  paste0("##fileDate=", format(Sys.time(), "%Y%m%d")),
  "##source=MiDAS",
  "##contig=<ID=chr6,length=170805979>",
  "##phasing=unphased",
  "##INFO=<ID=NS,Number=1,Type=Integer,Description=\"Number of Samples With Data\">",
  "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
  paste0('##ALT=<ID=INS,Description=\"', unlist(alleles_by_gene), '\">')
)

# create header
fix <- lapply(genes,
              FUN = function(x) {
                hlagene <- paste0("HLA-", x)
                alt <- alleles_by_gene[[x]] %>%
                  paste0("<", ., ">") %>%
                  paste(., collapse = ",")
                ns <- nrow(hla_calls) %>%
                  paste0("NS=", .)
                c("chr6", "31982024", hlagene, "A", alt, "50", "PASS", ns)
              }
)

if (length(genes) == 1) {
  fix <-
    matrix(fix,
           nrow = 1,
           ncol = 8,
           dimnames = list(
             NULL,
             c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO")
           ))
} else {
  fix <- do.call(rbind, fix)
  rownames(fix) <- NULL
  colnames(fix) <- c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO")
}

# create genotypes
gt <- lapply(genes, function(gene) {
  gene1 <- paste0(gene, "_1")
  allele1 <- hla_calls %>%
    select(starts_with(gene1)) %>%
    unlist() %>%
    match(alleles_by_gene[[gene]], nomatch = 0) %>%
    as.character()
  allele1[allele1 == "0"] <- "."

  gene2 <- paste0(gene, "_2")
  allele2 <- hla_calls %>%
    select(starts_with(gene2)) %>%
    unlist() %>%
    match(alleles_by_gene[[gene]], nomatch = 0) %>%
    as.character()
  allele2[allele2 == "0"] <- "."

  genotypes <- paste(allele1, allele2, sep = "/")
  c("GT", genotypes)
})

if (length(gt) == 1) {
  gt <-
    matrix(
      gt,
      nrow = 1,
      ncol = 1 + nrow(hla_calls),
      dimnames = list(NULL, c("FORMAT", hla_calls$ID))
    )
} else {
  gt <- do.call(rbind, gt)
  rownames(gt) <- NULL
  colnames(gt) <- c("FORMAT", hla_calls$ID)
}

vcf <- new(Class = "vcfR")
vcf@meta <- meta
vcf@fix <- fix
vcf@gt <- gt
check_keys(vcf)

write.vcf(vcf, file = vcf_out_file)
