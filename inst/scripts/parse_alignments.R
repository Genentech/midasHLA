#!/usr/bin/env R
# Script pre-parses alignment files for package use
# By Migdal
devtools::load_all()
library(stringi)

alignment_files <- list.files(path = "alignments", full.names = TRUE)
for (file in alignment_files) {
  # parse aln files without any processing
  alignment <-
    readHlaAlignments(file,
                      trim = FALSE,
                      unkchar = "*")
  # infer missing lower resolution alleles
  for (res in c(6, 4)) {
    allele_numbers <- reduceAlleleResolution(rownames(alignment), resolution = res)
    missing_alleles <- unique(allele_numbers[! allele_numbers %in% rownames(alignment)])
    missing_aln <- list()
    for (allele in missing_alleles) {
      i <- allele_numbers == allele
      missing_aln[[allele]] <- apply(alignment[i, , drop=FALSE], 2, function(col) {
        if (all(col == col[1])) col[1]
        else "*"
      })
    }
    missing_aln <- do.call(rbind, missing_aln)
    alignment <- rbind(alignment, missing_aln)
  }

  # readHlaAlignments numbers columns according to IPD-IMGT/HLA nomenclature,
  # columns where the reference allele has gaps (insertions) are not numbered
  first_codon_idx <- which(colnames(alignment) == "1")

  cached_aln_obj <- list(alignment, first_codon_idx)
  gene <- gsub(".*/([A-Z]+[0-9]*)_prot.txt", "\\1", file)
  saveRDS(cached_aln_obj,
          file = paste0(gene, "_prot.Rdata"))
}
