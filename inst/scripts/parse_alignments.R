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
                      unkchar = "*",
                      resolution = 8)
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

  # find first codon idx
  aln_raw <- stri_read_lines(file)
  aln <- stri_split_regex(aln_raw, "\\s+")
  nonempty_lines <- vapply(aln, length, integer(length = 1)) >= 2
  allele_numbers <- vapply(aln, `[`, character(length = 1), 2)
  allele_lines <- checkAlleleFormat(allele_numbers)
  aln_raw <- aln_raw[nonempty_lines]
  raw_first_codon_idx <- nchar(stri_subset_fixed(aln_raw, "Prot")[1])
  raw_alignment_line <- stri_sub(aln_raw[allele_lines][1],
                                 1,
                                 raw_first_codon_idx
  )
  raw_alignment_seq <- stri_split_regex(raw_alignment_line, "\\s+")
  raw_alignment_seq <- unlist(raw_alignment_seq)[-c(1, 2)]
  first_codon_idx <- nchar(stri_flatten(raw_alignment_seq))
  
  # find AA positions numbers
  if (first_codon_idx > 1) {
    aln_colnames <- c(seq(1 - first_codon_idx, -1, 1),
                      seq(1, ncol(alignment) + 1 - first_codon_idx, 1)
    )
  } else {
    aln_colnames <- seq(1, ncol(alignment) + 1 - first_codon_idx, 1)
  }
  colnames(alignment) <- aln_colnames

  cached_aln_obj <- list(alignment, first_codon_idx)
  gene <- gsub(".*/([A-Z]+[0-9]*)_prot.txt", "\\1", file)
  saveRDS(cached_aln_obj,
          file = paste0(gene, "_prot.Rdata"))
}
