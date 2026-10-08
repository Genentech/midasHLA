#!/usr/bin/env R
# Script pre-parses alignment files for package use
# By Migdal
devtools::load_all()

# prepareHlaAlignment parses alignment and infers sequences of lower resolution
# alleles; this is the format of alignments shipped with the package
alignment_files <- list.files(path = "alignments", full.names = TRUE)
for (file in alignment_files) {
  cached_aln_obj <- prepareHlaAlignment(file)
  gene <- gsub(".*/([A-Z]+[0-9]*)_prot.txt", "\\1", file)
  saveRDS(cached_aln_obj,
          file = paste0(gene, "_prot.Rdata"))
}
