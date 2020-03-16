#!/usr/bin/env R
# List of implemented MiDAS analysis types
midas_analysis_types <- data.frame(
  analysis_type = c("hla_allele", "aa_level", "expression_level",
                     "allele_g_group", "allele_supertype", "allele_group",
                     "kir_genes", "hla_kir_interactions", "none"
                     ),
  type = c("integer", "integer", "float", "integer", "integer", "integer",
            "integer", "integer", "unknown"
            ),
  stringsAsFactors = FALSE
)

usethis::use_data(midas_analysis_types, internal = TRUE, overwrite = TRUE)
