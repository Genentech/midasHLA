#!/usr/bin/env R
devtools::load_all()

MiDAS_tut_object <- prepareMiDAS(
  hla_calls = reduceHlaCalls(MiDAS_tut_HLA, resolution = 4),
  kir_calls = MiDAS_tut_KIR,
  colData = MiDAS_tut_pheno,
  experiment = c(
    "hla_alleles",
    "hla_aa",
    "hla_g_groups",
    "hla_supertypes",
    "hla_NK_ligands",
    "kir_genes",
    "kir_haplotypes",
    "hla_kir_interactions",
    "hla_divergence",
    "hla_het"
  )
)
usethis::use_data(MiDAS_tut_object, overwrite = TRUE)
