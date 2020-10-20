#!/usr/bin/env R
devtools::load_all()

MiDAS_tut_pheno <-
  read.table("inst/extdata/MiDAS_tut_pheno.txt",
             header = TRUE,
             stringsAsFactors = FALSE)

usethis::use_data(MiDAS_tut_pheno, overwrite = TRUE)
