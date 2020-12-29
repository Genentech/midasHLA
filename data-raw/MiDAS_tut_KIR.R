#!/usr/bin/env R
devtools::load_all()

MiDAS_tut_KIR <-
  read.table("inst/extdata/MiDAS_tut_KIR.txt",
             header = TRUE,
             stringsAsFactors = FALSE)

usethis::use_data(MiDAS_tut_KIR, overwrite = TRUE)
