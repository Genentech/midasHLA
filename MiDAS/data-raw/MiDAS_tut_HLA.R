#!/usr/bin/env R
devtools::load_all()

MiDAS_tut_HLA <- readHlaCalls("inst/extdata/MiDAS_tut_HLA.txt", resolution = 8)

usethis::use_data(MiDAS_tut_HLA, overwrite = TRUE)
