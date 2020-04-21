#!/usr/bin/env R
# Create Grantham distances vector
# each value gives Grantham distance for a pair of amino acids
# elements names are constructed by pasting the AA letters
# eg. "SS" or "SW"
path <- system.file("extdata", "grantham_distance_matrix.tsv", package = "MiDAS")
grantham_matrix <- read.table(path, header = TRUE, stringsAsFactors = FALSE)
dict_dist_grantham <- grantham_matrix %>%
  tidyr::gather(SECOND, SCORE, -FIRST)

# FIRST and SECOND contains all combinations of AA
# by using them as rownames efficient indexing will be a breeze
# and we can store it as a vector
idx <- paste(dict_dist_grantham[, "FIRST"], dict_dist_grantham[, "SECOND"], sep = "")
dict_dist_grantham <- setNames(dict_dist_grantham$SCORE, idx)

usethis::use_data(dict_dist_grantham, internal = TRUE, overwrite = TRUE)
